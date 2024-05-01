import logging

import pandas as pd
from alchemlyb.estimators import MBAR
from alchemlyb.parsing.gmx import extract_u_nk
from flow import aggregator

from martignac import config
from martignac.nomad.workflows import build_nomad_workflow, build_nomad_workflow_with_multiple_jobs
from martignac.utils.gromacs import gromacs_simulation_command
from martignac.utils.martini_flow_projects import (
    Job,
    MartiniFlowProject,
    fetched_from_nomad,
    import_job_from_other_flow,
    store_gromacs_log_to_doc_with_state_point,
    store_task_for_many_jobs,
    store_workflow_for_many_jobs,
    uploaded_to_nomad,
)
from martignac.utils.misc import func_name, sub_template_mdp, update_nested_dict
from martignac.workflows.system_generation.solute import SoluteGenFlow, get_solute_job
from martignac.workflows.system_generation.solute_in_solvent import SoluteInSolventFlow, get_solvation_job
from martignac.workflows.system_generation.solvent import SolventGenFlow, get_solvent_job

logger = logging.getLogger(__name__)

conf = config()["alchemical_transformation"]


class AlchemicalTransformationFlow(MartiniFlowProject):
    workspace_path: str = f"{MartiniFlowProject.workspaces_path}/{conf['relative_paths']['workspaces']}"
    itp_path = f"{MartiniFlowProject.input_files_path}/{conf['relative_paths']['itp_files']}"
    mdp_path = f"{MartiniFlowProject.input_files_path}/{conf['relative_paths']['mdp_files']}"
    itp_files = {k: v.get(str) for k, v in conf["itp_files"].items()}
    mdp_files = {k: v.get(str) for k, v in conf["mdp_files"].items()}
    simulation_settings = {
        "n_threads": conf["settings"]["n_threads"].get(int),
        "temperature": conf["settings"]["temperature"].get(float),
    }
    system_name = conf["output_names"]["system"].get(str)
    nomad_workflow: str = conf["output_names"]["nomad_workflow"].get(str)
    nomad_top_level_workflow: str = conf["output_names"]["nomad_top_level_workflow"].get(str)
    state_names = {k: v.get(str) for k, v in conf["output_names"]["states"].items()}

    @classmethod
    def get_system_run_name(cls, lambda_state: str) -> str:
        return AlchemicalTransformationFlow.get_state_name(state_name=f"production-{lambda_state}")


project_name = AlchemicalTransformationFlow.class_name()


def solvent_and_solute_aggregator(jobs):
    unique_combinations = {(job.sp["solvent_name"], job.sp["solute_name"]) for job in jobs}
    for solvent, solute in unique_combinations:
        yield [job for job in jobs if job.sp.solvent_name == solvent and job.sp.solute_name == solute]


def job_system_keys(job) -> [str, str]:
    return [job.sp.solvent_name, job.sp.solute_name]


def lowest_lambda_job(jobs) -> Job:
    return min(jobs, key=lambda j: j.sp.lambda_state)


@AlchemicalTransformationFlow.label
def system_prepared(job) -> bool:
    if project_name in job.document:
        return job.doc[project_name].get("system_prepared", False)
    return False


@AlchemicalTransformationFlow.label
def all_jobs_prepared(*jobs) -> bool:
    return all(system_prepared(j) for j in jobs)


@AlchemicalTransformationFlow.label
def nomad_workflow_built(*jobs) -> bool:
    return all(
        [
            any(job.isfile(AlchemicalTransformationFlow.nomad_workflow) for job in jobs),
            any(job.isfile(AlchemicalTransformationFlow.nomad_top_level_workflow) for job in jobs),
        ]
    )


@AlchemicalTransformationFlow.label
def any_uploaded_to_nomad(*jobs) -> bool:
    return any(uploaded_to_nomad(job) for job in jobs)


@AlchemicalTransformationFlow.pre(fetched_from_nomad)
@AlchemicalTransformationFlow.post(all_jobs_prepared, tag="system_prepared")
@AlchemicalTransformationFlow.operation_hooks.on_success(store_workflow_for_many_jobs)
@AlchemicalTransformationFlow.operation(aggregator=aggregator(aggregator_function=solvent_and_solute_aggregator))
def prepare_system(*jobs):
    job = lowest_lambda_job(jobs)
    logger.info(f"Preparing system for {job.id} @ AlchemicalTransformationFlow")
    solute_job: Job = get_solute_job(job)
    solute_keys = ["solute_itp"]
    import_job_from_other_flow(job, solute_gen_project, solute_job, solute_keys)
    solvent_job: Job = get_solvent_job(job)
    import_job_from_other_flow(job, solvent_gen_project, solvent_job, [])
    solvation_job: Job = get_solvation_job(job)
    solvation_keys = ["solute_solvent_gro", "solute_solvent_top"]
    project.operation_to_workflow[func_name()] = solute_solvation_project.class_name()
    import_job_from_other_flow(job, solute_solvation_project, solvation_job, solvation_keys)

    for other_job in [j for j in jobs if j != job]:
        logger.info(f"Importing data to job {other_job.id} @ AlchemicalTransformationFlow")
        import_job_from_other_flow(other_job, solute_gen_project, solute_job, solute_keys, run_child_job=False)
        import_job_from_other_flow(other_job, solvent_gen_project, solvent_job, [], run_child_job=False)
        import_job_from_other_flow(
            other_job, solute_solvation_project, solvation_job, solvation_keys, run_child_job=False
        )

    for job in jobs:
        job.doc = update_nested_dict(job.doc, {project_name: {"system_prepared": True}})


@AlchemicalTransformationFlow.label
def lambda_sampled(job):
    if project_name not in job.doc:
        return False
    log_file = job.doc[project_name].get("alchemical_log", None)
    if log_file and job.isfile(log_file):
        with open(job.fn(log_file)) as file:
            lines = file.read()
        return "Finished mdrun on rank" in lines
    return False


@AlchemicalTransformationFlow.label
def all_lambda_states_sampled(*jobs):
    result = all(
        [
            job.doc.get(project_name, False)
            and job.isfile(job.fn(job.doc[project_name].get("alchemical_xvg", default="")))
            for job in jobs
        ]
    )
    return result


@AlchemicalTransformationFlow.label
def free_energy_already_calculated(*jobs):
    return all(
        [
            job.doc.get(project_name, False)
            and job.doc[project_name].get("free_energy", False)
            and job.doc[project_name]["free_energy"].get("mean", False)
            for job in jobs
        ]
    )


@AlchemicalTransformationFlow.pre(system_prepared, tag="system_prepared")
@AlchemicalTransformationFlow.post(lambda_sampled, tag="lambda_sampled")
@AlchemicalTransformationFlow.operation_hooks.on_success(store_gromacs_log_to_doc_with_state_point)
@AlchemicalTransformationFlow.operation(cmd=True, with_job=True)
def production(job):
    lambda_state = str(job.sp.lambda_state)
    lambda_mdp = f"{AlchemicalTransformationFlow.get_system_run_name(job.sp.lambda_state)}.mdp"
    sub_template_mdp(AlchemicalTransformationFlow.mdp_files.get("production"), "sedstate", lambda_state, lambda_mdp)
    sub_template_mdp(lambda_mdp, "sedmolecule", job.sp["solute_name"], lambda_mdp)
    job.doc[project_name]["alchemical_xvg"] = f"{AlchemicalTransformationFlow.get_system_run_name(lambda_state)}.xvg"
    job.doc[project_name]["alchemical_log"] = f"{AlchemicalTransformationFlow.get_system_run_name(lambda_state)}.log"
    return gromacs_simulation_command(
        mdp=lambda_mdp,
        top=job.doc[solute_solvation_project.class_name()].get("solute_solvent_top"),
        gro=job.doc[solute_solvation_project.class_name()].get("solute_solvent_gro"),
        name=AlchemicalTransformationFlow.get_system_run_name(lambda_state),
        n_threads=AlchemicalTransformationFlow.simulation_settings.get("n_threads"),
    )


@AlchemicalTransformationFlow.pre(all_lambda_states_sampled, tag="lambda_sampled")
@AlchemicalTransformationFlow.post(free_energy_already_calculated, tag="mbar")
@AlchemicalTransformationFlow.operation_hooks.on_success(store_task_for_many_jobs)
@AlchemicalTransformationFlow.operation(
    aggregator=aggregator(aggregator_function=solvent_and_solute_aggregator, sort_by="lambda_state")
)
def compute_free_energy(*jobs):
    logger.info("calculating free energies using pyMBAR}")
    xvg_files = [job.fn(job.doc[project_name]["alchemical_xvg"]) for job in jobs]
    u_nk_list = [
        extract_u_nk(f, T=AlchemicalTransformationFlow.simulation_settings.get("temperature")) for f in xvg_files
    ]
    u_nk_combined = pd.concat(u_nk_list)
    mbar = MBAR().fit(u_nk_combined)
    free_energy = float(mbar.delta_f_.iloc[0, -1])
    d_free_energy = float(mbar.d_delta_f_.iloc[0, -1])
    solvent_name, solute_name = job_system_keys(jobs[0])
    logger.info(f"free energy {solvent_name}_{solute_name}: {free_energy:+6.2f} +/- {d_free_energy:5.2f} kT")
    for job in jobs:
        job.doc = update_nested_dict(
            job.doc, {project_name: {"free_energy": {"mean": free_energy, "std": d_free_energy}}}
        )


@AlchemicalTransformationFlow.pre(free_energy_already_calculated, tag="mbar")
@AlchemicalTransformationFlow.post(nomad_workflow_built)
@AlchemicalTransformationFlow.operation(aggregator=aggregator(aggregator_function=solvent_and_solute_aggregator))
def generate_nomad_workflow(*jobs):
    for job in jobs:
        with job:
            build_nomad_workflow(job, is_top_level=False, add_job_id=True)
    job = lowest_lambda_job(jobs)
    with job:
        build_nomad_workflow_with_multiple_jobs(project, list(jobs))


@AlchemicalTransformationFlow.pre(nomad_workflow_built)
@AlchemicalTransformationFlow.post(any_uploaded_to_nomad)
@AlchemicalTransformationFlow.operation(aggregator=aggregator(aggregator_function=solvent_and_solute_aggregator))
def upload_to_nomad(*jobs):
    return AlchemicalTransformationFlow.upload_to_nomad_multiple_jobs(list(jobs))


project = AlchemicalTransformationFlow.get_project(path=AlchemicalTransformationFlow.workspace_path)
solute_gen_project = SoluteGenFlow.get_project(path=SoluteGenFlow.workspace_path)
solvent_gen_project = SolventGenFlow.get_project(path=SolventGenFlow.workspace_path)
solute_solvation_project = SoluteInSolventFlow.get_project(path=SoluteInSolventFlow.workspace_path)

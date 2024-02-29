import logging

import networkx as nx
import pandas as pd
from alchemlyb.estimators import MBAR
from alchemlyb.parsing.gmx import extract_u_nk
from flow import aggregator

from martignac import config
from martignac.nomad.workflows import NomadTopLevelWorkflow, NomadWorkflow
from martignac.utils.gromacs import gromacs_simulation_command
from martignac.utils.martini_flow_projects import MartiniFlowProject
from martignac.utils.misc import copy_files_to, sub_template_mdp
from martignac.workflows.solute_solvation import SoluteSolvationFlow

logger = logging.getLogger(__name__)

conf = config()["alchemical_transformation"]


class AlchemicalTransformationFlow(MartiniFlowProject):
    run_mdp: str = conf["mdp_files"]["run"].get(str)
    n_threads: int = conf["settings"]["n_threads"].get(int)
    temperature: float = conf["settings"]["temperature"].get(float)
    system_name: str = conf["output_names"]["system"].get(str)
    nomad_workflow: str = conf["output_names"]["nomad_workflow"].get(str)
    nomad_top_level_workflow: str = conf["output_names"]["nomad_top_level_workflow"].get(str)
    run_name: str = conf["output_names"]["states"]["run"].get(str)

    @classmethod
    def get_system_run_name(cls, lambda_state: str) -> str:
        return f"{cls.system_name}_{cls.run_name}-{lambda_state}"

    @classmethod
    def get_system_run_gro(cls, lambda_state: str) -> str:
        return f"{cls.get_system_run_name(lambda_state)}.gro"


project = AlchemicalTransformationFlow().get_project()


def solvent_and_solute_aggregator(jobs):
    unique_combinations = {(job.sp["solvent_name"], job.sp["solute_name"]) for job in jobs}
    for solvent, solute in unique_combinations:
        yield [job for job in jobs if job.sp.solvent_name == solvent and job.sp.solute_name == solute]


def job_system_keys(job) -> [str, str]:
    return [job.sp.solvent_name, job.sp.solute_name]


def lowest_lambda_job(jobs) -> bool:
    return min(jobs, key=lambda j: j.sp.lambda_state)


@AlchemicalTransformationFlow.label
def system_prepared(job) -> bool:
    return job.document.get("system_prepared", False)


@AlchemicalTransformationFlow.label
def all_jobs_prepared(*jobs) -> bool:
    return all(system_prepared(j) for j in jobs)


@AlchemicalTransformationFlow.label
def all_nomad_workflows_built(*jobs) -> bool:
    return all(job.isfile(AlchemicalTransformationFlow.nomad_workflow) for job in jobs)


@AlchemicalTransformationFlow.label
def nomad_top_level_workflow_built(*jobs) -> bool:
    for job in jobs:
        if job == lowest_lambda_job(jobs):
            return job.isfile(AlchemicalTransformationFlow.nomad_top_level_workflow)


@AlchemicalTransformationFlow.post(all_jobs_prepared)
@AlchemicalTransformationFlow.operation(aggregator=aggregator(aggregator_function=solvent_and_solute_aggregator))
def prepare_system(*jobs):
    for job in jobs:
        if job == lowest_lambda_job(jobs):
            SoluteSolvationFlow().run(jobs=[job])
            for other_job in jobs:
                if "gromacs_logs" in job.document:
                    if "gromacs_logs" not in other_job.document:
                        other_job.document["gromacs_logs"] = {}
                    other_job.document["gromacs_logs"] = dict(other_job.document["gromacs_logs"]) | dict(
                        job.document["gromacs_logs"]
                    )
                if "nomad_workflow" in job.document:
                    if "nomad_workflow" not in other_job.document:
                        other_job.document["nomad_workflow"] = {}
                    other_job.document["nomad_workflow"] = dict(other_job.document["nomad_workflow"]) | dict(
                        job.document["nomad_workflow"]
                    )
                other_job.document["system_gro"] = SoluteSolvationFlow.get_system_equ_gro()
                other_job.document["system_top"] = SoluteSolvationFlow.get_top()
                other_job.document["solute_name"] = job.document["solute_name"]
                other_job.document["solute_itp"] = job.document["solute_itp"]
                other_job.document["system_prepared"] = True
                copy_files_to(
                    [
                        job.fn(other_job.document["system_gro"]),
                        job.fn(other_job.document["system_top"]),
                        job.fn(other_job.document["solute_itp"]),
                    ],
                    other_job.path,
                )


@AlchemicalTransformationFlow.label
def lambda_sampled(job):
    log_file = job.document.get("alchemical_log", None)
    if log_file and job.isfile(log_file):
        with open(job.fn(log_file)) as file:
            lines = file.read()
        return "Finished mdrun on rank" in lines
    return False


@AlchemicalTransformationFlow.label
def all_lambda_states_sampled(*jobs):
    result = all([job.isfile(job.fn(job.document.get("alchemical_xvg", default=""))) for job in jobs])
    return result


@AlchemicalTransformationFlow.label
def free_energy_already_calculated(*jobs):
    solvent_name, solute_name = job_system_keys(jobs[0])
    return any(
        [
            e
            for e in project.document.get("free_energies", [])
            if e["solute_name"] == solute_name and e["solvent_name"] == solvent_name
        ]
    )


@AlchemicalTransformationFlow.pre(system_prepared)
@AlchemicalTransformationFlow.post(lambda_sampled)
@AlchemicalTransformationFlow.operation(cmd=True, with_job=True)
@AlchemicalTransformationFlow.log_gromacs_simulation(AlchemicalTransformationFlow.get_system_run_name(""), True)
def sample_lambda(job):
    lambda_state = str(job.sp.lambda_state)
    lambda_mdp = job.fn("run_lambda.mdp")
    sub_template_mdp(AlchemicalTransformationFlow.run_mdp, "sedstate", lambda_state, lambda_mdp)
    sub_template_mdp(lambda_mdp, "sedmolecule", job.document["solute_name"], lambda_mdp)
    job.document["alchemical_xvg"] = f"{AlchemicalTransformationFlow.get_system_run_name(job.sp.lambda_state)}.xvg"
    job.document["alchemical_log"] = f"{AlchemicalTransformationFlow.get_system_run_name(job.sp.lambda_state)}.log"
    return gromacs_simulation_command(
        mdp=lambda_mdp,
        top=job.document.get("system_top"),
        gro=job.document.get("system_gro"),
        name=AlchemicalTransformationFlow.get_system_run_name(lambda_state),
        n_threads=AlchemicalTransformationFlow.n_threads,
    )


@AlchemicalTransformationFlow.pre(all_lambda_states_sampled)
@AlchemicalTransformationFlow.post(free_energy_already_calculated)
@AlchemicalTransformationFlow.operation(
    aggregator=aggregator(aggregator_function=solvent_and_solute_aggregator, sort_by="lambda_state")
)
def compute_free_energy(*jobs):
    xvg_files = [job.fn(job.document["alchemical_xvg"]) for job in jobs]
    u_nk_list = [extract_u_nk(f, T=AlchemicalTransformationFlow.temperature) for f in xvg_files]
    u_nk_combined = pd.concat(u_nk_list)
    mbar = MBAR().fit(u_nk_combined)
    free_energy = float(mbar.delta_f_.iloc[0, -1])
    d_free_energy = float(mbar.d_delta_f_.iloc[0, -1])
    solvent_name, solute_name = job_system_keys(jobs[0])
    logger.info(f"free energy {solvent_name}_{solute_name}: {free_energy:+6.2f} +- {d_free_energy:5.2f} kT")
    f_contrib = {
        "solvent_name": solvent_name,
        "solute_name": solute_name,
        "f_mean": free_energy,
        "f_std": d_free_energy,
    }
    if "free_energies" not in project.document:
        project.document["free_energies"] = []
    project.document["free_energies"].append(f_contrib)


@AlchemicalTransformationFlow.pre(free_energy_already_calculated)
@AlchemicalTransformationFlow.post(lambda job: job.isfile(AlchemicalTransformationFlow.nomad_workflow))
@AlchemicalTransformationFlow.operation(with_job=True)
def generate_nomad_workflow(job):
    workflow = NomadWorkflow(project, job)
    workflow.build_workflow_yaml(AlchemicalTransformationFlow.nomad_workflow)


@AlchemicalTransformationFlow.pre(all_nomad_workflows_built)
@AlchemicalTransformationFlow.post(nomad_top_level_workflow_built)
@AlchemicalTransformationFlow.operation(aggregator=aggregator(aggregator_function=solvent_and_solute_aggregator))
def generate_nomad_top_level_workflow(*jobs):
    for job in jobs:
        if job == lowest_lambda_job(jobs):
            graph = nx.DiGraph([("SoluteSolvationFlow", "AlchemicalTransformationFlow")])
            projects = {
                "SoluteSolvationFlow": SoluteSolvationFlow.get_project(),
                "AlchemicalTransformationFlow": project,
            }
            top_level_workflow = NomadTopLevelWorkflow(projects, job, graph)
            top_level_workflow.build_workflow_yaml(AlchemicalTransformationFlow.nomad_top_level_workflow)

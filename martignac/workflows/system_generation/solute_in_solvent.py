import logging

from martignac import config
from martignac.nomad.workflows import Job, build_nomad_workflow
from martignac.parsers.gromacs_topologies import Topology, append_all_includes_to_top
from martignac.utils.gromacs import (
    gromacs_simulation_command,
    solvate_solute_command,
)
from martignac.utils.martini_flow_projects import (
    MartiniFlowProject,
    fetched_from_nomad,
    flag_ready_for_upload,
    import_job_from_other_flow,
    is_ready_for_upload,
    store_gromacs_log_to_doc,
    store_task,
    store_workflow,
    uploaded_to_nomad,
)
from martignac.utils.misc import func_name
from martignac.workflows.system_generation.solute import SoluteGenFlow, get_solute_job
from martignac.workflows.system_generation.solvent import SolventGenFlow, get_solvent_job

logger = logging.getLogger(__name__)

conf = config()["solute_in_solvent"]


class SoluteInSolventFlow(MartiniFlowProject):
    workspace_path: str = f"{MartiniFlowProject.workspaces_path}/{conf['relative_paths']['workspaces']}"
    itp_files = {k: v.get(str) for k, v in conf["itp_files"].items()}
    mdp_path = f"{MartiniFlowProject.input_files_path}/{conf['relative_paths']['mdp_files']}"
    mdp_files = {k: v.get(str) for k, v in conf["mdp_files"].items()}
    simulation_settings = {"n_threads": conf["settings"]["n_threads"].get(int)}
    system_name = conf["output_names"]["system"].get(str)
    nomad_workflow: str = conf["output_names"]["nomad_workflow"].get(str)
    nomad_top_level_workflow: str = conf["output_names"]["nomad_top_level_workflow"].get(str)
    state_names = {k: v.get(str) for k, v in conf["output_names"]["states"].items()}


project_name = SoluteInSolventFlow.class_name()
solute_gen_name = SoluteGenFlow.class_name()
solvent_gen_name = SolventGenFlow.class_name()


@SoluteInSolventFlow.label
def solute_generated(job) -> bool:
    return (
        solute_gen_name in job.doc
        and job.doc[solute_gen_name].get("solute_gro")
        and job.doc[solute_gen_name].get("solute_top")
    )


@SoluteInSolventFlow.label
def solvent_generated(job) -> bool:
    return solvent_gen_name in job.doc and job.doc[solvent_gen_name].get("solvent_gro")


@SoluteInSolventFlow.pre(fetched_from_nomad)
@SoluteInSolventFlow.post(solute_generated, tag="solute_generated")
@SoluteInSolventFlow.operation_hooks.on_success(store_workflow)
@SoluteInSolventFlow.operation
def generate_solute(job: Job):
    solute_job: Job = get_solute_job(job)
    keys_for_files_to_copy = ["solute_gro", "solute_top", "solute_itp"]
    project.operation_to_workflow[func_name()] = solute_gen_name
    import_job_from_other_flow(job, solute_gen_project, solute_job, keys_for_files_to_copy)


@SoluteInSolventFlow.pre(fetched_from_nomad)
@SoluteInSolventFlow.post(solvent_generated, tag="solvent_generated")
@SoluteInSolventFlow.operation_hooks.on_success(store_workflow)
@SoluteInSolventFlow.operation
def generate_solvent(job):
    solvent_job: Job = get_solvent_job(job)
    keys_for_files_to_copy = ["solvent_gro", "solvent_top"]
    project.operation_to_workflow[func_name()] = solvent_gen_name
    import_job_from_other_flow(job, solvent_gen_project, solvent_job, keys_for_files_to_copy)


@SoluteInSolventFlow.label
def system_generated(job):
    return job.isfile(SoluteInSolventFlow.get_state_name("generate", "gro"))


@SoluteInSolventFlow.label
def system_minimized(job):
    return job.isfile(SoluteInSolventFlow.get_state_name("minimize", "gro"))


@SoluteInSolventFlow.label
def system_equilibrated(job):
    return job.isfile(SoluteInSolventFlow.get_state_name("equilibrate", "gro"))


@SoluteInSolventFlow.pre(solute_generated, tag="solute_generated")
@SoluteInSolventFlow.pre(solvent_generated, tag="solvent_generated")
@SoluteInSolventFlow.post(system_generated, tag="system_generated")
@SoluteInSolventFlow.operation_hooks.on_success(store_task)
@SoluteInSolventFlow.operation(cmd=True, with_job=True)
def solvate(job):
    return solvate_solute_command(
        gro_solute=job.doc[solute_gen_name]["solute_gro"],
        gro_solvent=job.doc[solvent_gen_name]["solvent_gro"],
        top_solute=job.doc[solute_gen_name]["solute_top"],
        top_output=SoluteInSolventFlow.get_state_name(extension="top"),
        output_name=SoluteInSolventFlow.get_state_name("generate"),
    )


@SoluteInSolventFlow.pre(system_generated, tag="system_generated")
@SoluteInSolventFlow.post(system_minimized, tag="system_minimized")
@SoluteInSolventFlow.operation_hooks.on_success(store_gromacs_log_to_doc)
@SoluteInSolventFlow.operation(cmd=True, with_job=True)
def minimize(job):
    solute_top = Topology.parse_top_file(job.doc[solute_gen_name]["solute_top"])
    solvent_top = Topology.parse_top_file(job.doc[solvent_gen_name]["solvent_top"])
    solute_solvent_top = Topology.parse_top_file(SoluteInSolventFlow.get_state_name(extension="top"))
    new_top = append_all_includes_to_top(solute_solvent_top, [solvent_top, solute_top])
    new_top.output_top(SoluteInSolventFlow.get_state_name(extension="top"))
    job.doc[project_name]["solute_solvent_top"] = SoluteInSolventFlow.get_state_name(extension="top")
    return gromacs_simulation_command(
        mdp=SoluteInSolventFlow.mdp_files.get("minimize"),
        top=SoluteInSolventFlow.get_state_name(extension="top"),
        gro=SoluteInSolventFlow.get_state_name("generate", "gro"),
        name=SoluteInSolventFlow.get_state_name("minimize"),
        n_threads=SoluteInSolventFlow.simulation_settings.get("n_threads"),
    )


@SoluteInSolventFlow.pre(system_minimized, tag="system_minimized")
@SoluteInSolventFlow.post(system_equilibrated, tag="system_equilibrated")
@SoluteInSolventFlow.operation_hooks.on_success(store_gromacs_log_to_doc)
@SoluteInSolventFlow.operation(cmd=True, with_job=True)
def equilibrate(job):
    job.doc[project_name]["solute_solvent_gro"] = SoluteInSolventFlow.get_state_name("equilibrate", "gro")
    return gromacs_simulation_command(
        mdp=SoluteInSolventFlow.mdp_files.get("equilibrate"),
        top=SoluteInSolventFlow.get_state_name(extension="top"),
        gro=SoluteInSolventFlow.get_state_name("minimize", "gro"),
        name=SoluteInSolventFlow.get_state_name("equilibrate"),
        n_threads=SoluteInSolventFlow.simulation_settings.get("n_threads"),
    )


@SoluteInSolventFlow.pre(system_equilibrated, tag="system_equilibrated")
@SoluteInSolventFlow.post(
    lambda job: all(
        [job.isfile(SoluteInSolventFlow.nomad_workflow), job.isfile(SoluteInSolventFlow.nomad_top_level_workflow)]
    ),
    tag="generated_nomad_workflow",
)
@SoluteInSolventFlow.operation_hooks.on_success(flag_ready_for_upload)
@SoluteInSolventFlow.operation(with_job=True)
def generate_nomad_workflow(job):
    build_nomad_workflow(job, is_top_level=False)
    build_nomad_workflow(job, is_top_level=True)


@SoluteInSolventFlow.pre(is_ready_for_upload, tag="generated_nomad_workflow")
@SoluteInSolventFlow.post(uploaded_to_nomad)
@SoluteInSolventFlow.operation(with_job=True)
def upload_to_nomad(job: Job):
    return SoluteInSolventFlow.upload_to_nomad(job)


def get_solvation_job(job: Job) -> Job:
    sp = {
        "type": "solute_solvation",
        "solvent_name": job.sp.get("solvent_name"),
        "solute_name": job.sp.get("solute_name"),
    }
    return project.open_job(sp).init()


project = SoluteInSolventFlow.get_project(path=SoluteInSolventFlow.workspace_path)
solute_gen_project = SoluteGenFlow.get_project(path=SoluteGenFlow.workspace_path)
solvent_gen_project = SolventGenFlow.get_project(path=SolventGenFlow.workspace_path)
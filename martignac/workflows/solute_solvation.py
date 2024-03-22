import logging
from typing import cast

import networkx as nx

from martignac import config
from martignac.nomad.workflows import Job, NomadTopLevelWorkflow, NomadWorkflow
from martignac.parsers.gromacs_topologies import Topology, append_all_includes_to_top
from martignac.utils.gromacs import (
    gromacs_simulation_command,
    solvate_solute_command,
)
from martignac.utils.martini_flow_projects import (
    MartiniFlowProject,
    fetched_from_nomad,
    import_job_from_other_flow,
    uploaded_to_nomad,
)
from martignac.workflows.solute_generation import SoluteGenFlow, get_solute_job
from martignac.workflows.solvent_generation import SolventGenFlow, get_solvent_job

logger = logging.getLogger(__name__)

conf = config()["solute_solvation"]


class SoluteSolvationFlow(MartiniFlowProject):
    workspace_path: str = f"{MartiniFlowProject.workspaces_path}/{conf['relative_paths']['workspaces']}"
    itp_path = f"{MartiniFlowProject.input_files_path}/{conf['relative_paths']['itp_files']}"
    itp_files = {k: v.get(str) for k, v in conf["itp_files"].items()}
    mdp_path = f"{MartiniFlowProject.input_files_path}/{conf['relative_paths']['mdp_files']}"
    mdp_files = {k: v.get(str) for k, v in conf["mdp_files"].items()}
    simulation_settings = {"n_threads": conf["settings"]["n_threads"].get(int)}
    system_name = conf["output_names"]["system"].get(str)
    nomad_workflow: str = conf["output_names"]["nomad_workflow"].get(str)
    nomad_top_level_workflow: str = conf["output_names"]["nomad_top_level_workflow"].get(str)
    state_names = {k: v.get(str) for k, v in conf["output_names"]["states"].items()}


@SoluteSolvationFlow.label
def solute_generated(job) -> bool:
    return job.document.get("solute_gro") and job.document.get("solute_top")


@SoluteSolvationFlow.label
def solvent_generated(job) -> bool:
    return job.document.get("solvent_gro")


@SoluteSolvationFlow.pre(fetched_from_nomad)
@SoluteSolvationFlow.post(solute_generated)
@SoluteSolvationFlow.operation
def generate_solute(job: Job):
    solute_job: Job = get_solute_job(job)
    keys_for_files_to_copy = ["solute_gro", "solute_top", "solute_itp"]
    update_keys = ["solute_gro", "solute_top", "solute_itp", "solute_name"]
    import_job_from_other_flow(job, solute_gen_project, solute_job, keys_for_files_to_copy, update_keys)


@SoluteSolvationFlow.pre(fetched_from_nomad)
@SoluteSolvationFlow.post(solvent_generated)
@SoluteSolvationFlow.operation
def generate_solvent(job):
    solvent_job: Job = get_solvent_job(job)
    keys_for_files_to_copy = ["solvent_gro", "solvent_top"]
    update_keys = ["solvent_gro", "solvent_top", "solvent_name"]
    import_job_from_other_flow(job, solvent_gen_project, solvent_job, keys_for_files_to_copy, update_keys)


@SoluteSolvationFlow.label
def system_generated(job):
    return job.isfile(SoluteSolvationFlow.get_state_name("generate", "gro"))


@SoluteSolvationFlow.label
def system_minimized(job):
    return job.isfile(SoluteSolvationFlow.get_state_name("minimize", "gro"))


@SoluteSolvationFlow.label
def system_equilibrated(job):
    return job.isfile(SoluteSolvationFlow.get_state_name("equilibrate", "gro"))


@SoluteSolvationFlow.pre(solute_generated)
@SoluteSolvationFlow.pre(solvent_generated)
@SoluteSolvationFlow.post(system_generated)
@SoluteSolvationFlow.operation(cmd=True, with_job=True)
def solvate(job):
    return solvate_solute_command(
        gro_solute=job.document["solute_gro"],
        gro_solvent=job.document["solvent_gro"],
        top_solute=job.document["solute_top"],
        top_output=SoluteSolvationFlow.get_state_name(extension="top"),
        output_name=SoluteSolvationFlow.get_state_name("generate"),
    )


@SoluteSolvationFlow.pre(system_generated)
@SoluteSolvationFlow.post(system_minimized)
@SoluteSolvationFlow.operation(cmd=True, with_job=True)
@SoluteSolvationFlow.log_gromacs_simulation(SoluteSolvationFlow.get_state_name("minimize"))
def minimize(job):
    solute_top = Topology.parse_top_file(job.document["solute_top"])
    solvent_top = Topology.parse_top_file(job.document["solvent_top"])
    solute_solvent_top = Topology.parse_top_file(SoluteSolvationFlow.get_state_name(extension="top"))
    new_top = append_all_includes_to_top(solute_solvent_top, [solvent_top, solute_top])
    new_top.output_top(SoluteSolvationFlow.get_state_name(extension="top"))
    job.document["solute_solvent_top"] = SoluteSolvationFlow.get_state_name(extension="top")
    return gromacs_simulation_command(
        mdp=SoluteSolvationFlow.mdp_files.get("minimize"),
        top=SoluteSolvationFlow.get_state_name(extension="top"),
        gro=SoluteSolvationFlow.get_state_name("generate", "gro"),
        name=SoluteSolvationFlow.get_state_name("minimize"),
        n_threads=SoluteSolvationFlow.simulation_settings.get("n_threads"),
    )


@SoluteSolvationFlow.pre(system_minimized)
@SoluteSolvationFlow.post(system_equilibrated)
@SoluteSolvationFlow.operation(cmd=True, with_job=True)
@SoluteSolvationFlow.log_gromacs_simulation(SoluteSolvationFlow.get_state_name("equilibrate"))
def equilibrate(job):
    job.document["solute_solvent_gro"] = SoluteSolvationFlow.get_state_name("equilibrate", "gro")
    return gromacs_simulation_command(
        mdp=SoluteSolvationFlow.mdp_files.get("equilibrate"),
        top=SoluteSolvationFlow.get_state_name(extension="top"),
        gro=SoluteSolvationFlow.get_state_name("minimize", "gro"),
        name=SoluteSolvationFlow.get_state_name("equilibrate"),
        n_threads=SoluteSolvationFlow.simulation_settings.get("n_threads"),
    )


@SoluteSolvationFlow.pre(system_equilibrated)
@SoluteSolvationFlow.post(
    lambda job: all(
        [job.isfile(SoluteSolvationFlow.nomad_workflow), job.isfile(SoluteSolvationFlow.nomad_top_level_workflow)]
    )
)
@SoluteSolvationFlow.operation(with_job=True)
def generate_nomad_workflow(job):
    projects = {
        "SoluteSolvationFlow": project,
        "SoluteGenFlow": solute_gen_project,
        "SolventGenFlow": solvent_gen_project,
    }
    jobs = {
        "SoluteSolvationFlow": job,
        "SoluteGenFlow": get_solute_job(job),
        "SolventGenFlow": get_solvent_job(job),
    }
    for workflow_name in projects:
        workflow = NomadWorkflow(projects[workflow_name], job)
        workflow_flow = cast(MartiniFlowProject, globals()[workflow_name])
        workflow.build_workflow_yaml(workflow_flow.nomad_workflow)
    graph = nx.DiGraph([("SoluteGenFlow", "SoluteSolvationFlow"), ("SolventGenFlow", "SoluteSolvationFlow")])
    top_level_workflow = NomadTopLevelWorkflow(projects, jobs, graph)
    top_level_workflow.build_workflow_yaml(SoluteSolvationFlow.nomad_top_level_workflow)


@SoluteSolvationFlow.pre(lambda job: job.isfile(SoluteSolvationFlow.nomad_top_level_workflow))
@SoluteSolvationFlow.post(lambda job: uploaded_to_nomad(job))
@SoluteSolvationFlow.operation(with_job=True)
def upload_to_nomad(job: Job):
    return SoluteSolvationFlow.upload_to_nomad(job)


def get_solvation_job(job: Job) -> Job:
    sp = {
        "type": "solute_solvation",
        "solvent_name": job.sp.get("solvent_name"),
        "solute_name": job.sp.get("solute_name"),
    }
    return project.open_job(sp).init()


project = SoluteSolvationFlow.get_project(path=SoluteSolvationFlow.workspace_path)
solute_gen_project = SoluteGenFlow.get_project(path=SoluteGenFlow.workspace_path)
solvent_gen_project = SolventGenFlow.get_project(path=SolventGenFlow.workspace_path)

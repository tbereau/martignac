import logging

import networkx as nx

from martignac import config
from martignac.nomad.workflows import NomadTopLevelWorkflow, NomadWorkflow
from martignac.parsers.gromacs_topologies import Topology, append_all_includes_to_top
from martignac.utils.gromacs import (
    gromacs_simulation_command,
    solvate_solute_command,
)
from martignac.utils.martini_flow_projects import Job, MartiniFlowProject
from martignac.workflows.solute_generation import SoluteGenFlow
from martignac.workflows.solvent_generation import SolventGenFlow

logger = logging.getLogger(__name__)

conf = config()["solute_solvation"]


class SoluteSolvationFlow(MartiniFlowProject):
    minimize_mdp: str = conf["mdp_files"]["minimize"].get(str)
    equilibrate_mdp: str = conf["mdp_files"]["equilibrate"].get(str)
    n_threads: int = conf["settings"]["n_threads"].get(int)
    system_name: str = conf["output_names"]["system"].get(str)
    nomad_workflow: str = conf["output_names"]["nomad_workflow"].get(str)
    nomad_top_level_workflow: str = conf["output_names"]["nomad_top_level_workflow"].get(str)
    generate_name: str = conf["output_names"]["states"]["generate"].get(str)
    minimize_name: str = conf["output_names"]["states"]["minimize"].get(str)
    equilibrate_name: str = conf["output_names"]["states"]["equilibrate"].get(str)

    @classmethod
    def get_system_gen_name(cls) -> str:
        return f"{cls.system_name}_{cls.generate_name}"

    @classmethod
    def get_system_gen_gro(cls) -> str:
        return f"{cls.get_system_gen_name()}.gro"

    @classmethod
    def get_system_min_name(cls) -> str:
        return f"{cls.system_name}_{cls.minimize_name}"

    @classmethod
    def get_system_min_gro(cls) -> str:
        return f"{cls.get_system_min_name()}.gro"

    @classmethod
    def get_system_equ_name(cls) -> str:
        return f"{cls.system_name}_{cls.equilibrate_name}"

    @classmethod
    def get_system_equ_gro(cls) -> str:
        return f"{cls.get_system_equ_name()}.gro"

    @classmethod
    def get_top(cls) -> str:
        return f"{cls.system_name}.top"


project = SoluteSolvationFlow().get_project()


@SoluteSolvationFlow.label
def solute_generated(job) -> bool:
    return job.document.get("solute_gro") and job.document.get("solute_top")


@SoluteSolvationFlow.label
def solvent_generated(job) -> bool:
    return job.document.get("solvent_gro")


@SoluteSolvationFlow.post(solute_generated)
@SoluteSolvationFlow.operation
def generate_solute(job):
    job.document["upload_to_nomad"] = False
    SoluteGenFlow().run(jobs=[job])
    job.document["upload_to_nomad"] = True


@SoluteSolvationFlow.post(solvent_generated)
@SoluteSolvationFlow.operation
def generate_solvent(job):
    job.document["upload_to_nomad"] = False
    SolventGenFlow().run(jobs=[job])
    job.document["upload_to_nomad"] = True


@SoluteSolvationFlow.label
def system_generated(job):
    return job.isfile(SoluteSolvationFlow.get_system_gen_gro())


@SoluteSolvationFlow.label
def system_minimized(job):
    return job.isfile(SoluteSolvationFlow.get_system_min_gro())


@SoluteSolvationFlow.label
def system_equilibrated(job):
    return job.isfile(SoluteSolvationFlow.get_system_equ_gro())


@SoluteSolvationFlow.pre(solute_generated)
@SoluteSolvationFlow.pre(solvent_generated)
@SoluteSolvationFlow.post(system_generated)
@SoluteSolvationFlow.operation(cmd=True, with_job=True)
def solvate(job):
    return solvate_solute_command(
        gro_solute=job.document["solute_gro"],
        gro_solvent=job.document["solvent_gro"],
        top_solute=job.document["solute_top"],
        top_output=SoluteSolvationFlow.get_top(),
        output_name=SoluteSolvationFlow.get_system_gen_name(),
    )


@SoluteSolvationFlow.pre(system_generated)
@SoluteSolvationFlow.post(system_minimized)
@SoluteSolvationFlow.operation(cmd=True, with_job=True)
@SoluteSolvationFlow.log_gromacs_simulation(SoluteSolvationFlow.get_system_min_name())
def minimize(job):
    solute_top = Topology.parse_top_file(job.document["solute_top"])
    solvent_top = Topology.parse_top_file(job.document["solvent_top"])
    solute_solvent_top = Topology.parse_top_file(SoluteSolvationFlow.get_top())
    new_top = append_all_includes_to_top(solute_solvent_top, [solvent_top, solute_top])
    new_top.output_top(SoluteSolvationFlow.get_top())
    job.document["solute_solvent_top"] = SoluteSolvationFlow.get_top()
    return gromacs_simulation_command(
        mdp=SoluteSolvationFlow.minimize_mdp,
        top=SoluteSolvationFlow.get_top(),
        gro=SoluteSolvationFlow.get_system_gen_gro(),
        name=SoluteSolvationFlow.get_system_min_name(),
        n_threads=SoluteSolvationFlow.n_threads,
    )


@SoluteSolvationFlow.pre(system_minimized)
@SoluteSolvationFlow.post(system_equilibrated)
@SoluteSolvationFlow.operation(cmd=True, with_job=True)
@SoluteSolvationFlow.log_gromacs_simulation(SoluteSolvationFlow.get_system_equ_name())
def equilibrate(job):
    job.document["solute_solvent_gro"] = SoluteSolvationFlow.get_system_equ_gro()
    return gromacs_simulation_command(
        mdp=SoluteSolvationFlow.equilibrate_mdp,
        top=SoluteSolvationFlow.get_top(),
        gro=SoluteSolvationFlow.get_system_min_gro(),
        name=SoluteSolvationFlow.get_system_equ_name(),
        n_threads=SoluteSolvationFlow.n_threads,
    )


@SoluteSolvationFlow.pre(system_equilibrated)
@SoluteSolvationFlow.post(
    lambda job: all(
        [job.isfile(SoluteSolvationFlow.nomad_workflow), job.isfile(SoluteSolvationFlow.nomad_top_level_workflow)]
    )
)
@SoluteSolvationFlow.operation(with_job=True)
def generate_nomad_workflow(job):
    workflow = NomadWorkflow(project, job)
    workflow.build_workflow_yaml(SoluteSolvationFlow.nomad_workflow)
    graph = nx.DiGraph([("SoluteGenFlow", "SoluteSolvationFlow"), ("SolventGenFlow", "SoluteSolvationFlow")])
    projects = {
        "SoluteSolvationFlow": project,
        "SoluteGenFlow": SoluteGenFlow.get_project(),
        "SolventGenFlow": SolventGenFlow.get_project(),
    }
    top_level_workflow = NomadTopLevelWorkflow(projects, job, graph)
    top_level_workflow.build_workflow_yaml(SoluteSolvationFlow.nomad_top_level_workflow)


@SoluteSolvationFlow.pre(lambda job: job.isfile(SoluteSolvationFlow.nomad_top_level_workflow))
@SoluteSolvationFlow.post(lambda job: job.document.get("nomad_upload_id", False))
@SoluteSolvationFlow.operation(with_job=True)
def upload_to_nomad(job: Job):
    return SoluteSolvationFlow.upload_to_nomad(job)

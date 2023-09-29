from flow import FlowProject

from martignac.utils.gromacs import (
    solvate_command,
    gromacs_simulation_command,
    copy_files_to,
)


ITP_FILES = [
    "../../resource_files/martini.itp",
    "../../resource_files/martini_na.itp",
    "../../resource_files/martini_solvents.itp"
]
SOLUTE_GRO = "../../resource_files/Na.gro"
SOLVENT_GRO = "../../resource_files/solvent.gro"
SOLUTE_TOP = "../../resource_files/na_topol.top"
MINIMIZE_MDP = "../../mdp_files/min.mdp"
EQUILIBRATE_MDP = "../../mdp_files/run.mdp"
N_THREADS = 1
TOPOLOGY_OUTPUT = "topology.top"
MINIMIZE_NAME = "minimize"
SOLVATE_NAME = "solvate"
EQUILIBRATE_NAME = "equilibrate"


class SolvationFlow(FlowProject):
    itp_files: list[str] = ITP_FILES
    solute_gro: str = SOLUTE_GRO
    solvent_gro: str = SOLVENT_GRO
    solute_top: str = SOLUTE_TOP
    minimize_mdp: str = MINIMIZE_MDP
    equilibrate_mdp: str = EQUILIBRATE_MDP
    n_threads: int = N_THREADS
    topology_output: str = TOPOLOGY_OUTPUT
    minimize_name: str = MINIMIZE_NAME
    solvate_name: str = SOLVATE_NAME
    equilibrate_name: str = EQUILIBRATE_NAME

    @classmethod
    def get_minimize_gro(cls) -> str:
        return f"{cls.minimize_name}.gro"

    @classmethod
    def get_solvate_gro(cls) -> str:
        return f"{cls.solvate_name}.gro"

    @classmethod
    def get_equilibrate_gro(cls) -> str:
        return f"{cls.equilibrate_name}.gro"


project = SolvationFlow().get_project()


@SolvationFlow.label
def solvated(job):
    return job.isfile(SolvationFlow.get_solvate_gro())


@SolvationFlow.label
def minimized(job):
    return job.isfile(SolvationFlow.get_minimize_gro())


@SolvationFlow.label
def equilibrated(job):
    return job.isfile(SolvationFlow.get_equilibrate_gro())


@SolvationFlow.post(solvated)
@SolvationFlow.operation(cmd=True, with_job=True)
def solvate(job):
    return solvate_command(
        gro_solute=SolvationFlow.solute_gro,
        gro_solvent=SolvationFlow.solvent_gro,
        top_solute=SolvationFlow.solute_top,
        top_output=SolvationFlow.topology_output,
        output_name=SolvationFlow.solvate_name
    )


@SolvationFlow.pre(solvated)
@SolvationFlow.post(minimized)
@SolvationFlow.operation(cmd=True, with_job=True)
def minimize(job):
    copy_files_to(SolvationFlow.itp_files, ".")
    return gromacs_simulation_command(
        mdp=SolvationFlow.minimize_mdp,
        top=SolvationFlow.topology_output,
        gro=SolvationFlow.get_solvate_gro(),
        name=SolvationFlow.minimize_name
    )


@SolvationFlow.pre(minimized)
@SolvationFlow.post(equilibrated)
@SolvationFlow.operation(cmd=True, with_job=True)
def equilibrate(job):
    copy_files_to(SolvationFlow.itp_files, ".")
    return gromacs_simulation_command(
        mdp=SolvationFlow.equilibrate_mdp,
        top=SolvationFlow.topology_output,
        gro=SolvationFlow.get_minimize_gro(),
        name=SolvationFlow.equilibrate_name
    )

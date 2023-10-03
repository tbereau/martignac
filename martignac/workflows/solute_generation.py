from flow import FlowProject

from martignac.utils.gromacs import (
    gromacs_simulation_command
)
from martignac.parsers.gromacs_particle_definitions import (
    find_particle_type_from_name,
    generate_gro_for_particle_type,
    generate_top_file_for_particle,
    generate_itp_file_for_particle,
)
from martignac import config

conf = config()["solute_generation"]


class SoluteGenFlow(FlowProject):
    particle_def_itp: str = conf["itp_files"]["particle_definition"].get(str)
    minimize_mdp: str = conf["mdp_files"]["minimize"].get(str)
    n_threads: int = conf["settings"]["n_threads"].get(int)
    system_name: str = conf["output_names"]["system"].get(str)
    generate_name: str = conf["output_names"]["states"]["generate"].get(str)
    minimize_name: str = conf["output_names"]["states"]["minimize"].get(str)

    @classmethod
    def get_solute_mol_name(cls) -> str:
        return f"{cls.system_name}_{cls.generate_name}"

    @classmethod
    def get_solute_mol_gro(cls) -> str:
        return f"{cls.get_solute_mol_name()}.gro"

    @classmethod
    def get_solute_min_name(cls) -> str:
        return f"{cls.system_name}_{cls.minimize_name}"

    @classmethod
    def get_solute_min_gro(cls) -> str:
        return f"{cls.get_solute_min_name()}.gro"

    @classmethod
    def get_solute_itp(cls) -> str:
        return f"{cls.system_name}.itp"

    @classmethod
    def get_solute_top(cls) -> str:
        return f"{cls.system_name}.top"


project = SoluteGenFlow().get_project()


@SoluteGenFlow.label
def generated(job):
    return job.isfile(SoluteGenFlow.get_solute_mol_gro())


@SoluteGenFlow.label
def minimized(job):
    return job.isfile(SoluteGenFlow.get_solute_min_gro())


@SoluteGenFlow.post(generated)
@SoluteGenFlow.operation(with_job=True)
def solvate(job) -> None:
    molecule = find_particle_type_from_name(
        SoluteGenFlow.particle_def_itp, job.sp.solute_name
    )
    generate_gro_for_particle_type(molecule, SoluteGenFlow.get_solute_mol_gro())
    generate_itp_file_for_particle(molecule, SoluteGenFlow.get_solute_itp())
    generate_top_file_for_particle(
        molecule,
        [SoluteGenFlow.particle_def_itp, SoluteGenFlow.get_solute_itp()],
        SoluteGenFlow.get_solute_top()
    )
    job.document["solute_itp"] = SoluteGenFlow.get_solute_itp()
    job.document["solute_top"] = SoluteGenFlow.get_solute_top()
    return None


@SoluteGenFlow.pre(generated)
@SoluteGenFlow.post(minimized)
@SoluteGenFlow.operation(cmd=True, with_job=True)
def minimize(job):
    job.document["solute_gro"] = SoluteGenFlow.get_solute_min_gro()
    return gromacs_simulation_command(
        mdp=SoluteGenFlow.minimize_mdp,
        top=SoluteGenFlow.get_solute_top(),
        gro=SoluteGenFlow.get_solute_mol_gro(),
        name=SoluteGenFlow.get_solute_min_name(),
        n_threads=SoluteGenFlow.n_threads
    )

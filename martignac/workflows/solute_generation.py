from flow import FlowProject

from martignac import config
from martignac.parsers.gromacs_forcefields import (
    generate_gro_file_for_molecule,
    generate_itp_file_for_molecule,
    generate_top_file_for_molecule,
    get_molecule_from_name,
)
from martignac.utils.gromacs import gromacs_simulation_command

conf = config()["solute_generation"]


class SoluteGenFlow(FlowProject):
    particle_def_itp: str = conf["itp_files"]["particle_definition"].get(str)
    minimize_mdp: str = conf["mdp_files"]["minimize"].get(str)
    n_threads: int = conf["settings"]["n_threads"].get(int)
    system_name: str = conf["output_names"]["system"].get(str)
    generate_name: str = conf["output_names"]["states"]["generate"].get(str)
    minimize_name: str = conf["output_names"]["states"]["minimize"].get(str)
    bond_length: float = conf["parameters"]["bond_length"].get(float)
    number_excl: int = conf["parameters"]["number_excl"].get(int)
    bond_constant: float = conf["parameters"]["bond_constant"].get(float)

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
    return (
        job.isfile(SoluteGenFlow.get_solute_mol_gro())
        and job.isfile(SoluteGenFlow.get_solute_itp())
        and job.isfile(SoluteGenFlow.get_solute_top())
    )


@SoluteGenFlow.label
def minimized(job):
    return job.isfile(SoluteGenFlow.get_solute_min_gro())


@SoluteGenFlow.post(generated)
@SoluteGenFlow.operation(with_job=True)
def solvate(job) -> None:
    molecule = get_molecule_from_name(
        job.sp.solute_name,
        bond_length=SoluteGenFlow.bond_length,
        bond_constant=SoluteGenFlow.bond_constant,
        number_excl=SoluteGenFlow.number_excl,
    )
    generate_gro_file_for_molecule(molecule, SoluteGenFlow.get_solute_mol_gro())
    generate_itp_file_for_molecule(molecule, SoluteGenFlow.get_solute_itp())
    generate_top_file_for_molecule(
        molecule,
        [SoluteGenFlow.particle_def_itp, SoluteGenFlow.get_solute_itp()],
        SoluteGenFlow.get_solute_top(),
    )
    job.document["solute_itp"] = SoluteGenFlow.get_solute_itp()
    job.document["solute_top"] = SoluteGenFlow.get_solute_top()
    job.document["solute_name"] = molecule.name
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
        n_threads=SoluteGenFlow.n_threads,
    )

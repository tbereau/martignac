from flow import FlowProject

from martignac.utils.gromacs import (
    generate_solvent_command,
    gromacs_simulation_command,
)
from martignac.parsers.gromacs_forcefields import (
    find_molecule_from_name,
    generate_gro_file_for_molecule,
    generate_top_file_for_molecule,
)
from martignac.parsers.gromacs_structures import (
    get_number_of_molecules_from_gro
)
from martignac import config

conf = config()["solvent_generation"]


class SolventGenFlow(FlowProject):
    itp_files: list[str] = [str(e) for e in conf["itp_files"].get()]
    n_threads: int = conf["settings"]["n_threads"].get(int)
    box_length: float = conf["settings"]["box_length"].get(float)
    minimize_mdp: str = conf["mdp_files"]["minimize"].get(str)
    equilibrate_mdp: str = conf["mdp_files"]["equilibrate"].get(str)
    production_mdp: str = conf["mdp_files"]["production"].get(str)
    system_name: str = conf["output_names"]["system"].get(str)
    gen_mol_name: str = conf["output_names"]["states"]["mol"].get(str)
    gen_box_name: str = conf["output_names"]["states"]["box"].get(str)
    min_name: str = conf["output_names"]["states"]["min"].get(str)
    equ_name: str = conf["output_names"]["states"]["equ"].get(str)
    prod_name: str = conf["output_names"]["states"]["prod"].get(str)

    @classmethod
    def get_solvent_mol_name(cls) -> str:
        return f"{cls.system_name}_{cls.gen_mol_name}"

    @classmethod
    def get_solvent_box_name(cls) -> str:
        return f"{cls.system_name}_{cls.gen_box_name}"

    @classmethod
    def get_solvent_min_name(cls) -> str:
        return f"{cls.system_name}_{cls.min_name}"

    @classmethod
    def get_solvent_equ_name(cls) -> str:
        return f"{cls.system_name}_{cls.equ_name}"

    @classmethod
    def get_solvent_prod_name(cls) -> str:
        return f"{cls.system_name}_{cls.prod_name}"

    @classmethod
    def get_solvent_mol_gro(cls) -> str:
        return f"{cls.get_solvent_mol_name()}.gro"

    @classmethod
    def get_solvent_box_gro(cls) -> str:
        return f"{cls.get_solvent_box_name()}.gro"

    @classmethod
    def get_solvent_min_gro(cls) -> str:
        return f"{cls.get_solvent_min_name()}.gro"

    @classmethod
    def get_solvent_equ_gro(cls) -> str:
        return f"{cls.get_solvent_equ_name()}.gro"

    @classmethod
    def get_solvent_prod_gro(cls) -> str:
        return f"{cls.get_solvent_prod_name()}.gro"

    @classmethod
    def get_solvent_mol_top(cls) -> str:
        return f"{cls.get_solvent_mol_name()}.top"

    @classmethod
    def get_solvent_box_top(cls) -> str:
        return f"{cls.get_solvent_box_name()}.top"


project = SolventGenFlow().get_project()


@SolventGenFlow.label
def generated_mol_gro(job):
    return job.isfile(SolventGenFlow.get_solvent_mol_gro())


@SolventGenFlow.label
def generated_box_gro(job):
    return job.isfile(SolventGenFlow.get_solvent_box_gro())


@SolventGenFlow.label
def generated_min_gro(job):
    return job.isfile(SolventGenFlow.get_solvent_min_gro())


@SolventGenFlow.label
def generated_equ_gro(job):
    return job.isfile(SolventGenFlow.get_solvent_equ_gro())


@SolventGenFlow.label
def generated_prod_gro(job):
    return job.isfile(SolventGenFlow.get_solvent_prod_gro())


@SolventGenFlow.post(generated_mol_gro)
@SolventGenFlow.operation(with_job=True)
def generate_solvent_molecule(job) -> None:
    molecule = find_molecule_from_name(
        SolventGenFlow.itp_files, job.sp.solvent_name
    )
    generate_gro_file_for_molecule(molecule, SolventGenFlow.get_solvent_mol_gro())
    generate_top_file_for_molecule(
        molecule, SolventGenFlow.itp_files, SolventGenFlow.get_solvent_mol_top()
    )
    return None


@SolventGenFlow.pre(generated_mol_gro)
@SolventGenFlow.post(generated_box_gro)
@SolventGenFlow.operation(cmd=True, with_job=True)
def solvate(job):
    return generate_solvent_command(
        gro_solvent_mol=SolventGenFlow.get_solvent_mol_gro(),
        box_length=SolventGenFlow.box_length,
        output_name=SolventGenFlow.get_solvent_box_name(),
    )


@SolventGenFlow.pre(generated_box_gro)
@SolventGenFlow.post(generated_min_gro)
@SolventGenFlow.operation(cmd=True, with_job=True)
def minimize(job):
    molecule = find_molecule_from_name(
        SolventGenFlow.itp_files, job.sp.solvent_name
    )
    generate_top_file_for_molecule(
        molecule,
        SolventGenFlow.itp_files,
        SolventGenFlow.get_solvent_box_top(),
        num_molecules=get_number_of_molecules_from_gro(
            SolventGenFlow.get_solvent_box_gro()
        )
    )
    job.document["solvent_top"] = SolventGenFlow.get_solvent_box_top()
    return gromacs_simulation_command(
        mdp=SolventGenFlow.minimize_mdp,
        top=SolventGenFlow.get_solvent_box_top(),
        gro=SolventGenFlow.get_solvent_box_gro(),
        name=SolventGenFlow.get_solvent_min_name(),
        n_threads=SolventGenFlow.n_threads,
    )


@SolventGenFlow.pre(generated_min_gro)
@SolventGenFlow.post(generated_equ_gro)
@SolventGenFlow.operation(cmd=True, with_job=True)
def equilibrate(job):
    return gromacs_simulation_command(
        mdp=SolventGenFlow.equilibrate_mdp,
        top=SolventGenFlow.get_solvent_box_top(),
        gro=SolventGenFlow.get_solvent_min_gro(),
        name=SolventGenFlow.get_solvent_equ_name(),
        n_threads=SolventGenFlow.n_threads
    )


@SolventGenFlow.pre(generated_equ_gro)
@SolventGenFlow.post(generated_prod_gro)
@SolventGenFlow.operation(cmd=True, with_job=True)
def run_production(job):
    job.document["solvent_gro"] = SolventGenFlow.get_solvent_prod_gro()
    return gromacs_simulation_command(
        mdp=SolventGenFlow.production_mdp,
        top=SolventGenFlow.get_solvent_box_top(),
        gro=SolventGenFlow.get_solvent_equ_gro(),
        name=SolventGenFlow.get_solvent_prod_name(),
        n_threads=SolventGenFlow.n_threads
    )

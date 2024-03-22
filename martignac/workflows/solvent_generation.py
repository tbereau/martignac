from martignac import config
from martignac.nomad.workflows import Job, NomadWorkflow
from martignac.parsers.gromacs_forcefields import (
    find_molecule_from_name,
    generate_gro_file_for_molecule,
    generate_top_file_for_molecule,
)
from martignac.parsers.gromacs_structures import get_number_of_molecules_from_gro
from martignac.utils.gromacs import (
    gromacs_simulation_command,
)
from martignac.utils.martini_flow_projects import MartiniFlowProject, fetched_from_nomad, uploaded_to_nomad
from martignac.utils.misc import convert_pdb_to_gro
from martignac.utils.packmol import generate_solvent_with_packmol

conf = config()["solvent_generation"]


class SolventGenFlow(MartiniFlowProject):
    workspace_path: str = f"{MartiniFlowProject.workspaces_path}/{conf['relative_paths']['workspaces']}"
    itp_path = f"{MartiniFlowProject.input_files_path}/{conf['relative_paths']['itp_files']}"
    mdp_path = f"{MartiniFlowProject.input_files_path}/{conf['relative_paths']['mdp_files']}"
    itp_files = {k: v.get(str) for k, v in conf["itp_files"].items()}
    mdp_files = {k: v.get(str) for k, v in conf["mdp_files"].items()}
    simulation_settings = {
        "n_threads": conf["settings"]["n_threads"].get(int),
        "box_length": conf["settings"]["box_length"].get(float),
    }
    system_name = conf["output_names"]["system"].get(str)
    nomad_workflow: str = conf["output_names"]["nomad_workflow"].get(str)
    state_names = {k: v.get(str) for k, v in conf["output_names"]["states"].items()}


@SolventGenFlow.label
def generated_box_pdb(job):
    return job.isfile(SolventGenFlow.get_state_name("box", "pdb"))


@SolventGenFlow.label
def generated_mol_gro(job):
    return job.isfile(SolventGenFlow.get_state_name("mol", "gro"))


@SolventGenFlow.label
def generated_box_gro(job):
    return job.isfile(SolventGenFlow.get_state_name("box", "gro"))


@SolventGenFlow.label
def generated_min_gro(job):
    return job.isfile(SolventGenFlow.get_state_name("min", "gro"))


@SolventGenFlow.label
def generated_equ_gro(job):
    return job.isfile(SolventGenFlow.get_state_name("equ", "gro"))


@SolventGenFlow.label
def generated_prod_gro(job):
    return job.isfile(SolventGenFlow.get_state_name("prod", "gro"))


@SolventGenFlow.pre(fetched_from_nomad)
@SolventGenFlow.post(generated_mol_gro)
@SolventGenFlow.operation(with_job=True)
def generate_solvent_molecule(job) -> None:
    molecule = find_molecule_from_name(list(SolventGenFlow.itp_files.values()), job.sp.solvent_name)
    generate_gro_file_for_molecule(molecule, SolventGenFlow.get_state_name("mol", "gro"))
    generate_top_file_for_molecule(
        molecule, list(SolventGenFlow.itp_files.values()), SolventGenFlow.get_state_name("mol", "top")
    )
    return None


@SolventGenFlow.pre(generated_mol_gro)
@SolventGenFlow.post(generated_box_pdb)
@SolventGenFlow.operation(cmd=True, with_job=True)
def solvate(_):
    return generate_solvent_with_packmol(
        gro_solvent_mol=SolventGenFlow.get_state_name("mol", "gro"),
        box_length=SolventGenFlow.simulation_settings.get("box_length"),
        output_pdb=SolventGenFlow.get_state_name("box", "pdb"),
    )


@SolventGenFlow.pre(generated_box_pdb)
@SolventGenFlow.post(generated_box_gro)
@SolventGenFlow.operation(with_job=True)
def solvate_convert_to_gro(_):
    convert_pdb_to_gro(
        SolventGenFlow.get_state_name("box", "pdb"),
        SolventGenFlow.get_state_name("box", "gro"),
        SolventGenFlow.simulation_settings.get("box_length"),
    )


@SolventGenFlow.pre(generated_box_gro)
@SolventGenFlow.post(generated_min_gro)
@SolventGenFlow.operation(cmd=True, with_job=True)
@SolventGenFlow.log_gromacs_simulation(SolventGenFlow.get_state_name("min"))
def minimize(job):
    molecule = find_molecule_from_name(list(SolventGenFlow.itp_files.values()), job.sp.solvent_name)
    generate_top_file_for_molecule(
        molecule,
        list(SolventGenFlow.itp_files.values()),
        SolventGenFlow.get_state_name("box", "top"),
        num_molecules=get_number_of_molecules_from_gro(SolventGenFlow.get_state_name("box", "gro")),
    )
    job.document["solvent_top"] = SolventGenFlow.get_state_name("box", "top")
    job.document["solvent_name"] = job.sp.solvent_name
    return gromacs_simulation_command(
        mdp=SolventGenFlow.mdp_files.get("minimize"),
        top=SolventGenFlow.get_state_name("box", "top"),
        gro=SolventGenFlow.get_state_name("box", "gro"),
        name=SolventGenFlow.get_state_name("min"),
        n_threads=SolventGenFlow.simulation_settings.get("n_threads"),
    )


@SolventGenFlow.pre(generated_min_gro)
@SolventGenFlow.post(generated_equ_gro)
@SolventGenFlow.operation(cmd=True, with_job=True)
@SolventGenFlow.log_gromacs_simulation(SolventGenFlow.get_state_name("equ"))
def equilibrate(_):
    return gromacs_simulation_command(
        mdp=SolventGenFlow.mdp_files.get("equilibrate"),
        top=SolventGenFlow.get_state_name("box", "top"),
        gro=SolventGenFlow.get_state_name("min", "gro"),
        name=SolventGenFlow.get_state_name("equ"),
        n_threads=SolventGenFlow.simulation_settings.get("n_threads"),
    )


@SolventGenFlow.pre(generated_equ_gro)
@SolventGenFlow.post(generated_prod_gro)
@SolventGenFlow.operation(cmd=True, with_job=True)
@SolventGenFlow.log_gromacs_simulation(SolventGenFlow.get_state_name("prod"))
def run_production(job):
    job.document["solvent_gro"] = SolventGenFlow.get_state_name("prod", "gro")
    return gromacs_simulation_command(
        mdp=SolventGenFlow.mdp_files.get("production"),
        top=SolventGenFlow.get_state_name("box", "top"),
        gro=SolventGenFlow.get_state_name("equ", "gro"),
        name=SolventGenFlow.get_state_name("prod"),
        n_threads=SolventGenFlow.simulation_settings.get("n_threads"),
    )


@SolventGenFlow.pre(generated_prod_gro)
@SolventGenFlow.post(lambda job: job.isfile(SolventGenFlow.nomad_workflow))
@SolventGenFlow.operation(with_job=True)
def generate_nomad_workflow(job):
    workflow = NomadWorkflow(project, job)
    workflow.build_workflow_yaml(SolventGenFlow.nomad_workflow)


@SolventGenFlow.pre(lambda job: job.isfile(SolventGenFlow.nomad_workflow))
@SolventGenFlow.post(lambda job: uploaded_to_nomad(job))
@SolventGenFlow.operation(with_job=True)
def upload_to_nomad(job: Job):
    return SolventGenFlow.upload_to_nomad(job)


def get_solvent_job(job: Job) -> Job:
    sp = {"type": "solvent", "solvent_name": job.sp.get("solvent_name")}
    return project.open_job(sp).init()


project = SolventGenFlow.get_project(path=SolventGenFlow.workspace_path)

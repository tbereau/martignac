import logging

from martignac import config
from martignac.nomad.workflows import Job, NomadWorkflow
from martignac.parsers.gromacs_forcefields import (
    generate_gro_file_for_molecule,
    generate_itp_file_for_molecule,
    generate_top_file_for_molecule,
    get_molecule_from_name,
)
from martignac.utils.gromacs import gromacs_simulation_command
from martignac.utils.martini_flow_projects import MartiniFlowProject, fetched_from_nomad, uploaded_to_nomad

conf = config()["solute_generation"]
logger = logging.getLogger(__name__)


class SoluteGenFlow(MartiniFlowProject):
    workspace_path: str = f"{MartiniFlowProject.workspaces_path}/{conf['relative_paths']['workspaces']}"
    itp_path = f"{MartiniFlowProject.input_files_path}/{conf['relative_paths']['itp_files']}"
    mdp_path = f"{MartiniFlowProject.input_files_path}/{conf['relative_paths']['mdp_files']}"
    itp_files = {k: v.get(str) for k, v in conf["itp_files"].items()}
    mdp_files = {k: v.get(str) for k, v in conf["mdp_files"].items()}
    simulation_settings = {"n_threads": conf["settings"]["n_threads"].get(int)}
    system_name = conf["output_names"]["system"].get(str)
    nomad_workflow: str = conf["output_names"]["nomad_workflow"].get(str)
    state_names = {k: v.get(str) for k, v in conf["output_names"]["states"].items()}
    ff_parameters = {
        "number_excl": conf["parameters"]["number_excl"].get(int),
        "bond_length": conf["parameters"]["bond_length"].get(float),
        "bond_constant": conf["parameters"]["bond_constant"].get(float),
    }


@SoluteGenFlow.label
def generated(job: Job):
    return (
        job.isfile(SoluteGenFlow.get_state_name("generate", "gro"))
        and job.isfile(SoluteGenFlow.get_state_name("", "itp"))
        and job.isfile(SoluteGenFlow.get_state_name("", "top"))
    )


@SoluteGenFlow.label
def minimized(job: Job):
    return job.isfile(SoluteGenFlow.get_state_name("minimize", "gro"))


@SoluteGenFlow.pre(fetched_from_nomad)
@SoluteGenFlow.post(generated)
@SoluteGenFlow.operation(with_job=True)
def solvate(job: Job) -> None:
    molecule = get_molecule_from_name(
        job.sp.solute_name,
        bond_length=SoluteGenFlow.ff_parameters["bond_length"],
        bond_constant=SoluteGenFlow.ff_parameters["bond_constant"],
        number_excl=SoluteGenFlow.ff_parameters["number_excl"],
    )
    generate_gro_file_for_molecule(molecule, SoluteGenFlow.get_state_name("generate", "gro"))
    generate_itp_file_for_molecule(molecule, SoluteGenFlow.get_state_name("", "itp"))
    generate_top_file_for_molecule(
        molecule,
        [*SoluteGenFlow.itp_files.values(), SoluteGenFlow.get_state_name("", "itp")],
        SoluteGenFlow.get_state_name("", "top"),
    )
    job.document["solute_itp"] = SoluteGenFlow.get_state_name("", "itp")
    job.document["solute_top"] = SoluteGenFlow.get_state_name("", "top")
    job.document["solute_name"] = molecule.name
    return None


@SoluteGenFlow.pre(generated)
@SoluteGenFlow.post(minimized)
@SoluteGenFlow.operation(cmd=True, with_job=True)
@SoluteGenFlow.log_gromacs_simulation(SoluteGenFlow.get_state_name("minimize"))
def minimize(job: Job):
    return gromacs_simulation_command(
        mdp=SoluteGenFlow.mdp_files["minimize"],
        top=SoluteGenFlow.get_state_name("", "top"),
        gro=SoluteGenFlow.get_state_name("generate", "gro"),
        name=SoluteGenFlow.get_state_name("minimize"),
        n_threads=SoluteGenFlow.simulation_settings["n_threads"],
    )


@SoluteGenFlow.label
def equilibrated(job: Job) -> bool:
    return job.isfile(SoluteGenFlow.get_state_name("equilibrate", "gro"))


@SoluteGenFlow.pre(minimized)
@SoluteGenFlow.post(equilibrated)
@SoluteGenFlow.operation(cmd=True, with_job=True)
@SoluteGenFlow.log_gromacs_simulation(SoluteGenFlow.get_state_name("equilibrate"))
def equilibrate(job: Job):
    job.document["solute_gro"] = SoluteGenFlow.get_state_name("equilibrate", "gro")
    return gromacs_simulation_command(
        mdp=SoluteGenFlow.mdp_files["equilibrate"],
        top=SoluteGenFlow.get_state_name("", "top"),
        gro=SoluteGenFlow.get_state_name("minimize", "gro"),
        name=SoluteGenFlow.get_state_name("equilibrate"),
        n_threads=SoluteGenFlow.simulation_settings["n_threads"],
    )


@SoluteGenFlow.pre(equilibrated)
@SoluteGenFlow.post(lambda job: job.isfile(SoluteGenFlow.nomad_workflow))
@SoluteGenFlow.operation(with_job=True)
def generate_nomad_workflow(job):
    workflow = NomadWorkflow(project, job)
    workflow.build_workflow_yaml(SoluteGenFlow.nomad_workflow)


@SoluteGenFlow.pre(lambda job: job.isfile(SoluteGenFlow.nomad_workflow))
@SoluteGenFlow.post(lambda job: uploaded_to_nomad(job))
@SoluteGenFlow.operation(with_job=True)
def upload_to_nomad(job: Job):
    return SoluteGenFlow.upload_to_nomad(job)


def get_solute_job(job: Job) -> Job:
    sp = {"type": "solute", "solute_name": job.sp.get("solute_name")}
    return project.open_job(sp).init()


project = SoluteGenFlow.get_project(path=SoluteGenFlow.workspace_path)

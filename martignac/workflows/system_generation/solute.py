from martignac import config
from martignac.nomad.workflows import Job, build_nomad_workflow, nomad_workflow_is_built
from martignac.parsers.gromacs_forcefields import (
    generate_gro_file_for_molecule,
    generate_itp_file_for_molecule,
    generate_top_file_for_molecule,
    get_molecule_from_name,
)
from martignac.utils.gromacs import gromacs_simulation_command
from martignac.utils.martini_flow_projects import (
    MartiniFlowProject,
    fetched_from_nomad,
    flag_ready_for_upload,
    store_gromacs_log_to_doc,
    store_task,
    uploaded_to_nomad,
)

conf = config()["solute_generation"]


class SoluteGenFlow(MartiniFlowProject):
    workspace_path: str = f"{MartiniFlowProject.workspaces_path}/{conf['relative_paths']['workspaces']}"
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


project_name = SoluteGenFlow.class_name()


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
@SoluteGenFlow.operation_hooks.on_success(store_task)
@SoluteGenFlow.operation(with_job=True)
def build(job: Job) -> None:
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
    job.doc[project_name]["solute_itp"] = SoluteGenFlow.get_state_name("", "itp")
    job.doc[project_name]["solute_top"] = SoluteGenFlow.get_state_name("", "top")
    job.doc[project_name]["solute_name"] = molecule.name
    return None


@SoluteGenFlow.pre(generated)
@SoluteGenFlow.post(minimized, tag="minimized")
@SoluteGenFlow.operation_hooks.on_success(store_gromacs_log_to_doc)
@SoluteGenFlow.operation(cmd=True, with_job=True)
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


@SoluteGenFlow.pre(minimized, tag="minimized")
@SoluteGenFlow.post(equilibrated, tag="equilibrated")
@SoluteGenFlow.operation_hooks.on_success(store_gromacs_log_to_doc)
@SoluteGenFlow.operation(cmd=True, with_job=True)
def equilibrate(job: Job):
    job.doc[project_name]["solute_gro"] = SoluteGenFlow.get_state_name("equilibrate", "gro")
    return gromacs_simulation_command(
        mdp=SoluteGenFlow.mdp_files["equilibrate"],
        top=SoluteGenFlow.get_state_name("", "top"),
        gro=SoluteGenFlow.get_state_name("minimize", "gro"),
        name=SoluteGenFlow.get_state_name("equilibrate"),
        n_threads=SoluteGenFlow.simulation_settings["n_threads"],
    )


@SoluteGenFlow.pre(equilibrated, tag="equilibrated")
@SoluteGenFlow.post(nomad_workflow_is_built)
@SoluteGenFlow.operation_hooks.on_success(flag_ready_for_upload)
@SoluteGenFlow.operation(with_job=True)
def generate_nomad_workflow(job):
    build_nomad_workflow(job, is_top_level=False)


@SoluteGenFlow.pre(nomad_workflow_is_built)
@SoluteGenFlow.post(uploaded_to_nomad)
@SoluteGenFlow.operation(with_job=True)
def upload_to_nomad(job: Job):
    return SoluteGenFlow.upload_to_nomad(job)


def get_solute_job(job: Job) -> Job:
    sp = {"type": "solute", "solute_name": job.sp.get("solute_name")}
    return project.open_job(sp).init()


project = SoluteGenFlow.get_project(path=SoluteGenFlow.workspace_path)

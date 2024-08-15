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
    system_equilibrated,
    system_minimized,
    uploaded_to_nomad,
)

conf = config()["solute_generation"]


class SoluteGenFlow(MartiniFlowProject):
    """
    Manages the workflow for generating, minimizing, and equilibrating solute systems for molecular dynamics simulations.

    This class extends `MartiniFlowProject` to provide specific functionalities for solute system preparation in the
    context of molecular dynamics simulations using GROMACS. It includes operations for generating the initial molecular
    structure, topology, and parameter files, minimizing the system's energy, and equilibrating the system under desired
    conditions.

    Attributes:
        workspace_path (str): Path to the workspace directory where simulation files are stored.
        mdp_path (str): Path to the directory containing MDP files for GROMACS simulations.
        itp_files (dict): Dictionary mapping solute names to their respective ITP file paths.
        mdp_files (dict): Dictionary mapping simulation types (e.g., 'minimize', 'equilibrate') to their MDP file paths.
        simulation_settings (dict): Settings for the simulation, such as the number of threads.
        system_name (str): The name of the system being simulated.
        nomad_workflow (str): The name of the NOMAD workflow associated with this project.
        state_names (dict): Dictionary mapping state names to their string representations.
        ff_parameters (dict): Force field parameters used in generating the molecular structure.

    The class provides methods for each step of the solute preparation process, including `build`, `minimize`, and
    `equilibrate`, as well as methods for managing the workflow's integration with the NOMAD database, such as
    `upload_to_nomad`.
    """

    workspace_path: str = (
        f"{MartiniFlowProject.workspaces_path}/{conf['relative_paths']['workspaces']}"
    )
    mdp_path = (
        f"{MartiniFlowProject.input_files_path}/{conf['relative_paths']['mdp_files']}"
    )
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
def system_generated(job: Job):
    return (
        job.isfile(SoluteGenFlow.get_state_name("generate", "gro"))
        and job.isfile(SoluteGenFlow.get_state_name("", "itp"))
        and job.isfile(SoluteGenFlow.get_state_name("", "top"))
    )


@SoluteGenFlow.pre(fetched_from_nomad)
@SoluteGenFlow.post(system_generated, tag="system_generated")
@SoluteGenFlow.operation_hooks.on_success(store_task)
@SoluteGenFlow.operation(with_job=True)
def build(job: Job) -> None:
    """
    Generate the initial structure, topology, and parameter files for a solute.

    This function is responsible for generating the initial molecular structure (GRO file),
    topology (TOP file), and parameter (ITP file) files for the solute specified in the job.
    It utilizes the `get_molecule_from_name` function to fetch the molecular data based on the
    solute name provided in the job's state point. The generated files are essential for
    conducting molecular dynamics simulations using GROMACS.

    The function updates the job document with paths to the generated files, marking the
    solute's initial system generation as complete. This step is crucial for preparing the
    system for subsequent energy minimization, equilibration, and production simulations.

    Args:
        job (Job): The job object representing the solute system for which the initial
                   structure, topology, and parameter files are to be generated.

    Returns:
        None: This function does not return a value but updates the job document with the
              paths to the generated files.
    """
    molecule = get_molecule_from_name(
        job.sp.solute_name,
        bond_length=SoluteGenFlow.ff_parameters["bond_length"],
        bond_constant=SoluteGenFlow.ff_parameters["bond_constant"],
        number_excl=SoluteGenFlow.ff_parameters["number_excl"],
    )
    generate_gro_file_for_molecule(
        molecule, SoluteGenFlow.get_state_name("generate", "gro")
    )
    generate_itp_file_for_molecule(molecule, SoluteGenFlow.get_state_name("", "itp"))
    generate_top_file_for_molecule(
        molecule,
        [*SoluteGenFlow.itp_files.values(), SoluteGenFlow.get_state_name("", "itp")],
        SoluteGenFlow.get_state_name("", "top"),
    )
    job.doc[project_name]["solute_itp"] = SoluteGenFlow.get_state_name("", "itp")
    job.doc[project_name]["solute_top"] = SoluteGenFlow.get_state_name("", "top")
    job.doc[project_name]["solute_name"] = molecule.name
    job.doc[project_name]["solute_has_charged_beads"] = molecule.has_charged_beads
    return None


@SoluteGenFlow.pre(system_generated, tag="system_generated")
@SoluteGenFlow.post(system_minimized, tag="system_minimized")
@SoluteGenFlow.operation_hooks.on_success(store_gromacs_log_to_doc)
@SoluteGenFlow.operation(cmd=True, with_job=True)
def minimize(job: Job) -> str:
    """
    Minimize the energy of the solute system.

    This function executes the energy minimization step for the solute using GROMACS. It utilizes the parameters
    specified in the `mdp_files["minimize"]` file for the minimization process. The minimized structure serves as a
    more stable starting point for subsequent equilibration and production simulations.

    The minimization process reduces the potential energy of the system, resolving any steric clashes or unrealistic
    geometries that may have arisen during the solute generation phase. This step is crucial for preparing the system
    for a stable simulation environment.

    Upon completion, the function updates the job document with the path to the minimized structure, indicating that
    the system has been successfully minimized.

    Args:
        job (Job): The job object representing the solute system to be minimized.

    Returns:
        str: The command to execute the GROMACS energy minimization.
    """
    return gromacs_simulation_command(
        mdp=SoluteGenFlow.mdp_files["minimize"],
        top=SoluteGenFlow.get_state_name("", "top"),
        gro=SoluteGenFlow.get_state_name("generate", "gro"),
        name=SoluteGenFlow.get_state_name("minimize"),
        n_threads=SoluteGenFlow.simulation_settings["n_threads"],
    )


@SoluteGenFlow.pre(system_minimized, tag="system_minimized")
@SoluteGenFlow.post(system_equilibrated, tag="system_equilibrated")
@SoluteGenFlow.operation_hooks.on_success(store_gromacs_log_to_doc)
@SoluteGenFlow.operation(cmd=True, with_job=True)
def equilibrate(job: Job) -> str:
    """
    Perform the equilibration process for the solute in the simulation.

    This function runs the GROMACS equilibration simulation using the parameters defined in the
    `mdp_files["equilibrate"]` file. It sets up the simulation with the minimized structure from the previous step as
    the input configuration and uses the topology file generated during the solute generation step. The number of
    threads for the simulation is determined by the `simulation_settings["n_threads"]` configuration.

    The equilibration process is crucial for stabilizing the system at the desired temperature and pressure before the
    main production run. It helps in relaxing the system and achieving a more realistic distribution of velocities and
    positions.

    Upon successful completion, the function updates the job document with the path to the generated GRO file, marking
    the state of the system as equilibrated.

    Args:
        job (Job): The job object for which the equilibration simulation is to be performed.

    Returns:
        str: The command to execute the GROMACS equilibration simulation.
    """
    job.doc[project_name]["solute_gro"] = SoluteGenFlow.get_state_name(
        "equilibrate", "gro"
    )
    return gromacs_simulation_command(
        mdp=SoluteGenFlow.mdp_files["equilibrate"],
        top=SoluteGenFlow.get_state_name("", "top"),
        gro=SoluteGenFlow.get_state_name("minimize", "gro"),
        name=SoluteGenFlow.get_state_name("equilibrate"),
        n_threads=SoluteGenFlow.simulation_settings["n_threads"],
    )


@SoluteGenFlow.pre(system_equilibrated, tag="system_equilibrated")
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

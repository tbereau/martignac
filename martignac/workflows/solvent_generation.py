import numpy as np

from martignac import config
from martignac.nomad.workflows import Job, build_nomad_workflow, nomad_workflow_is_built
from martignac.parsers.gromacs_forcefields import (
    find_molecule_from_name,
    generate_gro_file_for_molecule,
    generate_top_file_for_molecule,
)
from martignac.parsers.gromacs_structures import get_number_of_molecules_from_gro
from martignac.utils.gromacs import (
    gromacs_simulation_command,
)
from martignac.utils.martini_flow_projects import (
    MartiniFlowProject,
    fetched_from_nomad,
    flag_ready_for_upload,
    store_gromacs_log_to_doc,
    store_task,
    system_equilibrated,
    system_minimized,
    system_sampled,
    uploaded_to_nomad,
)
from martignac.utils.misc import convert_pdb_to_gro
from martignac.utils.packmol import generate_solvent_with_packmol

conf = config()["solvent_generation"]


class SolventGenFlow(MartiniFlowProject):
    """
    A workflow class for generating and preparing solvent systems for molecular dynamics simulations.

    This class extends `MartiniFlowProject` to manage the workflow specific to solvent system generation, including
    the creation of solvent molecules, building solvent boxes, and preparing these systems for simulations using
    GROMACS. It encompasses operations such as generating the initial molecular structure (GRO file), topology (TOP file),
    and parameter (ITP file) files for the solvent, followed by energy minimization, equilibration, and production
    simulations.

    Attributes:
        workspace_path (str): Path to the workspace directory where simulation files are stored.
        mdp_path (str): Path to the directory containing MDP files for GROMACS simulations.
        itp_files (dict): Dictionary mapping solvent names to their respective ITP file paths.
        mdp_files (dict): Dictionary mapping simulation types (e.g., 'minimize', 'equilibrate', 'production') to their MDP file paths.
        simulation_settings (dict): Settings for the simulation, such as the number of threads and box length.
        system_name (str): The name of the system being simulated.
        nomad_workflow (str): The name of the NOMAD workflow associated with this project.
        state_names (dict): Dictionary mapping state names to their string representations for tracking the simulation progress.

    The class provides methods for each step of the solvent preparation process, including generating solvent molecules,
    building solvent boxes, converting box formats, minimizing energy, equilibrating the system, running production
    simulations, and integrating with the NOMAD database for workflow management.
    """

    workspace_path: str = (
        f"{MartiniFlowProject.workspaces_path}/{conf['relative_paths']['workspaces']}"
    )
    mdp_path = (
        f"{MartiniFlowProject.input_files_path}/{conf['relative_paths']['mdp_files']}"
    )
    itp_files = {k: v.get(str) for k, v in conf["itp_files"].items()}
    mdp_files = {k: v.get(str) for k, v in conf["mdp_files"].items()}
    simulation_settings = {
        "n_threads": conf["settings"]["n_threads"].get(int),
        "box_length": conf["settings"]["box_length"].get(float),
    }
    system_name = conf["output_names"]["system"].get(str)
    nomad_workflow: str = conf["output_names"]["nomad_workflow"].get(str)
    state_names = {k: v.get(str) for k, v in conf["output_names"]["states"].items()}


project_name = SolventGenFlow.class_name()


@SolventGenFlow.label
def generated_box_pdb(job):
    return job.isfile(SolventGenFlow.get_state_name("box", "pdb"))


@SolventGenFlow.label
def generated_mol_gro(job):
    return job.isfile(SolventGenFlow.get_state_name("mol", "gro"))


@SolventGenFlow.label
def generated_box_gro(job):
    return job.isfile(SolventGenFlow.get_state_name("box", "gro"))


@SolventGenFlow.pre(fetched_from_nomad)
@SolventGenFlow.post(generated_mol_gro, tag="generated_mol_gro")
@SolventGenFlow.operation_hooks.on_success(store_task)
@SolventGenFlow.operation(with_job=True)
def generate_solvent_molecule(job) -> None:
    """
    Generates the initial molecular structure (GRO file) and topology (TOP file) for a specified solvent molecule.

    This function is a critical step in the solvent generation workflow, where it takes the name of the solvent
    from the job's state point and utilizes it to find the corresponding molecular data. It then generates the
    GRO and TOP files necessary for molecular dynamics simulations using GROMACS. These files represent the
    initial molecular structure and topology of the solvent molecule, respectively.

    The generation of these files marks the completion of the initial setup for the solvent molecule, allowing
    for further steps in the workflow, such as building the solvent box and performing energy minimization.

    Args:
        job (Job): The job object associated with the current workflow operation, containing the solvent name
                   and other relevant information in its state point.

    Returns:
        None: This function does not return a value but updates the job document with the paths to the generated
              GRO and TOP files, marking the 'generated_mol_gro' condition as complete.
    """
    molecule = find_molecule_from_name(
        list(SolventGenFlow.itp_files.values()), job.sp.solvent_name
    )
    generate_gro_file_for_molecule(
        molecule, SolventGenFlow.get_state_name("mol", "gro")
    )
    generate_top_file_for_molecule(
        molecule,
        list(SolventGenFlow.itp_files.values()),
        SolventGenFlow.get_state_name("mol", "top"),
    )
    return None


@SolventGenFlow.pre(generated_mol_gro, tag="generated_mol_gro")
@SolventGenFlow.post(generated_box_pdb, tag="generated_box_pdb")
@SolventGenFlow.operation_hooks.on_success(store_task)
@SolventGenFlow.operation(cmd=True, with_job=True)
def build_solvent_box(_):
    """
    Builds a solvent box using the generated molecular structure of the solvent.

    This operation is a key step in the solvent generation workflow, where it utilizes the previously generated
    GRO file of the solvent molecule to create a solvent box suitable for molecular dynamics simulations. The
    solvent box is generated using the Packmol tool, which arranges multiple copies of the solvent molecule
    within a defined box size to achieve a desired density.

    The function does not take any arguments directly; instead, it operates on the state of the workflow,
    specifically looking for the presence of a generated GRO file for the solvent molecule. Upon successful
    completion, it generates a PDB file representing the solvent box, marking the 'generated_box_pdb' condition
    as complete.

    This step is crucial for preparing the solvent system for subsequent energy minimization and equilibration
    processes, ensuring that the system is ready for molecular dynamics simulations.

    Note:
        This function is designed to be executed within the SolventGenFlow workflow and relies on the workflow's
        state management to control its execution and dependencies.
    """
    return generate_solvent_with_packmol(
        gro_solvent_mol=SolventGenFlow.get_state_name("mol", "gro"),
        box_length=SolventGenFlow.simulation_settings.get("box_length"),
        output_pdb=SolventGenFlow.get_state_name("box", "pdb"),
    )


@SolventGenFlow.pre(generated_box_pdb, tag="generated_box_pdb")
@SolventGenFlow.post(generated_box_gro, tag="generated_box_gro")
@SolventGenFlow.operation_hooks.on_success(store_task)
@SolventGenFlow.operation(with_job=True)
def convert_box_to_gro(_):
    """
    Converts the solvent box from PDB to GRO format and adjusts the box dimensions.

    This function is part of the solvent generation workflow in molecular dynamics simulations. It takes the solvent
    box in PDB format, generated by the `build_solvent_box` function, and converts it to GRO format, which is required
    for further simulation steps in GROMACS. Additionally, it adjusts the box dimensions based on the specified box
    length in the simulation settings, ensuring the correct size and orientation for the simulation.

    The conversion and adjustment of the box dimensions are crucial for accurate molecular dynamics simulations,
    as they directly affect the system's physical properties and behavior during the simulation.

    Note:
        This function does not take any arguments directly; instead, it operates on the state of the workflow,
        specifically looking for the presence of a generated PDB file for the solvent box. It updates the workflow's
        state by marking the 'generated_box_gro' condition as complete upon successful conversion.
    """
    box_length = SolventGenFlow.simulation_settings.get("box_length")
    box_vector = np.array(
        [box_length * 10.0, box_length * 10.0, box_length * 10.0, 90.0, 90.0, 90.0]
    )
    convert_pdb_to_gro(
        SolventGenFlow.get_state_name("box", "pdb"),
        SolventGenFlow.get_state_name("box", "gro"),
        box_vector,
    )


@SolventGenFlow.pre(generated_box_gro, tag="generated_box_gro")
@SolventGenFlow.post(system_minimized, tag="system_minimized")
@SolventGenFlow.operation_hooks.on_success(store_gromacs_log_to_doc)
@SolventGenFlow.operation(cmd=True, with_job=True)
def minimize(job):
    """
    Performs energy minimization on the solvent system.

    This function is a crucial step in the solvent generation workflow, aimed at minimizing the potential energy of the
    solvent system to ensure stability before proceeding to equilibration and production simulations. It generates the
    topology file for the solvent box with the correct number of molecules, updates the job document with the solvent
    topology and name, and then executes the GROMACS energy minimization command using the specified MDP file.

    The energy minimization process is essential for removing any physically unrealistic conformations and reducing the
    system's energy, which could otherwise lead to simulation artifacts or instabilities.

    Args:
        job (Job): The job object associated with the current workflow operation, containing the solvent name
                   and other relevant information in its state point.

    Returns:
        None: This function does not return a value but updates the job document with the paths to the generated
              topology file and executes the energy minimization command.
    """
    molecule = find_molecule_from_name(
        list(SolventGenFlow.itp_files.values()), job.sp.solvent_name
    )
    generate_top_file_for_molecule(
        molecule,
        list(SolventGenFlow.itp_files.values()),
        SolventGenFlow.get_state_name("box", "top"),
        num_molecules=get_number_of_molecules_from_gro(
            SolventGenFlow.get_state_name("box", "gro")
        ),
    )
    job.doc[project_name]["solvent_top"] = SolventGenFlow.get_state_name("box", "top")
    job.doc[project_name]["solvent_name"] = job.sp.solvent_name
    return gromacs_simulation_command(
        mdp=SolventGenFlow.mdp_files.get("minimize"),
        top=SolventGenFlow.get_state_name("box", "top"),
        gro=SolventGenFlow.get_state_name("box", "gro"),
        name=SolventGenFlow.get_state_name("minimize"),
        n_threads=SolventGenFlow.simulation_settings.get("n_threads"),
    )


@SolventGenFlow.pre(system_minimized, tag="system_minimized")
@SolventGenFlow.post(system_equilibrated, tag="system_equilibrated")
@SolventGenFlow.operation_hooks.on_success(store_gromacs_log_to_doc)
@SolventGenFlow.operation(cmd=True, with_job=True)
def equilibrate(_):
    """
    Performs the equilibration phase of the molecular dynamics simulation.

    This function is responsible for running the equilibration simulation step, which is crucial for stabilizing the
    system at the desired temperature and pressure before the production run. It uses the GROMACS simulation command
    configured with the 'equilibrate' MDP file, which contains the parameters for the equilibration process.

    The equilibration step ensures that the system's temperature, pressure, and density are stabilized, preparing it
    for the accurate and reliable collection of production data in the subsequent simulation phase.

    Note:
        This function does not take any arguments directly; instead, it operates on the state of the workflow,
        specifically looking for the presence of a minimized GRO file for the solvent box. It updates the workflow's
        state by marking the 'system_equilibrated' condition as complete upon successful equilibration.

    Returns:
        None: This function does not return a value but executes the equilibration simulation command and updates
              the workflow's state accordingly.
    """
    return gromacs_simulation_command(
        mdp=SolventGenFlow.mdp_files.get("equilibrate"),
        top=SolventGenFlow.get_state_name("box", "top"),
        gro=SolventGenFlow.get_state_name("minimize", "gro"),
        name=SolventGenFlow.get_state_name("equilibrate"),
        n_threads=SolventGenFlow.simulation_settings.get("n_threads"),
    )


@SolventGenFlow.pre(system_equilibrated, tag="system_equilibrated")
@SolventGenFlow.post(system_sampled, tag="system_sampled")
@SolventGenFlow.operation_hooks.on_success(store_gromacs_log_to_doc)
@SolventGenFlow.operation(cmd=True, with_job=True)
def production(job):
    """
    Executes the production phase of the molecular dynamics simulation.

    This function is the final step in the solvent generation workflow, where it performs the production simulation
    using GROMACS. The production simulation is crucial for generating the actual data for analysis, representing
    the behavior of the solvent system under specified conditions over time.

    The function configures and executes the GROMACS simulation command with the 'production' MDP file, which contains
    the parameters for the production phase, such as simulation time and temperature control settings. It updates the
    job document with the path to the final GRO file generated by this simulation, marking the 'system_sampled' condition
    as complete, indicating that the production phase has been successfully executed and data collection is finished.

    Args:
        job (Job): The job object associated with the current workflow operation, containing the solvent name
                   and other relevant information in its state point.

    Returns:
        None: This function does not return a value but updates the job document with the path to the production
              GRO file and executes the production simulation command.
    """
    job.doc[project_name]["solvent_gro"] = SolventGenFlow.get_state_name(
        "production", "gro"
    )
    return gromacs_simulation_command(
        mdp=SolventGenFlow.mdp_files.get("production"),
        top=SolventGenFlow.get_state_name("box", "top"),
        gro=SolventGenFlow.get_state_name("equilibrate", "gro"),
        name=SolventGenFlow.get_state_name("production"),
        n_threads=SolventGenFlow.simulation_settings.get("n_threads"),
    )


@SolventGenFlow.pre(system_sampled, tag="system_sampled")
@SolventGenFlow.post(nomad_workflow_is_built)
@SolventGenFlow.operation_hooks.on_success(flag_ready_for_upload)
@SolventGenFlow.operation(with_job=True)
def generate_nomad_workflow(job):
    build_nomad_workflow(job, is_top_level=False)


@SolventGenFlow.pre(nomad_workflow_is_built)
@SolventGenFlow.post(uploaded_to_nomad)
@SolventGenFlow.operation(with_job=True)
def upload_to_nomad(job: Job):
    return SolventGenFlow.upload_to_nomad(job)


def get_solvent_job(job: Job) -> Job:
    sp = {"type": "solvent", "solvent_name": job.sp.get("solvent_name")}
    return project.open_job(sp).init()


project = SolventGenFlow.get_project(path=SolventGenFlow.workspace_path)

from martignac import config
from martignac.liquid_models.mixtures import LiquidMixture
from martignac.nomad.workflows import Job, build_nomad_workflow, nomad_workflow_is_built
from martignac.parsers.gromacs_topologies import Topology
from martignac.utils.gromacs import gromacs_simulation_command
from martignac.utils.insane import generate_bilayer_with_insane
from martignac.utils.martini_flow_projects import (
    MartiniFlowProject,
    fetched_from_nomad,
    flag_ready_for_upload,
    store_gromacs_log_to_doc,
    system_equilibrated,
    system_generated,
    system_minimized,
    system_sampled,
    uploaded_to_nomad,
)
from martignac.utils.misc import sub_template_mdp

conf = config()["bilayer_generation"]


class BilayerGenFlow(MartiniFlowProject):
    """
    A workflow for generating, minimizing, equilibrating, and sampling molecular dynamics simulations of lipid bilayers.

    This class extends `MartiniFlowProject` to provide a specialized workflow for the generation and analysis of lipid
    bilayer systems using the MARTINI force field. It includes operations for the initial generation of the bilayer,
    energy minimization, equilibration, and production phase sampling. The workflow is designed to be used with the
    GROMACS simulation package and utilizes various utilities for file manipulation and simulation setup.

    Attributes:
        workspace_path (str): The path to the workspace directory where simulation files are stored.
        mdp_path (str): The path to the directory containing MDP (Molecular Dynamics Parameters) files.
        itp_files (dict): A dictionary mapping from descriptive names to ITP (Include Topology) file paths.
        mdp_files (dict): A dictionary mapping from simulation stages to their corresponding MDP file paths.
        simulation_settings (dict): A dictionary containing simulation settings such as the number of threads,
                                    and box dimensions.
        solvent (LiquidMixture): An instance of `LiquidMixture` representing the solvent in the simulation.
        system_name (str): The name of the system being simulated.
        nomad_workflow (str): The name of the NOMAD workflow associated with this project.
        state_names (dict): A dictionary mapping from state descriptions to state names used in file naming.

    The workflow is configured through a YAML configuration file, which specifies paths, file names, and simulation
    parameters. This class provides methods to generate the initial bilayer structure, perform energy minimization,
    equilibrate the system, and carry out the production phase of the simulation. It also includes methods for
    integrating with the NOMAD database for data storage and retrieval.
    """

    workspace_path: str = f"{MartiniFlowProject.workspaces_path}/{conf['relative_paths']['workspaces']}"
    mdp_path = f"{MartiniFlowProject.input_files_path}/{conf['relative_paths']['mdp_files']}"
    itp_files = {k: v.get(str) for k, v in conf["itp_files"].items()}
    mdp_files = {k: v.get(str) for k, v in conf["mdp_files"].items()}
    simulation_settings = {
        "n_threads": conf["settings"]["n_threads"].get(int),
        "box_length_xy": conf["settings"]["box_length_xy"].get(float),
        "box_length_z": conf["settings"]["box_length_z"].get(float),
    }
    solvent: LiquidMixture = LiquidMixture.from_list_of_dicts([s.get(dict) for s in conf["settings"]["solvent"]])
    system_name = conf["output_names"]["system"].get(str)
    nomad_workflow: str = conf["output_names"]["nomad_workflow"].get(str)
    state_names = {k: v.get(str) for k, v in conf["output_names"]["states"].items()}


project_name = BilayerGenFlow.class_name()


def lipid_names(job) -> list[str]:
    return [
        lip.get("name").removeprefix("M3.").removeprefix("M2.").removeprefix("M2o.") for lip in job.sp.get("lipids")
    ]


@BilayerGenFlow.pre(fetched_from_nomad)
@BilayerGenFlow.post(system_generated, tag="system_generated")
@BilayerGenFlow.operation(cmd=True, with_job=True)
def generate_initial_bilayer(job):
    """
    Generates the initial structure of a lipid bilayer for molecular dynamics simulations.

    This function utilizes the INSANE (INsertion of Solutes At specified Norms and Extensions) tool to generate
    a lipid bilayer structure based on the specifications provided in the job's state point. The generated structure
    includes the specified lipid mixture and solvent, arranged according to the given box dimensions. The resulting
    structure is saved in a GRO file, and the topology of the bilayer is updated to reflect the composition of the
    generated system.

    The function updates the job's state point with the solvent information in a format compatible with INSANE,
    ensuring that the solvent is correctly represented in the generated bilayer structure.

    Args:
        job (Job): The job object associated with the current workflow operation. It contains the state point
                   information specifying the lipid composition, solvent, and box dimensions for the bilayer
                   generation.

    Returns:
        str: The command executed to generate the bilayer with INSANE, primarily for logging or debugging purposes.
    """
    lipid_mixture = LiquidMixture.from_list_of_dicts(job.sp.lipids)
    return generate_bilayer_with_insane(
        lipids=lipid_mixture,
        solvent=BilayerGenFlow.solvent,
        box_length_xy=BilayerGenFlow.simulation_settings.get("box_length_xy"),
        box_length_z=BilayerGenFlow.simulation_settings.get("box_length_z"),
        gro_bilayer_gen=BilayerGenFlow.get_state_name("generate", "gro"),
        top_bilayer=BilayerGenFlow.get_state_name(extension="top"),
    )


@BilayerGenFlow.pre(system_generated, tag="system_generated")
@BilayerGenFlow.post(system_minimized, tag="system_minimized")
@BilayerGenFlow.operation_hooks.on_success(store_gromacs_log_to_doc)
@BilayerGenFlow.operation(cmd=True, with_job=True)
def minimize(_):
    """
    Performs energy minimization on the initial lipid bilayer structure.

    This function is a critical step in the molecular dynamics simulation workflow, aiming to relax the system by
    minimizing its potential energy. This process helps in stabilizing the initial structure before proceeding to
    more computationally intensive equilibration and production phases. The function modifies the topology file to
    include necessary ITP (Include Topology) files, ensuring that all molecular interactions are correctly defined
    for the minimization process.

    The GROMACS simulation tool is invoked with parameters specified in the corresponding MDP (Molecular Dynamics
    Parameters) file for minimization. The output includes a minimized structure in GRO format and an updated topology
    file, both of which are essential for subsequent simulation steps.

    Args:
        _ (Job): The job object associated with the current workflow operation. It contains the state point
                 information and provides context for the simulation, including paths to input and output files.
                 The underscore is used to indicate that this parameter is not directly used within the function,
                 but is required for compatibility with the workflow's operation structure.

    Returns:
        str: The command executed for the energy minimization process, primarily for logging or debugging purposes.
    """
    top = Topology.parse_top_file(BilayerGenFlow.get_state_name(extension="top"))
    top.includes = BilayerGenFlow.itp_files.values()
    top.output_top(BilayerGenFlow.get_state_name(extension="top"))
    return gromacs_simulation_command(
        mdp=BilayerGenFlow.mdp_files.get("minimize"),
        top=BilayerGenFlow.get_state_name(extension="top"),
        gro=BilayerGenFlow.get_state_name("generate", "gro"),
        name=BilayerGenFlow.get_state_name("minimize"),
        n_threads=BilayerGenFlow.simulation_settings.get("n_threads"),
        verbose=False,
    )


@BilayerGenFlow.pre(system_minimized, tag="system_minimized")
@BilayerGenFlow.post(system_equilibrated, tag="system_equilibrated")
@BilayerGenFlow.operation_hooks.on_success(store_gromacs_log_to_doc)
@BilayerGenFlow.operation(cmd=True, with_job=True)
def equilibrate(job):
    """
    Equilibrates the minimized lipid bilayer structure for molecular dynamics simulations.

    This function prepares and executes the equilibration phase of the molecular dynamics simulation workflow,
    following the energy minimization step. It involves modifying the MDP (Molecular Dynamics Parameters) file
    to include specific lipid and solvent names, which are then used to run the equilibration simulation with
    GROMACS. The equilibration process allows the system to reach a stable state before proceeding to the production
    phase of the simulation.

    The MDP file for equilibration is customized based on the job's state point, which includes the lipid composition
    and solvent information. This customization ensures that the simulation parameters are correctly set for the
    specific bilayer system being simulated.

    Args:
        job (Job): The job object associated with the current workflow operation. It contains the state point
                   information and provides context for the simulation, including paths to input and output files,
                   and the lipid and solvent composition.

    Returns:
        str: The command executed for the equilibration process, primarily for logging or debugging purposes.
    """
    mdp_file_template = BilayerGenFlow.mdp_files.get("equilibrate")
    mdp_file = job.fn(BilayerGenFlow.get_state_name("equilibrate", "mdp"))
    sub_template_mdp(mdp_file_template, "sedlipids", " ".join(lipid_names(job)), mdp_file)
    sub_template_mdp(mdp_file, "sedsolvent", " ".join(BilayerGenFlow.solvent.solvent_names), mdp_file)
    return gromacs_simulation_command(
        mdp=mdp_file,
        top=BilayerGenFlow.get_state_name(extension="top"),
        gro=BilayerGenFlow.get_state_name("minimize", "gro"),
        name=BilayerGenFlow.get_state_name("equilibrate"),
        n_threads=BilayerGenFlow.simulation_settings.get("n_threads"),
    )


@BilayerGenFlow.pre(system_equilibrated, tag="system_equilibrated")
@BilayerGenFlow.post(system_sampled, tag="system_sampled")
@BilayerGenFlow.operation_hooks.on_success(store_gromacs_log_to_doc)
@BilayerGenFlow.operation(cmd=True, with_job=True)
def production(job):
    """
    Executes the production phase of the molecular dynamics simulation for the lipid bilayer.

    This function is responsible for running the final, production phase of the molecular dynamics simulation
    workflow, following successful equilibration. It customizes the MDP (Molecular Dynamics Parameters) file for
    the production phase, incorporating specific lipid and solvent names. The simulation is executed using GROMACS,
    with the system now fully equilibrated and ready for detailed analysis and data collection.

    The MDP file for the production phase is tailored based on the job's state point, ensuring that simulation
    parameters are optimally set for the specific bilayer system under study. This phase is crucial for generating
    the simulation data that will be used for subsequent analysis and research findings.

    Args:
        job (Job): The job object associated with the current workflow operation. It contains the state point
                   information and provides context for the simulation, including paths to input and output files,
                   and the lipid and solvent composition.

    Returns:
        str: The command executed for the production phase, primarily for logging or debugging purposes. This includes
             the path to the modified MDP file, the output GRO file for the bilayer, and the updated topology file.
    """
    mdp_file_template = BilayerGenFlow.mdp_files.get("production")
    mdp_file = job.fn(BilayerGenFlow.get_state_name("production", "mdp"))
    sub_template_mdp(mdp_file_template, "sedlipids", " ".join(lipid_names(job)), mdp_file)
    sub_template_mdp(mdp_file, "sedsolvent", " ".join(BilayerGenFlow.solvent.solvent_names), mdp_file)
    job.doc[project_name]["bilayer_gro"] = BilayerGenFlow.get_state_name("production", "gro")
    job.doc[project_name]["bilayer_top"] = BilayerGenFlow.get_state_name(extension="top")
    job.doc[project_name]["lipid_names"] = lipid_names(job)
    return gromacs_simulation_command(
        mdp=mdp_file,
        top=BilayerGenFlow.get_state_name(extension="top"),
        gro=BilayerGenFlow.get_state_name("equilibrate", "gro"),
        name=BilayerGenFlow.get_state_name("production"),
        n_threads=BilayerGenFlow.simulation_settings.get("n_threads"),
    )


@BilayerGenFlow.pre(system_sampled, tag="system_sampled")
@BilayerGenFlow.post(nomad_workflow_is_built)
@BilayerGenFlow.operation_hooks.on_success(flag_ready_for_upload)
@BilayerGenFlow.operation(with_job=True)
def generate_nomad_workflow(job: Job):
    build_nomad_workflow(job, is_top_level=False)


@BilayerGenFlow.pre(nomad_workflow_is_built)
@BilayerGenFlow.post(uploaded_to_nomad)
@BilayerGenFlow.operation(with_job=True)
def upload_to_nomad(job: Job):
    print(f"comment {BilayerGenFlow.nomad_comment(job)}")
    return BilayerGenFlow.upload_to_nomad(job)


def get_solvent_job(job: Job) -> Job:
    sp = {"type": "bilayer", "lipids": job.sp.get("lipids")}
    return project.open_job(sp).init()


project = BilayerGenFlow.get_project(path=BilayerGenFlow.workspace_path)

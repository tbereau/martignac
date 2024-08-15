import logging

from martignac import config
from martignac.nomad.workflows import Job, build_nomad_workflow
from martignac.parsers.gromacs_topologies import Topology, append_all_includes_to_top
from martignac.utils.gromacs import (
    gromacs_simulation_command,
    solvate_solute_command,
)
from martignac.utils.martini_flow_projects import (
    MartiniFlowProject,
    fetched_from_nomad,
    flag_ready_for_upload,
    import_job_from_other_flow,
    is_ready_for_upload,
    store_gromacs_log_to_doc,
    store_task,
    store_workflow,
    system_equilibrated,
    system_generated,
    system_minimized,
    uploaded_to_nomad,
)
from martignac.utils.misc import func_name
from martignac.workflows.solute_generation import SoluteGenFlow, get_solute_job
from martignac.workflows.solvent_generation import SolventGenFlow, get_solvent_job

logger = logging.getLogger(__name__)

conf = config()["solute_in_solvent_generation"]


class SoluteInSolventGenFlow(MartiniFlowProject):
    """
    Represents the workflow for generating a solute in solvent system for molecular dynamics simulations.

    This class extends the MartiniFlowProject, incorporating specific configurations and operations
    for solute and solvent generation, solvation, minimization, equilibration, and preparation for
    NOMAD upload. It defines the workflow's file paths, input files, simulation settings, and state
    names based on the project configuration.

    Attributes:
        workspace_path (str): Path to the workspace directory for this workflow.
        itp_files (dict): Dictionary mapping from itp file identifiers to their paths.
        mdp_path (str): Path to the directory containing MDP files for GROMACS simulations.
        mdp_files (dict): Dictionary mapping from MDP file identifiers to their paths.
        simulation_settings (dict): Dictionary containing simulation settings, e.g., number of threads.
        system_name (str): Name of the system being generated.
        nomad_workflow (str): Name of the workflow file for NOMAD.
        nomad_top_level_workflow (str): Name of the top-level workflow file for NOMAD.
        state_names (dict): Dictionary mapping state identifiers to their names.

    The class also includes methods for each step of the workflow, including generating the solute and
    solvent, solvating the solute with the solvent, minimizing the energy of the system, equilibrating
    the system, and preparing and uploading the workflow to NOMAD.
    """

    workspace_path: str = (
        f"{MartiniFlowProject.workspaces_path}/{conf['relative_paths']['workspaces']}"
    )
    itp_files = {k: v.get(str) for k, v in conf["itp_files"].items()}
    mdp_path = (
        f"{MartiniFlowProject.input_files_path}/{conf['relative_paths']['mdp_files']}"
    )
    mdp_files = {k: v.get(str) for k, v in conf["mdp_files"].items()}
    simulation_settings = {"n_threads": conf["settings"]["n_threads"].get(int)}
    system_name = conf["output_names"]["system"].get(str)
    nomad_workflow: str = conf["output_names"]["nomad_workflow"].get(str)
    nomad_top_level_workflow: str = conf["output_names"][
        "nomad_top_level_workflow"
    ].get(str)
    state_names = {k: v.get(str) for k, v in conf["output_names"]["states"].items()}


project_name = SoluteInSolventGenFlow.class_name()
solute_gen_name = SoluteGenFlow.class_name()
solvent_gen_name = SolventGenFlow.class_name()


@SoluteInSolventGenFlow.label
def solute_generated(job) -> bool:
    return (
        solute_gen_name in job.doc
        and job.doc[solute_gen_name].get("solute_gro")
        and job.doc[solute_gen_name].get("solute_top")
    )


@SoluteInSolventGenFlow.label
def solvent_generated(job) -> bool:
    return solvent_gen_name in job.doc and job.doc[solvent_gen_name].get("solvent_gro")


@SoluteInSolventGenFlow.pre(fetched_from_nomad)
@SoluteInSolventGenFlow.post(solute_generated, tag="solute_generated")
@SoluteInSolventGenFlow.operation_hooks.on_success(store_workflow)
@SoluteInSolventGenFlow.operation
def generate_solute(job: Job):
    """
    Generates the solute for the molecular dynamics simulation.

    This operation is a critical step in the solute in solvent generation workflow. It imports the solute
    generation results from a previous workflow, specifically the solute's GRO, TOP, and ITP files, into the
    current job's document. This process ensures that the solute configuration is ready for subsequent steps,
    such as solvation and minimization.

    The function leverages the `import_job_from_other_flow` utility to copy the necessary files from the solute
    generation project to the current project, updating the job document with the paths to these files. It also
    sets the operation's name in the project's operation to workflow mapping, facilitating tracking and management
    of the workflow's progress.

    Args:
        job (Job): The job object associated with the current workflow operation. It contains the state point
                   and document for the job, which will be updated with the solute generation results.

    Returns:
        None: This function does not return a value but updates the job document with the paths to the solute's
              GRO, TOP, and ITP files, marking the solute generation step as complete.
    """
    solute_job: Job = get_solute_job(job)
    keys_for_files_to_copy = ["solute_gro", "solute_top", "solute_itp"]
    project.operation_to_workflow[func_name()] = solute_gen_name
    import_job_from_other_flow(
        job, solute_gen_project, solute_job, keys_for_files_to_copy
    )


@SoluteInSolventGenFlow.pre(fetched_from_nomad)
@SoluteInSolventGenFlow.post(solvent_generated, tag="solvent_generated")
@SoluteInSolventGenFlow.operation_hooks.on_success(store_workflow)
@SoluteInSolventGenFlow.operation
def generate_solvent(job):
    """
    Imports solvent generation results into the current job's document.

    This function is a key component of the solute in solvent generation workflow. It is responsible for
    importing the solvent's GRO and TOP files from a solvent generation project into the current job's
    document. This ensures that the solvent configuration is ready for the subsequent solvation and
    minimization steps.

    The function utilizes the `import_job_from_other_flow` utility to facilitate the copying of the
    necessary files from the solvent generation project to the current project. It updates the job
    document with the paths to these files and sets the operation's name in the project's operation
    to workflow mapping, aiding in the tracking and management of the workflow's progress.

    Args:
        job (Job): The job object associated with the current workflow operation. It contains the state point
                   and document for the job, which will be updated with the solvent generation results.

    Returns:
        None: This function does not return a value but updates the job document with the paths to the solvent's
              GRO and TOP files, marking the solvent generation step as complete.
    """
    solvent_job: Job = get_solvent_job(job)
    keys_for_files_to_copy = ["solvent_gro", "solvent_top"]
    project.operation_to_workflow[func_name()] = solvent_gen_name
    import_job_from_other_flow(
        job, solvent_gen_project, solvent_job, keys_for_files_to_copy
    )
    solvent_job: Job = get_solvent_job(job)
    keys_for_files_to_copy = ["solvent_gro", "solvent_top"]
    project.operation_to_workflow[func_name()] = solvent_gen_name
    import_job_from_other_flow(
        job, solvent_gen_project, solvent_job, keys_for_files_to_copy
    )


@SoluteInSolventGenFlow.pre(solute_generated, tag="solute_generated")
@SoluteInSolventGenFlow.pre(solvent_generated, tag="solvent_generated")
@SoluteInSolventGenFlow.post(system_generated, tag="system_generated")
@SoluteInSolventGenFlow.operation_hooks.on_success(store_task)
@SoluteInSolventGenFlow.operation(cmd=True, with_job=True)
def solvate(job):
    """
    Solvates the solute with the solvent in preparation for molecular dynamics simulation.

    This operation combines the solute and solvent geometries into a single system, creating a solvated
    system ready for energy minimization and further simulation steps. It uses the GROMACS tool to solvate
    the solute with the solvent, generating a new topology file and coordinate file for the combined system.

    The function updates the job document with the paths to the generated files, marking the solvation step
    as complete. This is a crucial step in the workflow, as it prepares the molecular system for the
    simulation environment.

    Args:
        job (Job): The job object associated with the current workflow operation. It contains the state point
                   and document for the job, which will be updated with the paths to the generated files for
                   the solvated system.

    Returns:
        str: The command to execute the solvation process using GROMACS, with parameters configured for the
             solute and solvent files specified in the job document.
    """
    return solvate_solute_command(
        gro_solute=job.doc[solute_gen_name]["solute_gro"],
        gro_solvent=job.doc[solvent_gen_name]["solvent_gro"],
        top_solute=job.doc[solute_gen_name]["solute_top"],
        top_output=SoluteInSolventGenFlow.get_state_name(extension="top"),
        output_name=SoluteInSolventGenFlow.get_state_name("generate"),
    )


@SoluteInSolventGenFlow.pre(system_generated, tag="system_generated")
@SoluteInSolventGenFlow.post(system_minimized, tag="system_minimized")
@SoluteInSolventGenFlow.operation_hooks.on_success(store_gromacs_log_to_doc)
@SoluteInSolventGenFlow.operation(cmd=True, with_job=True)
def minimize(job):
    """
    Minimizes the energy of the solvated system to prepare it for molecular dynamics simulation.

    This operation is essential for stabilizing the molecular system before proceeding to the equilibration
    and production phases of the simulation. It combines the topology files of the solute and solvent, and
    then uses GROMACS to minimize the energy of the combined system. The process ensures that any steric clashes
    or unrealistic bond lengths/angles are corrected before the system is subjected to dynamic simulation conditions.

    The function updates the job document with the path to the new combined topology file and executes the
    GROMACS energy minimization command. The output files, including the minimized coordinate file, are essential
    for the next steps in the simulation workflow.

    Args:
        job (Job): The job object associated with the current workflow operation. It contains the state point
                   and document for the job, which will be updated with the paths to the minimized system files.

    Returns:
        str: The command to execute the energy minimization process using GROMACS, with parameters configured
             for the combined solute and solvent system specified in the job document.
    """
    solute_top = Topology.parse_top_file(job.doc[solute_gen_name]["solute_top"])
    solvent_top = Topology.parse_top_file(job.doc[solvent_gen_name]["solvent_top"])
    solute_solvent_top = Topology.parse_top_file(
        SoluteInSolventGenFlow.get_state_name(extension="top")
    )
    new_top = append_all_includes_to_top(solute_solvent_top, [solvent_top, solute_top])
    new_top.system = f"{job.sp['solute_name']} in {job.sp['solvent_name']}"
    logger.info(f"setting topology system name to {new_top.system}")
    new_top.output_top(SoluteInSolventGenFlow.get_state_name(extension="top"))
    job.doc[project_name]["solute_solvent_top"] = SoluteInSolventGenFlow.get_state_name(
        extension="top"
    )
    return gromacs_simulation_command(
        mdp=SoluteInSolventGenFlow.mdp_files.get("minimize"),
        top=SoluteInSolventGenFlow.get_state_name(extension="top"),
        gro=SoluteInSolventGenFlow.get_state_name("generate", "gro"),
        name=SoluteInSolventGenFlow.get_state_name("minimize"),
        n_threads=SoluteInSolventGenFlow.simulation_settings.get("n_threads"),
    )


@SoluteInSolventGenFlow.pre(system_minimized, tag="system_minimized")
@SoluteInSolventGenFlow.post(system_equilibrated, tag="system_equilibrated")
@SoluteInSolventGenFlow.operation_hooks.on_success(store_gromacs_log_to_doc)
@SoluteInSolventGenFlow.operation(cmd=True, with_job=True)
def equilibrate(job):
    """
    Equilibrates the minimized molecular system to ensure stability under simulation conditions.

    This operation performs the equilibration phase of the molecular dynamics simulation process, using GROMACS.
    It aims to bring the system to a stable state where temperature, pressure, and densities are balanced and
    representative of the intended simulation environment. This step is crucial for achieving realistic simulation
    results and for preparing the system for the production run.

    The function updates the job document with the path to the equilibrated coordinate file, marking the
    equilibration step as complete. It constructs and returns the GROMACS command for equilibration, configured
    with the appropriate molecular dynamics parameters (MDP) file, topology file, and the minimized coordinate file.

    Args:
        job (Job): The job object associated with the current workflow operation. It contains the state point
                   and document for the job, which will be updated with the path to the equilibrated system files.

    Returns:
        str: The command to execute the equilibration process using GROMACS, with parameters configured for the
             system specified in the job document.
    """
    job.doc[project_name]["solute_solvent_gro"] = SoluteInSolventGenFlow.get_state_name(
        "equilibrate", "gro"
    )
    return gromacs_simulation_command(
        mdp=SoluteInSolventGenFlow.mdp_files.get("equilibrate"),
        top=SoluteInSolventGenFlow.get_state_name(extension="top"),
        gro=SoluteInSolventGenFlow.get_state_name("minimize", "gro"),
        name=SoluteInSolventGenFlow.get_state_name("equilibrate"),
        n_threads=SoluteInSolventGenFlow.simulation_settings.get("n_threads"),
    )


@SoluteInSolventGenFlow.pre(system_equilibrated, tag="system_equilibrated")
@SoluteInSolventGenFlow.post(
    lambda job: all(
        [
            job.isfile(SoluteInSolventGenFlow.nomad_workflow),
            job.isfile(SoluteInSolventGenFlow.nomad_top_level_workflow),
        ]
    ),
    tag="generated_nomad_workflow",
)
@SoluteInSolventGenFlow.operation_hooks.on_success(flag_ready_for_upload)
@SoluteInSolventGenFlow.operation(with_job=True)
def generate_nomad_workflow(job):
    build_nomad_workflow(job, is_top_level=False)
    build_nomad_workflow(job, is_top_level=True)


@SoluteInSolventGenFlow.pre(is_ready_for_upload, tag="generated_nomad_workflow")
@SoluteInSolventGenFlow.post(uploaded_to_nomad)
@SoluteInSolventGenFlow.operation(with_job=True)
def upload_to_nomad(job: Job):
    return SoluteInSolventGenFlow.upload_to_nomad(job)


def get_solvation_job(job: Job) -> Job:
    sp = {
        "type": "solute_solvation",
        "solvent_name": job.sp.get("solvent_name"),
        "solute_name": job.sp.get("solute_name"),
    }
    return project.open_job(sp).init()


project = SoluteInSolventGenFlow.init_and_get_project()
solute_gen_project = SoluteGenFlow.init_and_get_project()
solvent_gen_project = SolventGenFlow.init_and_get_project()

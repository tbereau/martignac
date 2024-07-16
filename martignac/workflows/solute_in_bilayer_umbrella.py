import logging

import gromacs
import numpy as np
from flow import aggregator
from MDAnalysis import Universe

from martignac import config
from martignac.nomad.workflows import Job, build_nomad_workflow
from martignac.parsers.gromacs_topologies import combine_multiple_topology_files
from martignac.utils.gromacs import gromacs_simulation_command, run_gmx_wham
from martignac.utils.martini_flow_projects import (
    MartiniFlowProject,
    fetched_from_nomad,
    flag_ready_for_upload,
    import_job_from_other_flow,
    is_ready_for_upload,
    store_gromacs_log_to_doc_with_depth_from_bilayer_core,
    store_task,
    store_task_for_many_jobs,
    store_workflow,
    system_equilibrated,
    system_generated,
    system_minimized,
    uploaded_to_nomad,
)
from martignac.utils.misc import (
    calculate_average_com,
    convert_pdb_to_gro,
    func_name,
    sub_template_mdp,
    translate_gro_by_vector,
    update_nested_dict,
)
from martignac.utils.packmol import place_solute_in_solvent_with_packmol
from martignac.workflows.bilayer_generation import BilayerGenFlow, get_solvent_job
from martignac.workflows.solute_generation import SoluteGenFlow, get_solute_job

logger = logging.getLogger(__name__)

conf = config()["solute_in_bilayer"]


class SoluteInBilayerUmbrellaFlow(MartiniFlowProject):
    """
    Defines the workflow for simulating a solute within a lipid bilayer using umbrella sampling.

    This class extends the MartiniFlowProject, incorporating specific configurations and operations
    necessary for setting up, running, and analyzing molecular dynamics simulations of solutes in
    lipid bilayers. It includes methods for generating the solute and bilayer, translating the solute
    to the desired depth within the bilayer, and performing umbrella sampling simulations to calculate
    the potential of mean force (PMF) across different bilayer depths.

    Attributes:
        workspace_path (str): Path to the workspace directory for this project.
        itp_files (dict): Dictionary mapping solute and lipid component names to their respective ITP file paths.
        mdp_path (str): Path to the directory containing MDP files for GROMACS simulations.
        mdp_files (dict): Dictionary mapping simulation types (e.g., 'minimize', 'equilibrate') to their MDP file paths.
        simulation_settings (dict): General settings for the simulation, such as the number of threads.
        system_name (str): Name of the system being simulated.
        nomad_workflow (str): Name of the workflow file for NOMAD integration.
        nomad_top_level_workflow (str): Name of the top-level workflow file for NOMAD integration.
        state_names (dict): Dictionary mapping state names (e.g., 'generated', 'minimized') to their string representations.
    """

    workspace_path: str = f"{MartiniFlowProject.workspaces_path}/{conf['relative_paths']['workspaces']}"
    itp_files = {k: v.get(str) for k, v in conf["itp_files"].items()}
    mdp_path = f"{MartiniFlowProject.input_files_path}/{conf['relative_paths']['mdp_files']}"
    mdp_files = {k: v.get(str) for k, v in conf["mdp_files"].items()}
    simulation_settings = {"n_threads": conf["settings"]["n_threads"].get(int)}
    system_name = conf["output_names"]["system"].get(str)
    nomad_workflow: str = conf["output_names"]["nomad_workflow"].get(str)
    nomad_top_level_workflow: str = conf["output_names"]["nomad_top_level_workflow"].get(str)
    state_names = {k: v.get(str) for k, v in conf["output_names"]["states"].items()}

    @classmethod
    def get_system_run_name(cls, depth_state: str) -> str:
        return SoluteInBilayerUmbrellaFlow.get_state_name(state_name=f"production-{depth_state}")


project_name = SoluteInBilayerUmbrellaFlow.class_name()
solute_gen_name = SoluteGenFlow.class_name()
bilayer_gen_name = BilayerGenFlow.class_name()


wham_tpr_summary = "tpr-files.dat"
wham_xvg_summary = "pullf-files.dat"
wham_profile_xvg = "wham_profile.xvg"
wham_bstrap_xvg = "wham_bstrap.xvg"
wham_bsprof_xvg = "wham_bsprof.xvg"
wham_hist_xvg = "wham_hist.xvg"
wham_profile = "wham_profile.npy"
wham_bstrap = "wham_bstrap.npy"
wham_bsprof = "wham_bsprof.npy"
wham_hist = "wham_hist.npy"


def solute_in_bilayer_aggregator(jobs):
    unique_combinations = {(job.sp["solute_name"], tuple(tuple(d.items()) for d in job.sp["lipids"])) for job in jobs}
    for solute, lipid_composition in unique_combinations:
        yield [
            job
            for job in jobs
            if job.sp.solute_name == solute and tuple(tuple(d.items()) for d in job.sp.lipids) == lipid_composition
        ]


def job_system_keys(job) -> [dict, str]:
    return [job.sp.lipids, job.sp.solute_name]


def lowest_depth_job(jobs) -> Job:
    return min(jobs, key=lambda j: j.sp.depth_from_bilayer_core)


@SoluteInBilayerUmbrellaFlow.label
def solute_generated(job) -> bool:
    return (
        solute_gen_name in job.doc
        and job.doc[solute_gen_name].get("solute_gro")
        and job.doc[solute_gen_name].get("solute_top")
    )


@SoluteInBilayerUmbrellaFlow.label
def solute_translated(job) -> bool:
    return (
        solute_gen_name in job.doc
        and job.doc[project_name].get("solute_translated_gro")
        and job.doc[solute_gen_name].get("solute_top")
    )


@SoluteInBilayerUmbrellaFlow.label
def generated_box_pdb(job):
    return job.isfile(SoluteInBilayerUmbrellaFlow.get_state_name("generate", "pdb"))


@SoluteInBilayerUmbrellaFlow.label
def bilayer_generated(job) -> bool:
    return bilayer_gen_name in job.doc and job.doc[bilayer_gen_name].get("bilayer_gro")


@SoluteInBilayerUmbrellaFlow.label
def topology_updated(job) -> bool:
    return project_name in job.doc and job.doc[project_name].get("solute_bilayer_top")


@SoluteInBilayerUmbrellaFlow.label
def bilayer_depth_sampled(job) -> bool:
    return job.doc.get(project_name, False) and job.isfile(
        job.fn(job.doc[project_name].get("umbrella_pullx_xvg", default=""))
    )


@SoluteInBilayerUmbrellaFlow.pre(fetched_from_nomad)
@SoluteInBilayerUmbrellaFlow.post(solute_generated, tag="solute_generated")
@SoluteInBilayerUmbrellaFlow.operation_hooks.on_success(store_workflow)
@SoluteInBilayerUmbrellaFlow.operation
def generate_solute(job: Job):
    """
    Generates the solute for the molecular dynamics simulation within a lipid bilayer.

    This operation is part of the SoluteInBilayerUmbrellaFlow workflow, which simulates a solute
    within a lipid bilayer using umbrella sampling. It imports the solute's GRO, TOP, and ITP files
    from a solute generation project into the current job's document. This step is crucial for
    preparing the solute for subsequent operations, such as translation to the desired depth within
    the bilayer and solvation.

    The function utilizes the `import_job_from_other_flow` utility to facilitate the copying of the
    necessary files from the solute generation project to the current project. It updates the job
    document with the paths to these files and sets the operation's name in the project's operation
    to workflow mapping, aiding in the tracking and management of the workflow's progress.

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
    import_job_from_other_flow(job, solute_gen_project, solute_job, keys_for_files_to_copy)


@SoluteInBilayerUmbrellaFlow.pre(solute_generated, tag="solute_generated")
@SoluteInBilayerUmbrellaFlow.post(solute_translated, tag="solute_translated")
@SoluteInBilayerUmbrellaFlow.operation_hooks.on_success(store_workflow)
@SoluteInBilayerUmbrellaFlow.operation(with_job=True)
def translate_solute(job: Job):
    """
    Translates the solute to the center of mass (COM) of the bilayer.

    This operation is crucial for positioning the solute at the correct depth within the lipid bilayer,
    ensuring that the simulation reflects the intended physical scenario. It calculates the COM of the bilayer
    and the solute separately, then translates the solute's position to align its COM with that of the bilayer.
    The translation vector is determined by the difference in the COMs of the solute and the bilayer, effectively
    centering the solute within the bilayer along the z-axis. The translated coordinates are saved in a new GRO file.

    Args:
        job (Job): The job object associated with the current workflow operation. It contains the state point
                   and document for the job, which will be updated with the path to the translated solute GRO file.

    Returns:
        None: This function does not return a value but updates the job document with the path to the translated
              solute GRO file, marking the solute translation step as complete.
    """
    lipid_names = job.doc[bilayer_gen_name].get("lipid_names")
    com_bilayer = calculate_average_com(job.doc[bilayer_gen_name].get("bilayer_gro"), lipid_names)
    logger.info(f"bilayer COM: {com_bilayer}")
    job.doc[project_name]["bilayer_z_mean"] = com_bilayer[2]
    solute_gro = job.doc[solute_gen_name].get("solute_gro")
    solute_translated_gro = SoluteInBilayerUmbrellaFlow.get_state_name("solute_translated", "gro")
    job.doc[project_name]["solute_translated_gro"] = solute_translated_gro
    com_solute = calculate_average_com(solute_gro, [])
    logger.debug(f"solute COM: {com_solute}")
    diff_solute_com = -com_solute
    logger.info(f"Translating the solute COM by {diff_solute_com}")
    translate_gro_by_vector(solute_gro, solute_translated_gro, diff_solute_com)
    project.operation_to_workflow[func_name()] = project_name


@SoluteInBilayerUmbrellaFlow.pre(fetched_from_nomad)
@SoluteInBilayerUmbrellaFlow.post(bilayer_generated, tag="bilayer_generated")
@SoluteInBilayerUmbrellaFlow.operation_hooks.on_success(store_workflow)
@SoluteInBilayerUmbrellaFlow.operation
def generate_bilayer(job):
    """
    Generates the lipid bilayer structure for the molecular dynamics simulation.

    This operation is a critical part of the SoluteInBilayerUmbrellaFlow workflow, focusing on the
    simulation of a solute within a lipid bilayer using umbrella sampling. It imports the bilayer's
    GRO and TOP files from a bilayer generation project into the current job's document. This step
    is essential for preparing the bilayer structure for subsequent operations, such as solute
    insertion, solvation, and the execution of molecular dynamics simulations.

    The function leverages the `import_job_from_other_flow` utility to facilitate the transfer of
    the necessary files from the bilayer generation project to the current project. It updates the
    job document with the paths to these files and sets the operation's name in the project's
    operation to workflow mapping, aiding in the tracking and management of the workflow's progress.

    Args:
        job (Job): The job object associated with the current workflow operation. It contains the
                   state point and document for the job, which will be updated with the bilayer
                   generation results.

    Returns:
        None: This function does not return a value but updates the job document with the paths to
              the bilayer's GRO and TOP files, marking the bilayer generation step as complete.
    """
    solvent_job: Job = get_solvent_job(job)
    keys_for_files_to_copy = ["bilayer_gro", "bilayer_top"]
    project.operation_to_workflow[func_name()] = bilayer_gen_name
    import_job_from_other_flow(job, bilayer_gen_project, solvent_job, keys_for_files_to_copy)


@SoluteInBilayerUmbrellaFlow.pre(solute_translated, tag="solute_translated")
@SoluteInBilayerUmbrellaFlow.pre(bilayer_generated, tag="bilayer_generated")
@SoluteInBilayerUmbrellaFlow.post(generated_box_pdb, tag="generated_box_pdb")
@SoluteInBilayerUmbrellaFlow.operation_hooks.on_success(store_task)
@SoluteInBilayerUmbrellaFlow.operation(cmd=True, with_job=True)
def insert_solute_in_box(job):
    """
    Inserts the solute into the simulation box with the lipid bilayer.

    This operation positions the solute within the lipid bilayer at a specified depth, using the Packmol tool to
    ensure the solute is correctly placed within the solvent environment. The depth is adjusted based on the
    job's state point information, specifically the 'depth_from_bilayer_core' parameter. The resulting structure
    is saved in a PDB file, which is later converted to GRO format for GROMACS simulations.

    The function calculates the restraint position along the z-axis for the solute, ensuring it is placed at the
    correct depth within the bilayer. This is crucial for simulating the interaction of the solute with the lipid
    bilayer at specific locations.

    Args:
        job (Job): The job object associated with the current workflow operation. It contains the state point
                   and document for the job, which will be updated with the path to the generated PDB file.

    Returns:
        str: The command to execute Packmol with the specified parameters, placing the solute in the solvent
             environment at the correct depth.
    """
    return place_solute_in_solvent_with_packmol(
        gro_solute=job.doc[project_name]["solute_translated_gro"],
        gro_solvent=job.doc[bilayer_gen_name]["bilayer_gro"],
        output_pdb=SoluteInBilayerUmbrellaFlow.get_state_name("generate", "pdb"),
        restraint_along_z=job.doc[project_name].get("bilayer_z_mean") + job.sp.depth_from_bilayer_core * 10.0,
    )


@SoluteInBilayerUmbrellaFlow.pre(generated_box_pdb, tag="generated_box_pdb")
@SoluteInBilayerUmbrellaFlow.post(system_generated, tag="system_generated")
@SoluteInBilayerUmbrellaFlow.operation_hooks.on_success(store_task)
@SoluteInBilayerUmbrellaFlow.operation(with_job=True)
def convert_box_to_gro(job: Job):
    """
    Converts the simulation box from PDB to GRO format after solute insertion.

    This operation is essential for preparing the simulation system in a format compatible with GROMACS simulations.
    It uses the MDAnalysis package to load the current state of the simulation box in PDB format, extracts the box
    dimensions, and then converts it to GRO format. The GRO format is required for subsequent simulation steps,
    including minimization, equilibration, and production runs in the GROMACS molecular dynamics package.

    Args:
        job (Job): The job object associated with the current workflow operation. It contains the state point
                   and document for the job, which will be updated with the path to the generated GRO file.

    Returns:
        None: This function does not return a value but updates the job document with the path to the converted
              GRO file, marking the box conversion step as complete.
    """
    universe = Universe(job.doc[bilayer_gen_name]["bilayer_gro"])
    box_dimensions = universe.dimensions
    convert_pdb_to_gro(
        SoluteInBilayerUmbrellaFlow.get_state_name("generate", "pdb"),
        SoluteInBilayerUmbrellaFlow.get_state_name("generate", "gro"),
        box_dimensions,
    )


@SoluteInBilayerUmbrellaFlow.pre(system_generated, tag="system_generated")
@SoluteInBilayerUmbrellaFlow.post(topology_updated, tag="topology_updated")
@SoluteInBilayerUmbrellaFlow.operation_hooks.on_success(store_task)
@SoluteInBilayerUmbrellaFlow.operation(with_job=True)
def update_topology_file(job: Job):
    """
    Updates the topology file by combining the solute and bilayer topology files.

    This operation is crucial for preparing the simulation system for GROMACS simulations. It combines the topology
    files of the solute and the bilayer into a single topology file that includes both components. This combined
    topology file is necessary for running simulations that involve both the solute and the bilayer, ensuring that
    the interactions between them are accurately represented.

    The function first gathers the paths to the solute and bilayer topology files from the job document. It then
    uses the `combine_multiple_topology_files` utility to merge these files into a single topology file. The combined
    topology file is then updated to ensure that the molecule counts match the actual system composition as defined
    by the GRO file generated in previous steps. Finally, the path to the combined topology file is stored in the
    job document for use in subsequent simulation steps.

    Args:
        job (Job): The job object associated with the current workflow operation. It contains the state point
                   and document for the job, which will be updated with the path to the combined topology file.

    Returns:
        None: This function does not return a value but updates the job document with the path to the combined
              topology file, marking the topology update step as complete.
    """
    comb_topology = combine_multiple_topology_files(
        [job.doc[bilayer_gen_name]["bilayer_top"], job.doc[solute_gen_name]["solute_top"]], "solute in bilayer"
    )
    top_file_name = SoluteInBilayerUmbrellaFlow.get_state_name(extension="top")
    job.doc[project_name]["solute_bilayer_top"] = top_file_name
    comb_topology.update_counts_against_gro(SoluteInBilayerUmbrellaFlow.get_state_name("generate", "gro"))
    comb_topology.output_top(top_file_name)


@SoluteInBilayerUmbrellaFlow.pre(topology_updated, tag="topology_updated")
@SoluteInBilayerUmbrellaFlow.post(system_minimized, tag="system_minimized")
@SoluteInBilayerUmbrellaFlow.operation_hooks.on_success(store_gromacs_log_to_doc_with_depth_from_bilayer_core)
@SoluteInBilayerUmbrellaFlow.operation(cmd=True, with_job=True)
def minimize(job: Job):
    """
    Executes the minimization step for the molecular dynamics simulation.

    This function prepares and runs the energy minimization step using GROMACS, which is essential for stabilizing
    the molecular system before further simulation steps. It utilizes the minimization parameters defined in the
    MDP file specified for the minimization process. The function constructs and executes the GROMACS command for
    energy minimization, ensuring the system is in a low-energy state and free of steric clashes or inappropriate
    geometries that could lead to simulation instabilities.

    Args:
        job (Job): The job object associated with the current workflow operation. It contains the state point
                   and document for the job, which will be updated with the results of the minimization step.

    Returns:
        str: The command executed for the minimization process, primarily for logging or debugging purposes.
    """
    return gromacs_simulation_command(
        mdp=SoluteInBilayerUmbrellaFlow.mdp_files.get("minimize"),
        top=SoluteInBilayerUmbrellaFlow.get_state_name(extension="top"),
        gro=SoluteInBilayerUmbrellaFlow.get_state_name("generate", "gro"),
        name=SoluteInBilayerUmbrellaFlow.get_state_name("minimize"),
        n_threads=SoluteInBilayerUmbrellaFlow.simulation_settings.get("n_threads"),
    )


@SoluteInBilayerUmbrellaFlow.pre(system_minimized, tag="system_minimized")
@SoluteInBilayerUmbrellaFlow.post(system_equilibrated, tag="system_equilibrated")
@SoluteInBilayerUmbrellaFlow.operation_hooks.on_success(store_gromacs_log_to_doc_with_depth_from_bilayer_core)
@SoluteInBilayerUmbrellaFlow.operation(cmd=True, with_job=True)
def equilibrate(job):
    """
    Prepares and executes the equilibration step of the molecular dynamics simulation.

    This function is responsible for equilibrating the system after the minimization step, ensuring that the system
    reaches a stable state before proceeding to the production phase of the simulation. It dynamically generates an
    MDP file tailored to the specific depth state of the solute within the bilayer, incorporating the solute name and
    the first lipid name from the job's state point information into the MDP file. This customization allows for
    precise control over the equilibration conditions based on the system's configuration.

    The equilibration process uses GROMACS for simulation, with parameters specified in the dynamically generated MDP
    file. The function constructs and executes the GROMACS command for equilibration, using the topology file and the
    GRO file generated from the minimization step as inputs.

    Args:
        job (Job): The job object associated with the current workflow operation. It contains the state point
                   and document for the job, which will be updated with the results of the equilibration step.

    Returns:
        str: The command executed for the equilibration process, primarily for logging or debugging purposes.
    """
    depth_state = str(job.sp.depth_from_bilayer_core)
    depth_mdp = f"{SoluteInBilayerUmbrellaFlow.get_system_run_name(depth_state)}.mdp"
    sub_template_mdp(SoluteInBilayerUmbrellaFlow.mdp_files.get("equilibrate"), "sedstate", depth_state, depth_mdp)
    sub_template_mdp(depth_mdp, "sedmolecule", job.sp["solute_name"], depth_mdp)
    sub_template_mdp(depth_mdp, "sedlipid", job.sp["lipids"][0]["name"], depth_mdp)
    return gromacs_simulation_command(
        mdp=depth_mdp,
        top=SoluteInBilayerUmbrellaFlow.get_state_name(extension="top"),
        gro=SoluteInBilayerUmbrellaFlow.get_state_name("minimize", "gro"),
        name=SoluteInBilayerUmbrellaFlow.get_state_name("equilibrate"),
        n_threads=SoluteInBilayerUmbrellaFlow.simulation_settings.get("n_threads"),
    )


@SoluteInBilayerUmbrellaFlow.pre(system_equilibrated, tag="system_equilibrated")
@SoluteInBilayerUmbrellaFlow.post(bilayer_depth_sampled, tag="system_sampled")
@SoluteInBilayerUmbrellaFlow.operation_hooks.on_success(store_gromacs_log_to_doc_with_depth_from_bilayer_core)
@SoluteInBilayerUmbrellaFlow.operation(cmd=True, with_job=True)
def production(job):
    """
    Executes the production phase of the molecular dynamics simulation.

    This function is responsible for running the production phase of the simulation, where the actual data collection
    occurs. It dynamically generates an MDP file tailored to the specific depth state of the solute within the bilayer,
    incorporating the solute name and the first lipid name from the job's state point information into the MDP file.
    This customization allows for precise control over the simulation conditions based on the system's configuration.

    The production process uses GROMACS for simulation, with parameters specified in the dynamically generated MDP
    file. The function constructs and executes the GROMACS command for the production run, using the topology file and
    the GRO file generated from the equilibration step as inputs. Additionally, it prepares files for umbrella sampling
    analysis by specifying output filenames for pull force and pull coordinate data, as well as the TPR and log files
    for the run.

    Args:
        job (Job): The job object associated with the current workflow operation. It contains the state point
                   and document for the job, which will be updated with the results of the production step and
                   the paths to the generated output files for umbrella sampling analysis.

    Returns:
        str: The command executed for the production process, primarily for logging or debugging purposes.
    """
    depth_state = str(job.sp.depth_from_bilayer_core)
    depth_mdp = f"{SoluteInBilayerUmbrellaFlow.get_system_run_name(depth_state)}.mdp"
    sub_template_mdp(SoluteInBilayerUmbrellaFlow.mdp_files.get("production"), "sedstate", depth_state, depth_mdp)
    sub_template_mdp(depth_mdp, "sedmolecule", job.sp["solute_name"], depth_mdp)
    sub_template_mdp(depth_mdp, "sedlipid", job.sp["lipids"][0]["name"], depth_mdp)
    job.doc[project_name][
        "umbrella_pullf_xvg"
    ] = f"{SoluteInBilayerUmbrellaFlow.get_system_run_name(depth_state)}_pullf.xvg"
    job.doc[project_name][
        "umbrella_pullx_xvg"
    ] = f"{SoluteInBilayerUmbrellaFlow.get_system_run_name(depth_state)}_pullx.xvg"
    job.doc[project_name]["tpr_file"] = f"{SoluteInBilayerUmbrellaFlow.get_system_run_name(depth_state)}.tpr"
    job.doc[project_name]["umbrella_log"] = f"{SoluteInBilayerUmbrellaFlow.get_system_run_name(depth_state)}.log"
    return gromacs_simulation_command(
        mdp=depth_mdp,
        top=SoluteInBilayerUmbrellaFlow.get_state_name(extension="top"),
        gro=SoluteInBilayerUmbrellaFlow.get_state_name("equilibrate", "gro"),
        name=SoluteInBilayerUmbrellaFlow.get_system_run_name(depth_state),
        n_threads=SoluteInBilayerUmbrellaFlow.simulation_settings.get("n_threads"),
    )


@SoluteInBilayerUmbrellaFlow.label
def wham_calculated(*jobs) -> bool:
    return any([job.isfile(wham_profile_xvg) and job.isfile(wham_hist_xvg) for job in jobs])


@SoluteInBilayerUmbrellaFlow.label
def free_energy_calculated(*jobs) -> bool:
    return all(
        [
            job.doc.get(project_name, False)
            and job.doc[project_name].get("free_energy", False)
            and job.doc[project_name]["free_energy"].get("profile", False)
            for job in jobs
        ]
    )


@SoluteInBilayerUmbrellaFlow.label
def all_depth_states_sampled(*jobs):
    result = all(
        [
            job.doc.get(project_name, False)
            and job.isfile(job.fn(job.doc[project_name].get("umbrella_pullx_xvg", default="")))
            for job in jobs
        ]
    )
    return result


@SoluteInBilayerUmbrellaFlow.pre(all_depth_states_sampled, tag="system_sampled")
@SoluteInBilayerUmbrellaFlow.post(wham_calculated, tag="wham_calculated")
@SoluteInBilayerUmbrellaFlow.operation_hooks.on_success(store_task_for_many_jobs)
@SoluteInBilayerUmbrellaFlow.operation(
    aggregator=aggregator(aggregator_function=solute_in_bilayer_aggregator, sort_by="depth_from_bilayer_core"),
    cmd=True,
)
def compute_wham(*jobs):
    """
    Calculates the Weighted Histogram Analysis Method (WHAM) for a set of jobs.

    This function is responsible for computing the potential of mean force (PMF) across different bilayer depths
    using the WHAM. It aggregates the pull force and TPR files from all jobs, generates summary files for these,
    and then runs the GROMACS WHAM tool. The WHAM analysis is crucial for understanding the free energy landscape
    of the solute across the bilayer, providing insights into its stability and behavior at various depths.

    The function first identifies the job with the lowest depth (as a representative job) to perform the WHAM
    calculation. It then generates summary files listing all the pull force and TPR files from the jobs. These
    summary files are used as input for the GROMACS WHAM tool, which calculates the PMF and generates several
    output files, including the PMF profile and histograms.

    Args:
        *jobs: A variable number of Job objects, each representing a simulation at a specific depth within the
               bilayer. These jobs should have completed the production phase of the simulation and contain the
               necessary files for WHAM analysis.

    Returns:
        str: The command executed for the WHAM analysis, primarily for logging or debugging purposes. This is
             the command that invokes the GROMACS WHAM tool with the appropriate parameters and input files.
    """
    logger.info("calculation of potential of mean force")
    xvg_files = [job.fn(job.doc[project_name]["umbrella_pullf_xvg"]) for job in jobs]
    tpr_files = [job.fn(job.doc[project_name]["tpr_file"]) for job in jobs]

    def generate_summary_file(output_filename: str, list_of_files: list[str]) -> None:
        with open(output_filename, "w") as file:
            for fil in list_of_files:
                file.write(f"{fil}\n")

    for job in jobs:
        if job == lowest_depth_job(jobs):
            with job:
                generate_summary_file(wham_tpr_summary, tpr_files)
                generate_summary_file(wham_xvg_summary, xvg_files)
                return run_gmx_wham(
                    job.fn(wham_tpr_summary),
                    job.fn(wham_xvg_summary),
                    job.fn(wham_profile_xvg),
                    job.fn(wham_hist_xvg),
                    job.fn(wham_bstrap_xvg),
                    job.fn(wham_bsprof_xvg),
                )


@SoluteInBilayerUmbrellaFlow.pre(wham_calculated, tag="wham_calculated")
@SoluteInBilayerUmbrellaFlow.post(free_energy_calculated, tag="free_energy_calculated")
@SoluteInBilayerUmbrellaFlow.operation_hooks.on_success(store_task_for_many_jobs)
@SoluteInBilayerUmbrellaFlow.operation(
    aggregator=aggregator(aggregator_function=solute_in_bilayer_aggregator, sort_by="depth_from_bilayer_core")
)
def analyze_wham(*jobs):
    """
    Analyzes the WHAM (Weighted Histogram Analysis Method) results to calculate free energy profiles.

    This function processes the output from the WHAM analysis, converting the XVG files generated by GROMACS WHAM tool
    into numpy arrays for easier analysis and visualization. It focuses on the lowest depth job as a representative
    for the analysis, assuming that the WHAM analysis has been completed and the relevant XVG files are available.

    The function reads the WHAM profile, histogram, and bootstrap files, converting each into a numpy array. These
    arrays are then saved to the job's filesystem for future use. Additionally, the job document is updated with
    the paths to these numpy arrays, facilitating access to the analysis results.

    Args:
        *jobs: A variable number of Job objects, each representing a simulation at a specific depth within the
               bilayer. These jobs should have completed the WHAM analysis phase and contain the necessary files
               for free energy calculation.

    Returns:
        None: This function does not return a value but updates the job document and saves numpy arrays of the
              WHAM analysis results to the job's filesystem.
    """
    for job in jobs:
        if job == lowest_depth_job(jobs):
            xvg_profile = gromacs.fileformats.XVG(job.fn(wham_profile_xvg))
            np.save(job.fn(wham_profile), xvg_profile.array)
            xvg_hist = gromacs.fileformats.XVG(job.fn(wham_hist_xvg))
            np.save(job.fn(wham_hist), xvg_hist.array)
            xvg_bstrap = gromacs.fileformats.XVG(job.fn(wham_bstrap_xvg))
            np.save(job.fn(wham_bstrap), xvg_bstrap.array)
        job.doc = update_nested_dict(
            job.doc,
            {
                project_name: {
                    "free_energy": {
                        "job_id": job.id,
                        "profile": wham_profile,
                        "hist": wham_hist,
                        "bootstrap": wham_bstrap,
                    }
                }
            },
        )


@SoluteInBilayerUmbrellaFlow.pre(free_energy_calculated, tag="free_energy_calculated")
@SoluteInBilayerUmbrellaFlow.post(
    lambda job: all(
        [
            job.isfile(SoluteInBilayerUmbrellaFlow.nomad_workflow),
            job.isfile(SoluteInBilayerUmbrellaFlow.nomad_top_level_workflow),
        ]
    ),
    tag="generated_nomad_workflow",
)
@SoluteInBilayerUmbrellaFlow.operation_hooks.on_success(flag_ready_for_upload)
@SoluteInBilayerUmbrellaFlow.operation(with_job=True)
def generate_nomad_workflow(job):
    build_nomad_workflow(job, is_top_level=False)
    build_nomad_workflow(job, is_top_level=True)


@SoluteInBilayerUmbrellaFlow.pre(is_ready_for_upload, tag="generated_nomad_workflow")
@SoluteInBilayerUmbrellaFlow.post(uploaded_to_nomad)
@SoluteInBilayerUmbrellaFlow.operation(with_job=True)
def upload_to_nomad(job: Job):
    return SoluteInBilayerUmbrellaFlow.upload_to_nomad(job)


def get_solvation_job(job: Job) -> Job:
    sp = {
        "type": "solute_solvation",
        "lipids": job.sp.get("lipids"),
        "solute_name": job.sp.get("solute_name"),
        "depth_from_bilayer_core": job.sp.get("depth_from_bilayer_core"),
    }
    return project.open_job(sp).init()


project = SoluteInBilayerUmbrellaFlow.get_project(path=SoluteInBilayerUmbrellaFlow.workspace_path)
solute_gen_project = SoluteGenFlow.get_project(path=SoluteGenFlow.workspace_path)
bilayer_gen_project = BilayerGenFlow.get_project(path=BilayerGenFlow.workspace_path)

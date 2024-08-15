import logging

import pandas as pd
from alchemlyb.estimators import MBAR
from alchemlyb.parsing.gmx import extract_u_nk
from flow import aggregator

from martignac import config
from martignac.nomad.workflows import (
    build_nomad_workflow,
    build_nomad_workflow_with_multiple_jobs,
)
from martignac.utils.gromacs import generate_lambdas, gromacs_simulation_command
from martignac.utils.martini_flow_projects import (
    Job,
    MartiniFlowProject,
    fetched_from_nomad,
    flag_ready_for_upload_multiple_jobs,
    import_job_from_other_flow,
    store_gromacs_log_to_doc_with_state_point,
    store_task_for_many_jobs,
    store_workflow_for_many_jobs,
    uploaded_to_nomad,
)
from martignac.utils.misc import func_name, sub_template_mdp, update_nested_dict
from martignac.workflows.solute_generation import SoluteGenFlow, get_solute_job
from martignac.workflows.solute_in_solvent_generation import (
    SoluteInSolventGenFlow,
    get_solvation_job,
)
from martignac.workflows.solvent_generation import SolventGenFlow, get_solvent_job

logger = logging.getLogger(__name__)

conf = config()["solute_in_solvent_alchemical"]


class SoluteInSolventAlchemicalFlow(MartiniFlowProject):
    """
    Defines the workflow for alchemical transformations of solutes in solvents within the Martini force field framework.

    This class extends the MartiniFlowProject, incorporating specific configurations and operations for conducting
    alchemical transformations. Alchemical transformations involve gradually changing the Hamiltonian of the system
    to transform a solute into another or modify its interactions with the solvent, which is a common technique in
    free energy calculations.

    Attributes:
        workspace_path (str): The path to the workspace directory where the project's files are stored.
        mdp_path (str): The path to the directory containing the MDP (Molecular Dynamics Parameters) files.
        itp_files (dict): A dictionary mapping the names of ITP (Include Topology) files to their paths.
        mdp_files (dict): A dictionary mapping the names of MDP files to their paths.
        simulation_settings (dict): A dictionary containing simulation settings such as the number of threads and temperature.
        system_name (str): The name of the system being simulated.
        nomad_workflow (str): The name of the workflow file for integration with the NOMAD framework.
        nomad_top_level_workflow (str): The name of the top-level workflow file for NOMAD.
        state_names (dict): A dictionary mapping state names to their string representations.

    The class provides methods for preparing the system, running production simulations, computing free energies,
    and generating and uploading workflows to NOMAD, facilitating an automated and scalable approach to molecular
    simulations.
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
        "temperature": conf["settings"]["temperature"].get(float),
    }
    system_name = conf["output_names"]["system"].get(str)
    nomad_workflow: str = conf["output_names"]["nomad_workflow"].get(str)
    nomad_top_level_workflow: str = conf["output_names"][
        "nomad_top_level_workflow"
    ].get(str)
    state_names = {k: v.get(str) for k, v in conf["output_names"]["states"].items()}

    @classmethod
    def get_system_run_name(cls, lambda_state: str) -> str:
        return SoluteInSolventAlchemicalFlow.get_state_name(
            state_name=f"production-{lambda_state}"
        )


project_name = SoluteInSolventAlchemicalFlow.class_name()


def solvent_and_solute_aggregator(jobs):
    unique_combinations = {
        (job.sp["solvent_name"], job.sp["solute_name"]) for job in jobs
    }
    for solvent, solute in unique_combinations:
        yield [
            job
            for job in jobs
            if job.sp.solvent_name == solvent and job.sp.solute_name == solute
        ]


def get_num_lambda_points(jobs):
    solute_name, solvent_name = jobs[0].sp["solute_name"], jobs[0].sp["solvent_name"]
    num_lambda_points = 0
    for job in jobs:
        if (
            job.sp["solute_name"] == solute_name
            or job.sp["solvent_name"] == solvent_name
        ):
            num_lambda_points += 1
    return num_lambda_points


def job_system_keys(job) -> [str, str]:
    return [job.sp.solvent_name, job.sp.solute_name]


def lowest_lambda_job(jobs) -> Job:
    return min(jobs, key=lambda j: j.sp.lambda_state)


@SoluteInSolventAlchemicalFlow.label
def system_prepared(job) -> bool:
    if project_name in job.document:
        return job.doc[project_name].get("system_prepared", False)
    return False


@SoluteInSolventAlchemicalFlow.label
def all_jobs_prepared(*jobs) -> bool:
    return all(system_prepared(j) for j in jobs)


@SoluteInSolventAlchemicalFlow.label
def nomad_workflow_built(*jobs) -> bool:
    return all(
        [
            any(
                job.isfile(SoluteInSolventAlchemicalFlow.nomad_workflow) for job in jobs
            ),
            any(
                job.isfile(SoluteInSolventAlchemicalFlow.nomad_top_level_workflow)
                for job in jobs
            ),
        ]
    )


@SoluteInSolventAlchemicalFlow.label
def any_uploaded_to_nomad(*jobs) -> bool:
    return any(uploaded_to_nomad(job) for job in jobs)


@SoluteInSolventAlchemicalFlow.pre(fetched_from_nomad)
@SoluteInSolventAlchemicalFlow.post(all_jobs_prepared, tag="system_prepared")
@SoluteInSolventAlchemicalFlow.operation_hooks.on_success(store_workflow_for_many_jobs)
@SoluteInSolventAlchemicalFlow.operation(
    aggregator=aggregator(aggregator_function=solvent_and_solute_aggregator)
)
def prepare_system(*jobs):
    """
    Prepares the system for alchemical transformation simulations by importing necessary data from related projects.

    This function is a critical part of the workflow in the SoluteInSolventAlchemicalFlow project. It ensures that
    all necessary input files and configurations are correctly set up for each job involved in the alchemical
    transformation process. The function performs several key operations:

    1. Identifies the job with the lowest lambda state to serve as the primary job for data import.
    2. Imports solute-related input files from the SoluteGenFlow project into the primary job.
    3. Imports solvent-related input files from the SolventGenFlow project into the primary job.
    4. Imports solvation-related input files from the SoluteInSolventGenFlow project into the primary job.
    5. Replicates the imported data across all other jobs involved in the transformation process.
    6. Updates the job document to mark the system as prepared for each job.

    The function leverages the `import_job_from_other_flow` utility to facilitate the import of data between
    different projects, ensuring that all necessary files are available for the simulation to proceed.

    Args:
        *jobs (Job): A variable number of job objects, each representing a distinct simulation job within the
                     alchemical transformation workflow.

    """
    job = lowest_lambda_job(jobs)
    logger.info(f"Preparing system for {job.id} @ AlchemicalTransformationFlow")
    solute_job: Job = get_solute_job(job)
    solute_keys = ["solute_itp"]
    import_job_from_other_flow(job, solute_gen_project, solute_job, solute_keys)
    solvent_job: Job = get_solvent_job(job)
    import_job_from_other_flow(job, solvent_gen_project, solvent_job, [])
    solvation_job: Job = get_solvation_job(job)
    solvation_keys = ["solute_solvent_gro", "solute_solvent_top"]
    project.operation_to_workflow[func_name()] = solute_solvation_project.class_name()
    import_job_from_other_flow(
        job, solute_solvation_project, solvation_job, solvation_keys
    )
    num_lambda_points = get_num_lambda_points(jobs)

    for other_job in [j for j in jobs if j != job]:
        logger.info(
            f"Importing data to job {other_job.id} @ AlchemicalTransformationFlow"
        )
        import_job_from_other_flow(
            other_job, solute_gen_project, solute_job, solute_keys, run_child_job=False
        )
        import_job_from_other_flow(
            other_job, solvent_gen_project, solvent_job, [], run_child_job=False
        )
        import_job_from_other_flow(
            other_job,
            solute_solvation_project,
            solvation_job,
            solvation_keys,
            run_child_job=False,
        )

    for job in jobs:
        job.doc = update_nested_dict(job.doc, {project_name: {"system_prepared": True}})
        job.doc[project_name]["num_lambda_points"] = num_lambda_points


@SoluteInSolventAlchemicalFlow.label
def lambda_sampled(job):
    if project_name not in job.doc:
        return False
    log_file = job.doc[project_name].get("alchemical_log", None)
    if log_file and job.isfile(log_file):
        with open(job.fn(log_file)) as file:
            lines = file.read()
        return "Finished mdrun on rank" in lines
    return False


@SoluteInSolventAlchemicalFlow.label
def all_lambda_states_sampled(*jobs):
    result = all(
        job.doc.get(project_name, False)
        and job.isfile(job.fn(job.doc[project_name].get("alchemical_xvg", default="")))
        for job in jobs
    )
    return result


@SoluteInSolventAlchemicalFlow.label
def free_energy_already_calculated(*jobs):
    return all(
        job.doc.get(project_name, False)
        and job.doc[project_name].get("free_energy", False)
        and job.doc[project_name]["free_energy"].get("mean", False)
        for job in jobs
    )


@SoluteInSolventAlchemicalFlow.pre(system_prepared, tag="system_prepared")
@SoluteInSolventAlchemicalFlow.post(lambda_sampled, tag="lambda_sampled")
@SoluteInSolventAlchemicalFlow.operation_hooks.on_success(
    store_gromacs_log_to_doc_with_state_point
)
@SoluteInSolventAlchemicalFlow.operation(cmd=True, with_job=True)
def production(job):
    """
    Executes the production phase of the molecular dynamics simulation for a given job.

    This function configures and runs the production phase of the alchemical transformation simulation. It involves
    setting up the Molecular Dynamics Parameters (MDP) file specific to the lambda state of the job, customizing it
    for the solute involved, and executing the GROMACS simulation. The function generates and stores the paths to
    the output XVG and log files in the job document for later analysis.

    The MDP file is first templated with the lambda state and solute name, then used to run the simulation with the
    specified number of threads. This phase is critical for generating the data necessary for free energy calculation.

    Args:
        job (Job): The job object representing a single simulation job within the alchemical transformation workflow.
                   It contains the state point information, including the lambda state and solute name, and provides
                   context for the simulation, including paths to input and output files.

    Returns:
        str: The command executed for the production phase, primarily for logging or debugging purposes. This includes
             the path to the modified MDP file, the output XVG file for the lambda state, and the log file.
    """
    solute_has_charged_beads = job.doc[solute_gen_project.class_name()].get(
        "solute_has_charged_beads"
    )
    num_lambda_points = job.doc[project_name]["num_lambda_points"]
    vdw_lambdas, coul_lambdas = generate_lambdas(
        num_lambda_points, turn_off_coulomb=not solute_has_charged_beads
    )
    vdw_lambda_str = " ".join(map(str, vdw_lambdas))
    coul_lambda_str = " ".join(map(str, coul_lambdas))
    lambda_state = str(job.sp.lambda_state)
    lambda_mdp = (
        f"{SoluteInSolventAlchemicalFlow.get_system_run_name(job.sp.lambda_state)}.mdp"
    )
    sub_template_mdp(
        SoluteInSolventAlchemicalFlow.mdp_files.get("production"),
        "sedstate",
        lambda_state,
        lambda_mdp,
    )
    sub_template_mdp(lambda_mdp, "sedmolecule", job.sp["solute_name"], lambda_mdp)
    sub_template_mdp(lambda_mdp, "sedvdwlambdas", vdw_lambda_str, lambda_mdp)
    sub_template_mdp(lambda_mdp, "sedcoullambdas", coul_lambda_str, lambda_mdp)
    job.doc[project_name][
        "alchemical_xvg"
    ] = f"{SoluteInSolventAlchemicalFlow.get_system_run_name(lambda_state)}.xvg"
    job.doc[project_name][
        "alchemical_log"
    ] = f"{SoluteInSolventAlchemicalFlow.get_system_run_name(lambda_state)}.log"
    return gromacs_simulation_command(
        mdp=lambda_mdp,
        top=job.doc[solute_solvation_project.class_name()].get("solute_solvent_top"),
        gro=job.doc[solute_solvation_project.class_name()].get("solute_solvent_gro"),
        name=SoluteInSolventAlchemicalFlow.get_system_run_name(lambda_state),
        n_threads=SoluteInSolventAlchemicalFlow.simulation_settings.get("n_threads"),
    )


@SoluteInSolventAlchemicalFlow.pre(all_lambda_states_sampled, tag="lambda_sampled")
@SoluteInSolventAlchemicalFlow.post(free_energy_already_calculated, tag="mbar")
@SoluteInSolventAlchemicalFlow.operation_hooks.on_success(store_task_for_many_jobs)
@SoluteInSolventAlchemicalFlow.operation(
    aggregator=aggregator(
        aggregator_function=solvent_and_solute_aggregator, sort_by="lambda_state"
    )
)
def compute_free_energy(*jobs):
    """
    Computes the free energy difference for alchemical transformations using the MBAR method.

    This function aggregates the simulation data from multiple jobs, each representing a different lambda state in
    the alchemical transformation process. It uses the alchemlyb library's MBAR estimator to calculate the free energy
    difference between the initial and final states of the transformation. The results, including the mean free energy
    difference and its standard deviation, are stored in the job document for each job.

    The function assumes that the necessary simulation data (XVG files) are already generated and stored by the
    production phase of the workflow. It leverages the extract_u_nk function from alchemlyb to parse these files and
    extract the reduced potential energies, which are then combined and passed to the MBAR estimator.

    Args:
        *jobs (Job): A variable number of job objects, each representing a distinct lambda state in the alchemical
                     transformation process. These jobs should have already completed the production phase, with XVG
                     files available for analysis.

    Note:
        This function is part of the SoluteInSolventAlchemicalFlow project and is designed to be used within a
        signac-flow workflow. It relies on specific project and job structure defined in the MartiniFlowProject class
        and its derivatives.
    """
    logger.info("calculating free energies using pyMBAR")
    xvg_files = [job.fn(job.doc[project_name]["alchemical_xvg"]) for job in jobs]
    u_nk_list = [
        extract_u_nk(
            f, T=SoluteInSolventAlchemicalFlow.simulation_settings.get("temperature")
        )
        for f in xvg_files
    ]
    u_nk_combined = pd.concat(u_nk_list)
    mbar = MBAR().fit(u_nk_combined)
    free_energy = float(mbar.delta_f_.iloc[0, -1])
    d_free_energy = float(mbar.d_delta_f_.iloc[0, -1])
    solvent_name, solute_name = job_system_keys(jobs[0])
    logger.info(
        f"free energy {solvent_name}_{solute_name}: {free_energy:+6.2f} +/- {d_free_energy:5.2f} kT"
    )
    for job in jobs:
        job.doc = update_nested_dict(
            job.doc,
            {
                project_name: {
                    "free_energy": {"mean": free_energy, "std": d_free_energy}
                }
            },
        )


@SoluteInSolventAlchemicalFlow.pre(free_energy_already_calculated, tag="mbar")
@SoluteInSolventAlchemicalFlow.post(nomad_workflow_built)
@SoluteInSolventAlchemicalFlow.operation_hooks.on_success(
    flag_ready_for_upload_multiple_jobs
)
@SoluteInSolventAlchemicalFlow.operation(
    aggregator=aggregator(aggregator_function=solvent_and_solute_aggregator)
)
def generate_nomad_workflow(*jobs):
    for job in jobs:
        with job:
            build_nomad_workflow(job, is_top_level=False, add_job_id=True)
    job = lowest_lambda_job(jobs)
    with job:
        build_nomad_workflow_with_multiple_jobs(project, list(jobs))


@SoluteInSolventAlchemicalFlow.pre(nomad_workflow_built)
@SoluteInSolventAlchemicalFlow.post(any_uploaded_to_nomad)
@SoluteInSolventAlchemicalFlow.operation(
    aggregator=aggregator(aggregator_function=solvent_and_solute_aggregator)
)
def upload_to_nomad(*jobs):
    return SoluteInSolventAlchemicalFlow.upload_to_nomad_multiple_jobs(list(jobs))


project = SoluteInSolventAlchemicalFlow.init_and_get_project()
solute_gen_project = SoluteGenFlow.init_and_get_project()
solvent_gen_project = SolventGenFlow.init_and_get_project()
solute_solvation_project = SoluteInSolventGenFlow.init_and_get_project()

import json
import logging
import os
import shutil
from hashlib import md5
from time import sleep
from typing import Any, Dict, TypeVar, cast

import signac
from flow import FlowProject
from signac.job import Job

from martignac import config
from martignac.nomad.datasets import get_dataset_by_id
from martignac.nomad.uploads import (
    get_upload_by_id,
    upload_files_to_nomad,
)
from martignac.utils.misc import update_nested_dict, zip_directories
from martignac.utils.nomad import generate_user_metadata

logger = logging.getLogger(__name__)

MartiniTypeFlow = TypeVar("MartiniTypeFlow", bound="MartiniFlowProject")
UPLOAD_TO_NOMAD: bool = config()["nomad"]["upload_to_nomad"].get(bool)
PUBLISH_TO_NOMAD: bool = config()["nomad"]["publish_uploads"].get(bool)


class MartiniFlowProject(FlowProject):
    """
    A specialized FlowProject for managing and executing Martini-based molecular dynamics simulations.

    This class extends the FlowProject class from the signac framework, providing functionalities specific to
    the setup, execution, and analysis of simulations using the Martini force field. It includes methods for
    registering simulation input files, uploading results to the NOMAD repository, and managing the simulation
    workflow states.

    Attributes:
        workspaces_path (str): Absolute path to the directory where simulation workspaces are stored.
        input_files_path (str): Absolute path to the directory containing input files for simulations.
        itp_path (str): Absolute path to the directory containing ITP (GROMACS topology) files.
        nomad_use_prod_database (bool): Flag indicating whether to use the production NOMAD database.
        nomad_dataset_id (str): Identifier for the dataset within the NOMAD repository.
        nomad_coauthors (list[str]): List of co-authors to be included in the NOMAD metadata.
        workspace_path (str): Path to the current job's workspace directory.
        itp_files (Dict[str, str]): Dictionary mapping ITP file names to their paths.
        mdp_path (str): Path to the directory containing MDP (GROMACS parameters) files.
        mdp_files (Dict[str, str]): Dictionary mapping MDP file names to their paths.
        simulation_settings (Dict[str, Any]): Dictionary containing settings for the simulation.
        system_name (str): Name of the simulation system.
        output_names (Dict[str, str]): Dictionary mapping output file types to their names.
        nomad_workflow (str): Identifier for the workflow within the NOMAD repository.
        nomad_top_level_workflow (str): Identifier for the top-level workflow within NOMAD.
        ff_parameters (Dict[str, Any]): Dictionary containing force field parameters.
        operation_to_workflow (Dict[str, str]): Mapping of operation names to workflow identifiers.
    """

    workspaces_path: str = config()["local"]["workspaces"]["absolute_path"].get(str)
    input_files_path: str = config()["local"]["input_files"]["absolute_path"].get(str)
    itp_path: str = config()["local"]["input_files"]["itp_files"].get(str)
    nomad_use_prod_database: bool = config()["nomad"]["use_prod"].get(bool)
    nomad_dataset_id: str = config()["nomad"]["dataset"]["id"].get(str)
    try:
        nomad_coauthors: list[str] = [
            c.get(str) for c in config()["nomad"]["coauthors"]
        ]
    except KeyError:
        nomad_coauthors: list[str] = []
    allow_symlinks: bool = config()["local"]["allow_symlinks"].get(bool)
    workspace_path: str = ""
    itp_files: Dict[str, str] = {}
    mdp_path: str
    mdp_files: Dict[str, str] = {}
    simulation_settings: Dict[str, Any] = {}
    system_name: str
    output_names: Dict[str, str] = {}
    nomad_workflow: str = ""
    nomad_top_level_workflow: str = ""
    ff_parameters: Dict[str, Any] = {}
    operation_to_workflow: Dict[str, str] = {}

    @classmethod
    def class_name(cls) -> str:
        return cls.__name__

    @classmethod
    def get_hash_for_files(cls, job: Job, files: list[str]) -> str:
        readouts = []
        for filename in files:
            with open(job.fn(filename)) as fp:
                readouts.append(fp.read())
        return md5(" ".join(readouts).encode("utf-8")).hexdigest()

    @classmethod
    def register_mdp_files(cls, job: Job) -> None:
        job.doc = update_nested_dict(
            job.doc,
            {
                cls.class_name(): {
                    "mdp_files": cls.get_hash_for_files(
                        job, list(cls.mdp_files.values())
                    )
                }
            },
        )

    @classmethod
    def register_itp_files(cls, job: Job) -> None:
        job.doc = update_nested_dict(
            job.doc,
            {
                cls.class_name(): {
                    "itp_files": cls.get_hash_for_files(
                        job, list(cls.itp_files.values())
                    )
                }
            },
        )

    @classmethod
    def get_state_name(cls, state_name: str = "", extension: str = "") -> str:
        state_ext = f"_{state_name}" if state_name else ""
        file_ext = f".{extension}" if extension else ""
        return f"{cls.system_name}{state_ext}{file_ext}"

    @classmethod
    def nomad_comment(cls, job: Job) -> dict:
        return {
            "job_id": job.id,
            "workflow_name": cls.class_name(),
            "state_point": job.sp(),
            "mdp_files": job.doc[cls.class_name()]["mdp_files"],
            "itp_files": job.doc[cls.class_name()]["itp_files"],
        }

    @classmethod
    def upload_to_nomad(
        cls,
        job: Job,
        nomad_upload_flag: bool = UPLOAD_TO_NOMAD,
        publish_flag: bool = PUBLISH_TO_NOMAD,
    ) -> None:
        if not nomad_upload_flag:
            logger.info("NOMAD upload turned off")
            return None
        if uploaded_to_nomad(job):
            project_name = cast("MartiniFlowProject", job.project).class_name()
            upload_id = job.doc[project_name].get("nomad_upload_id", "")
            logger.info(f"Workflow {cls.class_name()} already uploaded to {upload_id}")
            return None
        dataset = get_dataset_by_id(
            cls.nomad_dataset_id, use_prod=cls.nomad_use_prod_database
        )
        logger.info(f"made connection to dataset {dataset.dataset_id}")
        generate_user_metadata(
            file_name=job.fn("nomad.yaml"),
            comment=json.dumps(cls.nomad_comment(job)),
            coauthors=cls.nomad_coauthors,
            datasets=[dataset.dataset_id],
        )
        job.doc = update_nested_dict(job.doc, {"nomad_dataset_id": dataset.dataset_id})
        zip_file = zip_directories([job.path], job.id)
        upload_id = upload_files_to_nomad(zip_file, cls.nomad_use_prod_database)
        job.doc = update_nested_dict(
            job.doc, {cls.class_name(): {"nomad_upload_id": upload_id}}
        )
        if publish_flag:
            sleep(1)
            nomad_upload = get_upload_by_id(
                upload_id=upload_id, use_prod=cls.nomad_use_prod_database
            )
            nomad_upload.safe_publish()
        os.remove(zip_file)
        return None

    @classmethod
    def upload_to_nomad_multiple_jobs(
        cls,
        jobs: list[Job],
        nomad_upload_flag: bool = UPLOAD_TO_NOMAD,
        publish_flag: bool = PUBLISH_TO_NOMAD,
    ) -> None:
        if not nomad_upload_flag:
            logger.info("NOMAD upload turned off")
            return None
        dataset = get_dataset_by_id(
            cls.nomad_dataset_id, use_prod=cls.nomad_use_prod_database
        )
        logger.info(f"made connection to dataset {dataset.dataset_id}")
        for job in jobs:
            if uploaded_to_nomad(job):
                project_name = cast("MartiniFlowProject", job.project).class_name()
                upload_id = job.doc[project_name].get("nomad_upload_id", "")
                logger.info(
                    f"Workflow {cls.class_name()} already uploaded to {upload_id}"
                )
                return None
            job.doc = update_nested_dict(
                job.doc, {"nomad_dataset_id": dataset.dataset_id}
            )
            generate_user_metadata(
                file_name=job.fn("nomad.yaml"),
                comment=json.dumps(cls.nomad_comment(job)),
                coauthors=cls.nomad_coauthors,
                datasets=[dataset.dataset_id],
            )
        zip_file = zip_directories([job.path for job in jobs], jobs[0].id)
        upload_id = upload_files_to_nomad(zip_file, cls.nomad_use_prod_database)
        os.remove(zip_file)
        for job in jobs:
            job.doc = update_nested_dict(
                job.doc, {cls.class_name(): {"nomad_upload_id": upload_id}}
            )
        if publish_flag:
            sleep(1)
            nomad_upload = get_upload_by_id(
                upload_id=upload_id, use_prod=cls.nomad_use_prod_database
            )
            nomad_upload.safe_publish()
        return None

    @classmethod
    def unlink_itp_and_mdp_files(cls, job: Job) -> None:
        project = cast("MartiniFlowProject", job.project)
        logger.info(f"removing symbolic links for {job.id} @ {project.class_name()}")
        for root, _, files in os.walk(job.path):
            for file in filter(
                lambda f: f.endswith((".itp", ".mdp"))
                and os.path.islink(os.path.join(root, f)),
                files,
            ):
                if project.allow_symlinks:
                    os.unlink(os.path.join(root, file))
                else:
                    os.remove(os.path.join(root, file))
        if "files_symlinked" not in job.document:
            job.document["files_symlinked"] = {}
        job.document["files_symlinked"][project.class_name()] = False

    @classmethod
    def init_and_get_project(cls) -> MartiniTypeFlow:
        if not os.path.isdir(cls.workspace_path):
            os.makedirs(cls.workspace_path)
        signac.init_project(path=cls.workspace_path)
        return cls.get_project(path=cls.workspace_path)


def store_gromacs_log_to_doc(operation_name: str, job: Job):
    """
    Stores the GROMACS log file information in the job document without state point differentiation.

    This function is a simplified wrapper around `_store_gromacs_log_to_doc_flexible`, specifically designed
    for cases where differentiation based on state points is not required. It logs the GROMACS simulation
    details into the job's document, using only the operation name for identification.

    Parameters:
        operation_name (str): The name of the operation being logged. This typically corresponds to the
                              GROMACS command being executed.
        job (Job): The signac job instance for which the log is being stored. This job contains the document
                   where the log information is stored.

    Returns:
        None: This function does not return a value but updates the job's document with the log file information.
    """
    _store_gromacs_log_to_doc_flexible(operation_name, job, False)


def store_gromacs_log_to_doc_with_state_point(operation_name: str, job: Job):
    """
    Stores the GROMACS log file information in the job document with state point differentiation.

    This function is designed to log GROMACS simulation details into the job's document, incorporating
    state point information for differentiation. It utilizes `_store_gromacs_log_to_doc_flexible` with
    the `with_state_point` flag set to True, allowing for the inclusion of state point details in the
    log file's identification.

    Parameters:
        operation_name (str): The name of the operation being logged. This typically corresponds to the
                              GROMACS command being executed.
        job (Job): The signac job instance for which the log is being stored. This job contains the document
                   where the log information is stored.

    Returns:
        None: This function does not return a value but updates the job's document with the log file information,
              including state point differentiation.
    """
    _store_gromacs_log_to_doc_flexible(operation_name, job, True)


def store_gromacs_log_to_doc_with_depth_from_bilayer_core(
    operation_name: str, job: Job
):
    """
    Stores the GROMACS log file information in the job document, including differentiation based on the depth from the bilayer core.

    This function is specifically designed for simulations where the depth from the bilayer core is a relevant parameter. It utilizes
    `_store_gromacs_log_to_doc_flexible` with the `with_state_point` flag set to True and specifies "depth_from_bilayer_core" as the
    state point key for differentiation. This allows for the inclusion of depth-related details in the log file's identification,
    facilitating more granular analysis of simulation results based on their proximity to the bilayer core.

    Parameters:
        operation_name (str): The name of the operation being logged. This typically corresponds to the
                              GROMACS command being executed.
        job (Job): The signac job instance for which the log is being stored. This job contains the document
                   where the log information is stored.

    Returns:
        None: This function does not return a value but updates the job's document with the log file information,
              including differentiation based on the depth from the bilayer core.
    """
    _store_gromacs_log_to_doc_flexible(
        operation_name, job, True, state_point_key="depth_from_bilayer_core"
    )


def _store_gromacs_log_to_doc_flexible(
    operation_name: str,
    job: Job,
    with_state_point: bool,
    state_point_key: str = "lambda_state",
):
    logger.info(f"logging log file for gromacs simulation {operation_name} @ {job.id}")
    project_ = cast(MartiniFlowProject, job.project)
    file_name = (
        f"{operation_name}-{job.sp.get(state_point_key)}"
        if with_state_point
        else operation_name
    )
    job.doc = update_nested_dict(
        job.doc,
        {
            project_.class_name(): {
                "gromacs_logs": {
                    operation_name: project_.get_state_name(file_name, extension="log")
                }
            }
        },
    )
    project_.register_mdp_files(job)
    project_.register_itp_files(job)


def store_task(operation_name: str, job: Job):
    """
    Logs a workflow task for a given job in the project's documentation.

    This function updates the job's document to log the execution of a specific task within the workflow. It is
    part of the project's utilities for tracking the progress and state of jobs, specifically by recording which
    tasks have been run. This aids in the management and review of the job's lifecycle within the workflow.

    Parameters:
        operation_name (str): The name of the operation or task being logged. This name should correspond to a
                              defined operation within the workflow.
        job (Job): The signac job instance for which the task is being logged. This job's document is updated with
                   the task information.

    Returns:
        None: This function does not return a value but updates the job's document with the task information.
    """
    logger.info(f"logging workflow task {operation_name} @ {job.id}")
    project_ = cast(MartiniFlowProject, job.project)
    job.doc = update_nested_dict(
        job.doc, {project_.class_name(): {"tasks": {operation_name: "run"}}}
    )


def store_task_for_many_jobs(operation_name: str, *jobs):
    """
    Logs a workflow task for multiple jobs in the project's documentation.

    This function iterates over a collection of jobs, logging the execution of a specified task within the workflow
    for each job. It leverages the `store_task` function to individually update each job's document. This is useful
    for batch processing or when the same task is executed across multiple jobs, ensuring consistent documentation
    and tracking of workflow tasks across the project.

    Parameters:
        operation_name (str): The name of the operation or task being logged. This name should correspond to a
                              defined operation within the workflow.
        *jobs (Job): A variable number of signac job instances for which the task is being logged. Each job's
                     document is updated with the task information.

    Returns:
        None: This function does not return a value but updates each job's document with the task information.
    """
    for job in jobs:
        store_task(operation_name, job)


def store_workflow(operation_name: str, job: Job):
    """
    Logs the association of a workflow with a given job in the project's documentation.

    This function updates the job's document to log the association of a specific workflow, identified by the
    operation name, with the job. It is part of the project's utilities for tracking the progress and state of
    jobs, specifically by recording which workflows have been associated with them. This aids in the management
    and review of the job's lifecycle within the workflow system.

    Parameters:
        operation_name (str): The name of the operation or workflow being logged. This name should correspond to a
                              defined operation within the project's workflow system.
        job (Job): The signac job instance for which the workflow association is being logged. This job's document
                   is updated with the workflow association information.

    Raises:
        ValueError: If the operation name is not registered within the project's operation to workflow mapping,
                    indicating that the operation is not recognized by the project.

    Returns:
        None: This function does not return a value but updates the job's document with the workflow association
              information.
    """
    project_ = cast(MartiniFlowProject, job.project)
    logger.info(
        f"logging workflow {project_.class_name()} for {operation_name} @ {job.id}"
    )
    if operation_name not in project_.operation_to_workflow:
        raise ValueError(
            f"operation {operation_name} @ {job.id} has not been registered"
        )
    job.doc = update_nested_dict(
        job.doc,
        {
            project_.class_name(): {
                "workflows": {
                    operation_name: project_.operation_to_workflow[operation_name]
                }
            }
        },
    )


def store_workflow_for_many_jobs(operation_name: str, *jobs):
    """
    Logs the association of a specified workflow with multiple jobs in the project's documentation.

    This function iterates over a collection of jobs, logging the association of a specified workflow, identified by the
    operation name, with each job. It leverages the `store_workflow` function to individually update each job's document.
    This is useful for batch processing or when the same workflow is associated across multiple jobs, ensuring consistent
    documentation and tracking of workflow associations across the project.

    Parameters:
        operation_name (str): The name of the operation or workflow being logged. This name should correspond to a
                              defined operation within the project's workflow system.
        *jobs (Job): A variable number of signac job instances for which the workflow association is being logged.
                     Each job's document is updated with the workflow association information.

    Returns:
        None: This function does not return a value but updates each job's document with the workflow association
              information.
    """
    for job in jobs:
        store_workflow(operation_name, job)


def flag_ready_for_upload(_: str, job: Job):
    project_name = cast("MartiniFlowProject", job.project).class_name()
    job.doc = update_nested_dict(
        job.doc, {project_name: {"ready_for_nomad_upload": True}}
    )


def flag_ready_for_upload_multiple_jobs(_: str, *jobs):
    for job in jobs:
        flag_ready_for_upload(_, job)


@MartiniFlowProject.label
def is_ready_for_upload(job: Job) -> bool:
    project_name = cast("MartiniFlowProject", job.project).class_name()
    if project_name not in job.doc:
        return False
    return job.doc[project_name].get("ready_for_nomad_upload", False)


@MartiniFlowProject.label
def uploaded_to_nomad(job: Job) -> bool:
    from martignac.nomad.entries import get_entries_of_upload

    project_name = cast("MartiniFlowProject", job.project).class_name()
    if not UPLOAD_TO_NOMAD:
        return True
    project = cast("MartiniFlowProject", job.project)
    if job.doc.get("nomad_dataset_id", "") != project.nomad_dataset_id:
        logger.info(
            f"inconsistent dataset_id ({job.doc.get('nomad_dataset_id', '')} vs {project.nomad_dataset_id}"
        )
        return False
    if project_name not in job.doc:
        logger.info(f"missing project name {project_name} in job.doc")
        return False
    nomad_upload_id = job.doc[project_name].get("nomad_upload_id", "")
    if not nomad_upload_id:
        logger.info(f"missing nomad_upload_id in job.doc[{project_name}]")
        return False
    try:
        nomad_entries = get_entries_of_upload(
            nomad_upload_id, project.nomad_use_prod_database, with_authentication=True
        )
    except ValueError:
        logger.info(f"retrieving entries of upload {nomad_upload_id} on NOMAD failed")
        return False
    return any(e.job_id == job.id for e in nomad_entries)


@MartiniFlowProject.label
def fetched_from_nomad(*jobs):
    output_flag = True
    for job in jobs:
        project_name = cast("MartiniFlowProject", job.project).class_name()
        if project_name not in job.doc:
            return False
        output_flag = job.doc[project_name].get("fetched_nomad", False)
        if not output_flag:
            return False
    return output_flag


@MartiniFlowProject.label
def itp_mdp_files_symlinked(job: Job):
    project_name = cast("MartiniFlowProject", job.project).class_name()
    if project_name not in job.doc:
        return False
    return job.doc[project_name].get("files_symlinked", False)


@MartiniFlowProject.pre(itp_mdp_files_symlinked)
@MartiniFlowProject.post(fetched_from_nomad)
@MartiniFlowProject.operation(with_job=True)
def fetch_from_nomad(*jobs) -> None:
    """
    Fetches data for a given set of jobs from the NOMAD repository.

    This operation attempts to download raw data associated with each job from the NOMAD repository. It logs
    the success or failure of data retrieval for each job. If data is successfully fetched, it updates the job's
    document to reflect this. This function is useful for synchronizing local job data with data stored in NOMAD,
    ensuring that the local project state accurately reflects the data stored in the repository.

    Parameters:
        *jobs (Job): A variable number of signac job instances for which data is being fetched from NOMAD.

    Returns:
        None: This function does not return a value but updates each job's document with the fetch status.
    """
    from martignac.nomad.entries import download_raw_data_of_job

    for job in jobs:
        logger.info(f"Attempting to fetch job {job.id} from NOMAD")
        initialize_job_doc(job)
        result = download_raw_data_of_job(job)
        if result:
            logger.info(f"Remote data was found on NOMAD for {job.id}")
        else:
            logger.info(f"No remote data was found on NOMAD for {job.id}")
        project_name = cast("MartiniFlowProject", job.project).class_name()
        job.doc[project_name]["fetched_nomad"] = True


@MartiniFlowProject.post(itp_mdp_files_symlinked)
@MartiniFlowProject.operation(with_job=True)
def symlink_itp_and_mdp_files(job: Job) -> None:
    """
    Creates symbolic links for ITP and MDP files in the job's directory.

    This function generates symbolic links for each ITP (GROMACS topology) and MDP (GROMACS parameters) file
    within the job's directory, facilitating their use in simulation workflows. It ensures that the necessary
    simulation input files are accessible in the job's working directory without duplicating data. If a file
    already exists in the job's directory, it will not create a duplicate symlink.

    Parameters:
        job (Job): The signac job instance for which the symbolic links are being created.

    Returns:
        None: This function does not return a value but updates the job's directory with symbolic links.
    """
    project = cast("MartiniFlowProject", job.project)
    logger.info(f"generating symbolic links for {job.id} @ {project.class_name()}")

    def symlink_or_copy(files, _path):
        for file in files:
            if not job.isfile(os.path.basename(file)):
                if project.allow_symlinks:
                    os.symlink(f"{_path}/{file}", job.fn(os.path.basename(file)))
                else:
                    shutil.copy(f"{_path}/{file}", job.fn(os.path.basename(file)))

    symlink_or_copy(project.itp_files.values(), project.itp_path)
    symlink_or_copy(project.mdp_files.values(), project.mdp_path)

    job.doc = update_nested_dict(
        job.doc, {project.class_name(): {"files_symlinked": True}}
    )


def import_job_from_other_flow(
    job: Job,
    child_project: MartiniTypeFlow,
    child_job: Job,
    keys_for_files_to_copy: list[str],
    run_child_job: bool = True,
) -> None:
    """
    Imports data from a job in a child project into the current job's context.

    This function is designed to facilitate the transfer of data between jobs across different projects within
    the MartiniFlow framework. It optionally runs the child job before copying specified files from the child
    job's directory to the current job's directory. This is particularly useful for workflows that require
    data generated in one project to be used in another.

    Parameters:
        job (Job): The current signac job instance into which data is being imported.
        child_project (MartiniTypeFlow): The child project instance from which data is being imported.
        child_job (Job): The child project's signac job instance from which data is being imported.
        keys_for_files_to_copy (list[str]): A list of keys identifying the files in the child job's document
                                            that should be copied to the current job's directory.
        run_child_job (bool, optional): A flag indicating whether the child job should be executed before
                                        data is imported. Defaults to True.

    Returns:
        None: This function does not return a value but updates the current job's directory with the imported
              files and updates the job's document with information from the child job's document.
    """
    if run_child_job:
        logger.info(f"Running job {child_job.id} @  {child_project.class_name()}")
        child_project.run(jobs=[child_job])
        logger.info(
            f"Finished running job {child_job.id} @ {child_project.class_name()}"
        )
    for key in keys_for_files_to_copy:
        shutil.copy(
            child_job.fn(child_job.doc[child_project.class_name()].get(key)), job.path
        )
    job.doc = update_nested_dict(
        job.doc, {child_project.class_name(): child_job.doc[child_project.class_name()]}
    )


@MartiniFlowProject.label
def system_generated(job):
    project = cast("MartiniFlowProject", job.project)
    return job.isfile(project.get_state_name("generate", "gro"))


@MartiniFlowProject.label
def system_minimized(job):
    project = cast("MartiniFlowProject", job.project)
    return job.isfile(project.get_state_name("minimize", "gro"))


@MartiniFlowProject.label
def system_equilibrated(job):
    project = cast("MartiniFlowProject", job.project)
    return job.isfile(project.get_state_name("equilibrate", "gro"))


@MartiniFlowProject.label
def system_sampled(job):
    project = cast("MartiniFlowProject", job.project)
    return job.isfile(project.get_state_name("production", "gro"))


def initialize_job_doc(job: Job) -> None:
    project_name = cast("MartiniFlowProject", job.project).class_name()
    job.doc[project_name]["fetched_nomad"] = False
    job.doc[project_name]["gromacs_logs"] = {}
    job.doc[project_name]["itp_files"] = ""
    job.doc[project_name]["mdp_files"] = ""
    job.doc[project_name]["nomad_upload_id"] = ""
    job.doc[project_name]["nomad_workflow"] = ""
    job.doc[project_name]["ready_for_nomad_upload"] = False
    job.doc[project_name]["tasks"] = {}

import json
import logging
import os
import shutil
from hashlib import md5
from typing import Any, Dict, TypeVar, cast

from flow import FlowProject
from signac.job import Job

from martignac import config
from martignac.nomad.datasets import get_dataset_by_id
from martignac.nomad.uploads import upload_files_to_nomad
from martignac.utils.misc import update_nested_dict, zip_directories
from martignac.utils.nomad import generate_user_metadata

logger = logging.getLogger(__name__)

MartiniTypeFlow = TypeVar("MartiniTypeFlow", bound="MartiniFlowProject")


class MartiniFlowProject(FlowProject):
    workspaces_path: str = config()["local"]["workspaces"]["absolute_path"].get(str)
    input_files_path: str = config()["local"]["input_files"]["absolute_path"].get(str)
    nomad_use_prod_database: bool = config()["nomad"]["use_prod"].get(bool)
    nomad_dataset_id: str = config()["nomad"]["dataset"]["id"].get(str)
    nomad_coauthors: list[str] = [c.get(str) for c in config()["nomad"]["coauthors"]]
    workspace_path: str = ""
    itp_path: str
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
            job.doc, {cls.class_name(): {"mdp_files": cls.get_hash_for_files(job, list(cls.mdp_files.values()))}}
        )

    @classmethod
    def register_itp_files(cls, job: Job) -> None:
        job.doc = update_nested_dict(
            job.doc, {cls.class_name(): {"itp_files": cls.get_hash_for_files(job, list(cls.itp_files.values()))}}
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
            "state_point": dict(job.sp),
            "mdp_files": job.doc[cls.class_name()]["mdp_files"],
            "itp_files": job.doc[cls.class_name()]["itp_files"],
        }

    @classmethod
    def upload_to_nomad(cls, job: Job) -> None:
        if upload_id := job.doc[cls.class_name()].get("nomad_upload_id"):
            logger.info(f"Workflow {cls.class_name()} already uploaded to {upload_id}")
            return None
        dataset = get_dataset_by_id(cls.nomad_dataset_id, use_prod=cls.nomad_use_prod_database)
        logger.info(f"made connection to dataset {dataset.dataset_id}")
        generate_user_metadata(
            file_name=job.fn("nomad.yaml"),
            comment=json.dumps(cls.nomad_comment(job)),
            coauthors=cls.nomad_coauthors,
            datasets=[dataset.dataset_id],
        )
        zip_file = zip_directories([job.path], job.id)
        upload_id = upload_files_to_nomad(zip_file, cls.nomad_use_prod_database)
        job.doc = update_nested_dict(
            job.doc, {"nomad_dataset_id": dataset.dataset_id, cls.class_name(): {"nomad_upload_id": upload_id}}
        )
        os.remove(zip_file)

    @classmethod
    def upload_to_nomad_multiple_jobs(cls, jobs: list[Job]) -> None:
        for job in jobs:
            if upload_id := job.doc[cls.class_name()].get("nomad_upload_id"):
                logger.info(f"Workflow {cls.class_name()} already uploaded to {upload_id}")
                return None
        dataset = get_dataset_by_id(cls.nomad_dataset_id, use_prod=cls.nomad_use_prod_database)
        logger.info(f"made connection to dataset {dataset.dataset_id}")
        generate_user_metadata(
            file_name=job.fn("nomad.yaml"),
            comment=json.dumps(cls.nomad_comment(job)),
            coauthors=cls.nomad_coauthors,
            datasets=[dataset.dataset_id],
        )
        zip_file = zip_directories([job.path for job in jobs], job.id)
        upload_id = upload_files_to_nomad(zip_file, cls.nomad_use_prod_database)
        job.doc = update_nested_dict(
            job.doc, {"nomad_dataset_id": dataset.dataset_id, cls.class_name(): {"nomad_upload_id": upload_id}}
        )
        os.remove(zip_file)

    @classmethod
    def unlink_itp_and_mdp_files(cls, job: Job) -> None:
        project = cast("MartiniFlowProject", job.project)
        logger.info(f"removing symbolic links for {job.id} @ {project.class_name()}")
        for root, _, files in os.walk(job.path):
            for file in filter(lambda f: f.endswith((".itp", ".mdp")) and os.path.islink(os.path.join(root, f)), files):
                os.unlink(os.path.join(root, file))
        if "files_symlinked" not in job.document:
            job.document["files_symlinked"] = {}
        job.document["files_symlinked"][project.class_name()] = False


def store_gromacs_log_to_doc(operation_name: str, job: Job):
    _store_gromacs_log_to_doc_flexible(operation_name, job, False)


def store_gromacs_log_to_doc_with_state_point(operation_name: str, job: Job):
    _store_gromacs_log_to_doc_flexible(operation_name, job, True)


def _store_gromacs_log_to_doc_flexible(operation_name: str, job: Job, with_state_point):
    logger.info(f"logging log file for gromacs simulation {operation_name} @ {job.id}")
    project_ = cast(MartiniFlowProject, job.project)
    file_name = f"{operation_name}-{job.sp.lambda_state}" if with_state_point else operation_name
    job.doc = update_nested_dict(
        job.doc,
        {
            project_.class_name(): {
                "gromacs_logs": {operation_name: project_.get_state_name(file_name, extension="log")}
            }
        },
    )
    project_.register_mdp_files(job)
    project_.register_itp_files(job)


def store_task(operation_name: str, job: Job):
    logger.info(f"logging workflow task {operation_name} @ {job.id}")
    project_ = cast(MartiniFlowProject, job.project)
    job.doc = update_nested_dict(job.doc, {project_.class_name(): {"tasks": {operation_name: "run"}}})


def store_task_for_many_jobs(operation_name: str, *jobs):
    for job in jobs:
        store_task(operation_name, job)


def store_workflow(operation_name: str, job: Job):
    project_ = cast(MartiniFlowProject, job.project)
    logger.info(f"logging workflow {project_.class_name()} for {operation_name} @ {job.id}")
    if operation_name not in project_.operation_to_workflow:
        raise ValueError(f"operation {operation_name} @ {job.id} has not been registered")
    job.doc = update_nested_dict(
        job.doc,
        {project_.class_name(): {"workflows": {operation_name: project_.operation_to_workflow[operation_name]}}},
    )


def store_workflow_for_many_jobs(operation_name: str, *jobs):
    for job in jobs:
        store_workflow(operation_name, job)


def flag_ready_for_upload(_: str, job: Job):
    project_name = cast("MartiniFlowProject", job.project).class_name()
    job.doc = update_nested_dict(job.doc, {project_name: {"ready_for_nomad_upload": True}})


@MartiniFlowProject.label
def is_ready_for_upload(job: Job) -> bool:
    project_name = cast("MartiniFlowProject", job.project).class_name()
    if project_name not in job.doc:
        return False
    return job.doc[project_name].get("ready_for_nomad_upload", False)


@MartiniFlowProject.label
def uploaded_to_nomad(job: Job) -> bool:
    project_name = cast("MartiniFlowProject", job.project).class_name()
    return job.doc[project_name].get("nomad_upload_id", False)


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
    from martignac.nomad.entries import download_raw_data_of_job

    for job in jobs:
        logger.info(f"Attempting to fetch job {job.id} from NOMAD")
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
    project = cast("MartiniFlowProject", job.project)
    logger.info(f"generating symbolic links for {job.id} @ {project.class_name()}")
    for itp_file in project.itp_files.values():
        if not job.isfile(os.path.basename(itp_file)):
            os.symlink(f"{project.itp_path}/{itp_file}", job.fn(os.path.basename(itp_file)))
    for mdp_file in project.mdp_files.values():
        os.symlink(f"{project.mdp_path}/{mdp_file}", job.fn(os.path.basename(mdp_file)))
    job.doc = update_nested_dict(job.doc, {project.class_name(): {"files_symlinked": True}})


def import_job_from_other_flow(
    job: Job,
    child_project: MartiniTypeFlow,
    child_job: Job,
    keys_for_files_to_copy: list[str],
    run_child_job: bool = True,
) -> None:
    if run_child_job:
        logger.info(f"Running job {child_job.id} @  {child_project.class_name()}")
        child_project.run(jobs=[child_job])
        logger.info(f"Finished running job {child_job.id} @ {child_project.class_name()}")
    for key in keys_for_files_to_copy:
        shutil.copy(child_job.fn(child_job.doc[child_project.class_name()].get(key)), job.path)
    job.doc = update_nested_dict(job.doc, {child_project.class_name(): child_job.doc[child_project.class_name()]})

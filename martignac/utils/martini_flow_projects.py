import json
import logging
import os
import shutil
from functools import wraps
from hashlib import md5
from typing import Any, Dict, Optional, TypeVar, cast

from flow import FlowProject
from signac.job import Job

from martignac import config
from martignac.nomad.datasets import get_dataset_by_id
from martignac.nomad.uploads import upload_files_to_nomad
from martignac.utils.misc import update_nested_dict, zip_directory
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
    ff_parameters: Dict[str, Any] = {}

    @classmethod
    def get_mdp_hash(cls, job: Job) -> str:
        readouts = []
        for filename in cls.mdp_files.values():
            with open(job.fn(filename)) as fp:
                readouts.append(fp.read())
        return md5(" ".join(readouts).encode("utf-8")).hexdigest()

    @classmethod
    def register_mdp_files(cls, job: Job) -> None:
        if "mdp_files" in job.document:
            return
        if "mdp_files" not in job.document:
            job.document["mdp_files"] = {}
        job.document["mdp_files"][cls.__name__] = cls.get_mdp_hash(job)

    @classmethod
    def get_state_name(cls, state_name: str = "", extension: str = "") -> str:
        state_ext = f"_{state_name}" if state_name else ""
        file_ext = f".{extension}" if extension else ""
        return f"{cls.system_name}{state_ext}{file_ext}"

    @classmethod
    def log_gromacs_simulation(cls, simulation_name=None, add_lambda_state: bool = False):
        def decorator(func):
            @wraps(func)
            def wrapper(job, *args, **kwargs):
                if job:
                    class_name = cls.__name__
                    if "gromacs_logs" not in job.document:
                        job.document["gromacs_logs"] = {}
                    if class_name not in job.document["gromacs_logs"]:
                        job.document["gromacs_logs"][class_name] = {}
                    task_name = func.__name__
                    sim_name = simulation_name if simulation_name is not None else task_name
                    file_name = f"{sim_name}{job.sp.lambda_state}.log" if add_lambda_state else f"{sim_name}.log"
                    job.document["gromacs_logs"][class_name][task_name] = file_name
                    cls.register_mdp_files(job)
                    return func(job, *args, **kwargs)
                else:
                    raise ValueError("Job argument is required.")

            return wrapper

        return decorator

    @classmethod
    def nomad_comment(cls, job: Job) -> dict:
        return {
            "job_id": job.id,
            "workflow_name": cls.__name__,
            "state_point": dict(job.sp),
            "mdp_files": job.document["mdp_files"][cls.__name__],
        }

    @classmethod
    def upload_to_nomad(cls, job: Job) -> None:
        if not job.document.get("upload_to_nomad", True):
            return None
        dataset = get_dataset_by_id(cls.nomad_dataset_id, use_prod=cls.nomad_use_prod_database)
        logger.info(f"made connection to dataset {dataset.dataset_id}")
        generate_user_metadata(
            file_name=job.fn("nomad.yaml"),
            comment=json.dumps(cls.nomad_comment(job)),
            coauthors=cls.nomad_coauthors,
            datasets=[dataset.dataset_id],
        )
        zip_file = zip_directory(job.path, job.id)
        upload_id = upload_files_to_nomad(zip_file, cls.nomad_use_prod_database)
        job.document["nomad_dataset_id"] = dataset.dataset_id
        if "nomad_upload_id" not in job.document:
            job.document["nomad_upload_id"] = {}
        job.document["nomad_upload_id"][cls.__name__] = upload_id
        os.remove(zip_file)

    @classmethod
    def unlink_itp_and_mdp_files(cls, job: Job) -> None:
        project = cast("MartiniFlowProject", job.project)
        logger.info(f"removing symbolic links for {job.id} @ {project.__class__.__name__}")
        for root, _, files in os.walk(job.path):
            for file in filter(lambda f: f.endswith((".itp", ".mdp")) and os.path.islink(os.path.join(root, f)), files):
                os.unlink(os.path.join(root, file))
        if "files_symlinked" not in job.document:
            job.document["files_symlinked"] = {}
        job.document["files_symlinked"][project.__class__.__name__] = False


@MartiniFlowProject.label
def uploaded_to_nomad(job: Job) -> bool:
    if job.document.get("nomad_upload_id"):
        return job.document["nomad_upload_id"].get(job.project.__class__.__name__, False)
    return False


@MartiniFlowProject.label
def fetched_from_nomad(job: Job):
    return job.document.get("fetched_nomad", False)


@MartiniFlowProject.label
def itp_mdp_files_symlinked(job: Job):
    project = cast("MartiniFlowProject", job.project)
    if "files_symlinked" in job.document:
        return job.document["files_symlinked"].get(project.__class__.__name__, False)
    return False


@MartiniFlowProject.pre(itp_mdp_files_symlinked)
@MartiniFlowProject.post(fetched_from_nomad)
@MartiniFlowProject.operation(with_job=True)
def fetch_from_nomad(job: Job) -> None:
    from martignac.nomad.entries import download_raw_data_of_job

    logger.info(f"Attempting to fetch job {job.id} from NOMAD")
    result = download_raw_data_of_job(job)
    if result:
        logger.info(f"Remote data was found on NOMAD for {job.id}")
    else:
        logger.info(f"No remote data was found on NOMAD for {job.id}")
    job.document["fetched_nomad"] = True


@MartiniFlowProject.post(itp_mdp_files_symlinked)
@MartiniFlowProject.operation(with_job=True)
def symlink_itp_and_mdp_files(job: Job) -> None:
    project = cast("MartiniFlowProject", job.project)
    logger.info(f"generating symbolic links for {job.id} @ {project.__class__.__name__}")
    for itp_file in project.itp_files.values():
        if not job.isfile(os.path.basename(itp_file)):
            os.symlink(f"{project.itp_path}/{itp_file}", job.fn(os.path.basename(itp_file)))
    for mdp_file in project.mdp_files.values():
        os.symlink(f"{project.mdp_path}/{mdp_file}", job.fn(os.path.basename(mdp_file)))
    if "files_symlinked" not in job.document:
        job.document["files_symlinked"] = {}
    job.document["files_symlinked"][project.__class__.__name__] = True


def import_job_from_other_flow(
    job: Job,
    child_project: MartiniTypeFlow,
    child_job: Job,
    keys_for_files_to_copy: list[str],
    update_keys: Optional[list[str]] = None,
    run_child_job: bool = True,
) -> None:
    job.document["upload_to_nomad"] = False
    if run_child_job:
        logger.info(f"Running job {child_job.id} @  {child_project.__class__.__name__}")
        child_project.run(jobs=[child_job])
        logger.info(f"Finished running job {child_job.id} @ {child_project.__class__.__name__}")
    for key in keys_for_files_to_copy:
        shutil.copy(child_job.fn(child_job.document.get(key)), job.path)
    log_files = []
    for cls_name in child_job.doc["gromacs_logs"]:
        log_files = [*log_files, *list(child_job.doc["gromacs_logs"][cls_name].values())]
    for file in log_files:
        shutil.copy(child_job.fn(file), job.path)
    update_keys = ["gromacs_logs", "nomad_workflow", "nomad_upload_id", *update_keys]
    job.document = update_nested_dict(
        job.document, {k: child_job.document.get(k) for k in update_keys if k in child_job.document}
    )
    job.document["upload_to_nomad"] = True

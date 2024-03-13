import json
import logging
import os
from functools import wraps
from hashlib import md5
from typing import Any, Dict

from flow import FlowProject
from signac.job import Job

from martignac import config
from martignac.nomad.datasets import get_dataset_by_id
from martignac.nomad.uploads import upload_files_to_nomad
from martignac.utils.misc import zip_directory
from martignac.utils.nomad import generate_user_metadata

logger = logging.getLogger(__name__)


class MartiniFlowProject(FlowProject):
    nomad_use_prod_database: bool = config()["nomad"]["use_prod"].get(bool)
    nomad_dataset_id: str = config()["nomad"]["dataset"]["id"].get(str)
    nomad_coauthors: list[str] = [c.get(str) for c in config()["nomad"]["coauthors"]]
    mdp_files: Dict[str, str]
    itp_files: Dict[str, str]
    simulation_settings: Dict[str, Any]
    system_name: str
    output_names: Dict[str, str]
    ff_parameters: Dict[str, Any]

    @classmethod
    def register_mdp_files(cls, job: Job) -> None:
        if "mdp_files" not in job.document:
            job.document["mdp_files"] = {}
        for mdp_name, mdp_file in cls.mdp_files.items():
            if mdp_name not in job.document["mdp_files"]:
                with open(mdp_file) as f:
                    job.document["mdp_files"][mdp_name] = md5(f.read().encode("utf-8")).hexdigest()

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
            "mdp_files": {k: v for k, v in job.document["mdp_files"].items()},
        }

    @classmethod
    def upload_to_nomad(cls, job: Job) -> None:
        if not job.document.get("upload_to_nomad", True):
            return None
        dataset = get_dataset_by_id(cls.nomad_dataset_id, use_prod=cls.nomad_use_prod_database)
        logger.info(f"made connection to dataset {dataset}")
        generate_user_metadata(
            file_name=job.fn("nomad.yaml"),
            comment=json.dumps(cls.nomad_comment(job)),
            coauthors=cls.nomad_coauthors,
            datasets=[dataset.dataset_id],
        )
        zip_file = zip_directory(job.path, job.id)
        upload_id = upload_files_to_nomad(zip_file, cls.nomad_use_prod_database)
        job.document["nomad_dataset_id"] = dataset.dataset_id
        job.document["nomad_upload_id"] = upload_id
        os.remove(zip_file)

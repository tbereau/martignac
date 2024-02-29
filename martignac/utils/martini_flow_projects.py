import logging
import os
from functools import wraps

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
                    return func(job, *args, **kwargs)
                else:
                    raise ValueError("Job argument is required.")

            return wrapper

        return decorator

    @classmethod
    def upload_to_nomad(cls, job: Job) -> None:
        if not job.document.get("upload_to_nomad", True):
            return None
        dataset = get_dataset_by_id(cls.nomad_dataset_id, use_prod=cls.nomad_use_prod_database)
        logger.info(f"made connection to dataset {dataset}")
        generate_user_metadata(
            file_name=job.fn("nomad.yaml"),
            comment=job.id,
            coauthors=cls.nomad_coauthors,
            datasets=[dataset.dataset_id],
        )
        zip_file = zip_directory(job.path, job.id)
        upload_id = upload_files_to_nomad(zip_file, cls.nomad_use_prod_database)
        job.document["nomad_dataset_id"] = dataset.dataset_id
        job.document["nomad_upload_id"] = upload_id
        os.remove(zip_file)

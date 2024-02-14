import datetime as dt
import logging
from typing import Optional

from marshmallow import Schema, pre_load
from marshmallow_dataclass import class_schema, dataclass

from martignac.nomad.users import NomadUser, get_user_by_id
from martignac.nomad.utils import (
    delete_nomad_request,
    get_nomad_request,
    post_nomad_request,
)

logger = logging.getLogger(__name__)


class NomadDatasetSchema(Schema):
    @pre_load
    def convert_users(self, data, **kwargs):
        data["user"] = get_user_by_id(user_id=data["user_id"]).as_dict()
        del data["user_id"]
        return data


@dataclass(frozen=True)
class NomadDataset:
    dataset_id: str
    dataset_create_time: dt.datetime
    dataset_name: str
    dataset_type: Optional[str] = None
    dataset_modified_time: Optional[dt.datetime] = None
    user: Optional[NomadUser] = None
    doi: Optional[str] = None
    pid: Optional[int] = None
    m_annotations: Optional[dict] = None


def retrieve_datasets(
    dataset_id: str = None,
    dataset_name: str = None,
    user_id: str = None,
    page_size: int = 10,
    max_datasets: int = 50,
    use_prod: bool = True,
) -> list[NomadDataset]:
    parameters = []
    if dataset_id:
        parameters.append(f"dataset_id={dataset_id}")
    if dataset_name:
        parameters.append(f"dataset_name={dataset_name}")
    if user_id:
        parameters.append(f"user_id={user_id}")
    parameters.append(f"page_size={page_size}")
    url_suffix = "/datasets/"
    if len(parameters) > 0:
        url_suffix += f"?{parameters[0]}"
    for i in range(1, len(parameters)):
        url_suffix += f"&{parameters[i]}"
    headers = {"Accept": "application/json"}
    nomad_entry_schema = class_schema(NomadDataset, base_schema=NomadDatasetSchema)
    datasets = []
    page_after_value = None
    while (max_datasets > 0 and len(datasets) <= max_datasets) or (max_datasets < 0):
        url = f"{url_suffix}&page_after_value={page_after_value}" if page_after_value else url_suffix
        response = get_nomad_request(url, headers=headers, use_prod=use_prod)
        if len(response["data"]) == 0:
            break
        datasets.extend([nomad_entry_schema().load(d) for d in response["data"]])
        if response["pagination"]["page"] == response["pagination"]["total"]:
            break
        page_after_value = response["pagination"]["next_page_after_value"]
    return datasets


def create_dataset(dataset_name: str, use_prod: bool = False, timeout_in_sec: int = 10) -> str:
    logger.info(f"creating dataset name {dataset_name} on {'prod' if use_prod else 'test'} server")
    json_dict = {"dataset_name": dataset_name}
    response = post_nomad_request(
        "/datasets/",
        with_authentication=True,
        json_dict=json_dict,
        use_prod=use_prod,
        timeout_in_sec=timeout_in_sec,
    )
    return response.get("dataset_id")


def delete_dataset(dataset_id: str, use_prod: bool = False, timeout_in_sec: int = 10) -> None:
    logger.info(f"deleting dataset {dataset_id} on {'prod' if use_prod else 'test'} server")
    response = delete_nomad_request(
        f"/datasets/{dataset_id}",
        with_authentication=True,
        use_prod=use_prod,
        timeout_in_sec=timeout_in_sec,
    )
    if response.get("dataset_id"):
        logger.info(f"successfully deleted dataset {dataset_id}")
    else:
        logger.error("no dataset deleted")

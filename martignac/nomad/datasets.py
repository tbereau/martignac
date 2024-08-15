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
    """
    Represents a dataset within the NOMAD system.

    This class defines the structure of a dataset object used in the NOMAD application, including its metadata and
    associated user information. It is designed to be immutable, with all attributes set at the time of instantiation.

    Attributes:
        dataset_id (str): Unique identifier for the dataset.
        dataset_create_time (dt.datetime): The creation time of the dataset.
        dataset_name (str): The name of the dataset.
        dataset_type (Optional[str]): The type of the dataset, if specified.
        dataset_modified_time (Optional[dt.datetime]): The last modification time of the dataset, if any.
        user (Optional[NomadUser]): The user associated with the dataset, if any.
        doi (Optional[str]): The Digital Object Identifier (DOI) of the dataset, if any.
        pid (Optional[int]): The persistent identifier (PID) of the dataset, if any.
        m_annotations (Optional[dict]): A dictionary of metadata annotations associated with the dataset, if any.
    """

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
    dataset_id: Optional[str] = None,
    dataset_name: Optional[str] = None,
    user_id: Optional[str] = None,
    page_size: int = 10,
    max_datasets: int = 50,
    use_prod: bool = True,
) -> list[NomadDataset]:
    """
    Retrieves a list of NomadDataset objects based on the provided filters.

    This function queries the NOMAD system for datasets, optionally filtering by dataset ID, dataset name, or user ID.
    It supports pagination through the `page_size` parameter and allows limiting the total number of datasets returned
    with `max_datasets`. The `use_prod` flag determines whether to query the production or test environment.

    Args:
        dataset_id (str, optional): The unique identifier of the dataset to retrieve. Defaults to None.
        dataset_name (str, optional): The name of the dataset to filter by. Defaults to None.
        user_id (str, optional): The user ID to filter datasets by. Defaults to None.
        page_size (int, optional): The number of datasets to return per page. Defaults to 10.
        max_datasets (int, optional): The maximum number of datasets to retrieve. Defaults to 50.
        use_prod (bool, optional): Flag to use the production environment. Defaults to True.

    Returns:
        list[NomadDataset]: A list of NomadDataset objects matching the query.
    """
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
        url = (
            f"{url_suffix}&page_after_value={page_after_value}"
            if page_after_value
            else url_suffix
        )
        response = get_nomad_request(url, headers=headers, use_prod=use_prod)
        if len(response["data"]) == 0:
            break
        datasets.extend([nomad_entry_schema().load(d) for d in response["data"]])
        if response["pagination"]["page"] == response["pagination"]["total"]:
            break
        page_after_value = response["pagination"]["next_page_after_value"]
    return datasets


def get_dataset_by_id(dataset_id: str, use_prod: bool = True) -> NomadDataset:
    """
    Retrieves a single NomadDataset object by its dataset ID.

    This function queries the NOMAD system for a dataset with the specified ID. It leverages the `retrieve_datasets`
    function to perform the query, ensuring that only one dataset is returned. If the query returns more or fewer
    than one dataset, it raises a ValueError indicating an issue with the retrieval process.

    Args:
        dataset_id (str): The unique identifier of the dataset to retrieve.
        use_prod (bool, optional): Flag to use the production environment. Defaults to True.

    Returns:
        NomadDataset: The dataset object corresponding to the provided ID.

    Raises:
        ValueError: If no dataset is found with the provided ID, or if multiple datasets are returned.
    """
    datasets = retrieve_datasets(dataset_id=dataset_id, use_prod=use_prod)
    if len(datasets) != 1:
        raise ValueError(f"Problem retrieving dataset {dataset_id}: {datasets}")
    return datasets[0]


def create_dataset(
    dataset_name: str, use_prod: bool = False, timeout_in_sec: int = 10
) -> str:
    """
    Creates a new dataset in the NOMAD system with the specified name.

    This function sends a POST request to the NOMAD system to create a new dataset. The request includes the dataset
    name and is sent to either the production or test environment based on the `use_prod` flag. The function waits for
    a response for a specified timeout period.

    Args:
        dataset_name (str): The name of the dataset to be created.
        use_prod (bool, optional): Flag indicating whether to use the production environment. Defaults to False,
                                   indicating that the test environment is used by default.
        timeout_in_sec (int, optional): The maximum time in seconds to wait for a response from the server.
                                        Defaults to 10 seconds.

    Returns:
        str: The unique identifier of the newly created dataset, as returned by the NOMAD system.

    Raises:
        HTTPError: If the request fails or the NOMAD system returns an error response.
    """
    logger.info(
        f"creating dataset name {dataset_name} on {'prod' if use_prod else 'test'} server"
    )
    json_dict = {"dataset_name": dataset_name}
    response = post_nomad_request(
        "/datasets/",
        with_authentication=True,
        json_dict=json_dict,
        use_prod=use_prod,
        timeout_in_sec=timeout_in_sec,
    )
    return response.get("dataset_id")


def delete_dataset(
    dataset_id: str, use_prod: bool = False, timeout_in_sec: int = 10
) -> None:
    """
    Deletes a dataset from the NOMAD system by its dataset ID.

    This function sends a DELETE request to the NOMAD system to remove a dataset identified by its unique ID. The
    operation can be directed to either the production or test environment, as specified by the `use_prod` flag. The
    function allows specifying a timeout for the request.

    Args:
        dataset_id (str): The unique identifier of the dataset to be deleted.
        use_prod (bool, optional): Flag indicating whether to use the production environment. Defaults to False,
                                   indicating that the test environment is used by default.
        timeout_in_sec (int, optional): The maximum time in seconds to wait for a response from the server.
                                        Defaults to 10 seconds.

    Note:
        This function logs the outcome of the deletion operation, reporting success or failure through the logging
        system.
    """
    logger.info(
        f"deleting dataset {dataset_id} on {'prod' if use_prod else 'test'} server"
    )
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

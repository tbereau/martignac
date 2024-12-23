import datetime as dt
import logging
from collections.abc import ByteString
from os.path import basename
from time import sleep
from typing import Any, Optional, Union

from cachetools.func import ttl_cache
from marshmallow import Schema, pre_load
from marshmallow_dataclass import class_schema, dataclass

from martignac.nomad.users import NomadUser, get_user_by_id
from martignac.nomad.utils import (
    delete_nomad_request,
    get_nomad_base_url,
    get_nomad_request,
    post_nomad_request,
    put_nomad_request,
)

logger = logging.getLogger(__name__)


class NomadUploadSchema(Schema):
    @pre_load
    def convert_users(self, data, **kwargs):
        data["main_author"] = get_user_by_id(user_id=data["main_author"]).as_dict()
        data["writers"] = [get_user_by_id(user_id=w).as_dict() for w in data["writers"]]
        data["reviewers"] = [
            get_user_by_id(user_id=r).as_dict() for r in data["reviewers"]
        ]
        data["viewers"] = [get_user_by_id(user_id=v).as_dict() for v in data["viewers"]]
        return data


@dataclass
class NomadUpload:
    """
    Represents an upload in the NOMAD system, encapsulating all relevant metadata and state information.

    Attributes:
        upload_id (str): Unique identifier for the upload.
        upload_create_time (datetime.datetime): The creation time of the upload.
        main_author (NomadUser): The main author of the upload.
        process_running (bool): Flag indicating if a process is currently running for this upload.
        current_process (str): The name of the current process running.
        process_status (str): The status of the current process.
        last_status_message (str): The last status message received for the current process.
        errors (list[Any]): A list of errors associated with the upload.
        warnings (list[Any]): A list of warnings associated with the upload.
        coauthors (list[str]): List of coauthor identifiers.
        coauthor_groups (list[Any]): List of coauthor groups.
        reviewers (list[NomadUser]): List of reviewers.
        reviewer_groups (list[Any]): List of reviewer groups.
        writers (list[NomadUser]): List of writers with access to the upload.
        writer_groups (list[Any]): List of writer groups.
        viewers (list[NomadUser]): List of viewers with access to the upload.
        viewer_groups (list[Any]): List of viewer groups.
        published (bool): Flag indicating if the upload is published.
        published_to (list[Any]): List of platforms or locations the upload is published to.
        with_embargo (bool): Flag indicating if the upload is under embargo.
        embargo_length (float): The length of the embargo in days.
        license (str): The license associated with the upload.
        entries (int): The number of entries in the upload.
        n_entries (Optional[int]): The number of entries, if known.
        upload_files_server_path (Optional[str]): The server path to the uploaded files.
        publish_time (Optional[datetime.datetime]): The time the upload was published.
        references (Optional[list[str]]): List of references associated with the upload.
        datasets (Optional[list[str]]): List of dataset identifiers associated with the upload.
        external_db (Optional[str]): External database identifier.
        upload_name (Optional[str]): The name of the upload.
        comment (Optional[str]): A comment or description of the upload.
        use_prod (Optional[bool]): Flag indicating if the production environment is used.
        complete_time (Optional[datetime.datetime]): The time the upload was completed.
    """

    upload_id: str
    upload_create_time: dt.datetime
    main_author: NomadUser
    process_running: bool
    current_process: str
    process_status: str
    last_status_message: str
    errors: list[Any]
    warnings: list[Any]
    coauthors: list[str]
    coauthor_groups: list[Any]
    reviewers: list[NomadUser]
    reviewer_groups: list[Any]
    writers: list[NomadUser]
    writer_groups: list[Any]
    viewers: list[NomadUser]
    viewer_groups: list[Any]
    published: bool
    published_to: list[Any]
    with_embargo: bool
    embargo_length: float
    license: str
    entries: int
    n_entries: Optional[int] = None
    upload_files_server_path: Optional[str] = None
    publish_time: Optional[dt.datetime] = None
    references: Optional[list[str]] = None
    datasets: Optional[list[str]] = None
    external_db: Optional[str] = None
    upload_name: Optional[str] = None
    comment: Optional[str] = None
    use_prod: Optional[bool] = None
    complete_time: Optional[dt.datetime] = None

    @property
    def base_url(self) -> Optional[str]:
        if self.use_prod is not None:
            return get_nomad_base_url(self.use_prod)
        return None

    @property
    def nomad_gui_url(self) -> str:
        if self.base_url is None:
            raise ValueError(f"missing attribute 'use_prod' for upload {self}")
        return f"{self.base_url}/gui/user/uploads/upload/id/{self.upload_id}"

    def fetch(self, with_authentication: bool = True) -> None:
        updated_upload = get_upload_by_id(
            self.upload_id,
            use_prod=self.use_prod,
            with_authentication=with_authentication,
        )
        self.__dict__.update(updated_upload.__dict__)

    def safe_publish(self) -> None:
        if not self.published:
            wait_counter = 600  # 10 minutes overall
            while self.process_running and wait_counter > 0:
                logger.info("upload still processing, waiting before new attempt")
                sleep(10)
                self.fetch()
                wait_counter -= 1
            if self.process_running:
                raise ValueError(f"upload {self.upload_id} is still processing")
            publish_upload(upload_id=self.upload_id, use_prod=self.use_prod)


@ttl_cache(maxsize=128, ttl=180)
def get_all_my_uploads(
    use_prod: bool = False, timeout_in_sec: int = 10
) -> list[NomadUpload]:
    """
    Retrieves all uploads associated with the authenticated user from the NOMAD system.

    This function fetches a list of all uploads made by the currently authenticated user. It allows the user to specify
    whether to interact with the production or test environment of the NOMAD system. The function also supports specifying
    a timeout for the network request.

    Args:
        use_prod (bool, optional): Flag indicating whether to use the production environment of the NOMAD system.
                                   Defaults to False, which means the test environment is used by default.
        timeout_in_sec (int, optional): The maximum time in seconds to wait for a response from the NOMAD system.
                                        Defaults to 10 seconds.

    Returns:
        list[NomadUpload]: A list of `NomadUpload` objects, each representing an upload made by the user.
    """
    logger.info(f"retrieving all uploads on {'prod' if use_prod else 'test'} server")
    response = get_nomad_request(
        "/uploads",
        use_prod=use_prod,
        with_authentication=True,
        timeout_in_sec=timeout_in_sec,
    )
    upload_class_schema = class_schema(NomadUpload, base_schema=NomadUploadSchema)
    return [
        upload_class_schema().load({**r, "use_prod": use_prod})
        for r in response["data"]
    ]


def get_upload_by_id(
    upload_id: str,
    use_prod: bool = False,
    with_authentication: bool = True,
    timeout_in_sec: int = 10,
) -> NomadUpload:
    """
    Retrieves a specific upload by its ID from the NOMAD system.

    This function fetches the details of a single upload identified by its unique ID. It allows specifying whether to
    access the production or test environment of the NOMAD system. Additionally, a timeout for the network request can
    be defined.

    Args:
        upload_id (str): The unique identifier of the upload to retrieve.
        use_prod (bool, optional): Flag indicating whether to use the production environment of the NOMAD system.
                                   Defaults to False, which means the test environment is used by default.
        with_authentication: Flag indicating whether to include authentication headers in the request.
        timeout_in_sec (int, optional): The maximum time in seconds to wait for a response from the NOMAD system.
                                        Defaults to 10 seconds.

    Returns:
        NomadUpload: An instance of `NomadUpload` containing all the metadata and state information of the upload.
    """
    logger.info(
        f"retrieving upload {upload_id} on {'prod' if use_prod else 'test'} server"
    )
    response = get_nomad_request(
        f"/uploads/{upload_id}",
        use_prod=use_prod,
        with_authentication=with_authentication,
        timeout_in_sec=timeout_in_sec,
    )
    upload_class_schema = class_schema(NomadUpload, base_schema=NomadUploadSchema)
    return upload_class_schema().load({**response["data"], "use_prod": use_prod})


def delete_upload(
    upload_id: str, use_prod: bool = False, timeout_in_sec: int = 10
) -> NomadUpload:
    """
    Deletes a specific upload from the NOMAD system based on its unique ID.

    This function sends a request to the NOMAD system to delete an upload identified by its unique ID. It allows
    specifying whether to interact with the production or test environment of the NOMAD system. Additionally, a timeout
    for the network request can be defined to avoid indefinite waiting periods.

    Args:
        upload_id (str): The unique identifier of the upload to be deleted.
        use_prod (bool, optional): Flag indicating whether to use the production environment of the NOMAD system.
                                   Defaults to False, which means the test environment is used by default.
        timeout_in_sec (int, optional): The maximum time in seconds to wait for a response from the NOMAD system.
                                        Defaults to 10 seconds.

    Returns:
        NomadUpload: An instance of `NomadUpload` containing the metadata of the deleted upload. This can be useful
                     for logging or confirmation purposes.

    Raises:
        ValueError: If the deletion request fails or the response from the NOMAD system is unexpected.
    """
    logger.info(
        f"deleting upload {upload_id} on {'prod' if use_prod else 'test'} server"
    )
    response = delete_nomad_request(
        f"/uploads/{upload_id}",
        with_authentication=True,
        timeout_in_sec=timeout_in_sec,
    )
    upload_class_schema = class_schema(NomadUpload, base_schema=NomadUploadSchema)
    return upload_class_schema().load({**response["data"], "use_prod": use_prod})


def upload_files_to_nomad(
    filename: str, use_prod: bool = False, timeout_in_sec: int = 30
) -> str:
    """
    Uploads a file to the NOMAD system.

    This function uploads a specified file to the NOMAD system, allowing the user to choose between the production
    and test environments. It also supports setting a custom timeout for the upload process to prevent indefinite
    waiting periods. The function logs the outcome of the upload process, including success or failure messages.

    Args:
        filename (str): The path to the file to be uploaded.
        use_prod (bool, optional): Flag indicating whether to use the production environment of the NOMAD system.
                                   Defaults to False, which means the test environment is used by default.
        timeout_in_sec (int, optional): The maximum time in seconds to wait for the upload to complete.
                                        Defaults to 30 seconds.

    Returns:
        str: The unique identifier of the upload if successful, otherwise logs an error.

    Raises:
        IOError: If the file cannot be opened or read.
        ValueError: If the response from the NOMAD system does not contain an upload ID.
    """
    logger.info(f"uploading file {filename} on {'prod' if use_prod else 'test'} server")
    with open(filename, "rb") as f:
        response = post_nomad_request(
            "/uploads",
            with_authentication=True,
            data=f,
            use_prod=use_prod,
            timeout_in_sec=timeout_in_sec,
        )
    upload_id = response.get("upload_id")
    if upload_id:
        logger.info(f"successful upload to {upload_id}")
        return upload_id
    else:
        logger.error(f"could not upload {filename}. Response {response}")


def publish_upload(
    upload_id: str, use_prod: bool = False, timeout_in_sec: int = 10
) -> dict:
    """
    Publishes a specified upload in the NOMAD system.

    This function sends a request to the NOMAD system to publish an upload identified by its unique ID. It allows
    specifying whether to interact with the production or test environment of the NOMAD system. Additionally, a timeout
    for the network request can be defined to manage how long the function waits for a response from the NOMAD system.

    Args:
        upload_id (str): The unique identifier of the upload to be published.
        use_prod (bool, optional): Flag indicating whether to use the production environment of the NOMAD system.
                                   Defaults to False, which means the test environment is used by default.
        timeout_in_sec (int, optional): The maximum time in seconds to wait for a response from the NOMAD system.
                                        Defaults to 10 seconds.

    Returns:
        dict: A dictionary containing the response from the NOMAD system regarding the publish action.
    """
    logger.info(
        f"publishing upload {upload_id} on {'prod' if use_prod else 'test'} server"
    )
    response = post_nomad_request(
        f"/uploads/{upload_id}/action/publish",
        with_authentication=True,
        timeout_in_sec=timeout_in_sec,
    )
    logger.info(
        f"published upload {upload_id} on {'prod' if use_prod else 'test'} server"
    )
    return response


def edit_upload_metadata(
    upload_id: str,
    upload_name: Optional[str] = None,
    references: Optional[list[str]] = None,
    dataset_id: Optional[str] = None,
    embargo_length: Optional[float] = None,
    coauthors_ids: Optional[list[str]] = None,
    comment: Optional[str] = None,
    use_prod: bool = False,
    timeout_in_sec: int = 10,
) -> dict:
    """
    Edits the metadata of a specific upload in the NOMAD system.

    This function allows for the modification of various metadata fields of an existing upload, identified by its unique ID.
    It supports updating the upload's name, references, dataset ID, embargo length, coauthors, and additional comments.
    The function also allows specifying whether to interact with the production or test environment of the NOMAD system,
    along with a custom timeout for the network request.

    Args:
        upload_id (str): The unique identifier of the upload to be edited.
        upload_name (Optional[str], optional): The new name for the upload.
        references (Optional[list[str]], optional): A list of new references associated with the upload.
        dataset_id (Optional[str], optional): The new dataset ID associated with the upload.
        embargo_length (Optional[float], optional): The new embargo length in days.
        coauthors_ids (Optional[list[str]], optional): A list of new coauthor identifiers.
        comment (Optional[str], optional): A new comment or description for the upload.
        use_prod (bool, optional): Flag indicating whether to use the production environment of the NOMAD system.
                                   Defaults to False, which means the test environment is used by default.
        timeout_in_sec (int, optional): The maximum time in seconds to wait for a response from the NOMAD system.
                                        Defaults to 10 seconds.

    Returns:
        dict: A dictionary containing the response from the NOMAD system regarding the edit action.
    """
    logger.info(
        f"editing the metadata for upload {upload_id} on {'prod' if use_prod else 'test'} server"
    )
    metadata = {"metadata": {}}
    if upload_name:
        metadata["metadata"]["upload_name"] = upload_name
    if references:
        metadata["metadata"]["references"] = references
    if dataset_id:
        metadata["metadata"]["datasets"] = dataset_id
    if embargo_length:
        metadata["metadata"]["embargo_length"] = embargo_length
    if coauthors_ids:
        metadata["metadata"]["coauthors"] = coauthors_ids
    if comment:
        metadata["metadata"]["comment"] = comment
    response = post_nomad_request(
        f"/uploads/{upload_id}/edit",
        use_prod=use_prod,
        with_authentication=True,
        json_dict=metadata,
        timeout_in_sec=timeout_in_sec,
    )
    return response


@ttl_cache(maxsize=128, ttl=180)
def get_specific_file_from_upload(
    upload_id: str,
    path_to_file: str,
    use_prod: bool = False,
    with_authentication: bool = False,
    return_json: bool = True,
    timeout_in_sec: int = 10,
) -> Union[ByteString, dict]:
    """
    Downloads a file associated with a specific upload from the NOMAD system.

    This function retrieves a file associated with a specific upload from the NOMAD system and saves it to the specified
    path on the local filesystem. It allows specifying whether to interact with the production or test environment of
    the NOMAD system. Additionally, a timeout for the network request can be defined to manage how long the function
    waits for the file to be downloaded.

    Args:
        upload_id (str): The unique identifier of the upload containing the file to download.
        path_to_file (str): The path where the downloaded file should be saved on the local filesystem.
        use_prod (bool, optional): Flag indicating whether to use the production environment of the NOMAD system.
                                   Defaults to False, which means the test environment is used by default.
        with_authentication: Flag indicating whether to include authentication headers in the request.
        timeout_in_sec (int, optional): The maximum time in seconds to wait for the file to be downloaded.
                                        Defaults to 10 seconds.
        return_json: Flag indicating whether to return the response as JSON or as a byte string.

    Raises:
        ValueError: If the file cannot be downloaded or saved to the specified path.
    """
    logger.info(
        f"downloading file {path_to_file} from upload {upload_id} on {'prod' if use_prod else 'test'} server"
    )
    response = get_nomad_request(
        f"/uploads/{upload_id}/raw/{path_to_file}",
        use_prod=use_prod,
        with_authentication=with_authentication,
        timeout_in_sec=timeout_in_sec,
        return_json=return_json,
    )
    if response:
        return response
    else:
        raise ValueError(f"could not download file from upload {upload_id}")


def upload_file_to_specified_path(
    upload_id: str,
    remote_path: str,
    local_file: str,
    use_prod: bool = False,
    timeout_in_sec: int = 10,
) -> dict:
    """
    Uploads a file to a specific path within an upload in the NOMAD system.

    This function uploads a specified file to a specific path within an existing upload in the NOMAD system. It allows
    the user to choose between the production and test environments. It also supports setting a custom timeout for the
    upload process to prevent indefinite waiting periods. The function logs the outcome of the upload process, including
    success or failure messages.

    Args:
        upload_id (str): The unique identifier of the upload to which the file should be uploaded.
        path (str): The path where the file should be saved within the upload.
        file: The file to the file to be uploaded.
        use_prod (bool, optional): Flag indicating whether to use the production environment of the NOMAD system.
                                   Defaults to False, which means the test environment is used by default.
        timeout_in_sec (int, optional): The maximum time in seconds to wait for the upload to complete.
                                        Defaults to 10 seconds.

    Returns:
        dict: A dictionary containing the response from the NOMAD system regarding the upload action.
    """
    logger.info(
        f"uploading file {local_file} to upload {upload_id} on {'prod' if use_prod else 'test'} server"
    )
    remote_file = f"{remote_path}{basename(local_file)}"
    response = put_nomad_request(
        f"/uploads/{upload_id}/raw/{remote_path}?file_name={remote_file}",
        with_authentication=True,
        file=local_file,
        remote_path_to_file=remote_file,
        use_prod=use_prod,
        timeout_in_sec=timeout_in_sec,
    )
    upload_class_schema = class_schema(NomadUpload, base_schema=NomadUploadSchema)
    return upload_class_schema().load({**response["data"], "use_prod": use_prod})


def delete_file_to_specified_path(
    upload_id: str,
    remote_path: str,
    local_file: str,
    use_prod: bool = False,
    timeout_in_sec: int = 10,
) -> dict:
    """
    Deletes a file from a specific path within an upload in the NOMAD system.

    This function deletes a specified file from a specific path within an existing upload in the NOMAD system. It allows
    the user to choose between the production and test environments. It also supports setting a custom timeout for the
    deletion process to prevent indefinite waiting periods. The function logs the outcome of the deletion process,
    including success or failure messages.

    Args:
        upload_id (str): The unique identifier of the upload from which the file should be deleted.
        path (str): The path where the file should be deleted within the upload.
        file: The file to the file to be deleted.
        use_prod (bool, optional): Flag indicating whether to use the production environment of the NOMAD system.
                                   Defaults to False, which means the test environment is used by default.
        timeout_in_sec (int, optional): The maximum time in seconds to wait for the deletion to complete.
                                        Defaults to 10 seconds.

    Returns:
        dict: A dictionary containing the response from the NOMAD system regarding the deletion action.
    """
    logger.info(
        f"deleting file {local_file} from upload {upload_id} on {'prod' if use_prod else 'test'} server"
    )
    remote_file = f"{remote_path}{basename(local_file)}"
    response = delete_nomad_request(
        f"/uploads/{upload_id}/raw/{remote_file}",
        with_authentication=True,
        use_prod=use_prod,
        timeout_in_sec=timeout_in_sec,
    )
    upload_class_schema = class_schema(NomadUpload, base_schema=NomadUploadSchema)
    return upload_class_schema().load({**response["data"], "use_prod": use_prod})


def _get_raw_data_of_upload_by_id(
    upload_id: str,
    use_prod: bool = False,
    timeout_in_sec: int = 10,
    with_authentication: bool = False,
) -> ByteString:
    logger.info(
        f"retrieving raw data for upload {upload_id} on {'prod' if use_prod else 'test'} server"
    )
    response = get_nomad_request(
        f"/uploads/{upload_id}/raw",
        use_prod=use_prod,
        with_authentication=with_authentication,
        timeout_in_sec=timeout_in_sec,
        return_json=False,
    )
    if response:
        return response
    else:
        raise ValueError(f"could not retrieve raw data for upload {upload_id}")

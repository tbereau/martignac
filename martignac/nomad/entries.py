import datetime as dt
import json
import logging
import os
import shutil
from collections.abc import ByteString
from dataclasses import field
from pathlib import Path
from typing import Any, Optional, cast
from zipfile import ZipFile

from cachetools.func import ttl_cache
from marshmallow import Schema, pre_load
from marshmallow_dataclass import class_schema, dataclass
from signac.job import Job

from martignac import config
from martignac.nomad.datasets import NomadDataset
from martignac.nomad.uploads import _get_raw_data_of_upload_by_id, get_all_my_uploads
from martignac.nomad.users import NomadUser, get_user_by_id
from martignac.nomad.utils import (
    get_nomad_base_url,
    get_nomad_request,
    post_nomad_request,
)
from martignac.utils.martini_flow_projects import MartiniFlowProject
from martignac.utils.misc import update_nested_dict

logger = logging.getLogger(__name__)

DEFAULT_DATABASE = config()["nomad"]["dataset"]["id"].get(str)
DEFAULT_USE_PROD = config()["nomad"]["use_prod"].get(bool)


@dataclass(frozen=True)
class NomadSectionDefinition:
    used_directly: bool
    definition_id: str
    definition_qualified_name: str


class NomadEntrySchema(Schema):
    @pre_load
    def convert_users(self, data, **kwargs):
        data["main_author"] = get_user_by_id(
            user_id=data["main_author"]["user_id"]
        ).as_dict()
        data["writers"] = [
            get_user_by_id(user_id=w["user_id"]).as_dict() for w in data["writers"]
        ]
        data["authors"] = [
            get_user_by_id(user_id=a["user_id"]).as_dict() for a in data["authors"]
        ]
        data["viewers"] = [
            get_user_by_id(user_id=v["user_id"]).as_dict() for v in data["viewers"]
        ]
        return data


@dataclass(frozen=True)
class NomadEntry:
    """
    Represents an entry in the NOMAD system.

    This class encapsulates the data and metadata associated with an entry in the NOMAD database, including references,
    quantities, datasets, and processing information. It provides properties to access computed attributes such as
    URLs and job-related information derived from comments.

    Attributes:
        entry_id (str): Unique identifier for the entry.
        upload_id (str): Identifier of the upload this entry belongs to.
        references (list[str]): List of reference identifiers associated with this entry.
        origin (str): The origin of the entry data.
        quantities (list[str]): List of quantities measured or calculated for this entry.
        datasets (list[NomadDataset]): List of datasets associated with this entry.
        n_quantities (int): Number of quantities associated with this entry.
        nomad_version (str): Version of the NOMAD software used.
        upload_create_time (dt.datetime): Creation time of the upload this entry is part of.
        nomad_commit (str): Specific commit of the NOMAD software used.
        section_defs (list[NomadSectionDefinition]): Definitions of sections used in this entry.
        processing_errors (list[Any]): List of errors encountered during processing of this entry.
        last_processing_time (dt.datetime): Timestamp of the last processing attempt.
        parser_name (str): Name of the parser used to process this entry.
        calc_id (str): Calculation identifier associated with this entry.
        published (bool): Flag indicating whether this entry is published.
        writers (list[NomadUser]): List of users with write access to this entry.
        sections (list[str]): List of section identifiers associated with this entry.
        processed (bool): Flag indicating whether this entry has been processed.
        mainfile (str): Name of the main file for this entry.
        main_author (NomadUser): The main author of this entry.
        viewers (list[NomadUser]): List of users with view access to this entry.
        entry_create_time (dt.datetime): Creation time of this entry.
        with_embargo (bool): Flag indicating whether this entry is under embargo.
        files (list[str]): List of file names associated with this entry.
        authors (list[NomadUser]): List of authors associated with this entry.
        license (str): License under which this entry is published.
        results (Optional[dict]): Results associated with this entry, if any.
        entry_name (Optional[str]): Name of this entry, if any.
        entry_type (Optional[str]): Type of this entry, if any.
        domain (Optional[str]): Domain of this entry, if any.
        optimade (Optional[dict]): OPTIMADE data associated with this entry, if any.
        comment (Optional[str]): Additional comments associated with this entry, if any.
        upload_name (Optional[str]): Name of the upload this entry belongs to, if any.
        viewer_groups (Optional[list[Any]]): Groups of users with view access, if any.
        writer_groups (Optional[list[Any]]): Groups of users with write access, if any.
        text_search_contents (Optional[list[str]]): Contents available for text search, if any.
        publish_time (Optional[dt.datetime]): Time when this entry was published, if any.
        entry_references (Optional[list[dict]]): References associated with this entry, if any.
        use_prod (Optional[bool]): Flag indicating whether this entry is from the production environment.
    """

    entry_id: str
    upload_id: str
    references: list[str]
    origin: str
    quantities: list[str] = field(repr=False)
    datasets: list[NomadDataset] = field(repr=False)
    n_quantities: int
    nomad_version: str
    upload_create_time: dt.datetime
    nomad_commit: str
    section_defs: list[NomadSectionDefinition] = field(repr=False)
    processing_errors: list[Any]
    last_processing_time: dt.datetime
    parser_name: str
    calc_id: str
    published: bool
    writers: list[NomadUser]
    sections: list[str] = field(repr=False)
    processed: bool
    mainfile: str
    main_author: NomadUser
    viewers: list[NomadUser] = field(repr=False)
    entry_create_time: dt.datetime
    with_embargo: bool
    files: list[str] = field(repr=False)
    authors: list[NomadUser] = field(repr=False)
    license: str
    results: Optional[dict] = field(repr=False, default=None)
    entry_name: Optional[str] = None
    entry_type: Optional[str] = None
    domain: Optional[str] = None
    optimade: Optional[dict] = field(repr=False, default=None)
    comment: Optional[str] = None
    upload_name: Optional[str] = None
    viewer_groups: Optional[list[Any]] = field(repr=False, default=None)
    writer_groups: Optional[list[Any]] = field(repr=False, default=None)
    text_search_contents: Optional[list[str]] = None
    publish_time: Optional[dt.datetime] = None
    entry_references: Optional[list[dict]] = None
    use_prod: Optional[bool] = None

    @property
    def base_url(self) -> Optional[str]:
        if self.use_prod is not None:
            return get_nomad_base_url(self.use_prod)
        return None

    @property
    def nomad_gui_url(self) -> str:
        if self.base_url is None:
            raise ValueError(f"missing attribute 'use_prod' for entry {self}")
        return f"{self.base_url}/gui/user/uploads/upload/id/{self.upload_id}/entry/id/{self.entry_id}"

    @property
    def job_id(self) -> Optional[str]:
        return self._comment_dict.get("job_id", None)

    @property
    def workflow_name(self) -> Optional[str]:
        return self._comment_dict.get("workflow_name", None)

    @property
    def state_point(self) -> dict:
        return self._comment_dict.get("state_point", {})

    @property
    def mdp_files(self) -> Optional[str]:
        return self._comment_dict.get("mdp_files", None)

    @property
    def _comment_dict(self) -> dict:
        return json.loads(self.comment or "{}")


@ttl_cache(maxsize=128, ttl=180)
def get_entry_by_id(
    entry_id: str,
    use_prod: bool = True,
    with_authentication: bool = False,
    timeout_in_sec: int = 10,
) -> NomadEntry:
    """
    Retrieves a NomadEntry object by its unique entry ID.

    This function sends a GET request to the NOMAD system to retrieve the details of a specific entry identified by its
    unique ID. It allows specifying whether to use the production or test environment, whether authentication is
    required for the request, and a timeout for the request.

    Args:
        entry_id (str): The unique identifier of the entry to retrieve.
        use_prod (bool, optional): Flag indicating whether to use the production environment. Defaults to True.
        with_authentication (bool, optional): Flag indicating whether the request should include authentication.
                                               Defaults to False.
        timeout_in_sec (int, optional): The maximum time in seconds to wait for a response from the server.
                                        Defaults to 10 seconds.

    Returns:
        NomadEntry: The NomadEntry object corresponding to the provided entry ID.

    Raises:
        HTTPError: If the request fails or the NOMAD system returns an error response.
    """
    logger.info(
        f"retrieving entry {entry_id} on {'prod' if use_prod else 'test'} server"
    )
    response = get_nomad_request(
        f"/entries/{entry_id}",
        with_authentication=with_authentication,
        use_prod=use_prod,
        timeout_in_sec=timeout_in_sec,
    )
    nomad_entry_schema = class_schema(NomadEntry, base_schema=NomadEntrySchema)
    return nomad_entry_schema().load({**response["data"], "use_prod": use_prod})


@ttl_cache(maxsize=128, ttl=180)
def get_entries_of_upload(
    upload_id: str,
    use_prod: bool = False,
    with_authentication: bool = False,
    timeout_in_sec: int = 10,
) -> list[NomadEntry]:
    """
    Retrieves a list of NomadEntry objects associated with a specific upload ID.

    This function sends a GET request to the NOMAD system to retrieve all entries associated with a given upload ID.
    It allows specifying whether to use the production or test environment, whether authentication is required for
    the request, and a timeout for the request.

    Args:
        upload_id (str): The unique identifier of the upload whose entries are to be retrieved.
        use_prod (bool, optional): Flag indicating whether to use the production environment. Defaults to False.
        with_authentication (bool, optional): Flag indicating whether the request should include authentication.
                                               Defaults to False.
        timeout_in_sec (int, optional): The maximum time in seconds to wait for a response from the server.
                                        Defaults to 10 seconds.

    Returns:
        list[NomadEntry]: A list of NomadEntry objects corresponding to the entries of the specified upload ID.

    Raises:
        HTTPError: If the request fails or the NOMAD system returns an error response.
    """
    logger.info(
        f"retrieving entries for upload {upload_id} on {'prod' if use_prod else 'test'} server"
    )
    response = get_nomad_request(
        f"/uploads/{upload_id}/entries",
        with_authentication=with_authentication,
        use_prod=use_prod,
        timeout_in_sec=timeout_in_sec,
    )
    nomad_entry_schema = class_schema(NomadEntry, base_schema=NomadEntrySchema)
    return [
        nomad_entry_schema().load({**r["entry_metadata"], "use_prod": use_prod})
        for r in response["data"]
    ]


@ttl_cache(maxsize=128, ttl=180)
def get_entries_of_my_uploads(
    use_prod: bool = False, timeout_in_sec: int = 10
) -> list[NomadEntry]:
    """
    Retrieves a list of NomadEntry objects associated with the uploads owned by the current user.

    This function iterates over all uploads owned by the current user, identified through their unique upload IDs,
    and aggregates all entries associated with these uploads. It allows specifying whether to use the production
    or test environment for the NOMAD system, and a timeout for the request.

    Args:
        use_prod (bool, optional): Flag indicating whether to use the production environment. Defaults to False,
                                   indicating that the test environment is used by default.
        timeout_in_sec (int, optional): The maximum time in seconds to wait for a response from the server.
                                        Defaults to 10 seconds.

    Returns:
        list[NomadEntry]: A list of NomadEntry objects corresponding to the entries associated with the uploads
                          owned by the current user.
    """
    return [
        upload_entry
        for u in get_all_my_uploads(use_prod=use_prod, timeout_in_sec=timeout_in_sec)
        for upload_entry in get_entries_of_upload(
            u.upload_id, with_authentication=True, use_prod=use_prod
        )
    ]


@ttl_cache(maxsize=128, ttl=180)
def get_entries_in_database(
    database_id: str = DEFAULT_DATABASE, use_prod: bool = DEFAULT_USE_PROD
) -> list[NomadEntry]:
    """
    Retrieves a list of NomadEntry objects from a specified database.

    This function queries the NOMAD system for entries within a specified database, using the provided database ID.
    It allows specifying whether to use the production or test environment for the query. The function leverages
    the `query_entries` function to perform the actual query based on the provided parameters.

    Args:
        database_id (str, optional): The unique identifier of the database from which to retrieve entries.
                                     Defaults to the value of `DEFAULT_DATABASE` from the configuration.
        use_prod (bool, optional): Flag indicating whether to use the production environment. Defaults to the value
                                   of `DEFAULT_USE_PROD` from the configuration.

    Returns:
        list[NomadEntry]: A list of NomadEntry objects corresponding to the entries found in the specified database.
    """
    return query_entries(dataset_id=database_id, use_prod=use_prod)


@ttl_cache(maxsize=128, ttl=180)
def query_entries(
    worfklow_name: Optional[str] = None,
    program_name: Optional[str] = None,
    dataset_id: Optional[str] = None,
    origin: Optional[str] = None,
    page_size: int = 10,
    max_entries: int = 50,
    use_prod: bool = True,
) -> list[NomadEntry]:
    """
    Queries the NOMAD system for entries based on various filters and returns a list of NomadEntry objects.

    This function constructs a query to the NOMAD system, allowing for filtering based on workflow name, program name,
    dataset ID, and origin. It supports pagination through `page_size` and limits the number of entries returned with
    `max_entries`. The environment (production or test) can be specified with `use_prod`.

    Args:
        worfklow_name (str, optional): Filter entries by the name of the workflow. Defaults to None.
        program_name (str, optional): Filter entries by the program name. Defaults to None.
        dataset_id (str, optional): Filter entries by the dataset ID. Defaults to None.
        origin (str, optional): Filter entries by their origin. Defaults to None.
        page_size (int, optional): Number of entries to return per page. Defaults to 10.
        max_entries (int, optional): Maximum number of entries to return. Defaults to 50.
        use_prod (bool, optional): Flag indicating whether to query the production environment. Defaults to True.

    Returns:
        list[NomadEntry]: A list of NomadEntry objects that match the query criteria.
    """
    json_dict = {
        "query": {},
        "pagination": {"page_size": page_size},
        "required": {"include": ["entry_id"]},
    }
    entries = []
    while (max_entries > 0 and len(entries) <= max_entries) or (max_entries < 0):
        if dataset_id:
            json_dict["query"]["datasets"] = {"dataset_id": dataset_id}
        if worfklow_name:
            json_dict["query"]["results.method"] = {"workflow_name": worfklow_name}
        if program_name:
            json_dict["query"]["results.method"] = {
                "simulation": {"program_name": program_name}
            }
        if origin:
            json_dict["query"]["origin"] = origin
        query = post_nomad_request(
            "/entries/query", json_dict=json_dict, use_prod=use_prod
        )
        entries.extend([q["entry_id"] for q in query["data"]])
        next_page_after_value = query["pagination"].get("next_page_after_value", None)
        if next_page_after_value:
            json_dict["pagination"]["page_after_value"] = next_page_after_value
        else:
            break
    if max_entries > 0:
        entries = entries[:max_entries]
    return [get_entry_by_id(e, use_prod=use_prod) for e in entries]


@ttl_cache(maxsize=128, ttl=180)
def download_raw_data_of_job(job: Job, timeout_in_sec: int = 10) -> bool:
    """
    Downloads the raw data associated with a given job from the NOMAD system and stores it in the job's directory.

    This function attempts to find NOMAD entries corresponding to the specified job, downloads the raw data for the
    first matching entry, and extracts it into the job's directory. It handles the creation of a temporary ZIP file
    for the raw data, extracts the relevant files, and cleans up the temporary files. Additionally, it updates the
    job's document with NOMAD dataset and upload IDs, and specific workflow information if available.

    Args:
        job (Job): The signac job object for which to download the raw data.
        timeout_in_sec (int, optional): The maximum time in seconds to wait for a response from the server when
                                        downloading the raw data. Defaults to 10 seconds.

    Returns:
        bool: True if the raw data was successfully downloaded and processed, False if no corresponding NOMAD entries
              were found for the job.

    Raises:
        TypeError: If the job's project does not derive from MartiniFlowProject.
        ValueError: If the found entries have inconsistent upload IDs, indicating a data retrieval or processing error.
    """
    entries = find_entries_corresponding_to_job(job)
    if len(entries) == 0:
        return False
    entry = entries[0]
    logger.info(
        f"found {'un' if not entry.published else ''}published entry {entry.upload_id}"
    )
    if entry.published:
        zip_content = _get_raw_data_of_upload_by_id(
            entry.upload_id,
            use_prod=entry.use_prod,
            timeout_in_sec=timeout_in_sec,
            with_authentication=not entry.published,
        )
    else:
        zip_content = _get_raw_data_of_entry_by_id(
            entry.entry_id,
            use_prod=entry.use_prod,
            timeout_in_sec=timeout_in_sec,
            with_authentication=not entry.published,
        )
    zip_file_name = "nomad_archive.zip"
    with open(job.fn(zip_file_name), "wb") as f:
        f.write(bytes(zip_content))
    zip_file = ZipFile(job.fn(zip_file_name))
    archive_path = job.path + "/" + entry.upload_id if not entry.published else job.path
    name_list = zip_file.namelist()
    logger.info(f"zip content: {name_list}")
    if entry.published:
        name_list = [
            name for name in name_list if name.startswith(f"{entry.upload_id}/")
        ]
    zip_file.extractall(path=job.path, members=name_list)
    for file_name in name_list:
        if Path(file_name).name not in os.listdir(job.path):
            shutil.move(job.path + "/" + file_name, job.path)
    os.remove(job.fn(zip_file_name))
    if "signac_job_document.json" in os.listdir(archive_path):
        with open(archive_path + "/signac_job_document.json") as fp:
            json_data = json.load(fp)
            job.doc = update_nested_dict(job.doc, dict(json_data))
    if not entry.published:
        for file_name in os.listdir(archive_path):
            os.remove(archive_path + "/" + file_name)
        os.removedirs(archive_path)
    job.document["nomad_dataset_id"] = MartiniFlowProject.nomad_dataset_id
    if entry.workflow_name not in job.document:
        job.document[entry.workflow_name] = {}
    job.document[entry.workflow_name]["nomad_upload_id"] = entry.upload_id
    return True


@ttl_cache(maxsize=128, ttl=180)
def find_entries_corresponding_to_job(job: Job) -> list[NomadEntry]:
    """
    Finds NOMAD entries that correspond to a given job.

    This function searches for NOMAD entries that are associated with the specified job. It checks if the job's project
    is a subclass of MartiniFlowProject. If not, it raises a TypeError. The function then proceeds to match entries
    based on the job's ID and the hash of MDP files associated with the job. It combines entries found through querying
    the NOMAD dataset with entries associated with uploads owned by the current user. If entries with inconsistent
    upload IDs are found, a ValueError is raised.

    Args:
        job (Job): The signac job object to find corresponding NOMAD entries for.

    Returns:
        list[NomadEntry]: A list of NomadEntry objects that correspond to the given job.

    Raises:
        TypeError: If the job's project does not derive from MartiniFlowProject.
        ValueError: If found entries have inconsistent upload IDs, indicating a data retrieval or processing error.
    """
    if not issubclass(type(job.project), MartiniFlowProject):
        raise TypeError(
            f"job project {type(job.project)} does not derive from MartiniFlowProject"
        )
    project = cast("MartiniFlowProject", job.project)

    def associate_entry_to_job(entry_: NomadEntry) -> Optional[NomadEntry]:
        if (
            entry_.comment is not None
            and entry_.job_id == job.id
            and entry_.mdp_files
            == project.get_hash_for_files(job, list(project.mdp_files.values()))
        ):
            return entry_

    match_entries = []
    nomad_entries = query_entries(
        dataset_id=project.nomad_dataset_id, use_prod=project.nomad_use_prod_database
    )
    for entry in nomad_entries:
        if found_entry := associate_entry_to_job(entry):
            match_entries.append(found_entry)
    my_uploads = get_all_my_uploads(use_prod=project.nomad_use_prod_database)
    for upload in my_uploads:
        upload_entries = get_entries_of_upload(
            upload.upload_id,
            with_authentication=True,
            use_prod=project.nomad_use_prod_database,
        )
        for entry in upload_entries:
            if found_entry := associate_entry_to_job(entry):
                match_entries.append(found_entry)
    if len(match_entries) > 0 and not all(
        entry.upload_id == match_entries[0].upload_id for entry in match_entries
    ):
        raise ValueError(f"Inconsistent upload IDs in entries:\n{match_entries}")
    return match_entries


def _get_raw_data_of_entry_by_id(
    entry_id: str,
    use_prod: bool = False,
    timeout_in_sec: int = 10,
    with_authentication: bool = False,
) -> ByteString:
    logger.info(
        f"retrieving raw data of entry ID {entry_id} on {'prod' if use_prod else 'test'} server"
    )
    response = get_nomad_request(
        f"/entries/{entry_id}/raw?compress=true",
        with_authentication=with_authentication,
        use_prod=use_prod,
        timeout_in_sec=timeout_in_sec,
        return_json=False,
        accept_field="application/zip",
    )
    return response

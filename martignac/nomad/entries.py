import datetime as dt
import json
import logging
import os
import shutil
from dataclasses import field
from pathlib import Path
from typing import Any, ByteString, Optional, cast
from zipfile import ZipFile

from cachetools.func import ttl_cache
from marshmallow import Schema, pre_load
from marshmallow_dataclass import class_schema, dataclass
from signac.job import Job

from martignac import config
from martignac.nomad.datasets import NomadDataset
from martignac.nomad.uploads import get_all_my_uploads
from martignac.nomad.users import NomadUser, get_user_by_id
from martignac.nomad.utils import get_nomad_base_url, get_nomad_request, post_nomad_request
from martignac.utils.martini_flow_projects import MartiniFlowProject

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
        data["main_author"] = get_user_by_id(user_id=data["main_author"]["user_id"]).as_dict()
        data["writers"] = [get_user_by_id(user_id=w["user_id"]).as_dict() for w in data["writers"]]
        data["authors"] = [get_user_by_id(user_id=a["user_id"]).as_dict() for a in data["authors"]]
        data["viewers"] = [get_user_by_id(user_id=v["user_id"]).as_dict() for v in data["viewers"]]
        return data


@dataclass(frozen=True)
class NomadEntry:
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
    results: dict = field(repr=False)
    entry_name: str
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
    entry_type: str
    authors: list[NomadUser] = field(repr=False)
    license: str
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
    entry_id: str, use_prod: bool = True, with_authentication: bool = False, timeout_in_sec: int = 10
) -> NomadEntry:
    logger.info(f"retrieving entry {entry_id} on {'prod' if use_prod else 'test'} server")
    response = get_nomad_request(
        f"/entries/{entry_id}",
        with_authentication=with_authentication,
        use_prod=use_prod,
        timeout_in_sec=timeout_in_sec,
    )
    nomad_entry_schema = class_schema(NomadEntry, base_schema=NomadEntrySchema)
    return nomad_entry_schema().load({**response["data"], "use_prod": use_prod})


@ttl_cache(maxsize=128, ttl=180)
def get_entries_of_upload(upload_id: str, use_prod: bool = False, timeout_in_sec: int = 10) -> list[NomadEntry]:
    logger.info(f"retrieving entries for upload {upload_id} on {'prod' if use_prod else 'test'} server")
    response = get_nomad_request(
        f"/uploads/{upload_id}/entries",
        with_authentication=True,
        timeout_in_sec=timeout_in_sec,
    )
    nomad_entry_schema = class_schema(NomadEntry, base_schema=NomadEntrySchema)
    return [nomad_entry_schema().load({**r["entry_metadata"], "use_prod": use_prod}) for r in response["data"]]


def get_entries_of_my_uploads(use_prod: bool = False, timeout_in_sec: int = 10) -> list[NomadEntry]:
    return [
        upload_entry
        for u in get_all_my_uploads(use_prod=use_prod, timeout_in_sec=timeout_in_sec)
        for upload_entry in get_entries_of_upload(u.upload_id)
    ]


def get_entries_in_database(database_id: str = DEFAULT_DATABASE, use_prod: bool = DEFAULT_USE_PROD) -> list[NomadEntry]:
    return query_entries(dataset_id=database_id, use_prod=use_prod)


@ttl_cache(maxsize=128, ttl=180)
def query_entries(
    worfklow_name: str = None,
    program_name: str = None,
    dataset_id: str = None,
    origin: str = None,
    page_size: int = 10,
    max_entries: int = 50,
    use_prod: bool = True,
) -> list[NomadEntry]:
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
            json_dict["query"]["results.method"] = {"simulation": {"program_name": program_name}}
        if origin:
            json_dict["query"]["origin"] = origin
        query = post_nomad_request("/entries/query", json_dict=json_dict, use_prod=use_prod)
        entries.extend([q["entry_id"] for q in query["data"]])
        next_page_after_value = query["pagination"].get("next_page_after_value", None)
        if next_page_after_value:
            json_dict["pagination"]["page_after_value"] = next_page_after_value
        else:
            break
    if max_entries > 0:
        entries = entries[:max_entries]
    return [get_entry_by_id(e, use_prod=use_prod) for e in entries]


def download_raw_data_of_job(job: Job, timeout_in_sec: int = 10) -> bool:
    entries = find_entries_corresponding_to_job(job)
    if len(entries) == 0:
        return False
    entry = entries[0]
    zip_content = _get_raw_data_of_entry_by_id(
        entry.entry_id, use_prod=entry.use_prod, timeout_in_sec=timeout_in_sec, with_authentication=not entry.published
    )
    zip_file_name = "nomad_archive.zip"
    with open(job.fn(zip_file_name), "wb") as f:
        f.write(bytes(zip_content))
    zip_file = ZipFile(job.fn(zip_file_name))
    archive_path = job.path + "/" + entry.upload_id
    name_list = zip_file.namelist()
    name_list = [name for name in name_list if name.startswith(f"{entry.upload_id}/")]
    zip_file.extractall(path=job.path, members=name_list)
    for file_name in name_list:
        if Path(file_name).name not in os.listdir(job.path):
            shutil.move(job.path + "/" + file_name, job.path)
    os.remove(job.fn(zip_file_name))
    if "signac_job_document.json" in os.listdir(archive_path):
        with open(archive_path + "/signac_job_document.json") as fp:
            json_data = json.load(fp)
            job.document.update(json_data)
    for file_name in os.listdir(archive_path):
        os.remove(archive_path + "/" + file_name)
    os.removedirs(archive_path)
    job.document["nomad_dataset_id"] = MartiniFlowProject.nomad_dataset_id
    if "nomad_upload_id" not in job.document:
        job.document["nomad_upload_id"] = {}
    job.document["nomad_upload_id"][entry.workflow_name] = entry.upload_id
    return True


def find_entries_corresponding_to_job(job: Job) -> list[NomadEntry]:
    if not issubclass(type(job.project), MartiniFlowProject):
        raise TypeError(f"job project {type(job.project)} does not derive from MartiniFlowProject")
    project = cast("MartiniFlowProject", job.project)

    def associate_entry_to_job(entry_: NomadEntry) -> Optional[NomadEntry]:
        if (
            entry_.comment is not None
            and entry_.job_id == job.id
            and entry_.mdp_files == project.get_hash_for_files(job, list(project.mdp_files.values()))
        ):
            return entry_

    match_entries = []
    nomad_entries = query_entries(dataset_id=project.nomad_dataset_id, use_prod=project.nomad_use_prod_database)
    for entry in nomad_entries:
        if found_entry := associate_entry_to_job(entry):
            match_entries.append(found_entry)
    my_uploads = get_all_my_uploads(use_prod=project.nomad_use_prod_database)
    for upload in my_uploads:
        upload_entries = get_entries_of_upload(upload.upload_id)
        for entry in upload_entries:
            if found_entry := associate_entry_to_job(entry):
                match_entries.append(found_entry)
    if len(match_entries) > 0 and not all([entry.upload_id == match_entries[0].upload_id for entry in match_entries]):
        raise ValueError(f"Inconsistent upload IDs in entries:\n{match_entries}")
    return match_entries


def _get_raw_data_of_entry_by_id(
    entry_id: str, use_prod: bool = False, timeout_in_sec: int = 10, with_authentication: bool = False
) -> ByteString:
    logger.info(f"retrieving raw data of entry ID {entry_id} on {'prod' if use_prod else 'test'} server")
    response = get_nomad_request(
        f"/entries/{entry_id}/raw?compress=true",
        with_authentication=with_authentication,
        use_prod=use_prod,
        timeout_in_sec=timeout_in_sec,
        return_json=False,
        accept_field="application/zip",
    )
    return response

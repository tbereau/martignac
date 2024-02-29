import datetime as dt
import logging
from dataclasses import field
from typing import Any, Optional

from marshmallow import Schema, pre_load
from marshmallow_dataclass import class_schema, dataclass

from martignac.nomad.datasets import NomadDataset
from martignac.nomad.users import NomadUser, get_user_by_id
from martignac.nomad.utils import get_nomad_request, post_nomad_request

logger = logging.getLogger(__name__)


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


def get_entry_by_id(entry_id: str, use_prod: bool = True, with_authentication: bool = False) -> NomadEntry:
    logger.info(f"retrieving entry {entry_id} on {'prod' if use_prod else 'test'} server")
    response = get_nomad_request(
        f"/entries/{entry_id}",
        with_authentication=with_authentication,
        use_prod=use_prod,
    )
    nomad_entry_schema = class_schema(NomadEntry, base_schema=NomadEntrySchema)
    return nomad_entry_schema().load(response["data"])


def get_entries_of_upload(upload_id: str, use_prod: bool = False, timeout_in_sec: int = 10) -> list[NomadEntry]:
    logger.info(f"retrieving entries for upload {upload_id} on {'prod' if use_prod else 'test'} server")
    response = get_nomad_request(
        f"/uploads/{upload_id}/entries",
        with_authentication=True,
        timeout_in_sec=timeout_in_sec,
    )
    nomad_entry_schema = class_schema(NomadEntry, base_schema=NomadEntrySchema)
    return [nomad_entry_schema().load(r["entry_metadata"]) for r in response["data"]]


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

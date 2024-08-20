import json
import logging
from typing import Optional, cast

from cachetools.func import ttl_cache
from marshmallow import Schema, pre_load
from marshmallow_dataclass import class_schema, dataclass
from signac.job import Job

from martignac import config
from martignac.nomad.utils import post_nomad_request
from martignac.utils.martini_flow_projects import MartiniFlowProject

logger = logging.getLogger(__name__)


class NomadMiniQuerySchema(Schema):
    @pre_load
    def convert_comments(self, data, **kwargs):
        try:
            data["comment"] = json.loads(data.get("comment"))
            data["workflow_name"] = data["comment"].get("workflow_name")
        except (json.JSONDecodeError, TypeError):
            data["comment"] = {}
        return data


@dataclass
class NomadMiniQuery:
    entry_id: str
    upload_id: str
    comment: dict
    published: bool
    workflow_name: Optional[str] = None
    use_prod: Optional[bool] = None

    @property
    def is_martignac_entry(self) -> bool:
        return self.comment != {} and "mdp_files" in self.comment

    def matches_with_job(self, job: Job) -> bool:
        if not self.is_martignac_entry:
            return False
        project = cast("MartiniFlowProject", job.project)
        return (
            self.is_martignac_entry
            and self.comment.get("mdp_files")
            == project.get_hash_for_files(job, list(project.mdp_files.values()))
            and self.comment.get("job_id") == job.id
        )


@ttl_cache(maxsize=128, ttl=180)
def query_entries_in_mini_format(
    use_prod: bool = False,
    with_authentication: bool = False,
    page_size: int = 10,
    timeout_in_sec: int = 10,
) -> list[NomadMiniQuery]:
    logger.info(f"querying comments on {'prod' if use_prod else 'test'} server")
    nomad_query_authors: list[str] = [
        c.get(str) for c in config()["nomad"]["query_authors"]
    ]
    query = {
        "query": {
            "authors.name": {
                "any": nomad_query_authors,
            },
        },
        "pagination": {
            "page_size": page_size,
        },
        "required": {
            "include": [
                "entry_id",
                "upload_id",
                "comment",
                "published",
            ],
        },
    }
    raw_entries = []
    while True:
        response = post_nomad_request(
            "entries/query",
            use_prod=use_prod,
            data=json.dumps(query),
            with_authentication=with_authentication,
            timeout_in_sec=timeout_in_sec,
        )
        raw_entries.extend(response["data"])
        next_value = response["pagination"].get("next_page_after_value")
        if not next_value:
            break
        query["pagination"]["page_after_value"] = next_value
    query_class_schema = class_schema(NomadMiniQuery, base_schema=NomadMiniQuerySchema)
    processed_entries: list[NomadMiniQuery] = [
        query_class_schema().load({**entry, "use_prod": use_prod})
        for entry in raw_entries
    ]
    processed_entries = [
        entry for entry in processed_entries if entry.is_martignac_entry
    ]
    logger.info(f"found {len(processed_entries)} Martignac entries online")
    return processed_entries


def find_mini_queries_corresponding_to_job(job: Job) -> list[NomadMiniQuery]:
    if not issubclass(type(job.project), MartiniFlowProject):
        raise TypeError(
            f"job project {type(job.project)} does not derive from MartiniFlowProject"
        )
    found_entries = [
        query for query in query_entries_in_mini_format() if query.matches_with_job(job)
    ]
    logger.info(f"found {len(found_entries)} entries online for job {job.id}")
    return found_entries

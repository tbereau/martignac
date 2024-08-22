import json
import logging
from typing import Optional, cast

from cachetools.func import ttl_cache
from marshmallow import Schema, pre_load
from marshmallow_dataclass import class_schema, dataclass
from signac.job import Job

from martignac import config
from martignac.nomad.datasets import NomadDataset
from martignac.nomad.queries import query_nomad_entries
from martignac.utils.martini_flow_projects import MartiniFlowProject, MartiniTypeFlow

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
    datasets: Optional[list[NomadDataset]] = None
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

    def matches_with_workflow(self, workflow_name: str) -> bool:
        if not self.is_martignac_entry:
            return False
        return self.is_martignac_entry and self.workflow_name == workflow_name


@ttl_cache(maxsize=128, ttl=180)
def query_entries_in_mini_format(
    use_prod: bool = False,
    with_authentication: bool = False,
    page_size: int = 100,
    timeout_in_sec: int = 10,
    only_workflows: bool = True,
    dataset_ids: Optional[tuple[str]] = None,
    query_authors: Optional[tuple[str]] = None,
) -> list[NomadMiniQuery]:
    if query_authors is None:
        try:
            query_authors = tuple(
                str(c.get(str)) for c in config()["nomad"]["query_authors"]
            )
        except KeyError:
            query_authors = ()
    dataset_ids = dataset_ids or ()
    required = ("entry_id", "upload_id", "comment", "published", "datasets.*")
    raw_entries = query_nomad_entries(
        use_prod=use_prod,
        with_authentication=with_authentication,
        page_size=page_size,
        timeout_in_sec=timeout_in_sec,
        only_workflows=only_workflows,
        dataset_ids=dataset_ids,
        query_authors=query_authors,
        required=required,
    )
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
    project = cast("MartiniFlowProject", job.project)
    found_entries = [
        query
        for query in query_entries_in_mini_format(
            use_prod=project.nomad_use_prod_database
        )
        if query.matches_with_job(job)
    ]
    logger.info(f"found {len(found_entries)} entries online for job {job.id}")
    return found_entries


def find_mini_queries_corresponding_to_workflow(
    workflow: MartiniTypeFlow,
    dataset_ids: Optional[list[str]] = None,
    use_prod: Optional[bool] = None,
) -> list[NomadMiniQuery]:
    if dataset_ids is None:
        dataset_ids = []
    if use_prod is None:
        use_prod = workflow.nomad_use_prod_database
    found_entries = [
        query
        for query in query_entries_in_mini_format(
            use_prod=use_prod, dataset_ids=tuple(dataset_ids)
        )
        if query.matches_with_workflow(workflow.class_name())
    ]
    logger.info(
        f"found {len(found_entries)} entries online for workflow {workflow.class_name()}"
    )
    return found_entries

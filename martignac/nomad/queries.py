import json
import logging
from typing import Optional

from cachetools.func import ttl_cache

from martignac.nomad.utils import post_nomad_request

logger = logging.getLogger(__name__)


@ttl_cache(maxsize=128, ttl=180)
def query_nomad_entries(
    use_prod: bool = False,
    with_authentication: bool = False,
    page_size: int = 100,
    timeout_in_sec: int = 10,
    only_workflows: bool = True,
    dataset_ids: Optional[tuple[str, ...]] = None,
    query_authors: Optional[tuple[str, ...]] = None,
    required: Optional[tuple[str, ...]] = None,
    entry_ids: Optional[tuple[str, ...]] = None,
) -> list[dict]:
    dataset_ids = dataset_ids or ()
    query_authors = query_authors or ()
    required = required or ()
    entry_ids = entry_ids or ()
    logger.info(f"querying comments on {'prod' if use_prod else 'test'} server")
    query = {
        "query": {
            "authors.name": {
                "any": query_authors,
            },
            "datasets.dataset_id": {
                "any": dataset_ids,
            },
            "entry_id": {
                "any": entry_ids,
            },
        },
        "pagination": {
            "page_size": page_size,
        },
    }
    if required:
        query["required"] = {"include": required}
    if only_workflows:
        query["query"]["entry_type"] = "Workflow"
    logger.info(f"querying entries with query {query}")
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
    return raw_entries

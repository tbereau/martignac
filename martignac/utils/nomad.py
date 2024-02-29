import json
from typing import Optional


def generate_user_metadata(
    file_name: str,
    comment: Optional[str] = None,
    references: list[str] = None,
    coauthors: list[str] = None,
    datasets: list[str] = None,
) -> None:
    user_metadata = {"comment": comment, "references": references, "coauthors": coauthors, "datasets": datasets}
    with open(file_name, "w") as f:
        json.dump(user_metadata, f)

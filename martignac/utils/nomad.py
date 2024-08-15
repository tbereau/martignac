import json
from typing import Optional


def generate_user_metadata(
    file_name: str,
    comment: Optional[str] = None,
    references: Optional[list[str]] = None,
    coauthors: Optional[list[str]] = None,
    datasets: Optional[list[str]] = None,
) -> None:
    """
    Generates a JSON file containing user metadata for a project.

    This function creates a JSON file with metadata information such as comments, references, coauthors, and datasets
    related to a project. It's useful for documenting project details, collaborations, and data sources.

    Parameters:
        file_name (str): The name of the file where the metadata will be saved.
        comment (Optional[str]): An optional comment or description about the project.
        references (list[str]): An optional list of references or citations related to the project.
        coauthors (list[str]): An optional list of coauthors involved in the project.
        datasets (list[str]): An optional list of datasets used or generated during the project.

    Returns:
        None: This function does not return a value but writes the metadata to a file.
    """
    user_metadata = {
        "comment": comment,
        "references": references,
        "coauthors": coauthors,
        "datasets": datasets,
    }
    with open(file_name, "w") as f:
        json.dump(user_metadata, f)

from abc import abstractmethod
from dataclasses import dataclass
from functools import cached_property
from typing import Optional, TypeVar

from marshmallow import Schema
from marshmallow_dataclass import class_schema

from martignac.nomad.entries import (
    NomadEntry,
    get_entries_of_upload,
    get_specific_file_from_entry,
)
from martignac.nomad.uploads import get_specific_file_from_upload

Interface = TypeVar("Interface", bound="GenericInterface")


@dataclass(frozen=True)
class GenericInterface:
    """
    A base class representing a generic interface for handling NOMAD uploads and workflows.

    Attributes:
        fetched_nomad (bool): Indicates if the NOMAD data has been fetched.
        files_symlinked (bool): Indicates if the files have been symlinked.
        gromacs_logs (dict): A dictionary containing GROMACS log information.
        itp_files (str): A string representing the ITP files.
        mdp_files (str): A string representing the MDP files.
        nomad_upload_id (str): The ID of the NOMAD upload.
        nomad_workflow (str): The workflow associated with the NOMAD upload.
        ready_for_nomad_upload (bool): Indicates if the data is ready for NOMAD upload.
        tasks (dict): A dictionary containing task information.
    """

    fetched_nomad: bool
    files_symlinked: bool
    gromacs_logs: dict
    itp_files: str
    mdp_files: str
    nomad_upload_id: str
    nomad_workflow: str
    ready_for_nomad_upload: bool
    tasks: dict
    upload_id: str
    job_id: str
    use_prod: bool
    with_authentication: bool

    @classmethod
    @abstractmethod
    def from_upload(
        cls,
        upload_id: str,
        job_id: Optional[str] = None,
        use_prod: bool = False,
        with_authentication: bool = True,
    ) -> type["GenericInterface"]:
        """
        Abstract method to create an instance of the interface from a NOMAD upload.

        Args:
            upload_id (str): The ID of the NOMAD upload.
            job_id (Optional[str]): The ID of the job. Defaults to None.
            use_prod (bool): Flag to indicate if production environment should be used. Defaults to False.
            with_authentication (bool): Flag to indicate if authentication is required. Defaults to True.

        Returns:
            type[GenericInterface]: An instance of the GenericInterface class.

        Raises:
            NotImplementedError: If the method is not implemented in a subclass.
        """
        raise NotImplementedError

    @classmethod
    @abstractmethod
    def from_entry(
        cls, entry: NomadEntry, with_authentication: bool = False
    ) -> type["GenericInterface"]:
        raise NotImplementedError

    @cached_property
    def state_point(self) -> dict:
        file_name = "signac_statepoint.json"
        path_to_file = f"{self.job_id}/{file_name}" if self.job_id else file_name
        sp_doc = get_specific_file_from_upload(
            self.upload_id,
            path_to_file,
            use_prod=self.use_prod,
            with_authentication=self.with_authentication,
        )
        return sp_doc


def get_interface_for_entry(
    entry: NomadEntry,
    interface_class: type,
    workflow_class_name: str,
    find_first_job_id: bool = False,
    file_name: str = "signac_job_document.json",
    with_authentication: bool = True,
    base_schema: Optional[type(Schema)] = None,
) -> Interface:
    job_id = None
    if find_first_job_id:
        job_id = entry.job_id
    path_to_file = f"{job_id}/{file_name}" if job_id else file_name
    job_doc = get_specific_file_from_entry(
        entry.entry_id,
        path_to_file,
        use_prod=entry.use_prod,
        with_authentication=with_authentication,
    )
    interface_class_schema = class_schema(interface_class, base_schema=base_schema)
    info_doc = {
        **job_doc[workflow_class_name],
        "upload_id": entry.upload_id,
        "job_id": job_id or "",
        "use_prod": entry.use_prod,
        "with_authentication": with_authentication,
    }
    interface = interface_class_schema().load(info_doc)
    return interface


def get_interface_for_upload_id_and_job_id(
    upload_id: str,
    interface_class: type,
    workflow_class_name: str,
    job_id: Optional[str] = None,
    find_first_job_id: bool = False,
    file_name: str = "signac_job_document.json",
    use_prod: bool = False,
    with_authentication: bool = True,
    base_schema: Optional[type(Schema)] = None,
) -> Interface:
    """
    Retrieve an interface instance for a given upload ID and job ID.

    Args:
        upload_id (str): The ID of the NOMAD upload.
        interface_class (type): The class type of the interface to be instantiated.
        workflow_class_name (str): The name of the workflow class.
        job_id (Optional[str]): The ID of the job. Defaults to None.
        find_first_job_id: bool: Flag to indicate if the first job ID should be used. Defaults to False.
        file_name (str): The name of the file to retrieve. Defaults to "signac_job_document.json".
        use_prod (bool): Flag to indicate if production environment should be used. Defaults to False.
        with_authentication (bool): Flag to indicate if authentication is required. Defaults to True.
        base_schema: Optional[Schema]: The base schema to use for the interface class. Defaults to None.

    Returns:
        Interface: An instance of the specified interface class.

    Raises:
        ValueError: If the job ID cannot be determined from the upload entries.
    """
    if not job_id and find_first_job_id:
        upload_entries = get_entries_of_upload(upload_id, use_prod, with_authentication)
        # Take the first one
        job_id = next(iter({u.job_id for u in upload_entries}))
    path_to_file = f"{job_id}/{file_name}" if job_id else file_name
    job_doc = get_specific_file_from_upload(
        upload_id,
        path_to_file,
        use_prod=use_prod,
        with_authentication=with_authentication,
    )
    interface_class_schema = class_schema(interface_class, base_schema=base_schema)
    info_doc = {
        **job_doc[workflow_class_name],
        "upload_id": upload_id,
        "job_id": job_id or "",
        "use_prod": use_prod,
        "with_authentication": with_authentication,
    }
    interface = interface_class_schema().load(info_doc)
    return interface

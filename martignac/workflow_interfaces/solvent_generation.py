from dataclasses import dataclass
from typing import Optional

from martignac.nomad.entries import NomadEntry
from martignac.workflow_interfaces.generic import (
    GenericInterface,
    get_interface_for_entry,
    get_interface_for_upload_id_and_job_id,
)


@dataclass(frozen=True)
class SolventGenerationInterface(GenericInterface):
    solvent_gro: str
    solvent_name: str
    solvent_top: str

    @classmethod
    def from_upload(
        cls,
        upload_id: str,
        job_id: Optional[str] = None,
        use_prod: bool = False,
        with_authentication: bool = True,
    ) -> "SolventGenerationInterface":
        return get_interface_for_upload_id_and_job_id(
            upload_id,
            cls,
            "SolventGenFlow",
            job_id,
            use_prod=use_prod,
            with_authentication=with_authentication,
            find_first_job_id=False,
        )

    @classmethod
    def from_entry(
        cls, entry: NomadEntry, with_authentication: bool = False
    ) -> "SolventGenerationInterface":
        return get_interface_for_entry(
            entry,
            cls,
            "SolventGenFlow",
            with_authentication=with_authentication,
            find_first_job_id=False,
        )

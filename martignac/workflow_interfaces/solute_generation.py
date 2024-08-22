from dataclasses import dataclass
from typing import Optional

from martignac.nomad.entries import NomadEntry
from martignac.workflow_interfaces.generic import (
    GenericInterface,
    get_interface_for_entry,
    get_interface_for_upload_id_and_job_id,
)


@dataclass(frozen=True)
class SoluteGenerationInterface(GenericInterface):
    solute_has_charged_beads: bool
    solute_gro: str
    solute_itp: str
    solute_name: str
    solute_top: str

    @classmethod
    def from_upload(
        cls,
        upload_id: str,
        job_id: Optional[str] = None,
        use_prod: bool = False,
        with_authentication: bool = True,
    ) -> "SoluteGenerationInterface":
        return get_interface_for_upload_id_and_job_id(
            upload_id,
            cls,
            "SoluteGenFlow",
            job_id,
            use_prod=use_prod,
            with_authentication=with_authentication,
            find_first_job_id=False,
        )

    @classmethod
    def from_entry(
        cls, entry: NomadEntry, with_authentication: bool = False
    ) -> "SoluteGenerationInterface":
        raise get_interface_for_entry(
            entry,
            cls,
            "SoluteGenFlow",
            with_authentication=with_authentication,
            find_first_job_id=False,
        )

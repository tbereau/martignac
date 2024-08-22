from dataclasses import dataclass
from typing import Optional

from martignac.nomad.entries import NomadEntry
from martignac.workflow_interfaces.generic import (
    GenericInterface,
    get_interface_for_entry,
    get_interface_for_upload_id_and_job_id,
)


@dataclass(frozen=True)
class SoluteInSolventGenerationInterface(GenericInterface):
    solute_solvent_gro: str
    solute_solvent_top: str
    workflows: dict

    @property
    def solvent_name(self) -> str:
        return self.state_point["solvent_name"]

    @property
    def solute_name(self) -> str:
        return self.state_point["solute_name"]

    @classmethod
    def from_upload(
        cls,
        upload_id: str,
        job_id: Optional[str] = None,
        use_prod: bool = False,
        with_authentication: bool = True,
    ) -> "SoluteInSolventGenerationInterface":
        return get_interface_for_upload_id_and_job_id(
            upload_id,
            cls,
            "SoluteInSolventGenFlow",
            job_id,
            use_prod=use_prod,
            with_authentication=with_authentication,
            find_first_job_id=False,
        )

    @classmethod
    def from_entry(
        cls, entry: NomadEntry, with_authentication: bool = False
    ) -> "SoluteInSolventGenerationInterface":
        return get_interface_for_entry(
            entry,
            cls,
            "SoluteInSolventGenFlow",
            with_authentication=with_authentication,
            find_first_job_id=False,
        )

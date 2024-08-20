from dataclasses import dataclass
from typing import Optional

from martignac.workflow_interfaces.generic import (
    GenericInterface,
    get_interface_for_upload_id_and_job_id,
)


@dataclass(frozen=True)
class BilayerGenerationInterface(GenericInterface):
    lipid_names: list[str]
    bilayer_gro: str
    bilayer_top: str

    @property
    def solvent(self) -> str:
        return self.state_point["solvent"]

    @property
    def lipids(self) -> list[dict]:
        return self.state_point["lipids"]

    @classmethod
    def from_upload(
        cls,
        upload_id: str,
        job_id: Optional[str] = None,
        use_prod: bool = False,
        with_authentication: bool = True,
    ) -> "BilayerGenerationInterface":
        return get_interface_for_upload_id_and_job_id(
            upload_id,
            cls,
            "BilayerGenFlow",
            job_id,
            use_prod=use_prod,
            with_authentication=with_authentication,
            find_first_job_id=False,
        )

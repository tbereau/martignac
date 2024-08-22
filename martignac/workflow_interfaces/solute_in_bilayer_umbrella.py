import io
from dataclasses import dataclass
from typing import Optional

import numpy as np

from martignac.nomad.entries import NomadEntry
from martignac.nomad.uploads import get_specific_file_from_upload
from martignac.workflow_interfaces.generic import (
    GenericInterface,
    get_interface_for_entry,
    get_interface_for_upload_id_and_job_id,
)


@dataclass(frozen=True)
class SoluteInBilayerInterface(GenericInterface):
    """
    Interface for handling solute in bilayer umbrella sampling workflows.

    Attributes:
    -----------
    umbrella_pullx_xvg : str
        Path to the umbrella pullx xvg file.
    umbrella_pullf_xvg : str
        Path to the umbrella pullf xvg file.
    umbrella_log : str
        Path to the umbrella log file.
    free_energy : dict
        Dictionary containing free energy data.
    bilayer_z_mean : float
        Mean position of the bilayer along the z-axis.
    workflows : dict
        Dictionary containing workflow information.
    solute_translated_gro : str
        Path to the solute translated gro file.
    solute_bilayer_top : str
        Path to the solute bilayer topology file.
    tpr_file : str
        Path to the tpr file.
    """

    umbrella_pullx_xvg: str
    umbrella_pullf_xvg: str
    umbrella_log: str
    free_energy: dict
    bilayer_z_mean: float
    workflows: dict
    solute_translated_gro: str
    solute_bilayer_top: str
    tpr_file: str

    @property
    def lipids(self) -> list[dict]:
        return self.state_point["lipids"]

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
    ) -> "SoluteInBilayerInterface":
        return get_interface_for_upload_id_and_job_id(
            upload_id,
            cls,
            "SoluteInBilayerUmbrellaFlow",
            job_id,
            use_prod=use_prod,
            with_authentication=with_authentication,
            find_first_job_id=True,
        )

    @classmethod
    def from_entry(
        cls, entry: NomadEntry, with_authentication: bool = False
    ) -> "SoluteInBilayerInterface":
        return get_interface_for_entry(
            entry,
            cls,
            "SoluteInBilayerUmbrellaFlow",
            with_authentication=with_authentication,
            find_first_job_id=True,
        )

    def get_wham_npy(self, use_bootstrap: bool = True) -> np.ndarray:
        """
        Retrieve and load the WHAM (Weighted Histogram Analysis Method) numpy array.

        This method constructs the path to the WHAM file based on the `free_energy` attribute
        and retrieves the file from the upload. The file is then loaded into a numpy array.

        Parameters:
        -----------
        use_bootstrap : bool, optional
            If True, use the bootstrap key to construct the file path. If False, use the profile key.
            Default is True.

        Returns:
        --------
        np.ndarray
            The WHAM data loaded into a numpy array.
        """
        bootstrap_key = "bootstrap" if use_bootstrap else "profile"
        path_to_file = f"{self.free_energy['job_id']}/{self.free_energy[bootstrap_key]}"
        response = get_specific_file_from_upload(
            self.upload_id,
            path_to_file=path_to_file,
            use_prod=self.use_prod,
            with_authentication=self.with_authentication,
            return_json=False,
        )
        return np.load(io.BytesIO(response))

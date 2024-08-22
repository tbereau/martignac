from dataclasses import dataclass
from typing import Optional

from marshmallow import Schema, fields, pre_load

from martignac.nomad.entries import NomadEntry
from martignac.workflow_interfaces.generic import (
    GenericInterface,
    get_interface_for_entry,
    get_interface_for_upload_id_and_job_id,
)


@dataclass
class FreeEnergy:
    mean: float
    std: float


class FreeEnergySchema(Schema):
    mean = fields.Float(required=True)
    std = fields.Float(required=True)


class AlchemicalSchema(Schema):
    @pre_load
    def convert_free_energy(self, data, **kwargs):
        data["_free_energy"] = data["free_energy"]
        del data["free_energy"]
        return data


@dataclass(frozen=True)
class SoluteInSolventAlchemicalInterface(GenericInterface):
    """
    Interface for handling solute in solvent alchemical workflows.

    Attributes:
    -----------
    alchemical_log : str
        Path to the alchemical log file.
    alchemical_xvg : str
        Path to the alchemical xvg file.
    _free_energy : dict
        Dictionary containing free energy data.
    num_lambda_points : int
        Number of lambda points used in the alchemical calculation.
    system_prepared : bool
        Indicates whether the system has been prepared.
    workflows : dict
        Dictionary containing workflow information.
    """

    alchemical_log: str
    alchemical_xvg: str
    _free_energy: dict
    num_lambda_points: int
    system_prepared: bool
    workflows: dict

    @property
    def free_energy(self) -> FreeEnergy:
        """
        Retrieve the free energy as a FreeEnergy object.

        Returns:
        --------
        FreeEnergy
            The free energy data.
        """
        return FreeEnergy(**self._free_energy)

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
    ) -> "SoluteInSolventAlchemicalInterface":
        """
        Create an instance of SoluteInSolventAlchemicalInterface from an upload.

        Parameters:
        -----------
        upload_id : str
            The ID of the upload.
        job_id : Optional[str], optional
            The ID of the job, by default None.
        use_prod : bool, optional
            Whether to use production settings, by default False.
        with_authentication : bool, optional
            Whether to use authentication, by default True.
        find_first_job_id : bool, optional
            Whether to find the first job ID, by default False.

        Returns:
        --------
        SoluteInSolventAlchemicalInterface
            An instance of the interface.
        """
        return get_interface_for_upload_id_and_job_id(
            upload_id,
            cls,
            "SoluteInSolventAlchemicalFlow",
            job_id,
            use_prod=use_prod,
            with_authentication=with_authentication,
            base_schema=AlchemicalSchema,
            find_first_job_id=True,
        )

    @classmethod
    def from_entry(
        cls, entry: NomadEntry, with_authentication: bool = False
    ) -> "SoluteInSolventAlchemicalInterface":
        return get_interface_for_entry(
            entry,
            cls,
            "SoluteInSolventAlchemicalFlow",
            with_authentication=with_authentication,
            find_first_job_id=True,
        )

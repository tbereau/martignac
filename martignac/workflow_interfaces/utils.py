import json

from martignac.nomad.entries import NomadEntry
from martignac.workflow_interfaces.bilayer_generation import BilayerGenerationInterface
from martignac.workflow_interfaces.generic import Interface
from martignac.workflow_interfaces.solute_generation import SoluteGenerationInterface
from martignac.workflow_interfaces.solute_in_bilayer_umbrella import (
    SoluteInBilayerInterface,
)
from martignac.workflow_interfaces.solute_in_solvent_alchemical import (
    SoluteInSolventAlchemicalInterface,
)
from martignac.workflow_interfaces.solute_in_solvent_generation import (
    SoluteInSolventGenerationInterface,
)
from martignac.workflow_interfaces.solvent_generation import SolventGenerationInterface


def convert_entry_to_specific_interface(
    entry: NomadEntry, use_prod: bool = False, with_authentication: bool = True
) -> Interface:
    if not entry.comment:
        raise ValueError("missing entry comment")
    workflow_name = json.loads(entry.comment).get("workflow_name")
    workflows = [
        "SoluteGenFlow",
        "SolventGenFlow",
        "SoluteInSolventGenFlow",
        "SoluteInSolventAlchemicalFlow",
        "BilayerGenFlow",
        "SoluteInBilayerUmbrellaFlow",
    ]
    interfaces = [
        SoluteGenerationInterface,
        SolventGenerationInterface,
        SoluteInSolventGenerationInterface,
        SoluteInSolventAlchemicalInterface,
        BilayerGenerationInterface,
        SoluteInBilayerInterface,
    ]
    for workflow, interface in zip(workflows, interfaces):
        if workflow_name == workflow:
            return interface.from_upload(
                entry.upload_id,
                use_prod=use_prod,
                with_authentication=with_authentication,
            )
    raise ValueError(f"could not find specific interface for entry {entry.entry_id}")

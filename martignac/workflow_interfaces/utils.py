import json

from martignac.nomad.entries import (
    NomadEntry,
    get_entry_by_id,
    get_multiple_entries_by_id,
)
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
                upload_id=entry.upload_id,
                with_authentication=with_authentication,
                use_prod=use_prod,
            )
    raise ValueError(f"could not find specific interface for entry {entry.entry_id}")


def convert_entry_id_to_specific_interface(
    entry_id: str, use_prod: bool = False, with_authentication: bool = True
) -> Interface:
    entry = get_entry_by_id(
        entry_id, use_prod=use_prod, with_authentication=with_authentication
    )
    return convert_entry_to_specific_interface(
        entry, use_prod=use_prod, with_authentication=with_authentication
    )


def convert_multiple_entry_ids_to_specific_interfaces(
    entry_ids: tuple[str, ...], use_prod: bool = False, with_authentication: bool = True
) -> list[Interface]:
    entries = get_multiple_entries_by_id(
        tuple(entry_ids),
        use_prod=use_prod,
        with_authentication=with_authentication,
    )
    return [
        convert_entry_to_specific_interface(
            entry, use_prod=use_prod, with_authentication=with_authentication
        )
        for entry in entries
    ]

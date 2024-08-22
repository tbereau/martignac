import secrets

import pytest

from martignac.nomad.datasets import (
    create_dataset,
    delete_dataset,
    get_dataset_by_id,
    retrieve_datasets,
)


@pytest.mark.order(1)
def test_create_nomad_dataset(global_state: dict):
    # Generate a random string of hexadecimal characters
    hex_string = secrets.token_hex(16)

    print("Random Hexadecimal String:", hex_string)

    # Create the dataset
    dataset_name = f"Test dataset {hex_string}"
    global_state["dataset_name"] = dataset_name
    dataset_id = create_dataset(dataset_name, use_prod=False)

    # Verify the dataset was created
    assert dataset_id is not None

    global_state["dataset_id"] = dataset_id


def test_get_nomad_dataset_by_id(global_state: dict):
    assert global_state["dataset_id"] is not None
    dataset_id = global_state["dataset_id"]

    # Retrieve the dataset
    dataset = get_dataset_by_id(dataset_id, use_prod=False)

    # Verify the dataset was retrieved
    assert dataset.dataset_id == dataset_id
    assert dataset.dataset_name == global_state["dataset_name"]


def test_get_nomad_dataset_by_name(global_state: dict):
    assert global_state["dataset_name"] is not None
    dataset_name = global_state["dataset_name"]

    # Retrieve the dataset
    datasets = retrieve_datasets(dataset_name=dataset_name, use_prod=False)
    assert len(datasets) > 0


@pytest.mark.order("last")
def test_delete_nomad_dataset(global_state: dict):
    assert global_state["dataset_id"] is not None
    dataset_id = global_state["dataset_id"]

    # Delete the dataset
    delete_dataset(dataset_id, use_prod=False)

    # Verify the dataset was deleted
    # Attempt to retrieve the dataset to ensure it no longer exists
    with pytest.raises(ValueError, match=f"Problem retrieving dataset {dataset_id}"):
        get_dataset_by_id.cache.clear()
        get_dataset_by_id(dataset_id, use_prod=False)

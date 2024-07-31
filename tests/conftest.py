from typing import Optional

import pytest
from confuse.core import Configuration

from martignac import config


@pytest.fixture(scope="session")
def global_state() -> dict[str, Optional[str]]:
    return {"dataset_id": None, "dataset_name": None}


@pytest.fixture(scope="session")
def conf() -> Configuration:
    conf = config()
    conf.set_file("tests/config_tests.yaml")
    return conf


@pytest.fixture(scope="session")
def solute_gen_workspace_path(conf) -> str:
    return f"{conf['local']['workspaces']['absolute_path']}/{conf['solute_generation']['relative_paths']['workspaces']}"

from os import makedirs
from os.path import isdir

from martignac import config
from martignac.workflows.bilayer_generation import BilayerGenFlow
from martignac.workflows.solute_generation import SoluteGenFlow
from martignac.workflows.solute_in_bilayer_umbrella import SoluteInBilayerUmbrellaFlow
from martignac.workflows.solute_in_solvent_alchemical import (
    SoluteInSolventAlchemicalFlow,
)
from martignac.workflows.solute_in_solvent_generation import SoluteInSolventGenFlow
from martignac.workflows.solvent_generation import SolventGenFlow


def init_for_streamlit():
    conf = config()
    conf.set_file("/mount/src/martignac/dashboard/config_st.yaml")

    workflows = [
        SoluteGenFlow,
        SolventGenFlow,
        SoluteInSolventGenFlow,
        BilayerGenFlow,
        SoluteInSolventAlchemicalFlow,
        SoluteInBilayerUmbrellaFlow,
    ]

    for workflow in workflows:
        if not isdir(workflow.workspace_path):
            makedirs(workflow.workspace_path)
        workflow.get_project(workflow.workspace_path)

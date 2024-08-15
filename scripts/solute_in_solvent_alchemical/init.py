import itertools
from os import makedirs
from os.path import isdir

import signac

from martignac.workflows.solute_in_solvent_alchemical import (
    SoluteInSolventAlchemicalFlow,
)

if not isdir(SoluteInSolventAlchemicalFlow.workspace_path):
    makedirs(SoluteInSolventAlchemicalFlow.workspace_path)

project = signac.init_project(path=SoluteInSolventAlchemicalFlow.workspace_path)

solvent_names = ["HD"]  # ["HD", "OCO", "CLF", "ETH", "BENZ", "CHEX", "W"]
solute_names = ["P6"]  # ["P6", "Q1", "D", "N3a"]
lambda_states = range(3)  # 11

triplets = list(itertools.product(solvent_names, solute_names, lambda_states))

for solvent_name, solute_name, lambda_state in triplets:
    sp = {
        "type": "solute_in_solvent_alchemical",
        "solvent_name": solvent_name,
        "solute_name": solute_name,
        "lambda_state": lambda_state,
    }
    job = project.open_job(sp).init()

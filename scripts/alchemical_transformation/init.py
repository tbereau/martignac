import itertools

import signac

from martignac.workflows.free_energy_calculations.alchemical_transformation import AlchemicalTransformationFlow

project = signac.init_project(path=AlchemicalTransformationFlow.workspace_path)

solvent_names = [
    # "W",
    "OCO"
]
solute_names = ["P6"]
lambda_states = range(4)  # 11

triplets = list(itertools.product(solvent_names, solute_names, lambda_states))

for solvent_name, solute_name, lambda_state in triplets:
    sp = {
        "type": "alchemical_transformation",
        "solvent_name": solvent_name,
        "solute_name": solute_name,
        "lambda_state": lambda_state,
    }
    job = project.open_job(sp).init()

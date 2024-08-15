import itertools

import signac

from martignac.workflows.solute_in_solvent_generation import SoluteInSolventGenFlow

project = signac.init_project(path=SoluteInSolventGenFlow.workspace_path)

solvent_names = ["HD"]  # "ETH", "BENZ"]
solute_names = ["P6"]

pairs = list(itertools.product(solute_names, solvent_names))

for solute_name, solvent_name in pairs:
    sp = {
        "type": "solute_in_solvent_generation",
        "solvent_name": solvent_name,
        "solute_name": solute_name,
    }
    job = project.open_job(sp).init()

import signac

from martignac.workflows.system_generation.solute_in_solvent import SoluteInSolventFlow

project = signac.init_project(path=SoluteInSolventFlow.workspace_path)

solvent_names = ["HD"]
solute_names = ["P6"]

for solute_name, solvent_name in zip(solute_names, solvent_names):
    sp = {"type": "solute_solvation", "solvent_name": solvent_name, "solute_name": solute_name}
    job = project.open_job(sp).init()

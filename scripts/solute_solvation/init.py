import signac

from martignac.workflows.solute_solvation import SoluteSolvationFlow

project = signac.init_project(path=SoluteSolvationFlow.workspace_path)

solvent_names = ["HD"]
solute_names = ["P6"]

for solute_name, solvent_name in zip(solute_names, solvent_names):
    sp = {"type": "solute_solvation", "solvent_name": solvent_name, "solute_name": solute_name}
    job = project.open_job(sp).init()

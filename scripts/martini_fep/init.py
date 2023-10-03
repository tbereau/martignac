import signac
import itertools

project = signac.init_project()

solvent_names = ["HD"]
solute_names = ["P6"]
lambda_states = range(2)

triplets = list(itertools.product(solvent_names, solute_names, lambda_states))

for solvent_name, solute_name, lambda_state in triplets:
    sp = {"solvent_name": solvent_name, "solute_name": solute_name, "lambda_state": lambda_state}
    job = project.open_job(sp).init()

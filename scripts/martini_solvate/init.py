import signac

project = signac.init_project()

solvent_names = ["HD"]
solute_names = ["P6"]

for solute_name, solvent_name in zip(solute_names, solvent_names):
    sp = {"solvent_name": solvent_name, "solute_name": solute_name}
    job = project.open_job(sp).init()

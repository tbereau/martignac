import signac

project = signac.init_project()

for solvent_name in ["HD"]:
    sp = {"type": "solvent", "solvent_name": solvent_name}
    job = project.open_job(sp).init()

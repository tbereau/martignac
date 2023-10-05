import signac

project = signac.init_project()

for solvent_name in ["OCO", "W"]:
    sp = {"type": "solvent", "solvent_name": solvent_name}
    job = project.open_job(sp).init()

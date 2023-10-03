import signac

project = signac.init_project()

for solute_name in ["P6"]:
    sp = {"type": "solute", "solute_name": solute_name}
    job = project.open_job(sp).init()

import signac

project = signac.init_project()

for lambda_state in range(11):
    sp = {"lambda_state": lambda_state}
    job = project.open_job(sp).init()

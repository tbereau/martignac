import signac

from martignac.workflows.solute_generation import SoluteGenFlow

project = signac.init_project(path=SoluteGenFlow.workspace_path)

for solute_name in ["P6"]:
    sp = {"type": "solute", "solute_name": solute_name}
    job = project.open_job(sp).init()

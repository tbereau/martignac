import signac

from martignac.workflows.system_generation.solute import SoluteGenFlow

project = signac.init_project(path=SoluteGenFlow.workspace_path)

for solute_name in ["C4"]:
    sp = {"type": "solute", "solute_name": solute_name}
    job = project.open_job(sp).init()

import signac

from martignac.workflows.system_generation.solvent import SolventGenFlow

project = signac.init_project(path=SolventGenFlow.workspace_path)

for solvent_name in ["HD"]:
    sp = {"type": "solvent", "solvent_name": solvent_name}
    job = project.open_job(sp).init()

import signac

from martignac.workflows.solvent_generation import SolventGenFlow

project = signac.init_project(path=SolventGenFlow.workspace_path)

for solvent_name in ["HD", "BENZ"]:  # ["HD", "OCO", "CLF", "ETH", "BENZ", "CHEX", "W"]:
    sp = {"type": "solvent", "solvent_name": solvent_name}
    job = project.open_job(sp).init()

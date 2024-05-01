import signac

from martignac.liquid_models.mixtures import LiquidComponent, LiquidMixture
from martignac.workflows.system_generation.solute_in_bilayer import SoluteInBilayerFlow

project = signac.init_project(path=SoluteInBilayerFlow.workspace_path)

lipids = [LiquidMixture([LiquidComponent("POPC", 1.0)])]
lipid_names = []
for lipid in lipids:
    lipid_names.append([{"name": c.name, "fraction": c.fraction} for c in lipid.components])
solute_names = ["P6"]

for solute, lipid in zip(solute_names, lipid_names):
    sp = {"type": "solute_in_bilayer", "lipids": lipid, "solute_name": solute}
    job = project.open_job(sp).init()

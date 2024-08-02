import signac

from martignac.liquid_models.mixtures import LiquidComponent, LiquidMixture
from martignac.workflows.bilayer_generation import BilayerGenFlow

project = signac.init_project(path=BilayerGenFlow.workspace_path)

for lipid_name in [LiquidMixture([LiquidComponent("M3.POPC", 1.0)])]:
    sp = {
        "type": "bilayer",
        "lipids": [{"name": c.name, "fraction": c.fraction} for c in lipid_name.components],
        "solvent": BilayerGenFlow.solvent.to_insane_format(),
    }
    job = project.open_job(sp).init()

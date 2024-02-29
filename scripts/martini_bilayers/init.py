import signac

from martignac.liquid_models.mixtures import LiquidComponent, LiquidMixture

project = signac.init_project()

for lipid_name in [LiquidMixture([LiquidComponent("POPC", 1.0)])]:
    sp = {
        "type": "bilayer",
        "lipids": [{"name": c.name, "fraction": c.fraction} for c in lipid_name.components],
    }
    job = project.open_job(sp).init()

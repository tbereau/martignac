import itertools

import signac

from martignac.liquid_models.mixtures import LiquidComponent, LiquidMixture
from martignac.workflows.system_generation.solute_in_bilayer import SoluteInBilayerFlow

project = signac.init_project(path=SoluteInBilayerFlow.workspace_path)

lipids = [LiquidMixture([LiquidComponent("POPC", 1.0)])]
lipid_names = []
for lipid in lipids:
    lipid_names.append([{"name": c.name, "fraction": c.fraction} for c in lipid.components])
solute_names = ["P6"]
depths_from_bilayer_core = [0.0]  # in nm

triplets = list(itertools.product(solute_names, lipid_names, depths_from_bilayer_core))

for solute, lipid, depth in triplets:
    sp = {
        "type": "solute_in_bilayer",
        "lipids": lipid,
        "solute_name": solute,
        "depth_from_bilayer_core": depth,
    }
    job = project.open_job(sp).init()
import itertools
from os import makedirs
from os.path import isdir

import signac

from martignac.liquid_models.mixtures import LiquidComponent, LiquidMixture
from martignac.workflows.solute_in_bilayer_umbrella import SoluteInBilayerUmbrellaFlow

if not isdir(SoluteInBilayerUmbrellaFlow.workspace_path):
    makedirs(SoluteInBilayerUmbrellaFlow.workspace_path)

project = signac.init_project(path=SoluteInBilayerUmbrellaFlow.workspace_path)

lipids = [LiquidMixture([LiquidComponent("M3.POPC", 0.8), LiquidComponent("M3.DPPC", 0.2)])]
lipid_names = []
for lipid in lipids:
    lipid_names.append([{"name": c.name, "fraction": c.fraction} for c in lipid.components])
solute_names = [
    "P6",
    # "C1 P3, 0-1",
]
depths_from_bilayer_core = [1.0, 1.05, 1.10]  # np.arange(0.0, 4.05, 0.05)  # in nm

triplets = list(itertools.product(solute_names, lipid_names, depths_from_bilayer_core))

for solute, lipid, depth in triplets:
    sp = {
        "type": "solute_in_bilayer",
        "lipids": lipid,
        "solute_name": solute,
        "depth_from_bilayer_core": float(depth),
    }
    job = project.open_job(sp).init()

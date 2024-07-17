# Running workflows

!!! warning "Safeguards for first-time users"

    The [installation guide](install.md) provides a detailed explanation of the safeguards that are in place to prevent
    accidental data upload to NOMAD. Please make sure to read it before running Martignac for the first time.
    This assumes you've also made it through the [Quick start](quickstart.md).

Here's a list of workflows currently implemented in Martignac:

| Workflow                     | Description                                        | Class name                                           |
|------------------------------|----------------------------------------------------|------------------------------------------------------|
| Solute generation            | generate a solute in the gas phase                 | [martignac.workflows.solute_generation][]            |
| Solvent generation           | generate a solvent: homogeneous liquid in a box    | [martignac.workflows.solvent_generation][]           |
| Solute-in-solvent generation | generate a solute in a box of homogeneous liquid   | [martignac.workflows.solute_in_solvent_generation][] |
| Solute-in-solvent alchemical | alchemical transformation of a solute in a solvent | [martignac.workflows.solute_in_solvent_alchemical][] |
| Bilayer generation           | generate a bilayer                                 | [martignac.workflows.bilayer_generation][]           |
| Solute-in-bilayer umbrella   | umbrella sampling of a solute in a bilayer         | [martignac.workflows.solute_in_bilayer_umbrella][]   |

and a corresponding set of example `init.py` files to get started with screening each workflow:

### Solute generation

specify the solute chemistry

``` python
import signac

from martignac.workflows.solute_generation import SoluteGenFlow

project = signac.init_project(path=SoluteGenFlow.workspace_path)

for solute_name in ["C4", "P6 C3, 0-1", "P6 C3 N1, 0_1 1_2 0-2"]:
    sp = {"type": "solute", "solute_name": solute_name}
    job = project.open_job(sp).init()
```

### Solvent generation

specify the solvent chemistry

``` python
import signac

from martignac.workflows.solvent_generation import SolventGenFlow

project = signac.init_project(path=SolventGenFlow.workspace_path)

for solvent_name in ["HD"]:
    sp = {"type": "solvent", "solvent_name": solvent_name}
    job = project.open_job(sp).init()
```

### Solute-in-solvent generation

specify the solute and solvent chemistry

``` python
import signac

from martignac.workflows.solute_in_solvent_generation import SoluteInSolventGenFlow

project = signac.init_project(path=SoluteInSolventGenFlow.workspace_path)

solvent_names = ["HD"]
solute_names = ["P6"]

for solute_name, solvent_name in zip(solute_names, solvent_names):
    sp = {"type": "solute_solvation", "solvent_name": solvent_name, "solute_name": solute_name}
    job = project.open_job(sp).init()
```

### Solute-in-solvent alchemical

specify the solute and solvent chemistry, and the lambda states for the alchemical transformation

``` python
import itertools

import signac

from martignac.workflows.solute_in_solvent_alchemical import SoluteInSolventAlchemicalFlow

project = signac.init_project(path=SoluteInSolventAlchemicalFlow.workspace_path)

solvent_names = ["W", "OCO"]
solute_names = ["P6"]
lambda_states = range(11)

triplets = list(itertools.product(solvent_names, solute_names, lambda_states))

for solvent_name, solute_name, lambda_state in triplets:
    sp = {
        "type": "alchemical_transformation",
        "solvent_name": solvent_name,
        "solute_name": solute_name,
        "lambda_state": lambda_state,
    }
    job = project.open_job(sp).init()
```

### Bilayer generation

specify the phospholipid chemistry (solvent, e.g., water, is specific in the `config.yaml` file)

``` python
import signac

from martignac.liquid_models.mixtures import LiquidComponent, LiquidMixture
from martignac.workflows.bilayer_generation import BilayerGenFlow

project = signac.init_project(path=BilayerGenFlow.workspace_path)

for lipid_name in [LiquidMixture([LiquidComponent("POPC", 1.0)])]:
    sp = {
        "type": "bilayer",
        "lipids": [{"name": c.name, "fraction": c.fraction} for c in lipid_name.components],
    }
    job = project.open_job(sp).init()
```

### Solute-in-bilayer umbrella

specify the solute and bilayer chemistry, as well as the set of umbrella restraints

``` python
import itertools

import signac
import numpy as np

from martignac.liquid_models.mixtures import LiquidComponent, LiquidMixture
from martignac.workflows.solute_in_bilayer_umbrella import SoluteInBilayerUmbrellaFlow

project = signac.init_project(path=SoluteInBilayerUmbrellaFlow.workspace_path)

lipids = [LiquidMixture([LiquidComponent("POPC", 1.0)])]
lipid_names = []
for lipid in lipids:
    lipid_names.append([{"name": c.name, "fraction": c.fraction} for c in lipid.components])
solute_names = ["C4"]
depths_from_bilayer_core = np.linspace(0.0, 4.0, 41)  # in nm

triplets = list(itertools.product(solute_names, lipid_names, depths_from_bilayer_core))

for solute, lipid, depth in triplets:
    sp = {
        "type": "solute_in_bilayer",
        "lipids": lipid,
        "solute_name": solute,
        "depth_from_bilayer_core": depth,
    }
    job = project.open_job(sp).init()
```

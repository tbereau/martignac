# Quick start

Make sure you first install and set up Martignac as described in the [installation guide](install.md).

## Hello, Martignac

!!! warning "Safeguards for first-time users"

    The [installation guide](install.md) provides a detailed explanation of the safeguards that are in place to prevent
    accidental data upload to NOMAD. Please make sure to read it before running Martignac for the first time.

You are now ready to run the workflows. 
For instance, to generate solutes in bilayers, go to the input absolute path (e.g., `scripts/solute_in_bilayer`).
There should be two files 

`init.py` 

:   initializes the computational `signac-workflow`; defines the solutes chemistries to screen (here one single
    C4 bead, one P6-C3 molecule, and one P6-C3-N1 molecule). This solely tells `signac` which computational jobs are 
    to be _defined_. The actual computation is defined in the `project.py` script.

``` python
import signac

from martignac.workflows.solute_generation import SoluteGenFlow

project = signac.init_project(path=SoluteGenFlow.workspace_path)

for solute_name in ["C4", "P6 C3, 0-1", "P6 C3 N1, 0_1 1_2 0-2"]:
    sp = {"type": "solute", "solute_name": solute_name}
    job = project.open_job(sp).init()
```

`project.py`

:   defines the workflow to run by importing the class and setting the desired workspace path. All computations 
    will occur from this script, and the results will be stored in the `workspace/` directory.

``` python
from martignac.workflows.solute_generation import SoluteGenFlow

if __name__ == "__main__":
    SoluteGenFlow(path=SoluteGenFlow.workspace_path).main()
```

You can simply run the two scripts by the following commands:
```bash
python init.py
python project.py run
```

The results can be found in your output `workspace/` directory. 
If you've decided to push the data to NOMAD, you can check the content of your uploads on the NOMAD webserver,
by logging in and navigating to your Datasets. For instance, the url for the test database would be 
<https://nomad-lab.eu/prod/v1/test/gui/user/datasets>.

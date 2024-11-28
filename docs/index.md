# Martignac: Coarse-grained Martini simulation workflows

Martignac is a Python-based toolkit designed to streamline the preparation and analysis of coarse-grained Martini 
molecular dynamics simulations. 
It provides a suite of utilities for generating initial configurations---including solutes, solvent boxes, and 
phospholipid bilayers---runs the simulations, and performs analysis, 
such as alchemical transformations and umbrella sampling . 
Martignac directly connects to the NOMAD database to pull existing simulations and push missing ones to the database,
thereby constantly enriching the database with new data.

The toolkit efficiently sets up and runs simulations, and perform analysis. 
By automating routine tasks and offering a structured approach to simulation setup, Martignac aims to enhance 
productivity and reproducibility in molecular dynamics studies.

<div class="grid cards" markdown>

- [:material-arrow-down-box: __Install__ `martignac` to get started](install.md)
- [:material-clock-fast: __Quick start__ to run your first workflow](quickstart.md)
- [:material-run-fast: __Workflows__ explore all workflows](running_workflows.md)
- [:material-file-document: __Reference__ check out the docs](nomad/datasets.md)

</div>

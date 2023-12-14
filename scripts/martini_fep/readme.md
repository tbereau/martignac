# Martignac alchemical free-energy calculations 

Tristan Bereau. 
Dated: Thu Dec. 14, 2023

## Project structure

Computation of solvation free energies for a solute inside a box of a solvent. Everything computed at the CG Martini 3 resolution. `signac-flow` is used to handle the workflow.

The present example contains the following systems:

- P6 in water
  - Solute: a single P6 bead
  - Solvent: Water
- P6 in octanol
  - Solute: a single P6 bead
  - Solvent: Octanol

For each system, alchemical free energies with 11 values of the coupling parameter are used:
```
vdw-lambdas              = 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 
```

## File structure

Each coupling parameter of each system (i.e., **no hierarchy**) will be its own `signac` statepoint. The relevant `signac-flow` files are:
```
❯ ll *.py
-rw-r--r--@ 1 bereau  staff   259B Oct  3 16:24 dashboard.py
-rw-r--r--@ 1 bereau  staff   416B Dec 14 09:27 init.py
-rw-r--r--@ 1 bereau  staff   158B Dec 14 09:29 project.py
```
where

- `init.py` initializes the statepoints
- `project.py` is used to run and query the status of the workflow calculations
- `dashboard.py` is used to visualize in the browser the status of the workflow.

Statepoint calculations are stored in the subolder `/workspace`. A typical statepoint will consist of the following files:
```
-rw-r--r--@ 1 bereau  staff    39K Dec 14 09:29 fep_run-6.cpt
-rw-r--r--@ 1 bereau  staff   4.0K Dec 14 09:29 fep_run-6.edr
-rw-r--r--@ 1 bereau  staff   105K Dec 14 09:29 fep_run-6.gro
-rw-r--r--@ 1 bereau  staff    24K Dec 14 09:29 fep_run-6.log
-rw-r--r--@ 1 bereau  staff    42K Dec 14 09:29 fep_run-6.tpr
-rw-r--r--@ 1 bereau  staff    42K Dec 14 09:29 fep_run-6.xtc
-rw-r--r--@ 1 bereau  staff   758K Dec 14 09:29 fep_run-6.xvg
-rw-r--r--@ 1 bereau  staff    11K Dec 14 09:29 mdout.mdp
-rw-r--r--@ 1 bereau  staff   2.8K Dec 14 09:29 run_lambda.mdp
-rw-r--r--@ 1 bereau  staff   220B Dec 14 09:29 signac_job_document.json
-rw-r--r--@ 1 bereau  staff    63B Dec 14 09:28 signac_statepoint.json
-rw-r--r--@ 1 bereau  staff   144B Dec 14 09:28 solute.itp
-rw-r--r--@ 1 bereau  staff   230B Dec 14 09:28 solute_solvent.top
-rw-r--r--@ 1 bereau  staff   105K Dec 14 09:28 solute_solvent_equ.gro
```
where 

- this is the 6th value of the coupling parameter; `fep_run-*` files consist of the free-energy perturbation simulations
- `run_lambda.mdp` is the gromacs simulation input file
- `signac_job_document.json` points variables used in the workflow setup to the correct files (for internal use of the workflow)
- `signac_statepoint.json` is a `signac-flow` file that specifies the statepoint (e.g., `{"solvent_name": "OCO", "solute_name": "P6", "lambda_state": 6}`)
- `solute.itp` is the gromacs itp file for the solute
- `solute_solvent.top` is the gromacs top file for the soluve-solvent system
- `solute_solvent_equ.gro` is the equilibrated gromacs gro file for the solute-solvent system.

On top of the free-energy calculations, for each system there is a system preparation. This system preparation is **identical** to every values of the coupling parameter. The workflow uses an aggregator to ensure that only a single system preparation is performed, and blocks all free-energy simulations until then. System preparation consists of several steps:

- Solute
  - Generation
  - Minimization
- Solvent
  - Molecule generation
  - Box generation (using `packmol`)
  - Box minimization
  - Box equilibration
  - Box production
- Solute-Solvent
  - System generation (from the minimized solute and production solvent)
  - System minimization
  - System equilibration

The last step is the input to all alchemical free-energy simulations. The overall file structure looks like this:
```
-rw-r--r--@ 1 bereau  staff    39K Dec 14 09:29 fep_run-0.cpt
-rw-r--r--@ 1 bereau  staff   3.8K Dec 14 09:29 fep_run-0.edr
-rw-r--r--@ 1 bereau  staff   105K Dec 14 09:29 fep_run-0.gro
-rw-r--r--@ 1 bereau  staff    24K Dec 14 09:29 fep_run-0.log
-rw-r--r--@ 1 bereau  staff    42K Dec 14 09:29 fep_run-0.tpr
-rw-r--r--@ 1 bereau  staff    46K Dec 14 09:29 fep_run-0.xtc
-rw-r--r--@ 1 bereau  staff   733K Dec 14 09:29 fep_run-0.xvg
-rw-r--r--@ 1 bereau  staff    11K Dec 14 09:29 mdout.mdp
-rw-r--r--@ 1 bereau  staff   150B Dec 14 09:28 packmol.inp
-rw-r--r--@ 1 bereau  staff   2.8K Dec 14 09:29 run_lambda.mdp
-rw-r--r--@ 1 bereau  staff   466B Dec 14 09:29 signac_job_document.json
-rw-r--r--@ 1 bereau  staff    61B Dec 14 09:28 signac_statepoint.json
-rw-r--r--@ 1 bereau  staff   144B Dec 14 09:28 solute.itp
-rw-r--r--@ 1 bereau  staff   133B Dec 14 09:28 solute.top
-rw-r--r--@ 1 bereau  staff   104B Dec 14 09:28 solute_gen_mol.gro
-rw-r--r--@ 1 bereau  staff   728B Dec 14 09:28 solute_min.edr
-rw-r--r--@ 1 bereau  staff    92B Dec 14 09:28 solute_min.gro
-rw-r--r--@ 1 bereau  staff    15K Dec 14 09:28 solute_min.log
-rw-r--r--@ 1 bereau  staff   2.7K Dec 14 09:28 solute_min.tpr
-rw-r--r--@ 1 bereau  staff   132B Dec 14 09:28 solute_min.trr
-rw-r--r--@ 1 bereau  staff   231B Dec 14 09:29 solute_solvent.top
-rw-r--r--@ 1 bereau  staff    38K Dec 14 09:29 solute_solvent_equ.cpt
-rw-r--r--@ 1 bereau  staff   3.6K Dec 14 09:29 solute_solvent_equ.edr
-rw-r--r--@ 1 bereau  staff   105K Dec 14 09:29 solute_solvent_equ.gro
-rw-r--r--@ 1 bereau  staff    21K Dec 14 09:29 solute_solvent_equ.log
-rw-r--r--@ 1 bereau  staff    41K Dec 14 09:29 solute_solvent_equ.tpr
-rw-r--r--@ 1 bereau  staff    46K Dec 14 09:29 solute_solvent_equ.xtc
-rw-r--r--@ 1 bereau  staff    69K Dec 14 09:29 solute_solvent_gen.gro
-rw-r--r--@ 1 bereau  staff    66K Dec 14 09:29 solute_solvent_min.edr
-rw-r--r--@ 1 bereau  staff    69K Dec 14 09:29 solute_solvent_min.gro
-rw-r--r--@ 1 bereau  staff   102K Dec 14 09:29 solute_solvent_min.log
-rw-r--r--@ 1 bereau  staff    23K Dec 14 09:29 solute_solvent_min.tpr
-rw-r--r--@ 1 bereau  staff    18K Dec 14 09:29 solute_solvent_min.trr
-rw-r--r--@ 1 bereau  staff    38K Dec 14 09:29 solvent_equ.cpt
-rw-r--r--@ 1 bereau  staff   3.5K Dec 14 09:29 solvent_equ.edr
-rw-r--r--@ 1 bereau  staff   105K Dec 14 09:29 solvent_equ.gro
-rw-r--r--@ 1 bereau  staff    21K Dec 14 09:29 solvent_equ.log
-rw-r--r--@ 1 bereau  staff    41K Dec 14 09:28 solvent_equ.tpr
-rw-r--r--@ 1 bereau  staff    46K Dec 14 09:29 solvent_equ.xtc
-rw-r--r--@ 1 bereau  staff    69K Dec 14 09:28 solvent_gen_box.gro
-rw-r--r--@ 1 bereau  staff   124K Dec 14 09:28 solvent_gen_box.pdb
-rw-r--r--@ 1 bereau  staff   174B Dec 14 09:28 solvent_gen_box.top
-rw-r--r--@ 1 bereau  staff   104B Dec 14 09:28 solvent_gen_mol.gro
-rw-r--r--@ 1 bereau  staff   207B Dec 14 09:28 solvent_gen_mol.pdb
-rw-r--r--@ 1 bereau  staff   174B Dec 14 09:28 solvent_gen_mol.top
-rw-r--r--@ 1 bereau  staff    66K Dec 14 09:28 solvent_min.edr
-rw-r--r--@ 1 bereau  staff    69K Dec 14 09:28 solvent_min.gro
-rw-r--r--@ 1 bereau  staff   102K Dec 14 09:28 solvent_min.log
-rw-r--r--@ 1 bereau  staff    22K Dec 14 09:28 solvent_min.tpr
-rw-r--r--@ 1 bereau  staff    18K Dec 14 09:28 solvent_min.trr
-rw-r--r--@ 1 bereau  staff    39K Dec 14 09:29 solvent_prod.cpt
-rw-r--r--@ 1 bereau  staff   3.8K Dec 14 09:29 solvent_prod.edr
-rw-r--r--@ 1 bereau  staff   105K Dec 14 09:29 solvent_prod.gro
-rw-r--r--@ 1 bereau  staff    21K Dec 14 09:29 solvent_prod.log
-rw-r--r--@ 1 bereau  staff    41K Dec 14 09:29 solvent_prod.tpr
-rw-r--r--@ 1 bereau  staff    46K Dec 14 09:29 solvent_prod.xtc
```

## Free-energy analysis

At the end of the workflow, an aggregator function computes the free energy by using the `MBAR` package for Python. The results are stored in the overall project's JSON document (`/signac_project_document.json`):
```
❯ cat signac_project_document.json
{"free_energies": [{"solvent_name": "W", "solute_name": "P6", "f_mean": 7.403748589036349, "f_std": 0.03610783576548312}, {"solvent_name": "OCO", "solute_name": "P6", "f_mean": 3.1498282167097624, "f_std": 0.02864011196314877}]}%
```
which provides both the mean and the standard deviation of the solvation free energy.

An alternative overview of the workflow statepoints and calculation of the free-energy calculations is provided in the jupyter notebook `/analysis.ipynb`. It makes use of `signac-flow` methods to filter workflow directories by the relevant statepoints, and recomputes the free energies of solvation, as well as plot the results.

# Martignac

!!! note ""

    Martignac: Coarse-grained Martini simulation worfklows

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

## Installation


First, clone the repository::
``` bash
git clone git@lin0.thphys.uni-heidelberg.de:bereau/martignac.git
```
Setup a virtual environment (Python 3.9 recommended) and install the package
```bash
virtualenv --python="/path/to/python/version" .venv
pip install --upgrade pip
python setup.py install
pip install -e .
```

## Martignac dependencies installation 

Update `requirements.txt` from `pyproject.toml` with piptools, and then install dependencies with pip:
```bash
pip install pip-tools
python -m piptools compile -o requirements.txt pyproject.toml
pip install -r requirements.txt
```
For development:
```bash
pip install -r test_requirements.txt
```
If you don't already have Gromacs installed, you can quickly install with conda via:
```bash
conda config --add channels bioconda
conda install conda install ocl-icd-system==1.0.0
conda install gromacs==2018.6
```

alternatively, to install `conda`:
```bash
wget -c https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
# on MacOS replace with:
# wget -c https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -O miniconda.sh
bash ./miniconda.sh -b -f -p /path/to/conda
export PATH="/path/to/conda/bin:$PATH"
source activate base
```

and gromacs
```bash
conda install -c conda-forge gromacs
```

## NOMAD

### Credentials

Martignac will need to connect to NOMAD. Create an account on <https://nomad-lab.eu/>. 
Store your credentials in a `.env` file at the root of the `martignac` directory, with the following content
```bash
NOMAD_USERNAME="MyLogin"
NOMAD_PASSWORD="MyPassword"
```
and insert your username and password. 

!!! warning "`.env` file"
    
    Never push your `.env` file to a repository. This would expose your password.

### Dataset

You will need to create your own dataset to push to NOMAD. 
See the [NOMAD datasets](nomad/datasets.md) page for more information. 
The function [martignac.nomad.datasets.create_dataset][] will help you create a dataset. 

### Production vs. test databases

There are two main databases in NOMAD: production and test. Make sure to run all tests on the **test database**.

!!! warning "Production database"
    
    Never push test data to the production database.

## Config file

Martignac is controlled through a `config.yaml`. You can find an example file in `martignac/config_default.yaml`. 
All aspects of simulation input files and parameters can be controlled from the yaml file. 

## Workflow paths 

You can set the directories for input and output of your workflows in the `config.yaml`. 

- Output directory: set the key `local > workspaces > absolute_path`
- Input directories: set the keys 
  - `local > input_files > absolute_path` (default: `scripts/`), containing python scripts and mdp files 
  - `local > input_files > itp_files`, containing force-field itp files

In your workspace directory `${WORKSPACE}`, create the necessary subdirectories and initialize `signac`:
```bash
for dir in alchemical_transformation solute_generation solute_in_solvent bilayer_generation \
  solute_in_bilayer solvent_generation; do
  mkdir -p ${dir}
  cd ${dir}
  signac init
  cd ..
done
```
you should read `Initialized project` for each directory.

## Running Martignac

You are now ready to run the workflows. 
For instance, to generate solutes in bilayers, go to the input absolute path (e.g., `scripts/solute_in_bilayer`) 
and execute the following commands:
```bash
python init.py
python project.py run
``` 

The results can be found in your output `workspace/` directory.

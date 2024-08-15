# Martignac

Martignac: Coarse-grained Martini simulation worfklows

- **GitHub repository**: <https://github.com/tbereau/martignac/>

## Getting started with your project

First, create a repository on GitHub with the same name as this project, and then run the following commands:

``` bash
git init -b main
git add .
git commit -m "init commit"
git remote add origin git@github.com:tbereau/martignac.git
git push -u origin main
```

Finally, setup a virtual environment (Python 3.9 recommended) and install the package 

```bash
virtualenv --python="/path/to/python/version" .venv
pip install --upgrade pip
python setup.py install
pip install -e .
```

## Martignac dependencies installation 

For development:

```bash
pip install -r test_requirements.txt
```

If you don't already have Gromacs installed, you can quickly install with conda via:

```bash
conda config --add channels bioconda
conda install conda install ocl-icd-system==1.0.0
conda install gromacs==2023.1
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

Martignac will need to connect to NOMAD. Create an account on https://nomad-lab.eu/. 
Store your credentials in a `.env` file at the root of the `martignac` directory, with the following content
```bash
NOMAD_USERNAME="MyLogin"
NOMAD_PASSWORD="MyPassWord"
```
and insert your username and password. 

> [!CAUTION]
> Never push your `.env` file to a repository. This would expose your password.

## Config file

Martignac is controlled through a `config.yaml`. You can find an example file in `martignac/config_default.yaml`. 
All aspects of simulation input files and parameters can be controlled from the yaml file. 

## Workflow paths 

You can set the directories for input and output of your workflows in the `config.yaml`. 

- Output directory: set the key `local > workspaces > absolute_path`
- Input directories: set the keys 
  - `local > input_files > absolute_path` (default: `scripts/`), containing python scripts and mdp files 
  - `local > input_files > itp_files`, containing force-field itp files

For each workflow, run the python scripts. 
For instance, to generate solutes in bilayers, go to the input absolute path (e.g., `scripts/solute_in_bilayer`) 
and execute the following commands:
```bash
python init.py
python project.py run
``` 

The results can be found in your output `workspace/` directory.

## Documentation

To generate the documentation, run the command:
```bash
mkdocs serve
```

for a general NOMAD API documentation, see https://nomad-lab.eu/prod/v1/api/v1/extensions/docs


[![Build and Deploy](https://github.com/tbereau/martignac/actions/workflows/main.yml/badge.svg)](https://github.com/tbereau/martignac/actions/workflows/main.yml)
[![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://martignac.streamlit.app/)
[![Docs site](https://img.shields.io/badge/docs-GitHub_Pages-blue)](https://tbereau.github.io/martignac/)


> Martignac: Coarse-grained Martini simulation worfklows

See our [preprint](https://arxiv.org/abs/2409.15478) for more details on the Martignac workflows:
```
@article{bereau2024martignac,
  title={Martignac: Computational workflows for reproducible, traceable, and composable coarse-grained Martini simulations},
  author={Bereau, Tristan and Walter, Luis J and Rudzinski, Joseph F},
  journal={arXiv preprint arXiv:2409.15478},
  year={2024}
}
```

## Getting started with your project

First, clone (or fork and then clone) the repository::
``` bash
git clone git@github.com:tbereau/martignac.git
```

Setup a virtual environment (Python 3.9 recommended) within the root directory of the local repo:

```bash
python -m venv .pyenv python==3.9
. .pyenv/bin/activate
pip install --upgrade pip
```

We recommend installing the `uv` package for efficient installations with pip:

```bash
pip install uv
```

Now, install the package dependencies (`test_requirements.txt` is only needed for development):

``` bash
uv pip install -r requirements.txt
uv pip install -r test_requirements.txt
```

Finally, run the setup and install the package:
```bash
python setup.py install
pip install -e .
```

The `-e` option installs the package in "editable" mode, which allows on the fly changes to the code without re-installing during testing or debugging phases.

## Gomacs installation

If you don't already have Gromacs installed, you can install with `apt-get` on linux (tested on Ubuntu 22.04):

```bash
sudo apt-get update
sudo apt-get install -y build-essential cmake git libfftw3-dev libgsl0-dev
sudo apt-get install -y libboost-all-dev libopenmpi-dev
sudo apt-get install -y gromacs
```

Alternatively, you can quickly install with conda via:

```bash
conda config --add channels bioconda conda-forge
conda install conda install ocl-icd-system==1.0.0
conda install gromacs==2023.1
```

In either case, you will need to make the gromacs executables available to your environment. If you use conda, make sure not to intertwine your conda and venv (i.e., only work within 1 environment at a time).

To install `conda`:
```bash
wget -c https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
# on MacOS replace with:
# wget -c https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -O miniconda.sh
bash ./miniconda.sh -b -f -p /path/to/conda
export PATH="/path/to/conda/bin:$PATH"
source activate base
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

Martignac is controlled through a `config.yaml` within the `martigac/` folder of the root directory. There you will find an example file called `config_default.yaml`. We recommend you copy this default config to `config.yaml`. You may want to initially delete some options including the dataset id and coauthor nomad ids, and add them back in later. You can add your own name under `query_authors`:

```yaml
nomad:
  upload_to_nomad: false
  publish_uploads: false
  use_prod: false
  dataset:
    id: "<dataset_id>"
  coauthors:
    - "<coauthor_nomad_id>" # coauthor #1
  query_authors:
    - "<Your_Name>"
```

You will need to add the appropriate paths to your workspace:

```yaml
local:
  workspaces:
    absolute_path: "<path_to_martignac_root>/martignac/workspaces"
  input_files:
    absolute_path: "<path_to_martignac_root>/scripts"
    itp_files: "<path_to_martignac_root>/scripts/martini_v300"
  allow_symlinks: true
```

The rest of the yaml file consists of a variety of other options which control all aspects of simulation input files and parameters.

Note: The variable `MARTIGNACDIR` ensures that your Martignac config file can be found by the code. It should be set automatically, but there have been instances of import errors associated with config variables, which can be easily fixed by explicitly setting this variable with `export MARTIGNACDIR=<path_to_martignac_root>/martignac/`.

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

## Tutorial

A tutorial in the form of a Jupyter notebook is available in the [scripts/tutorial.ipynb](https://github.com/tbereau/martignac/blob/main/scripts/tutorial.ipynb) file.


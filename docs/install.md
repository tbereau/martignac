# Installation

First, clone (or fork and then clone) the repository::
``` bash
git clone git@github.com:tbereau/martignac.git
```

To begin your Martignac project, fork the repository and then create a clone in your local workspace.

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


<!-- Update `requirements.txt` from `pyproject.toml` with piptools, and then install dependencies with pip: -->
<!-- ```bash
pip install pip-tools
python -m piptools compile -o requirements.txt pyproject.toml
pip install -r requirements.txt
``` -->

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

### Push to NOMAD

Martignac will push all data to NOMAD. For testing purposes, there are safeguards in place:

| Method            | Description                                    |
|-------------------|------------------------------------------------|
| `upload_to_nomad` | :arrow_up:     upload data to NOMAD            |
| `use_prod`        | :star: use the prod/test database              |
| `publish`         | :construction:       publish the NOMAD entries |

!!! warning "`publish` toggle"

    :construction: `publish` is not yet supported.

These safeguards can be changed in the `config.yaml` file. Note that any simulation that is published on NOMAD
**canot** be deleted.

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

!!! Note "Martignac workspace path"

    The variable `MARTIGNACDIR` ensures that your Martignac config file can be found by the code. It should be set automatically, but there have been instances of import errors associated with config variables, which can be easily fixed by explicitly setting this variable with `export MARTIGNACDIR=<path_to_martignac_root>/martignac/`.

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


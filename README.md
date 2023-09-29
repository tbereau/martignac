# martignac

This is a template repository for Python projects that use Poetry for their dependency management.

- **GitLab repository**: <https://lin0.thphys.uni-heidelberg.de:4443/bereau/martignac/>

## Getting started with your project

First, create a repository on GitHub with the same name as this project, and then run the following commands:

``` bash
git init -b main
git add .
git commit -m "init commit"
git remote add origin git@lin0.thphys.uni-heidelberg.de:bereau/martignac.git
git push -u origin main
```

Finally, install the environment and the pre-commit hooks with 

```bash
make install
```

You are now ready to start development on your project! The CI/CD
pipeline will be triggered when you open a pull request, merge to main,
or when you create a new release.

To finalize the set-up for publishing to PyPi or Artifactory, see
[here](https://fpgmaas.github.io/cookiecutter-poetry/features/publishing/#set-up-for-pypi).
For activating the automatic documentation with MkDocs, see
[here](https://fpgmaas.github.io/cookiecutter-poetry/features/mkdocs/#enabling-the-documentation-on-github).
To enable the code coverage reports, see [here](https://fpgmaas.github.io/cookiecutter-poetry/features/codecov/).

## Martignac dependencies installation 

Update `requirements.txt` from `pyproject.toml` with piptools, and then install dependencies with pip:

```bash
pip install pip-tools
python -m piptools compile -o requirements.txt pyproject.toml
pip install -r requirements.txt
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

## Examples 

Examples are found in `martini_fep/` and  `martini_solvate/`. In each case execute the following commands:
```bash
python init.py
python project.py run
``` 

The results can be found in `workspace/`

## Releasing a new version



---

Repository initiated with [fpgmaas/cookiecutter-poetry](https://github.com/fpgmaas/cookiecutter-poetry).

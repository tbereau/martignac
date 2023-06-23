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

## Martignac dependencies installation with Conda
add channels to conda config, and then install dependencies in requirements.txt

```bash
conda config --add channels conda-forge bioconda
conda install --files requirements.txt
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

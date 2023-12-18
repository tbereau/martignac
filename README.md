# martignac

## Martignac dependencies installation 

It is recommended to use [Anaconda](https://www.anaconda.com/) for installing and managing dependencies. Create a new conda environment with *Python 3.11* and install the package installer *pip*:
```
conda create -n conda-martinac python=3.11 pip
source activate conda-martinac
```
Then, install all required dependencies and the program *packmol*.
```bash
pip install -r requirements.txt
conda install -c conda-forge packmol
pip install -e .
```
We also need the Gromacs simulation package. Check the [Gromacs Website](https://manual.gromacs.org/documentation/current/install-guide/index.html) for instructions on how to install the latest version of Gromacs. Older versions of Gromacs are available via Anaconda:
```bash
conda install -c bioconda gromacs
```

## Examples 

Examples are found in `scripts/martini_fep/` and  `scripts/martini_solvate/`. In each case execute the following commands:
```bash
python init.py
python project.py run
``` 

The results can be found in `workspace/`


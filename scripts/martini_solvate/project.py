import numpy as np
from flow import FlowProject


RESOURCE_DIR = "../../resource_files/"
MDP_FILES = "../../mdp_files/"
NUM_ALCH_STATES = 11


class MartiniProject(FlowProject):
    pass


def solvate_command(
        gro_solute: str, gro_solvent: str, box_array: np.ndarray, max_sol: int
) -> str:
    box_size = " ".join([str(b) for b in box_array])
    cmd = (
        f'gmx solvate -cp {gro_solute} -cs {gro_solvent} -o solvated.gro -box {box_size} '
        f'-maxsol {max_sol}'
    )
    return cmd


def gromacs_simulation_command(
        mdp: str,
        top: str,
        gro: str,
        name: str,
        n_maxwarn: int = 10,
        n_threads: int = 1
) -> str:
    grompp_cmd = f'gmx grompp -f {mdp} -p {top} -c {gro} -o {name}.tpr -maxwarn {n_maxwarn}'
    mdrun_cmd = f'gmx mdrun -nt {n_threads} -v -deffnm {name}'
    return f"{grompp_cmd} && {mdrun_cmd}"


@MartiniProject.label
def solvated(job):
    return job.isfile('solvated.gro')


@MartiniProject.label
def minimized(job):
    return job.isfile('minimized.gro')


@MartiniProject.post(solvated)
@MartiniProject.operation(cmd=True, with_job=True)
def solvate(job):
    # TODO: Make it job dependent
    return solvate_command(
        gro_solute=RESOURCE_DIR + "Na.gro",
        gro_solvent=RESOURCE_DIR + "solvent.gro",
        box_array=np.array([5.5, 5.5, 5.5]),
        max_sol=320,
    )


@MartiniProject.pre(solvated)
@MartiniProject.post(minimized)
@MartiniProject.operation(cmd=True, with_job=True)
def minimize(job):
    # TODO: Make it job dependent
    # TODO: Generate top file
    return gromacs_simulation_command(
        mdp=MDP_FILES + "min.mdp",
        top=RESOURCE_DIR + "topol.top",
        gro="solvated.gro",
        name="min"
    )


if __name__ == "__main__":
    MartiniProject().main()

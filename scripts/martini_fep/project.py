import numpy as np
import regex
from flow import FlowProject


RESOURCE_DIR = "../../resource_files/"
MDP_FILES = "../../mdp_files/"


class MartiniProject(FlowProject):
    pass


def get_run_name(job) -> str:
    lambda_state = str(job.sp.lambda_state)
    return f"run_{lambda_state}"


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
def sampled(job):
    log_file = f"run_{job.sp.lambda_state}.log"
    if job.isfile(log_file):
        with open(job.fn(log_file)) as file:
            lines = file.read()
        return "Finished mdrun on rank" in lines
    return False


@MartiniProject.post(sampled)
@MartiniProject.operation(cmd=True, with_job=True)
def sample(job):
    lambda_state = str(job.sp.lambda_state)
    mdp_file_orig = "run2.mdp"
    run_name = get_run_name(job)
    mdp_file = f"{run_name}.mdp"
    with open(MDP_FILES + mdp_file_orig) as pipe:
        mdp_content = pipe.read()
    mdp_content = regex.sub("sedstate", lambda_state, mdp_content)
    with open(job.fn(mdp_file), "w") as pipe:
        pipe.write(mdp_content)
    return gromacs_simulation_command(
        mdp=job.fn(mdp_file),
        top=RESOURCE_DIR + "topol.top",
        gro=RESOURCE_DIR + "min.gro",
        name=run_name,
        n_threads=4,
    )


if __name__ == "__main__":
    MartiniProject().main()

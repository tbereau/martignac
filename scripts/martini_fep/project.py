import pandas as pd
import regex
import os
from os.path import join

from alchemlyb.parsing.gmx import extract_u_nk
from alchemlyb.estimators import MBAR
from flow import FlowProject, aggregator


RESOURCE_DIR = "../../resource_files/"
MDP_FILES = "../../mdp_files/"
TEMPERATURE = 298.0


class MartiniProject(FlowProject):
    pass


project = MartiniProject().get_project()


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
    grompp_cmd = (
        f'gmx grompp -f {mdp} -p {top} -c {gro} -o {name}.tpr -maxwarn {n_maxwarn}'
    )
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


@MartiniProject.label
def all_sampled(*jobs):
    result = all([
        any([".xvg" in f for f in os.listdir(job.path)])
        for job in project.find_jobs()
    ])
    return result


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


@MartiniProject.pre(all_sampled)
@MartiniProject.post(lambda *jobs: project.doc.get("free_energy", False))
@MartiniProject.operation(aggregator=aggregator(sort_by="lambda_state"))
def compute_free_energy(*jobs):
    xvg_files = [
        join(job.path, f) for job in jobs
        for f in os.listdir(job.path)
        if ".xvg" in f
    ]
    u_nk_list = [extract_u_nk(f, T=TEMPERATURE) for f in xvg_files]
    u_nk_combined = pd.concat(u_nk_list)
    mbar = MBAR().fit(u_nk_combined)
    # mbar.delta_f_.to_csv("delta_f.csv")
    # mbar.d_delta_f_.to_csv("d_delta_f.csv")
    free_energy = float(mbar.delta_f_.iloc[0, -1])
    print(f"free energy: {free_energy:+.3f} kT")
    project.document["free_energy"] = free_energy


if __name__ == "__main__":
    MartiniProject().main()

import shutil
import logging
from os.path import join, basename


__all__ = [
    "solvate_command",
    "gromacs_simulation_command",
    "copy_files_to",
]

logger = logging.getLogger(__name__)


def solvate_command(
        gro_solute: str,
        gro_solvent: str,
        top_solute: str,
        top_output: str = "topol.top",
        output_name: str = "solvate"
) -> str:
    def extract_box_size_from_gro(gro_file):
        with open(gro_file, 'r') as f:
            lines = f.readlines()
            box_dimensions = lines[-1].split()
            return " ".join([str(dim) for dim in box_dimensions])

    box_size = extract_box_size_from_gro(gro_solvent)
    shutil.copy(top_solute, top_output)
    cmd = (
        f'gmx solvate -cp {gro_solute} -cs {gro_solvent} '
        f'-box {box_size} -p {top_output} -o {output_name}.gro '
    )
    return cmd


def gromacs_simulation_command(
        mdp: str,
        top: str,
        gro: str,
        name: str,
        n_max_warn: int = 10,
        n_threads: int = 1
) -> str:
    grompp_cmd = (
        f'gmx grompp -f {mdp} -p {top} -c {gro} -o {name}.tpr -maxwarn {n_max_warn}'
    )
    mdrun_cmd = f'gmx mdrun -nt {n_threads} -v -deffnm {name}'
    return f"{grompp_cmd} && {mdrun_cmd}"


def copy_files_to(
        files: list[str],
        destination_dir: str
) -> None:
    for file in files:
        logger.info(
            f"copying {file} to {destination_dir}"
        )
        shutil.copy(file, join(destination_dir, basename(file)))

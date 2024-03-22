import logging
import shutil

VOLUME_PER_CG_BEAD_IN_NM3 = 0.08

__all__ = [
    "generate_solvent_with_gromacs",
    "solvate_solute_command",
    "gromacs_simulation_command",
]

logger = logging.getLogger(__name__)


def generate_solvent_with_gromacs(
    gro_solvent_mol: str, box_length: float, output_name: str = "solvent", scale: float = 0.1
) -> str:
    box_size = " ".join([str(box_length) for _ in range(3)])
    return f"gmx solvate -cs {gro_solvent_mol} -box {box_size} -o {output_name}.gro -scale {scale}"


def solvate_solute_command(
    gro_solute: str, gro_solvent: str, top_solute: str, top_output: str = "topol.top", output_name: str = "solvate"
) -> str:
    def extract_box_size_from_gro(gro_file):
        with open(gro_file) as f:
            lines = f.readlines()
            box_dimensions = lines[-1].split()
            return " ".join([str(dim) for dim in box_dimensions])

    box_size = extract_box_size_from_gro(gro_solvent)
    shutil.copy(top_solute, top_output)
    cmd = f"gmx solvate -cp {gro_solute} -cs {gro_solvent} -box {box_size} -p {top_output} -o {output_name}.gro "
    return cmd


def gromacs_simulation_command(
    mdp: str, top: str, gro: str, name: str, n_max_warn: int = 10, n_threads: int = 1, verbose: bool = True
) -> str:
    grompp_cmd = f"gmx grompp -f {mdp} -p {top} -c {gro} -o {name}.tpr -po {name}_out.mdp -maxwarn {n_max_warn}"
    mdrun_cmd = f"gmx mdrun -nt {n_threads} -deffnm {name}"
    if verbose:
        mdrun_cmd += " -v"
    return f"{grompp_cmd} && {mdrun_cmd}"

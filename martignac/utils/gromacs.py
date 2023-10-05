import shutil
import logging
import regex
from MDAnalysis import Universe
from os.path import join, basename


VOLUME_PER_CG_BEAD_IN_NM3 = 0.08

__all__ = [
    "generate_solvent_with_gromacs",
    "generate_solvent_with_packmol",
    "convert_pdb_to_gro",
    "solvate_solute_command",
    "gromacs_simulation_command",
    "copy_files_to",
    "generate_top_file_for_generic_molecule",
    "sub_template_mdp"
]

logger = logging.getLogger(__name__)


def generate_solvent_with_gromacs(
        gro_solvent_mol: str,
        box_length: float,
        output_name: str = "solvent",
        scale: float = 0.1
) -> str:
    box_size = " ".join([str(box_length) for _ in range(3)])
    return (
        f'gmx solvate -cs {gro_solvent_mol} -box {box_size} '
        f'-o {output_name}.gro -scale {scale}'
    )


def generate_solvent_with_packmol(
        gro_solvent_mol: str,
        box_length: float,
        output_pdb: str,
        packmol_input_file: str = "packmol.inp"
) -> str:
    pdb_solvent_mol = f"{gro_solvent_mol.rstrip('.gro')}.pdb"
    universe = Universe(gro_solvent_mol)
    universe.atoms.write(pdb_solvent_mol)
    length_in_a = box_length * 10.
    num_beads_per_mol = len(universe.residues[0].atoms)
    volume_box = box_length ** 3
    num_beads_in_box = int(
            volume_box / (num_beads_per_mol * VOLUME_PER_CG_BEAD_IN_NM3)
    )
    logger.info(f"calling packmol to generate box with {num_beads_in_box} molecules")
    packmol_inp = f"""tolerance 2.0
filetype pdb
output {output_pdb}

structure {pdb_solvent_mol}
  number {num_beads_in_box}
  inside box 0. 0. 0. {length_in_a} {length_in_a} {length_in_a}
end structure
"""
    with open(packmol_input_file, "w") as pipe:
        pipe.write(packmol_inp)
    return "packmol < packmol.inp"


def convert_pdb_to_gro(pdb_file: str, output_gro: str, box_length: float) -> None:
    universe = Universe(pdb_file)
    box_size = [box_length * 10., box_length * 10., box_length * 10., 90., 90., 90.]
    universe.dimensions = box_size
    universe.atoms.write(output_gro)


def solvate_solute_command(
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
        try:
            shutil.copy(file, join(destination_dir, basename(file)))
        except shutil.SameFileError:
            logger.info(f"file {file} already located at {destination_dir}")
            pass


def generate_top_file_for_generic_molecule(
        molecule_name: str,
        force_field_filenames: list[str],
        top_filename: str,
        num_molecules: int = 1
) -> None:
    with open(top_filename, 'w') as f:
        # Include force field files
        for ff_file in force_field_filenames:
            f.write(f"#include \"{ff_file}\"\n")

        f.write("\n")
        f.write("[ system ]\n")
        f.write(f"{molecule_name} system\n\n")
        f.write("[ molecules ]\n")
        f.write(f"{molecule_name:4s}            {num_molecules:5d}\n")


def sub_template_mdp(mdp_source: str, template: str, new_entry: str, mdp_dest: str) -> None:
    with open(mdp_source) as pipe:
        mdp_content = pipe.read()
    mdp_content = regex.sub(template, new_entry, mdp_content)
    with open(mdp_dest, "w") as pipe:
        pipe.write(mdp_content)

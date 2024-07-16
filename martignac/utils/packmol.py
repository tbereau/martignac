import logging
from typing import Optional

from MDAnalysis import Universe

from martignac.utils.gromacs import VOLUME_PER_CG_BEAD_IN_NM3

__all__ = ["generate_solvent_with_packmol", "place_solute_in_solvent_with_packmol"]

logger = logging.getLogger(__name__)


def generate_solvent_with_packmol(
    gro_solvent_mol: str, box_length: float, output_pdb: str, packmol_input_file: str = "packmol.inp"
) -> str:
    """
    Generates a solvent box with specified dimensions using PACKMOL.

    This function converts a GROMACS .gro file of a solvent molecule to a .pdb file, calculates the number of molecules
    needed to fill a box of the given volume, and generates a PACKMOL input file to create a solvent box. The function
    then returns the command to run PACKMOL with the generated input file.

    Parameters:
        gro_solvent_mol (str): Path to the .gro file of the solvent molecule.
        box_length (float): Length of the cubic box in nm to be filled with solvent molecules.
        output_pdb (str): Path where the output .pdb file containing the solvent box should be saved.
        packmol_input_file (str, optional): Path to the PACKMOL input file to be generated. Defaults to "packmol.inp".

    Returns:
        str: The command to run PACKMOL with the generated input file.

    Note:
        PACKMOL must be installed and accessible from the command line to execute the returned command successfully.
    """
    pdb_solvent_mol = f"{gro_solvent_mol.rstrip('.gro')}.pdb"
    universe = Universe(gro_solvent_mol)
    universe.atoms.write(pdb_solvent_mol)
    length_in_a = box_length * 10.0
    num_beads_per_mol = len(universe.residues[0].atoms)
    volume_box = box_length**3
    num_beads_in_box = int(volume_box / (num_beads_per_mol * VOLUME_PER_CG_BEAD_IN_NM3))
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
    return f"packmol < {packmol_input_file}"


def place_solute_in_solvent_with_packmol(
    gro_solute: str,
    gro_solvent: str,
    output_pdb: str,
    restraint_along_x: Optional[float] = None,
    restraint_along_y: Optional[float] = None,
    restraint_along_z: Optional[float] = None,
    restraint_std_dev: float = 1.0,
    packmol_input_file: str = "packmol.inp",
) -> str:
    """
    Places a solute molecule within a solvent box using PACKMOL with optional positional restraints.

    This function prepares a PACKMOL input file to place a solute molecule within a pre-existing solvent box,
    allowing for the specification of positional restraints along the x, y, and z axes. This can be particularly
    useful for simulations where the solute needs to be positioned at a specific location within the solvent box,
    such as near a surface or interface. The function returns the command to run PACKMOL with the generated input file.

    Parameters:
        gro_solute (str): Path to the .gro file of the solute molecule.
        gro_solvent (str): Path to the .gro file of the solvent box.
        output_pdb (str): Path where the output .pdb file, containing both solute and solvent, should be saved.
        restraint_along_x (Optional[float]): Optional x-coordinate to restrain the solute molecule along the x-axis.
        restraint_along_y (Optional[float]): Optional y-coordinate to restrain the solute molecule along the y-axis.
        restraint_along_z (Optional[float]): Optional z-coordinate to restrain the solute molecule along the z-axis.
        restraint_std_dev (float): Standard deviation for the positional restraint, defining the allowed fluctuation range.
        packmol_input_file (str): Path to the PACKMOL input file to be generated. Defaults to "packmol.inp".

    Returns:
        str: The command to run PACKMOL with the generated input file.

    Note:
        PACKMOL must be installed and accessible from the command line to execute the returned command successfully.
    """
    logger.info("calling packmol to place solute in solvent")
    pdb_solute = f"{gro_solute.rstrip('.gro')}.pdb"
    universe = Universe(gro_solute)
    universe.atoms.write(pdb_solute)
    pdb_solvent = f"{gro_solvent.rstrip('.gro')}.pdb"
    universe = Universe(gro_solvent)
    universe.atoms.write(pdb_solvent)
    box_dimensions = universe.dimensions
    logger.info(f"box dimensions: {box_dimensions}")
    x_min, x_max = 0.0, box_dimensions[0]
    if restraint_along_x is not None:
        x_min = restraint_along_x - restraint_std_dev
        x_max = restraint_along_x + restraint_std_dev
    y_min, y_max = 0.0, box_dimensions[1]
    if restraint_along_y is not None:
        y_min = restraint_along_y - restraint_std_dev
        y_max = restraint_along_y + restraint_std_dev
    z_min, z_max = 0.0, box_dimensions[2]
    if restraint_along_z is not None:
        z_min = restraint_along_z - restraint_std_dev
        z_max = restraint_along_z + restraint_std_dev
    packmol_inp = f"""tolerance 2.0

filetype pdb

output {output_pdb}

structure {pdb_solvent}
  number 1
  fixed 0. 0. 0. 0. 0. 0.
end structure

structure {pdb_solute}
  number 1
  inside box {x_min} {y_min} {z_min} {x_max} {y_max} {z_max}
end structure
"""
    with open(packmol_input_file, "w") as pipe:
        pipe.write(packmol_inp)
    return f"packmol < {packmol_input_file}"

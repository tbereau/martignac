import logging

from MDAnalysis import Universe

from martignac.utils.gromacs import VOLUME_PER_CG_BEAD_IN_NM3

logger = logging.getLogger(__name__)


def generate_solvent_with_packmol(
    gro_solvent_mol: str, box_length: float, output_pdb: str, packmol_input_file: str = "packmol.inp"
) -> str:
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
    return "packmol < packmol.inp"

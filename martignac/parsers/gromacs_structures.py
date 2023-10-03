import MDAnalysis as mda

__all__ = ["get_number_of_molecules_from_gro"]


def get_number_of_molecules_from_gro(gro_filename: str) -> int:
    u = mda.Universe(gro_filename)
    return len(u.residues)
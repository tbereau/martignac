import MDAnalysis as mda

__all__ = ["get_number_of_molecules_from_gro"]


def get_number_of_molecules_from_gro(gro_filename: str) -> int:
    """
    Calculates the number of molecules in a GROMACS .gro file.

    This function utilizes MDAnalysis to load a .gro file and counts the number of residues, which corresponds
    to the number of molecules present in the file. It's useful for quickly assessing the composition of a
    simulation box without manually parsing the file.

    Parameters:
        gro_filename (str): The path to the .gro file from which to count molecules.

    Returns:
        int: The number of molecules found in the specified .gro file.
    """
    u = mda.Universe(gro_filename)
    return len(u.residues)

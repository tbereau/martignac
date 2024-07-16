import logging
from dataclasses import dataclass

import MDAnalysis as mda

from martignac.utils.misc import generate_top_file_for_generic_molecule

logger = logging.getLogger(__name__)


__all__ = [
    "ParticleType",
    "generate_top_file_for_particle",
    "generate_itp_file_for_particle",
    "generate_gro_for_particle_type",
    "find_particle_type_from_name",
]


@dataclass
class ParticleType:
    """
    Represents a particle type in a molecular dynamics simulation.

    Attributes:
        name (str): The name of the particle type.
        mass (float): The mass of the particle type.
        charge (float): The charge of the particle type.
        ptype (str): The type of the particle (e.g., atomistic, coarse-grained).
        c6 (float): The Lennard-Jones C6 parameter for the particle type.
        c12 (float): The Lennard-Jones C12 parameter for the particle type.
    """

    name: str
    mass: float
    charge: float
    ptype: str
    c6: float
    c12: float


def parse_particle_types_from_itp(filename: str) -> list[ParticleType]:
    """
    Parses particle types from a GROMACS .itp file.

    This function reads a GROMACS topology (.itp) file and extracts the definitions of particle types
    listed under the [ atomtypes ] section. Each particle type is represented by a `ParticleType` object,
    which includes properties such as name, mass, charge, particle type (ptype), and Lennard-Jones parameters (c6, c12).

    Parameters:
        filename (str): The path to the .itp file to be parsed.

    Returns:
        list[ParticleType]: A list of `ParticleType` objects extracted from the .itp file.
    """
    particle_types = []
    with open(filename) as file:
        lines = file.readlines()
        in_atomtypes = False
        for line in lines:
            line = line.strip()
            if line == "[ atomtypes ]":
                in_atomtypes = True
            elif in_atomtypes:
                if line.startswith(";") or not line:  # Skip comments and empty lines
                    continue
                try:
                    parts = line.split()
                    name = parts[0]
                    mass = float(parts[1])
                    charge = float(parts[2])
                    ptype = parts[3]
                    c6 = float(parts[4])
                    c12 = float(parts[5])
                    bead = ParticleType(name, mass, charge, ptype, c6, c12)
                    particle_types.append(bead)
                except (ValueError, IndexError):  # Exit the loop if line format doesn't match
                    break

    return particle_types


def find_particle_type_from_name(itp_filename: str, particle_name: str) -> ParticleType:
    """
    Finds and returns a ParticleType object from a GROMACS .itp file based on the particle name.

    This function searches through the particle types defined in a specified .itp file, looking for a particle
    type that matches the given name. It utilizes the `parse_particle_types_from_itp` function to read and parse
    the .itp file, extracting all defined particle types into a list of ParticleType objects. It then iterates
    through this list to find a particle type with a name that matches the `particle_name` parameter.

    Parameters:
        itp_filename (str): The path to the .itp file to be searched.
        particle_name (str): The name of the particle type to find.

    Returns:
        ParticleType: The ParticleType object corresponding to the specified particle name.

    Raises:
        StopIteration: If no particle type with the specified name is found in the .itp file.
    """
    particle_types = parse_particle_types_from_itp(itp_filename)
    return next(p for p in particle_types if p.name == particle_name)


def generate_gro_for_particle_type(particle_type: ParticleType, gro_filename: str, box_length: float = 100.0) -> None:
    """
    Generates a GROMACS .gro file for a single particle type with specified box dimensions.

    This function creates a .gro file representing a simulation box containing a single particle of the specified type.
    The particle is placed at the origin of the box. The function uses MDAnalysis to create an empty universe, sets the
    particle's attributes, and then writes the .gro file with the specified box dimensions.

    Parameters:
        particle_type (ParticleType): The particle type for which to generate the .gro file. This should be an instance
                                      of the ParticleType class, containing the necessary attributes like name, mass,
                                      and charge.
        gro_filename (str): The path and filename where the .gro file will be saved. If the file already exists, it will
                            be overwritten.
        box_length (float): The length of the sides of the cubic simulation box in which the particle is placed. The
                            default value is 100.0 Ångströms.

    Returns:
        None: This function does not return a value but writes directly to a file specified by `gro_filename`.
    """
    u = mda.Universe.empty(n_atoms=1, trajectory=True)

    # Setting properties
    u.add_TopologyAttr("name", values=[particle_type.name])
    u.add_TopologyAttr("type", values=[particle_type.ptype])
    u.add_TopologyAttr("resid", values=[1])
    u.add_TopologyAttr("resname", values=[particle_type.name])
    u.add_TopologyAttr("masses", values=[particle_type.mass])
    u.add_TopologyAttr("charges", values=[particle_type.charge])
    u.atoms.positions = [[0, 0, 0]]  # Setting it to the origin for simplicity
    box_size = [box_length, box_length, box_length, 90.0, 90.0, 90.0]
    u.dimensions = box_size

    logger.info(f"generating {gro_filename} for molecule {particle_type.name}")
    u.atoms.write(gro_filename)


def generate_top_file_for_particle(
    particle: ParticleType,
    force_field_filenames: list[str],
    top_filename: str,
    num_molecules: int = 1,
) -> None:
    """
    Generates a GROMACS topology (.top) file for a single particle type.

    This function creates a .top file for a specified particle type, incorporating the specified force field files.
    The generated .top file includes the definition of the particle type and the number of molecules of that type
    to be included in the simulation. It delegates the actual writing of the file to the
    `generate_top_file_for_generic_molecule` function, which handles the formatting and inclusion of force field
    parameters.

    Parameters:
        particle (ParticleType): The particle type for which to generate the .top file. This should be an instance
                                 of the ParticleType class, containing the necessary attributes like name, mass,
                                 and charge.
        force_field_filenames (list[str]): A list of paths to the force field files to be included in the .top file.
                                           These files contain the parameters necessary for the simulation, such as
                                           bond lengths, angles, and Lennard-Jones parameters.
        top_filename (str): The path and filename where the .top file will be saved. If the file already exists, it
                            will be overwritten.
        num_molecules (int): The number of molecules of the specified particle type to include in the simulation.
                             Defaults to 1.

    Returns:
        None: This function does not return a value but writes directly to a file specified by `top_filename`.
    """
    return generate_top_file_for_generic_molecule(particle.name, force_field_filenames, top_filename, num_molecules)


def generate_itp_file_for_particle(
    particle: ParticleType,
    itp_filename: str,
) -> None:
    """
    Generates a GROMACS .itp file for a single particle type.

    This function writes to a .itp file specific information about a particle type, including its name, type, and charge.
    The .itp (include topology) file format is used by GROMACS to define molecule types in simulations. This function
    creates a minimal .itp file containing a single [moleculetype] entry and one [atoms] entry for the specified particle.

    Parameters:
        particle (ParticleType): The particle type for which to generate the .itp file. This should be an instance
                                 of the ParticleType class, containing the necessary attributes like name, mass,
                                 and charge.
        itp_filename (str): The path and filename where the .itp file will be saved. If the file already exists, it
                            will be overwritten.

    Returns:
        None: This function does not return a value but writes directly to a file specified by `itp_filename`.
    """
    with open(itp_filename, "w") as f:
        f.write(f";;;;;; {particle.name} bead\n")
        f.write("[moleculetype]\n")
        f.write("; molname       nrexcl\n")
        f.write(f"  {particle.name}\t1\n\n")

        f.write("[atoms]\n")
        f.write("; id    type    resnr   residu  atom    cgnr    charge\n")
        f.write(f"  1\t{particle.name}\t1\t{particle.name}\t{particle.ptype}\t1\t{particle.charge}\n")

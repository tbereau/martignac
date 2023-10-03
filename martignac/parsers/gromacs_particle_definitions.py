from dataclasses import dataclass
import MDAnalysis as mda
import logging

from martignac.utils.gromacs import generate_top_file_for_generic_molecule

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
    name: str
    mass: float
    charge: float
    ptype: str
    c6: float
    c12: float


def parse_particle_types_from_itp(filename: str) -> list[ParticleType]:
    particle_types = []
    with open(filename, 'r') as file:
        lines = file.readlines()
        in_atomtypes = False
        for line in lines:
            line = line.strip()
            if line == "[ atomtypes ]":
                in_atomtypes = True
            elif in_atomtypes:
                if line.startswith(';') or not line:  # Skip comments and empty lines
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
    particle_types = parse_particle_types_from_itp(itp_filename)
    return next(p for p in particle_types if p.name == particle_name)


def generate_gro_for_particle_type(
        particle_type: ParticleType, gro_filename: str, box_length: float = 100.0
) -> None:
    u = mda.Universe.empty(n_atoms=1, trajectory=True)

    # Setting properties
    u.add_TopologyAttr('name', values=[particle_type.name])
    u.add_TopologyAttr('type', values=[particle_type.ptype])
    u.add_TopologyAttr('resid', values=[1])
    u.add_TopologyAttr('resname', values=[particle_type.name])
    u.add_TopologyAttr('masses', values=[particle_type.mass])
    u.add_TopologyAttr('charges', values=[particle_type.charge])
    u.atoms.positions = [[0, 0, 0]]  # Setting it to the origin for simplicity
    box_size = [box_length, box_length, box_length, 90., 90., 90.]
    u.dimensions = box_size

    logger.info(f"generating {gro_filename} for molecule {particle_type.name}")
    u.atoms.write(gro_filename)


def generate_top_file_for_particle(
        particle: ParticleType,
        force_field_filenames: list[str],
        top_filename: str,
        num_molecules: int = 1,
) -> None:
    return generate_top_file_for_generic_molecule(
        particle.name, force_field_filenames, top_filename, num_molecules
    )


def generate_itp_file_for_particle(
        particle: ParticleType,
        itp_filename: str,
) -> None:
    with open(itp_filename, 'w') as f:
        f.write(f";;;;;; {particle.name} bead\n")
        f.write("[moleculetype]\n")
        f.write("; molname       nrexcl\n")
        f.write(f"  {particle.name}\t1\n\n")

        f.write("[atoms]\n")
        f.write("; id    type    resnr   residu  atom    cgnr    charge\n")
        f.write(
            f"  1\t{particle.name}\t1\t{particle.name}\t{particle.ptype}\t1\t{particle.charge}\n"
        )

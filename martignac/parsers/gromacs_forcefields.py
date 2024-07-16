import logging
from contextlib import suppress
from dataclasses import dataclass
from random import uniform
from typing import Optional

import numpy as np
from MDAnalysis import Universe

from martignac.utils.misc import generate_top_file_for_generic_molecule

__all__ = [
    "Atom",
    "Bond",
    "Angle",
    "Constraint",
    "Molecule",
    "parse_molecules_from_itp",
    "find_molecule_from_name",
    "generate_gro_file_for_molecule",
    "generate_top_file_for_molecule",
    "get_molecule_from_name",
    "generate_itp_file_for_molecule",
]

logger = logging.getLogger(__name__)


@dataclass
class Atom:
    id: int
    type: str
    residue_number: int
    residue: str
    atom: str
    charge_number: int
    charge: int

    @classmethod
    def parse_from_itp_entry(cls, entry: list) -> "Atom":
        return Atom(
            int(entry[0]), str(entry[1]), int(entry[2]), str(entry[3]), str(entry[4]), int(entry[5]), int(entry[6])
        )


@dataclass
class Constraint:
    id_i: int
    id_j: int
    funct: int
    length: float

    @classmethod
    def parse_from_itp_entry(cls, entry: list) -> "Constraint":
        return Constraint(int(entry[0]), int(entry[1]), int(entry[2]), float(entry[3]))


@dataclass
class Bond(Constraint):
    force_constant: float

    @classmethod
    def parse_from_itp_entry(cls, entry: list) -> "Bond":
        return Bond(int(entry[0]), int(entry[1]), int(entry[2]), float(entry[3]), float(entry[4]))


@dataclass
class Angle:
    id_i: int
    id_j: int
    id_k: int
    funct: int
    angle: float
    force_constant: float

    @classmethod
    def parse_from_itp_entry(cls, entry: list) -> "Angle":
        return Angle(int(entry[0]), int(entry[1]), int(entry[2]), int(entry[3]), float(entry[4]), float(entry[5]))


@dataclass
class Molecule:
    name: str
    number_excl: int
    atoms: list[Atom]
    bonds: Optional[list[Bond]] = None
    angles: Optional[list[Angle]] = None
    constraints: Optional[list[Constraint]] = None

    @property
    def num_atoms(self) -> int:
        return len(self.atoms)

    def generate_coordinates(self, jitter: float = 0.1) -> np.ndarray:
        coordinates = np.array([[0.0, 0.0, 0.0] for _ in range(self.num_atoms)])
        placed_atoms = set()

        def add_jitter(coord):
            return [val + uniform(-jitter, jitter) for val in coord]

        def get_next_coordinate(id_ref, length):
            if id_ref - 1 not in placed_atoms:
                return add_jitter([length, 0, 0])
            else:
                return add_jitter([coordinates[id_ref - 1][0] + length, 0, 0])

        if self.constraints:
            for constraint in self.constraints:
                coordinates[constraint.id_j - 1] = get_next_coordinate(constraint.id_i, constraint.length)
                placed_atoms.add(constraint.id_j - 1)

        if self.bonds:
            for bond in self.bonds:
                if bond.id_j - 1 not in placed_atoms:
                    coordinates[bond.id_j - 1] = get_next_coordinate(bond.id_i, bond.length)
                    placed_atoms.add(bond.id_j - 1)

        return coordinates

    @classmethod
    def parse_from_itp_entry(cls, entry: dict) -> "Molecule":
        name = entry["moleculetype"][0][0]
        number_excl = entry["moleculetype"][0][1]
        atoms = [Atom.parse_from_itp_entry(e) for e in entry["atoms"]]
        bonds, angles, constraints = [], [], []
        with suppress(KeyError):
            bonds = [Bond.parse_from_itp_entry(e) for e in entry["bonds"]]
        with suppress(KeyError):
            angles = [Angle.parse_from_itp_entry(e) for e in entry["angles"]]
        with suppress(KeyError):
            constraints = [Constraint.parse_from_itp_entry(e) for e in entry["constraints"]]
        return Molecule(
            name=name, number_excl=number_excl, atoms=atoms, bonds=bonds, angles=angles, constraints=constraints
        )


def parse_molecules_from_itp(itp_filename: str) -> list[Molecule]:
    with open(itp_filename) as f:
        lines = f.readlines()

    molecules = []
    molecule_data = {}
    section = None

    for line in lines:
        line = line.strip()

        # Skip comments and empty lines
        if line.startswith(";") or not line:
            continue

        # Detect section headers
        if line.startswith("["):
            section = line.strip("[] ").lower()

            # If a new molecule starts
            if section == "moleculetype" and molecule_data:
                molecules.append(molecule_data)
                molecule_data = {}

            molecule_data[section] = []
            continue

        # Add data for the current section
        if section:
            molecule_data[section].append(line.split())

    # Add the last molecule
    if molecule_data:
        molecules.append(molecule_data)

    return [Molecule.parse_from_itp_entry(m) for m in molecules]


def find_molecule_from_name(itp_filenames: list[str], molecule_name: str) -> Molecule:
    molecules = []
    for itp_filename in itp_filenames:
        with suppress(KeyError):
            molecules = [*molecules, *parse_molecules_from_itp(itp_filename)]
    return next(m for m in molecules if m.name == molecule_name)


def generate_gro_file_for_molecule(molecule: Molecule, gro_filename: str, box_length: float = 100.0) -> None:
    n_atoms = len(molecule.atoms)
    dtype = [
        ("name", np.dtype("<U4")),
        ("type", np.dtype("<U4")),
        ("resid", int),
        ("resname", np.dtype("<U4")),
        ("charge", float),
        ("id", int),
    ]

    atom_data = np.zeros(n_atoms, dtype=dtype)

    for i, atom in enumerate(molecule.atoms):
        atom_data[i] = (atom.atom, atom.type, atom.residue_number, atom.residue, atom.charge, atom.id)

    assert len(set(atom_data["resid"])) == 1, "only one resid supported"
    assert len(set(atom_data["resname"])) == 1, "only one resname supported"
    residue_indices = [0] * n_atoms
    u = Universe.empty(n_atoms=n_atoms, n_residues=1, atom_resindex=residue_indices, trajectory=True)
    u.add_TopologyAttr("name", atom_data["name"])
    u.add_TopologyAttr("type", atom_data["type"])
    resid = atom_data["resid"][0] if len(atom_data["resid"]) > 1 else atom_data["resid"]
    resname = atom_data["resname"][0] if len(atom_data["resname"]) > 1 else atom_data["resname"]
    resid = [resid] if type(resid) not in [list, np.ndarray] else resid
    resname = [resname] if type(resname) not in [list, np.ndarray] else resname
    u.add_TopologyAttr("resid", resid)
    u.add_TopologyAttr("resname", resname)
    u.add_TopologyAttr("charge", atom_data["charge"])
    u.add_TopologyAttr("id", atom_data["id"])

    u.atoms.positions = [[x * 10 for x in coord] for coord in molecule.generate_coordinates()]
    if molecule.bonds:
        bond_tuples = [(bond.id_i - 1, bond.id_j - 1) for bond in molecule.bonds]  # adjust index by -1
        u.add_TopologyAttr("bonds", bond_tuples)

    box_size = [box_length, box_length, box_length, 90.0, 90.0, 90.0]
    u.dimensions = box_size

    logger.info(f"Generating {gro_filename} for molecule {molecule.name}")
    u.atoms.write(gro_filename)


def generate_top_file_for_molecule(
    molecule: Molecule, force_field_filenames: list[str], top_filename: str, num_molecules: int = 1
) -> None:
    return generate_top_file_for_generic_molecule(molecule.name, force_field_filenames, top_filename, num_molecules)


def _get_atom_from_string(atom_string: str, i: int, mol_name: str) -> Atom:
    charge = 0
    if "+" in atom_string:
        particle_info = atom_string.split("+")
        charge = int(particle_info[1]) if particle_info[1].isdigit() else 1
        atom_string = particle_info[0]
    elif "-" in atom_string:
        particle_info = atom_string.split("-")
        charge = -1 * int(particle_info[1]) if particle_info[1].isdigit() else -1
        atom_string = particle_info[0]
    atom = Atom(
        i + 1,  # id
        atom_string,  # type
        1,  # residue_number
        mol_name,  # residue
        atom_string,  # atom
        i + 1,  # charge_number
        charge,  # charge
    )
    return atom


def get_molecule_from_name(
    molecule_name: str,
    bond_length: float,
    bond_constant: Optional[float] = None,
    number_excl: int = 3,
    molecule_label: Optional[str] = None,
) -> Molecule:
    particle_names = molecule_name.split(",")[0].split()
    if molecule_label is None:
        molecule_label = "".join(particle_names)[:5]
    # Construct list of atoms
    atoms = [_get_atom_from_string(name, i, molecule_label) for i, name in enumerate(particle_names)]
    molecule = Molecule(molecule_label, number_excl, atoms)
    # Add bonds and constraints to the molecule
    if "," in molecule_name:
        # Obtain bonds and constraints from molecule string
        bond_constraints_info = molecule_name.split(",")[1].split()
        bond_info = [[int(idx) + 1 for idx in b.split("-")] for b in bond_constraints_info if "-" in b]
        constraint_info = [[int(idx) + 1 for idx in c.split("_")] for c in bond_constraints_info if "_" in c]
        # Add bonds and constraints to the molecule
        if bond_constant is None and len(bond_info) > 0:
            raise ValueError("bond_constant must be specified if the molecule contains bonds")
        if len(constraint_info) > 0:
            constraints = []
            for id_i, id_j in constraint_info:
                constraint = Constraint(id_i, id_j, 1, bond_length)
                constraints.append(constraint)
            molecule.constraints = constraints
        if len(bond_info) > 0:
            bonds = []
            for id_i, id_j in bond_info:
                bond = Bond(id_i, id_j, 1, bond_length, bond_constant)
                bonds.append(bond)
            molecule.bonds = bonds
    return molecule


def generate_itp_file_for_molecule(
    molecule: Molecule,
    itp_filename: str,
) -> None:
    with open(itp_filename, "w") as f:
        f.write(f"\n;;;;;; {molecule.name} molecule\n")
        f.write("\n[moleculetype]\n")
        f.write("; molname       nrexcl\n")
        f.write(f" {molecule.name}\t{molecule.number_excl}\n")
        f.write("\n[atoms]\n")
        f.write("; id    type    resnr   residu  atom    cgnr    charge\n")
        for atom in molecule.atoms:
            f.write(
                f" {atom.id}\t{atom.type}\t{atom.residue_number}\t{atom.residue}\t"
                f"{atom.atom}\t{atom.charge_number}\t{atom.charge}\n"
            )
        if molecule.bonds is not None:
            f.write("\n[bonds]\n")
            f.write("; i     j       funct   length  force_constant\n")
            for bond in molecule.bonds:
                f.write(f" {bond.id_i}\t{bond.id_j}\t{bond.funct}\t{bond.length}\t{bond.force_constant}\n")
        if molecule.angles is not None:
            f.write("\n[angles]\n")
            f.write("; i     j       k       funct   angle   force_constant\n")
            for angle in molecule.angles:
                f.write(
                    f" {angle.id_i}\t{angle.id_j}\t{angle.id_k}\t{angle.funct}\t{angle.angle}\t{angle.force_constant}\n"
                )
        if molecule.constraints is not None:
            f.write("\n[constraints]\n")
            f.write("; i     j       funct   length\n")
            for constraint in molecule.constraints:
                f.write(f" {constraint.id_i}\t{constraint.id_j}\t{constraint.funct}\t{constraint.length}\n")

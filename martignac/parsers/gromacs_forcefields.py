import logging
import re
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
    """
    Represents an atom within a molecule.

    This class models an atom's properties, including its unique identifier, type, residue information, and charge.
    It also provides a class method to parse atom information from a topology file entry.

    Attributes:
        id (int): Unique identifier of the atom within the molecule.
        type (str): Type of the atom, defining its chemical characteristics.
        residue_number (int): Identifier of the residue the atom belongs to.
        residue (str): Name of the residue the atom belongs to.
        atom (str): Name of the atom.
        charge_number (int): Identifier for the charge group the atom belongs to.
        charge (int): The charge of the atom.

    Class Methods:
        parse_from_itp_entry(cls, entry: list) -> "Atom": Parses atom information from a given topology file entry and returns an Atom instance.
    """

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
            int(entry[0]),
            str(entry[1]),
            int(entry[2]),
            str(entry[3]),
            str(entry[4]),
            int(entry[5]),
            int(entry[6]),
        )

    @property
    def is_charged(self) -> bool:
        return self.charge != 0


@dataclass
class Constraint:
    """
    Represents a constraint between two atoms in a molecule.

    This class models a constraint, which is a fixed distance between two atoms, typically used to maintain
    a certain structure within the molecule. Constraints are often used in molecular dynamics simulations
    to simplify the model or enforce certain conditions.

    Attributes:
        id_i (int): The identifier of the first atom in the constraint.
        id_j (int): The identifier of the second atom in the constraint.
        funct (int): The function type of the constraint. In GROMACS, this typically refers to the type of constraint algorithm used.
        length (float): The length of the constraint, usually in nanometers.

    Class Methods:
        parse_from_itp_entry(cls, entry: list) -> "Constraint": Parses a constraint from a given entry in a topology file.
    """

    id_i: int
    id_j: int
    funct: int
    length: float

    @classmethod
    def parse_from_itp_entry(cls, entry: list) -> "Constraint":
        return Constraint(int(entry[0]), int(entry[1]), int(entry[2]), float(entry[3]))


@dataclass
class Bond(Constraint):
    """
    Represents a bond between two atoms in a molecule, extending the Constraint class.

    This class models a bond as a specialized type of constraint with an additional attribute for the force constant.
    Bonds are used to define the fixed distance between two atoms along with the force constant that describes the
    strength of the bond. This is particularly useful in molecular dynamics simulations for defining the interactions
    and energy calculations between atoms.

    Attributes:
        id_i (int): The identifier of the first atom in the bond.
        id_j (int): The identifier of the second atom in the bond.
        funct (int): The function type of the bond. In GROMACS, this typically refers to the type of bond algorithm used.
        length (float): The length of the bond, usually in nanometers.
        force_constant (float): The force constant of the bond, describing its strength.

    Class Methods:
        parse_from_itp_entry(cls, entry: list) -> "Bond": Parses a bond from a given entry in a topology file, extending the Constraint class method with an additional parameter for the force constant.
    """

    force_constant: float

    @classmethod
    def parse_from_itp_entry(cls, entry: list) -> "Bond":
        return Bond(
            int(entry[0]),
            int(entry[1]),
            int(entry[2]),
            float(entry[3]),
            float(entry[4]),
        )


@dataclass
class Angle:
    """
    Represents an angle formed by three atoms in a molecule.

    This class models an angle, which is defined by three atoms, typically used to represent the geometric
    structure of a molecule. Angles are crucial in molecular dynamics simulations for defining the spatial
    arrangement of atoms and calculating potential energy based on geometric constraints.

    Attributes:
        id_i (int): The identifier of the first atom forming the angle.
        id_j (int): The identifier of the second atom (vertex) forming the angle.
        id_k (int): The identifier of the third atom forming the angle.
        funct (int): The function type of the angle. In GROMACS, this typically refers to the type of angle potential used.
        angle (float): The value of the angle, usually in degrees.
        force_constant (float): The force constant of the angle, describing its stiffness.

    Class Methods:
        parse_from_itp_entry(cls, entry: list) -> "Angle": Parses angle information from a given entry in a topology file.
    """

    id_i: int
    id_j: int
    id_k: int
    funct: int
    angle: float
    force_constant: float

    @classmethod
    def parse_from_itp_entry(cls, entry: list) -> "Angle":
        return Angle(
            int(entry[0]),
            int(entry[1]),
            int(entry[2]),
            int(entry[3]),
            float(entry[4]),
            float(entry[5]),
        )


@dataclass
class Molecule:
    """
    Represents a molecule with its constituent atoms, bonds, angles, and constraints.

    This class encapsulates a molecule's structure, including its atoms and the relationships between them
    (bonds, angles, and constraints). It provides methods for generating molecular coordinates, parsing molecule
    data from topology files, and generating files for molecular dynamics simulations.

    Attributes:
        name (str): The name of the molecule.
        number_excl (int): The number of exclusions for the molecule, used in simulations to define non-bonded interactions.
        atoms (list[Atom]): A list of `Atom` instances representing the atoms in the molecule.
        bonds (Optional[list[Bond]]): A list of `Bond` instances representing the bonds in the molecule. Default is None.
        angles (Optional[list[Angle]]): A list of `Angle` instances representing the angles in the molecule. Default is None.
        constraints (Optional[list[Constraint]]): A list of `Constraint` instances representing the constraints in the molecule. Default is None.

    Methods:
        num_atoms (property): Returns the number of atoms in the molecule.
        generate_coordinates(self, jitter: float = 0.1) -> np.ndarray: Generates 3D coordinates for the atoms in the molecule.
        parse_from_itp_entry(cls, entry: dict) -> "Molecule": Class method to parse a molecule from a topology file entry.
    """

    name: str
    number_excl: int
    atoms: list[Atom]
    bonds: Optional[list[Bond]] = None
    angles: Optional[list[Angle]] = None
    constraints: Optional[list[Constraint]] = None

    @property
    def num_atoms(self) -> int:
        return len(self.atoms)

    @property
    def has_charged_beads(self) -> bool:
        return any(atom.is_charged for atom in self.atoms)

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
                coordinates[constraint.id_j - 1] = get_next_coordinate(
                    constraint.id_i, constraint.length
                )
                placed_atoms.add(constraint.id_j - 1)

        if self.bonds:
            for bond in self.bonds:
                if bond.id_j - 1 not in placed_atoms:
                    coordinates[bond.id_j - 1] = get_next_coordinate(
                        bond.id_i, bond.length
                    )
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
            constraints = [
                Constraint.parse_from_itp_entry(e) for e in entry["constraints"]
            ]
        return Molecule(
            name=name,
            number_excl=number_excl,
            atoms=atoms,
            bonds=bonds,
            angles=angles,
            constraints=constraints,
        )


def parse_molecules_from_itp(itp_filename: str) -> list[Molecule]:
    """
    Parses molecules from a GROMACS topology (.itp) file.

    This function reads a .itp file, extracts molecule information, and creates a list of Molecule instances
    based on the data found in the file. It handles multiple molecules within the same file, separating them
    based on the [moleculetype] directive in the .itp file format. Each molecule's data is parsed and used to
    instantiate a Molecule object, which is then added to the list of molecules to be returned.

    Parameters:
        itp_filename (str): The path to the .itp file to be parsed.

    Returns:
        list[Molecule]: A list of Molecule instances representing each molecule found in the .itp file.
    """
    with open(itp_filename) as f:
        content = f.read()

    # Pattern to find and process conditional compilation blocks
    conditional_block_pattern = re.compile(
        r"#ifdef FLEXIBLE\s*\[bonds\]\s*(.*?)#else\s*\[constraints\]\s*(.*?)#endif",
        re.DOTALL,
    )

    # Function to keep the 'bonds' part and discard the 'constraints' part
    def choose_bonds_over_constraints(match):
        bonds_part = match.group(1)  # The 'bonds' section of the match
        return f"[bonds]{bonds_part}"

    # Replace the matched sections with only the 'bonds' part
    modified_content = re.sub(
        conditional_block_pattern, choose_bonds_over_constraints, content
    )

    # Split the modified content into lines for further processing
    lines = modified_content.splitlines()

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
    """
    Searches for and returns a Molecule instance matching the specified name from a list of .itp files.

    This function iterates over a list of GROMACS topology (.itp) files, parsing each to find a molecule
    that matches the given name. It leverages the `parse_molecules_from_itp` function to parse the .itp files
    and then searches through the resulting list of Molecule instances for a name match. The first matching
    Molecule instance found is returned. If no match is found after all files have been searched, a StopIteration
    exception is raised.

    Parameters:
        itp_filenames (list[str]): A list of paths to .itp files to be searched.
        molecule_name (str): The name of the molecule to find.

    Returns:
        Molecule: The first Molecule instance found with a name matching `molecule_name`.

    Raises:
        StopIteration: If no molecule with the specified name is found in the provided .itp files.
    """
    molecules = []
    for itp_filename in itp_filenames:
        with suppress(KeyError):
            molecules = [*molecules, *parse_molecules_from_itp(itp_filename)]
    return next(m for m in molecules if m.name == molecule_name)


def generate_gro_file_for_molecule(
    molecule: Molecule, gro_filename: str, box_length: float = 100.0
) -> None:
    """
    Generates a GROMACS .gro file for a given molecule.

    This function creates a .gro file, which is a GROMACS file format used to describe the positions of atoms in a molecule.
    It sets up a simulation box with a specified box length and places the atoms of the molecule within this box. The positions
    of the atoms are determined by the `generate_coordinates` method of the `Molecule` class, and additional molecular properties
    such as bonds are also considered if present.

    Parameters:
        molecule (Molecule): The molecule for which to generate the .gro file.
        gro_filename (str): The path and name of the .gro file to be generated.
        box_length (float): The length of the sides of the cubic simulation box in which the molecule is placed. Default is 100.0.

    Returns:
        None: This function does not return a value but writes directly to a file.
    """
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
        atom_data[i] = (
            atom.atom,
            atom.type,
            atom.residue_number,
            atom.residue,
            atom.charge,
            atom.id,
        )

    assert len(set(atom_data["resid"])) == 1, "only one resid supported"
    assert len(set(atom_data["resname"])) == 1, "only one resname supported"
    residue_indices = [0] * n_atoms
    u = Universe.empty(
        n_atoms=n_atoms, n_residues=1, atom_resindex=residue_indices, trajectory=True
    )
    u.add_TopologyAttr("name", atom_data["name"])
    u.add_TopologyAttr("type", atom_data["type"])
    resid = atom_data["resid"][0] if len(atom_data["resid"]) > 1 else atom_data["resid"]
    resname = (
        atom_data["resname"][0]
        if len(atom_data["resname"]) > 1
        else atom_data["resname"]
    )
    resid = [resid] if type(resid) not in [list, np.ndarray] else resid
    resname = [resname] if type(resname) not in [list, np.ndarray] else resname
    u.add_TopologyAttr("resid", resid)
    u.add_TopologyAttr("resname", resname)
    u.add_TopologyAttr("charge", atom_data["charge"])
    u.add_TopologyAttr("id", atom_data["id"])

    u.atoms.positions = [
        [x * 10 for x in coord] for coord in molecule.generate_coordinates()
    ]
    if molecule.bonds:
        bond_tuples = [
            (bond.id_i - 1, bond.id_j - 1) for bond in molecule.bonds
        ]  # adjust index by -1
        u.add_TopologyAttr("bonds", bond_tuples)

    box_size = [box_length, box_length, box_length, 90.0, 90.0, 90.0]
    u.dimensions = box_size

    logger.info(f"Generating {gro_filename} for molecule {molecule.name}")
    u.atoms.write(gro_filename)


def generate_top_file_for_molecule(
    molecule: Molecule,
    force_field_filenames: list[str],
    top_filename: str,
    num_molecules: int = 1,
) -> None:
    """
    Generates a GROMACS topology (.top) file for a given molecule using specified force fields.

    This function delegates the task of generating a .top file to the `generate_top_file_for_generic_molecule`
    utility function. It prepares the necessary parameters, including the molecule's name, the list of force field
    filenames, the target .top file name, and the number of molecules to be included in the simulation.

    Parameters:
        molecule (Molecule): The molecule instance for which to generate the .top file.
        force_field_filenames (list[str]): A list of strings representing the paths to the force field files to be used.
        top_filename (str): The path and name of the .top file to be generated.
        num_molecules (int): The number of molecules to be included in the .top file. Default is 1.

    Returns:
        None: This function does not return a value but writes directly to a file.
    """
    return generate_top_file_for_generic_molecule(
        molecule.name, force_field_filenames, top_filename, num_molecules
    )


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
    """
    Constructs a Molecule instance from a simplified molecule name string.

    This function allows for the creation of a Molecule object by specifying a string that represents the molecule's
    composition, bond lengths, and optionally, bond constants. It supports the definition of atoms, bonds, and constraints
    within the molecule based on the provided string format. The molecule name string format should include atom types,
    followed by bond information (if any), separated by commas.

    Parameters:
        molecule_name (str): A string representing the molecule's composition and structure. Atom types should be
                             separated by spaces, and bond information (if any) should follow after a comma. Bonds are
                             indicated by indices of atoms (starting from 1) connected by a dash ("-") for bonds or an
                             underscore ("_") for constraints.
        bond_length (float): The default length for bonds and constraints within the molecule, in nanometers.
        bond_constant (Optional[float]): The force constant for the bonds within the molecule. This parameter is required
                                         if the molecule_name string includes bond information. Default is None.
        number_excl (int): The number of exclusions for the molecule, used in simulations to define non-bonded interactions.
                           Default is 3.
        molecule_label (Optional[str]): An optional label for the molecule. If not provided, a label is generated from the
                                        first five characters of the concatenated atom types. Default is None.

    Returns:
        Molecule: An instance of the Molecule class, constructed based on the provided parameters.

    Raises:
        ValueError: If bond_constant is None but the molecule_name string includes bond information.
    """
    particle_names = molecule_name.split(",")[0].split()
    if molecule_label is None:
        molecule_label = "".join(particle_names)[:4]
    # Construct list of atoms
    atoms = [
        _get_atom_from_string(name, i, molecule_label)
        for i, name in enumerate(particle_names)
    ]
    molecule = Molecule(molecule_label, number_excl, atoms)
    # Add bonds and constraints to the molecule
    if "," in molecule_name:
        # Obtain bonds and constraints from molecule string
        bond_constraints_info = molecule_name.split(",")[1].split()
        bond_info = [
            [int(idx) + 1 for idx in b.split("-")]
            for b in bond_constraints_info
            if "-" in b
        ]
        constraint_info = [
            [int(idx) + 1 for idx in c.split("_")]
            for c in bond_constraints_info
            if "_" in c
        ]
        # Add bonds and constraints to the molecule
        if bond_constant is None and len(bond_info) > 0:
            raise ValueError(
                "bond_constant must be specified if the molecule contains bonds"
            )
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
    """
    Generates a GROMACS topology (.itp) file for a specified molecule.

    This function writes a .itp file containing the definition of a molecule, including atoms, bonds, angles,
    and constraints. The file format adheres to the GROMACS topology file standards, making it suitable for
    use in molecular dynamics simulations. The function outputs sections for moleculetype, atoms, bonds (if any),
    angles (if any), and constraints (if any), providing a comprehensive description of the molecule's structure.

    Parameters:
        molecule (Molecule): The molecule instance for which to generate the .itp file. This object should contain
                             all necessary information about the molecule, including its atoms, bonds, angles, and
                             constraints.
        itp_filename (str): The path and name of the .itp file to be generated. If the file already exists, it will
                            be overwritten.

    Returns:
        None: This function does not return a value but writes directly to a file specified by `itp_filename`.
    """
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
                f.write(
                    f" {bond.id_i}\t{bond.id_j}\t{bond.funct}\t{bond.length}\t{bond.force_constant}\n"
                )
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
                f.write(
                    f" {constraint.id_i}\t{constraint.id_j}\t{constraint.funct}\t{constraint.length}\n"
                )

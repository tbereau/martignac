import logging
from copy import copy
from dataclasses import dataclass, field

import MDAnalysis as mda

logger = logging.getLogger(__name__)


@dataclass
class MoleculeTopology:
    """
    Represents the topology of a molecule within a GROMACS simulation.

    Attributes:
        name (str): The name of the molecule.
        number_elements (int): The number of elements (atoms) in the molecule.
    """

    name: str
    number_elements: int


@dataclass
class Topology:
    """
    Represents the overall topology of a molecular system in a GROMACS simulation.

    This class encapsulates the system's name, the molecules present within the system, and any include files
    necessary for the simulation. It provides methods to parse topology files (.top), update molecule counts based
    on a given .gro file, and output a new topology file.

    Attributes:
        system (str): A descriptive name for the system.
        molecules (list[MoleculeTopology]): A list of `MoleculeTopology` instances representing the molecules
                                            present in the system.
        includes (list[str]): A list of paths to additional topology files to be included in the simulation.

    Methods:
        parse_top_file(cls, top_filename: str) -> "Topology":
            Class method to parse a GROMACS topology file (.top) and return a `Topology` instance representing
            the parsed data.

        update_counts_against_gro(self, gro_filename: str) -> None:
            Updates the molecule counts in the topology based on the actual counts found in a .gro file.

        output_top(self, filename: str) -> None:
            Writes the topology to a new .top file, including all molecules, system name, and include files.
    """

    system: str = ""
    molecules: list[MoleculeTopology] = field(default_factory=list)
    includes: list[str] = field(default_factory=list)

    @classmethod
    def parse_top_file(cls, top_filename: str) -> "Topology":
        with open(top_filename) as file:
            lines = file.readlines()

        includes = []
        system = ""
        molecules = []

        # Identify sections
        in_system_section = False
        in_molecules_section = False

        for line in lines:
            line = line.strip()

            # Skip comments and empty lines
            if not line or line.startswith(";"):
                continue

            # Identify the start of a section
            if line.startswith("["):
                section = line.strip("[]").strip().lower()
                if section == "system":
                    in_system_section = True
                    in_molecules_section = False
                elif section == "molecules":
                    in_molecules_section = True
                    in_system_section = False
                continue

            # Process the sections
            if in_system_section:
                system = line
            elif in_molecules_section:
                molecule_name, molecule_count = line.split()
                molecules.append(MoleculeTopology(molecule_name, int(molecule_count)))
            elif line.startswith("#include"):
                includes.append(line.split('"')[1])

        return Topology(includes=includes, system=system, molecules=molecules)

    def update_counts_against_gro(self, gro_filename: str) -> None:
        u = mda.Universe(gro_filename)
        existing_mol_names = []
        final_list_of_molecules = []
        for molecule in self.molecules:
            if molecule.name not in existing_mol_names:
                atom_selection = u.select_atoms(f"resname {molecule.name}")
                molecule.number_elements = atom_selection.residues.n_residues
                logger.info(
                    f"found {molecule.number_elements} {molecule.name} in {gro_filename}"
                )
                existing_mol_names.append(molecule.name)
                final_list_of_molecules.append(molecule)
        self.molecules = final_list_of_molecules

    def output_top(self, filename: str) -> None:
        with open(filename, "w") as file:
            # Write includes
            for include in self.includes:
                file.write(f'#include "{include}"\n')
            file.write("\n")

            # Write system
            file.write("[ system ]\n")
            file.write(f"{self.system}\n")
            file.write("\n")

            # Write molecules
            file.write("[ molecules ]\n")
            for molecule in self.molecules:
                file.write(f"{molecule.name:<20}{molecule.number_elements}\n")


def append_all_includes_to_top(main_top: Topology, others: list[Topology]) -> Topology:
    """
    Combines the include files from multiple Topology instances into a single Topology instance.

    This function takes a main Topology instance and a list of other Topology instances, then combines all unique
    include file paths from these instances into the main Topology instance. The main Topology's include list is
    first cleared to ensure no duplicates from its original includes, and then it is populated with unique includes
    from all provided Topology instances.

    Parameters:
        main_top (Topology): The main Topology instance to which all unique include file paths from other Topology
                             instances will be added.
        others (list[Topology]): A list of other Topology instances from which to gather all unique include file paths.

    Returns:
        Topology: The updated main Topology instance with a combined list of unique include file paths.
    """
    new_top = copy(main_top)
    new_top.includes = []
    for top in others:
        for i in top.includes:
            if i not in new_top.includes:
                new_top.includes.append(i)
    return new_top


def combine_multiple_topology_files(
    topology_files: list[str], system_name: str
) -> Topology:
    """
    Combines multiple GROMACS topology files into a single Topology instance.

    This function parses multiple topology (.top) files, each representing a part of a molecular system, and combines
    them into a single Topology instance. It aggregates all molecules and include files from each topology file,
    ensuring no duplicates in the include files. The system name for the combined topology is specified by the user.

    Parameters:
        topology_files (list[str]): A list of paths to the topology files to be combined.
        system_name (str): The name to assign to the combined system.

    Returns:
        Topology: An instance of the Topology class representing the combined topology of all input files.
    """
    comb_topology = Topology(system_name, molecules=[], includes=[])
    for top_file in topology_files:
        top = Topology.parse_top_file(top_file)
        comb_topology.molecules.extend(top.molecules)
        for include in top.includes:
            if include not in comb_topology.includes:
                comb_topology.includes.append(include)
    return comb_topology

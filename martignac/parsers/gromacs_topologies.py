from dataclasses import dataclass, field
from copy import copy


@dataclass
class MoleculeTopology:
    name: str
    number_elements: int


@dataclass
class Topology:
    system: str = ""
    molecules: list[MoleculeTopology] = field(default_factory=list)
    includes: list[str] = field(default_factory=list)

    @classmethod
    def parse_top_file(cls, top_filename: str) -> "Topology":
        with open(top_filename, 'r') as file:
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

    def output_top(self, filename: str) -> None:
        with open(filename, 'w') as file:
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
    new_top = copy(main_top)
    new_top.includes = []
    for top in others:
        for i in top.includes:
            if i not in new_top.includes:
                new_top.includes.append(i)
    return new_top

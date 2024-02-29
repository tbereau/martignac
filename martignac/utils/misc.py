import logging
import os
import shutil
from os.path import abspath, basename, join

import regex
from MDAnalysis import Universe

logger = logging.getLogger(__name__)


def copy_files_to(files: list[str], destination_dir: str) -> None:
    for file in files:
        logger.info(f"copying {file} to {destination_dir}")
        try:
            shutil.copy(file, join(destination_dir, basename(file)))
        except shutil.SameFileError:
            logger.info(f"file {file} already located at {destination_dir}")
            pass


def generate_top_file_for_generic_molecule(
    molecule_name: str, force_field_filenames: list[str], top_filename: str, num_molecules: int = 1
) -> None:
    with open(top_filename, "w") as f:
        for ff_file in force_field_filenames:
            f.write(f'#include "{ff_file}"\n')

        f.write("\n")
        f.write("[ system ]\n")
        f.write(f"{molecule_name} system\n\n")
        f.write("[ molecules ]\n")
        f.write(f"{molecule_name:4s}            {num_molecules:5d}\n")


def sub_template_mdp(mdp_source: str, template: str, new_entry: str, mdp_dest: str) -> None:
    with open(mdp_source) as pipe:
        mdp_content = pipe.read()
    mdp_content = regex.sub(template, new_entry, mdp_content)
    with open(mdp_dest, "w") as pipe:
        pipe.write(mdp_content)


def convert_pdb_to_gro(pdb_file: str, output_gro: str, box_length: float) -> None:
    universe = Universe(pdb_file)
    box_size = [box_length * 10.0, box_length * 10.0, box_length * 10.0, 90.0, 90.0, 90.0]
    universe.dimensions = box_size
    universe.atoms.write(output_gro)


def zip_directory(directory_name: str, output_file_name: str) -> str:
    # Determine the parent directory to initially save the archive
    parent_dir = abspath(join(directory_name, os.pardir))
    intermediate_output_path = os.path.join(parent_dir, os.path.basename(output_file_name))

    logger.info(f"Zipping {directory_name} to {intermediate_output_path}.zip")
    archive_path = shutil.make_archive(intermediate_output_path, "zip", directory_name)

    # Move the archive to the original directory
    final_output_path = os.path.join(directory_name, os.path.basename(output_file_name) + ".zip")
    shutil.move(archive_path, final_output_path)

    logger.info(f"Moved archive to {final_output_path}")
    return final_output_path

import inspect
import logging
import os
import shutil
from collections.abc import Mapping
from os.path import abspath, basename, islink, join
from zipfile import ZIP_DEFLATED, ZipFile

import numpy as np
import regex
from MDAnalysis import Universe

logger = logging.getLogger(__name__)


def func_name():
    return inspect.stack()[1][3]


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


def zip_directories(directory_names: list, output_file_name: str) -> str:
    if not directory_names:
        raise ValueError("directory_names list cannot be empty")

    # Assume all directories are at the same parent level for the output path
    parent_dir = abspath(join(directory_names[0], os.pardir))
    intermediate_output_path = join(parent_dir, basename(output_file_name)) + ".zip"

    logger.info(f"Zipping directories {directory_names} to {intermediate_output_path}")

    with ZipFile(intermediate_output_path, "w", ZIP_DEFLATED) as zipf:
        for directory_name in directory_names:
            for root, _, files in os.walk(directory_name):
                for file in files:
                    file_path = join(root, file)
                    # Exclude symlinks
                    if not islink(file_path):
                        # Calculate archive name based on the file's relative path to the directory being zipped
                        arcname = os.path.relpath(file_path, start=os.path.commonpath(directory_names))
                        zipf.write(file_path, arcname=arcname)

    # The final output path could be in the same directory as the script or a specific directory
    final_output_path = join(os.getcwd(), basename(output_file_name) + ".zip")
    os.rename(intermediate_output_path, final_output_path)

    logger.info(f"Moved archive to {final_output_path}")
    return final_output_path


def update_nested_dict(d, u):
    for k, v in u.items():
        if isinstance(v, Mapping):
            d[k] = update_nested_dict(d.get(k, {}), v)
        else:
            d[k] = v
    return d


def calculate_average_com(gro_file_path: str, molecule_names: list[str]) -> np.ndarray:
    u = Universe(gro_file_path)
    if not molecule_names:
        molecule_atoms = u.atoms
    else:
        selection_str = " or ".join([f"resname {mol}" for mol in molecule_names])
        molecule_atoms = u.select_atoms(selection_str)

    if len(molecule_atoms) > 0:
        return molecule_atoms.center_of_mass()
    else:
        raise ValueError(f"No atoms found for the types {molecule_names} in {gro_file_path}")


def translate_gro_by_vector(gro_file_path, output_gro_file_path, com_diff_vector):
    u = Universe(gro_file_path)
    u.atoms.positions += com_diff_vector
    u.atoms.write(output_gro_file_path)

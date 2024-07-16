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
    """
    Copies a list of files to a specified destination directory.

    This function iterates over a list of file paths, copying each file to the destination directory specified.
    If a file already exists at the destination, a message is logged, and the file is not copied again to avoid
    overwriting. This is useful for organizing files or preparing directories with necessary files without manual
    copying.

    Parameters:
        files (list[str]): A list of file paths to be copied.
        destination_dir (str): The directory path where the files should be copied to.

    Returns:
        None: This function does not return a value but copies files and logs actions.
    """
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
    """
    Generates a GROMACS topology file for a generic molecule.

    This function creates a topology file (.top) for a specified molecule using a list of force field filenames.
    It writes the necessary include statements for the force field files, defines the system, and specifies the
    number of molecules present. This is particularly useful for setting up molecular dynamics simulations where
    custom molecules or force fields are involved.

    Parameters:
        molecule_name (str): The name of the molecule for which the topology is being generated.
        force_field_filenames (list[str]): A list of filenames for the force field files to be included in the topology.
        top_filename (str): The path and filename where the topology file will be saved.
        num_molecules (int, optional): The number of molecules of this type present in the system. Defaults to 1.

    Returns:
        None: This function does not return a value but writes a .top file at the specified location.
    """
    with open(top_filename, "w") as f:
        for ff_file in force_field_filenames:
            f.write(f'#include "{ff_file}"\n')

        f.write("\n")
        f.write("[ system ]\n")
        f.write(f"{molecule_name} system\n\n")
        f.write("[ molecules ]\n")
        f.write(f"{molecule_name:4s}            {num_molecules:5d}\n")


def sub_template_mdp(mdp_source: str, template: str, new_entry: str, mdp_dest: str) -> None:
    """
    Substitutes a template placeholder in a GROMACS MDP file with a new entry.

    This function reads a GROMACS MDP (molecular dynamics parameters) file, searches for a specified template
    placeholder, and replaces it with a new entry. The modified MDP content is then written to a new destination file.
    This is useful for dynamically adjusting simulation parameters without manually editing MDP files.

    Parameters:
        mdp_source (str): The path to the source MDP file containing the template placeholder.
        template (str): The template placeholder to be replaced. This should be a regex pattern.
        new_entry (str): The new entry to replace the template placeholder with.
        mdp_dest (str): The path where the modified MDP file will be saved.

    Returns:
        None: This function does not return a value but writes the modified MDP content to a new file.
    """
    with open(mdp_source) as pipe:
        mdp_content = pipe.read()
    mdp_content = regex.sub(template, new_entry, mdp_content)
    with open(mdp_dest, "w") as pipe:
        pipe.write(mdp_content)


def convert_pdb_to_gro(pdb_file: str, output_gro: str, box_vector: np.ndarray) -> None:
    """
    Converts a PDB file to a GROMACS GRO file with specified box dimensions.

    This function utilizes MDAnalysis to load a PDB file, set the unit cell dimensions to the provided box vector,
    and then save the structure in GROMACS GRO format. This is useful for preparing molecular dynamics simulations
    where a specific box size is required.

    Parameters:
        pdb_file (str): The path to the input PDB file.
        output_gro (str): The path where the output GRO file will be saved.
        box_vector (np.ndarray): A 1D array containing the box dimensions (x, y, z, alpha, beta, gamma) in nanometers.

    Returns:
        None: This function does not return a value but writes the converted file to `output_gro`.
    """
    universe = Universe(pdb_file)
    universe.dimensions = box_vector
    universe.atoms.write(output_gro)


def zip_directories(directory_names: list, output_file_name: str) -> str:
    """
    Zips multiple directories into a single archive and moves it to a specified output path.

    This function takes a list of directory paths and an output file name, creates a ZIP archive of these directories,
    and then moves the archive to the current working directory or a specified directory. It excludes symlinks from
    the archive and ensures that the archive's internal structure reflects the relative paths of files to the
    directories being zipped. This is useful for bundling multiple directories into a single, compressed file for
    easy distribution or storage.

    Parameters:
        directory_names (list): A list of strings representing the paths to the directories to be zipped.
        output_file_name (str): The name for the output ZIP file. The '.zip' extension is appended automatically.

    Returns:
        str: The path to the final output ZIP file.

    Raises:
        ValueError: If the `directory_names` list is empty.
    """
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
    """
    Recursively updates a nested dictionary with another dictionary.

    This function takes two dictionaries, `d` and `u`. It iterates through `u`, updating `d` with `u`'s key-value pairs.
    If a value in `u` is a dictionary and the corresponding key exists in `d` with a dictionary as its value, this
    function recursively updates the nested dictionary. Otherwise, it directly updates the value. This is useful for
    merging configurations or settings where nested structures are involved without losing the nested specificity.

    Parameters:
        d (dict): The original dictionary to be updated.
        u (dict): The dictionary with updates to apply to `d`.

    Returns:
        dict: The updated dictionary `d` after applying updates from `u`.
    """
    for k, v in u.items():
        if isinstance(v, Mapping):
            d[k] = update_nested_dict(d.get(k, {}), v)
        else:
            d[k] = v
    return d


def calculate_average_com(gro_file_path: str, molecule_names: list[str]) -> np.ndarray:
    """
    Calculates the average center of mass (COM) of specified molecules in a GRO file.

    This function loads a molecular structure from a GRO file and calculates the average center of mass for the
    specified molecules. If no molecule names are provided, it calculates the COM for all atoms in the file. This
    can be useful for determining the central point of a group of molecules or the entire system for further
    analysis or manipulation.

    Parameters:
        gro_file_path (str): The path to the input GRO file containing the molecular structure.
        molecule_names (list[str]): A list of molecule names (residue names) for which to calculate the COM. If empty,
                                    the COM of all atoms in the file is calculated.

    Returns:
        np.ndarray: A 1D numpy array containing the x, y, z coordinates of the average center of mass.

    Raises:
        ValueError: If no atoms are found for the specified molecule names in the GRO file.
    """
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


def translate_gro_by_vector(gro_file_path: str, output_gro_file_path: str, com_diff_vector: np.ndarray) -> None:
    """
    Translates the positions of all atoms in a GRO file by a specified vector and saves the result to a new file.

    This function loads a molecular structure from a GRO file, applies a translation to all atom positions by adding
    the specified center of mass difference vector, and then writes the modified structure to a new GRO file. This
    can be useful for adjusting the position of molecules within a simulation box, for example, to center a molecule
    or to move it to a specific location within the box.

    Parameters:
        gro_file_path (str): The path to the input GRO file containing the molecular structure to be translated.
        output_gro_file_path (str): The path where the translated GRO file will be saved.
        com_diff_vector (np.ndarray): A 1D numpy array containing the x, y, z components of the translation vector.

    Returns:
        None: This function does not return a value but writes the translated structure to `output_gro_file_path`.
    """
    u = Universe(gro_file_path)
    u.atoms.positions += com_diff_vector
    u.atoms.write(output_gro_file_path)

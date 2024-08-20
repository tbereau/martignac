import logging
import shutil

VOLUME_PER_CG_BEAD_IN_NM3 = 0.08

__all__ = [
    "generate_solvent_with_gromacs",
    "solvate_solute_command",
    "gromacs_simulation_command",
    "run_gmx_wham",
    "VOLUME_PER_CG_BEAD_IN_NM3",
]

logger = logging.getLogger(__name__)


def generate_solvent_with_gromacs(
    gro_solvent_mol: str,
    box_length: float,
    output_name: str = "solvent",
    scale: float = 0.1,
) -> str:
    """
    Generates a command to solvate a given solvent molecule file using GROMACS.

    This function constructs a command for the GROMACS 'solvate' tool, which is used to solvate a given
    solvent molecule file within a cubic box of a specified length. The density of the solvent can be
    adjusted using the scale parameter.

    Parameters:
        gro_solvent_mol (str): The path to the GRO file of the solvent molecules.
        box_length (float): The length of each side of the cubic box in which the solvent molecules will be placed, in nm.
        output_name (str): The base name for the output GRO file. Defaults to "solvent".
        scale (float): The scaling factor for adjusting the density of the solvent. Defaults to 0.1.

    Returns:
        str: The command to execute in GROMACS to generate the solvated system.
    """
    box_size = " ".join([str(box_length) for _ in range(3)])
    return f"gmx solvate -cs {gro_solvent_mol} -box {box_size} -o {output_name}.gro -scale {scale}"


def solvate_solute_command(
    gro_solute: str,
    gro_solvent: str,
    top_solute: str,
    top_output: str = "topol.top",
    output_name: str = "solvate",
) -> str:
    """
    Generates a command to solvate a solute with a given solvent using GROMACS.

    This function constructs a command for the GROMACS 'solvate' tool, which is used to combine a solute and solvent
    in a simulation box. The dimensions of the box are extracted from the solvent .gro file. It also updates the topology
    file to include the solvent molecules.

    Parameters:
        gro_solute (str): The path to the GRO file of the solute.
        gro_solvent (str): The path to the GRO file of the solvent molecules.
        top_solute (str): The path to the topology file of the solute.
        top_output (str, optional): The path for the output topology file. Defaults to "topol.top".
        output_name (str, optional): The base name for the output GRO file. Defaults to "solvate".

    Returns:
        str: The command to execute in GROMACS to solvate the solute with the solvent.
    """

    def extract_box_size_from_gro(gro_file):
        with open(gro_file) as f:
            lines = f.readlines()
            box_dimensions = lines[-1].split()
            return " ".join([str(dim) for dim in box_dimensions])

    box_size = extract_box_size_from_gro(gro_solvent)
    shutil.copy(top_solute, top_output)
    cmd = f"gmx solvate -cp {gro_solute} -cs {gro_solvent} -box {box_size} -p {top_output} -o {output_name}.gro "
    return cmd


def run_gmx_wham(
    tpr_files: str,
    pullf_files: str,
    output_profile: str,
    output_hist: str,
    output_bstrap: str,
    output_bsprof: str,
    num_boostrap: int = 100,
    unit: str = "kT",
    z_min: float = 0.0,
    z_max: float = 5.0,
) -> str:
    """
    Generates a command to perform Weighted Histogram Analysis Method (WHAM) analysis using GROMACS.

    This function constructs a command for the GROMACS 'wham' tool, which is used to analyze the results of
    simulations, particularly for calculating free energy landscapes from umbrella sampling simulations. It
    supports specifying multiple input files for tpr and pullf files, outputting various analysis results,
    and setting the number of bootstrap analyses to improve statistical reliability.

    Parameters:
        tpr_files (str): The path to the text file containing a list of TPR files for the analysis.
        pullf_files (str): The path to the text file containing a list of pull force files for the analysis.
        output_profile (str): The path for the output free energy profile.
        output_hist (str): The path for the output histogram file.
        output_bstrap (str): The path for the output bootstrap results file.
        output_bsprof (str): The path for the output bootstrap profile file.
        num_boostrap (int): The number of bootstrap analyses to perform. Defaults to 100.
        unit (str): The unit for the free energy calculation. Defaults to "kT".

    Returns:
        str: The command to execute in GROMACS for WHAM analysis.
    """
    cmd = (
        f"gmx wham -it {tpr_files} -if {pullf_files} -o {output_profile} "
        f"-hist {output_hist} -bsres {output_bstrap} -bsprof {output_bsprof} "
        f"-min {z_min} -max {z_max} -zprof0 {z_max} "
        f"-nBootstrap {num_boostrap} -unit {unit}"
    )
    return cmd


def gromacs_simulation_command(
    mdp: str,
    top: str,
    gro: str,
    name: str,
    n_max_warn: int = 10,
    n_threads: int = 1,
    verbose: bool = True,
) -> str:
    """
    Prepares and executes a molecular dynamics simulation using GROMACS.

    This function generates the necessary commands to prepare (grompp) and run (mdrun) a simulation in GROMACS. It allows
    for specifying the input files, the number of threads, and whether to run in verbose mode. The function combines the
    preparation and execution steps into a single command, facilitating streamlined execution of simulations.

    Parameters:
        mdp (str): The path to the MDP file containing the simulation parameters.
        top (str): The path to the topology file for the system.
        gro (str): The path to the GRO file containing the system's initial structure.
        name (str): The base name for the output files.
        n_max_warn (int, optional): The maximum number of allowed warnings during the preparation step. Defaults to 10.
        n_threads (int, optional): The number of threads to use for the simulation. Defaults to 1.
        verbose (bool, optional): If True, runs the simulation in verbose mode. Defaults to True.

    Returns:
        str: The command to execute in GROMACS that combines the preparation and execution steps.
    """
    grompp_cmd = f"gmx grompp -f {mdp} -p {top} -c {gro} -o {name}.tpr -po {name}_out.mdp -maxwarn {n_max_warn}"
    mdrun_cmd = f"gmx mdrun -nt {n_threads} -deffnm {name}"
    if verbose:
        mdrun_cmd += " -v"
    return f"{grompp_cmd} && {mdrun_cmd}"


def generate_lambdas(num_sim: int, turn_off_coulomb: bool = False):
    """
    Generates van der Waals and Coulomb lambda schedules for Gromacs alchemical transformations.

    Parameters:
    - num_sim (int): Number of state points.
    - turn_off_coulomb (bool): If True, Coulomb interactions are turned off and the schedule is dedicated to vdw.

    Returns:
    - Tuple[List[float], List[float]]: Lambda schedules for van der Waals and Coulomb interactions.
    """
    if num_sim < 2:
        raise ValueError("Number of state points must be at least 2.")

    # Adjust the calculation for vdw_points to ensure a linear increase
    vdw_points = num_sim if turn_off_coulomb else (num_sim + 1) // 2

    # Generate linearly increasing lambda values for vdw
    vdw_lambdas = [i / (vdw_points - 1) for i in range(vdw_points)]

    if turn_off_coulomb:
        coul_lambdas = [0.0] * num_sim
    else:
        # Adjust the lambda values for Coulomb, starting from 0 after vdw is fully coupled
        coul_lambdas = [0.0] * (vdw_points - 1) + [
            (i - (vdw_points - 1)) / (num_sim - vdw_points)
            for i in range(vdw_points - 1, num_sim)
        ]

    # Ensure lambdas list length is num_sim for vdw_lambdas when turn_off_coulomb is False
    if not turn_off_coulomb:
        vdw_lambdas += [1.0] * (num_sim - vdw_points)

    return vdw_lambdas, coul_lambdas

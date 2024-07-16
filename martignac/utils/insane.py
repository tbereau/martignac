import logging

from martignac.liquid_models.mixtures import LiquidMixture

logger = logging.getLogger(__name__)


def generate_bilayer_with_insane(
    lipids: LiquidMixture,
    solvent: LiquidMixture,
    box_length_xy: float,
    box_length_z: float,
    gro_bilayer_gen: str,
    top_bilayer: str,
) -> str:
    """
    Generates a command for creating a lipid bilayer system using the INSANE tool.

    This function constructs a command line instruction for the INSANE (INitiation of Solvated membrANEs) tool,
    which is used to generate a lipid bilayer system. The command includes specifications for the lipid and solvent
    compositions, the dimensions of the simulation box, and the output files for the generated system.

    Parameters:
        lipids (LiquidMixture): A LiquidMixture object representing the lipid composition of the bilayer.
        solvent (LiquidMixture): A LiquidMixture object representing the solvent composition.
        box_length_xy (float): The length of the simulation box in the x and y dimensions, in nanometers.
        box_length_z (float): The height of the simulation box in the z dimension, in nanometers.
        gro_bilayer_gen (str): The filename for the generated GRO file of the bilayer system.
        top_bilayer (str): The filename for the generated topology (TOP) file of the bilayer system.

    Returns:
        str: The command to be executed for generating the bilayer system with INSANE.
    """
    return (
        f"insane -l {lipids.to_insane_format()} -d {box_length_xy} -dz {box_length_z} "
        f"-sol {solvent.to_insane_format()} -o {gro_bilayer_gen} -p {top_bilayer}"
    )

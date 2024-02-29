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
    return (
        f"insane -l {lipids.to_insane_format()} -d {box_length_xy} -dz {box_length_z} "
        f"-sol {solvent.to_insane_format()} -o {gro_bilayer_gen} -p {top_bilayer}"
    )

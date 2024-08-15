from dataclasses import dataclass
from typing import Optional

import numpy as np


@dataclass
class LiquidComponent:
    """
    Represents a single component of a liquid mixture.

    This class models a component within a liquid mixture, including its name, the fraction of the total mixture it
    represents, and an optional integer component value. It provides methods for initialization validation, creating
    instances from dictionaries, and formatting component information for specific external systems.

    Attributes:
        name (str): The name of the liquid component.
        fraction (float): The fraction of the total mixture this component represents, must be between 0 and 1.
        integer_component (Optional[int]): An optional integer value associated with the component, default is None.

    Methods:
        __post_init__(self): Validates the fraction attribute to ensure it is within the expected range (0 to 1).
        from_dict(cls, input_dict: dict) -> "LiquidComponent": Class method to create an instance from a dictionary.
        to_insane_format(self) -> str: Formats the component information for use in an external system, returning a string.
    """

    name: str
    fraction: float
    integer_component: Optional[int] = None

    def __post_init__(self):
        if not 0 <= self.fraction <= 1:
            raise ValueError(
                f"unexpected component fraction {self.fraction} for {self.name}"
            )

    @classmethod
    def from_dict(cls, input_dict: dict) -> "LiquidComponent":
        if "name" not in input_dict:
            raise KeyError(f"Missing key 'name' in {input_dict}")
        if "fraction" not in input_dict:
            raise KeyError(f"Missing key 'fraction' in {input_dict}")
        return LiquidComponent(
            name=input_dict.get("name"), fraction=input_dict.get("fraction")
        )

    def to_insane_format(self) -> str:
        return f"{self.name}:{self.integer_component}"


@dataclass
class LiquidMixture:
    """
    Represents a mixture of liquid components.

    This class models a mixture of various liquid components, allowing for operations such as normalization checks,
    setting integer contributions based on a scaling factor, and formatting the mixture information for specific
    external systems. It supports initialization from a list of dictionaries representing individual components,
    ensuring that the total fraction of components equals 1.0 for a normalized mixture.

    Attributes:
        components (list[LiquidComponent]): A list of `LiquidComponent` instances representing the components of the mixture.
        scaling_factor (float): A scaling factor used to calculate integer contributions of each component, default is 10.0.

    Methods:
        __post_init__(self): Validates the mixture to ensure it is normalized (i.e., the sum of component fractions equals 1.0).
        from_list_of_dicts(cls, list_of_dicts: list[dict]) -> "LiquidMixture": Class method to create an instance from a list of dictionaries.
        set_integer_contributions(self): Sets the integer contribution for each component based on the scaling factor.
        to_insane_format(self) -> str: Formats the mixture information for use in an external system, returning a string.
        solvent_names (property): Returns a list of names of the components in the mixture.
    """

    components: list[LiquidComponent]
    scaling_factor: float = 10.0

    def __post_init__(self):
        if not np.isclose(sum([c.fraction for c in self.components]), 1.0):
            raise ValueError(f"liquid mixture is not normalized: {self.components}")

    @classmethod
    def from_list_of_dicts(cls, list_of_dicts: list[dict]) -> "LiquidMixture":
        components = [LiquidComponent.from_dict(d) for d in list_of_dicts]
        return LiquidMixture(components)

    def set_integer_contributions(self) -> None:
        for c in self.components:
            c.integer_component = int(c.fraction * self.scaling_factor)

    def to_insane_format(self, separator: str = " ") -> str:
        self.set_integer_contributions()
        return separator.join([c.to_insane_format() for c in self.components])

    @property
    def solvent_names(self) -> list[str]:
        return [c.name for c in self.components]

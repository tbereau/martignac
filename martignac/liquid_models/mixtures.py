from dataclasses import dataclass
from typing import Optional

import numpy as np


@dataclass
class LiquidComponent:
    name: str
    fraction: float
    integer_component: Optional[int] = None

    def __post_init__(self):
        if not 0 <= self.fraction <= 1:
            raise ValueError(f"unexpected component fraction {self.fraction} for {self.name}")

    @classmethod
    def from_dict(cls, input_dict: dict) -> "LiquidComponent":
        if "name" not in input_dict.keys():
            raise KeyError(f"Missing key 'name' in {input_dict}")
        if "fraction" not in input_dict.keys():
            raise KeyError(f"Missing key 'fraction' in {input_dict}")
        return LiquidComponent(name=input_dict.get("name"), fraction=input_dict.get("fraction"))

    def to_insane_format(self) -> str:
        return f"{self.name}:{self.integer_component}"


@dataclass
class LiquidMixture:
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

    def to_insane_format(self) -> str:
        self.set_integer_contributions()
        return " ".join([c.to_insane_format() for c in self.components])

    @property
    def solvent_names(self) -> list[str]:
        return [c.name for c in self.components]

"""Top-level package for Martignac."""

__author__ = """Tristan Bereau"""
__email__ = "bereau@uni-heidelberg.de"

import confuse
from confuse.core import Configuration


__all__ = ["config"]

CONFIG = confuse.Configuration('martignac', __name__)


def config() -> Configuration:
    return CONFIG


"""Top-level package for Martignac."""

__author__ = """Tristan Bereau"""
__email__ = "bereau@uni-heidelberg.de"

import confuse
from confuse.core import Configuration

from martignac.logging import logger

__all__ = ["config"]

CONFIG = confuse.Configuration("martignac", __name__)
LOGGING_LEVEL = "INFO"
logger.setLevel(LOGGING_LEVEL)


def config() -> Configuration:
    return CONFIG

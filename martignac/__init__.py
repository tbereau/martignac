"""Top-level package for Martignac."""

__author__ = """Tristan Bereau"""
__email__ = "bereau@uni-heidelberg.de"

import os
from pathlib import Path
import confuse
from confuse.core import Configuration

from martignac.logging import logger

__all__ = ["config"]

if os.environ.get('MARTIGNACDIR') is None:
    os.environ['MARTIGNACDIR'] = str(Path.cwd())
CONFIG = confuse.Configuration("martignac", __name__)
LOGGING_LEVEL = "INFO"
logger.setLevel(LOGGING_LEVEL)

if Path(CONFIG.config_dir()).resolve() == Path(__file__).parent.resolve():
    logger.warning('You are currently using the default configuration of Martignac.')

def config() -> Configuration:
    return CONFIG

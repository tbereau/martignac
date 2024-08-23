from os import makedirs
from os.path import isdir

import signac

from martignac.workflows.solute_generation import SoluteGenFlow

if not isdir(SoluteGenFlow.workspace_path):
    makedirs(SoluteGenFlow.workspace_path)

project = signac.init_project(path=SoluteGenFlow.workspace_path)

for solute_name in [
    "P6",
    "P5",
    "P4",
    "C4",
    "C3",
    "C2",
    "C1",
    "C5",
]:  # , "P6 C3, 0-1", "P6 C3 N1, 0_1 1_2 0-2", "P1"]:
    sp = {"type": "solute", "solute_name": solute_name}
    job = project.open_job(sp).init()

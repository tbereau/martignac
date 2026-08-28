import importlib
import importlib.resources
import importlib.util
import sys
import types


def _install_pkg_resources_shim() -> None:
    """
    Supplies the sliver of ``pkg_resources`` that gravis 0.1.0 relies on.

    gravis reads its bundled HTML and JavaScript templates through
    ``pkg_resources.resource_string``. setuptools removed ``pkg_resources`` in
    version 82.0.0, and gravis has had no release since 0.1.0, so on a recent
    environment ``import gravis`` fails outright. Rather than pinning setuptools
    below the removal, back the single function gravis calls with
    ``importlib.resources`` and register it under the expected module name.

    The shim is installed only when ``pkg_resources`` is genuinely absent, and
    is deliberately minimal: it covers gravis's usage, not the wider API.
    """
    if "pkg_resources" in sys.modules or importlib.util.find_spec("pkg_resources"):
        return

    def resource_string(package: str, resource_path: str) -> bytes:
        anchor = sys.modules.get(package) or importlib.import_module(package)
        if getattr(anchor, "__path__", None) is None:
            # gravis anchors on a module rather than a package; its resources
            # sit alongside that module, in the parent package's directory.
            package = anchor.__package__ or package.rpartition(".")[0]
        return importlib.resources.files(package).joinpath(resource_path).read_bytes()

    shim = types.ModuleType("pkg_resources")
    shim.resource_string = resource_string
    sys.modules["pkg_resources"] = shim


_install_pkg_resources_shim()

import gravis as gv
import networkx as nx
import numpy as np

from martignac.utils.martini_flow_projects import MartiniFlowProject


def generate_gravis_network(project: MartiniFlowProject):
    """
    Generates a directed graph visualization for a given MartiniFlowProject using Gravis.

    This function takes a MartiniFlowProject instance, extracts its operations and the operation graph,
    and generates a directed graph visualization. The visualization is created using the Gravis library,
    which is built on top of D3.js for interactive and dynamic graph visualizations in Jupyter notebooks
    and other web-based environments.

    Parameters:
        project (MartiniFlowProject): An instance of MartiniFlowProject containing operations and their relationships.

    Returns:
        gv.d3: A Gravis D3 graph object that can be displayed in Jupyter notebooks or web environments.
    """
    ops = project.operations.keys()
    adj = np.asarray(project.detect_operation_graph())

    g = nx.DiGraph(adj)
    for i, op in enumerate(ops):
        g.nodes[i]["label"] = op
    return gv.d3(g, node_label_data_source="label", zoom_factor=1.5)

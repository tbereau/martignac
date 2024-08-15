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
    ops = project._operations.keys()
    adj = np.asarray(project.detect_operation_graph())

    g = nx.DiGraph(adj)
    for i, op in enumerate(ops):
        g.nodes[i]["label"] = op
    return gv.d3(g, node_label_data_source="label", zoom_factor=1.5)

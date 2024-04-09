import gravis as gv
import numpy as np
import networkx as nx

from martignac.utils.martini_flow_projects import MartiniFlowProject

def generate_gravis_network(project: MartiniFlowProject):
    ops = project.operations.keys()
    adj = np.asarray(project.detect_operation_graph())

    g = nx.DiGraph(adj)
    for i, op in enumerate(ops):
        g.nodes[i]["label"] = op
    return gv.d3(g, node_label_data_source="label", zoom_factor=1.5)

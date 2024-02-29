from dataclasses import dataclass
from functools import cached_property
from typing import Dict

import networkx as nx
import yaml
from signac.job import Job

from martignac.utils.martini_flow_projects import MartiniFlowProject


@dataclass
class NomadWorkflow:
    project: MartiniFlowProject
    job: Job

    @property
    def project_name(self) -> str:
        return self.project.__class__.__name__

    @property
    def gromacs_logs(self) -> dict:
        return self.job.document["gromacs_logs"].get(self.project_name, default={})

    @cached_property
    def graph(self) -> nx.DiGraph:
        graph = nx.DiGraph()
        operations = self.project._collect_operations()
        if "gromacs_logs" not in self.job.document:
            return nx.DiGraph()
        for op_name, _op_func in operations:
            if op_name in self.gromacs_logs:
                graph.add_node(op_name, label=self.gromacs_logs[op_name])
        for op_name, _op_func in operations:
            for f_i, f_j in self._forward_link_dict.items():
                if self._co_func_to_name[f_i] == op_name and (
                    self._co_func_to_name[f_i] in self.gromacs_logs and self._co_func_to_name[f_j] in self.gromacs_logs
                ):
                    graph.add_edge(self._co_func_to_name[f_i], self._co_func_to_name[f_j])
        return graph

    def build_workflow_yaml(self, destination_filename: str) -> None:
        graph = self.graph
        workflow_output = {"workflow2": {}}
        sub_dict = workflow_output["workflow2"]
        if len(graph.nodes) == 0:
            return None
        node_name, node_dict = list(graph.nodes.items())[0]
        sub_dict["inputs"] = [_construct_yaml_section("input system", node_dict["label"], 0)]
        node_name, node_dict = list(graph.nodes.items())[-1]
        sub_dict["outputs"] = [_construct_yaml_section("final output", node_dict["label"], -1)]
        sub_dict["tasks"] = []
        for i, (node_name, node_dict) in enumerate(list(graph.nodes.items())):
            node_label = node_dict["label"]
            input_node_label = node_label if i == 0 else list(graph.nodes.items())[i - 1][1]["label"]
            snapshot_number = 0 if i == 0 else -1
            sub_dict["tasks"].append(
                {
                    "m_def": "nomad.datamodel.metainfo.workflow.TaskReference",
                    "task": f"../upload/archive/mainfile/{node_label}#/workflow2",
                    "name": node_name,
                    "inputs": [_construct_yaml_section("input", input_node_label, snapshot_number)],
                    "outputs": [_construct_yaml_section("output", node_label, -1)],
                }
            )
        self.register_workflow_yaml(destination_filename)
        with open(destination_filename, "w") as f:
            yaml.dump(workflow_output, f)

    def register_workflow_yaml(self, filename: str) -> None:
        class_name = self.project.__class__.__name__
        if "nomad_workflow" not in self.job.document:
            self.job.document["nomad_workflow"] = {}
        self.job.document["nomad_workflow"][class_name] = filename

    @property
    def _forward_link_dict(self) -> dict:
        forward_link_dict = {}
        for f, pre_f in list(self.project._collect_preconditions().items()):
            for f_2, post_f2 in list(self.project._collect_postconditions().items()):
                if pre_f == post_f2:
                    forward_link_dict[f_2] = f
        return forward_link_dict

    @property
    def _co_func_to_name(self) -> dict:
        return {co_func: co_name for co_name, co_func in self.project._collect_operations()}


@dataclass
class NomadTopLevelWorkflow:
    projects: Dict[str, MartiniFlowProject]
    job: Job
    graph: nx.DiGraph

    def build_workflow_yaml(self, destination_filename: str) -> None:
        workflow_output = {"workflow2": {}}
        sub_dict = workflow_output["workflow2"]
        if len(self.graph.nodes) == 0:
            return None
        start_nodes = [n for n, d in self.graph.in_degree if d == 0]
        sub_dict["inputs"] = []
        for node in start_nodes:
            project = self.projects[node]
            input_graph = NomadWorkflow(project, self.job).graph
            if len(input_graph) == 0:
                return None
            _, node_dict = list(input_graph.nodes.items())[0]
            sub_dict["inputs"].append(_construct_yaml_section(f"Input for {node}", node_dict["label"], 0))

        end_nodes = [n for n, d in self.graph.out_degree if d == 0]
        sub_dict["outputs"] = []
        for node in end_nodes:
            project = self.projects[node]
            input_graph = NomadWorkflow(project, self.job).graph
            if len(input_graph) == 0:
                return None
            _, node_dict = list(input_graph.nodes.items())[-1]
            sub_dict["outputs"].append(_construct_yaml_section(f"Output for {node}", node_dict["label"], -1))

        sub_dict["tasks"] = []
        for node_first in self.graph.nodes:
            other_nodes = (
                [node_first]
                if len(self.graph.edges(node_first)) == 0
                else [n2 for _, n2 in self.graph.edges(node_first)]
            )
            for node_next in other_nodes:
                workflow_name = self.job.document["nomad_workflow"][node_first]
                input_graph = NomadWorkflow(self.projects[node_first], self.job).graph
                input_file = list(input_graph.nodes.values())[0]["label"]
                output_graph = NomadWorkflow(self.projects[node_next], self.job).graph
                snapshot_number = -1 if node_first == node_next else 0
                output_file = list(output_graph.nodes.values())[snapshot_number]["label"]
                snapshot_number = 0 if self.graph.in_degree[node_first] == 0 else -1
                sub_dict["tasks"].append(
                    {
                        "m_def": "nomad.datamodel.metainfo.workflow.TaskReference",
                        "task": f"../upload/archive/mainfile/{workflow_name}#/workflow2",
                        "name": node_first,
                        "inputs": [_construct_yaml_section("input", input_file, snapshot_number)],
                        "outputs": [_construct_yaml_section("output", output_file, -1)],
                    }
                )
        with open(destination_filename, "w") as f:
            yaml.dump(workflow_output, f)


def _construct_yaml_section(name: str, log_file: str, number: int) -> dict:
    return {
        "name": name,
        "section": f"../upload/archive/mainfile/{log_file}#/run/0/calculation/{number}",
    }

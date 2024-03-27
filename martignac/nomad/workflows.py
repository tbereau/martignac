import logging
from copy import copy
from dataclasses import dataclass, field
from functools import cached_property
from typing import Dict, Literal, Optional, cast

import networkx as nx
import numpy as np
import yaml
from signac.job import Job

from martignac.utils.martini_flow_projects import MartiniFlowProject
from martignac.utils.misc import update_nested_dict

logger = logging.getLogger(__name__)
DEFAULT_M_DEF = "nomad.datamodel.metainfo.workflow.TaskReference"


SectionType = Literal["task", "workflow", "default"]


@dataclass
class NomadSection:
    name: str
    section_type: SectionType
    label: str
    snapshot_number: int = 0
    run_number: int = 0
    upload_id: Optional[str] = None

    @property
    def is_task(self) -> bool:
        return self.section_type == "task"

    @property
    def is_workflow(self) -> bool:
        return self.section_type == "workflow"

    @property
    def log_file(self) -> Optional[str]:
        return self.label if not self.is_task else None

    @property
    def section(self) -> str:
        return f"{self.upload_prefix}#/run/{self.run_number}/calculation/{self.snapshot_number}"

    @property
    def upload_prefix(self) -> str:
        upload_prefix = f"/uploads/{self.upload_id}" if self.upload_id else "../upload"
        return f"{upload_prefix}/archive/mainfile/{self.log_file}" if self.log_file else ""

    def to_dict(self) -> dict:
        return {"name": self.name, "section": self.section}

    def add_job_id(self, job_: Job) -> "NomadSection":
        if job_.id not in self.label and self.label != "run" and not self.upload_id:
            self.label = f"{job_.id}/{self.label}"
        return self


@dataclass
class NomadTask:
    name: str
    m_def: str = DEFAULT_M_DEF
    inputs: list[NomadSection] = field(default_factory=list)
    outputs: list[NomadSection] = field(default_factory=list)
    task_section: Optional[NomadSection] = None

    def __post_init__(self):
        for i in range(len(self.inputs)):
            self.inputs[i].name = "input"
        for i in range(len(self.outputs)):
            self.outputs[i].name = "output"

    @property
    def task(self) -> Optional[str]:
        if self.task_section.is_task:
            return None
        return self.task_section.upload_prefix + "#/workflow2"

    def to_dict(self) -> dict:
        output_dict = {
            "name": self.name,
            "m_def": self.m_def,
            "inputs": [i.to_dict() for i in self.inputs],
            "outputs": [o.to_dict() for o in self.outputs],
        }
        if self.task:
            output_dict["task"] = self.task
        return output_dict


@dataclass
class NomadWorkflowArchive:
    name: str = "workflow2"
    inputs: list[NomadSection] = field(default_factory=list)
    outputs: list[NomadSection] = field(default_factory=list)
    tasks: list[NomadTask] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {
            self.name: {
                "inputs": [i.to_dict() for i in self.inputs],
                "outputs": [o.to_dict() for o in self.outputs],
                "tasks": [t.to_dict() for t in self.tasks],
            },
        }

    def to_yaml(self, destination_filename: str) -> None:
        with open(destination_filename, "w") as f:
            yaml.dump(self.to_dict(), f)

    @classmethod
    def from_multiple_jobs(
        cls, project: MartiniFlowProject, jobs: list[Job], aggregate_same_task_names: bool = True
    ) -> "NomadWorkflowArchive":
        def filter_unique(ele):
            final_inputs = []
            for inp in ele:
                if inp not in final_inputs:
                    final_inputs.append(inp)
            return final_inputs

        archive = NomadWorkflowArchive()
        for job in jobs:
            workflow = NomadWorkflow(project, job, is_top_level=True)
            job_inputs = [inp.add_job_id(job) for inp in copy(workflow.generate_archive().inputs)]
            archive.inputs.extend(job_inputs)
            job_outputs = [out.add_job_id(job) for out in copy(workflow.generate_archive().outputs)]
            archive.outputs.extend(job_outputs)
            job_tasks = copy(workflow.generate_archive().tasks)

            for task in job_tasks:
                task.inputs = [inp.add_job_id(job) for inp in task.inputs]
                task.outputs = [out.add_job_id(job) for out in task.outputs]
                task.task_section = task.task_section.add_job_id(job)
            archive.tasks.extend(job_tasks)

        if aggregate_same_task_names:
            final_tasks = []
            for task in archive.tasks:
                if task.name not in [t.name for t in final_tasks]:
                    final_tasks.append(task)
                else:
                    dest_task = next(t for t in final_tasks if t.name == task.name)
                    dest_task.inputs = filter_unique(dest_task.inputs)
                    dest_task.outputs = filter_unique(dest_task.outputs)
            archive.tasks = final_tasks

        archive.inputs = filter_unique(archive.inputs)
        archive.outputs = filter_unique(archive.outputs)
        archive.tasks = filter_unique(archive.tasks)
        return archive


@dataclass
class NomadWorkflow:
    project: MartiniFlowProject
    job: Job
    is_top_level: bool
    add_job_id: bool = False

    def __post_init__(self):
        self.task_elements: Dict[str, NomadSection] = {}
        self.task_counter: int = 0

    @property
    def project_name(self) -> str:
        return self.project.__class__.__name__

    @property
    def gromacs_logs(self) -> dict:
        return self.job.doc[self.project_name].get("gromacs_logs", {})

    @property
    def tasks(self) -> dict:
        return self.job.doc[self.project_name].get("tasks", {})

    @property
    def workflows(self) -> dict:
        return self.job.doc[self.project_name].get("workflows", {}) if self.is_top_level else {}

    @property
    def all_tasks(self) -> dict:
        return dict(self.gromacs_logs) | dict(self.tasks) | dict(self.workflows)

    def register_section(self, operation_name: str) -> None:
        section_type = self._section_type(operation_name)
        label = self.all_tasks[operation_name]
        if section_type == "workflow":
            label = self.job.doc[self.all_tasks[operation_name]].get("nomad_workflow", self.all_tasks[operation_name])
        upload_id = (
            self.job.doc[self.all_tasks[operation_name]].get("nomad_upload_id", None)
            if section_type == "workflow"
            else None
        )
        section = NomadSection(
            name=operation_name,
            section_type=section_type,
            label=label,
            run_number=self.task_counter if section_type == "task" else 0,
            upload_id=upload_id,
        )
        if section.is_task:
            self.task_counter += 1
        if self.add_job_id:
            section.add_job_id(self.job)
        self.task_elements[operation_name] = section

    @cached_property
    def graph(self) -> nx.DiGraph:
        operations = list(self.project.operations.keys())
        adjacency_matrix = np.asarray(self.project.detect_operation_graph())
        signac_graph = nx.DiGraph(adjacency_matrix)
        graph = nx.DiGraph()
        all_tasks = dict(self.gromacs_logs) | dict(self.tasks) | dict(self.workflows)
        for node_index in signac_graph.nodes:
            op_name = operations[node_index]
            if op_name in all_tasks:
                graph.add_node(op_name, label=all_tasks[op_name], is_task=op_name in self.tasks)
                self.register_section(op_name)
        for node_1, node_2 in signac_graph.edges:
            if (op_name_1 := operations[node_1]) in graph.nodes and (op_name_2 := operations[node_2]) in graph.nodes:
                graph.add_edge(op_name_1, op_name_2)
        return graph

    def build_workflow_yaml(self, destination_filename: str) -> None:
        archive = self.generate_archive()
        archive.to_yaml(destination_filename)
        project_name = self.project.class_name()
        self.job.doc = update_nested_dict(self.job.doc, {project_name: {"nomad_workflow": destination_filename}})

    def generate_archive(self) -> NomadWorkflowArchive:
        archive = NomadWorkflowArchive()
        archive.inputs = []
        for node in [n for n, d in self.graph.in_degree if d == 0]:
            element = self.task_elements[node]
            archive.inputs.append(element)

        for node in [n for n, d in self.graph.out_degree if d == 0]:
            element = self.task_elements[node]
            element.snapshot_number = -1
            archive.outputs.append(element)
        node_names = list(self.graph.nodes)
        for _i, node_name in enumerate(node_names):
            in_nodes = [n for n, _ in self.graph.in_edges([node_name])]
            out_nodes = [n for _, n in self.graph.out_edges([node_name])]
            if not in_nodes:
                input_node = copy(self.task_elements[node_name])
                task_inputs = [input_node]
            else:
                if all([self.task_elements[n].is_task for n in in_nodes]):
                    task_inputs = [copy(self.task_elements[node_name])]
                else:
                    task_inputs = [copy(self.task_elements[n]) for n in in_nodes]
                    for _input in task_inputs:
                        _input.snapshot_number = -1
            if not out_nodes:
                output_node = copy(self.task_elements[node_name])
                task_outputs = [output_node]
            else:
                if self.task_elements[node_name].is_task:
                    task_outputs = [copy(self.task_elements[n]) for n in out_nodes]
                else:
                    task_outputs = [copy(self.task_elements[node_name])]
                    for _output in task_outputs:
                        _output.snapshot_number = -1
            archive.tasks.append(
                NomadTask(
                    node_name,
                    inputs=task_inputs,
                    outputs=task_outputs,
                    task_section=self.task_elements[node_name],
                )
            )
        return archive

    def _section_type(self, operation_name: str) -> SectionType:
        if operation_name in self.tasks:
            return "task"
        elif operation_name in self.workflows:
            return "workflow"
        return "default"


@MartiniFlowProject.label
def nomad_workflow_is_built(job: Job) -> bool:
    project = cast(MartiniFlowProject, job.project)
    return job.isfile(project.nomad_workflow)


def build_nomad_workflow(job, is_top_level: bool = False, add_job_id: bool = False):
    project = cast(MartiniFlowProject, job.project)
    workflow = NomadWorkflow(project, job, is_top_level, add_job_id=add_job_id)
    destination_filename = project.nomad_top_level_workflow if is_top_level else project.nomad_workflow
    workflow.build_workflow_yaml(destination_filename)


def build_nomad_workflow_with_multiple_jobs(project: MartiniFlowProject, jobs: list[Job]):
    archive = NomadWorkflowArchive.from_multiple_jobs(project, jobs)
    destination_filename = project.nomad_top_level_workflow
    archive.to_yaml(destination_filename)

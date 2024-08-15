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
    """
    Represents a section within a NOMAD workflow or task.

    This class is used to define and manage the properties of a section, which can be a task, a workflow, or a default
    type within the NOMAD system. It includes attributes to specify the section's name, type, label, snapshot number,
    run number, and an optional upload ID. Additional properties provide convenience methods to determine if the
    section is a task or workflow, to generate a log file path, to construct the section's URL, and to convert the
    section information into a dictionary format.

    Attributes:
        name (str): The name of the section.
        section_type (SectionType): The type of the section, which can be 'task', 'workflow', or 'default'.
        label (str): A label for the section, used in generating paths and URLs.
        snapshot_number (int): The snapshot number associated with this section, defaulting to 0.
        run_number (int): The run number associated with this section, defaulting to 0.
        upload_id (Optional[str]): An optional upload ID for the section, defaulting to None.

    Properties:
        is_task (bool): Returns True if the section is a task, False otherwise.
        is_workflow (bool): Returns True if the section is a workflow, False otherwise.
        log_file (Optional[str]): Returns the log file path if the section is not a task, None otherwise.
        section (str): Constructs and returns the section's URL based on its properties.
        upload_prefix (str): Constructs and returns the upload prefix URL for the section.

    Methods:
        to_dict() -> dict: Converts the section's attributes into a dictionary.
        add_job_id(job_: Job) -> 'NomadSection': Adds a job ID to the section's label if applicable.
    """

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
        return (
            f"{upload_prefix}/archive/mainfile/{self.log_file}" if self.log_file else ""
        )

    def to_dict(self) -> dict:
        return {"name": self.name, "section": self.section}

    def add_job_id(self, job_: Job) -> "NomadSection":
        if job_.id not in self.label and self.label != "run" and not self.upload_id:
            self.label = f"{job_.id}/{self.label}"
        return self


@dataclass
class NomadTask:
    """
    Represents a task within a NOMAD workflow.

    This class is designed to encapsulate the details of a task in the NOMAD system, including its name, definition,
    inputs, outputs, and an optional task-specific section. It provides a structured way to manage tasks within
    workflows, facilitating the creation, manipulation, and serialization of task-related information.

    Attributes:
        name (str): The name of the task, serving as a unique identifier within the workflow.
        m_def (str): The metainfo definition for the task, defaulting to a predefined value.
        inputs (list[NomadSection]): A list of `NomadSection` instances representing the inputs to the task.
        outputs (list[NomadSection]): A list of `NomadSection` instances representing the outputs from the task.
        task_section (Optional[NomadSection]): An optional `NomadSection` instance representing a task-specific section.

    Methods:
        __post_init__(self): A post-initialization method to set default names for inputs and outputs.
        task(self) -> Optional[str]: A property that returns the task's URL if it is not a task-specific section.
        to_dict(self) -> dict: Serializes the task's attributes into a dictionary for easy export and manipulation.
    """

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
    """
    Represents the archival structure of a NOMAD workflow.

    This class encapsulates the structure necessary for representing a complete NOMAD workflow, including its inputs,
    outputs, and the tasks that comprise the workflow. It provides methods for serializing the workflow to a dictionary
    and saving it as a YAML file, which can be used for workflow reconstruction or analysis. Additionally, it offers a
    class method for creating a workflow archive from multiple jobs, allowing for the aggregation of tasks and data
    across different computational jobs within the same project.

    Attributes:
        name (str): The name of the workflow, used as a key in the serialized output.
        inputs (list[NomadSection]): A list of sections representing the inputs to the workflow.
        outputs (list[NomadSection]): A list of sections representing the outputs from the workflow.
        tasks (list[NomadTask]): A list of tasks that are part of the workflow.

    Methods:
        to_dict(self) -> dict: Serializes the workflow archive to a dictionary.
        to_yaml(self, destination_filename: str) -> None: Saves the serialized workflow archive as a YAML file to the specified location.
        from_multiple_jobs(cls, project: MartiniFlowProject, jobs: list[Job], aggregate_same_task_names: bool = True) -> 'NomadWorkflowArchive':
            Class method to create a workflow archive from multiple jobs, optionally aggregating tasks with the same name.
    """

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
        cls,
        project: MartiniFlowProject,
        jobs: list[Job],
        aggregate_same_task_names: bool = True,
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
            job_inputs = [
                inp.add_job_id(job) for inp in copy(workflow.generate_archive().inputs)
            ]
            archive.inputs.extend(job_inputs)
            job_outputs = [
                out.add_job_id(job) for out in copy(workflow.generate_archive().outputs)
            ]
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
    """
    Manages the construction and serialization of a NOMAD workflow.

    This class is responsible for creating a representation of a workflow within the NOMAD system, which includes
    tasks, inputs, and outputs as defined by the user's project and job configurations. It utilizes the project's
    operations and the job's documentation to construct a directed graph representing the workflow, which can then
    be serialized to a YAML file for use within the NOMAD system.

    Attributes:
        project (MartiniFlowProject): The project instance containing operations and configurations for the workflow.
        job (Job): The job instance containing specific execution details and documentation for the workflow.
        is_top_level (bool): Flag indicating if the workflow is at the top level of the project hierarchy.
        add_job_id (bool): Flag indicating if the job ID should be added to section labels for uniqueness.

    Properties:
        project_name (str): Returns the name of the project class.
        gromacs_logs (dict): Retrieves GROMACS log entries from the job's documentation.
        tasks (dict): Retrieves task entries from the job's documentation.
        workflows (dict): Retrieves workflow entries from the job's documentation if `is_top_level` is True.
        all_tasks (dict): Aggregates all tasks, workflows, and GROMACS logs into a single dictionary.
        graph (nx.DiGraph): Constructs and returns a directed graph representing the workflow structure.

    Methods:
        register_section(operation_name: str): Registers a section (task, workflow, or default) based on the operation name.
        build_workflow_yaml(destination_filename: str): Serializes the workflow to a YAML file.
        generate_archive() -> NomadWorkflowArchive: Generates a `NomadWorkflowArchive` instance representing the workflow.
        _section_type(operation_name: str) -> SectionType: Determines the section type (task, workflow, default) for a given operation name.
    """

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
        return (
            self.job.doc[self.project_name].get("workflows", {})
            if self.is_top_level
            else {}
        )

    @property
    def all_tasks(self) -> dict:
        return dict(self.gromacs_logs) | dict(self.tasks) | dict(self.workflows)

    def register_section(self, operation_name: str) -> None:
        section_type = self._section_type(operation_name)
        label = self.all_tasks[operation_name]
        if section_type == "workflow":
            label = self.job.doc[self.all_tasks[operation_name]].get(
                "nomad_workflow", self.all_tasks[operation_name]
            )
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
        operations = list(self.project._operations.keys())
        adjacency_matrix = np.asarray(self.project.detect_operation_graph())
        signac_graph = nx.DiGraph(adjacency_matrix)
        graph = nx.DiGraph()
        all_tasks = dict(self.gromacs_logs) | dict(self.tasks) | dict(self.workflows)
        for node_index in signac_graph.nodes:
            op_name = operations[node_index]
            if op_name in all_tasks:
                graph.add_node(
                    op_name, label=all_tasks[op_name], is_task=op_name in self.tasks
                )
                self.register_section(op_name)
        for node_1, node_2 in signac_graph.edges:
            if (op_name_1 := operations[node_1]) in graph.nodes and (
                op_name_2 := operations[node_2]
            ) in graph.nodes:
                graph.add_edge(op_name_1, op_name_2)
        return graph

    def build_workflow_yaml(self, destination_filename: str) -> None:
        archive = self.generate_archive()
        archive.to_yaml(destination_filename)
        project_name = self.project.class_name()
        self.job.doc = update_nested_dict(
            self.job.doc, {project_name: {"nomad_workflow": destination_filename}}
        )

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
                if all(self.task_elements[n].is_task for n in in_nodes):
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
    """
    Builds and serializes a NOMAD workflow for a given job.

    This function initializes a `NomadWorkflow` instance with the specified job and configuration flags,
    then serializes the constructed workflow into a YAML file. The YAML file is saved in a location determined
    by whether the workflow is considered top-level or not. If `is_top_level` is True, the workflow is saved
    using the `nomad_top_level_workflow` path from the project; otherwise, it uses the `nomad_workflow` path.
    The `add_job_id` flag determines whether job IDs should be added to section labels to ensure uniqueness.

    Args:
        job (Job): The job instance for which the workflow is being built.
        is_top_level (bool, optional): Flag indicating if the workflow is at the top level of the project hierarchy.
                                       Defaults to False.
        add_job_id (bool, optional): Flag indicating if the job ID should be added to section labels for uniqueness.
                                     Defaults to False.
    """
    project = cast(MartiniFlowProject, job.project)
    workflow = NomadWorkflow(project, job, is_top_level, add_job_id=add_job_id)
    destination_filename = (
        project.nomad_top_level_workflow if is_top_level else project.nomad_workflow
    )
    workflow.build_workflow_yaml(destination_filename)


def build_nomad_workflow_with_multiple_jobs(
    project: MartiniFlowProject, jobs: list[Job]
):
    """
    Builds and serializes a NOMAD workflow archive for multiple jobs within a project.

    This function aggregates the workflows of multiple jobs into a single NOMAD workflow archive. It utilizes the
    `NomadWorkflowArchive.from_multiple_jobs` class method to create an archive that represents the combined workflow
    of the specified jobs. The resulting workflow archive is then serialized into a YAML file, which is saved using
    the `nomad_top_level_workflow` path from the project. This allows for the analysis and reconstruction of complex
    workflows that span multiple computational jobs within the same project.

    Args:
        project (MartiniFlowProject): The project instance containing operations and configurations for the workflow.
        jobs (list[Job]): A list of job instances to be included in the workflow archive.
    """
    archive = NomadWorkflowArchive.from_multiple_jobs(project, jobs)
    destination_filename = project.nomad_top_level_workflow
    archive.to_yaml(destination_filename)

from martignac import config
from martignac.liquid_models.mixtures import LiquidMixture
from martignac.nomad.workflows import Job, build_nomad_workflow
from martignac.parsers.gromacs_topologies import Topology
from martignac.utils.gromacs import gromacs_simulation_command
from martignac.utils.insane import generate_bilayer_with_insane
from martignac.utils.martini_flow_projects import (
    MartiniFlowProject,
    fetched_from_nomad,
    uploaded_to_nomad,
    store_task,
    store_gromacs_log_to_doc,
)
from martignac.utils.misc import sub_template_mdp

conf = config()["bilayer_generation"]


class BilayerGenFlow(MartiniFlowProject):
    workspace_path: str = f"{MartiniFlowProject.workspaces_path}/{conf['relative_paths']['workspaces']}"
    itp_path = f"{MartiniFlowProject.input_files_path}/{conf['relative_paths']['itp_files']}"
    mdp_path = f"{MartiniFlowProject.input_files_path}/{conf['relative_paths']['mdp_files']}"
    itp_files = {k: v.get(str) for k, v in conf["itp_files"].items()}
    mdp_files = {k: v.get(str) for k, v in conf["mdp_files"].items()}
    simulation_settings = {
        "n_threads": conf["settings"]["n_threads"].get(int),
        "box_length_xy": conf["settings"]["box_length_xy"].get(float),
        "box_length_z": conf["settings"]["box_length_z"].get(float),
    }
    solvent: LiquidMixture = LiquidMixture.from_list_of_dicts([s.get(dict) for s in conf["settings"]["solvent"]])
    system_name = conf["output_names"]["system"].get(str)
    nomad_workflow: str = conf["output_names"]["nomad_workflow"].get(str)
    state_names = {k: v.get(str) for k, v in conf["output_names"]["states"].items()}


def lipid_names(job) -> list[str]:
    return [lip.get("name") for lip in job.sp.get("lipids")]


@BilayerGenFlow.label
def generated_gen_gro(job):
    return job.isfile(BilayerGenFlow.get_state_name("generate", "gro"))


@BilayerGenFlow.label
def generated_min_gro(job):
    return job.isfile(BilayerGenFlow.get_state_name("minimize", "gro"))


@BilayerGenFlow.label
def generated_equ_gro(job):
    return job.isfile(BilayerGenFlow.get_state_name("equilibrate", "gro"))


@BilayerGenFlow.label
def generated_prod_gro(job):
    return job.isfile(BilayerGenFlow.get_state_name("production", "gro"))


@BilayerGenFlow.pre(fetched_from_nomad)
@BilayerGenFlow.post(generated_gen_gro, tag="generated_gen_gro")
@BilayerGenFlow.operation_hooks.on_success(store_task)
@BilayerGenFlow.operation(cmd=True, with_job=True)
def generate_initial_bilayer(job):
    lipid_mixture = LiquidMixture.from_list_of_dicts(job.sp.lipids)
    return generate_bilayer_with_insane(
        lipids=lipid_mixture,
        solvent=BilayerGenFlow.solvent,
        box_length_xy=BilayerGenFlow.simulation_settings.get("box_length_xy"),
        box_length_z=BilayerGenFlow.simulation_settings.get("box_length_z"),
        gro_bilayer_gen=BilayerGenFlow.get_state_name("generate", "gro"),
        top_bilayer=BilayerGenFlow.get_state_name(extension="top"),
    )


@BilayerGenFlow.pre(generated_gen_gro, tag="generated_gen_gro")
@BilayerGenFlow.post(generated_min_gro, tag="generated_min_gro")
@BilayerGenFlow.operation_hooks.on_success(store_gromacs_log_to_doc)
@BilayerGenFlow.operation(cmd=True, with_job=True)
def minimize(job):
    top = Topology.parse_top_file(BilayerGenFlow.get_state_name(extension="top"))
    top.includes = BilayerGenFlow.itp_files.values()
    top.output_top(BilayerGenFlow.get_state_name(extension="top"))
    return gromacs_simulation_command(
        mdp=BilayerGenFlow.mdp_files.get("minimize"),
        top=BilayerGenFlow.get_state_name(extension="top"),
        gro=BilayerGenFlow.get_state_name("generate", "gro"),
        name=BilayerGenFlow.get_state_name("minimize"),
        n_threads=BilayerGenFlow.simulation_settings.get("n_threads"),
        verbose=False,
    )


@BilayerGenFlow.pre(generated_min_gro, tag="generated_min_gro")
@BilayerGenFlow.post(generated_equ_gro, tag="generated_equ_gro")
@BilayerGenFlow.operation_hooks.on_success(store_gromacs_log_to_doc)
@BilayerGenFlow.operation(cmd=True, with_job=True)
def equilibrate(job):
    mdp_file_template = BilayerGenFlow.mdp_files.get("equilibrate")
    mdp_file = job.fn(BilayerGenFlow.get_state_name("equilibrate", "mdp"))
    sub_template_mdp(mdp_file_template, "sedlipids", " ".join(lipid_names(job)), mdp_file)
    sub_template_mdp(mdp_file, "sedsolvent", " ".join(BilayerGenFlow.solvent.solvent_names), mdp_file)
    return gromacs_simulation_command(
        mdp=mdp_file,
        top=BilayerGenFlow.get_state_name(extension="top"),
        gro=BilayerGenFlow.get_state_name("minimize", "gro"),
        name=BilayerGenFlow.get_state_name("equilibrate"),
        n_threads=BilayerGenFlow.simulation_settings.get("n_threads"),
    )


@BilayerGenFlow.pre(generated_equ_gro, tag="generated_equ_gro")
@BilayerGenFlow.post(generated_prod_gro, tag="generated_prod_gro")
@BilayerGenFlow.operation_hooks.on_success(store_gromacs_log_to_doc)
@BilayerGenFlow.operation(cmd=True, with_job=True)
def production(job):
    mdp_file_template = BilayerGenFlow.mdp_files.get("production")
    mdp_file = job.fn(BilayerGenFlow.get_state_name("production", "mdp"))
    sub_template_mdp(mdp_file_template, "sedlipids", " ".join(lipid_names(job)), mdp_file)
    sub_template_mdp(mdp_file, "sedsolvent", " ".join(BilayerGenFlow.solvent.solvent_names), mdp_file)
    return gromacs_simulation_command(
        mdp=mdp_file,
        top=BilayerGenFlow.get_state_name(extension="top"),
        gro=BilayerGenFlow.get_state_name("equilibrate", "gro"),
        name=BilayerGenFlow.get_state_name("production"),
        n_threads=BilayerGenFlow.simulation_settings.get("n_threads"),
    )


@BilayerGenFlow.pre(generated_prod_gro, tag="generated_prod_gro")
@BilayerGenFlow.post(lambda job: job.isfile(BilayerGenFlow.nomad_workflow), tag="generated_nomad_workflow")
@BilayerGenFlow.operation(with_job=True)
def generate_nomad_workflow(job: Job):
    build_nomad_workflow(job, is_top_level=False)
    build_nomad_workflow(job, is_top_level=True)


@BilayerGenFlow.pre(lambda job: job.isfile(BilayerGenFlow.nomad_workflow), tag="generated_nomad_workflow")
@BilayerGenFlow.post(uploaded_to_nomad)
@BilayerGenFlow.operation(with_job=True)
def upload_to_nomad(job: Job):
    return BilayerGenFlow.upload_to_nomad(job)


project = BilayerGenFlow.get_project(path=BilayerGenFlow.workspace_path)

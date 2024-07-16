import logging

import gromacs
import numpy as np
from flow import aggregator
from MDAnalysis import Universe

from martignac import config
from martignac.nomad.workflows import Job, build_nomad_workflow
from martignac.parsers.gromacs_topologies import combine_multiple_topology_files
from martignac.utils.gromacs import gromacs_simulation_command, run_gmx_wham
from martignac.utils.martini_flow_projects import (
    MartiniFlowProject,
    fetched_from_nomad,
    flag_ready_for_upload,
    import_job_from_other_flow,
    is_ready_for_upload,
    store_gromacs_log_to_doc_with_depth_from_bilayer_core,
    store_task,
    store_task_for_many_jobs,
    store_workflow,
    system_equilibrated,
    system_generated,
    system_minimized,
    uploaded_to_nomad,
)
from martignac.utils.misc import (
    calculate_average_com,
    convert_pdb_to_gro,
    func_name,
    sub_template_mdp,
    translate_gro_by_vector,
    update_nested_dict,
)
from martignac.utils.packmol import place_solute_in_solvent_with_packmol
from martignac.workflows.system_generation.bilayer import BilayerGenFlow, get_solvent_job
from martignac.workflows.system_generation.solute import SoluteGenFlow, get_solute_job

logger = logging.getLogger(__name__)

conf = config()["solute_in_bilayer"]


class SoluteInBilayerFlow(MartiniFlowProject):
    workspace_path: str = f"{MartiniFlowProject.workspaces_path}/{conf['relative_paths']['workspaces']}"
    itp_files = {k: v.get(str) for k, v in conf["itp_files"].items()}
    mdp_path = f"{MartiniFlowProject.input_files_path}/{conf['relative_paths']['mdp_files']}"
    mdp_files = {k: v.get(str) for k, v in conf["mdp_files"].items()}
    simulation_settings = {"n_threads": conf["settings"]["n_threads"].get(int)}
    system_name = conf["output_names"]["system"].get(str)
    nomad_workflow: str = conf["output_names"]["nomad_workflow"].get(str)
    nomad_top_level_workflow: str = conf["output_names"]["nomad_top_level_workflow"].get(str)
    state_names = {k: v.get(str) for k, v in conf["output_names"]["states"].items()}

    @classmethod
    def get_system_run_name(cls, depth_state: str) -> str:
        return SoluteInBilayerFlow.get_state_name(state_name=f"production-{depth_state}")


project_name = SoluteInBilayerFlow.class_name()
solute_gen_name = SoluteGenFlow.class_name()
bilayer_gen_name = BilayerGenFlow.class_name()


wham_tpr_summary = "tpr-files.dat"
wham_xvg_summary = "pullf-files.dat"
wham_profile_xvg = "wham_profile.xvg"
wham_bstrap_xvg = "wham_bstrap.xvg"
wham_bsprof_xvg = "wham_bsprof.xvg"
wham_hist_xvg = "wham_hist.xvg"
wham_profile = "wham_profile.npy"
wham_bstrap = "wham_bstrap.npy"
wham_bsprof = "wham_bsprof.npy"
wham_hist = "wham_hist.npy"


def solute_in_bilayer_aggregator(jobs):
    unique_combinations = {(job.sp["solute_name"], tuple(tuple(d.items()) for d in job.sp["lipids"])) for job in jobs}
    for solute, lipid_composition in unique_combinations:
        yield [
            job
            for job in jobs
            if job.sp.solute_name == solute and tuple(tuple(d.items()) for d in job.sp.lipids) == lipid_composition
        ]


def job_system_keys(job) -> [dict, str]:
    return [job.sp.lipids, job.sp.solute_name]


def lowest_depth_job(jobs) -> Job:
    return min(jobs, key=lambda j: j.sp.depth_from_bilayer_core)


@SoluteInBilayerFlow.label
def solute_generated(job) -> bool:
    return (
        solute_gen_name in job.doc
        and job.doc[solute_gen_name].get("solute_gro")
        and job.doc[solute_gen_name].get("solute_top")
    )


@SoluteInBilayerFlow.label
def solute_translated(job) -> bool:
    return (
        solute_gen_name in job.doc
        and job.doc[project_name].get("solute_translated_gro")
        and job.doc[solute_gen_name].get("solute_top")
    )


@SoluteInBilayerFlow.label
def generated_box_pdb(job):
    return job.isfile(SoluteInBilayerFlow.get_state_name("generate", "pdb"))


@SoluteInBilayerFlow.label
def bilayer_generated(job) -> bool:
    return bilayer_gen_name in job.doc and job.doc[bilayer_gen_name].get("bilayer_gro")


@SoluteInBilayerFlow.label
def topology_updated(job) -> bool:
    return project_name in job.doc and job.doc[project_name].get("solute_bilayer_top")


@SoluteInBilayerFlow.label
def bilayer_depth_sampled(job) -> bool:
    return job.doc.get(project_name, False) and job.isfile(
        job.fn(job.doc[project_name].get("umbrella_pullx_xvg", default=""))
    )


@SoluteInBilayerFlow.pre(fetched_from_nomad)
@SoluteInBilayerFlow.post(solute_generated, tag="solute_generated")
@SoluteInBilayerFlow.operation_hooks.on_success(store_workflow)
@SoluteInBilayerFlow.operation
def generate_solute(job: Job):
    solute_job: Job = get_solute_job(job)
    keys_for_files_to_copy = ["solute_gro", "solute_top", "solute_itp"]
    project.operation_to_workflow[func_name()] = solute_gen_name
    import_job_from_other_flow(job, solute_gen_project, solute_job, keys_for_files_to_copy)


@SoluteInBilayerFlow.pre(solute_generated, tag="solute_generated")
@SoluteInBilayerFlow.post(solute_translated, tag="solute_translated")
@SoluteInBilayerFlow.operation_hooks.on_success(store_workflow)
@SoluteInBilayerFlow.operation(with_job=True)
def translate_solute(job: Job):
    lipid_names = job.doc[bilayer_gen_name].get("lipid_names")
    com_bilayer = calculate_average_com(job.doc[bilayer_gen_name].get("bilayer_gro"), lipid_names)
    logger.info(f"bilayer COM: {com_bilayer}")
    job.doc[project_name]["bilayer_z_mean"] = com_bilayer[2]
    solute_gro = job.doc[solute_gen_name].get("solute_gro")
    solute_translated_gro = SoluteInBilayerFlow.get_state_name("solute_translated", "gro")
    job.doc[project_name]["solute_translated_gro"] = solute_translated_gro
    com_solute = calculate_average_com(solute_gro, [])
    logger.debug(f"solute COM: {com_solute}")
    diff_solute_com = -com_solute
    logger.info(f"Translating the solute COM by {diff_solute_com}")
    translate_gro_by_vector(solute_gro, solute_translated_gro, diff_solute_com)
    project.operation_to_workflow[func_name()] = project_name


@SoluteInBilayerFlow.pre(fetched_from_nomad)
@SoluteInBilayerFlow.post(bilayer_generated, tag="bilayer_generated")
@SoluteInBilayerFlow.operation_hooks.on_success(store_workflow)
@SoluteInBilayerFlow.operation
def generate_bilayer(job):
    solvent_job: Job = get_solvent_job(job)
    keys_for_files_to_copy = ["bilayer_gro", "bilayer_top"]
    project.operation_to_workflow[func_name()] = bilayer_gen_name
    import_job_from_other_flow(job, bilayer_gen_project, solvent_job, keys_for_files_to_copy)


@SoluteInBilayerFlow.pre(solute_translated, tag="solute_translated")
@SoluteInBilayerFlow.pre(bilayer_generated, tag="bilayer_generated")
@SoluteInBilayerFlow.post(generated_box_pdb, tag="generated_box_pdb")
@SoluteInBilayerFlow.operation_hooks.on_success(store_task)
@SoluteInBilayerFlow.operation(cmd=True, with_job=True)
def insert_solute_in_box(job):
    return place_solute_in_solvent_with_packmol(
        gro_solute=job.doc[project_name]["solute_translated_gro"],
        gro_solvent=job.doc[bilayer_gen_name]["bilayer_gro"],
        output_pdb=SoluteInBilayerFlow.get_state_name("generate", "pdb"),
        restraint_along_z=job.doc[project_name].get("bilayer_z_mean") + job.sp.depth_from_bilayer_core * 10.0,
    )


@SoluteInBilayerFlow.pre(generated_box_pdb, tag="generated_box_pdb")
@SoluteInBilayerFlow.post(system_generated, tag="system_generated")
@SoluteInBilayerFlow.operation_hooks.on_success(store_task)
@SoluteInBilayerFlow.operation(with_job=True)
def convert_box_to_gro(job: Job):
    universe = Universe(job.doc[bilayer_gen_name]["bilayer_gro"])
    box_dimensions = universe.dimensions
    convert_pdb_to_gro(
        SoluteInBilayerFlow.get_state_name("generate", "pdb"),
        SoluteInBilayerFlow.get_state_name("generate", "gro"),
        box_dimensions,
    )


@SoluteInBilayerFlow.pre(system_generated, tag="system_generated")
@SoluteInBilayerFlow.post(topology_updated, tag="topology_updated")
@SoluteInBilayerFlow.operation_hooks.on_success(store_task)
@SoluteInBilayerFlow.operation(with_job=True)
def update_topology_file(job):
    comb_topology = combine_multiple_topology_files(
        [job.doc[bilayer_gen_name]["bilayer_top"], job.doc[solute_gen_name]["solute_top"]], "solute in bilayer"
    )
    top_file_name = SoluteInBilayerFlow.get_state_name(extension="top")
    job.doc[project_name]["solute_bilayer_top"] = top_file_name
    comb_topology.update_counts_against_gro(SoluteInBilayerFlow.get_state_name("generate", "gro"))
    comb_topology.output_top(top_file_name)


@SoluteInBilayerFlow.pre(topology_updated, tag="topology_updated")
@SoluteInBilayerFlow.post(system_minimized, tag="system_minimized")
@SoluteInBilayerFlow.operation_hooks.on_success(store_gromacs_log_to_doc_with_depth_from_bilayer_core)
@SoluteInBilayerFlow.operation(cmd=True, with_job=True)
def minimize(job):
    return gromacs_simulation_command(
        mdp=SoluteInBilayerFlow.mdp_files.get("minimize"),
        top=SoluteInBilayerFlow.get_state_name(extension="top"),
        gro=SoluteInBilayerFlow.get_state_name("generate", "gro"),
        name=SoluteInBilayerFlow.get_state_name("minimize"),
        n_threads=SoluteInBilayerFlow.simulation_settings.get("n_threads"),
    )


@SoluteInBilayerFlow.pre(system_minimized, tag="system_minimized")
@SoluteInBilayerFlow.post(system_equilibrated, tag="system_equilibrated")
@SoluteInBilayerFlow.operation_hooks.on_success(store_gromacs_log_to_doc_with_depth_from_bilayer_core)
@SoluteInBilayerFlow.operation(cmd=True, with_job=True)
def equilibrate(job):
    depth_state = str(job.sp.depth_from_bilayer_core)
    depth_mdp = f"{SoluteInBilayerFlow.get_system_run_name(depth_state)}.mdp"
    sub_template_mdp(SoluteInBilayerFlow.mdp_files.get("equilibrate"), "sedstate", depth_state, depth_mdp)
    sub_template_mdp(depth_mdp, "sedmolecule", job.sp["solute_name"], depth_mdp)
    sub_template_mdp(depth_mdp, "sedlipid", job.sp["lipids"][0]["name"], depth_mdp)
    return gromacs_simulation_command(
        mdp=depth_mdp,
        top=SoluteInBilayerFlow.get_state_name(extension="top"),
        gro=SoluteInBilayerFlow.get_state_name("minimize", "gro"),
        name=SoluteInBilayerFlow.get_state_name("equilibrate"),
        n_threads=SoluteInBilayerFlow.simulation_settings.get("n_threads"),
    )


@SoluteInBilayerFlow.pre(system_equilibrated, tag="system_equilibrated")
@SoluteInBilayerFlow.post(bilayer_depth_sampled, tag="system_sampled")
@SoluteInBilayerFlow.operation_hooks.on_success(store_gromacs_log_to_doc_with_depth_from_bilayer_core)
@SoluteInBilayerFlow.operation(cmd=True, with_job=True)
def production(job):
    depth_state = str(job.sp.depth_from_bilayer_core)
    depth_mdp = f"{SoluteInBilayerFlow.get_system_run_name(depth_state)}.mdp"
    sub_template_mdp(SoluteInBilayerFlow.mdp_files.get("production"), "sedstate", depth_state, depth_mdp)
    sub_template_mdp(depth_mdp, "sedmolecule", job.sp["solute_name"], depth_mdp)
    sub_template_mdp(depth_mdp, "sedlipid", job.sp["lipids"][0]["name"], depth_mdp)
    job.doc[project_name]["umbrella_pullf_xvg"] = f"{SoluteInBilayerFlow.get_system_run_name(depth_state)}_pullf.xvg"
    job.doc[project_name]["umbrella_pullx_xvg"] = f"{SoluteInBilayerFlow.get_system_run_name(depth_state)}_pullx.xvg"
    job.doc[project_name]["tpr_file"] = f"{SoluteInBilayerFlow.get_system_run_name(depth_state)}.tpr"
    job.doc[project_name]["umbrella_log"] = f"{SoluteInBilayerFlow.get_system_run_name(depth_state)}.log"
    return gromacs_simulation_command(
        mdp=depth_mdp,
        top=SoluteInBilayerFlow.get_state_name(extension="top"),
        gro=SoluteInBilayerFlow.get_state_name("equilibrate", "gro"),
        name=SoluteInBilayerFlow.get_system_run_name(depth_state),
        n_threads=SoluteInBilayerFlow.simulation_settings.get("n_threads"),
    )


@SoluteInBilayerFlow.label
def wham_calculated(*jobs) -> bool:
    return any([job.isfile(wham_profile_xvg) and job.isfile(wham_hist_xvg) for job in jobs])


@SoluteInBilayerFlow.label
def free_energy_calculated(*jobs) -> bool:
    return all(
        [
            job.doc.get(project_name, False)
            and job.doc[project_name].get("free_energy", False)
            and job.doc[project_name]["free_energy"].get("profile", False)
            for job in jobs
        ]
    )


@SoluteInBilayerFlow.label
def all_depth_states_sampled(*jobs):
    result = all(
        [
            job.doc.get(project_name, False)
            and job.isfile(job.fn(job.doc[project_name].get("umbrella_pullx_xvg", default="")))
            for job in jobs
        ]
    )
    return result


@SoluteInBilayerFlow.pre(all_depth_states_sampled, tag="system_sampled")
@SoluteInBilayerFlow.post(wham_calculated, tag="wham_calculated")
@SoluteInBilayerFlow.operation_hooks.on_success(store_task_for_many_jobs)
@SoluteInBilayerFlow.operation(
    aggregator=aggregator(aggregator_function=solute_in_bilayer_aggregator, sort_by="depth_from_bilayer_core"),
    cmd=True,
)
def compute_wham(*jobs):
    logger.info("calculation of potential of mean force")
    xvg_files = [job.fn(job.doc[project_name]["umbrella_pullf_xvg"]) for job in jobs]
    tpr_files = [job.fn(job.doc[project_name]["tpr_file"]) for job in jobs]

    def generate_summary_file(output_filename: str, list_of_files: list[str]) -> None:
        with open(output_filename, "w") as file:
            for fil in list_of_files:
                file.write(f"{fil}\n")

    for job in jobs:
        if job == lowest_depth_job(jobs):
            with job:
                generate_summary_file(wham_tpr_summary, tpr_files)
                generate_summary_file(wham_xvg_summary, xvg_files)
                return run_gmx_wham(
                    job.fn(wham_tpr_summary),
                    job.fn(wham_xvg_summary),
                    job.fn(wham_profile_xvg),
                    job.fn(wham_hist_xvg),
                    job.fn(wham_bstrap_xvg),
                    job.fn(wham_bsprof_xvg),
                )


@SoluteInBilayerFlow.pre(wham_calculated, tag="wham_calculated")
@SoluteInBilayerFlow.post(free_energy_calculated, tag="free_energy_calculated")
@SoluteInBilayerFlow.operation_hooks.on_success(store_task_for_many_jobs)
@SoluteInBilayerFlow.operation(
    aggregator=aggregator(aggregator_function=solute_in_bilayer_aggregator, sort_by="depth_from_bilayer_core")
)
def analyze_wham(*jobs):
    for job in jobs:
        if job == lowest_depth_job(jobs):
            xvg_profile = gromacs.fileformats.XVG(job.fn(wham_profile_xvg))
            np.save(job.fn(wham_profile), xvg_profile.array)
            xvg_hist = gromacs.fileformats.XVG(job.fn(wham_hist_xvg))
            np.save(job.fn(wham_hist), xvg_hist.array)
            xvg_bstrap = gromacs.fileformats.XVG(job.fn(wham_bstrap_xvg))
            np.save(job.fn(wham_bstrap), xvg_bstrap.array)
        job.doc = update_nested_dict(
            job.doc,
            {
                project_name: {
                    "free_energy": {
                        "job_id": job.id,
                        "profile": wham_profile,
                        "hist": wham_hist,
                        "bootstrap": wham_bstrap,
                    }
                }
            },
        )


@SoluteInBilayerFlow.pre(free_energy_calculated, tag="free_energy_calculated")
@SoluteInBilayerFlow.post(
    lambda job: all(
        [job.isfile(SoluteInBilayerFlow.nomad_workflow), job.isfile(SoluteInBilayerFlow.nomad_top_level_workflow)]
    ),
    tag="generated_nomad_workflow",
)
@SoluteInBilayerFlow.operation_hooks.on_success(flag_ready_for_upload)
@SoluteInBilayerFlow.operation(with_job=True)
def generate_nomad_workflow(job):
    build_nomad_workflow(job, is_top_level=False)
    build_nomad_workflow(job, is_top_level=True)


@SoluteInBilayerFlow.pre(is_ready_for_upload, tag="generated_nomad_workflow")
@SoluteInBilayerFlow.post(uploaded_to_nomad)
@SoluteInBilayerFlow.operation(with_job=True)
def upload_to_nomad(job: Job):
    return SoluteInBilayerFlow.upload_to_nomad(job)


def get_solvation_job(job: Job) -> Job:
    sp = {
        "type": "solute_solvation",
        "lipids": job.sp.get("lipids"),
        "solute_name": job.sp.get("solute_name"),
        "depth_from_bilayer_core": job.sp.get("depth_from_bilayer_core"),
    }
    return project.open_job(sp).init()


project = SoluteInBilayerFlow.get_project(path=SoluteInBilayerFlow.workspace_path)
solute_gen_project = SoluteGenFlow.get_project(path=SoluteGenFlow.workspace_path)
bilayer_gen_project = BilayerGenFlow.get_project(path=BilayerGenFlow.workspace_path)

import os
import shutil
from time import sleep

import pytest
import signac

from martignac.nomad.uploads import delete_upload, get_upload_by_id


@pytest.mark.order(2)
def test_initialize_workspaces(conf, solute_gen_workspace_path):
    os.makedirs(solute_gen_workspace_path, exist_ok=True)
    shutil.rmtree(solute_gen_workspace_path)
    project = signac.init_project(solute_gen_workspace_path)
    assert project is not None


@pytest.mark.order(3)
def test_run_solute_gen_workspace(conf, solute_gen_workspace_path):
    from martignac.workflows.solute_generation import SoluteGenFlow

    project = SoluteGenFlow.get_project(solute_gen_workspace_path)
    sp = {"type": "solute", "solute_name": "C4"}
    project.open_job(sp).init()
    project.run()
    assert len(project.find_jobs()) == 1


@pytest.mark.order(4)
def test_analyze_solute_gen_workspace(conf, solute_gen_workspace_path):
    from martignac.workflows.solute_generation import SoluteGenFlow

    project = SoluteGenFlow.get_project(solute_gen_workspace_path)
    # There should only be one job
    job = list(project.find_jobs())[0]
    assert job.isfile(job.doc["SoluteGenFlow"]["solute_itp"])
    assert job.isfile(job.doc["SoluteGenFlow"]["solute_top"])
    assert job.isfile(job.doc["SoluteGenFlow"]["solute_gro"])


@pytest.mark.order(5)
def test_send_solute_gen_workspace_to_nomad(conf, global_state, solute_gen_workspace_path):
    from martignac.utils.martini_flow_projects import UPLOAD_TO_NOMAD
    from martignac.workflows.solute_generation import SoluteGenFlow

    project = SoluteGenFlow.get_project(solute_gen_workspace_path)
    print(f"upload to nomad: {UPLOAD_TO_NOMAD}")
    print(f"global state: {global_state}")
    SoluteGenFlow.nomad_dataset_id = global_state["dataset_id"]
    SoluteGenFlow.nomad_use_prod_database = False
    print(f"NOMAD dataset ID: {SoluteGenFlow.nomad_dataset_id}")
    job = list(project.find_jobs())[0]
    SoluteGenFlow.upload_to_nomad(job, True, publish_flag=False)
    sleep(2.0)  # wait for processing to finish
    upload_id = job.doc["SoluteGenFlow"]["nomad_upload_id"]
    nomad_upload = get_upload_by_id(upload_id, use_prod=SoluteGenFlow.nomad_use_prod_database)
    assert nomad_upload is not None
    delete_upload(upload_id, use_prod=SoluteGenFlow.nomad_use_prod_database)

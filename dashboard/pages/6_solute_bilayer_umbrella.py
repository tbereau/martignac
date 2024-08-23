import os
import sys

current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)

from dataclasses import asdict

import pandas as pd
import streamlit as st
import streamlit.components.v1 as components
from init import init_for_streamlit

st.set_page_config(
    page_title="Martignac Solute-in-bilayer umbrella sampling",
    page_icon="🧪",
    layout="wide",
)

init_for_streamlit()

from init2 import paths_for_streamlit

from martignac.liquid_models.mixtures import LiquidMixture
from martignac.nomad.entries import get_entry_by_id
from martignac.nomad.mini_entries import find_mini_queries_corresponding_to_workflow
from martignac.utils.dashboard import generate_gravis_network
from martignac.workflow_interfaces.utils import (
    convert_multiple_entry_ids_to_specific_interfaces,
)
from martignac.workflows.solute_in_bilayer_umbrella import project

paths_for_streamlit()

st.title("Martignac")

st.subheader("Solute-in-bilayer umbrella sampling")

components.html(generate_gravis_network(project).to_html(), height=400)

database_id = st.text_input("Database ID", "HJdEI1q4SV-c5Di43BTT_Q")
prod_db = st.toggle("Production database", False)
st.button("Re-run")

with st.spinner("Querying NOMAD..."):
    df = pd.json_normalize(
        [
            asdict(e)
            for e in find_mini_queries_corresponding_to_workflow(
                project, dataset_ids=[database_id], use_prod=prod_db
            )
        ]
    )

if len(df) == 0:
    st.warning(
        f"No solute-in-bilayer umbrella entry found in database ID '{database_id}' on {'prod' if prod_db else 'test'} server",
        icon="⚠️",
    )
else:
    df = df.rename(
        columns={
            "comment.job_id": "job_id",
            "comment.state_point.lipids": "lipids",
            "comment.state_point.depth_from_bilayer_core": "depth_from_bilayer_core",
            "comment.state_point.solute_name": "solute_name",
            "comment.mdp_files": "mdp_files",
            "comment.itp_files": "itp_files",
        }
    )
    df["lipids"] = df["lipids"].apply(
        lambda x: LiquidMixture.from_list_of_dicts(x).to_insane_format()
    )
    df = (
        df.sort_values("depth_from_bilayer_core")
        .groupby(by=["solute_name", "lipids"])
        .first()
    )
    df = df.drop(
        [
            "workflow_name",
            "comment.workflow_name",
            "datasets",
            "comment.state_point.type",
            "depth_from_bilayer_core",
        ],
        axis=1,
    )
    tuple_of_entry_ids = tuple(df["entry_id"].values)
    with st.spinner("Querying NOMAD..."):
        interfaces = convert_multiple_entry_ids_to_specific_interfaces(
            tuple_of_entry_ids, use_prod=project.nomad_use_prod_database
        )
        df["pmf"] = [i.get_wham_npy(use_bootstrap=False)[1] for i in interfaces]
        df["nomad_url"] = df["entry_id"].apply(
            lambda x: get_entry_by_id(
                x, use_prod=project.nomad_use_prod_database
            ).nomad_gui_url
        )
    column_nomad_url = df.pop("nomad_url")
    column_pmf = df.pop("pmf")
    df.insert(0, "pmf", column_pmf)
    df.insert(1, "nomad_url", column_nomad_url)

    st.dataframe(
        df,
        column_config={
            "nomad_url": st.column_config.LinkColumn("nomad_url", display_text="link"),
            "pmf": st.column_config.LineChartColumn(
                "pmf",
                width="medium",
                help="Potential of mean force (PMF) in kT",
            ),
        },
    )

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
    page_title="Martignac Solute generation", page_icon="🧪", layout="wide"
)

init_for_streamlit()

from init2 import paths_for_streamlit

from martignac.nomad.entries import get_entry_by_id
from martignac.nomad.mini_entries import find_mini_queries_corresponding_to_workflow
from martignac.utils.dashboard import generate_gravis_network
from martignac.workflows.solute_generation import project

paths_for_streamlit()

st.title("Martignac")

st.subheader("Solute generation")

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
        f"No solute generation entry found in database ID '{database_id}' on {'prod' if prod_db else 'test'} server",
        icon="⚠️",
    )
else:
    df = df.rename(
        columns={
            "comment.state_point.solute_name": "solute_name",
            "comment.job_id": "job_id",
            "comment.mdp_files": "mdp_files",
            "comment.itp_files": "itp_files",
        }
    )
    df = df.sort_values("solute_name", ignore_index=True).drop(
        [
            "workflow_name",
            "comment.workflow_name",
            "datasets",
            "comment.state_point.type",
        ],
        axis=1,
    )
    with st.spinner("Querying NOMAD..."):
        df["nomad_url"] = df["entry_id"].apply(
            lambda x: get_entry_by_id(x, use_prod=prod_db).nomad_gui_url
        )
    df = df.set_index("solute_name")
    column_nomad_url = df.pop("nomad_url")
    df.insert(0, "nomad_url", column_nomad_url)

    st.dataframe(
        df,
        column_config={
            "nomad_url": st.column_config.LinkColumn("nomad_url", display_text="link")
        },
    )

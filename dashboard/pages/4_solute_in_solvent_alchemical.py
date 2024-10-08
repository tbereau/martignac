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
    page_title="Martignac Solute-in-solvent alchemical transformation",
    page_icon="🧪",
    layout="wide",
)

init_for_streamlit()

from init2 import paths_for_streamlit

from martignac.nomad.entries import get_multiple_entries_by_id
from martignac.nomad.mini_entries import find_mini_queries_corresponding_to_workflow
from martignac.utils.dashboard import generate_gravis_network
from martignac.workflow_interfaces.utils import (
    convert_multiple_entry_ids_to_specific_interfaces,
)
from martignac.workflows.solute_in_solvent_alchemical import project

paths_for_streamlit()

st.title("Martignac")

st.subheader("Solute-in-solvent alchemical transformation")

components.html(generate_gravis_network(project).to_html(), height=400)

database_id = st.text_input("Database ID", "Q6CrHWQzS2CYi8NwQDw5ng")
prod_db = st.toggle("Production database", True)
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
        f"No alchemical entry found in database ID '{database_id}' on {'prod' if prod_db else 'test'} server",
        icon="⚠️",
    )
else:
    df = df.rename(
        columns={
            "comment.job_id": "job_id",
            "comment.state_point.lambda_state": "lambda_state",
            "comment.state_point.solute_name": "solute_name",
            "comment.state_point.solvent_name": "solvent_name",
            "comment.mdp_files": "mdp_files",
            "comment.itp_files": "itp_files",
        }
    )
    df = (
        df.sort_values("lambda_state")
        .groupby(by=["solute_name", "solvent_name"])
        .first()
    )
    df = df.drop(
        [
            "workflow_name",
            "comment.workflow_name",
            "datasets",
            "comment.state_point.type",
            "lambda_state",
        ],
        axis=1,
    )
    with st.spinner("Querying NOMAD..."):
        tuple_of_entry_ids = tuple(df["entry_id"].values)
        interfaces = convert_multiple_entry_ids_to_specific_interfaces(
            tuple_of_entry_ids,
            use_prod=prod_db,
        )
        df["free_energy_in_kt"] = [i.free_energy.mean for i in interfaces]
        nomad_entries = get_multiple_entries_by_id(
            tuple_of_entry_ids,
            use_prod=prod_db,
        )
        df["nomad_url"] = [e.nomad_gui_url for e in nomad_entries]
    column_free_energy = df.pop("free_energy_in_kt")
    column_nomad_url = df.pop("nomad_url")
    df.insert(0, "free_energy_in_kt", column_free_energy)
    df.insert(1, "nomad_url", column_nomad_url)

    st.dataframe(
        df,
        column_config={
            "nomad_url": st.column_config.LinkColumn("nomad_url", display_text="link")
        },
    )

    # df = df.reset_index(level=[0, 1])
    #
    # st.data_editor(
    #     styled_df,
    #     column_config={
    #         "nomad_url": st.column_config.LinkColumn(
    #             "nomad_url",
    #             help="Hyperlink to NOMAD GUI",
    #             display_text="url",
    #         )
    #     },
    # )
    st.info(
        "Column `free_energy_in_kt` displays the solvation free energy of `solute_name` "
        r"in `solvent_name` in units of $k_\mathrm{B}T$.",
    )

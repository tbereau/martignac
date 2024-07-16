from dataclasses import asdict

import pandas as pd
import streamlit as st
import streamlit.components.v1 as components

from martignac.nomad.entries import (
    NomadEntry,
    get_entries_of_my_uploads,
)
from martignac.utils.dashboard import generate_gravis_network
from martignac.workflows.system_generation.solute_in_bilayer import project as solute_bilayer_project

st.set_page_config(page_title="Solute-in-bilayer umbrella sampling", page_icon="📊")

st.title("Martignac: Solute-in-bilayer umbrella sampling")

components.html(generate_gravis_network(solute_bilayer_project).to_html(), height=400)

# uploads
upload_entries: list[NomadEntry] = get_entries_of_my_uploads(use_prod=False)
df_u = pd.json_normalize(
    [
        {
            **asdict(e),
            "workflow_name": e.workflow_name,
            "state_point": e.state_point,
            "mdp_files": e.mdp_files,
            "url": e.nomad_gui_url,
        }
        for e in upload_entries
        if e.comment and e._comment_dict.get("state_point") and e._comment_dict.get("workflow_name")
    ]
)
sp_columns = [c for c in df_u.columns if "state_point" in c]
df_u = df_u[
    [
        "entry_id",
        "upload_id",
        "entry_name",
        "entry_type",
        "workflow_name",
        "mdp_files",
        *sp_columns,
        "entry_create_time",
        "url",
    ]
]
df_u["entry_create_time"] = pd.to_datetime(df_u["entry_create_time"]).dt.date

df_u = (
    df_u[(df_u["workflow_name"] == "SoluteInBilayerFlow") & (df_u["entry_type"] == "Workflow")]
    .drop(
        [
            "workflow_name",
            "state_point.type",
            "state_point.solvent_name",
            "state_point.lambda_state",
        ],
        axis=1,
    )
    .sort_values("state_point.lipids")
)


# Expand the 'state_point.lipids' column
def expand_lipids(lipids):
    names = [lipid["name"] for lipid in lipids]
    fractions = [lipid["fraction"] for lipid in lipids]
    return pd.Series([names, fractions])


# Apply the function to the 'state_point.lipids' column
df_expanded = df_u["state_point.lipids"].apply(expand_lipids)
df_expanded.columns = ["state_point.lipids.names", "state_point.lipids.fractions"]

# Concatenate the expanded DataFrame with the original DataFrame (if you want to keep other columns)
df_u = pd.concat([df_u.drop(columns=["state_point.lipids"]), df_expanded], axis=1)

st.dataframe(df_u, column_config={"url": st.column_config.LinkColumn("URL", display_text="link")})

st.button("Re-run")

from dataclasses import asdict

import pandas as pd
import streamlit as st

from martignac.nomad.entries import (
    NomadEntry,
    get_entries_of_my_uploads,
)

st.set_page_config(page_title="Solute generation", page_icon="📊")

st.title("Martignac: Solute generation")

# uploads
upload_entries: list[NomadEntry] = get_entries_of_my_uploads()
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
    df_u[(df_u["workflow_name"] == "SoluteGenFlow") & (df_u["entry_type"] == "Workflow")]
    .drop(
        [
            "workflow_name",
            "state_point.type",
            "state_point.solvent_name",
            # "state_point.lambda_state"
        ],
        axis=1,
    )
    .sort_values("state_point.solute_name")
)

st.dataframe(df_u, column_config={"url": st.column_config.LinkColumn("URL", display_text="link")})

st.button("Re-run")

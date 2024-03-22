from dataclasses import asdict

import pandas as pd
import streamlit as st

from martignac.nomad.entries import (
    get_entries_in_database,
    get_entries_of_my_uploads,
)

st.set_page_config(
    page_title="Martignac Home",
    page_icon="👋",
)

st.title("Martignac dataset")

entries = get_entries_in_database()
df_e = pd.json_normalize([{**asdict(e), "url": e.nomad_gui_url} for e in entries])
df_e = df_e[["entry_id", "upload_id", "entry_type", "entry_create_time", "url"]]

# uploads
upload_entries = get_entries_of_my_uploads()
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

st.subheader("Entries")
st.dataframe(df_e, column_config={"url": st.column_config.LinkColumn("URL", display_text="link")})

st.subheader("Uploads")
st.dataframe(df_u, column_config={"url": st.column_config.LinkColumn("URL", display_text="link")})

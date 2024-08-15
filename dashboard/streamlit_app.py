from dataclasses import asdict

import pandas as pd
import streamlit as st
from init import init_for_streamlit

st.set_page_config(
    page_title="Martignac Home",
    page_icon="👋",
)

init_for_streamlit()

from init2 import paths_for_streamlit

from martignac.nomad.entries import (
    get_entries_in_database,
    get_entries_of_my_uploads,
)

paths_for_streamlit()

st.title("Martignac dataset")

database_id = st.text_input("Database ID", "HJdEI1q4SV-c5Di43BTT_Q")
prod_db = st.toggle("Production database", False)

df_e = pd.DataFrame()
entries = get_entries_in_database(database_id=database_id, use_prod=prod_db)
if not (entries):
    st.warning(
        f"No entry found in database ID '{database_id}' on {'prod' if prod_db else 'test'} server",
        icon="⚠️",
    )
else:
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
        if e.comment
        and e._comment_dict.get("state_point")
        and e._comment_dict.get("workflow_name")
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
st.dataframe(
    df_e, column_config={"url": st.column_config.LinkColumn("URL", display_text="link")}
)

st.subheader("Uploads")
st.dataframe(
    df_u, column_config={"url": st.column_config.LinkColumn("URL", display_text="link")}
)

from dataclasses import asdict

import pandas as pd
import streamlit as st
from init import init_for_streamlit

st.set_page_config(page_title="Martignac Home", page_icon="🧪", layout="wide")

init_for_streamlit()

from init2 import paths_for_streamlit

from martignac.nomad.datasets import get_dataset_by_id
from martignac.nomad.queries import query_entries_in_mini_format

paths_for_streamlit()

st.title("Martignac")

"Welcome to the Martignac streamlit app! Here you can explore the datasets and entries in the Martignac database."

st.link_button("Martignac's Github repository", "https://github.com/tbereau/martignac/")

st.subheader("Martignac datasets")

"Toggle the switch to see the Martignac datasets in the NOMAD production or test databases."

prod_db = st.toggle("Production database", False)


with st.spinner("Querying NOMAD..."):
    df_e = pd.json_normalize(
        [asdict(e) for e in query_entries_in_mini_format(use_prod=prod_db)]
    )

    dataset_ids = (
        list({e["dataset_id"] for d in df_e["datasets"] if d for e in d})
        if len(df_e) > 0
        else []
    )
    df = pd.json_normalize(
        [asdict(get_dataset_by_id(d, use_prod=prod_db)) for d in dataset_ids]
    )

if len(df_e) > 0:
    df = df.drop(
        [
            "dataset_type",
            "dataset_modified_time",
            "pid",
            "m_annotations",
            "user.first_name",
            "user.last_name",
            "user.username",
            "user.affiliation_address",
            "user.email",
            "user.is_oasis_admin",
            "user.repo_user_id",
            "user.created",
            "user.is_admin",
        ],
        axis=1,
    )
    df["nomad_url"] = df["dataset_id"].apply(
        lambda x: get_dataset_by_id(x, use_prod=prod_db).nomad_gui_url
    )
    df = df.set_index("dataset_id")
    column_nomad_url = df.pop("nomad_url")
    column_dataset_name = df.pop("dataset_name")
    df.insert(0, "nomad_url", column_nomad_url)
    df.insert(1, "dataset_name", column_dataset_name)

col1, col2 = st.columns(2)
col1.metric("Number of NOMAD datasets", len(df))
col2.metric("Number of NOMAD workflow entries", len(df_e))

if len(df_e) == 0:
    st.warning(
        f"No Martignac dataset found on {'prod' if prod_db else 'test'} server",
        icon="⚠️",
    )
else:
    st.dataframe(
        df,
        column_config={
            "nomad_url": st.column_config.LinkColumn("URL", display_text="link")
        },
    )

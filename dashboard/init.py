from martignac import config


def init_for_streamlit():
    conf = config()
    conf.set_file("/mount/src/martignac/dashboard/config_st.yaml")

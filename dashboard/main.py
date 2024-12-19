import streamlit as st

def streamlit_main():
    st.set_page_config(
        page_title="Chromatography Simulation",
        page_icon="🧊",
        layout="wide",
        initial_sidebar_state="expanded",
        menu_items={
            "About": "### Chromatography Simulation Tool v0.1\n Created by Benedict Hasenauer, 2020"
        },
    )



if __name__=="__main__":
    streamlit_main()
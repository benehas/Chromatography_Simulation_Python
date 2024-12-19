import streamlit as st
import pandas as pd
st.title("Results")

if 'results' in st.session_state:
    x = dict()
    for c in range(st.session_state['n_components']):
        x[f'Component {c}'] = st.session_state.results[c]
    df = pd.DataFrame(x)
    st.line_chart(df,use_container_width=True)

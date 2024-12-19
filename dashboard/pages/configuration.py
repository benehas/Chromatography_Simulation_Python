import streamlit as st
import pandas as pd
import sys
import os

# Add the parent directory to the system path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import simulation_core


def create_config_dicts(df):
    x = dict()
    for i in range(len(df)):
        x[df.Parameter[i]] = df.Value[i]
    return x

def run_simulation():
    match st.session_state['isotherme']:
        case 'Henry':
            _isotherme = simulation_core.IsothermeType.Henry
        case 'Langmuir':
            _isotherme = simulation_core.IsothermeType.Lagmuir
        case 'SMA':
            _isotherme = simulation_core.IsothermeType.SMA
    st.session_state['results'] = simulation_core.ChromatographySimulator(st.session_state['column_config_dict'],
                                                                          st.session_state['numeric_config_dict'],
                                                                          _isotherme,
                                                                          st.session_state['n_components'],
                                                                          st.session_state['component_config_dicts']).simulate()


st.title("Configuration")
col1, col2, col3 = st.columns(3)
with col1:
    st.subheader("Column Parameters")
    df_column_config = pd.DataFrame(
        [
            {"Parameter": 'L', "Value": 15},
            {"Parameter": 'tmax', "Value": 30},
            {"Parameter": 'u_inter', "Value": 1.5},
            {"Parameter": 'epsilon', "Value": 0.72},
            {"Parameter": 't_slug', "Value": 2},
            {"Parameter": 'keff', "Value": 0.01},
            {"Parameter": 'Rp', "Value": 0.0025},
            {"Parameter": 'Lambda', "Value": 2},

        ]
    )
    if 'column_config_df' in st.session_state:
        df_column_config = st.session_state['column_config_df']
    edited_df_column_config = st.data_editor(df_column_config,disabled=["Parameter"],use_container_width=True, hide_index=True)
    st.session_state['column_config_df'] = edited_df_column_config
    st.session_state['column_config_dict'] = create_config_dicts(edited_df_column_config)

with col2:
    st.subheader("Numeric Parameters")
    df_numeric_config = pd.DataFrame(
        [
            {"Parameter": 'dt', "Value": 0.005},
            {"Parameter": 'dx', "Value": 0.1},

        ]
    )
    if 'edited_df_numeric_config' in st.session_state:
        df_numeric_config = st.session_state['edited_df_numeric_config']

    edited_df_numeric_config = st.data_editor(df_numeric_config,disabled=["Parameter"],
                                              use_container_width=True, hide_index=True)

    st.session_state['edited_df_numeric_config'] = edited_df_numeric_config
    st.session_state['numeric_config_dict'] = create_config_dicts(edited_df_numeric_config)

with col3:
    st.subheader("Component Parameters")
    isotherme = st.selectbox('Isotherme Model',('Henry','Langmuir','SMA'))
    st.session_state['isotherme'] = isotherme

    if 'n_components' in st.session_state:
        n_components_default = st.session_state['n_components']
    else:
        n_components_default = 2
    C = st.number_input('Number of Components', min_value=2, value=n_components_default, step=1)
    st.session_state['n_components'] = C

    edited_df_component_config = list()
    if 'component_config_dicts' not in st.session_state:
        st.session_state['component_config_dicts'] = list()
    for i in range(C):
        st.text(f'Component {i+1}')
        if isotherme == 'Henry':
            df_component_config = pd.DataFrame(
                [
                    {"Parameter": 'Kh', "Value": 1.0},
                    {"Parameter": 'c_in', "Value": 1.0}
                ]
            )
        if isotherme == 'Langmuir':
            df_component_config = pd.DataFrame(
                [
                    {"Parameter": 'KL', "Value": 1.0},
                    {"Parameter": 'qmax', "Value": 0.5},
                    {"Parameter": 'c_in', "Value": 1.0},
                ]
            )
        if isotherme == 'SMA':
            df_component_config = pd.DataFrame(
                [
                    {"Parameter": 'Ksma', "Value": 3.0},
                    {"Parameter": 'sigma', "Value": 2.0},
                    {"Parameter": 'ny', "Value": 4.0},
                    {"Parameter": 'c_in', "Value": 1.0},
                ]
            )
        if 'component_configs' in st.session_state:
            if len(st.session_state['component_configs']) - 1 > i:
               df_component_config = st.session_state['component_configs'][i]
        edited_df_component_config.append(st.data_editor(df_component_config,disabled=["Parameter"],
                                                         key=f'Component_{i}',
                                                         use_container_width=True,
                                                         hide_index=True))

    st.session_state['component_config_dicts'] = list()
    st.session_state['component_configs'] = list()

    for i in range(len(edited_df_component_config)):
        st.session_state['component_config_dicts'].append(create_config_dicts(edited_df_component_config[i]))
        st.session_state['component_configs'].append(edited_df_component_config[i])

st.button('Start Simulation', use_container_width=True,on_click=run_simulation)
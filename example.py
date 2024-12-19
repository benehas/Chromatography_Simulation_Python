from dashboard.simulation_core import *
# column Parameters
column_para = {}
column_para['L'] = 15  # cm
column_para['tmax'] = 30  # s
column_para['u_inter'] = 1.5  # cm/s
column_para['epsilon'] = 0.72
column_para['t_slug'] = 2
column_para['keff'] = 0.01
column_para['Rp'] = 0.0025
column_para['Lambda'] = 2

numeric_para = dict()
# Numeric Parameters
numeric_para['dt'] = 0.005  # s
numeric_para['dx'] = 0.1  # cm

# Components
isotherme = 'Henry'
C = 2  # Number of Components
params = np.ndarray(C, dict)

if isotherme == 'Langmuir':
    params[0] = {'KL': 1, 'qmax': 0.5, 'c_in': 0.2}
    params[1] = {'KL': 0.1, 'qmax': 0.5, 'c_in': 0.2}

if isotherme == 'Henry':
    params[0] = {'Kh': 1, 'c_in': 1}
    params[1] = {'Kh': 0.1, 'c_in': 1}

if isotherme == 'SMA':
    params[0] = {'Ksma': 3, 'sigma': 2, 'ny': 2, 'c_in': 2}
    params[1] = {'Ksma': 1, 'sigma': 4, 'ny': 3, 'c_in': 2}



c = ChromatographySimulator(column_para,numeric_para,IsothermeType.Henry,2,params)
c.simulate()
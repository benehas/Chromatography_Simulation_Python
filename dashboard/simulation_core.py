import numpy as np
import matplotlib.pyplot as plt
import copy as cp
import time
from enum import Enum


class IsothermeType(Enum):
    Henry = 'Henry'
    Langmuir = 'Langmuir'
    SMA = 'SMA'

class ChromatographySimulator:
    def __init__(self, column_parameters: dict, numeric_parameters: dict, isotherme: IsothermeType, n_components: int,
                 isotherme_parameters: dict):
        """

        :param column_parameters:
        :type column_parameters:
        :param numeric_parameters:
        :type numeric_parameters:
        :param isotherme:
        :type isotherme:
        :param n_components:
        :type n_components:
        :param isotherme_parameters:
        :type isotherme_parameters:
        """
        self.column_parameters = column_parameters
        self.numeric_parameters = numeric_parameters
        self.isotherme = isotherme
        self.n_components = n_components
        self.isotherme_parameters = isotherme_parameters

        self.const_1 = (1 - column_parameters['epsilon']) / column_parameters['epsilon']  # constants to reduce computations
        self.const_2 = 3 / column_parameters['Rp']  # constants to reduce computations
        self.D_ax = 0.05  # cm^2/s

        self.N = int(column_parameters['L'] / numeric_parameters['dx'])
        self.M = int(column_parameters['tmax'] / numeric_parameters['dt'])

    def compute_dc_dx(self, c_i, i):
        return (c_i[i + 1] - c_i[i - 1]) / (2 * self.numeric_parameters['dx'])

    def compute_dc2_dx2(self, c_i, i):
        return (c_i[i + 1] + c_i[i - 1] - 2 * c_i[i]) / (self.numeric_parameters['dx'] * self.numeric_parameters['dx'])

    def compute_cp_LG(self, c, q, i, params):
        A = np.ndarray((self.n_components, self.n_components))
        D = np.ndarray(self.n_components)

        for k in range(self.n_components):
            for l in range(self.n_components):
                if k == l:
                    A[k][l] = params[l]['KL'] * (q[k][i] - params[l]['qmax'])
                else:
                    A[k][l] = q[k][i] * params[l]['KL']
            D[k] = -q[k][i]
        return D.dot(np.linalg.inv(A))

    def compute_cp_SMA(self,q, k, i, params, csalt):
        tmp = 0
        for r in range(self.n_components):
            tmp += q[r][i] * (params[r]['ny'] + params[r]['sigma'])
        return (q[k][i] / params[k]['Ksma']) * (csalt / (self.column_parameters['Lambda'] - tmp)) ** params[k]['ny']

    def simulate(self):
        concentration = np.zeros((self.n_components, self.N))
        q = np.zeros((self.n_components, self.N))
        if self.isotherme == IsothermeType.SMA:
            c_salt = np.zeros(self.N)
            c_salt_save = np.zeros(self.N)
        concentration_s = np.zeros((self.n_components, self.N))
        q_s = np.zeros((self.n_components, self.N))
        plot = np.zeros((self.n_components, self.M))
        t = [i * self.numeric_parameters['dt'] for i in range(self.M)]

        start = time.time()

        # Initial conditions
        print(self.isotherme_parameters)
        for k in range(self.n_components):
            concentration[k][0] = self.isotherme_parameters[k]['c_in']
        if self.isotherme == IsothermeType.SMA:
            c_salt[0] = 1

        # Main Looper
        for i in range(self.M):
            concentration_s = cp.deepcopy(concentration)
            q_s = cp.deepcopy(q)

            if self.isotherme == IsothermeType.SMA:
                c_salt_save = cp.deepcopy(c_salt)

            if self.isotherme == IsothermeType.Langmuir:
                for l in range(1, self.N - 1):
                    tmp = self.compute_cp_LG(concentration_s, q_s, l, self.isotherme_parameters)
                    for k in range(self.n_components):
                        concentration[k][l] += (self.D_ax * self.compute_dc2_dx2(concentration_s[k], l) -
                                                self.column_parameters['u_inter'] * self.compute_dc_dx(
                            concentration_s[k], l) - self.const_1 * self.const_2 * self.column_parameters['keff'] *
                                                (concentration_s[k][l] - tmp[k])) * self.numeric_parameters['dt']

                        q[k][l] += (self.const_1 * self.const_2 * self.column_parameters['keff'] *
                                    (concentration_s[k][l] - tmp[k]) * self.numeric_parameters['dt'])

            if self.isotherme == IsothermeType.Henry:
                for l in range(1, self.N - 1):
                    for k in range(self.n_components):
                        concentration[k][l] += (self.D_ax * self.compute_dc2_dx2(concentration_s[k], l) - self.column_parameters['u_inter'] * self.compute_dc_dx(
                            concentration_s[k], l) - self.const_1 * self.const_2 * self.column_parameters['keff'] * (
                                                        concentration_s[k][l] - q[k][l] / self.isotherme_parameters[k]['Kh'])) * self.numeric_parameters['dt']
                        q[k][l] += self.const_1 * self.const_2 * self.column_parameters['keff'] * (concentration_s[k][l] - q[k][l] / self.isotherme_parameters[k]['Kh']) * self.numeric_parameters['dt']

            if self.isotherme == IsothermeType.SMA:
                for l in range(1, self.N - 1):
                    c_salt[l] += (self.D_ax * self.compute_dc2_dx2(c_salt_save, l) - self.column_parameters['u_inter'] *
                                  self.compute_dc_dx(c_salt_save, l)) * self.numeric_parameters['dt']

                    for k in range(self.n_components):
                        concentration[k][l] += (self.D_ax * self.compute_dc2_dx2(concentration_s[k], l) -
                                                self.column_parameters['u_inter'] * self.compute_dc_dx(
                            concentration_s[k], l) - self.const_1 * self.const_2 * self.column_parameters['keff'] * (
                                                        concentration_s[k][l] - self.compute_cp_SMA(q, k, l, self.isotherme_parameters,
                                                                                               c_salt[l]))) * self.numeric_parameters['dt']
                        q[k][l] += self.const_1 * self.const_2 * self.column_parameters['keff'] * (
                                concentration_s[k][l] - self.compute_cp_SMA(q, k, l, self.isotherme_parameters, c_salt[l])) * self.numeric_parameters['dt']

            # boundary conditions
            for k in range(self.n_components):
                if i * self.numeric_parameters['dt'] < self.column_parameters['t_slug']:
                    concentration[k][0] = self.isotherme_parameters[k]['c_in']
                else:
                    concentration[k][0] = 0
                concentration[k][-1] = concentration[k][-2]
            if self.isotherme == 'SMA':
                c_salt[-1] = c_salt[-2]

            # add detector output for visualization

            for k in range(self.n_components):
                plot[k][i] = concentration[k][-2]
            if 100 * i / self.M % 10 == 0:
                print('Simulation progress: ' + str(100 * i / self.M) + '%')
                print('Simulation time: ' + str(i * self.numeric_parameters['dt']) + 's')
                print('Phys time: ' + str(time.time() - start) + 's')



        print('Total Computaion time: ' + str(time.time() - start))
        # # visualization
        # plt.xlabel('time in s')
        # plt.ylabel('C')
        # for k in range(self.n_components):
        #     plt.plot(t, plot[k, :])
        # plt.show()
        return plot



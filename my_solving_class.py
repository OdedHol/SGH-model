import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve
import seaborn as sns
sns.set(style="ticks", font_scale=1.5)
import os
import json
from datetime import datetime
import inspect
import os
from my_model_class import Model


class Solving:
    def __init__(self, model):
        self.model = model
        # inst1=inst()
        pass
    
    def solve_beta_S(self, S_values=np.linspace(0, 1, 100)):
        # s_values = np.linspace(0, 1, 100)
        beta_S_values = []
        for s in S_values:
            self.model.y = (self.model.y[0], s)
            beta_S_values.append(self.model.beta_S())
        return beta_S_values
        
    def solve_treewater(self, TC_values, S_border=np.linspace(0.3,0.6,100)):
        # s_values = np.linspace(0, 1, 100)
        s_values = S_border
        tree_water_values = []
        self.model.params['TC'] = TC_values
        for s in s_values:
            # if s < 0.3:
            #     tree_water_values.append(self.model.params['lamda_min']) 
            # elif s > 0.6:
            #     tree_water_values.append(self.model.params['lamda_max'])
            # else:
            self.model.y = (self.model.y[0], s)
            tree_water_values.append(self.model.tree_water())
        return tree_water_values

    def S_of_P(self, P_values, TC_values):
        # params_backup = self.params.copy()
        S = np.zeros_like(P_values)
        self.model.params['TC'] = TC_values
        for i, p in enumerate(P_values):
            self.model.params['P'] = p
            _, S[i] = self.model.solutions_fsolve()
        # self.params = params_backup
        return S

    def B_of_P(self, P_values, TC_values, betaA=[], betaT=[], K_values=[]):
        B = np.zeros_like(P_values)
        self.model.params['TC'] = TC_values
        for i, p in enumerate(P_values):
            self.model.params['P'] = p
            if K_values:
                self.model.params['K'] = K_values
            if betaA and betaT:
                for index_a, ba in enumerate(betaA):
                    for index_t, bt in enumerate(betaT):
                        self.model.params['beta_A'] = ba
                        self.model.params['beta_T'] = bt
                        B[i], _ = self.model.solutions_fsolve()
                        # B[i], _ = self.model.solutions_sympy()
                   
            else:
                B[i], _ = self.model.solutions_fsolve()

        return B
    
    def solve_sh_of_TC(self, betaA=[], betaT=[]):
        tc = np.linspace(0, 1, 100)
        self.model.params['TC'] = tc
        # if betaA and betaT:
        #         for index_a, ba in enumerate(betaA):
        #             for index_t, bt in enumerate(betaT):
        #                 self.model.params['beta_A'] = ba
        #                 self.model.params['beta_T'] = bt
        #                 assimilation = self.model.shading_assimilation()
        #                 transpiration = self.model.shading_transpiration()
        # else:
        assimilation = self.model.shading_assimilation()
        transpiration = self.model.shading_transpiration()
        return assimilation, transpiration

    def B_of_P_Tc(self, K_values=[], betaA=[], betaT=[]):
        P_values = np.linspace (0,800/180,50)
        TC_values = np.linspace (0,1,50)
        if isinstance(K_values, list):
            self.model.params['K'] = K_values[0]
        else:
            self.model.params['K'] = K_values

        if betaA and betaT:
            for index_a, ba in enumerate(betaA):
                for index_t, bt in enumerate(betaT):
                    self.model.params['beta_A'] = ba
                    self.model.params['beta_T'] = bt
                    P_mesh, TC_mesh = np.meshgrid(P_values, TC_values)
                    B = np.zeros_like(P_mesh)
                    for i in range(len(P_values)):
                        for j in range(len(TC_values)):
                            self.model.params['P'] = P_values[i]
                            self.model.params['TC'] = TC_values[j]
                            B[j,i], _ = self.model.solutions_fsolve()
        else:
            P_mesh, TC_mesh = np.meshgrid(P_values, TC_values)
            B = np.zeros_like(P_mesh)
            for i in range(len(P_values)):
                for j in range(len(TC_values)):
                    self.model.params['P'] = P_values[i]
                    self.model.params['TC'] = TC_values[j]
                    B[j,i], _ = self.model.solutions_fsolve()

        return P_mesh, TC_mesh, B
    
    def solve_stream(self, borderB=np.linspace(0,1,100), borderS=np.linspace(0,1,100)):
        # params_backup = self.params.copy()
        B = borderB
        s = borderS
        B_mesh, s_mesh = np.meshgrid(B, s)

        solution_fsolve = self.model.solutions_fsolve()
        # solution_fsolve = self.model.solutions_sympy()        
        equation_values = np.array([self.model.equ_wrapper((B_mesh,s_mesh)) for B_mesh, s_mesh in zip(B_mesh.flatten(), s_mesh.flatten())])
        equation_values = equation_values.T.reshape(2, B_mesh.shape[0], B_mesh.shape[1])
        
        # self.params = params_backup
        return B_mesh, s_mesh, equation_values, solution_fsolve
    
    def B_of_TC(self, TC_values, P_values, betaA=[], betaT=[]):
        B = np.zeros_like(TC_values)
        self.model.params['P'] = P_values
        for i, tc in enumerate(TC_values):
            self.model.params['TC'] = tc
            for index_a, ba in enumerate(betaA):
                for index_t, bt in enumerate(betaT):
                    self.model.params['beta_A'] = ba
                    self.model.params['beta_T'] = bt
                    B[i], _ = self.model.solutions_fsolve()
            else:
                B[i], _ = self.model.solutions_fsolve()

        return B
    
    def solve_logistic(self, B_values):
        # params_backup = self.params.copy()
        dbdb = self.model.logistic(B_values)
        return dbdb
   
    def solve_interaction_intensity(self, TC_values, P_values, is_part_of_grid=False):
        p = P_values
        sol = []
        intensity_real_list = []
        intensity_log_list = []
        
        # Always calculate sol[0] (reference) first
        self.model.params['TC'] = TC_values[0]
        sol0 = self.B_of_P(p, TC_values[0])

        for tc in TC_values[1:]:  # start from the second TC
            self.model.params['TC'] = tc
            sol_current = self.B_of_P(p, tc)
            sol.append(sol_current)

            intensity_real = sol_current - sol0
            intensity_real_list.append(intensity_real)

            intensity_log = np.log(sol_current / sol0)
            intensity_log_list.append(intensity_log)

        return intensity_real_list, intensity_log_list
    
    def calculate_equilibrium_times(self):
        _, _, _, solutionsbs = self.solve_stream(borderB=np.linspace(0, 1, 100), borderS=np.linspace(0, 1, 100))
        b_star, s_star = solutionsbs[0] / 1000, solutionsbs[1]
        a_a = self.model.params['A'] * self.model.shading_assimilation() * s_star * (1 - (2 * b_star) / self.model.params['K']) - self.model.params['D']
        b_b = self.model.params['A'] * self.model.shading_assimilation() * b_star * (1 - b_star / self.model.params['K'])
        c_c = (-1 / (self.model.params['n'] * self.model.params['Z'])) * (self.model.params['E0'] * self.model.shading_transpiration()) * s_star 
        d_d = (-1 / self.model.params['n'] * self.model.params['Z']) * (self.model.params['E0'] * self.model.shading_transpiration() * b_star + self.model.params['Ksat'] * self.model.params['c'] * s_star ** (self.model.params['c'] - 1))
        matrix_j = np.array([[a_a, b_b], [c_c, d_d]])
        eigenvalues, _ = np.linalg.eig(matrix_j)
        teq_1 = -1 / eigenvalues[0]
        teq_2 = -1 / eigenvalues[1]
        return eigenvalues
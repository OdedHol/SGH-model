
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve
from scipy import optimize
import numpy as np


class Model:
    def __init__(self, A=0.05, K_1=0.2, TC=1, K=1, D=0.01,P=1.5, E0=5, n=0.4, Z=300, Ksat=800, c=10, beta_A=3, beta_T=4, lamda_min=-2, lamda_max=0, r=0.25, c_0=1, c_1=1.25, c_2=-3.5, c_3=1.25, S_h=0.3, S_sc=0.45, S_s=0.6, tree=False, Shading=1, is_logistic=1, tree_is='both', menten=False, new_S=False, linear=False, dark_ET=0,
                 y=(0.745811, 0.50321), fsolve_iterations=20, **kwargs):
        self.y = y         
        self.params = {
                  'A':A,      # Assimilation (1/d)
                  'TC':TC,    # Tree Cover (kg)
                  'K':K,      # Max Biomass (kg)
                  'D':D,      # Death (kg/m^2*d)
                  'P':P,      # Precipitation Rate (mm/day)
                  'E0':E0,    # Evaportaion Rate (mm*kg/d*m^2 )
                  'n':n,      # Porosity (-)
                  'Z':Z,      # Root depth (mm) 
                  'Ksat':Ksat, # Saturated Hydraulic Conductivity (mm/day)
                  'c':c,      # exponent of leakage
                  'beta_A':beta_A, # Assimilation exp 
                  'beta_T':beta_T, # Transpiration exp
                  'lamda_min':lamda_min, # Tree water constant
                  'lamda_max':lamda_max, # Tree water constant
                  'r': r, # Transpiration constant
                  'c_0': c_0, # Assimilation constant
                  'c_1': c_1, # Assimilation constant
                  'c_2': c_2, # Assimilation constant
                  'c_3': c_3, # Assimilation constant
                  'S_h': S_h, # Hygroscopic point
                  'S_sc': S_sc, # Field capacity
                  'S_s': S_s, # Saturation point
                  'tree': tree, # Tree or no tree
                  'Shading': Shading, # Shading or no shading
                  'is_logistic': is_logistic, # Carrying capacity or no carrying capacity
                  'tree_is': tree_is, # Tree is present or not
                  'menten': menten, # Menten or no menten
                  'new_S': new_S, # New S or old S
                  'dark_ET': dark_ET, # Dark evapotranspiration
                  'linear': linear, # Linear or non-linear
                  'initial_guess':(0.9, 0.5),  # initial guess for fsolve
                  'fsolve_iterations':fsolve_iterations,
                  }
        self.params.update(kwargs)

    def backup_params(self):
        """Backs up the current model parameters."""
        self._backup_params = self.params.copy()

    def restore_params(self):
        """Restores model parameters to the backed up state."""
        if hasattr(self, "_backup_params"):
            self.params = self._backup_params.copy()
        else:
            raise AttributeError("No backup found. Please backup parameters before restoration.")

    def leakage(self):
        B, s = self.y
        return self.params['Ksat'] * s ** self.params['c']

    def transpiration(self):
        B, s = self.y                           
        return self.params['E0'] * B * s

    def shading_assimi(self):
        return (1 - self.params['TC']) ** (1/self.params['beta_A'])
    
    def shading_transp(self):
        return (1 - self.params['TC']) ** self.params['beta_T']
    
    def shading_assimilation(self):
        if self.params['linear']:
            return (1-self.params['TC'])
            # return (1 - self.params['TC']) ** (self.params['beta_A'])
        else:
            return (1 - self.params['Shading']*self.params['TC']) ** (1/float(self.params['beta_A']))

        # return (self.params['c_0'] + self.params['c_1'] * self.params['TC'] + self.params['c_2'] * self.params['TC'] ** 2 + self.params['c_3'] * self.params['TC'] ** 3)

    def shading_transpiration(self):
        if self.params['linear']:
            return (1-self.params['TC'])
            # return (1 - self.params['Shading']*self.params['TC']) ** (1/self.params['beta_T'])
        elif self.params['dark_ET']==1:
            e_dark = 0#0.3
            return self.params['Shading']* ((1 - e_dark) * (1 - self.params['TC'])**self.params['beta_T']  + e_dark)
        else:
            return (1 - self.params['Shading']*self.params['TC']) ** float(self.params['beta_T'])
               
        # return ((1-self.params['r'])*(1-self.params['TC'])**self.params['beta_T'] + self.params['r'])

    def beta_S(self):
        B, s = self.y
        if self.params['new_S']:
            if self.params['menten']:
                if s <= self.params['S_h']:
                    return 0
                else:
                    return (s - self.params['S_h']) / (0.1 + (s - self.params['S_h']))
            else:
                return np.piecewise(s, 
                            [s < self.params['S_h'], (self.params['S_h'] <= s) & (s <= self.params['S_sc']), s > self.params['S_sc']], 
                            [0, lambda s: (s - self.params['S_h']) / (self.params['S_sc'] - self.params['S_h']), 1])
        else:
            return s

    def tree_water(self):
        B, s = self.y
        ds = self.params['S_s'] - self.params['S_h']
        dlamda = self.params['lamda_max'] - self.params['lamda_min']
        if self.params['tree']:
            if self.params['tree_is'] == 'pro':
                self.params['lamda_min'] = -2
                self.params['lamda_max'] = 0
            elif self.params['tree_is'] == 'use':
                self.params['lamda_min'] = 0
                self.params['lamda_max'] = 10
            else:
                self.params['lamda_min'] = -2
                self.params['lamda_max'] = 10
            if s < self.params['S_h']:
                return (self.params['TC'] * dlamda * self.params['S_h']) / ds + self.params['TC'] * self.params['lamda_min'] - (self.params['TC'] * dlamda * self.params['S_h']) / ds
                # return 0
            if s > self.params['S_s']:
                return (self.params['TC'] * dlamda * self.params['S_s']) / ds + self.params['TC'] * self.params['lamda_min'] - (self.params['TC'] * dlamda * self.params['S_h']) / ds
            else:
                return (self.params['TC'] * dlamda * s) / ds + self.params['TC'] * self.params['lamda_min'] - (self.params['TC'] * dlamda * self.params['S_h']) / ds
        else:
            return 0

    def carrying_capacity(self):
        B, s = self.y
        return (1 - self.params['is_logistic']*B / self.params['K'])
           
    def growth_biomass(self, is_tree=False):
        B, s = self.y
        if self.params['new_S']:
            return self.params['A'] * self.shading_assimilation() * self.beta_S() * B * self.carrying_capacity()    
        else:
            return self.params['A'] * self.shading_assimilation() * s * B * self.carrying_capacity()
        
    def death_biomass(self):
        B, s = self.y 
        return B * self.params['D']

    def evapotranspiration(self, is_tree=False):
        B ,s = self.y 
        if self.params['new_S']:
            return self.params['E0'] * self.shading_transpiration() * B * self.beta_S()
        else:
            return self.params['E0'] * self.shading_transpiration() * B * s
        # return 0

    def logistic(self, B):
            return B * (1 - (B / self.params['K']))
       
    def equ(self):
        dydt = np.array([
            self.growth_biomass() - self.death_biomass(),
            # (self.params['P'] - self.evapotranspiration() - self.leakage()) / (self.params['n'] * self.params['Z'])
            (self.params['P'] - self.tree_water() - self.evapotranspiration() - self.leakage()) / (self.params['n'] * self.params['Z'])
        ])
        return dydt
    
    def equ_wrapper(self, y):
        self.y = y
        return self.equ()

    def equ_wrapper_ivp(self, t, y):
        self.y = y
        return self.equ()

    def solutions_fsolve(self):
        init_guess = self.params['initial_guess']
        t_span = (0, 50000)
        t_eval = np.linspace(t_span[0], t_span[1], 300)
        for i in range(self.params['fsolve_iterations']):
            sol_fsolve = fsolve(self.equ_wrapper, init_guess)
            B = sol_fsolve[0] * 1000
            s = sol_fsolve[1]
            init_guess = (np.random.uniform(0.5, 1.5), np.random.uniform(0.5, 0.6))
            if (B > 10) and (0 < s < 1):
                return [B, s]
            else:
                init_guess = (np.random.uniform(0.5, 1.5), np.random.uniform(0.5, 0.6))
        # Integrate time integration, using solve ode; return the last value
        sol = solve_ivp(self.equ_wrapper_ivp, t_span, init_guess, t_eval=t_eval) # Implement a threshold for the change of B that gets to equilbrium
        B = sol.y[0, -1] * 1000
        s = sol.y[1, -1]
        return [B,s]

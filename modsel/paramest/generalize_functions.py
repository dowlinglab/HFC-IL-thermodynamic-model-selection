'''
Generalized Functions

Setup access to all models, all parameters, all data for all systems

Created by Bridgette Befort, Jialu Wang, Alejandro Garciadiego, Alexander Dowling
University of Notre Dame (2022)
'''

# Import objects from pyomo package
from pyomo.environ import (ConcreteModel,
                           SolverFactory,
                           units as pyunits)

import idaes.logger as idaeslog
import pyomo.contrib.parmest.parmest as parmest
import pandas as pd
import pytest

from idaes.core import FlowsheetBlock
from idaes.generic_models.unit_models import Flash
from idaes.generic_models.properties.core.generic.generic_property import (
        GenericParameterBlock)

from idaes.core.util.model_statistics import degrees_of_freedom

# Setup class for PR models
class PRModels:
    def __init__(self, theta, configuration, comp_1, comp_2, x_comp_1, x_comp_2, mode='constant'):
        '''
        To run a PR model, you need:
        
        theta: parameters dataframe (to initialize or already fit)
        configuration: file with information on the system and model
        comp_1: component 1
        comp_2: component 2
        x_comp_1: name of component 1 mole fraction column in csv file
        x_comp_2: name of component 1 mole fraction column in csv file
        mode (str): type of PR model that is being used
        '''
        
        self.configuration = configuration
        
        self.comp_1 = comp_1
        self.comp_2 = comp_2
        self.x_comp_1 = x_comp_1
        self.x_comp_2 = x_comp_2
        
#         if mode == '1param_Opt1':
#             self.__1param_Opt1_parse_theta()
            
#         if mode == '1param_Opt2':
#             self.__1param_Opt2_parse_theta()
        
#         if mode == 'constant':
#             self.__constant_parse_theta()
            
#         if mode == '3params_Opt1':
#             self.__3params_Opt1_parse_theta()
            
#         if mode == '3params_Opt2':
#             self.__3params_Opt2_parse_theta()
            
#         if mode == 'linear':
#             self.__linear_parse_theta()
            
#         if mode == 'quadratic':
#             self.__quadratic_parse_theta()
            
#         if mode == 'polynomial':
#             self.__polynomial_parse_theta()
            
            
    def __parse_theta(self): #need to change this in access files
        self.PR_kappa_A_comp_1_comp_2 = theta[0]
        self.PR_kappa_A_comp_2_comp_1 = theta[1]
        self.PR_kappa_B_comp_1_comp_2 = theta[2]
        self.PR_kappa_B_comp_2_comp_1 = theta[3]
        self.PR_kappa_C_comp_1_comp_2 = theta[4]
        self.PR_kappa_C_comp_2_comp_1 = theta[5]
        self.PR_kappa_D_comp_1_comp_2 = theta[6]
        self.PR_kappa_D_comp_2_comp_1 = theta[7]
        
#     def __1param_Opt2_parse_theta(self):
#         self.PR_kappa_A_comp_1_comp_2 = 0
#         self.PR_kappa_A_comp_2_comp_1 = theta[0]
#         self.PR_kappa_B_comp_1_comp_2 = 0
#         self.PR_kappa_B_comp_2_comp_1 = 0
#         self.PR_kappa_C_comp_1_comp_2 = 0
#         self.PR_kappa_C_comp_2_comp_1 = 0
#         self.PR_kappa_D_comp_1_comp_2 = 0
#         self.PR_kappa_D_comp_2_comp_1 = 0
    
#     def __constant_parse_theta(self):
#         self.PR_kappa_A_comp_1_comp_2 = theta[]
#         self.PR_kappa_A_comp_2_comp_1 = theta[]
#         self.PR_kappa_B_comp_1_comp_2 = 0
#         self.PR_kappa_B_comp_2_comp_1 = 0
#         self.PR_kappa_C_comp_1_comp_2 = 0
#         self.PR_kappa_C_comp_2_comp_1 = 0
#         self.PR_kappa_D_comp_1_comp_2 = 0
#         self.PR_kappa_D_comp_2_comp_1 = 0
        
#     def __3params_Opt1_parse_theta(self):
#         self.PR_kappa_A_comp_1_comp_2 = theta[]
#         self.PR_kappa_A_comp_2_comp_1 = theta[]
#         self.PR_kappa_B_comp_1_comp_2 = theta[]
#         self.PR_kappa_B_comp_2_comp_1 = 0
#         self.PR_kappa_C_comp_1_comp_2 = 0
#         self.PR_kappa_C_comp_2_comp_1 = 0
#         self.PR_kappa_D_comp_1_comp_2 = 0
#         self.PR_kappa_D_comp_2_comp_1 = 0
        
#     def __3params_Opt2_parse_theta(self):
#         self.PR_kappa_A_comp_1_comp_2 = theta[]
#         self.PR_kappa_A_comp_2_comp_1 = theta[]
#         self.PR_kappa_B_comp_1_comp_2 = 0
#         self.PR_kappa_B_comp_2_comp_1 = theta[]
#         self.PR_kappa_C_comp_1_comp_2 = 0
#         self.PR_kappa_C_comp_2_comp_1 = 0
#         self.PR_kappa_D_comp_1_comp_2 = 0
#         self.PR_kappa_D_comp_2_comp_1 = 0
        
#     def __linear_parse_theta(self):
#         self.PR_kappa_A_comp_1_comp_2 = theta[]
#         self.PR_kappa_A_comp_2_comp_1 = theta[]
#         self.PR_kappa_B_comp_1_comp_2 = theta[]
#         self.PR_kappa_B_comp_2_comp_1 = theta[]
#         self.PR_kappa_C_comp_1_comp_2 = 0
#         self.PR_kappa_C_comp_2_comp_1 = 0
#         self.PR_kappa_D_comp_1_comp_2 = 0
#         self.PR_kappa_D_comp_2_comp_1 = 0
        
#     def __quadratic_parse_theta(self):
#         self.PR_kappa_A_comp_1_comp_2 = theta[]
#         self.PR_kappa_A_comp_2_comp_1 = theta[]
#         self.PR_kappa_B_comp_1_comp_2 = theta[]
#         self.PR_kappa_B_comp_2_comp_1 = theta[]
#         self.PR_kappa_C_comp_1_comp_2 = theta[]
#         self.PR_kappa_C_comp_2_comp_1 = theta[]
#         self.PR_kappa_D_comp_1_comp_2 = 0
#         self.PR_kappa_D_comp_2_comp_1 = 0
        
#     def __polynomial_parse_theta(self):
#         self.PR_kappa_A_comp_1_comp_2 = theta[]
#         self.PR_kappa_A_comp_2_comp_1 = theta[]
#         self.PR_kappa_B_comp_1_comp_2 = theta[]
#         self.PR_kappa_B_comp_2_comp_1 = theta[]
#         self.PR_kappa_C_comp_1_comp_2 = theta[]
#         self.PR_kappa_C_comp_2_comp_1 = theta[]
#         self.PR_kappa_D_comp_1_comp_2 = theta[]
#         self.PR_kappa_D_comp_2_comp_1 = theta[]
        
    def create_model(data):
        '''
        
        '''
        
        m = ConcreteModel()

        m.fs = FlowsheetBlock(default={"dynamic": False})

        m.fs.properties = GenericParameterBlock(default=configuration)

        m.fs.state_block = m.fs.properties.state_block_class(
            default={"parameters": m.fs.properties,
                     "defined_state": True})
        x = float(data[x_comp_1])+eps
        m.fs.state_block.flow_mol.fix(1)
        m.fs.state_block.temperature.fix(float(data["temperature"]))
        m.fs.state_block.pressure.fix(float(data["pressure"]))
        m.fs.state_block.mole_frac_comp[comp_2].fix(1-x)
        m.fs.state_block.mole_frac_comp[comp_1].fix(x)

        # parameters
        m.fs.properties.PR_kappa_A[comp_2, comp_2].fix(0)
        m.fs.properties.PR_kappa_A[comp_2, comp_1].fix(PR_kappa_A_comp_2_comp_1)
        m.fs.properties.PR_kappa_A[comp_1, comp_1].fix(0)
        m.fs.properties.PR_kappa_A[comp_1, comp_2].fix(PR_kappa_A_comp_1_comp_2)
        m.fs.properties.PR_kappa_B[comp_2, comp_2].fix(0)
        m.fs.properties.PR_kappa_B[comp_2, comp_1].fix(PR_kappa_B_comp_2_comp_1)
        m.fs.properties.PR_kappa_B[comp_1, comp_1].fix(0)
        m.fs.properties.PR_kappa_B[comp_1, comp_2].fix(PR_kappa_B_comp_1_comp_2)
        m.fs.properties.PR_kappa_C[comp_2, comp_2].fix(0)
        m.fs.properties.PR_kappa_C[comp_2, comp_1].fix(PR_kappa_C_comp_2_comp_1)
        m.fs.properties.PR_kappa_C[comp_1, comp_1].fix(0)
        m.fs.properties.PR_kappa_C[comp_1, comp_2].fix(PR_kappa_C_comp_1_comp_2)
        m.fs.properties.PR_kappa_D[comp_2, comp_2].fix(0)
        m.fs.properties.PR_kappa_D[comp_2, comp_1].fix(PR_kappa_D_comp_2_comp_1)
        m.fs.properties.PR_kappa_D[comp_1, comp_1].fix(0)
        m.fs.properties.PR_kappa_D[comp_1, comp_2].fix(PR_kappa_D_comp_1_comp_2)

        # Initialize the flash unit
#         m.fs.state_block.initialize(outlvl=idaeslog.CRITICAL)
        m.fs.state_block.initialize(outlvl=50)

        # Fix the state variables on the state block
        m.fs.state_block.pressure.unfix()
        m.fs.state_block.mole_frac_comp[comp_2].unfix()
        m.fs.state_block.mole_frac_comp[comp_1].unfix()
        m.fs.state_block.temperature.fix(float(data["temperature"]))
        m.fs.state_block.mole_frac_phase_comp['Liq', comp_1].fix(float(data[x_comp_1]))
        m.fs.state_block.mole_frac_phase_comp['Liq', comp_2].fix(float(data[x_comp_2]))
        m.fs.state_block.mole_frac_comp[comp_1].fix(float(data[x_comp_1])+eps)
        m.fs.state_block.mole_frac_comp[comp_2].unfix()


        # Set bounds on variables to be estimated
        m.fs.properties.PR_kappa_A[comp_2, comp_1].setlb(-20)
        m.fs.properties.PR_kappa_A[comp_2, comp_1].setub(20)

        m.fs.properties.PR_kappa_A[comp_1, comp_2].setlb(-20)
        m.fs.properties.PR_kappa_A[comp_1, comp_2].setub(20)

        m.fs.properties.PR_kappa_B[comp_2, comp_1].setlb(-20)
        m.fs.properties.PR_kappa_B[comp_2, comp_1].setub(20)

        m.fs.properties.PR_kappa_B[comp_1, comp_2].setlb(-20)
        m.fs.properties.PR_kappa_B[comp_1, comp_2].setub(20)

        m.fs.properties.PR_kappa_C[comp_2, comp_1].setlb(-20)
        m.fs.properties.PR_kappa_C[comp_2, comp_1].setub(20)

        m.fs.properties.PR_kappa_C[comp_1, comp_2].setlb(-20)
        m.fs.properties.PR_kappa_C[comp_1, comp_2].setub(20)
        
        m.fs.properties.PR_kappa_D[comp_2, comp_1].setlb(-20)
        m.fs.properties.PR_kappa_D[comp_2, comp_1].setub(20)

        m.fs.properties.PR_kappa_D[comp_1, comp_2].setlb(-20)
        m.fs.properties.PR_kappa_D[comp_1, comp_2].setub(20)

        # Return initialized flash model
        return m
        
        
'''
To Add:

-estimate parameters function
-FIM function
'''

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

# Setup class for PR, SRK models
class CEOSModels:
    def __init__(self, theta, configuration, comp_1, comp_2, x_comp_1, x_comp_2, 
                  Model_type = 'PR', Tdep_type='1Param_Opt1'):
        '''
        To run a CEOS model, you need:
        
        
        theta: either file or parameters dataframe (to initialize or already fit)
        Model_type: a string of PR or SRK model type
        Tdep_type: a string of Temperature dependence/model type
        configuration: file with information on the system and model
        comp_1: component 1
        comp_2: component 2
        x_comp_1: name of component 1 mole fraction column in csv file
        x_comp_2: name of component 1 mole fraction column in csv file
        init_temp: temperature to initialize the model 
        init_pressure: pressure to initialize the model
        init_x_c1: mole fraction to initialize the model
        '''
        
        self.configuration = configuration
        
        self.comp_1 = comp_1
        self.comp_2 = comp_2
        self.x_comp_1 = x_comp_1
        self.x_comp_2 = x_comp_2
        
        #self.init_temp = init_temp
        #self.init_pressure = init_pressure
        #self.init_x_c1 = init_x_c1
        
        self.Model_type = Model_type
        
        self.Tdep_type = Tdep_type
        
        # check file or vector here.
        if type(theta) is str: 
            self.__parse_theta_csv(theta)
        else:
            self.__parse_theta(theta)
        
        self.__get_param_names()
        
    def __parse_theta_csv(self, file_name):
        
        params = pd.read_csv(file_name,header=None)

        # all parameters for initialization
        params_list = list(params[0])
        
        self.__parse_theta(params_list)
        
        
    def __parse_theta(self, theta): 
        
        if self.Model_type == "PR":
            self.PR_kappa_A_comp_1_comp_2 = theta[0]
            self.PR_kappa_A_comp_2_comp_1 = theta[1]
            self.PR_kappa_B_comp_1_comp_2 = theta[2]
            self.PR_kappa_B_comp_2_comp_1 = theta[3]
            self.PR_kappa_C_comp_1_comp_2 = theta[4]
            self.PR_kappa_C_comp_2_comp_1 = theta[5]
            self.PR_kappa_D_comp_1_comp_2 = theta[6]
            self.PR_kappa_D_comp_2_comp_1 = theta[7]
            
        elif self.Model_type == "SRK":
            self.SRK_kappa_A_comp_1_comp_2 = theta[0]
            self.SRK_kappa_A_comp_2_comp_1 = theta[1]
            self.SRK_kappa_B_comp_1_comp_2 = theta[2]
            self.SRK_kappa_B_comp_2_comp_1 = theta[3]
            self.SRK_kappa_C_comp_1_comp_2 = theta[4]
            self.SRK_kappa_C_comp_2_comp_1 = theta[5]
            self.SRK_kappa_D_comp_1_comp_2 = theta[6]
            self.SRK_kappa_D_comp_2_comp_1 = theta[7]
            
        else:
            print("Wrong model type.")
        
        
    def __get_param_names(self):
        if self.Model_type == 'PR':
        
            if self.Tdep_type == '1Param_Opt1':
                param_name_dict = { 'fs.properties.PR_kappa_A["R32", "emimTf2N"]': self.PR_kappa_A_comp_1_comp_2}

            elif self.Tdep_type == '1Param_Opt2':
                param_name_dict = { 'fs.properties.PR_kappa_A["emimTf2N", "R32"]': self.PR_kappa_A_comp_2_comp_1}

            elif self.Tdep_type == 'No':
                param_name_dict = { 'fs.properties.PR_kappa_A["R32", "emimTf2N"]': self.PR_kappa_A_comp_1_comp_2, 
                                    'fs.properties.PR_kappa_A["emimTf2N", "R32"]': self.PR_kappa_A_comp_2_comp_1}

            elif self.Tdep_type == '3Params_Opt1':
                param_name_dict = { 'fs.properties.PR_kappa_A["R32", "emimTf2N"]': self.PR_kappa_A_comp_1_comp_2, 
                                    'fs.properties.PR_kappa_A["emimTf2N", "R32"]': self.PR_kappa_A_comp_2_comp_1,
                                    'fs.properties.PR_kappa_B["R32", "emimTf2N"]': self.PR_kappa_B_comp_1_comp_2}

            elif self.Tdep_type == '3Params_Opt2':
                param_name_dict = { 'fs.properties.PR_kappa_A["R32", "emimTf2N"]': self.PR_kappa_A_comp_1_comp_2, 
                                    'fs.properties.PR_kappa_A["emimTf2N", "R32"]': self.PR_kappa_A_comp_2_comp_1,
                                    'fs.properties.PR_kappa_B["emimTf2N", "R32"]': self.PR_kappa_B_comp_2_comp_1}

            elif self.Tdep_type == 'Linear':
                param_name_dict = { 'fs.properties.PR_kappa_A["R32", "emimTf2N"]': self.PR_kappa_A_comp_1_comp_2, 
                                    'fs.properties.PR_kappa_A["emimTf2N", "R32"]': self.PR_kappa_A_comp_2_comp_1,
                                    'fs.properties.PR_kappa_B["R32", "emimTf2N"]': self.PR_kappa_B_comp_1_comp_2,
                                    'fs.properties.PR_kappa_B["emimTf2N", "R32"]': self.PR_kappa_B_comp_2_comp_1}

            elif self.Tdep_type == 'Quadratic':
                param_name_dict = { 'fs.properties.PR_kappa_A["R32", "emimTf2N"]': self.PR_kappa_A_comp_1_comp_2, 
                                    'fs.properties.PR_kappa_A["emimTf2N", "R32"]': self.PR_kappa_A_comp_2_comp_1,
                                    'fs.properties.PR_kappa_B["R32", "emimTf2N"]': self.PR_kappa_B_comp_1_comp_2,
                                    'fs.properties.PR_kappa_B["emimTf2N", "R32"]': self.PR_kappa_B_comp_2_comp_1,
                                    'fs.properties.PR_kappa_C["R32", "emimTf2N"]': self.PR_kappa_C_comp_1_comp_2,
                                    'fs.properties.PR_kappa_C["emimTf2N", "R32"]': self.PR_kappa_C_comp_2_comp_1}

            elif self.Tdep_type == 'Polynomial':
                param_name_dict = { 'fs.properties.PR_kappa_A["R32", "emimTf2N"]': self.PR_kappa_A_comp_1_comp_2, 
                                    'fs.properties.PR_kappa_A["emimTf2N", "R32"]': self.PR_kappa_A_comp_2_comp_1,
                                    'fs.properties.PR_kappa_B["R32", "emimTf2N"]': self.PR_kappa_B_comp_1_comp_2,
                                    'fs.properties.PR_kappa_B["emimTf2N", "R32"]': self.PR_kappa_B_comp_2_comp_1,
                                    'fs.properties.PR_kappa_C["R32", "emimTf2N"]': self.PR_kappa_C_comp_1_comp_2,
                                    'fs.properties.PR_kappa_C["emimTf2N", "R32"]': self.PR_kappa_C_comp_2_comp_1,
                                    'fs.properties.PR_kappa_D["R32", "emimTf2N"]': self.PR_kappa_D_comp_1_comp_2,
                                    'fs.properties.PR_kappa_D["emimTf2N", "R32"]': self.PR_kappa_D_comp_2_comp_1}
                
        elif self.Model_type == 'SRK':
            
            if self.Tdep_type == '1Param_Opt1':
                param_name_dict = { 'fs.properties.SRK_kappa_A["R32", "emimTf2N"]': self.SRK_kappa_A_comp_1_comp_2}

            elif self.Tdep_type == '1Param_Opt2':
                param_name_dict = { 'fs.properties.SRK_kappa_A["emimTf2N", "R32"]': self.SRK_kappa_A_comp_2_comp_1}

            elif self.Tdep_type == 'No':
                param_name_dict = { 'fs.properties.SRK_kappa_A["R32", "emimTf2N"]': self.SRK_kappa_A_comp_1_comp_2, 
                                    'fs.properties.SRK_kappa_A["emimTf2N", "R32"]': self.SRK_kappa_A_comp_2_comp_1}

            elif self.Tdep_type == '3Params_Opt1':
                param_name_dict = { 'fs.properties.SRK_kappa_A["R32", "emimTf2N"]': self.SRK_kappa_A_comp_1_comp_2, 
                                    'fs.properties.SRK_kappa_A["emimTf2N", "R32"]': self.SRK_kappa_A_comp_2_comp_1,
                                    'fs.properties.SRK_kappa_B["R32", "emimTf2N"]': self.SRK_kappa_B_comp_1_comp_2}

            elif self.Tdep_type == '3Params_Opt2':
                param_name_dict = { 'fs.properties.SRK_kappa_A["R32", "emimTf2N"]': self.SRK_kappa_A_comp_1_comp_2, 
                                    'fs.properties.SRK_kappa_A["emimTf2N", "R32"]': self.SRK_kappa_A_comp_2_comp_1,
                                    'fs.properties.SRK_kappa_B["emimTf2N", "R32"]': self.SRK_kappa_B_comp_2_comp_1}

            elif self.Tdep_type == 'Linear':
                param_name_dict = { 'fs.properties.SRK_kappa_A["R32", "emimTf2N"]': self.SRK_kappa_A_comp_1_comp_2, 
                                    'fs.properties.SRK_kappa_A["emimTf2N", "R32"]': self.SRK_kappa_A_comp_2_comp_1,
                                    'fs.properties.SRK_kappa_B["R32", "emimTf2N"]': self.SRK_kappa_B_comp_1_comp_2,
                                    'fs.properties.SRK_kappa_B["emimTf2N", "R32"]': self.SRK_kappa_B_comp_2_comp_1}

            elif self.Tdep_type == 'Quadratic':
                param_name_dict = { 'fs.properties.SRK_kappa_A["R32", "emimTf2N"]': self.SRK_kappa_A_comp_1_comp_2, 
                                    'fs.properties.SRK_kappa_A["emimTf2N", "R32"]': self.SRK_kappa_A_comp_2_comp_1,
                                    'fs.properties.SRK_kappa_B["R32", "emimTf2N"]': self.SRK_kappa_B_comp_1_comp_2,
                                    'fs.properties.SRK_kappa_B["emimTf2N", "R32"]': self.SRK_kappa_B_comp_2_comp_1,
                                    'fs.properties.SRK_kappa_C["R32", "emimTf2N"]': self.SRK_kappa_C_comp_1_comp_2,
                                    'fs.properties.SRK_kappa_C["emimTf2N", "R32"]': self.SRK_kappa_C_comp_2_comp_1}

            elif self.Tdep_type == 'Polynomial':
                param_name_dict = { 'fs.properties.SRK_kappa_A["R32", "emimTf2N"]': self.SRK_kappa_A_comp_1_comp_2, 
                                    'fs.properties.SRK_kappa_A["emimTf2N", "R32"]': self.SRK_kappa_A_comp_2_comp_1,
                                    'fs.properties.SRK_kappa_B["R32", "emimTf2N"]': self.SRK_kappa_B_comp_1_comp_2,
                                    'fs.properties.SRK_kappa_B["emimTf2N", "R32"]': self.SRK_kappa_B_comp_2_comp_1,
                                    'fs.properties.SRK_kappa_C["R32", "emimTf2N"]': self.SRK_kappa_C_comp_1_comp_2,
                                    'fs.properties.SRK_kappa_C["emimTf2N", "R32"]': self.SRK_kappa_C_comp_2_comp_1,
                                    'fs.properties.SRK_kappa_D["R32", "emimTf2N"]': self.SRK_kappa_D_comp_1_comp_2,
                                    'fs.properties.SRK_kappa_D["emimTf2N", "R32"]': self.SRK_kappa_D_comp_2_comp_1}
                
            
        self.param_name_dict = param_name_dict

        
    def create_model(self, data, eps=0.0, 
                     init_temp = 283.1, init_pressure = 399300, init_x_c1 = 0.45, polynomial = False):
        '''
        
        '''
        
        m = ConcreteModel()

        m.fs = FlowsheetBlock(default={"dynamic": False})

        m.fs.properties = GenericParameterBlock(default=self.configuration)

        m.fs.state_block = m.fs.properties.state_block_class(
            default={"parameters": m.fs.properties,
                     "defined_state": True})
        
        #if not init:
        #    init_x_c1 = self.init_x_c1
            
        #    init_temp = self.init_temp
            
        #    init_pressure = self.init_pressure
        
        x = float(init_x_c1)+eps
        m.fs.state_block.flow_mol.fix(1)
        m.fs.state_block.temperature.fix(init_temp)
        m.fs.state_block.pressure.fix(init_pressure)
        m.fs.state_block.mole_frac_comp[self.comp_2].fix(1-x)
        m.fs.state_block.mole_frac_comp[self.comp_1].fix(x)
        
        if self.Model_type == 'PR':

            if polynomial == False:
                # parameters
                m.fs.properties.PR_kappa_A[self.comp_2, self.comp_2].fix(0)
                m.fs.properties.PR_kappa_A[self.comp_2, self.comp_1].fix(self.PR_kappa_A_comp_2_comp_1)
                m.fs.properties.PR_kappa_A[self.comp_1, self.comp_1].fix(0)
                m.fs.properties.PR_kappa_A[self.comp_1, self.comp_2].fix(self.PR_kappa_A_comp_1_comp_2)
                m.fs.properties.PR_kappa_B[self.comp_2, self.comp_2].fix(0)
                m.fs.properties.PR_kappa_B[self.comp_2, self.comp_1].fix(self.PR_kappa_B_comp_2_comp_1)
                m.fs.properties.PR_kappa_B[self.comp_1, self.comp_1].fix(0)
                m.fs.properties.PR_kappa_B[self.comp_1, self.comp_2].fix(self.PR_kappa_B_comp_1_comp_2)
                m.fs.properties.PR_kappa_C[self.comp_2, self.comp_2].fix(0)
                m.fs.properties.PR_kappa_C[self.comp_2, self.comp_1].fix(self.PR_kappa_C_comp_2_comp_1)
                m.fs.properties.PR_kappa_C[self.comp_1, self.comp_1].fix(0)
                m.fs.properties.PR_kappa_C[self.comp_1, self.comp_2].fix(self.PR_kappa_C_comp_1_comp_2)

            else:
                # parameters
                m.fs.properties.PR_kappa_A[self.comp_2, self.comp_2].fix(0)
                m.fs.properties.PR_kappa_A[self.comp_2, self.comp_1].fix(self.PR_kappa_A_comp_2_comp_1)
                m.fs.properties.PR_kappa_A[self.comp_1, self.comp_1].fix(0)
                m.fs.properties.PR_kappa_A[self.comp_1, self.comp_2].fix(self.PR_kappa_A_comp_1_comp_2)
                m.fs.properties.PR_kappa_B[self.comp_2, self.comp_2].fix(0)
                m.fs.properties.PR_kappa_B[self.comp_2, self.comp_1].fix(self.PR_kappa_B_comp_2_comp_1)
                m.fs.properties.PR_kappa_B[self.comp_1, self.comp_1].fix(0)
                m.fs.properties.PR_kappa_B[self.comp_1, self.comp_2].fix(self.PR_kappa_B_comp_1_comp_2)
                m.fs.properties.PR_kappa_C[self.comp_2, self.comp_2].fix(0)
                m.fs.properties.PR_kappa_C[self.comp_2, self.comp_1].fix(self.PR_kappa_C_comp_2_comp_1)
                m.fs.properties.PR_kappa_C[self.comp_1, self.comp_1].fix(0)
                m.fs.properties.PR_kappa_C[self.comp_1, self.comp_2].fix(self.PR_kappa_C_comp_1_comp_2)
                m.fs.properties.PR_kappa_D[self.comp_2, self.comp_2].fix(0)
                m.fs.properties.PR_kappa_D[self.comp_2, self.comp_1].fix(self.PR_kappa_D_comp_2_comp_1)
                m.fs.properties.PR_kappa_D[self.comp_1, self.comp_1].fix(0)
                m.fs.properties.PR_kappa_D[self.comp_1, self.comp_2].fix(self.PR_kappa_D_comp_1_comp_2)
                
        elif self.Model_type == 'SRK':

            if polynomial == False:
                # parameters
                m.fs.properties.SRK_kappa_A[self.comp_2, self.comp_2].fix(0)
                m.fs.properties.SRK_kappa_A[self.comp_2, self.comp_1].fix(self.SRK_kappa_A_comp_2_comp_1)
                m.fs.properties.SRK_kappa_A[self.comp_1, self.comp_1].fix(0)
                m.fs.properties.SRK_kappa_A[self.comp_1, self.comp_2].fix(self.SRK_kappa_A_comp_1_comp_2)
                m.fs.properties.SRK_kappa_B[self.comp_2, self.comp_2].fix(0)
                m.fs.properties.SRK_kappa_B[self.comp_2, self.comp_1].fix(self.SRK_kappa_B_comp_2_comp_1)
                m.fs.properties.SRK_kappa_B[self.comp_1, self.comp_1].fix(0)
                m.fs.properties.SRK_kappa_B[self.comp_1, self.comp_2].fix(self.SRK_kappa_B_comp_1_comp_2)
                m.fs.properties.SRK_kappa_C[self.comp_2, self.comp_2].fix(0)
                m.fs.properties.SRK_kappa_C[self.comp_2, self.comp_1].fix(self.SRK_kappa_C_comp_2_comp_1)
                m.fs.properties.SRK_kappa_C[self.comp_1, self.comp_1].fix(0)
                m.fs.properties.SRK_kappa_C[self.comp_1, self.comp_2].fix(self.SRK_kappa_C_comp_1_comp_2)

            else:
                # parameters
                m.fs.properties.SRK_kappa_A[self.comp_2, self.comp_2].fix(0)
                m.fs.properties.SRK_kappa_A[self.comp_2, self.comp_1].fix(self.SRK_kappa_A_comp_2_comp_1)
                m.fs.properties.SRK_kappa_A[self.comp_1, self.comp_1].fix(0)
                m.fs.properties.SRK_kappa_A[self.comp_1, self.comp_2].fix(self.SRK_kappa_A_comp_1_comp_2)
                m.fs.properties.SRK_kappa_B[self.comp_2, self.comp_2].fix(0)
                m.fs.properties.SRK_kappa_B[self.comp_2, self.comp_1].fix(self.SRK_kappa_B_comp_2_comp_1)
                m.fs.properties.SRK_kappa_B[self.comp_1, self.comp_1].fix(0)
                m.fs.properties.SRK_kappa_B[self.comp_1, self.comp_2].fix(self.SRK_kappa_B_comp_1_comp_2)
                m.fs.properties.SRK_kappa_C[self.comp_2, self.comp_2].fix(0)
                m.fs.properties.SRK_kappa_C[self.comp_2, self.comp_1].fix(self.SRK_kappa_C_comp_2_comp_1)
                m.fs.properties.SRK_kappa_C[self.comp_1, self.comp_1].fix(0)
                m.fs.properties.SRK_kappa_C[self.comp_1, self.comp_2].fix(self.SRK_kappa_C_comp_1_comp_2)
                m.fs.properties.SRK_kappa_D[self.comp_2, self.comp_2].fix(0)
                m.fs.properties.SRK_kappa_D[self.comp_2, self.comp_1].fix(self.SRK_kappa_D_comp_2_comp_1)
                m.fs.properties.SRK_kappa_D[self.comp_1, self.comp_1].fix(0)
                m.fs.properties.SRK_kappa_D[self.comp_1, self.comp_2].fix(self.SRK_kappa_D_comp_1_comp_2)

        # Initialize the flash unit
        #m.fs.state_block.initialize(outlvl=idaeslog.CRITICAL)
        m.fs.state_block.initialize(outlvl=50)

        # Fix the state variables on the state block
        m.fs.state_block.pressure.unfix()
        m.fs.state_block.mole_frac_comp[self.comp_2].unfix()
        m.fs.state_block.mole_frac_comp[self.comp_1].unfix()
        m.fs.state_block.temperature.fix(float(data["temperature"]))
        m.fs.state_block.mole_frac_phase_comp['Liq', self.comp_1].fix(float(data[self.x_comp_1]))
        m.fs.state_block.mole_frac_phase_comp['Liq', self.comp_2].fix(float(data[self.x_comp_2]))
        m.fs.state_block.mole_frac_comp[self.comp_1].fix(float(data[self.x_comp_1])+eps)
        m.fs.state_block.mole_frac_comp[self.comp_2].unfix()

        # Return initialized flash model
        return m
        
        
'''
To Add:

-estimate parameters function
-FIM function
'''

'''
Fitting function for estimating parameters for the two-parameter PR EOS

Bridgette Befort, Alejandro Garciadiego, Edward Maginn, Alexander Dowling
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

# Function to compute two-parameters for PR EoS

def constant(file, configuration, comp_1, comp_2, x_comp_1, x_comp_2,
    init_temp = 323.15, init_press = 399800, init_x_c1 = 0.5, init_x_c2 = 0.5,
    init_kappa_A_1_2 = -0.3, init_kappa_A_2_1 = 0.5, eps = 0.0, scaling_fac = 1e-4, filename='filename_test.csv'):
    
    """
    Estimates kappa_A parameters for Peng Robinson equation.
    Args:
        file: csv data file in Pa and K
        configuration: imported configuration dictionary
        comp_1: component 1
        comp_2: component 2
        x_comp_1: name of component 1 mole fraction column in csv file
        x_comp_2: name of component 1 mole fraction column in csv file
        init_temp = temperature to initialize model [K]
        init_press = pressure to initialize model [Pa]
        init_x_c1 = component 1 mole fraction to initialize model [mol/mol]
        init_x_c2 = component 1 mole fraction to initialize model [mol/mol]
        init_kappa_A_1_2 = initial guess for kappa_A parameter component 2-component 1
        init_kappa_A_2_1 = initial guess for kappa_A parameter component 1-component 2
        eps = provides variation/tolerance within the dataset
        scaling_fac = 1e-4 --> scales data
        filename: file name of the ipopt output text file
    Returns:
        binary interaction parameters, kappa_A, objective function value, covariance matrix
    """

    # Read in data
    data = file
    
    # Define/initialize PR model in IDAES framework

    def PR_model(data):
        
        """
        Function which returns initilized model.
        Args:
            data: pandas DataFrame with data
        Returns:
            initialized model
        """
        
        # Define model, read in configuration
        m = ConcreteModel()

        m.fs = FlowsheetBlock(default={"dynamic": False})

        m.fs.properties = GenericParameterBlock(default=configuration)

        m.fs.state_block = m.fs.properties.state_block_class(
            default={"parameters": m.fs.properties,
                     "defined_state": True})
        
        # Set up flash specifications
        x = float(data[x_comp_1])+eps
        m.fs.state_block.flow_mol.fix(1)
        m.fs.state_block.temperature.fix(float(data["temperature"]))
        m.fs.state_block.pressure.fix(float(data["pressure"]))
        m.fs.state_block.mole_frac_comp[comp_2].fix(1-x)
        m.fs.state_block.mole_frac_comp[comp_1].fix(x)

        # parameter - kappa_A_ij (set at initial value, 0 if i=j)
        m.fs.properties.PR_kappa_A[comp_2, comp_2].fix(0)
        m.fs.properties.PR_kappa_A[comp_2, comp_1].fix(init_kappa_A_2_1)
        m.fs.properties.PR_kappa_A[comp_1, comp_1].fix(0)
        m.fs.properties.PR_kappa_A[comp_1, comp_2].fix(init_kappa_A_1_2)
        m.fs.properties.PR_kappa_B[comp_2, comp_2].fix(0)
        m.fs.properties.PR_kappa_B[comp_2, comp_1].fix(0)
        m.fs.properties.PR_kappa_B[comp_1, comp_1].fix(0)
        m.fs.properties.PR_kappa_B[comp_1, comp_2].fix(0)
        m.fs.properties.PR_kappa_C[comp_2, comp_2].fix(0)
        m.fs.properties.PR_kappa_C[comp_2, comp_1].fix(0)
        m.fs.properties.PR_kappa_C[comp_1, comp_1].fix(0)
        m.fs.properties.PR_kappa_C[comp_1, comp_2].fix(0)
        
        # Initialize the flash unit
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

        # Return initialized flash model
        return m
    
    # Define objective function as squared error

    def SSE(m, data):
        """
        returns objective function expression.
        Args:
            m: model
            data: pandas DataFrame with data
        Returns:
            objective function scaled expresion
        """
        expr = ((float(data["pressure"]) - m.fs.state_block.pressure)**2)

        return expr * scaling_fac

    # Define fitting parameters    
    variable_name = ["fs.properties.PR_kappa_A['" + comp_2 + "','" + comp_1 + "']",
                 "fs.properties.PR_kappa_A['" + comp_1 + "','" + comp_2 + "']"]
    
    # Define additional solver options
    solver_opts_dict = {'output_file':filename}
    
    #Estimate Parameters via Parmest
    pest = parmest.Estimator(PR_model, data, variable_name, SSE, tee=False, solver_options=solver_opts_dict)

    obj_value, parameters, a = pest.theta_est(calc_cov=True)
    
    return parameters, obj_value, a

class Data_reg:

    def __init__(self, configuration, comp_1, comp_2, x_comp_1, x_comp_2):
        self.configuration = configuration
        self.comp_1 = comp_1
        self.comp_2 = comp_2
        self.x_comp_1 = x_comp_1
        self.x_comp_2 = x_comp_2
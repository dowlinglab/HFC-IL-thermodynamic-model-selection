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
from idaes.generic_models.properties.activity_coeff_models.activity_coeff_prop_pack import (
    ActivityCoeffParameterData,
)
from idaes.core.util.model_statistics import degrees_of_freedom

def taus(file, configuration, comp_1, comp_2, x_comp_1, x_comp_2, init_tau_1_2 = -0.3, init_tau_2_1 = 0.5, eps = 0.0, scaling_fac = 1e-4,filename='filename_test.csv'):

    """
    Estimates tau parameters for NRTL
    Args:
        file: csv data file in Pa and K
        configuration: imported configuration dictionary
        comp_1: component 1
        comp_2: component 2
        x_comp_1: name of component 1 mole fraction column in csv file
        x_comp_2: name of component 1 mole fraction column in csv file
        init_tau_1_2 = initial guess for tau parameter
        init_tau_2_1 = initial guess for tau parameter 
        eps = extra
        scaling_fac = 1e-4)
    Returns:
        printed parameters for binary interaction parameters
    """
    
    data = file
    
    def NRTL_model(data):

        #Todo: Create a ConcreteModel object
        m = ConcreteModel()

        #Todo: Create FlowsheetBlock object
        m.fs = FlowsheetBlock(default={"dynamic": False})


        #Todo: Create a properties parameter object with the following options:
        # "valid_phase": ('Liq', 'Vap')
        # "activity_coeff_model": 'NRTL'
        m.fs.properties = configuration(default={"valid_phase":
                                                     ('Liq', 'Vap'),
                                                     "activity_coeff_model":
                                                     'NRTL'})

        m.fs.state_block = m.fs.properties.state_block_class(
            default={"parameters": m.fs.properties,
                     "defined_state": True})


        # Fix the state variables on the state block
        # hint: state variables exist on the state block i.e. on m.fs.state_block

        x = float(data[x_comp_1])+eps

        m.fs.state_block.flow_mol.fix(1)
        m.fs.state_block.temperature.fix(float(data["temperature"]))
        m.fs.state_block.pressure.fix(float(data["pressure"]))
        m.fs.state_block.mole_frac_comp[comp_2].fix(1-x)
        m.fs.state_block.mole_frac_comp[comp_1].fix(x)

        # Fix NRTL specific parameters. 

        # non-randomness parameter - alpha_ij (set at 0.3, 0 if i=j)
        m.fs.properties.\
            alpha[comp_2, comp_2].fix(0)
        m.fs.properties.\
            alpha[comp_2, comp_1].fix(0.3)
        m.fs.properties.\
            alpha[comp_1, comp_1].fix(0)
        m.fs.properties.\
            alpha[comp_1, comp_2].fix(0.3)

        # binary interaction parameter - tau_ij (0 if i=j, else to be estimated later but fixing to initialize)
        m.fs.properties.\
            tau[comp_2, comp_2].fix(0)
        m.fs.properties.\
            tau[comp_2, comp_1].fix(init_tau_2_1)
        m.fs.properties.\
            tau[comp_1, comp_1].fix(0)
        m.fs.properties.\
            tau[comp_1, comp_2].fix(init_tau_1_2)

        # Initialize the flash unit
        print('Trying to Initialize')
        print(data["pressure"])
#         m.fs.state_block.initialize(outlvl=idaeslog.INFO_HIGH)
        m.fs.state_block.initialize(outlvl=50)
#         m.fs.state_block.initialize(outlvl=idaeslog.DEBUG)
        print('Finished with Init')

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
        m.fs.properties.\
            tau[comp_2, comp_1].setlb(-10)
        m.fs.properties.\
            tau[comp_2, comp_1].setub(10)

        m.fs.properties.\
            tau[comp_1, comp_2].setlb(-10)
        m.fs.properties.\
            tau[comp_1, comp_2].setub(10)

        # Return initialized flash model
        return m
    
    def SSE(m, data):
        """
        returns objective function expresion.
        Args:
            m: model
            data: pandas DataFrame with data
        Returns:
            objectuve function scaled expresion
        """
        expr = ((float(data["pressure"]) - m.fs.state_block.pressure)**2)
        
        return expr * scaling_fac

    #Create a list of vars to estimate
    variable_name = ["fs.properties.tau['" + comp_2 + "','" + comp_1 + "']",
                     "fs.properties.tau['" + comp_1 + "','" + comp_2 + "']"]
    
    solver_opts_dict = {'output_file':filename}
    
    pest = parmest.Estimator(NRTL_model, data, variable_name, SSE, tee=True, solver_options=solver_opts_dict)

    obj_value, parameters, a = pest.theta_est(calc_cov=True)

#     print_params(obj_value, parameters)
    
#     print("covariance_matrix",a)
    
    return parameters, obj_value, a


def print_params(obj_value, parameters):
    """
    returns printed parameters.
    Args:
        obj_value: objective function value
        parameters: estimated binary parameters
    Returns:
        printed parameters for binary interaction parameters
    """
    print("The SSE at the optimal solution is %0.6f" % obj_value)
    print()
    print("The values for the parameters are as follows:")
    for k,v in parameters.items():
        print(k, "=", v)
        

class Data_reg:

    def __init__(self, configuration, comp_1, comp_2, x_comp_1, x_comp_2):
        self.configuration = configuration
        self.comp_1 = comp_1
        self.comp_2 = comp_2
        self.x_comp_1 = x_comp_1
        self.x_comp_2 = x_comp_2
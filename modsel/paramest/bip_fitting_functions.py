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

def polynomial(file, configuration, comp_1, comp_2, x_comp_1, x_comp_2,
            init_temp = 323.15, init_press = 399800, init_x_c1 = 0.5, init_x_c2 = 0.5,
            init_kappa_1_2A = -0.3, init_kappa_2_1A = 0.5,
            init_kappa_1_2B = -0.3, init_kappa_2_1B = 0.5,
            init_kappa_1_2C = -0.3, init_kappa_2_1C = 0.5,
            init_kappa_1_2D = -0.3, init_kappa_2_1D = 0.5,
            eps = 0.0, scaling_fac = 1e-4, filename='filename_test.csv'):
    """
    Estimates kappa parameters for Peng Robinson equation.
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
        init_kappa_1_2 = initial guess for kappa parameter component 2-component 1
        init_kappa_2_1 = initial guess for kappa parameter component 1-component 2
        eps = extra
        scaling_fac = 1e-4)
    Returns:
        printed parameters for binary interaction parameters
    """
    data = file

    def PR_model(data):
        """
        Function which returns initilized model.
        Args:
            data: pandas DataFrame with data
        Returns:
            initialized model
        """
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

        # parameter - kappa_ij (set at 0.3, 0 if i=j)
        m.fs.properties.PR_kappa_A[comp_2, comp_2].fix(0)
        m.fs.properties.PR_kappa_A[comp_2, comp_1].fix(init_kappa_2_1A)
        m.fs.properties.PR_kappa_A[comp_1, comp_1].fix(0)
        m.fs.properties.PR_kappa_A[comp_1, comp_2].fix(init_kappa_1_2A)
        m.fs.properties.PR_kappa_B[comp_2, comp_2].fix(0)
        m.fs.properties.PR_kappa_B[comp_2, comp_1].fix(init_kappa_2_1B)
        m.fs.properties.PR_kappa_B[comp_1, comp_1].fix(0)
        m.fs.properties.PR_kappa_B[comp_1, comp_2].fix(init_kappa_1_2B)
        m.fs.properties.PR_kappa_C[comp_2, comp_2].fix(0)
        m.fs.properties.PR_kappa_C[comp_2, comp_1].fix(init_kappa_2_1C)
        m.fs.properties.PR_kappa_C[comp_1, comp_1].fix(0)
        m.fs.properties.PR_kappa_C[comp_1, comp_2].fix(init_kappa_1_2C)
        m.fs.properties.PR_kappa_D[comp_2, comp_2].fix(0)
        m.fs.properties.PR_kappa_D[comp_2, comp_1].fix(init_kappa_2_1D)
        m.fs.properties.PR_kappa_D[comp_1, comp_1].fix(0)
        m.fs.properties.PR_kappa_D[comp_1, comp_2].fix(init_kappa_1_2D)

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

    variable_name = ["fs.properties.PR_kappa_A['" + comp_2 + "','" + comp_1 + "']",
                     "fs.properties.PR_kappa_A['" + comp_1 + "','" + comp_2 + "']",
                     "fs.properties.PR_kappa_B['" + comp_2 + "','" + comp_1 + "']",
                     "fs.properties.PR_kappa_B['" + comp_1 + "','" + comp_2 + "']",
                     "fs.properties.PR_kappa_C['" + comp_2 + "','" + comp_1 + "']",
                     "fs.properties.PR_kappa_C['" + comp_1 + "','" + comp_2 + "']",
                     "fs.properties.PR_kappa_D['" + comp_2 + "','" + comp_1 + "']",
                     "fs.properties.PR_kappa_D['" + comp_1 + "','" + comp_2 + "']"]


    solver_opts_dict = {'output_file':filename}
    
    pest = parmest.Estimator(PR_model, data, variable_name, SSE, tee=False, solver_options=solver_opts_dict)

    obj_value, parameters, a= pest.theta_est(calc_cov=True)

#     print_params(obj_value, parameters)

#     print("covariance_matrix",a)
    
    return parameters, obj_value, a

def cuadratic(file, configuration, comp_1, comp_2, x_comp_1, x_comp_2,
            init_temp = 323.15, init_press = 399800, init_x_c1 = 0.5, init_x_c2 = 0.5,
            init_kappa_1_2A = -0.3, init_kappa_2_1A = 0.5,
            init_kappa_1_2B = -0.3, init_kappa_2_1B = 0.5,
            init_kappa_1_2C = -0.3, init_kappa_2_1C = 0.5,
            eps = 0.0, scaling_fac = 1e-4, filename='filename_test.csv'):#, terms=3):
    """
    Estimates kappa parameters for Peng Robinson equation.
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
        init_kappa_1_2 = initial guess for kappa parameter component 2-component 1
        init_kappa_2_1 = initial guess for kappa parameter component 1-component 2
        eps = extra
        scaling_fac = 1e-4)
    Returns:
        printed parameters for binary interaction parameters
    """
    data = file

    def PR_model(data):
        """
        Function which returns initilized model.
        Args:
            data: pandas DataFrame with data
        Returns:
            initialized model
        """
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

        # parameter - kappa_ij (set at 0.3, 0 if i=j)
        m.fs.properties.PR_kappa_A[comp_2, comp_2].fix(0)
        m.fs.properties.PR_kappa_A[comp_2, comp_1].fix(init_kappa_2_1A)
        m.fs.properties.PR_kappa_A[comp_1, comp_1].fix(0)
        m.fs.properties.PR_kappa_A[comp_1, comp_2].fix(init_kappa_1_2A)
        m.fs.properties.PR_kappa_B[comp_2, comp_2].fix(0)
        m.fs.properties.PR_kappa_B[comp_2, comp_1].fix(init_kappa_2_1B)
        m.fs.properties.PR_kappa_B[comp_1, comp_1].fix(0)
        m.fs.properties.PR_kappa_B[comp_1, comp_2].fix(init_kappa_1_2B)
        m.fs.properties.PR_kappa_C[comp_2, comp_2].fix(0)
        m.fs.properties.PR_kappa_C[comp_2, comp_1].fix(init_kappa_2_1C)
        m.fs.properties.PR_kappa_C[comp_1, comp_1].fix(0)
        m.fs.properties.PR_kappa_C[comp_1, comp_2].fix(init_kappa_1_2C)

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

    variable_name = ["fs.properties.PR_kappa_A['" + comp_2 + "','" + comp_1 + "']",
                     "fs.properties.PR_kappa_A['" + comp_1 + "','" + comp_2 + "']",
                     "fs.properties.PR_kappa_B['" + comp_2 + "','" + comp_1 + "']",
                     "fs.properties.PR_kappa_B['" + comp_1 + "','" + comp_2 + "']",
                     "fs.properties.PR_kappa_C['" + comp_2 + "','" + comp_1 + "']",
                     "fs.properties.PR_kappa_C['" + comp_1 + "','" + comp_2 + "']"]


    solver_opts_dict = {'output_file':filename}
    
    pest = parmest.Estimator(PR_model, data, variable_name, SSE, tee=False, solver_options=solver_opts_dict)

    obj_value, parameters, a= pest.theta_est(calc_cov=True)

#     print_params(obj_value, parameters)

#     print("covariance_matrix",a)
    
    return parameters, obj_value, a

def linear(file, configuration, comp_1, comp_2, x_comp_1, x_comp_2,
    init_temp = 323.15, init_press = 399800, init_x_c1 = 0.5, init_x_c2 = 0.5,
    init_kappa_1_2A = -0.3, init_kappa_2_1A = 0.5,
    init_kappa_1_2B = -0.3, init_kappa_2_1B = 0.5,
    eps = 0.1, scaling_fac = 1e-9, optional_params = 'Four', filename='filename_test.csv'):
    """
    Estimates kappa parameters for Peng Robinson equation.
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
        init_kappa_1_2 = initial guess for kappa parameter component 2-component 1
        init_kappa_2_1 = initial guess for kappa parameter component 1-component 2
        eps = extra
        scaling_fac = 1e-4)
    Returns:
        printed parameters for binary interaction parameters
    """
    data = file

    def PR_model(data):
        """
        Function which returns initilized model.
        Args:
            data: pandas DataFrame with data
        Returns:
            initialized model
        """
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

        # parameter - kappa_ij (set at 0.3, 0 if i=j)
        m.fs.properties.PR_kappa_A[comp_2, comp_2].fix(0)
        m.fs.properties.PR_kappa_A[comp_2, comp_1].fix(init_kappa_2_1A)
        m.fs.properties.PR_kappa_A[comp_1, comp_1].fix(0)
        m.fs.properties.PR_kappa_A[comp_1, comp_2].fix(init_kappa_1_2A)
        m.fs.properties.PR_kappa_B[comp_2, comp_2].fix(0)
        m.fs.properties.PR_kappa_B[comp_2, comp_1].fix(init_kappa_2_1B)
        m.fs.properties.PR_kappa_B[comp_1, comp_1].fix(0)
        m.fs.properties.PR_kappa_B[comp_1, comp_2].fix(init_kappa_1_2B)
        m.fs.properties.PR_kappa_C[comp_2, comp_2].fix(0)
        m.fs.properties.PR_kappa_C[comp_2, comp_1].fix(0)
        m.fs.properties.PR_kappa_C[comp_1, comp_1].fix(0)
        m.fs.properties.PR_kappa_C[comp_1, comp_2].fix(0)
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
        m.fs.properties.PR_kappa_A[comp_2, comp_1].setlb(-5)
        m.fs.properties.PR_kappa_A[comp_2, comp_1].setub(5)

        m.fs.properties.PR_kappa_A[comp_1, comp_2].setlb(-5)
        m.fs.properties.PR_kappa_A[comp_1, comp_2].setub(5)

        m.fs.properties.PR_kappa_B[comp_2, comp_1].setlb(-5)
        m.fs.properties.PR_kappa_B[comp_2, comp_1].setub(5)

        m.fs.properties.PR_kappa_B[comp_1, comp_2].setlb(-5)
        m.fs.properties.PR_kappa_B[comp_1, comp_2].setub(5)

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
    
    if optional_params == 'Four':

        variable_name = ["fs.properties.PR_kappa_A['" + comp_1 + "','" + comp_2 + "']",
                         "fs.properties.PR_kappa_A['" + comp_2 + "','" + comp_1 + "']",
                         "fs.properties.PR_kappa_B['" + comp_1 + "','" + comp_2 + "']",
                         "fs.properties.PR_kappa_B['" + comp_2 + "','" + comp_1 + "']"]
    
    elif optional_params == 'Opt1': #Fix kji_B or k[emimTF2N, R32]_B
        
        variable_name = ["fs.properties.PR_kappa_A['" + comp_1 + "','" + comp_2 + "']",
                         "fs.properties.PR_kappa_A['" + comp_2 + "','" + comp_1 + "']",
                         "fs.properties.PR_kappa_B['" + comp_1 + "','" + comp_2 + "']"]
    
    elif optional_params == 'Opt2': #Fix kij_B or k[R32, emimTF2N]_B
        
        variable_name = ["fs.properties.PR_kappa_A['" + comp_1 + "','" + comp_2 + "']",
                         "fs.properties.PR_kappa_A['" + comp_2 + "','" + comp_1 + "']",
                         "fs.properties.PR_kappa_B['" + comp_2 + "','" + comp_1 + "']"]

    solver_opts_dict = {'output_file':filename}
    
    pest = parmest.Estimator(PR_model, data, variable_name, SSE, tee=False, solver_options=solver_opts_dict)

    obj_value, parameters,a= pest.theta_est(calc_cov=True)

#     print_params(obj_value, parameters)

#     print("covariance_matrix",a)
    
    return parameters, obj_value, a

def constant(file, configuration, comp_1, comp_2, x_comp_1, x_comp_2,
    init_temp = 323.15, init_press = 399800, init_x_c1 = 0.5, init_x_c2 = 0.5,
    init_kappa_A_1_2 = -0.3, init_kappa_A_2_1 = 0.5, eps = 0.0, scaling_fac = 1e-4,read= True, filename='filename_test.csv',optional_params = 'Two'):
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
        eps = extra
        scaling_fac = 1e-4)
        read: determines if read csv file or not
    Returns:
        printed parameters for binary interaction parameters
    """
    if read:
        data = pd.read_csv(file)
    else:
        data = file

    def PR_model(data):
        """
        Function which returns initilized model.
        Args:
            data: pandas DataFrame with data
        Returns:
            initialized model
        """
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

        # parameter - kappa_A_ij (set at 0.3, 0 if i=j)
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
#         print(data["pressure"])
        m.fs.state_block.initialize(outlvl=50)#idaeslog.ERROR)#INFO_LOW)#(outlvl=idaeslog.CRITICAL)

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

    if optional_params == 'Two':
        
        variable_name = ["fs.properties.PR_kappa_A['" + comp_2 + "','" + comp_1 + "']",
                     "fs.properties.PR_kappa_A['" + comp_1 + "','" + comp_2 + "']"]
   
    elif optional_params == 'Opt1': #fit kij_A: k[R32,IL]
        
        variable_name = ["fs.properties.PR_kappa_A['" + comp_1 + "','" + comp_2 + "']"]
        
    elif optional_params == 'Opt2': #fit kji_A: k[IL,R32]
        
        variable_name = ["fs.properties.PR_kappa_A['" + comp_2 + "','" + comp_1 + "']"]
    
#     print(filename)
    
    solver_opts_dict = {'output_file':filename}
    
#     print(solver_opts_dict)
    
    pest = parmest.Estimator(PR_model, data, variable_name, SSE, tee=False, solver_options=solver_opts_dict)

    obj_value, parameters, a = pest.theta_est(calc_cov=True)

#     print_params(obj_value, parameters)
    
#     print("covariance_matrix",a)
    
    return parameters, obj_value, a

def constant_CPSA(file, configuration, comp_1, comp_2, x_comp_1, x_comp_2,
    init_temp = 323.15, init_press = 399800, init_x_c1 = 0.5, init_x_c2 = 0.5,
    init_kappa_A_1_2 = -0.3, init_kappa_A_2_1 = 0.5, init_IL_Pc = 2.92E6, init_IL_Tc = 906.91, init_omega = 0.4223, eps = 0.0, scaling_fac = 1e-4,read= True):
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
        eps = extra
        scaling_fac = 1e-4)
        read: determines if read csv file or not
    Returns:
        printed parameters for binary interaction parameters
    """
    if read:
        data = pd.read_csv(file)
    else:
        data = file

    def PR_model(data):
        """
        Function which returns initilized model.
        Args:
            data: pandas DataFrame with data
        Returns:
            initialized model
        """
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

        # CP values for IL - component 2
        m.fs.properties.emimTf2N.pressure_crit.fix(init_IL_Pc)
        m.fs.properties.emimTf2N.temperature_crit.fix(init_IL_Tc)
        m.fs.properties.emimTf2N.omega.fix(init_omega)
        
        # parameter - kappa_A_ij (set at 0.3, 0 if i=j)
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
        m.fs.state_block.initialize(outlvl=idaeslog.CRITICAL)

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
        m.fs.properties.PR_kappa_A[comp_2, comp_1].setlb(-5)
        m.fs.properties.PR_kappa_A[comp_2, comp_1].setub(5)

        m.fs.properties.PR_kappa_A[comp_1, comp_2].setlb(-5)
        m.fs.properties.PR_kappa_A[comp_1, comp_2].setub(5)
        
        m.fs.properties.emimTf2N.pressure_crit.setlb(0.5E6)
        m.fs.properties.emimTf2N.pressure_crit.setub(4E6)
        
        m.fs.properties.emimTf2N.temperature_crit.setlb(700)
        m.fs.properties.emimTf2N.temperature_crit.setub(1200)
        
        m.fs.properties.emimTf2N.omega.setlb(0.2)
        m.fs.properties.emimTf2N.omega.setlb(0.8)

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

        
    variable_name = ["fs.properties.PR_kappa_A['" + comp_2 + "','" + comp_1 + "']",
                 "fs.properties.PR_kappa_A['" + comp_1 + "','" + comp_2 + "']",
                 "fs.properties.emimTf2N.pressure_crit",
                 "fs.properties.emimTf2N.temperature_crit",
                 "fs.properties.emimTf2N.omega"]

    pest = parmest.Estimator(PR_model, data, variable_name, SSE, tee=True)

    obj_value, parameters, a = pest.theta_est(calc_cov=True)

    print_params(obj_value, parameters)
    
    print("covariance_matrix",a)
    
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
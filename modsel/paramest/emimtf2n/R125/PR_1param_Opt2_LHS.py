'''
Fit emimTF2N data

EoS: PR

Parameter T dependence: Constant, 1 parameter - emimTF2N, R125 pair

N (total fitting parameters): 1

Python file to run LHS multistart.

'''

## Import 
import idaes

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize

import sys
sys.path.append('../../')

from pyomo.environ import (Constraint,
                           Var,
                           ConcreteModel,
                           Expression,
                           Param,
                           Objective,
                           SolverFactory,
                           TransformationFactory,
                           value)
from pyomo.opt import TerminationCondition, SolverStatus

from idaes.core import FlowsheetBlock
import idaes.logger as idaeslog
# Import the Generic Parameter Block
from idaes.generic_models.properties.core.generic.generic_property import (
        GenericParameterBlock)
# Import unit models from the model library
from idaes.generic_models.unit_models import Flash
# Import degrees of freedom tool
from idaes.core.util.model_statistics import degrees_of_freedom

# parmest (binary_param2)
from bip_fitting_functions import constant

## Load Data
data_subset = pd.read_csv('r125_emimtf2n_subset.csv')

from hfc125_emimtf2n_PR import configuration 

m = ConcreteModel()
m.fs = FlowsheetBlock(default={"dynamic": False})
m.fs.properties = GenericParameterBlock(default=configuration)
m.fs.F101 = Flash(default={"property_package": m.fs.properties,
                           "has_heat_transfer": True,
                           "has_pressure_change": True})

## LHS Multistart
lhs = pd.read_csv('../../LHS_250x8.txt')
lhs = lhs[['Param1']]

#Bound - scaling
lhs['sc_param1'] = ((lhs['Param1']-lhs['Param1'].min())/(lhs['Param1'].max()-lhs['Param1'].min()))*2+(-1)

#Run through parameter sets
result = np.zeros((len(lhs),2))

for i in range(len(lhs)):

    try:
        print(i)
        parameters, obj_value, a = constant(data_subset, configuration, 'R125', 'emimTf2N', "x_R125", "x_emimTf2N", 
        init_temp =  283.1, init_press =   399300 , init_x_c1 =    0.448, init_x_c2 = 0.552,
        init_kappa_A_2_1 = lhs['sc_param1'][i], init_kappa_A_1_2 = 0, 
        eps = 0.1, scaling_fac = 1e-9 , read=False, optional_params='Opt2', filename='Data/Fits/LHS_Fits/Ipopt_Output/PR_1param_Opt2/PR_1param_Opt2_LHS_ipopt_params'+str(i)+'.txt')

#         print(parameters)

        result[i,:] = parameters[0], obj_value
    except:
        parameters, obj_value, a = np.nan, np.nan, np.nan

        result[i,:] = np.nan, np.nan

#Save result
lhs['PR_kappa_A[emimTf2N,R32]']=result[:,0]
lhs['SSR']=result[:,1]

lhs.to_csv('Data/Fits/LHS_Fits/PR_1param_Opt2_LHS.csv',sep=',')
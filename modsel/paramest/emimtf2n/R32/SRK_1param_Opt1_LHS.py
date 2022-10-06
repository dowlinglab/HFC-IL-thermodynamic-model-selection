'''
Fit emimTF2N data

EoS: SRK

Parameter T dependence: Constant, 1 parameter - R32,emimTF2N pair

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
from bip_fitting_functions_SRK import constant

## Load Data
data_subset = pd.read_csv('r32_emimtf2n_subset.csv')

from hfc32_emimtf2n_SRK import configuration 

m = ConcreteModel()
m.fs = FlowsheetBlock(default={"dynamic": False})
m.fs.properties = GenericParameterBlock(default=configuration)
m.fs.F101 = Flash(default={"property_package": m.fs.properties,
                           "has_heat_transfer": True,
                           "has_pressure_change": True})

## LHS Multistart
lhs = pd.read_csv('../../LHS_250x8.txt')
lhs = lhs[['Param2']]

#Bound - scaling
lhs['sc_param2'] = ((lhs['Param2']-lhs['Param2'].min())/(lhs['Param2'].max()-lhs['Param2'].min()))*2+(-1)

#Run through parameter sets
result = np.zeros((len(lhs),2))

for i in range(len(lhs)):
    
    print('Trying i=',i)

    try:

        parameters, obj_value, a = constant(data_subset, configuration, 'R32', 'emimTf2N', "x_R32", "x_emimTf2N", 
        init_temp =  283.1, init_press =   399300 , init_x_c1 =    0.448, init_x_c2 = 0.552,
        init_kappa_A_2_1 = 0, init_kappa_A_1_2 = lhs['sc_param2'][i],
        eps = 0.1, scaling_fac = 1e-9 , read=False, optional_params='Opt1', filename='Data/Fits/LHS_Fits/Ipopt_Output/SRK_1param_Opt1/SRK_1param_Opt1_LHS_ipopt_params'+str(i)+'.txt', )

#         print(parameters)

        result[i,:] = parameters[0], obj_value
    except:
        parameters, obj_value, a = np.nan, np.nan, np.nan

        result[i,:] = np.nan, np.nan

#Save result
lhs['SRK_kappa_A[R32,emimTf2N]']=result[:,0]
lhs['SSR']=result[:,1]

lhs.to_csv('Data/Fits/LHS_Fits/SRK_1param_Opt1_LHS.csv',sep=',')
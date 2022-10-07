'''
Fit bmimpf6 data

EoS: SRK

Parameter T dependence: Linear in [bmimpf6,R32] direction

N (total fitting parameters): 3

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
from bip_fitting_functions_SRK import linear

## Load Data
data_subset = pd.read_csv('r32_bmimpf6_subset.csv')

from hfc32_bmimpf6_SRK import configuration 

m = ConcreteModel()
m.fs = FlowsheetBlock(default={"dynamic": False})
m.fs.properties = GenericParameterBlock(default=configuration)
m.fs.F101 = Flash(default={"property_package": m.fs.properties,
                           "has_heat_transfer": True,
                           "has_pressure_change": True})

## LHS Multistart
lhs = pd.read_csv('../../LHS_250x8.txt')
lhs = lhs[['Param1','Param2','Param3']]

#Bound - scaling
lhs['sc_param1'] = ((lhs['Param1']-lhs['Param1'].min())/(lhs['Param1'].max()-lhs['Param1'].min()))*2+(-1)
lhs['sc_param2'] = ((lhs['Param2']-lhs['Param2'].min())/(lhs['Param2'].max()-lhs['Param2'].min()))*2+(-1)
lhs['sc_param3'] = ((lhs['Param3']-lhs['Param3'].min())/(lhs['Param3'].max()-lhs['Param3'].min()))*2+(-1)

#Run through parameter sets
result = np.zeros((len(lhs),4))

for i in range(len(lhs)):

    print(i)
    
    try:

        parameters, obj_value, a = linear(data_subset, configuration, 'R32', 'bmimpf6', "x_R32", "x_bmimpf6", 
        init_temp =  283.1, init_press =   399300 , init_x_c1 =    0.448, init_x_c2 = 0.552,
        init_kappa_2_1A = lhs['sc_param1'][i], init_kappa_1_2A = lhs['sc_param2'][i], 
        init_kappa_2_1B = lhs['sc_param3'][i], init_kappa_1_2B = 0, 
        eps = 0.1, scaling_fac = 1e-9 , optional_params = 'Opt2', filename='Data/Fits/LHS_Fits/Ipopt_Output/SRK_3params_Opt2/SRK_3params_Opt2_LHS_ipopt_params'+str(i)+'.txt')


        result[i,:] = parameters[0], parameters[1], parameters[2], obj_value
        
    except:
        parameters, obj_value, a = np.nan, np.nan, np.nan

        result[i,:] = np.nan, np.nan, np.nan, np.nan

#Save result
lhs['SRK_kappa_A[bmimpf6,R32]']=result[:,0]
lhs['SRK_kappa_A[R32,bmimpf6]']=result[:,1]
lhs['SRK_kappa_B[bmimpf6,R32]']=result[:,2]
lhs['SSR']=result[:,3]

lhs.to_csv('Data/Fits/LHS_Fits/SRK_3params_Opt2_LHS.csv',sep=',')
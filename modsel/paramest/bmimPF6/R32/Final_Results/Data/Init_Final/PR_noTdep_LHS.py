'''
Fit bmimpf6 data

EoS: PR

Parameter T dependence: Constant

N (total fitting parameters): 2

Python file to run LHS multistart.

'''

## Import 
import idaes

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize

import sys
sys.path.append('../../../../../')

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

from bip_fitting_functions import constant

## Load Data
data_subset = pd.read_csv('../../r32_bmimpf6_subset.csv')

from hfc32_bmimpf6_PR import configuration 

m = ConcreteModel()
m.fs = FlowsheetBlock(default={"dynamic": False})
m.fs.properties = GenericParameterBlock(default=configuration)
m.fs.F101 = Flash(default={"property_package": m.fs.properties,
                           "has_heat_transfer": True,
                           "has_pressure_change": True})

## LHS Multistart

lhs = pd.read_csv('../../../../../LHS_5000x8.txt')
lhs = lhs[['Param1','Param2']]

# #Bound - scaling
lhs['sc_param1'] = ((lhs['Param1']-lhs['Param1'].min())/(lhs['Param1'].max()-lhs['Param1'].min()))*4+(-2)
lhs['sc_param2'] = ((lhs['Param2']-lhs['Param2'].min())/(lhs['Param2'].max()-lhs['Param2'].min()))*4+(-2)

#Run through parameter sets
result = np.zeros((len(lhs),3))

for i in range(len(lhs)):
    
    print('Trying i:',i)

    try:
        parameters, obj_value, a = constant(data_subset, configuration, 'R32', 'bmimpf6', "x_R32", "x_bmimpf6", 
        init_temp =  283.1, init_press =   399300 , init_x_c1 =    0.448, init_x_c2 = 0.552,
        init_kappa_A_2_1 = lhs['sc_param1'][i], init_kappa_A_1_2 = lhs['sc_param2'][i], 
        eps = 0.1, scaling_fac = 1e-9 , read=False, filename='Ipopt_Output/PR_noTdep/PR_noTdep_LHS_ipopt_params_final_'+str(i)+'.txt')

        print(parameters)

        result[i,:] = parameters[0], parameters[1], obj_value
    except:
        parameters, obj_value, a = np.nan, np.nan, np.nan

        result[i,:] = np.nan, np.nan, np.nan

#Save result
lhs['PR_kappa_A[bmimpf6,R32]']=result[:,0]
lhs['PR_kappa_A[R32,bmimpf6]']=result[:,1]
lhs['SSR']=result[:,2]

lhs.to_csv('LHS_Output/PR_noTdep_LHS_final.csv',sep=',')

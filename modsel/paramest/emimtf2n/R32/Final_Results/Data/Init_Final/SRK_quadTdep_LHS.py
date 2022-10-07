'''
Fit emimTF2N data

EoS: SRK

Parameter T dependence: quadratic

N (total fitting parameters): 6

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

from bip_fitting_functions_SRK import cuadratic

## Load Data
data_subset = pd.read_csv('../../r32_emimtf2n_subset.csv')

from hfc32_emimtf2n_SRK import configuration 

m = ConcreteModel()
m.fs = FlowsheetBlock(default={"dynamic": False})
m.fs.properties = GenericParameterBlock(default=configuration)
m.fs.F101 = Flash(default={"property_package": m.fs.properties,
                           "has_heat_transfer": True,
                           "has_pressure_change": True})

## LHS Multistart

lhs = pd.read_csv('../../../../../LHS_5000x8.txt')
lhs = lhs[['Param1','Param2','Param3','Param4','Param5','Param6']]

#Bound - scaling
lhs['sc_param1'] = ((lhs['Param1']-lhs['Param1'].min())/(lhs['Param1'].max()-lhs['Param1'].min()))*4+(-2)
lhs['sc_param2'] = ((lhs['Param2']-lhs['Param2'].min())/(lhs['Param2'].max()-lhs['Param2'].min()))*4+(-2)
lhs['sc_param3'] = ((lhs['Param3']-lhs['Param3'].min())/(lhs['Param3'].max()-lhs['Param3'].min()))*4+(-2)
lhs['sc_param4'] = ((lhs['Param4']-lhs['Param4'].min())/(lhs['Param4'].max()-lhs['Param4'].min()))*4+(-2)
lhs['sc_param5'] = ((lhs['Param5']-lhs['Param5'].min())/(lhs['Param5'].max()-lhs['Param5'].min()))*4+(-2)
lhs['sc_param6'] = ((lhs['Param6']-lhs['Param6'].min())/(lhs['Param6'].max()-lhs['Param6'].min()))*4+(-2)

#Run through parameter sets
result = np.zeros((len(lhs),7))

for i in range(len(lhs)):

    print('Trying i:',i)
    
    try:

        parameters, obj_value, a = cuadratic(data_subset, configuration, 'R32', 'emimTf2N', "x_R32", "x_emimTf2N", 
        init_temp =  283.1, init_press =   399300 , init_x_c1 =    0.448, init_x_c2 = 0.552,
        init_kappa_2_1A = lhs['sc_param1'][i], init_kappa_1_2A = lhs['sc_param2'][i], 
        init_kappa_2_1B = lhs['sc_param3'][i], init_kappa_1_2B = lhs['sc_param4'][i],
        init_kappa_2_1C = lhs['sc_param5'][i], init_kappa_1_2C = lhs['sc_param6'][i],
        eps = 0.1, scaling_fac = 1e-7, filename='Ipopt_Output/SRK_quadTdep/SRK_quadTdep_LHS_ipopt_params_final_'+str(i)+'.txt')

        print(parameters)

        result[i,:] = parameters[0], parameters[1], parameters[2], parameters[3], parameters[4], parameters[5], obj_value
    
    except:
        parameters, obj_value, a = np.nan, np.nan, np.nan

        result[i,:] = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan

#Save result
lhs['SRK_kappa_A[emimTf2N,R32]']=result[:,0]
lhs['SRK_kappa_A[R32,emimTf2N]']=result[:,1]
lhs['SRK_kappa_B[emimTf2N,R32]']=result[:,2]
lhs['SRK_kappa_B[R32,emimTf2N]']=result[:,3]
lhs['SRK_kappa_C[emimTf2N,R32]']=result[:,4]
lhs['SRK_kappa_C[R32,emimTf2N]']=result[:,5]
lhs['SSR']=result[:,6]

lhs.to_csv('LHS_Output/SRK_quadTdep_LHS_final.csv',sep=',')

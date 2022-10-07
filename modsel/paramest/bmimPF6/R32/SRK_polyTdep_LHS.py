'''
Fit bmimpf6 data

EoS: SRK

Parameter T dependence: order 3 polynomial

N (total fitting parameters):8

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
from bip_fitting_functions_SRK import polynomial

## Load Data
data_subset = pd.read_csv('r32_bmimpf6_subset.csv')

from hfc32_bmimpf6_SRK_polynomial import configuration 

m = ConcreteModel()
m.fs = FlowsheetBlock(default={"dynamic": False})
m.fs.properties = GenericParameterBlock(default=configuration)
m.fs.F101 = Flash(default={"property_package": m.fs.properties,
                           "has_heat_transfer": True,
                           "has_pressure_change": True})

## LHS Multistart

#Try initialized from SRK quad using initial scaled parameters for Params 1-6, LHS for params 7-8
lhs = pd.read_csv('Data/Fits/LHS_Fits/SRK_polyTdep_fromquadLHS_5000_ext.csv')

# #Try initialized from PR quad using initial act parameters for Params 1-6, LHS for params 7-8
# lhs = pd.read_csv('Data/Fits/LHS_Fits/PR_polyTdep_actLHS.csv')

# lhs = pd.read_csv('../../LHS_5000x8.txt')
# lhs = lhs[['Param1','Param2','Param3','Param4','Param5','Param6','Param7','Param8']]

# #Bound - scaling
# lhs['sc_param1'] = ((lhs['Param1']-lhs['Param1'].min())/(lhs['Param1'].max()-lhs['Param1'].min()))*4+(-2)
# lhs['sc_param2'] = ((lhs['Param2']-lhs['Param2'].min())/(lhs['Param2'].max()-lhs['Param2'].min()))*4+(-2)
# lhs['sc_param3'] = ((lhs['Param3']-lhs['Param3'].min())/(lhs['Param3'].max()-lhs['Param3'].min()))*4+(-2)
# lhs['sc_param4'] = ((lhs['Param4']-lhs['Param4'].min())/(lhs['Param4'].max()-lhs['Param4'].min()))*4+(-2)
# lhs['sc_param5'] = ((lhs['Param5']-lhs['Param5'].min())/(lhs['Param5'].max()-lhs['Param5'].min()))*4+(-2)
# lhs['sc_param6'] = ((lhs['Param6']-lhs['Param6'].min())/(lhs['Param6'].max()-lhs['Param6'].min()))*4+(-2)
# lhs['sc_param7'] = ((lhs['Param7']-lhs['Param7'].min())/(lhs['Param7'].max()-lhs['Param7'].min()))*4+(-2)
# lhs['sc_param8'] = ((lhs['Param8']-lhs['Param8'].min())/(lhs['Param8'].max()-lhs['Param8'].min()))*4+(-2)

#Run through parameter sets
result = np.zeros((len(lhs),9))

for i in range(len(lhs)):

    try:
        print(i)
        parameters, obj_value, a = polynomial(data_subset, configuration, 'R32', 'bmimpf6', "x_R32", "x_bmimpf6", 
        init_temp =  283.1, init_press =   399300 , init_x_c1 =    0.448, init_x_c2 = 0.552,
        init_kappa_2_1A = lhs['sc_param1'][i], init_kappa_1_2A = lhs['sc_param2'][i],
        init_kappa_2_1B = lhs['sc_param3'][i], init_kappa_1_2B = lhs['sc_param4'][i],
        init_kappa_2_1C = lhs['sc_param5'][i], init_kappa_1_2C = lhs['sc_param6'][i],
        init_kappa_2_1D = lhs['sc_param7'][i], init_kappa_1_2D = lhs['sc_param8'][i],
        eps = 0.1, scaling_fac = 1e-7, filename='Data/Fits/LHS_Fits/Ipopt_Output/SRK_polyTdep/SRK_polyTdep_LHS_ipopt_params_fromquad_5000ext_'+str(i)+'.txt')

        print(parameters)

        result[i,:] = parameters[0], parameters[1], parameters[2], parameters[3], parameters[4], parameters[5], parameters[6], parameters[7], obj_value
    
    except:
        parameters, obj_value, a = np.nan, np.nan, np.nan

        result[i,:] = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan

#Save result
lhs['SRK_kappa_A[bmimpf6,R32]']=result[:,0]
lhs['SRK_kappa_A[R32,bmimpf6]']=result[:,1]
lhs['SRK_kappa_B[bmimpf6,R32]']=result[:,2]
lhs['SRK_kappa_B[R32,bmimpf6]']=result[:,3]
lhs['SRK_kappa_C[bmimpf6,R32]']=result[:,4]
lhs['SRK_kappa_C[R32,bmimpf6]']=result[:,5]
lhs['SRK_kappa_D[bmimpf6,R32]']=result[:,6]
lhs['SRK_kappa_D[R32,bmimpf6]']=result[:,7]
lhs['SSR']=result[:,8]

lhs.to_csv('Data/Fits/LHS_Fits/SRK_polyTdep_LHS_fromquad_5000ext.csv',sep=',')
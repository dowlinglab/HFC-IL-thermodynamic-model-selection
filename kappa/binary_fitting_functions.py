#!/usr/bin/env python
# coding: utf-8

# In[7]:

# import functions
import idaes
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from pyomo.environ import (Constraint,
                           Var,
                           ConcreteModel,
                           Expression,
                           Param,
                           Objective,
                           SolverFactory,
                           TransformationFactory,
                           value)

from idaes.core import FlowsheetBlock
import idaes.logger as idaeslog
from pyomo.opt import TerminationCondition, SolverStatus
# Import the Generic Parameter Block
from idaes.generic_models.properties.core.generic.generic_property import (
        GenericParameterBlock)
# Import unit models from the model library
from idaes.generic_models.unit_models import Flash
# Import degrees of freedom tool
from idaes.core.util.model_statistics import degrees_of_freedom
# parmest (binary_param2)
from binary_param2 import binary_params_peng


def calc_outlet_xP(data, m, HFC, IL, k21,k12,verbose=False):
    '''Inputs
            Data - csv file for x,T,P data HFC-IL
            m    - model
            HFC  - string with name of HFC (i.e."R125")
            IL   - string with name of IL (i.e. "emimTf2N")
            k21  - binary parameter IL-HFC
            k12  - binary parameter HFC-IL

            '''
    # create zeros arrays to be filled with x1 and P1
    HFC_x = np.zeros((len(data)))
    P = np.zeros((len(data)))

    # data labels
    x_IL= "x_"+IL
    x_HFC= "x_"+HFC

#     #print
#     print(x_IL)
#     print(x_HFC)

    # model constraints
    m.fs.liq = Param(mutable=True,default=0.040)
    m.fs.liquid = Constraint(expr=m.fs.F101.liq_outlet.mole_frac_comp[0, IL] == m.fs.liq)

    for i in range(len(data)):
        m.fs.liq = data[x_IL].iloc[i]
        m.fs.F101.inlet.flow_mol.fix(1)
        m.fs.F101.inlet.temperature.fix(float(data["temperature"].iloc[i]))
        m.fs.F101.inlet.pressure.fix(float(data["pressure"].iloc[i]))
        m.fs.F101.inlet.mole_frac_comp[0, HFC].fix(float(data[x_HFC].iloc[i])+0.1)
        m.fs.F101.inlet.mole_frac_comp[0, IL].fix(float(1-(data[x_HFC].iloc[i]+0.1)))
        m.fs.F101.vap_outlet.temperature.fix(float(data["temperature"].iloc[i]))
        m.fs.properties.PR_kappa_A[IL, HFC].fix(k21) # (-0.20093)
        m.fs.properties.PR_kappa_A[HFC, IL].fix(k12) # (-0.05619)

        if verbose:
            DOF_final = degrees_of_freedom(m)
            print("The final DOF is {0}".format(DOF_final))

        # solver
        m.fs.F101.initialize(outlvl=idaeslog.CRITICAL)
        solver = SolverFactory('ipopt')
        solver.options = {'tol': 1e-6}
        status = solver.solve(m, tee = False)
    #     m.fs.F101.report()
        if (status.solver.status == SolverStatus.ok) and (status.solver.termination_condition == TerminationCondition.optimal):
            HFC_x[i] = value(m.fs.F101.liq_outlet.mole_frac_comp[0,HFC])
            P[i] = value(m.fs.F101.vap_outlet.pressure[0])
        else:
            print('Infeasible.')

    return(HFC_x,P)

# sensitivity analysis plot function
def plot_sens_analysis(data,m,HFC,IL,kappa_A21,kappa_A12,T_label):
    x1,P1= calc_outlet_xP(data,m, HFC, IL, kappa_A21[0],kappa_A12[0])
    plt.plot(x1,P1,"b+",label="K21= "+str(kappa_A21[0])+ "K12= "+ str(kappa_A12[0]))
    x2,P2= calc_outlet_xP(data,m, HFC, IL, kappa_A21[0],kappa_A12[1])
    plt.plot(x2,P2,"g+",label="K21= "+str(kappa_A21[0])+ "K12= "+ str(kappa_A12[1]))
    x3,P3= calc_outlet_xP(data,m, HFC, IL,kappa_A21[0],kappa_A12[2])
    plt.plot(x3,P3,"m+",label="K21= "+str(kappa_A21[0])+ "K12= "+ str(kappa_A12[2]))
    x4,P4= calc_outlet_xP(data,m, HFC, IL,kappa_A21[0],kappa_A12[3])
    plt.plot(x4,P4,"c+",label="K21= "+str(kappa_A21[0])+ "K12= "+ str(kappa_A12[3]))
    x5,P5= calc_outlet_xP(data,m, HFC, IL,kappa_A21[1],kappa_A12[0])
    plt.plot(x5,P5,"b*",label="K21= "+str(kappa_A21[1])+ "K12= "+ str(kappa_A12[0]))
    x6,P6= calc_outlet_xP(data,m, HFC, IL,kappa_A21[1],kappa_A12[1])
    plt.plot(x6,P6,"g*-",label="K21= "+str(kappa_A21[1])+ "K12= "+ str(kappa_A12[1]))
    x7,P7= calc_outlet_xP(data,m, HFC, IL,kappa_A21[1],kappa_A12[2])
    plt.plot(x7,P7,"m*",label="K21= "+str(kappa_A21[1])+ "K12= "+ str(kappa_A12[2]))
    x8,P8= calc_outlet_xP(data,m, HFC, IL,kappa_A21[1],kappa_A12[3])
    plt.plot(x8,P8,"c*",label="K21= "+str(kappa_A21[1])+ "K12= "+ str(kappa_A12[3]))
    x9,P9= calc_outlet_xP(data,m, HFC, IL,kappa_A21[2],kappa_A12[0])
    plt.plot(x9,P9,"b^",label="K21= "+str(kappa_A21[2])+ "K12= "+ str(kappa_A12[0]))
    x10,P10= calc_outlet_xP(data,m, HFC, IL,kappa_A21[2],kappa_A12[1])
    plt.plot(x10,P10,"g^",label="K21= "+str(kappa_A21[2])+ "K12= "+ str(kappa_A12[1]))
    x11,P11= calc_outlet_xP(data,m, HFC, IL,kappa_A21[2],kappa_A12[2])
    plt.plot(x11,P11,"m^",label="K21= "+str(kappa_A21[2])+ "K12= "+ str(kappa_A12[2]))
    x12,P12= calc_outlet_xP(data,m, HFC, IL,kappa_A21[2],kappa_A12[3])
    plt.plot(x12,P12,"c^",label="K21= "+str(kappa_A21[2])+ "K12= "+ str(kappa_A12[3]))
    x13,P13= calc_outlet_xP(data,m, HFC, IL,kappa_A21[3],kappa_A12[0])
    plt.plot(x13,P13,"b>",label="K21= "+str(kappa_A21[3])+ "K12= "+ str(kappa_A12[0]))
    x14,P14= calc_outlet_xP(data,m, HFC, IL,kappa_A21[3],kappa_A12[1])
    plt.plot(x14,P14,"g>",label="K21= "+str(kappa_A21[3])+ "K12= "+ str(kappa_A12[1]))
    x15,P15= calc_outlet_xP(data,m, HFC, IL,kappa_A21[3],kappa_A12[2])
    plt.plot(x15,P15,"m>",label="K21= "+str(kappa_A21[3])+ "K12= "+ str(kappa_A12[2]))
    x16,P16= calc_outlet_xP(data,m, HFC, IL,kappa_A21[3],kappa_A12[3])
    plt.plot(x16,P16,"c>",label="K21= "+str(kappa_A21[3])+ "K12= "+ str(kappa_A12[3]))
    # plot data and configuration
    plt.plot(data["x_"+HFC],data["pressure"],"r.",label="data")
    plt.title("Isotherm "+ "["+HFC+"]["+IL+"] at "+ T_label)
    plt.ylabel('Pressure (Pa)')
    plt.xlabel("x_"+HFC)
    plt.ylim(0,1100000)
    plt.grid(True)
    plt.legend(bbox_to_anchor=(1.05, 1))
    plt.show()

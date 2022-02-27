"""
Add a descriptive sentence

Created by Bridgette Befort, University of Notre Dame, 2022
in collaboration with Alejandro Garciadiego, Alex Dowling, Ed Maginn

Functions defined in this file:

custom_ipopt - returns ipopt solver/executable
System Parameters (SP) - class with attributes of components, critical properties, EoS constants
Fitted Parameters (FP) - class with attributes of BIPs
cubic roots - returns solutions to a cubic equation (can be used with volume, compressibility)
fugacity - returns a dictionary of values related to solving for fugacity using vdW EoS, can be scaled (compressibility cubic root) or unscaled (V cubic root)
inexact newton - calculates root using newton's method
calc_P - calculates the Pressure using fugacity function and inexact newton's method
create_FP - makes a fitted parameters dictionary
regression func - returns fugacity residual
calculate_fugacity_detailed - calculates the fugacity for a given pressure and gives detailed results
plot_isotherms - plots isotherm using 1 set of parameters
plot_isotherms3 - plots isothmers for 2 different parameter sets
create model - builds Pyomo model
extract results - builds dictionary of Pyomo model values
compare dictionaries - compares Scipy and Pyomo model results
calc_P_pyomo - calculates pressure from Pyomo fitted parameters using pyomo model to solve
solve_Pyomo_model - solves the model using ipopt
plot_isotherms_Pyomo - uses BIPs fitted from Pyomo model and calc_P_pyomo function to plot pressures over a series of compositions
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import math
import scipy.optimize as optimize

def custom_ipopt():
    ''' Return compiled version of Ipopt
    
    Arguments:
        None
        
    Returns:
        solver
    
    '''

    import socket
    name = socket.gethostname()
    print(name)
   
    if name == "Bridgette’s MacBook Pro":
        print("Loading custom Ipopt...")
        solver = SolverFactory('ipopt',executable="/Users/bridgettebefort/src/CoinIpopt/build/bin/ipopt")
        solver.options['linear_solver'] = "ma57"
    else:
        print("Loading default Ipopt.")
        solver = SolverFactory('ipopt')
        
    return solver

class SystemParameters:

    ## Class attributes
    '''
    # list of chemical species
    components
    
    # dictionaries of critical properties
    Tc
    Pc
    '''

    def __init__(self, components, Tc, Pc):
        
        self.components = components
        self.Tc = Tc
        self.Pc = Pc

class FittedParameters:

    def __init__(self, k):
        self.k = k
        
        
def cubic_roots(a,b,c,phase="all",verbose=False):
    ''' Return cubic root corresponding to the specified phase
        for equation
        
        f(Z) = Z**3 + a*Z**2 + b*Z + c
        
        if three real roots exist,
        return the largest for vapor
        return the smallest for liquid
        
        if only one real root exists,
        return that one
    
    Args:
        a: dimensionless coefficient (float)
        b: dimensionless coefficient (float)
        c: dimensionless coefficient (float)
        phase: "L" for liquid, "V" for vapor (string), "all" for all three roots
        
    Returns:
        Z: root (float)
    '''
    
    # calculate intermediates
    Q = (a*a - 3*b)/9
    R = (2*a**3 - 9*a*b + 27*c) / 54
    M = R**2 - Q**3
    if verbose==True:
        print("Q",Q)
        print("R",R)
        print("M",M)
    if M < 0:
        # Three real roots exist
        thetaZ = math.acos(R / Q**(3/2))
        sQ = math.sqrt(Q)
        
        
        Z = [ -(2*sQ*math.cos(thetaZ/3)) - a/3 , # root 1
              -(2*sQ*math.cos(thetaZ/3 + 2*math.pi/3)) - a/3, # root 2
              -(2*sQ*math.cos(thetaZ/3 - 2*math.pi/3)) - a/3 # root 3
        ]
        
      
        
    else:
        # One real root exists
        # calculate sign of R
        if R < 0:
            r = -1
        else:
            r = 1
        
        S = -r*(math.fabs(R) + math.sqrt(M))**(1/3)
        
        if S == 0:
            T = 0
        else:
            T = Q/S
        
        Z = S + T - a/3
    
    if verbose==True:
        print("Z",Z)
    if phase == "L":
        # return the smallest root, which corresponds with the liquid phase
        return np.min(Z)
    elif phase == "V":
        # return the largest root, which corresponds with the vapor phase
        return np.max(Z)
        # return all of the returns
    elif phase == "" or phase =="all":
        return Z
    else:
        # give a warning message
        print("Only phase='L', phase='V' and phase='all' supported.")
        print("All real roots found:")
        print(Z)
        return Z
    
#Define fugacity function of which to find roots
def fugacity(P,T,x1,k,SP,verbose=False,scaled=True,unit_change=True,one_l=True):
    
    '''
    Compute liquid and vapor fugacity of CO2
    
    Input:
    P - guess pressure of system, MPa, to be determined (scalar)
    T - temperature, K (scalar)
    x_1 - mole fraction of component 1 (scalar)
    k - fitted parameters
    SP - SystemsParameters
    verbose - Boolean to toggle on/off extra print statements
    
    Output:
        Fugacity Residual (scalar)
    
    '''
    
#     P = max(P,1e-8)
    
    if verbose:
        print("\n\nP=",P)
        for i in SP.components:
            for j in SP.components:
                print("k =",k[i,j])
    
    #Constants
    if unit_change==True:
        R = 8.314e-3 #m3*MPa/K/kmol
    else:
        R = 8.314e-6 #m3*MPa/K/mol
    
    comps = SP.components
    
    #MoleFrac
    MoleFrac = {}
    MoleFrac[comps[0]] = x1 # solute
    MoleFrac[comps[1]] = 1 - x1 # IL
    
    if verbose:
        print("Mole Fractions = ",MoleFrac)
    
    #Find pure a - simple vdW eq. 1, 2/21/22
    a_pure = {}
    for i in SP.components:
        a_pure[i] = 0.421875*R**2*SP.Tc[i]**2/SP.Pc[i]
    if verbose:
        print("pure a = ",a_pure)
    
    #Find pure b - simple vdW eq. 2, 2/21/22
    b_pure = {}
    for i in SP.components:
        b_pure[i] = 0.125*R*SP.Tc[i]/SP.Pc[i]
    if verbose:
        print("pure b = ",b_pure)
        
    #Find mixture a - simple vdW eq. 3, 2/21/22
    a_mix = []
    a_mix = 0
    for i in SP.components:
        for j in SP.components:
            a_mix += math.sqrt(a_pure[i]*a_pure[j])*(1 - k[i,j])*MoleFrac[i]*MoleFrac[j]
    if verbose:
        print("mixture a =",a_mix)
        
    #Find mixture b - simple vdW eq. 4, 2/21/22
    b_mix = []
    b_mix = 0
    for i in SP.components:
        b_mix += MoleFrac[i]*b_pure[i]
    if verbose:
        print("mixture b =",b_mix)
    
    #Find Z - A, B different for liquid vs. vapor, simple vdW eq. 8, 2/21/22
    
    # Liquid A and B
    A_star_L = a_mix*P/R**2/T**2 
    B_star_L = b_mix*P/R/T 

    # Vapor A and B
    A_star_V = a_pure[comps[0]]*P/R**2/T**2
    B_star_V = b_pure[comps[0]]*P/R/T

    # Liquid phase (IL + HFC)
    alpha_l_comp = -(1+B_star_L)
    beta_l_comp = A_star_L
    gamma_l_comp = -A_star_L*B_star_L
    Z_compressibility_L = cubic_roots(alpha_l_comp,beta_l_comp,gamma_l_comp,phase="L",verbose=verbose)
    V_L = Z_compressibility_L*R*T/P

    # Vapor phase (HFC only)

    alpha_v_comp = -(1+B_star_V)
    beta_v_comp = A_star_V
    gamma_v_comp = -A_star_V*B_star_V
    Z_compressibility_V = cubic_roots(alpha_v_comp,beta_v_comp,gamma_v_comp,phase="V",verbose=verbose)
    V_V = Z_compressibility_V*R*T/P

    alpha_l_save = alpha_l_comp
    beta_l_save = beta_l_comp
    gamma_l_save = gamma_l_comp
    alpha_v_save = alpha_v_comp
    beta_v_save = beta_v_comp
    gamma_v_save = gamma_v_comp
        
    if verbose:
        print("A_star_L = ", A_star_L)
        print("B_star_L = ", B_star_L)
        print("V_L = ",V_L)
        print("Z_compressibility_L = ", Z_compressibility_L)
        
        print("A_star_V = ", A_star_V)
        print("B_star_V = ", B_star_V)
        print("V_V = ",V_V)
        print("Z_compressibility_V = ", Z_compressibility_V)
        
        print("alpha_l:",alpha_l_save)
        print("beta_l:",beta_l_save)
        print("gamma_l:",gamma_l_save)
        
        print("alpha_v:",alpha_v_save)
        print("beta_v:",beta_v_save)
        print("gamma_v:",gamma_v_save)
        
    #Find coefficients in log fugacity term
    B = {}
    for i in SP.components:
        B[i] = b_pure[i]*P/R/T
    if verbose:
        print("B = ",B)
    
    B_mix = P*b_mix/R/T #same as B_star_L
    
    A_mix = {}
    for i in SP.components:
        A_mix[i] = sum(math.sqrt(a_pure[i]*a_pure[j])*MoleFrac[j]*(1 - k[i,j])*P/R/T  for j in SP.components) #if i == j)
    
    #Find liquid fugacity
    lnphi_L = {}
    for i in SP.components:
        lnphi_L[i] = B[i]/(Z_compressibility_L - B_mix) - np.log(Z_compressibility_L - B_mix) - 2*A_mix[i]/Z_compressibility_L
        
    phi_L_HFC = np.exp(lnphi_L[comps[0]])
        
    #Find vapor fugacity
    lnphi_V = B[comps[0]]/(Z_compressibility_V - B[comps[0]]) - np.log(Z_compressibility_V - B[comps[0]]) - 2*a_pure[comps[0]]*P/R/T/Z_compressibility_V
    
    phi_V_HFC = np.exp(lnphi_V)

    #Check vdW pressure calc
    def check_P(V_check, a_check, b_check):
        return R*T/(V_check - b_check) - a_check / V_check**2
    
    if verbose:
        print("V_L = ",V_L)
        print("corresponding pressure:",check_P(V_L, a_mix, b_mix))
        
        print("\nV_V = ",V_V)
        print("corresponding pressure:",check_P(V_V, a_pure[comps[0]], b_pure[comps[0]]))
        
    def check_equil(x,lnphiL_check,lnphiV_check):
        return np.log(x) + lnphiL_check - lnphiV_check
        
    if verbose:
        print("ln(x)+lnphi_L - lnphi_V = ",check_equil(MoleFrac[comps[0]],lnphi_L[comps[0]],lnphi_V))
        
#     assert V_L >0.00001, "V_L is too small"
#     assert V_V <10000, "V_V is too large"
        
    return {"x1":MoleFrac[comps[0]], "a_pure":a_pure, "a_mix":a_mix, "b_pure":b_pure, "b_mix":b_mix, "alpha_l":alpha_l_save,"beta_l":beta_l_save,"gamma_l":gamma_l_save,"VL":V_L,"A_star_L":A_star_L,"B_star_L":B_star_L,"Z_compressibility_L":Z_compressibility_L, "ln_phi_L":lnphi_L[comps[0]],"phi_L_HFC":phi_L_HFC,"alpha_v":alpha_v_save,"beta_v":beta_v_save,"gamma_v":gamma_v_save, "VV":V_V, "A_star_V":A_star_V,"B_star_V":B_star_V, "Z_compressibility_V":Z_compressibility_V,"ln_phi_V":lnphi_V,"phi_V_HFC":phi_V_HFC,"P":P}
    

def inexact_newton(f,x0,delta = 1.0e-7, epsilon=1.0e-6,max_iter=50, LOUD=False,bnds=None):
    """Find the root of the function f via Newton-Raphson method
    Args:
        f: function to find root of
        x0: initial guess
        delta: finite difference parameter
        epsilon: tolerance
        
    Returns:
        estimate of root
    """
    x = x0
    if (LOUD):
        print("x0 =",x0)
    iterations = 0
    while (np.fabs(f(x)) > epsilon):
        fx = f(x)
        
        fxdelta = f(x + delta)
        slope = (fxdelta - fx)/delta
        
        if (LOUD):
            print("x_",iterations+1,"=",x,"-",fx,"/",slope,"=",x - fx/slope)
        x = x - fx/slope
        if bnds is not None:
            x = max(x,bnds[0])
            x = min(x,bnds[1])
        iterations += 1
        assert iterations <max_iter,"Maximum number of iterations exceeded"
    
    if LOUD:
        print("It took",iterations,"iterations")
    return x #return estimate of root

def inexact_newton2(f,x0,delta = 1.0e-5, epsilon=1.0e-6, LOUD=True,max_iter=50,bnds=[-np.Inf,np.Inf]):
    """Find the root of the function f via Newton-Raphson method with line search
    Args:
        f: function to find root of [function]
        x0: initial guess [float]
        delta: finite difference parameter [float]
        epsilon: tolerance [float]
        LOUD: toggle on/off print statements [boolean]
        max_iter: maximum number of iterations [int]
        bnds: bounds on the x variable, [list with 2 floats]
       
    Returns:
        estimate of root [float]
    """
    x0 = x0.item()
    assert callable(f), "Warning: 'f' should be a Python function"
    assert type(x0) is float or type(x0) is int, "Warning: 'x0' should be a float or integer"
    assert type(delta) is float, "Warning: 'delta' should be a float"
    assert type(epsilon) is float, "Warning: 'eps' should be a float"
    assert type(max_iter) is int, "Warning: 'max_iter' should be an integer"
    assert max_iter >= 0, "Warning: 'max_iter' should be non-negative"
   
    x = x0
    # evaluate function 'f' at 'x0'
    fx = f(x)
   
    if (LOUD):
        print("x0 =",x0)
    iterations = 0
    keep_going = True
    terminate_message = ''
   
    # Check if the residual is close enough to zero
    while keep_going:
       

        # calculate 'slope'
        fxdelta = f(x + delta)
        slope = (fxdelta - fx)/delta
               
        # check if slope is zero
        if np.abs(slope) < 1E-14:
            keep_going = False
            terminate_message = 'Terminated. |slope| is less than 1E-14.'
       
        else:
       
            # apply line search
            keep_going, terminate_message, x, fx = linesearch(f,x,fx,slope,iterations,LOUD,bnds)
       
            iterations += 1
       
            # check if converged
            if np.fabs(fx) < epsilon:
                terminate_message = 'Converged!'
                keep_going = False

            # check if number of iterations exceeded
            elif iterations >= max_iter:
                keep_going = False
                terminate_message = 'Terminated. Maximum number of iterations exceeded.'
           

    print(terminate_message)
    assert terminate_message == 'Converged!'
    print("It took",iterations,"iterations")

    return x #return estimate of root

def linesearch(f,x,fx,slope,iterations,LOUD,bnds):
    ''' Basic linesearch algorithm to augment Newton's method
   
    Arguments:
        f: function to find root of [function]
        x: curret value for root [float]
        fx: function value at x [float]
        slope: slope estimated at x [float]
        iterations: current iteration [int]
        LOUD: flag to toggle on detailed printing [boolean]
        bnds: bounds for x
   
    Returns:
        keep_going: flag if the search should continue [boolean]
        terminate_message: termination message [string]
        x_new: accepted step [float]
        f_new: f(x_new) [float]
       
    '''

    # specify constant for 'sufficient decrease'
    eta = 0.6
   
    # initialize line search
    alpha = 1.0
    alpha_min = 1E-6

    # calculate the step direction
    d = -fx/slope
   
    # flag for linesearch
    linesearch = True
   
    # allocate variable to store result
    xtemp = 0.0
    fxtemp = 0.0
     
    while(linesearch):
               
        # calculate new trial point
        xtemp = x + alpha*d
        
        print_string = "x_{0} = {1:.3e} - {2:.3e} * {3:.3e} / {4:.3e} = {5:0.6e}".format(
            iterations+1,x,alpha,fx,slope,xtemp)
        
        if xtemp < bnds[0]:
            print(print_string +" < lower bound of {0:.3e}".format(bnds[0]))

        elif xtemp > bnds[1]:
            print(print_string +" > upper bound of {0:.3e}".format(bnds[1]))

        else:
            # evaluate f at new trial point
            fxtemp = f(xtemp)

            if (LOUD):
                print(print_string+",  f(x_{0}) = {1:0.3e}".format(iterations+1,fxtemp))
   
            # check if there is a sufficient decrease
            if np.abs(fxtemp) < eta*np.abs(fx):
                # line search converged
                linesearch = False
                keep_going = True
                terminate_message = ' '

        if linesearch:
            # descrease alpha
            alpha = alpha * 0.9

            if alpha < alpha_min:
                linesearch = False
                keep_going = False
                terminate_message = 'Terminated. Line search failure with alpha < {0:0.1e}.'.format(alpha_min)
   
    return keep_going, terminate_message, xtemp, fxtemp

def calc_P(T,x1,FP,SP,NM,P0=1.0,verbose=False):
    '''Compute the pressure with EOS at given T and x1
    
    Arguments:
        T: temperature in K (float)
        x1: mole fraction CO2 (float)
        FP: fitted parameters (FittedParameters type)
        SP: class with system parameters
        P0: initial guess
        NM: newton method type
    
    Returns:
        P: calculated pressure in MPa
    '''
    
    def my_f(logP):
        ''' Compute residual to equilibrium equation at guess pressure
        
        Arguments:
           logP: guess pressure in [MPa] 
        
        Returns:
            r: residual
        
        '''
#         print(logP)
        P = np.exp(logP)
        
        r = fugacity(P,T,x1,FP.k,SP,verbose=verbose)
        return r["phi_L_HFC"]*r["x1"] - r["phi_V_HFC"] #return equilibrium calc
    
    if NM==1:
        logP = inexact_newton(my_f,np.log(P0))#,bnds=[np.log(0.01),np.log(10.1)])
    else:
        logP = inexact_newton2(my_f,np.log(P0),bnds=[np.log(0.001),np.log(100.1)])
        
    return np.exp(logP)

def create_FP(theta,comps):
    ''' Load values of paramters into FittedParameter object
    
    Argument:
        theta: parameter vector (l, l, m, tau)
        comps: list of component names. comps[0] is solute, comps[1] is IL
        
    Return:
        FP: FittedParameters object
        
    '''
    
    assert type(comps) is list
    
    # unpack theta
    k_ = {(comps[0],comps[0]): 0,
        (comps[0],comps[1]): theta[0],
        (comps[1],comps[0]): theta[1],
        (comps[1],comps[1]): 0}


    FP = FittedParameters(k_)
    return FP

def regression_func(theta, T, x1, P, SP,NM,verbose=False):
    '''
    Function to define regression function for least-squares fitting
    Arguments:
        theta: parameter vector
        x: independent variable vector
        y: dependent variable vector (measurements)
        SP: class containing system parameters
    Returns:
        e: residual vector
    '''
    if len(theta)==2:
        FP = create_FP(theta, SP.components)
    
    n = len(T)
    e = np.zeros(n)
    
    for i in range(n):
        Pexp = P[i]
        Pcalc = P[i]
        
        if i==0:
            Pcalc = 0.001

        e[i] = Pexp - calc_P(T[i],x1[i],FP,SP,NM,Pcalc,verbose=verbose)
    
    return e

def calculate_fugacity_detailed(FP, SP,data, equilibrium=True, scaling=True):
    ''' Calculate the fugacity at specified T-P-x1 with Scipy model.
        Return detailed calculations in a dictionary of dictionaries
        
    Arguments:
    
    Returns:
        R - dictionary of dictionaries (for each datapoint)
    '''
    
    R = dict()
    
    T = data["T"].to_dict()
    P = data["P"].to_dict()
    x1 = data["x1"].to_dict()
    for m in range(len(data)):
    
        if equilibrium: #calculate equilibrium pressure
            P_ = calc_P(T[m],x1[m],FP, SP, P[m])
        else: #otherwise print out guess pressure
            P_ = P[m]
        
        #Calculate fugacity
        R[m] = fugacity(P_,T[m],x1[m],FP.k,SP)

        R[m]["P"] = P_ #replace pressure fugacity with this the calculated pressure
            
    return R
    

def plot_isotherms(FP,SP,NM,data,num_sets,verbose=False):
    '''Plot the isotherms
    
    Arugment:
        FP: fitted parameters
        SP: system parameters
        data: Pandas dataframe containing P-T-x1
    
    Returns:
        Nothing
        
    Other:
        Creates plot
    '''
    
    colors = {'S1':'red','S2':'blue','S3':'green','S4':'purple'}
    markers = {'S1':'o','S2':'s','S3':'*','S4':'^'}
    N = 20*num_sets
    results = np.zeros((N,2))
    for s in data.Set.unique():
        temp = data[data.Set == s]
        temperature = temp['T'].median()
        plt.plot(temp.x1,temp.P,color=colors[s],label=str(temperature)+" K",
                 marker=markers[s],linestyle=' ')
        n = 20
        P_calculated = []
        x_calculated = []
        count = 0
        for x in np.linspace(temp.x1.min(),temp.x1.max(),n):
            P0 = np.interp(x, temp.x1, temp.P)
            
            print("\nSolving for x=",x,"using P0=",P0,"MPa at T=",temperature,"K")
            P = calc_P(temperature,x,FP,SP,NM,P0,one_l=one_l,verbose=verbose)  
            results[count,:] = [x,P]
            P_calculated.append(P)
            x_calculated.append(x)
            count += 1
        print(P_calculated)
        plt.plot(x_calculated,P_calculated,color=colors[s],linestyle='-')

    plt.legend()
    plt.grid(False)
    plt.xlabel('Mole Fraction HFC in Liquid')
    plt.ylabel('Pressure [MPa]')
#     plt.savefig('R125_bmimPF6_Scipy.pdf')
    plt.show()
    
    return results


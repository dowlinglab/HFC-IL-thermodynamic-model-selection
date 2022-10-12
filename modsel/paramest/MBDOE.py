import matplotlib.pyplot as plt
import numpy as np
import pandas as pd 
import idaes
import scipy.stats as stats

# import pyomo.doe
from fim_doe import *

from idaes.core import FlowsheetBlock
import idaes.logger as idaeslog
# Import the Generic Parameter Block
from idaes.generic_models.properties.core.generic.generic_property import (
        GenericParameterBlock)
# Import unit models from the model library
from idaes.generic_models.unit_models import Flash
# Import degrees of freedom tool
from idaes.core.util.model_statistics import degrees_of_freedom

# Import configuration
from hfc32_emimtf2n_PR import configuration 

from generalize_functions import PRModels


class MBDOE: 
    def __init__(self, data_file, create_model, verbose=True):
        '''
        data_file: csv file including all design variables
        create_model: model object
        '''
        
        self.data_exp = pd.read_csv(data_file)
        self.create_model_object = create_model
        self.verbose = verbose
        
        # get parameter information
        self.param_dict = self.create_model_object.param_name_dict
        # param names
        self.param_name = list(self.param_dict.keys())
        
        if self.verbose:
            print("parameter set:", self.param_dict)
        
        
    def sumDOE(self, exp_idx_set, scale_opt=False):
        
        num_param = len(self.param_name)
        
        totalFIM = [[0]*num_param for i in range(num_param)]
        
        failed_set = []
        
        for i in exp_idx_set:
            print("==========Experiment index:", i, "===============")
            try:
                res = self.doe(i, scale=scale_opt)
            
                totalFIM += res.FIM
            except:
                failed_set.append(i)
                print("Failure initialization!")
            
            
        print("Failed set:", failed_set)
        if self.verbose: 
            print('======Result summary======')
            print('Four design criteria log10() value:')
            print('A-optimality:', np.log10(np.trace(totalFIM)))
            print('D-optimality:', np.log10(np.linalg.det(totalFIM)))
            print('E-optimality:', np.log10(np.linalg.eigvals(totalFIM)))
            print('Modified E-optimality:', np.log10(np.linalg.cond(totalFIM)))
            
        return totalFIM

        
    def doe(self, exp_idx, scale=False):
        '''
        param_file: csv file including all parameter values 
        exp_idx: an integer, indicating which line in data_exp is the design vector for this model
        '''
        
        
        createmod = self.create_model_object.create_model(self.data_exp.iloc[exp_idx])
        
        # Control time set [h]
        t_control = [0]

        # design variable and its control time set
        # (What goes here does not matter because we tell Pyomo.DOE that it does not need to fix design variables)
        dv_pass = {'fs.F101.inlet.temperature': t_control,
                   'fs.F101.inlet.pressure': t_control,
                 'fs.F101.inlet.mole_frac_comp[0,"R32"]': t_control,
                  "fs.F101.inlet.mole_frac_comp[0,'emimTf2N']":t_control}

        # Create measurement object
        #measure_pass = {'fs.F101.control_volume.properties_out[0.0].pressure': t_control}
        measure_pass = {'fs.state_block.pressure': t_control}
        measure_class =  Measurements(measure_pass)


        # prior information: none
        prior_pass = [[0]]
        
        # pyomo doe mode option
        sensi_opt = 'direct_kaug'
    
        # Define experiments and design variable names 
        # This does not matter for this problem but needed for inputs
        exp1 = {'fs.F101.inlet.temperature': {0: 1},
                   'fs.F101.inlet.pressure': {0: 1},
                 'fs.F101.inlet.mole_frac_comp[0,"R32"]': {0: 1},
                  "fs.F101.inlet.mole_frac_comp[0,'emimTf2N']": {0: 1}}

        design_names = ['model.fs.F101.inlet.temperature[0]', 'model.fs.F101.inlet.pressure[0]',
                        'model.fs.F101.inlet.mole_frac_comp[0,"R32"]', 'fs.F101.inlet.mole_frac_comp[0,"emimTf2N"]']

        
        
        # create object
        doe_object = DesignOfExperiments(self.param_dict, dv_pass,
                                         measure_class, createmod,
                                        prior_FIM=prior_pass, fixed=True)


        # compute FIM for a square MBDOE problem
        # Note that I did not scale the Jacobian 
        result = doe_object.compute_FIM(exp1, mode=sensi_opt, FIM_store_name = 'thermo.csv', 
                                        scale_nominal_param_value=scale, 
                                        store_output = 'store_output', read_output=None,
                                        formula='central')

        # calculate FIM 
        result.calculate_FIM(doe_object.design_values)
        
        
        if self.verbose:
            print('======Result summary======')
            print('Four design criteria log10() value:')
            print('A-optimality:', np.log10(result.trace))
            print('D-optimality:', np.log10(result.det))
            print('E-optimality:', np.log10(result.min_eig))
            print('Modified E-optimality:', np.log10(result.cond))
            
            
        return result



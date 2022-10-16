import matplotlib.pyplot as plt
import numpy as np
import pandas as pd 
import idaes
import scipy.stats as stats
import json

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
 

from generalize_functions import CEOSModels


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
        
        
    def sumDOE(self, exp_idx_set, scale_opt=False, record_name=None,  init_temp_option = 283.1, init_pressure_option = 399300 , init_x_c1_option = 0.45, poly_option=False):
        
        num_param = len(self.param_name)
        
        totalFIM = [[0]*num_param for i in range(num_param)]
        
        failed_set = []
        
        record = {}
        record["model_name"] = record_name
        
        for i in exp_idx_set:
            print("==========Experiment index:", i, "===============")
            
            res = self.doe(i, scale=scale_opt, init_temp_opt = init_temp_option, 
                          init_pressure_opt = init_pressure_option, init_x_c1_opt = init_x_c1_option, poly_opt = poly_option)
            
            record[str(i)] = res.FIM.tolist()
            
            totalFIM += res.FIM
            #except:
                #failed_set.append(i)
                #print("Failure initialization!")
        
        
        # record
        record["Total"] = totalFIM.tolist()
        record = [record]
        print(record)
        file_name = json.dumps(record)
        f1 = open('./emimtf2n_FIM_info/'+record_name+".json", 'w')
        f1.write(file_name)
        f1.close()
            
        print("Failed set:", failed_set)
        if self.verbose: 
            print('======Result summary======')
            print('Four design criteria log10() value:')
            print('A-optimality:', np.log10(np.trace(totalFIM)))
            print('D-optimality:', np.log10(np.linalg.det(totalFIM)))
            print('E-optimality:', np.log10(min(np.linalg.eigvals(totalFIM))))
            print('Modified E-optimality:', np.log10(np.linalg.cond(totalFIM)))
            
        return totalFIM

        
    def doe(self, exp_idx, scale=False, init_temp_opt = 283.1, init_pressure_opt = 399300 , init_x_c1_opt = 0.45, poly_opt=False):
        '''
        param_file: csv file including all parameter values 
        exp_idx: an integer, indicating which line in data_exp is the design vector for this model
        ''' 
        
        createmod = self.create_model_object.create_model(self.data_exp.iloc[exp_idx], init_temp=init_temp_opt,
                                                         init_pressure = init_pressure_opt, init_x_c1 = init_x_c1_opt, polynomial=poly_opt)
        
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
    
    
    def run_grid_search(self, design_range_values):
        
        createmod = self.create_model_object.create_model(self.data_exp.iloc[exp_idx], init_temp=init_temp_opt,
                                                         init_pressure = init_pressure_opt, init_x_c1 = init_x_c1_opt, polynomial=poly_opt)
        
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
        
        
        # specify inputs 
        # Design variable ranges as lists 
        design_ranges = [list(np.linspace(1,5,5)), list(np.linspace(300,700,5))]

        # Design variable names 
        dv_apply_name = ['CA0','T']

        # Design variable should be fixed at these time points
        dv_apply_time = [[0],t_control]

        # Define experiments. This is a starting point of which the value does not matter
        exp1 = generate_exp(t_control, 5, [300, 300, 300, 300, 300, 300, 300, 300, 300])
        
        # add prior information
        prior_all = [[ 22.52943024 , 1.84034314, -70.23273336, -11.09432962],
         [   1.84034314 ,  18.09848116 ,  -5.73565034 , -109.15866135],
         [ -70.23273336 ,  -5.73565034 , 218.94192843 ,  34.57680848],
         [ -11.09432962 , -109.15866135 ,  34.57680848 ,  658.37644634]]

        print(np.shape(prior_all))

        prior_pass=np.asarray(prior_all)
        print(np.shape(prior_pass))

        print('The prior information FIM:', prior_pass)
        print('Prior Det:', np.linalg.det(prior_pass))
        
        result = doe_object.run_grid_search(design_values, design_ranges, design_dimension_names, design_control_time, mode='direct_kaug',fixed_model = True, tee_option=False, scale_nominal_param_value=False, scale_constant_value=1, store_name= None, read_name=None,filename=None, formula='central', step=0.001)



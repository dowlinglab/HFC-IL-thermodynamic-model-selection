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
        
        record_pressure = {}
        record_pressure["model_name"] = record_name+"_pressure"
        
        for i in exp_idx_set:
            print("==========Experiment index:", i, "===============")
            
            
            
            try: 
                res = self.doe(i, scale=scale_opt, init_temp_opt = init_temp_option, 
                              init_pressure_opt = init_pressure_option, init_x_c1_opt = init_x_c1_option, poly_opt = poly_option)

                record[str(i)] = res.FIM.tolist()

                record_pressure[str(i)] = res.pressure

                totalFIM += res.FIM
            except:
                failed_set.append(i)
                print("Failure initialization!")
        
        print("Failed set:", failed_set)
        # record
        record["Total"] = totalFIM.tolist()
        record = [record]
        #print(record)
        file_name = json.dumps(record)
        f1 = open('./32_bmimpf6_FIM_info/'+record_name+".json", 'w')
        f1.write(file_name)
        f1.close()
        
        # record pressure
        record_pressure = [record_pressure]
        file_name2 = json.dumps(record_pressure)
        f2 = open('./32_bmimpf6_FIM_info/'+record_name+"_pressure.json", 'w')
        f2.write(file_name2)
        f2.close()
            
        
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
                  # 'fs.F101.inlet.mole_frac_comp[0,"R125"]': t_control, 
                  # "fs.F101.inlet.mole_frac_comp[0,'emimTf2N']":t_control}
                    "fs.F101.inlet.mole_frac_comp[0,'bmimpf6']":t_control}
                   
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
                #'fs.F101.inlet.mole_frac_comp[0,"R125"]': {0: 1},   
                #"fs.F101.inlet.mole_frac_comp[0,'emimTf2N']": {0: 1}}
                "fs.F101.inlet.mole_frac_comp[0,'bmimpf6']": {0: 1}}

        design_names = ['model.fs.F101.inlet.temperature[0]', 'model.fs.F101.inlet.pressure[0]',
                        'model.fs.F101.inlet.mole_frac_comp[0,"R32"]', 
                        #'model.fs.F101.inlet.mole_frac_comp[0,"R125"]',
                        #'fs.F101.inlet.mole_frac_comp[0,"emimTf2N"]']
                        'fs.F101.inlet.mole_frac_comp[0,"bmimpf6"]']

        
        
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
    
    
    def run_grid_search(self, design_range_values, fixed_models=None, prior_FIM=None, scale_opt=False,
                        scale_pressure = 1, store_name="record"):
        
        #createmod = self.create_model_object.create_model(self.data_exp.iloc[exp_idx], init_temp=init_temp_opt,
        #                                                 init_pressure = init_pressure_opt, init_x_c1 = init_x_c1_opt, polynomial=poly_opt)
        
        # here this createmod only helps create the object; it is not the real model used for computation.
        createmod = fixed_models[0]
        
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
        
        
        num_param = len(self.param_name)
        
        if not prior_FIM:
            prior_FIM = [[0]*num_param for i in range(num_param)]
        
        prior_pass = np.asarray(prior_FIM)
        
        
        # create object
        doe_object = DesignOfExperiments(self.param_dict, dv_pass,
                                         measure_class, createmod,
                                        prior_FIM=prior_pass, fixed=True)
        
        
        # specify inputs 
        # Design variable ranges as lists 
        design_ranges = design_range_values

        # Design variable names 
        #dv_apply_name = ['Temperature', 'x_R32']
        dv_apply_name = ['Temperature', 'x_R125']
        # Design variable should be fixed at these time points
        dv_apply_time = [[0],[0]]

        # fixed_model_list
        all_fim = doe_object.run_grid_search(exp1, design_ranges, dv_apply_name, dv_apply_time,record_name=store_name, mode='direct_kaug',fixed_model_list = fixed_models, tee_option=False, scale_nominal_param_value=scale_opt,
                                            scale_constant_value=scale_pressure)

        
        return all_fim
        # draw heatmap
        #fixed = {}
        #all_fim.figure_drawing(fixed, ['CA0','T'], 'Reactor case','$C_{A0}$ [M]', 'T [K]' )

# Data Science for Thermodynamic Modeling: Case Study for Ionic Liquid and Hydrofluorocarbon Refrigerant Mixtures

## Data

HFC-32/[emim][TF2N]: modsel/paramest/emimtf2n/R32/Final_Results/r32_emimtf2n_subset.csv

HFC-125/[emim][TF2N]: modsel/paramest/emimtf2n/R125/Final_Results/r125_emimtf2n_subset.csv

HFC-32/[bmim][PF6]: modsel/paramest/bmimpf6/R32/Final_Results/r32_bmimpf6_subset.csv

## Fitting and Analysis File Locations

### HFC-32/[emim][TF2N] System

Position: modsel/paramest/emimtf2n/R32/Final_Results/

File List:
 - PR_1param_Opt1_final.ipynb: Fitting and anlysis notebook for model PR-1A
 - PR_1param_Opt2_final.ipynb: Fitting and anlysis notebook for model PR-1B
 - PR_3params_Opt1_final.ipynb: Fitting and anlysis notebook for model PR-3A
 - PR_3params_Opt2_final.ipynb: Fitting and anlysis notebook for model PR-3B
 - PR_noTdep_final.ipynb: Fitting and anlysis notebook for model PR-2
 - PR_linTdep_final.ipynb: Fitting and anlysis notebook for model PR-4
 - PR_quadTdep_final.ipynb: Fitting and anlysis notebook for model PR-6
 - PR_polyTdep_final.ipynb: Fitting and anlysis notebook for model PR-8
 - SRK_1param_Opt1_final.ipynb: Fitting and anlysis notebook for model SRK-1A
 - SRK_1param_Opt2_final.ipynb: Fitting and anlysis notebook for model SRK-1B
 - SRK_3params_Opt1_final.ipynb: Fitting and anlysis notebook for model SRK-3A
 - SRK_3params_Opt2_final.ipynb: Fitting and anlysis notebook for model SRK-3B
 - SRK_noTdep_final.ipynb: Fitting and anlysis notebook for model SRK-2
 - SRK_linTdep_final.ipynb: Fitting and anlysis notebook for model SRK-4
 - SRK_quadTdep_final.ipynb: Fitting and anlysis notebook for model SRK-6
 - SRK_polyTdep_final.ipynb: Fitting and anlysis notebook for model SRK-8

Initialization for each model: Data/Init_Final

Fits for each model: Data/Fits

AIC analysis for each model: Data/AIC

Covariance matrices for each model (these covariances were calculated using parmest and were not included in the paper): Data/Covariance

### HFC-125/[emim][TF2N] System

Position: modsel/paramest/emimtf2n/R125/Final_Results/

Same file structure as HFC-32/[emim][TF2N] system.

### HFC-32/[bmim][PF6] System

Position: modsel/paramest/bmimpf6/R32/Final_Results/

Same file structure as HFC-32/[emim][TF2N] system.

## MBDoE Analysis Instructions 

Position: modsel/paramest

File list:
 - fim_doe.py: A modified version of Pyomo.DoE for this project.
 - generalize_functions.py: Generalize all candidate models for easier switch. 
 - MBDOE.py: User interface for Pyomo.DoE for this project. 
 - Apply_MBDOE_to_models.ipynb: Model simulations and MBDoE analysis.
 - draw_heatmap.ipynb: MBDoE data analysis
 - /emimtf2n/R32/Final_Results/MBDoE/ModelDisc/: Model discrimination analysis for the HFC-32/emimTF2N system
 - /emimtf2n/R32/Final_Results/MBDoE/InfoContent/: Information content analysis for the HFC-32/emimTF2N system
 
Data list:
 - emimtf2n_FIM_info/: folder contains HFC-R32 (emimtf2n) system FIM prior information and heatmap solution
 - 125_emimtf2n_FIM_info/: folder contains HFC-R125 (emimtf2n) system FIM prior information 
 - 32_bmimpf6_FIM_info/: folder contains HFC-R32 (bmimpf6) system FIM prior information


 ## Number of experiments based on D-optimality

Position: modsel/paramest

File list:
- Apply_MBDOE_to_num_experiments.ipynb

Data list:
- PR_quadTdep_scaleopt.json
 
 
 
 
 
 
 
 
 

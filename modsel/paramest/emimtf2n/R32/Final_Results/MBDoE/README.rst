'''
Model Selection: MBDoE 

Contributors: Bridgette Befort, Alejandro Garciadiego, Jialu Wang, Ke Wang, Alexander Dowling
University of Notre Dame (2022)
'''

This zip file contains:

17 jupyter notebooks (1 for each model)
Params (folder): contains optimized parameters for each model
configuration files:
    hfc32_emimtf2n_PR.py
    hfc32_emimtf2n_PR_polynomial.py
    hfc32_emimtf2n_SRK.py
    hfc32_emimtf2n_SRK_polynomial.py
data file:
    r32_emimtf2n_subset.csv
idaes eos files:
    ceos_k_tempdep.py
    ceos_k_tempdep_polynomial.py
    
A notebook:
    -loads in data
    -loads in the optimized parameters
    -loads in the configuration
    -applies the model and optimized parameters to the individual data points
    

To run a notebook:

- add the idaes eos files to the following path in your idaes workflow: idaes-pse/idaes/generic_models/properties/core/eos
- open the notebooks, hit runall

To update the data:

-make a new csv of data, load in, runall

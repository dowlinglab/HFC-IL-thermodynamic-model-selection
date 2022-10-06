Instructions for Using IDAES to Estimate Parameters for the Two-Parameter PR EoS and Predict New Pressures

Bridgette Befort, Alejandro Garciadiego, Edward Maginn, Alexander Dowling
University of Notre Dame (2022)

This provides general instructions for using IDAES to estimate parameters and then use them to predict pressure using the two-parameter PR EOS for a given HFC/IL system. This readme contains information on the files needed and procedure. Within this document is an example for the R-32/emimTF2N system.

Required Files to Run:

-bip_fitting_functions.py
    -sets up the parameter estimation
    
-csv file of data containing temperature (K), pressure (Pa), x_HFC, x_IL
    -fitting won't work on x_HFC=0 (or near 0 values)

-configuration file of system of interest

-.ipynb for the parameter estimation and solubility prediction


Procedure:

Step 0: Set up IDAES

Download the IDAES development package following these instructions: https://idaes-pse.readthedocs.io/en/stable/tutorials/advanced_install/index.html

Conda activate your IDAES environment

Open a jupyter notebook

Step 1: Confirm necessary files are created/available

-fill out the configuration file with information for the given system of interest (the information to fillin is noted in the example configuration file)

-fill out the .ipynb with the information for the given system of interest (the information to fillin is noted in the example .ipynb)

Step 2: Run the .ipynb

-you may need to mess around with initial parameters to get a good fit

-you will need to update the predicted points with desired initial data and guesses
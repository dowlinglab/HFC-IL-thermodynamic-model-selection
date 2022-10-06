'''
HFC-32/emimTF2N Configuration File

Bridgette Befort, Alejandro Garciadiego, Edward Maginn, Alexander Dowling

University of Notre Dame, 2022
'''

#!/usr/bin/env python
# coding: utf-8

# Import Python libraries
import logging

# Import Pyomo units
from pyomo.environ import units as pyunits

# Import IDAES cores
from idaes.core import LiquidPhase, VaporPhase, Component
from idaes.core.phases import PhaseType as PT
from idaes.generic_models.properties.core.state_definitions import FTPx
from idaes.generic_models.properties.core.eos.ceos_k_tempdep import Cubic, CubicType
from idaes.generic_models.properties.core.phase_equil import SmoothVLE
from idaes.generic_models.properties.core.phase_equil.bubble_dew import LogBubbleDew
from idaes.generic_models.properties.core.phase_equil.forms import log_fugacity


from idaes.generic_models.properties.core.pure import RPP4
from idaes.generic_models.properties.core.pure import NIST

# Set up logger
_log = logging.getLogger(__name__)


configuration = {
    # Specifying components
    "components": {
        "emimTf2N": {"type": Component, #FILLIN: change name from 'emimTf2N' to desired IL
                    "enth_mol_ig_comp": RPP4,
                    "entr_mol_ig_comp": RPP4,
                    "valid_phase_types": PT.liquidPhase,
                    "phase_equilibrium_form": {("Vap", "Liq"): log_fugacity},
                    "parameter_data": {
                        "mw": (291.0E-3, pyunits.kg/pyunits.mol),  #FILLIN: molecular weight
                        "pressure_crit": (2.92E6, pyunits.Pa),  #FILLIN: critical pressure
                        "temperature_crit": (906.91, pyunits.K),  #FILLIN: critical temperature
                        "omega": 0.4223,  #FILLIN: acentric factor
                        "cp_mol_ig_comp_coeff": {
                            'A': (3259.5745, pyunits.J/pyunits.mol/pyunits.K),  
                            'B': (-28.5610, pyunits.J/pyunits.mol/pyunits.K**2),
                            'C': (0.09354, pyunits.J/pyunits.mol/pyunits.K**3),
                            'D': (-0.0001000673, pyunits.J/pyunits.mol/pyunits.K**4)},
                    "enth_mol_form_vap_comp_ref": (
                        0.145, pyunits.J/pyunits.mol),  
                    "entr_mol_form_vap_comp_ref": (
                        137.5, pyunits.J/pyunits.mol/pyunits.K)}}, 
        "R32": {"type": Component, #FILLIN: change name from 'R32" to desired HFC
                  "enth_mol_ig_comp": RPP4,
                  "entr_mol_ig_comp": RPP4,
                  "pressure_sat_comp": NIST,
                  "phase_equilibrium_form": {("Vap", "Liq"): log_fugacity},
                  "parameter_data": {
                      "mw": (52.02E-3, pyunits.kg/pyunits.mol),  #FILLIN: molecular weight
                      "pressure_crit": (57.85e5, pyunits.Pa),  #FILLIN: critical pressure
                      "temperature_crit": (351.3, pyunits.K),  #FILLIN: critical temperature
                      "omega": 0.277,  #FILLIN: acentric factor
                      "cp_mol_ig_comp_coeff": {  #FILLIN: heat capacity coefficients
                          "A": (4785.106, pyunits.J/pyunits.mol/pyunits.K),
                          "B": (-68.0944, pyunits.J/pyunits.mol/pyunits.K**2),
                          "C": (0.3278, pyunits.J/pyunits.mol/pyunits.K**3),
                          "D": (-0.000524, pyunits.J/pyunits.mol/pyunits.K**4)},
                  "enth_mol_form_vap_comp_ref": (
                      21.2e3, pyunits.J/pyunits.mol),  #FILLIN: enthalpy of formation
                  "entr_mol_form_vap_comp_ref": (
                      -269, pyunits.J/pyunits.mol/pyunits.K),  #FILLIN: entropy of formation reference
                  "pressure_sat_comp_coeff": {"A": (4.26224, None),  #FILLIN: Saturation pressure coefficients for Antoine
                                              "B": (821.092, pyunits.K),
                                              "C": (-28.554, pyunits.K)}}}},
    # Specifying phases
    "phases":  {'Liq': {"type": LiquidPhase,
                        "equation_of_state": Cubic,
                        "equation_of_state_options": {
                            "type": CubicType.PR}},
                'Vap': {"type": VaporPhase,
                        "equation_of_state": Cubic,
                        "equation_of_state_options": {
                            "type": CubicType.PR}}},

    # Set base units of measurement
    "base_units": {"time": pyunits.s,
                   "length": pyunits.m,
                   "mass": pyunits.kg,
                   "amount": pyunits.mol,
                   "temperature": pyunits.K},

    # Specifying state definition
    "state_definition": FTPx,
    "state_bounds": {"flow_mol": (0, 100, 1000, pyunits.mol/pyunits.s),
                     "temperature": (10, 300, 500, pyunits.K),
                     "pressure": (5e-4, 1e5, 1e10, pyunits.Pa)},
    "pressure_ref": (101325, pyunits.Pa),
    "temperature_ref": (298.15, pyunits.K),

    # Defining phase equilibria
    "phases_in_equilibrium": [("Vap", "Liq")],
    "phase_equilibrium_state": {("Vap", "Liq"): SmoothVLE},
    "bubble_dew_method": LogBubbleDew,
     #FILLIN: change components in parameter names to desired HFC and IL
    "parameter_data": {"PR_kappa_A": {("emimTf2N", "emimTf2N"): 0.000, 
                                    ("emimTf2N", "R32"): 0.000,
                                    ("R32", "R32"): 0.000,
                                    ("R32", "emimTf2N"): 0.000},
                     "PR_kappa_B": {("emimTf2N", "emimTf2N"): 0.000,
                                    ("emimTf2N", "R32"): 0.0,
                                    ("R32", "R32"): 0.000,
                                    ("R32", "emimTf2N"): 0.0},
                     "PR_kappa_C": {("emimTf2N", "emimTf2N"): 0.000,
                                    ("emimTf2N", "R32"): 0.0,
                                    ("R32", "R32"): 0.000,
                                    ("R32", "emimTf2N"): 0.0}}}
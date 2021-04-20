#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# Import Python libraries
import logging

# Import Pyomo units
from pyomo.environ import units as pyunits

# Import IDAES cores
from idaes.core import LiquidPhase, VaporPhase, Component
from idaes.core.phases import PhaseType as PT
from idaes.generic_models.properties.core.state_definitions import FTPx
from idaes.generic_models.properties.core.eos.ceos import Cubic, CubicType
from idaes.generic_models.properties.core.phase_equil import smooth_VLE
from idaes.generic_models.properties.core.phase_equil.bubble_dew import         LogBubbleDew
from idaes.generic_models.properties.core.phase_equil.forms import log_fugacity


import idaes.generic_models.properties.core.pure.RPP4 as RPP
import idaes.generic_models.properties.core.pure.NIST as NIST

# Set up logger
_log = logging.getLogger(__name__)


# ---------------------------------------------------------------------
# Configuration dictionary for an ideal Benzene-Toluene system

# Data Sources:

# [1] Critical Properties, Normal Boiling Temperature, and Acentric Factor of Another 200 Ionic Liquids.
# José O. Valderrama, Wilson W. Sanga, and Juan A. Lazzús. Industrial & Engineering Chemistry Research 2008 47 (4),
# 1318-1330. Industrial & Engineering Chemistry Research 2008 47 (4), 1318-1330. DOI: 10.1021/ie071055d

# [2] The Properties of Gases and Liquids (1987) 4th edition, Chemical Engineering Series - Robert C. Reid

# [4] Engineering Toolbox, https://www.engineeringtoolbox.com
#     Retrieved 1st December, 2019

# ---------------------------------------------------------------------


configuration = {
    # Specifying components
    "components": {
        "hmimFAP": {"type": Component,
                    "enth_mol_ig_comp": RPP,
                    "entr_mol_ig_comp": RPP,
                    "valid_phase_types": PT.liquidPhase,
                    "phase_equilibrium_form": {("Vap", "Liq"): log_fugacity},
                    "parameter_data": {
                        "mw": (6.12E-3, pyunits.kg/pyunits.mol),  # [1]
                        "pressure_crit": (2.7E06, pyunits.Pa),  # [1]
                        "temperature_crit": (800, pyunits.K),  # [1]
                        "omega": 0.600,  # [1]
                        "cp_mol_ig_comp_coeff": {
                            'A': (3259.5745, pyunits.J/pyunits.mol/pyunits.K), 
                            'B': (-28.5610, pyunits.J/pyunits.mol/pyunits.K**2),
                            'C': (0.09354, pyunits.J/pyunits.mol/pyunits.K**3),
                            'D': (-0.0001000673, pyunits.J/pyunits.mol/pyunits.K**4)},
                        "enth_mol_form_vap_comp_ref": (
                            0.145, pyunits.J/pyunits.mol),  
                        "entr_mol_form_vap_comp_ref": (
                            137.5, pyunits.J/pyunits.mol/pyunits.K)}},
        "R125": {"type": Component,
                  "enth_mol_ig_comp": RPP,
                  "entr_mol_ig_comp": RPP,
                  "pressure_sat_comp": NIST,
                  "phase_equilibrium_form": {("Vap", "Liq"): log_fugacity},
                  "parameter_data": {
                      "mw": (120.002E-3, pyunits.kg/pyunits.mol),  # [1]
                      "pressure_crit": (36.29e5, pyunits.Pa),  # [1]
                      "temperature_crit": (339.160, pyunits.K),  # [1]
                      "omega": 0.305,  # [1]
                      "cp_mol_ig_comp_coeff": {
                          "A": (38.7614, pyunits.J/pyunits.mol/pyunits.K),  # [1]
                          "B": (0.78561, pyunits.J/pyunits.mol/pyunits.K**2),
                          "C": (-0.00468, pyunits.J/pyunits.mol/pyunits.K**3),
                          "D": (9.09976, pyunits.J/pyunits.mol/pyunits.K**4)},
                      "enth_mol_form_vap_comp_ref": (
                          22.8e3, pyunits.J/pyunits.mol),  # [ ]
                      "entr_mol_form_vap_comp_ref": (
                          -269, pyunits.J/pyunits.mol/pyunits.K),  # [ ]
                      "pressure_sat_comp_coeff": {"A": (4.7694, None),  # [ ]
                                                  "B": (1153.2136, pyunits.K),
                                                  "C": (19.5064, pyunits.K)}}}},
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
    "phase_equilibrium_state": {("Vap", "Liq"): smooth_VLE},
    "bubble_dew_method": LogBubbleDew,
    "parameter_data": {"PR_kappa": {("hmimFAP", "hmimFAP"): 0.000,
                                    ("hmimFAP", "R125"): -0.01408977877696036,
                                    ("R125", "R125"): 0.000,
                                    ("R125", "hmimFAP"): -0.018059739008510666}}}


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
# [1] Separation of CO2 and H2S using room-temperature ionic liquid [bmim][PF6] Mark B.Shiflett, A.Yokozeki, 2010
# [2] The Properties of Gases and Liquids (1987) 4th edition, Chemical Engineering Series - Robert C. Reid
# [3] Critical Properties, Normal Boiling Temperatures, and Acentric Factors of Fifty Ionic Liquids.
#     J. O. Valderrama and P. A. Robles. Industrial & Engineering Chemistry Research 2007 46 (4), 1338-1344.
#     DOI: 10.1021/ie0603058
# [4] Engineering Toolbox, https://www.engineeringtoolbox.com
#     Retrieved 1st December, 2019

# ---------------------------------------------------------------------


configuration = {
    # Specifying components
    "components": {
        "bmimPF6": {"type": Component,
                    "enth_mol_ig_comp": RPP,
                    "entr_mol_ig_comp": RPP,
                    "valid_phase_types": PT.liquidPhase,
                    "phase_equilibrium_form": {("Vap", "Liq"): log_fugacity},
                    "parameter_data": {
                        "mw": (284.0E-3, pyunits.kg/pyunits.mol),  # [3]
                        "pressure_crit": (2.40E06, pyunits.Pa),  # [1]
                        "temperature_crit": (860, pyunits.K),  # [3]
                        "omega": 0.7917,  # [3]
                        "cp_mol_ig_comp_coeff": {
                            'A': (3259.5745, pyunits.J/pyunits.mol/pyunits.K), 
                            'B': (-28.5610, pyunits.J/pyunits.mol/pyunits.K**2),
                            'C': (0.0934, pyunits.J/pyunits.mol/pyunits.K**3),
                            'D': (-0.0001000673, pyunits.J/pyunits.mol/pyunits.K**4)},
                        "enth_mol_form_vap_comp_ref": (
                            0.145, pyunits.J/pyunits.mol),  # [3] for bmimPF6
                        "entr_mol_form_vap_comp_ref": (
                            137.5, pyunits.J/pyunits.mol/pyunits.K)}},
        "R32": {"type": Component,
                  "enth_mol_ig_comp": RPP,
                  "entr_mol_ig_comp": RPP,
                  "pressure_sat_comp": NIST,
                  "phase_equilibrium_form": {("Vap", "Liq"): log_fugacity},
                  "parameter_data": {
                      "mw": (52.02E-3, pyunits.kg/pyunits.mol),  # [2]
                      "pressure_crit": (57.85e5, pyunits.Pa),  # [2]
                      "temperature_crit": (351.3, pyunits.K),  # [2]
                      "omega": 0.277,  # [1]
                      "cp_mol_ig_comp_coeff": {
                          "A": (4785.106, pyunits.J/pyunits.mol/pyunits.K),  # [2]
                          "B": (-68.0944, pyunits.J/pyunits.mol/pyunits.K**2),
                          "C": (0.3278, pyunits.J/pyunits.mol/pyunits.K**3),
                          "D": (-0.000524, pyunits.J/pyunits.mol/pyunits.K**4)},
                      "enth_mol_form_vap_comp_ref": (
                          21.2e3, pyunits.J/pyunits.mol),  # [4]
                      "entr_mol_form_vap_comp_ref": (
                          -269, pyunits.J/pyunits.mol/pyunits.K),  # [4]
                      "pressure_sat_comp_coeff": {"A": (4.26224, None),  # [4]
                                                  "B": (821.092, pyunits.K),
                                                  "C": (-28.554, pyunits.K)}}}},
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
    "parameter_data": {"PR_kappa": {("bmimPF6", "bmimPF6"): 0.000,
                                    ("bmimPF6", "R32"): -0.01408977877696036,
                                    ("R32", "R32"): 0.000,
                                    ("R32", "bmimPF6"): -0.018059739008510666}}}





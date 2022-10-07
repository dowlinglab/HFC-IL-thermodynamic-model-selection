#!/usr/bin/env python
# coding: utf-8

# In[2]:


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


# ---------------------------------------------------------------------
# Configuration dictionary for an ideal Benzene-Toluene system

# Data Sources:
# [5] Mohammad Shokouhi, Amir Hossein Saali, Mehdi Vahidi, Ali Taghi Zoghi, Amir Hossein Jalili,
#     Diffusivity and solubility of carbonyl sulfide and sulfur dioxide in 1-ethyl-3-methylimidazolium
#     bis (trifluoromethyl) sulfonylimide ([emim][Tf2N]): Experimental measurement and modelling, The
#     Journal of Chemical Thermodynamics,Vol 132, 2019, Pages 411-422, SSN 0021-9614,
#     https://doi.org/10.1016/j.jct.2019.01.019. (https://www.sciencedirect.com/science/article/pii/S0021961419300412)
# [1] The Properties of Gases and Liquids (1987)
#     4th edition, Chemical Engineering Series - Robert C. Reid
# [3] Engineering Toolbox, https://www.engineeringtoolbox.com
#     Retrieved 1st December, 2019
# ---------------------------------------------------------------------


configuration = {
    # Specifying components
    "components": {
        "emimTf2N": {"type": Component,
                    "enth_mol_ig_comp": RPP4,
                    "entr_mol_ig_comp": RPP4,
                    "valid_phase_types": PT.liquidPhase,
                    "phase_equilibrium_form": {("Vap", "Liq"): log_fugacity},
                    "parameter_data": {
                        "mw": (291.0E-3, pyunits.kg/pyunits.mol),  # [5]
                        "pressure_crit": (2.92E6, pyunits.Pa),  # [5]
                        "temperature_crit": (906.91, pyunits.K),  # [5]
                        "omega": 0.4223,  # 51]
                        "cp_mol_ig_comp_coeff": {
                            'A': (3259.5745, pyunits.J/pyunits.mol/pyunits.K),  # [1] for bmimPF6
                            'B': (-28.5610, pyunits.J/pyunits.mol/pyunits.K**2),
                            'C': (0.09354, pyunits.J/pyunits.mol/pyunits.K**3),
                            'D': (-0.0001000673, pyunits.J/pyunits.mol/pyunits.K**4)},
                        "enth_mol_form_vap_comp_ref": (
                            0.145, pyunits.J/pyunits.mol),  # [3] for bmimPF6
                        "entr_mol_form_vap_comp_ref": (
                            137.5, pyunits.J/pyunits.mol/pyunits.K)}},
        "R125": {"type": Component,
                  "enth_mol_ig_comp": RPP4,
                  "entr_mol_ig_comp": RPP4,
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
                          22.8e3, pyunits.J/pyunits.mol),  # [2]
                      "entr_mol_form_vap_comp_ref": (
                          -269, pyunits.J/pyunits.mol/pyunits.K),  # [2]
                      "pressure_sat_comp_coeff": {"A": (4.7694, None),  # [2]
                                                  "B": (1153.2136, pyunits.K),
                                                  "C": (19.5064, pyunits.K)}}}},
    # Specifying phases
    "phases":  {'Liq': {"type": LiquidPhase,
                        "equation_of_state": Cubic,
                        "equation_of_state_options": {
                            "type": CubicType.SRK}},
                'Vap': {"type": VaporPhase,
                        "equation_of_state": Cubic,
                        "equation_of_state_options": {
                            "type": CubicType.SRK}}},

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
    "parameter_data": {"SRK_kappa_A": {("emimTf2N", "emimTf2N"): 0.000,
                                    ("emimTf2N", "R125"): -0.200,
                                    ("R125", "R125"): 0.000,
                                    ("R125", "emimTf2N"): -0.05619},
                     "SRK_kappa_B": {("emimTf2N", "emimTf2N"): 0.000,
                                    ("emimTf2N", "R125"): 0.0,
                                    ("R125", "R125"): 0.000,
                                    ("R125", "emimTf2N"): 0.0},
                     "SRK_kappa_C": {("emimTf2N", "emimTf2N"): 0.000,
                                    ("emimTf2N", "R125"): 0.0,
                                    ("R125", "R125"): 0.000,
                                    ("R125", "emimTf2N"): 0.0}}}
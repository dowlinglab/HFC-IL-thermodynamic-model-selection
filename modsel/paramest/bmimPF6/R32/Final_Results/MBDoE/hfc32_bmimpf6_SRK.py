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
        "bmimpf6": {"type": Component,
                    "enth_mol_ig_comp": RPP4,
                    "entr_mol_ig_comp": RPP4,
                    "valid_phase_types": PT.liquidPhase,
                    "phase_equilibrium_form": {("Vap", "Liq"): log_fugacity},
                    "parameter_data": {
                        "mw": (284.18E-3, pyunits.kg/pyunits.mol),  # [1]
                        "pressure_crit": (24e5, pyunits.Pa),  # [1]
                        "temperature_crit": (860, pyunits.K),  # [1]
                        "omega": 0.7917,  # [1]
                        "cp_mol_ig_comp_coeff": {
                            'A': (3259.5745, pyunits.J/pyunits.mol/pyunits.K),  # [1]
                            'B': (-28.5610, pyunits.J/pyunits.mol/pyunits.K**2),
                            'C': (0.09354, pyunits.J/pyunits.mol/pyunits.K**3),
                            'D': (-0.0001000673, pyunits.J/pyunits.mol/pyunits.K**4)},
                        "enth_mol_form_vap_comp_ref": (
                            0.145, pyunits.J/pyunits.mol),  # [3]
                        "entr_mol_form_vap_comp_ref": (
                            137.5, pyunits.J/pyunits.mol/pyunits.K)}},
        "R32": {"type": Component,
                  "enth_mol_ig_comp": RPP4,
                  "entr_mol_ig_comp": RPP4,
                  "pressure_sat_comp": NIST,
                  "phase_equilibrium_form": {("Vap", "Liq"): log_fugacity},
                  "parameter_data": {
                      "mw": (52.02E-3, pyunits.kg/pyunits.mol),  # [1]
                      "pressure_crit": (57.85e5, pyunits.Pa),  # [1]
                      "temperature_crit": (351.3, pyunits.K),  # [1]
                      "omega": 0.277,  # [1]
                      "cp_mol_ig_comp_coeff": {
                          "A": (4785.106, pyunits.J/pyunits.mol/pyunits.K),  # [1]
                          "B": (-68.0944, pyunits.J/pyunits.mol/pyunits.K**2),
                          "C": (0.3278, pyunits.J/pyunits.mol/pyunits.K**3),
                          "D": (-0.000524, pyunits.J/pyunits.mol/pyunits.K**4)},
                      "enth_mol_form_vap_comp_ref": (
                          21.2e3, pyunits.J/pyunits.mol),  # [2]
                      "entr_mol_form_vap_comp_ref": (
                          -269, pyunits.J/pyunits.mol/pyunits.K),  # [2]
                      "pressure_sat_comp_coeff": {"A": (4.26224, None),  # [2]
                                                  "B": (821.092, pyunits.K),
                                                  "C": (-28.554, pyunits.K)}}}},
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
    "parameter_data": {"SRK_kappa_A": {("bmimpf6", "bmimpf6"): 0.000,
                                    ("bmimpf6", "R32"): -0.200,
                                    ("R32", "R32"): 0.000,
                                    ("R32", "bmimpf6"): -0.05619},
                     "SRK_kappa_B": {("bmimpf6", "bmimpf6"): 0.000,
                                    ("bmimpf6", "R32"): 0.0,
                                    ("R32", "R32"): 0.000,
                                    ("R32", "bmimpf6"): 0.0},
                     "SRK_kappa_C": {("bmimpf6", "bmimpf6"): 0.000,
                                    ("bmimpf6", "R32"): 0.0,
                                    ("R32", "R32"): 0.000,
                                    ("R32", "bmimpf6"): 0.0}}}
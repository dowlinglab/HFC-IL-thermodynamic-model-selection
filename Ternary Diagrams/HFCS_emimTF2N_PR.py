#!/usr/bin/env python
# coding: utf-8

# In[ ]:


##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2020, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################
"""
Benzene-Toluene phase equilibrium package using ideal liquid and vapor.
Example property package using the Generic Property Package Framework.
This exmample shows how to set up a property package to do benzene-toluene
phase equilibrium in the generic framework using ideal liquid and vapor
assumptions along with methods drawn from the pre-built IDAES property
libraries.
"""
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
# [1] The Properties of Gases and Liquids (1987)
#     4th edition, Chemical Engineering Series - Robert C. Reid
# [3] Engineering Toolbox, https://www.engineeringtoolbox.com
#     Retrieved 1st December, 2019

configuration = {
    # Specifying components
    "components": {
        "emimTf2N": {"type": Component,
                    "enth_mol_ig_comp": RPP,
                    "entr_mol_ig_comp": RPP,
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

        "R32": {"type": Component,
                  "enth_mol_ig_comp": RPP,
                  "entr_mol_ig_comp": RPP,
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
                                                  "C": (-28.554, pyunits.K)}}},

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
                     "temperature": (10, 300, 1000, pyunits.K),
                     "pressure": (5e-4, 1e5, 1e10, pyunits.Pa)},
    "pressure_ref": (101325, pyunits.Pa),
    "temperature_ref": (298.15, pyunits.K),

    # Defining phase equilibria
    "phases_in_equilibrium": [("Vap", "Liq")],
    "phase_equilibrium_state": {("Vap", "Liq"): smooth_VLE},
    "bubble_dew_method": LogBubbleDew,
    "parameter_data": {"PR_kappa": {("R32", "R32"): 0.000,
                                    ("R32", "R125"): 0.003831, #-0.002535,
                                    ("R32", "emimTf2N"): -0.056,
                                    ("R125", "R125"): 0.000,
                                    ("R125", "R32"): 0.00271, #-0.002535,
                                    ("R125", "emimTf2N"): -0.016,
                                    ("emimTf2N", "R32"): -0.2,
                                    ("emimTf2N", "R125"): -0.035,
                                    ("emimTf2N", "emimTf2N"): 0.000}}}


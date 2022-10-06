#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021
# by the software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia University
# Research Corporation, et al.  All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and
# license information.
#################################################################################
"""
Property package for the VLE calucations for a R-32/emimTF2N 
system. If using the activity coefficient models (NRTL or Wilson), the user is
expected to provide the paramters necessary for these models. Please note that
these parameters are declared as variables here to allow for use in a parameter
estimation problem if the VLE data is available.
"""

#See: https://github.com/IDAES/idaes-pse/blob/main/idaes/models/properties/activity_coeff_models/BTX_activity_coeff_VLE.py
#See: https://github.com/IDAES/examples-pse/blob/main/src/Tutorials/Advanced/ParamEst/parameter_estimation_NRTL_using_state_block_solution_testing.ipynb

# Import Pyomo libraries
from pyomo.environ import Param, NonNegativeReals, Set, units as pyunits

# Import IDAES cores
from idaes.core import declare_process_block_class, Component
from idaes.core.util.misc import extract_data

from idaes.generic_models.properties.activity_coeff_models.activity_coeff_prop_pack import (
    ActivityCoeffParameterData,
)
from idaes.logger import getIdaesLogger


# Set up logger
_log = getIdaesLogger(__name__)


@declare_process_block_class("HFCILParameterBlock")
class HFCILParameterData(ActivityCoeffParameterData):
    def build(self):
        """
        Callable method for Block construction.
        """
        self.component_list_master = Set(initialize=["R32","emimTf2N"])

        # Create component objects
        # NOTE: User needs to update this list; can be a subset or
        # equal to the master component list
        self.emimTf2N = Component()
        self.R32 = Component()

        super(HFCILParameterData, self).build()

        # List of phase equilibrium index
        self.phase_equilibrium_idx_master = Set(initialize=[1, 2])

        self.phase_equilibrium_idx = Set(initialize=[1, 2])

        self.phase_equilibrium_list_master = {
            1: ["R32", ("Vap", "Liq")],
            2: ["emimTf2N", "Liq"],
        }

        self.phase_equilibrium_list = {
            1: ["R32", ("Vap", "Liq")],
            2: ["emimTf2N", "Liq"],
        }

        # Thermodynamic reference state
        self.pressure_reference = Param(
            mutable=True,
            default=101325,
            doc="Reference pressure [Pa]",
            units=pyunits.Pa,
        )
        self.temperature_reference = Param(
            mutable=True,
            default=298.15,
            doc="Reference temperature [K]",
            units=pyunits.K,
        )

        pressure_critical_data = {
            "emimTf2N": 2.92E6,
            "R32": 57.85e5
        }

        self.pressure_critical = Param(
            self.component_list,
            within=NonNegativeReals,
            mutable=False,
            initialize=extract_data(pressure_critical_data),
            doc="Critical pressure [Pa]",
            units=pyunits.Pa,
        )

        temperature_critical_data = {
            "emimTf2N": 906.91,
            "R32": 351.3
        }

        self.temperature_critical = Param(
            self.component_list,
            within=NonNegativeReals,
            mutable=False,
            initialize=extract_data(temperature_critical_data),
            doc="Critical temperature [K]",
            units=pyunits.K,
        )

        # Source: The Properties of Gases and Liquids (1987)
        # 4th edition, Chemical Engineering Series - Robert C. Reid
        mw_comp_data = {
            "emimTf2N": 291.0E-3,
            "R32": 52.02E-3
        }

        self.mw_comp = Param(
            self.component_list,
            mutable=False,
            initialize=extract_data(mw_comp_data),
            doc="molecular weight kg/mol",
            units=pyunits.kg / pyunits.mol,
        )
        
        # Source: The Properties of Gases and Liquids (1987)
        # 4th edition, Chemical Engineering Series - Robert C. Reid
        pressure_sat_coeff_data = {
            ("R32", "A"): 4.26224,
            ("R32", "B"): 821.092,
            ("R32", "C"):-28.554,
            ("R32", "D"):0,
            ("emimTf2N", "A"):1,
            ("emimTf2N", "B"): 1,
            ("emimTf2N", "C"): 1,
            ("emimTf2N", "D"): 0,
        }

        self.pressure_sat_coeff = Param(
            self.component_list,
            ["A", "B", "C","D"],
            mutable=False,
            initialize=extract_data(pressure_sat_coeff_data),
            doc="parameters to compute P_sat",
        )


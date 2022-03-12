## @file toroidal_bfield.py
#  @brief containing a problem class with magnetic fields in two cyclical coordinates for pyoculus ODE solver
#  @author Zhisong Qu (zhisong.qu@anu.edu.au)
#

from .toroidal_problem import ToroidalProblem
from .bfield_problem import BfieldProblem
import numpy as np


class ToroidalBfield(ToroidalProblem, BfieldProblem):
    def __init__(self):
        """! Set up the problem with two cyclical coordinates, e.g. \f[ (s, \theta, \zeta)  \f]"""
        super().__init__()

    def f(self, zeta, st, args=None):
        """! Returns ODE RHS
        @param zeta cylindrical angle in ODE
        @param st \f$(s, \theta)\f$ in ODE
        @param arg1 parameter for the ODE
        @returns the RHS of the ODE
        """
        stz = np.array([st[0], st[1], zeta])
        B = self.B(stz, args)
        f = np.array([B[0] / B[2], B[1] / B[2]])
        return f

    def f_tangent(self, zeta, st, args=None):
        """! Returns ODE RHS, with tangent
        @param zeta cylindrical angle in ODE
        @param st \f$(s, \theta, ds_1, d\theta_1, ds_2, d\theta_2)\f$ in ODE
        @param arg1 parameter for the ODE
        @returns the RHS of the ODE, with tangent
        """
        stz = np.array([st[0], st[1], zeta])
        Bu, dBu = self.dBdX(stz, args)

        deltax = np.reshape(st[2:], [2,2])
        gBzeta = Bu[2]

        M = dBu[0:2, 0:2] * gBzeta - dBu[0:2, 2, np.newaxis] * Bu[0:2]

        deltax = deltax @ M / gBzeta**2

        df = np.zeros([6])
        df[0:2] = Bu[0:1] / Bu[2]
        df[2:6] = deltax.flatten()

        return df
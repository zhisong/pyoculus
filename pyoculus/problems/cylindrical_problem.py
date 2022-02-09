## @file cylindrical_problem.py
#  @brief containing cylindrical problem class for pyoculus ODE solver
#  @author Zhisong Qu (zhisong.qu@anu.edu.au)
#

from .base_problem import BaseProblem
import numpy as np


class CylindricalProblem(BaseProblem):
    def __init__(self, R0, Z0, Nfp=1):
        """! Set up cylindrical system \f[ (R, \varphi, Z)  \f]"""
        self.problem_size = 2
        self.poincare_plot_type = "RZ"
        self.poincare_plot_xlabel = "R(m)"
        self.poincare_plot_ylabel = "Z(m)"
        self._R0 = R0
        self._Z0 = Z0
        self.Nfp = Nfp

    def f_RZ(self, phi, RZ, args=None):
        """! Returns ODE RHS
        @param zeta cylindrical angle in ODE
        @param RZ \f$(R, Z)\f$ in ODE
        @param arg1 parameter for the ODE
        @returns the RHS of the ODE
        """
        raise NotImplementedError(
            "A CylindricalProblem class should implement member function f"
        )

    def f_RZ_tangent(self, phi, RZ, args=None):
        """! Returns ODE RHS, with tangent
        @param zeta cylindrical angle in ODE
        @param RZ \f$(R, Z, dR_1, dZ_1, dR_2, dZ_2)\f$ in ODE
        @param arg1 parameter for the ODE
        @returns the RHS of the ODE, with tangent
        """
        raise NotImplementedError(
            "A CylindricalProblem class should implement member function f_tangent"
        )

    def f(self, zeta, y, args=None):
        """! Returns ODE RHS
        @param zeta cylindrical angle in ODE
        @param y \f$(R, Z, R_0, Z_0, \theta)
        @returns the RHS of the ODE
        """
        R = y[0]
        Z = y[1]
        R0 = y[2]
        Z0 = y[3]

        dRZ = self.f_RZ(zeta, np.array([R, Z]))
        dRZ0 = self.f_RZ(zeta, np.array([R0, Z0]))

        dR = dRZ[0]
        dZ = dRZ[1]
        dR0 = dRZ0[0]
        dZ0 = dRZ0[1]

        deltaR = R - R0
        deltaZ = Z - Z0

        dtheta = (deltaR * (dZ - dZ0) - deltaZ * (dR - dR0)) / (
            deltaR ** 2 + deltaZ ** 2
        )

        return np.array([dR, dZ, dR0, dZ0, dtheta])

    def f_tangent(self, zeta, y, args=None):
        """! Returns ODE RHS
        @param zeta cylindrical angle in ODE
        @param y \f$(R, Z, R_0, Z_0, \theta, dR_1, dZ_1, dR_2, dZ_2)
        @returns the RHS of the ODE
        """
        R = y[0]
        Z = y[1]
        R0 = y[2]
        Z0 = y[3]

        dRZ = self.f_RZ_tangent(zeta, np.array([R, Z, y[5], y[6], y[7], y[8]]))
        dRZ0 = self.f_RZ(zeta, np.array([R0, Z0]))

        dR = dRZ[0]
        dZ = dRZ[1]
        dR0 = dRZ0[0]
        dZ0 = dRZ0[1]

        deltaR = R - R0
        deltaZ = Z - Z0

        dtheta = (deltaR * (dZ - dZ0) - deltaZ * (dR - dR0)) / (
            deltaR ** 2 + deltaZ ** 2
        )

        return np.array([dR, dZ, dR0, dZ0, dtheta, dRZ[2], dRZ[3], dRZ[4], dRZ[5]])

    def set_axis(self, R0, Z0):
        """! Set the axis
        @param R0 \f$R_0\f$
        @param Z0 \f$Z_0\f$
        """
        self._R0 = R0
        self._Z0 = Z0

## @file simsgeo_biot_savart.py
#  @brief containing a class for pyoculus ODE solver that deals with Simsgeo B field
#  @author Zhisong Qu (zhisong.qu@anu.edu.au)
#

from pyoculus.problems import CartesianBfield
import numpy as np


class SimsgeoBiotSavart(CartesianBfield):
    def __init__(self, bs, R0, Z0, Nfp=1):
        """! Set up the problem for a simsopt.geo.BiotSavart
        @param bs an instance of simsopt.geo.BiotSavart, use to compute the magnetic field
        """
        from simsopt.geo.biotsavart import BiotSavart

        super().__init__(R0, Z0, Nfp)

        if not isinstance(bs, BiotSavart):
            raise TypeError("bs should be an instance of simsopt.geo.BiotSavart")

        self._bs = bs

    def B(self, xyz, args=None):

        point = [xyz]
        self._bs.set_points(point)
        Bfield=self._bs.B()
        return Bfield[0]

    def dBdX(self, xyz, args=None):

        point = [xyz]
        self._bs.set_points(point)
        dB=self._bs.dB_by_dX()
        return dB[0]
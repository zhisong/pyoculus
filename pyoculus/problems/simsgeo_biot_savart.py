## @file simsgeo_biot_savart.py
#  @brief containing a class for pyoculus ODE solver that deals with Simsgeo B field
#  @author Zhisong Qu (zhisong.qu@anu.edu.au)
#

from pyoculus.problems import CartesianBfield
import numpy as np

class SimsgeoBiotSavart(CartesianBfield):
    def __init__(self, bs, R0, Z0, Nfp=1):
        """! Set up the problem to compute the magnetic field for simsopt.geo.BiotSavart
        @param bs an instance of simsopt.geo.BiotSavart, use to compute the magnetic field
        @param R0 the magnetic axis R coordinate at phi=0 plane
        @param Z0 the magnetic axis Z coordinate at phi=0 plane
        """
        from simsopt.geo.biotsavart import BiotSavart

        super().__init__(R0, Z0, Nfp)

        if not isinstance(bs, BiotSavart):
            raise TypeError("bs should be an instance of simsopt.geo.BiotSavart")

        self._bs = bs

    def B(self, coords, args=None):
        """! The magnetic field, being used by parent class CartesianBfield
        @param coords array with coordinates \f$(x, y z)\f$
        @returns \f$(B_x, B_y, B_z)\f$
        """
        point = [coords]
        self._bs.set_points(point)
        Bfield=self._bs.B()
        return Bfield[0]

    def dBdX(self, coords, args=None):
        """! The derivative of the magnetic field, being used by parent class CartesianBfield
        @param coords array with coordinates \f$(x, y z)\f$
        @returns B, dBdX, B and \f$\partial (B_x, B_y, B_z)/\partial (x,y,z)\f$
        """
        point = [coords]
        self._bs.set_points(point)
        dB=self._bs.dB_by_dX()
        B = self.B(coords)
        return B, dB[0]
## @file qfm_bfield.py
#  @brief Setup the magnetic field system - using the Quadratic minimising surfaces coordinates
#  @author Zhisong Qu (zhisong.qu@anu.edu.au)
#
from .toroidal_bfield import ToroidalBfield
from pyoculus.solvers import SurfacesToroidal
import numpy as np

## Class that used to setup the bfield problem based on QFM
#
class QFMBfield(ToroidalBfield):

    def __init__(self, pb : ToroidalBfield, surfaces : SurfacesToroidal, lvol):
        """! Set up the problems based on known magnetic field and the QFM surfaces
        @param pb  the TorodialBfield problem we will used as the original coordinates
        @param surfaces  the QFM surfaces
        """
        self.pb = pb
        self.surfaces = surfaces
        super().__init__()
    
    def B(self, coords, args=None):
        """! Returns magnetic fields
        @param coordinates \f$(\rho,\vartheta,\zeta)\f$
        @param arg1 parameter
        @returns the contravariant magnetic fields
        """
        coords2d = np.atleast_2d(coords)
        r = coords2d[:,0]
        t = coords2d[:,1]
        z = coords2d[:,2]

        coordsout = self.surfaces.get_coords(r, t, z, derivative=1)
        coordsin = np.array([coordsout.s, coordsout.t, z]).T
        B = self.pb.B(coordsin)

        

    def dBdX(self, coords, args=None):
        """! Returns magnetic fields
        @param coordinates \f$(s,\theta,\zeta)\f$
        @param arg1 parameter
        @returns B, dBdX, the contravariant magnetic fields, the derivatives of them
        """
        B, dB = self.fortran_module.specbfield.get_bfield_tangent(coords)
        return B, dB.T

    def convert_coords(self, stz):
        """! Python wrapper for getting the xyz coordinates from stz
        @param stz  the stz coordinate
        @returns the xyz coordinates
        """
        xyz = self.fortran_module.speccoords.get_xyz(stz)

        # depending on the geometry, return RZ or yx
        if self.Igeometry == 1:
            # for a slab, return x=R, y=theta, z=zeta
            return np.array(
                [
                    xyz[0],
                    np.mod(stz[1], 2 * np.pi) * self.rpol,
                    np.mod(stz[2], 2 * np.pi) * self.rtor,
                ],
                dtype=np.float64,
            )
        if self.Igeometry == 2:
            # for cylinderical geometry, return x=r*cos theta, y=zeta*rtor, z=sin theta
            return np.array(
                [xyz[0] * np.cos(stz[1]), stz[2] * self.rtor, xyz[0] * np.sin(stz[1])],
                dtype=np.float64,
            )
        if self.Igeometry == 3:
            # for toroidal geometry, return x=R, y=zeta, z=Z
            return xyz

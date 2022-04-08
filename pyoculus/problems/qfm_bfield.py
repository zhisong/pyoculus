## @file qfm_bfield.py
#  @brief Setup the magnetic field system - using the Quadratic minimising surfaces coordinates
#  @author Zhisong Qu (zhisong.qu@anu.edu.au)
#
from .toroidal_bfield import ToroidalBfield
from .interpolate_coordinates import SurfacesToroidal
import numpy as np

## Class that used to setup the bfield problem based on QFM
#
class QFMBfield(ToroidalBfield):

    def __init__(self, pb : ToroidalBfield, surfaces : SurfacesToroidal):
        """! Set up the problems based on a known magnetic field and the QFM surfaces
        @param pb  the TorodialBfield problem we will used as the original coordinates and fields
        @param surfaces  the QFM surfaces
        """
        self.pb = pb
        self.surfaces = surfaces
        self.poincare_plot_type = pb.poincare_plot_type
        self.poincare_plot_xlabel = pb.poincare_plot_xlabel
        self.poincare_plot_ylabel = pb.poincare_plot_ylabel
        self.Nfp = pb.Nfp
        self.has_jacobian = pb.has_jacobian
        super().__init__()
    
    def B(self, coords, args=None):
        """! Returns magnetic fields
        @param coords \f$(\rho,\vartheta,\zeta)\f$
        @param arg1 parameter
        @returns the contravariant magnetic fields
        """
 
        r = np.atleast_1d(coords[0])
        t = np.atleast_1d(coords[1])
        z = np.atleast_1d(coords[2])

        coordsout = self.surfaces.get_coords(r, t, z, derivative=1, input1D=True)
        coordsin = np.array([coordsout.s[0], coordsout.t[0], coordsout.z[0]])
        B = np.atleast_2d(self.pb.B(coordsin))
        B = self.surfaces.contra_vector_transform(B, coordsout, has_jacobian=self.pb.has_jacobian, derivative=False).flatten()

        return B

    def dBdX(self, coords, args=None):
        """! Returns magnetic fields
        @param coords \f$(s,\theta,\zeta)\f$
        @param arg1 parameter
        @returns B, dBdX, the contravariant magnetic fields, the derivatives of them
        """

        r = np.atleast_1d(coords[0])
        t = np.atleast_1d(coords[1])
        z = np.atleast_1d(coords[2])

        coordsout = self.surfaces.get_coords(r, t, z, derivative=2, input1D=True)
        coordsin = np.array([coordsout.s[0], coordsout.t[0], coordsout.z[0]])
        B, dBdX = self.pb.dBdX(coordsin)
        B, dBdX = self.surfaces.contra_vector_transform(np.atleast_2d(B), coordsout, has_jacobian=self.has_jacobian, derivative=True, dv=np.atleast_3d(dBdX))

        return B[0], dBdX[0]

    def B_many(self, coords, args=None):
        """! Returns magnetic fields given multiple inputs
        @param coords \f$(\rho,\vartheta,\zeta)\f$
        @param args the precomputed output from self.surfaces.get_coords. If None is given then it will be computed
        @returns the contravariant magnetic fields
        """

        r = coords[:,0]
        t = coords[:,1]
        z = coords[:,2]

        if args is None:
            coordsout = self.surfaces.get_coords(r, t, z, derivative=1, input1D=True)
        else:
            coordsout = args

        B = self.pb.B_many(np.stack([coordsout.s, coordsout.t, z], -1))
        B = self.surfaces.contra_vector_transform(B, coordsout, has_jacobian=self.has_jacobian)

        return B

    def dBdX_many(self, coords, args=None):
        """! Returns magnetic fields
        @param coords \f$(s,\theta,\zeta)\f$
        @param args the precomputed output from self.surfaces.get_coords. If None is given then it will be computed
        @returns B, dBdX, the contravariant magnetic fields, the derivatives of them
        """
        r = coords[:,0]
        t = coords[:,1]
        z = coords[:,2]

        if args is None:
            coordsout = self.surfaces.get_coords(r, t, z, derivative=2, input1D=True)
        else:
            coordsout = args

        B, dBdX = self.pb.dBdX_many(np.stack([coordsout.s, coordsout.t, z], -1))
        B, dBdX = self.surfaces.contra_vector_transform(B, coordsout, has_jacobian=self.pb.has_jacobian, derivative=True, dv=dBdX)

        return B

    def convert_coords(self, stz):
        """! Python wrapper for getting the xyz coordinates from stz
        @param stz  the stz coordinate
        @returns the xyz coordinates
        """

        r = np.atleast_1d(stz[0])
        t = np.atleast_1d(stz[1])
        z = np.atleast_1d(stz[2])

        coordsout = self.surfaces.get_coords(r, t, z, derivative=0, input1D=True)

        return self.pb.convert_coords(np.array([coordsout.s[0], coordsout.t[0], coordsout.z[0]]))


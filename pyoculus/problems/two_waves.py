## @file two_waves.py
#  @brief the perturbed slab problem
#  @author Zhisong Qu (zhisong.qu@anu.edu.au)

from .toroidal_bfield import ToroidalBfield
import numpy as np

##
# A very simple (but still hard) problem for testing the tools and methods
#
# Details in
# S.R. Hudson, Phys. Plasmas 11, 677 (2004).
#
# The magnetic field is given by
# \f[
#    \mathbf{B} = \nabla s \times \nabla \theta - \nabla \chi(s, \theta, \zeta) \times \nabla \zeta
# \f]
#
# The Hamiltonian of the magnetic field is given by
# \f[ \chi(s, \theta,\zeta) = \frac{s^2}{2} - k \left[ \frac{1}{2} \cos (2\theta - \zeta) + \frac{1}{3} \cos(3\theta - 2\zeta) \right] \f]
#
# The magnetic fields are given by
#
#  \f[ B^s = - k \left[ \sin (2\theta - \zeta) + \sin(3\theta - 2\zeta)\right], \quad B^\theta = s , \quad B^\zeta = 1
# \f]
#
# To use the class:
#
#     ps = TwoWaves(k=0.002)
#
class TwoWaves(ToroidalBfield):


    def __init__(self, k=0.002):
        """! Set up the problem
        @param k the value used in the Hamiltonian
        """
        super().__init__()
        self.k = k
        self.Nfp = 1
        self.has_jacobian = True

        ## the problem size, 2 for 1.5D/2D Hamiltonian system
        self.problem_size = 2
        ## by default plotting the yx plane
        self.poincare_plot_type = "yx"
        ## by default x axis has label q
        self.poincare_plot_xlabel = "q"
        ## by default y axis has label p
        self.poincare_plot_ylabel = "p"

    def set_k(self, k):
        """! Set the value of k
        @param k the value used in the Hamiltonian
        """
        self.k = k

    def B(self, coords, *args):
        """! Returns magnetic fields
        @param coordinates
        @param *args extra parameters
        @returns the contravariant magnetic fields
        """
        q = coords[1]
        p = coords[0]
        t = coords[2]

        dqdt = p
        dpdt = -self.k * (np.sin(2 * q - t) + np.sin(3 * q - 2 * t))
        dtdt = 1

        return np.array([dpdt, dqdt, dtdt], dtype=np.float64)

    def dBdX(self, coords, *args):
        """! Returns magnetic fields
        @param coordinates
        @param *args extra parameters
        @returns B, dBdX, the contravariant magnetic fields, the derivatives of them
        """
        q = coords[1]
        p = coords[0]
        t = coords[2]

        dqdt = p
        dpdt = -self.k * (np.sin(2 * q - t) + np.sin(3 * q - 2 * t))
        dtdt = 1

        dpdq = -self.k * (2.0 * np.cos(2 * q - t) + 3.0 * np.cos(3 * q - 2 * t))
        dqdp = 1.0

        dBu = np.zeros([3, 3], dtype=np.float64)

        dBu[0, 1] = dqdp
        dBu[1, 0] = dpdq

        Bu = np.array([dpdt, dqdt, dtdt], dtype=np.float64)

        return Bu, dBu

    def convert_coords(self, incoords):
        """! Convert coordinates for Poincare plots
        @param incoords \f$(p,q,t)\f$
        @returns \f$(p,q \ \mbox{mod} \ 2\pi,t \ \mbox{mod} \ 2\pi)\f$
        """
        return np.array(
            [
                incoords[0],
                np.mod(incoords[1], 2.0 * np.pi),
                np.mod(incoords[2], 2.0 * np.pi),
            ],
            dtype=np.float64,
        )
    
    def B_many(self, x1arr, x2arr, x3arr, input1D=True, *args):
        """! Returns magnetic fields, with multipy coordinate inputs
        @param x1arr the first coordinates. Should have the same length as the other two if input1D=True.
        @param x2arr the second coordinates. Should have the same length as the other two if input1D=True.
        @param x3arr the third coordinates. Should have the same length as the other two if input1D=True.
        @param input1D if False, create a meshgrid with sarr, tarr and zarr, if True, treat them as a list of points
        @param *args extra parameters
        @returns the contravariant magnetic fields
        """
        x1arr = np.atleast_1d(x1arr)
        x2arr = np.atleast_1d(x2arr)
        x3arr = np.atleast_1d(x3arr)

        if not input1D:
            size = (x1arr.size, x2arr.size, x3arr.size)
            p = np.broadcast_to(x1arr[:, np.newaxis, np.newaxis], size).flatten()
            q = np.broadcast_to(x2arr[np.newaxis, :, np.newaxis], size).flatten()
            t = np.broadcast_to(x3arr[np.newaxis, np.newaxis, :], size).flatten()
            n = p.size

        else:
            p = x1arr
            q = x2arr
            t = x3arr
            n = p.size

        dqdt = p
        dpdt = -self.k * (np.sin(2 * q - t) + np.sin(3 * q - 2 * t))
        dtdt = np.ones_like(p)

        Blist = np.stack([dpdt, dqdt, dtdt], -1)

        return Blist

    def dBdX_many(self, x1arr, x2arr, x3arr, input1D=True, *args):
        """! Returns magnetic fields
        @param x1arr the first coordinates. Should have the same length as the other two if input1D=True.
        @param x2arr the second coordinates. Should have the same length as the other two if input1D=True.
        @param x3arr the third coordinates. Should have the same length as the other two if input1D=True.
        @param input1D if False, create a meshgrid with sarr, tarr and zarr, if True, treat them as a list of points
        @param *args extra parameters
        @returns B, dBdX, the contravariant magnetic fields, the derivatives of them
        """
        x1arr = np.atleast_1d(x1arr)
        x2arr = np.atleast_1d(x2arr)
        x3arr = np.atleast_1d(x3arr)

        if not input1D:
            size = (x1arr.size, x2arr.size, x3arr.size)
            p = np.broadcast_to(x1arr[:, np.newaxis, np.newaxis], size).flatten()
            q = np.broadcast_to(x2arr[np.newaxis, :, np.newaxis], size).flatten()
            t = np.broadcast_to(x3arr[np.newaxis, np.newaxis, :], size).flatten()
            n = p.size

        else:
            p = x1arr
            q = x2arr
            t = x3arr
            n = p.size

        dqdt = p
        dpdt = -self.k * (np.sin(2 * q - t) + np.sin(3 * q - 2 * t))
        dtdt = np.ones_like(p)

        Blist = np.stack([dpdt, dqdt, dtdt], -1)

        dpdq = -self.k * (2.0 * np.cos(2 * q - t) + 3.0 * np.cos(3 * q - 2 * t))
        dqdp = np.ones_like(dpdq)

        zeros = np.zeros_like(dpdq)

        dBplist = np.stack([zeros, dpdq, zeros], -1)
        dBqlist = np.stack([dqdp, zeros, zeros], -1)
        dBtlist = np.stack([zeros, zeros, zeros], -1)

        dBlist = np.stack([dBplist, dBqlist, dBtlist], -1)

        return Blist, dBlist
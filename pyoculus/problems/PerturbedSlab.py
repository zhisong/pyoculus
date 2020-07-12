########################################
# PerturbedSlab.py: the perturbed slab problem
# A very simple (but still hard) problem for testing the tools and methods
#
# Details in
# S.R. Hudson, Phys. Plasmas 11, 677 (2004).
#
# The Hamiltonian of the system is given by
# H(q,p,t) = p^2/2 - k[1/2 cos (2q - t) + 1/3 cos(3q - 2t)]
#
# written by @zhisong (zhisong.qu@anu.edu.au)
#

from .BaseProblem import BaseProblem
import numpy as np

class PerturbedSlab(BaseProblem):
    """
    Class that used to setup the perturbed slab problem used in ODE solver.
    A very simple (but still hard) problem for testing the tools and methods
    Details in
    S.R. Hudson, Phys. Plasmas 11, 677 (2004).
    The Hamiltonian of the system is given by
    H(q,p,t) = p^2/2 - k[1/2 cos (2q - t) + 1/3 cos(3q - 2t)]

    Call signature:
        ps = PerturbedSlab(k=0.002) 

    Contains:
        f - function to compute the RHS of the ODE
        f_tangent - function to compute the RHS of the ODE, with tangent
        coords_convert - function that converts curvilinear coordinates to real coordinates
    """

    problem_size = 2

    def __init__(self, k=0.002):
        '''Set up the problem
        parameters:
            k -- the value used in the Hamiltonian
        '''
        self.k = k

        super().__init__()

    def set_k (self, k):
        '''Set the value of k
        parameters:
            k -- the value used in the Hamiltonian
        '''
        self.k = k

    def f(self, t, qp, arg1=None):
        '''The RHS of the Hamilton's equations 
        parameters:
            t -- the zeta coordinate
            qp -- array size 2, the (q,p) coordinate
            arg1 -- parameter for the ODE, not used here

        return:
            array size 2, the RHS of the ODE
        '''

        q = qp[1]
        p = qp[0]

        dqdt = p
        dpdt = - self.k * (np.sin(2*q - t) + np.sin(3*q - 2*t))

        return np.array([dpdt, dqdt], dtype=np.float64)

    def f_tangent(self, t, qp, arg1=None):
        '''The RHS of the Hamilton's equations, with tangent
        parameters:
            t -- the zeta coordinate
            qp -- array size 2, the (q,p) coordinate
            arg1 -- parameter for the ODE, not used here

        return:
            array size 6, the RHS of the ODE, with tangent
        '''
    #return self.fortran_module.bfield.get_bfield_tangent(zeta, st)

    def convert_coords(self, incoords):
        '''We need to mod y and z coordinates by 2 * pi
        '''
        return np.array([incoords[0], np.mod(incoords[1], 2.0*np.pi), np.mod(incoords[2], 2.0*np.pi)], dtype=np.float64)

    
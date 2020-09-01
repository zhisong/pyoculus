## @file PerturbedSlab.py
#  @brief the perturbed slab problem
#  @author Zhisong Qu (zhisong.qu@anu.edu.au)

from .BaseProblem import BaseProblem
import numpy as np

##
# A very simple (but still hard) problem for testing the tools and methods
#
# Details in
# S.R. Hudson, Phys. Plasmas 11, 677 (2004).
#
# The Hamiltonian of the system is given by
# \f[ H(q,p,t) = \frac{p^2}{2} - k \left[ \frac{1}{2} \cos (2q - t) + \frac{1}{3} \cos(3q - 2t) \right] \f]
#
# The ODEs are given by
#
#  \f[ \frac{dq}{dt} = p, \quad \frac{dp}{dt} = - k \left[ \sin (2q - t) + \sin(3q - 2t) \right] \f]
#
# To use the class:
#
#     ps = PerturbedSlab(k=0.002)
#
class PerturbedSlab(BaseProblem):

    ## the problem size, 2 for 1.5D/2D Hamiltonian system
    problem_size = 2
    ## by default plotting the yx plane
    poincare_plot_type = 'yx'
    ## by default x axis has label q
    poincare_plot_xlabel = 'q'
    ## by default y axis has label p
    poincare_plot_ylabel = 'p'

    ## Set up the problem
    # @param k the value used in the Hamiltonian
    def __init__(self, k=0.002):
        self.k = k

        super().__init__()

    ## Set the value of k
    # @param k the value used in the Hamiltonian
    def set_k (self, k):

        self.k = k

    ## The RHS of the Hamilton's equations
    # @param t the zeta coordinate
    # @param qp array size 2, the \f$(p,q)\f$ coordinate
    # @param arg1 parameter for the ODE, not used here, can be set to anything
    # @returns array size 2, the RHS of the ODE
    def f(self, t, qp, arg1=None):

        q = qp[1]
        p = qp[0]

        dqdt = p
        dpdt = - self.k * (np.sin(2*q - t) + np.sin(3*q - 2*t))

        return np.array([dpdt, dqdt], dtype=np.float64)

    ## The RHS of the Hamilton's equations, with tangent
    # @param t the zeta coordinate
    # @param qp array size 6, the \f$(p,q,\Delta p_1, \Delta q_1, \Delta p_2, \Delta q_2 )\f$ coordinate
    # @param arg1 parameter for the ODE, not used here, can be set to anything
    # @returns array size 6, the RHS of the ODE, with tangent
    def f_tangent(self, t, qp, arg1=None):

        q = qp[1]
        p = qp[0]

        dpq = np.array([[qp[2],qp[4]],[qp[3],qp[5]]], dtype=np.float64)
        M = np.zeros([2,2], dtype=np.float64)

        dqdt = p
        dpdt = - self.k * (np.sin(2*q - t) + np.sin(3*q - 2*t))

        M[0,0] = 0.0
        M[0,1] = - self.k * (2.0*np.cos(2*q - t) + 3.0*np.cos(3*q - 2*t))
        M[1,0] = 1.0
        M[1,1] = 0.0

        dqp = np.matmul(M, dpq)

        return np.array([dpdt, dqdt, dqp[0,0], dqp[1,0], dqp[0,1], dqp[1,1]], dtype=np.float64)

    ## Convert coordinates for Poincare plots
    ## @param incoords \f$(p,q,t)\f$
    ## @returns \f$(p,q \ \mbox{mod} \ 2\pi,t \ \mbox{mod} \ 2\pi)\f$
    def convert_coords(self, incoords):
        return np.array([incoords[0], np.mod(incoords[1], 2.0*np.pi), np.mod(incoords[2], 2.0*np.pi)], dtype=np.float64)

    
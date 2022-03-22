## @file qfm.py: class for generating the (weighted) Quadratic Flux Minimising (QFM) surfaces
#  @brief class for generating the QFMs
#  @author Zhisong Qu (zhisong.qu@anu.edu.au)
#

from .base_solver import BaseSolver
from pyoculus.problems import ToroidalBfield
import numpy as np

nax = np.newaxis

class QFM(BaseSolver):
    def __init__(
        self, problem : ToroidalBfield, params=dict(), integrator=None, integrator_params=dict()
    ):
        """! Set up the class of the fixed point finder
        @param problem must inherit pyoculus.problems.BaseProblem, the problem to solve
        @param params dict, the parameters for the solver
        @param integrator the integrator to use, must inherit \pyoculus.integrators.BaseIntegrator, if set to None by default using RKIntegrator (not used here)
        @param integrator_params dict, the parmaters passed to the integrator (not used here)

        <code> params['ntheta']=100 </code> -- the number of theta points for theta integration
        <code> params['nfft_multiplier']=8 </code> -- the extended (multiplier) resolution for FFT
        <code> params['stellar_sym']=True </code> -- if stellarator symmetry is assumed
        """

        if "ntheta" not in params.keys():
            params["ntheta"] = 100

        if "fft_multiplier" not in params.keys():
            params["nfft_multiplier"] = 8

        if "stellar_sym" not in params.keys():
            params["stellar_sym"] = True

        self.sym = params["stellar_sym"]
        self._MM = params["nfft_multiplier"]

        if not isinstance(problem, ToroidalBfield):
            raise ValueError('Currently only support the problem class ToroidalBfield')
            
        super().__init__(problem, params, integrator, integrator_params)

    def action(self, pp, qq):
        # shorthand
        MM = self._MM
        pqNtor = self._pqNtor
        iota = pp / qq

        qN = qq * pqNtor
        Nfft = MM * qN
        dz = 2 * np.pi / ( MM * pqNtor )

        self._nlist = np.arange(0, qN+1)
        self._zeta = np.arange(0, Nfft) * dz
        self._nzq = self._nlist[:, nax] * self._zeta[nax, :] / qq
        self._cnzq = np.cos(self._nzq)
        self._snzq = np.sin(self._nzq)

    def action_integral(self, xx, pp, qq, a):
        """! Computes the action integral, being used in root finding
        @param xx  the packed degrees of freedom. If self.sym == True, it should contain rcn, tsn, nv. Otherwise it should contain rcn, tsn, rsn, tcn, nv.
        @param pp  the poloidal periodicity of the island, should be an integer
        @param pp  the toroidal periodicity of the island, should be an integer
        @param a   the target area
        @returns ff  the equtions to find zeros, see below.

        Construct the Fourier transform of \f$B^\vartheta_i / B^\zeta_i\f$ and \f$B^\rho_i / B^\zeta_i + \bar \nu / (J B^\zeta_i)\f$,
        \f[
        B^\t / B^\z & = & f^c_0 + \sum_{n=1}^{qN} \left[ f^c_n \cos(n\z/q) + f^s_n \sin(n\z/q) \right], \label{eqn:f}
        \f] \f[
        B^\rho / B^\z + \bar \nu / \sqrt g B^\z & = & g^c_0 + \sum_{n=1}^{qN} \left[ g^c_n \cos(n\z/q) + g^s_n \sin(n\z/q) \right], \label{eqn:g} 
        \f]
        
        \item The Fourier harmonics of $\dot\rho$ and $\dot\t$ are given directly by, 
        \begin{eqnarray}
        \dot \t(\z) & = & p/q + \sum_{n=1}^{qN} \left[ - \t^c_n \sin(n\z/q) + \t^s_n\cos(n\z/q) \right](n/q), \label{eqn:tdot} \\
        \dot \rho(\z) & = & \sum_{n=1}^{qN} \left[ - \rho^c_n \sin(n\z/q) + \rho^s_n\cos(n\z/q) \right] (n/q), \label{eqn:rdot}
        \label{eqn:dtrialcurve}
        \end{eqnarray}
        """
        # shorthand
        MM = self._MM
        pqNtor = self._pqNtor
        iota = pp / qq

        qN = qq * pqNtor
        Nfft = MM * qN

        # unpack dof
        nv, rcn, tsn, rsn, tcn = self._unpack_dof(xx, qN)

        r = np.sum(rcn[:,nax] * self._cnzq, axis=0) + np.sum(rsn[:,nax] * self._snzq, axis=0)
        t = np.sum(tcn[:,nax] * self._cnzq, axis=0) + np.sum(tsn[:,nax] * self._snzq, axis=0)

        




        

    def _unpack_dof(self, xx, qN):
        """! Unpack the degrees of freedom into Fourier harmonics
        @param xx  the packed degrees of freedom
        @param qN  Fourier resolution
        @returns nv, rcn, tsn, rsn, tcn
        """
        nv = xx[0]
        rcn = xx[1:qN+2]
        tsn = np.concatenate([0, xx[qN+2:2*qN+2]])
        rsn = np.concatenate([0, xx[2*qN+2:3*qN+3]])
        tcn = xx[3*qN+3:]
        
        return nv, rcn, tsn, rsn, tcn

    def _pack_dof(self, nv, rcn, tsn, rsn=None, tcn=None):
        """! Unpack the degrees of freedom into Fourier harmonics
        @param nv
        @param rcn
        @param tsn
        @param rsn
        @param tcn
        """
        xx = np.concatenate([nv, rcn, tsn[1:], rsn[1:], tcn])
        
        return xx




        

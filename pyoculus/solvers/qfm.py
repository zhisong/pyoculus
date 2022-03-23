## @file qfm.py: class for generating the (weighted) Quadratic Flux Minimising (QFM) surfaces
#  @brief class for generating the QFMs
#  @author Zhisong Qu (zhisong.qu@anu.edu.au)
#

from .base_solver import BaseSolver
from pyoculus.problems import QFMBfield
import numpy as np

nax = np.newaxis

class QFM(BaseSolver):
    def __init__(
        self, problem : QFMBfield, params=dict(), integrator=None, integrator_params=dict()
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

        if "nfft_multiplier" not in params.keys():
            params["nfft_multiplier"] = 8

        if "pqNtor" not in params.keys():
            params["pqNtor"] = 1

        if "stellar_sym" not in params.keys():
            params["stellar_sym"] = True

        self.sym = params["stellar_sym"]
        self._MM = params["nfft_multiplier"]
        self._pqNtor = params["pqNtor"]

        integrator_params["ode"] = problem.f
            
        super().__init__(problem, params, integrator, integrator_params)

    def action(self, pp, qq):
        from scipy.optimize import root
        # shorthand
        MM = self._MM
        pqNtor = self._pqNtor
        iota = pp / qq

        qN = qq * pqNtor
        Nfft = MM * qN
        dz = 2 * np.pi / ( MM * pqNtor )
        self.dz = dz

        self._nlist = np.arange(0, qN+1)
        self._zeta = np.arange(0, Nfft) * dz
        self._nzq = self._nlist[:, nax] * self._zeta[nax, :] / qq
        self._cnzq = np.cos(self._nzq)
        self._snzq = np.sin(self._nzq)

        nv0 = 0
        rcn0=np.zeros(qN+1)
        tcn0=np.zeros(qN+1)
        rsn0=np.zeros(qN+1)
        tsn0=np.zeros(qN+1)

        rcn0[0] = 0.5
        xx0 = self._pack_dof(nv0, rcn0, tsn0, rsn0, tcn0)
        sol = root(self.action_gradient, xx0, args=(pp,qq,0), method='hybr')
        success = sol.success
        if success:
            nv, rcn, tsn, rsn, tcn = self._unpack_dof(sol.x, qN)
        else:
            raise RuntimeError("QFM orbit for pp=" + str(pp) + ",qq=" + str(qq) + ",a=" + str(a) + " not found.")

        print(nv, rcn, tsn, rsn, tcn)

    def action_gradient(self, xx, pp, qq, a):
        """! Computes the action gradient, being used in root finding
        @param xx  the packed degrees of freedom. It should contain rcn, tsn, rsn, tcn, nv.
        @param pp  the poloidal periodicity of the island, should be an integer
        @param qq  the toroidal periodicity of the island, should be an integer
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
        z = self._zeta
        t += iota * z

        area = ( np.sum( t ) + np.pi * pp) * self.dz / (qq*2*np.pi) - pp * np.pi

        B = self._problem.B_many(np.stack([r, t, z], -1))
        
        gBr = B[:,0]
        gBt = B[:,1]
        gBz = B[:,2]

        rhs_tdot = gBt / gBz
        rhs_rdot = gBr / gBz - nv / gBz

        rhs_tdot_fft = np.fft.rfft(rhs_tdot)
        rhs_rdot_fft = np.fft.rfft(rhs_rdot)

        rhs_tdot_fft_cos = np.real(rhs_tdot_fft) / Nfft * 2
        rhs_tdot_fft_cos[0] /= 2

        rhs_tdot_fft_sin = -np.imag(rhs_tdot_fft) / Nfft * 2
        rhs_tdot_fft_sin[0] = 0

        rhs_rdot_fft_cos = np.real(rhs_rdot_fft) / Nfft * 2
        rhs_rdot_fft_cos[0] /= 2

        rhs_rdot_fft_sin = -np.imag(rhs_rdot_fft) / Nfft * 2
        rhs_rdot_fft_sin[0] = 0

        # now pack the function values
        ff = np.zeros_like(xx)

        ff[0] = area - a
        ff[1:qN+2] = rsn * self._nlist /qq - rhs_rdot_fft_cos[0:qN+1]
        ff[qN+2:2*qN+2] = (- rcn * self._nlist / qq - rhs_rdot_fft_sin[0:qN+1])[1:]
        ff[2*qN+2: 3*qN+3] = tsn * self._nlist / qq - rhs_tdot_fft_cos[0:qN+1]
        ff[2*qN+2] += iota
        ff[3*qN+3:] =  (- tcn * self._nlist / qq - rhs_tdot_fft_sin[0:qN+1])[1:]

        return ff

    def _unpack_dof(self, xx, qN):
        """! Unpack the degrees of freedom into Fourier harmonics
        @param xx  the packed degrees of freedom
        @param qN  Fourier resolution
        @returns nv, rcn, tsn, rsn, tcn
        """
        nv = xx[0]
        rcn = xx[1:qN+2]
        tsn = np.concatenate([[0], xx[qN+2:2*qN+2]])
        rsn = np.concatenate([[0], xx[2*qN+2:3*qN+2]])
        tcn = xx[3*qN+2:]
        
        return nv, rcn, tsn, rsn, tcn

    def _pack_dof(self, nv, rcn, tsn, rsn=None, tcn=None):
        """! Unpack the degrees of freedom into Fourier harmonics
        @param nv
        @param rcn
        @param tsn
        @param rsn
        @param tcn
        """
        xx = np.concatenate([[nv], rcn, tsn[1:], rsn[1:], tcn])
        
        return xx




        

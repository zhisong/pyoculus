########################################
# LyapunovExponent: class for computing the LyapunovExponent
# written by @zhisong (zhisong.qu@anu.edu.au)
#
from pyoculus.integrators.BaseIntegrator import BaseIntegrator
from pyoculus.integrators.RKIntegrator import RKIntegrator
from pyoculus.problems.BaseProblem import BaseProblem
import numpy as np

class LyapunovExponent:
    def __init__(self, problem,  params=dict(), integrator=RKIntegrator, integrator_params=dict()):
        '''Set up the class
        parameters:
            system -- BaseProblem class, the problem to solve
            integrator -- the integrator to use, by default using RKIntegrator
            params -- dict, the parameters used in the ODE solver
            integrator_params -- dict, the parmaters passed to the integrator

            params['nPpts']=2000 -- the number of iterations
            params['nsave']=100 -- save the Lyapunov Exponent each nsave iteration
            params['Nfp']=1 -- period in zeta direction
        '''
        if 'nPpts' not in params.keys():
            params['nPpts'] = 2000

        if 'nsave' not in params.keys():
            params['nsave'] = 100

        if 'Nfp' not in params.keys():
            params['Nfp'] = 1
            print('Nfp not specified in params, using 1 as default, please confirm')

        self._params = params
        self.nPpts = params['nPpts']
        self.nsave = params['nsave']
        self.Nfp = params['Nfp']

        integrator_params['ode'] = problem.f_tangent
        self._integrator = integrator(integrator_params)
        self._problem = problem

        # check the integrator and the coordinate fun
        if not isinstance(self._integrator, BaseIntegrator):
            raise Exception('The Integrator is not a derived type of BaseIntegrator class')
        if not isinstance(self._problem, BaseProblem):
            raise Exception('The problem is not a derived type of BaseProblem class')

        self.ile = np.arange([0,self.nPpts+1,self.nsave], dtype=np.float64)
        self.ile = self.ile[2:]
        self.le = np.zeros_like(self.ile)

    def compute(self, t0, ic, dic=[1.0, 0.0]):
        '''Set up the class
        parameters:
            t0 -- the start time (or angle)
            ic -- the initial conidition
            dic -- the initial perturbation direction (non-zero, for most of the cases it doesn't matter)
        '''

        # short hand for the size of the problem
        n = self._problem.problem_size

        # the zeta interval for each integration
        self.dt = 2*np.pi / float(self._params['Nfp'])

        # normalize dic
        dic_norm = np.array(dic, dtype=np.float64).flatten()
        dic_norm = dic_norm / np.sqrt(np.sum(dic_norm**2))

        icnp = np.array(ic, dtype=np.float64).flatten()

        # check the size of the input
        assert len(icnp) == n
        assert len(dic_norm) == n

        self.t0 = t0
        self.ic = icnp.copy()
        self.dic = dic_norm.copy()

        # put together all the initial conditions
        ic_all = np.concatenate((icnp, np.tile(dic_norm, n)))

        # initialize di to save the data
        self.di = np.zeros((self.nPpts), dtype=np.float64)

        t = t0
        dt = self.dt

        for ii in range(self.nPpts):

            # initialize the integrator
            self._integrator.set_initial_value(t, ic_all)

            # run the integrator
            try:
                st = self._integrator.integrate(t + dt)
            except Exception:
                # integration failed, abort for this orbit
                print('Integration failed for s=', ic[0])
                break

            # extract the result
            ic_new = st[0:n]
            dic_new = st[n:2*n]

            # normalize and save di
            self.di[ii] = np.sqrt(np.sum(dic_new**2))
            dic_new = dic_new / self.di[ii]

            # set the new initial condition
            ic_all = np.concatenate((ic_new, np.tile(dic_new, n)))

            # advance in time
            t = t + dt





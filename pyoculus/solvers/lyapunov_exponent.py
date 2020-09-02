## @file lyapunov_exponent.py
# Containing a class for computing the Lyapunov Exponent
#  @author Zhisong Qu (zhisong.qu@anu.edu.au)
#
from pyoculus.problems import BaseProblem
from .base_solver import BaseSolver
import numpy as np

## Class that used to setup the Lyapunov exponent computation
class LyapunovExponent(BaseSolver):
    def __init__(
        self, problem, params=dict(), integrator=None, integrator_params=dict()
    ):
        """! Sets up the class that compute the Lyapunov exponent
        @param problem must inherit pyoculus.problems.BaseProblem, the problem to solve
        @param params dict, the parameters for the solver
        @param integrator the integrator to use, must inherit \pyoculus.integrators.BaseIntegrator, if set to None by default using RKIntegrator
        @param integrator_params dict, the parmaters passed to the integrator

        <code> params['nPpts']=2000 </code> -- the number of iterations

        <code> params['nsave']=100 </code> -- save the Lyapunov Exponent each nsave iteration

        <code> params['Nfp']=1 </code> -- period in zeta direction
        """
        if "nPpts" not in params.keys():
            params["nPpts"] = 2000

        if "nsave" not in params.keys():
            params["nsave"] = 100

        self._params = params
        self.nPpts = params["nPpts"]
        self.nsave = params["nsave"]
        self.Nfp = problem.Nfp

        integrator_params["ode"] = problem.f_tangent

        super().__init__(
            problem=problem,
            params=params,
            integrator=integrator,
            integrator_params=integrator_params,
        )

        self.ile = np.arange(0, self.nPpts + 1, self.nsave, dtype=np.int)
        self.ile = self.ile[2:]
        self.le = np.zeros(self.ile.shape, dtype=np.float64)

    def compute(self, t0, ic, dic=[1.0, 0.0]):
        """! Compute the maximal Lyapunov exponents
        @param t0  the start time (or angle)
        @param ic  the initial conidition
        @param dic the initial perturbation direction (non-zero, for most of the cases it doesn't matter)

        @returns a class with results

        `result.le` -- the computed maximal Lyapunov Exponent (as a function of number of map iterations)

        `result.ile` -- the number of iterations
        """

        # short hand for the size of the problem
        n = self._problem.problem_size

        # the zeta interval for each integration
        self.dt = 2 * np.pi / self.Nfp

        # normalize dic
        dic_norm = np.array(dic, dtype=np.float64).flatten()
        dic_norm = dic_norm / np.sqrt(np.sum(dic_norm ** 2))

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
        self.di = np.ones((self.nPpts), dtype=np.float64)

        t = t0
        dt = self.dt

        for ii in range(self.nPpts):

            # initialize the integrator
            self._integrator.set_initial_value(t, ic_all)

            # run the integrator
            try:
                st = self._integrator.integrate(t + dt)
            except:
                # integration failed, abort for this orbit
                print("Integration failed for s=", ic[0])
                break

            # extract the result
            ic_new = st[0:n]
            dic_new = st[n : 2 * n]

            # normalize and save di
            self.di[ii] = np.sqrt(np.sum(dic_new ** 2))
            dic_new = dic_new / self.di[ii]

            # set the new initial condition
            ic_all = np.concatenate((ic_new, np.tile(dic_new, n)))

            # advance in time
            t = t + dt

        for ii in range(len(self.le)):
            self.le[ii] = np.sum(np.log(self.di[0 : self.ile[ii]])) / self.ile[ii] / dt

        result = LyapunovExponent.OutputData()
        result.ile = self.ile.copy()
        result.le = self.le.copy()

        self.successful = True

        return result

    def plot(self, **kwargs):
        import matplotlib.pyplot as plt

        """! Generates the plot of maximal Lyapunov exponent vs number of iterations
        @param **kwargs passed to the plotting routine "plot"
        """

        if not self.successful:
            raise Exception("A successful call of compute() is needed")

        if plt.get_fignums():
            fig = plt.gcf()
            ax = plt.gca()
        else:
            fig, ax = plt.subplots()

        leplot = ax.plot(
            np.log10(self.ile.astype(np.float64)), np.log10(self.le), **kwargs
        )

        plt.xlabel("Log10 Num of iters", fontsize=20)
        plt.ylabel("Log10 Maximal LE", fontsize=20)
        plt.xticks(fontsize=16)
        plt.yticks(fontsize=16)

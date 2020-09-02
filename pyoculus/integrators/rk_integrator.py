########################################
## @file rk_integrator.py
#  @brief Contains the class of RK ODE integrator
#  @author Zhisong Qu (zhisong.qu@anu.edu.au)
#

from .base_integrator import BaseIntegrator
from scipy.integrate import ode
import numpy as np

## RKIntegrator wraps the explicit Runge-Kutta implimented in scipy.integrate.ode. for use of pyoculus
#
# Default integrator for pyoculus. Not very fast but versatile and robust.
#
# See __init__ for how to set up the integrator
class RKIntegrator(BaseIntegrator):
    def __init__(self, params):
        """! Sets up the ODE solver
        @param params dict, the parameters used in the ODE solver

        <code>params['ode']</code> -- callable f: rhs=f(t,x,arg1), must provide

        <code>params['args']=None</code> -- the argment that will be used to call f

        <code>params['rtol']=1e-7</code> -- relative tolerance

        <code>params['type']='dopri5'</code> -- the type of integrator, 'dopri5' for RK45, 'dop853' for RK853
        """

        # check if the ode is provided. If not, raise an error

        if "ode" not in params.keys():
            raise ValueError("Please specify the ODE to solve for the Integrator class")
        else:
            self.rhs = params["ode"]

        if "type" not in params.keys():
            params["type"] = "dopri5"  # set default to RK45

        if params["type"] not in ["dopri5", "dop853"]:
            raise ValueError(
                "Please specify the correct type of RK solver, dopri5 for RK45, dop853 for RK853"
            )

        if "rtol" not in params.keys():
            params["rtol"] = 1e-7  # set to default value
        self.rtol = params["rtol"]

        if "args" not in params.keys():
            params["args"] = None
        self.args = params["args"]

        # set up the integrator
        self.integrator = ode(self.rhs).set_integrator(
            params["type"], rtol=params["rtol"]
        )

        super().__init__(params)

    def set_initial_value(self, t, x):
        """! Sets up the initial value for the ODE solver
        @param t the start of time
        @param x the start of coordinates
        """

        self.integrator.set_initial_value(x, t).set_f_params(self._params["args"])
        try:
            testoutput = self.rhs(t, x, self.args)
        except:
            print("ODE function not callable")
            raise

        super().set_initial_value(t, x)

    def integrate(self, tend):
        """! Integrates the ODE until `tend`
        @param tend the target end time
        @returns the new value of x
        """
        x_new = self.integrator.integrate(tend)

        if not self.integrator.successful():
            raise Exception("Integration failed")

        self.x = x_new
        self.t = tend
        return x_new

    def copy(self):
        """! Returns a copy of self, to use if want to compute in parallel
        @returns a copy of self
        """

        # set up a new integrator
        return RKIntegrator(self._params)

    @staticmethod
    def _test_fun(t, y, args):
        return [0.1 * np.cos(y[1]), -y[0]]

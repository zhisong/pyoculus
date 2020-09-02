## @file base_integrator.py
#  @brief Contains base class of ODE integrator
#  @author Zhisong Qu (zhisong.qu@anu.edu.au)
#

## Class that used to setup the ODE integrator.
#
# This is an abstract class, should never be used as an instance.
#
# All integrators derived from BaseIntegrator should contain the following member functions
#
#   - set_initial_value -- Set up initial value for the ODE solver
#   - integrate -- Solve the ODE until a given time
#   .
#
# Optional:
#   - copy -- make a copy of the integrator as a new instance
class BaseIntegrator:
    def __init__(self, params):
        """! Set up the ODE solver
        @param params dict, the parameters used in the ODE solver
        """
        self._params = dict(params)

    def set_initial_value(self, t, x):
        """! Set up the initial value for the ODE solver
        @param t the start of time
        @param x the start of coordinates
        """
        self.t = t
        self.x = x

    def integrate(self, tend):
        """! Integrate the ODE until `tend`
        @param tend the target end time
        @returns the new value of x
        """
        raise NotImplementedError("ERROR: Integrator has to implement integrate")

    def get_solution(self):
        """! Get the solution at current time
        @returns the solution
        """
        return self.x

    def copy(self):
        """! Return a copy of self, to use if want to compute in parallel
        @returns a copy of self
        """
        raise NotImplementedError("ERROR: Integrator has to implement copy")

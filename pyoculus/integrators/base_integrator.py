## @file base_integrator.py
#  @brief Contains base class of ODE integrator
#  @author Zhisong Qu (zhisong.qu@anu.edu.au)
#

## Class that used to setup the ODE integrator.
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

    ## Set up the ODE solver
    # @param params dict, the parameters used in the ODE solver
    def __init__(self, params):

        self._params = dict(params)

    ## Set up the initial value for the ODE solver
    # @param t the start of time
    # @param x the start of coordinates
    def set_initial_value(self, t, x):

        self.t = t
        self.x = x

    ## Integrate the ODE until `tend`
    # @param tend the target end time
    # @returns the new value of x
    def integrate(self, tend):

        raise NotImplementedError("ERROR: Integrator has to implement integrate")

    ## Get the solution at current time
    # @returns the solution
    def get_solution(self):

        return self.x

    ## Return a copy of self, to use if want to compute in parallel
    # @returns a copy of self
    def copy(self):

        raise NotImplementedError("ERROR: Integrator has to implement copy")

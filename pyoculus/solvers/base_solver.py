## @file base_solver.py
#  @brief Contains base class for solvers
#  @author Zhisong Qu (zhisong.qu@anu.edu.au)
#

from pyoculus.integrators import BaseIntegrator
from pyoculus.integrators import RKIntegrator
from pyoculus.problems import BaseProblem

## Abstract class that used to setup all other solvers.
class BaseSolver:

    ## Used to return the output data
    class OutputData:
        def __init__(self):
            pass

    def __init__(
        self, problem, params=dict(), integrator=None, integrator_params=dict()
    ):
        """! Sets up the solver
        @param problem must inherit pyoculus.problems.BaseProblem, the problem to solve
        @param params dict, the parameters for the solver
        @param integrator the integrator to use, must inherit \pyoculus.integrators.BaseIntegrator, if set to None by default using RKIntegrator
        @param integrator_params dict, the parmaters passed to the integrator
        """
        ## flagging if the computation is done and successful
        self.successful = False

        # check the integrator
        if integrator is None:
            self._integrator_type = RKIntegrator
        else:
            # check the integrator
            if not issubclass(integrator, BaseIntegrator):
                raise ValueError(
                    "The Integrator is not a derived type of BaseIntegrator class"
                )
            self._integrator_type = integrator

        # check the problem
        if not isinstance(problem, BaseProblem):
            raise ValueError("The problem is not a derived type of BaseProblem class")

        self._params = dict(params)
        self._integrator = self._integrator_type(integrator_params)
        self._problem = problem

        self._integrator_params = dict(integrator_params)

    def is_successful(self):
        """! Returns True if the computation is successfully completed
        @returns successful -- True if the computation is successfully completed, False otherwise
        """
        return self.successful

########################################
# BaseSolver.py: base class for solvers
# written by @zhisong (zhisong.qu@anu.edu.au)
#

from pyoculus.integrators.BaseIntegrator import BaseIntegrator
from pyoculus.integrators.RKIntegrator import RKIntegrator
from pyoculus.problems.BaseProblem import BaseProblem

class BaseSolver:
    """
    Abstract class that used to setup all other solvers.

    Call signature:
        my_solver = BaseSolver(problem, params, integrator, integrator_params) 
    """

    # flagging if the computation is done and successful
    successful = False

    class OutputData:
        """Used to return the output data
        """
        def __init__(self):
            pass

    def __init__(self, problem,  params=dict(), integrator=None, integrator_params=dict()):
        '''Set up the Poincare plot 
        parameters:

            problem -- BaseEquilibrium.BaseProblem class, the problem to solve
            integrator -- the integrator to use, if set to None by default using RKIntegrator
            params -- dict, the parameters used in the ODE solver
            integrator_params -- dict, the parmaters passed to the integrator
        '''
        
        # check the integrator
        if integrator is None:
            self._integrator_type = RKIntegrator
        else:
            # check the integrator
            if not issubclass(integrator, BaseIntegrator):
                raise ValueError('The Integrator is not a derived type of BaseIntegrator class')
            self._integrator_type = integrator

        # check the problem
        if not isinstance(problem, BaseProblem):
            raise ValueError('The problem is not a derived type of BaseProblem class')

        self._params = dict(params)
        self._integrator = self._integrator_type(integrator_params)
        self._problem = problem

        self._integrator_params = dict(integrator_params)

    def is_successful(self):
        '''Return if the computation is successfully completed
        Returns:
            successful -- True if the computation is successfully completed, False otherwise
        '''
        return self.successful
    

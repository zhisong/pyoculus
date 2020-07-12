########################################
# BaseIntegrator.py: base class of ODE integrator
# written by @zhisong (zhisong.qu@anu.edu.au)
#

class BaseIntegrator:
    """
    Class that used to setup the ODE integrator.
    This is an abstract class, should never be used as an instance

    Call signature:
        my_integrator = BaseIntegrator(params) 

    Contains:
        set_initial_value -- Set up initial value for the ODE solver
        integrate -- Solve the ODE until a given time
        get_solution -- Get the solution at current time
        copy -- make a copy of the integrator as a new instance
    """

    def __init__(self, params):
        '''Set up the ODE solver 
        parameters:
            params -- the parameters used in the ODE solver
        '''
        self._params = dict(params)

    def set_initial_value(self, t, x):
        '''Set up the initial value for the ODE solver
        parameters:
            t -- the start of time
            x -- the start of coordinates

        Returns:
            None
        '''
        self.t = t
        self.x = x

    def integrate(self, tend):
        '''Integrate the ODE until tend
        parameters:
            tend -- the target end time

        Returns:
            None
        '''
        raise NotImplementedError('ERROR: Integrator has to implement integrate')

    def get_solution(self):
        '''Get the solution at current time
        Returns:
            x -- the current solution
        '''
        return self.x

    def copy(self):
        '''Return a copy of self
        Returns:
            integrator -- a copy of self
        '''
        raise NotImplementedError('ERROR: Integrator has to implement copy')
########################################
# RKIntegrator.py: class of RK ODE integrator
# uses scipy.integrate.ode
# written by @zhisong (zhisong.qu@anu.edu.au)
#

from .BaseIntegrator import BaseIntegrator
from scipy.integrate import ode
import numpy as np

class RKIntegrator(BaseIntegrator):
    """
    Class that used to setup the RK45 ODE integrator.

    Call signature:
        my_integrator = RKIntegrator(params) 

    Contains:
        set_initial_value -- Set up initial value for the ODE solver
        integrate -- Solve the ODE until a given time
        copy -- make a copy of the integrator as a new instance
    """

    def __init__(self, params):
        '''Set up the ODE solver 
        parameters:
            params -- dict, the parameters used in the ODE solver

            params['ode'] -- callable f: rhs=f(t,x,arg1)
            params['args']=None -- the argment that will be used to call f
            params['rtol']=1e-7 -- relative tolerance
            params['type']='dopri5' -- the type of integrator, 'dopri5' for RK45, 'dop853' for RK853
        '''
        # check if the ode is provided. If not, raise an error

        if 'ode' not in params.keys():
            raise ValueError('Please specify the ODE to solve for the Integrator class')
        else:
            self.rhs = params['ode']

        if 'type' not in params.keys():
            params['type'] = 'dopri5' # set default to RK45
        
        if params['type'] not in ['dopri5', 'dop853']:
            raise ValueError('Please specify the correct type of RK solver, dopri5 for RK45, dop853 for RK853')

        if 'rtol' not in params.keys():
            params['rtol'] = 1e-7   # set to default value
        self.rtol = params['rtol']

        if 'args' not in params.keys():
            params['args'] = None
        self.args = params['args']

        # set up the integrator
        self.integrator = ode(self.rhs).set_integrator(params['type'], rtol=params['rtol'])
        
        super().__init__(params)

    def set_initial_value(self, t, x):
        '''Set up the initial value for the ODE solver
        parameters:
            t -- the start of time
            x -- the start of coordinates

        Returns:
            None
        '''
        self.integrator.set_initial_value(x, t).set_f_params(self._params['args'])
        try:
            testoutput = self.rhs(t,x,self.args)
        except:
            print('ODE function not callable')
            raise

        super().set_initial_value(t,x)

    def integrate(self, tend):
        '''Integrate the ODE until tend
        parameters:
            tend -- the target end time

        Returns:
            x_new -- the integration result
        '''
        x_new = self.integrator.integrate(tend)

        if not self.integrator.successful():
            raise Exception('Integration failed')
        
        self.x = x_new
        self.t = tend
        return x_new

    def copy(self):
        '''Return a copy of self
        Returns:
            integrator -- a copy of self
        '''
        # set up a new integrator
        return RKIntegrator(self._params)

    @staticmethod
    def _test_fun(t, y, args):
        return [0.1*np.cos(y[1]), -y[0]]
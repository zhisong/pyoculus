from pyoculus.problems.BaseProblem import BaseProblem
from .BaseSolver import BaseSolver
from .FixedPoint import FixedPoint
import numpy as np

class Manifolds(BaseSolver):
    """
    Computes the stable and unstable manifolds at the fixed point. 
    We create the tangent map at each fixed point and extract the eigenvalues and eigenvector. If eigenvalues are found
    negative , we then integrate the magnetic fields line flow equation to the unstable eigenvector coreesponds to
    negative eigenvalues to form unstable manifold and vice-versa.
    """

    def __int__(self, problem,  params=dict(), integrator=None, integrator_params=dict()):

        '''
        system -- BaseProblem class, the problem to solve
            integrator -- the integrator to use, by default using RKIntegrator
            params -- dict, the parameters used in the ODE solver
            integrator_params -- dict, the parmaters passed to the integrator

            params['nPpts']=2000 -- the number of iterations
            params['nsave']=100 -- save the Lyapunov Exponent each nsave iteration
            params['Nfp']=1 -- period in zeta direction
        '''

        if 'theta' not in params.keys():
            params['theta'] = 0.0
            
        if 'zeta' not in params.keys():
            params['zeta'] = 0.0

        if 'nPpts' not in params.keys():
            params['nPpts'] = 500

        if 'nPtrj' not in params.keys():
            params['nPtrj'] = 10

        if 'sbegin' not in params.keys():
            params['sbegin'] = -1.0

        if 'send' not in params.keys():
            params['send'] = 1.0

        if 'Nfp' not in params.keys():
            params['Nfp'] = 1
            print('Nfp not specified in params, using 1 as default, please confirm')

        if 'nthreads' not in params.keys():
            params['nthreads'] = 1

        integrator_params['ode'] = problem.f
        super().__init__(problem=problem,  params=params, integrator=integrator, integrator_params=integrator_params)

        # set up the result array
        self.x = np.zeros([self._params['nPtrj']+1,self._params['nPpts']+1],dtype=np.float64)
        self.y = np.zeros_like(self.x)
        self.z = np.zeros_like(self.x)
        self.s = np.zeros_like(self.x)
        self.theta = np.zeros_like(self.x)
        self.zeta = np.zeros_like(self.x)

        def compute(self,veclam_1,veclam_2):
        '''Integration of untable and stable manifolds 
        Returns:
            pdata -- a class that contains the results
            pdata.x,pdata.y,pdata,z -- the Poincare data in xyz coordinates
            pdata.s,pdata,theta,pdata,zeta -- the Poincare data in s,theta,zeta coordinates
        '''

        self.successful = False


        fxdata=self.x
        fydata=self.y

        dUS= 1E-02

        mdata=np.array([fxdata + dUS + veclam_1,fydata + dUS + veclam_2])

        # the zeta interval for each integration
        self.dt = 2*np.pi / float(self._params['Nfp'])

        for ii in range(self._params['nPtrj']+1):
            # find which s to start with
            swindow = self._params['send']-self._params['sbegin']
            ds = swindow / float(self._params['nPtrj'])
            self.s[ii,0] = self._params['sbegin'] + ds * float(ii)

            # put in the initial conditions
            self.theta[ii,0] = self._params['theta']
            self.zeta[ii,0] = self._params['zeta']

        if self._params['nthreads'] == 1: # single thread, do it straight away

            for ii in range(self._params['nPtrj']+1):
                # set up the initial value
                ic = [self.s[ii,0], self.theta[ii,0]]
                t0 = self.zeta[ii,0]

                # initialize the integrator
                self._integrator.set_initial_value(t0, ic)

                t = t0
                dt = self.dt

                for jj in range(1,self._params['nPpts']+1):

                    # run the integrator
                    try:
                        st = self._integrator.integrate(t + dt)
                    except Exception:
                        # integration failed, abort for this orbit
                        print('Integration failed for s=', ic[0])
                        break

                    # extract the result to s theta zeta
                    self.s[ii,jj] = st[0]
                    self.theta[ii,jj] = st[1]
                    self.zeta[ii,jj] = t + dt

                    # advance in time
                    t = t + dt
        
        self.successful = True

        pdata = PoincarePlot.OutputData()
        pdata.x=self.x.copy()
        pdata.y=self.y.copy()
        pdata.z=self.z.copy()
        pdata.s=self.s.copy()
        pdata.theta=self.theta.copy()
        pdata.zeta=self.zeta.copy()

        return pdata




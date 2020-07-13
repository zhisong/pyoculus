########################################
# FixPoint.py: class for finding fixed points
# written by @zhisong (zhisong.qu@anu.edu.au)
#

from .BaseSolver import BaseSolver
import numpy as np

class FixedPoint(BaseSolver):
    """
    Class that used to setup the fixed point finder.

    Call signature:
        my_fixedpoint = FixPoint(problem, params, integrator, integrator_params) 

    Contains:
        compute -- find the fixed point
        plot -- plot the fixed point
    """

    successful = False
    is_theta_fixed = True

    def __init__(self, problem,  params=dict(), integrator=None, integrator_params=dict()):
        '''Set up the class of the fixed point finder.
        parameters:
            system -- BaseProblem class, the problem to solve
            integrator -- the integrator to use, by default using RKIntegrator
            params -- dict, the parameters used in the ODE solver
            integrator_params -- dict, the parmaters passed to the integrator

            params['niter']=100 -- the maximum number of Newton iterations
            params['theta']=None -- if we look for fixed point on some symmetry line
                                    =None : theta is also a free variable to look for
                                    =somenumber : only look for theta with this number
            params['zeta']=0.0 -- the toroidal plane we are after
            params['nrestart']=1 -- if search failed, the number of time to restart (randomly within the domain)
            params['Nfp']=1 -- the periodicity
        '''

        if 'niter' not in params.keys():
            params['niter'] = 100
        
        if 'theta' not in params.keys():
            params['theta'] = None

        if 'zeta' not in params.keys():
            params['zeta'] = 0.0

        if 'nrestart' not in params.keys():
            params['nrestart'] = 1

        if 'Nfp' not in params.keys():
            params['Nfp'] = 1
            print('Nfp not specified in params, using 1 as default, please confirm')

        integrator_params['ode'] = problem.f_tangent

        super().__init__(problem=problem,  params=params, integrator=integrator, integrator_params=integrator_params)
        
        self.Nfp = params['Nfp']
        self.niter = params['niter']
        self.nrestart = params['nrestart'] 
        if params['theta'] is None:
            self.is_theta_fixed = False
        else:
            self.is_fixed_theta = True

    def compute(self, guess, pp, qq, sbegin=-1.0, send=1.0, tol=None):
        '''Look for the fixed point with rotation number pp/qq
            guess -- the initial guess 
                    [s, theta], if params['theta'] == None
                    [s], if params['theta'] == somevalue
            pp,qq -- integers, the numerator and denominator of the rotation number
            sbegin=-1.0 -- the allowed minimum s
            send=1.0 -- the allowed maximum s
            tol=self._integrator_params['rtol']*qq -- the tolerance of the fixed point
        '''
        
        if not isinstance(pp,int) or not isinstance(qq,int):
            raise Exception('pp and qq should be integers')

        if tol is None:
            tol = self._integrator_params['rtol']*qq

        if pp * qq >= 0:
            pp = int(np.abs(pp))
            qq = int(np.abs(qq))
        else:
            pp = -int(np.abs(pp))
            qq = int(np.abs(qq))

        self.pp = pp
        self.qq = qq
        self.dzeta = 2 * np.pi / self.Nfp

        # arrays that save the data
        self.s = np.zeros([qq+1], dtype=np.float64)
        self.theta = np.zeros([qq+1], dtype=np.float64)
        self.zeta = np.zeros([qq+1], dtype=np.float64)
        self.x = np.zeros([qq+1], dtype=np.float64)
        self.y = np.zeros([qq+1], dtype=np.float64)
        self.z = np.zeros([qq+1], dtype=np.float64)

        self.history = []

        guess = np.array(guess, dtype=np.float64)

        # set up the guess
        s_guess = guess[0]
        if self._params['theta'] is None:
            theta_guess = guess[1]
        else:
            theta_guess = self._params['theta']

        # run the Newton's method
        for ii in range(self._params['nrestart']+1):
            try: # run the solver, if failed, try a different random initial condition
                if self.is_theta_fixed:
                    result = self._newton_method_1(pp,qq,s_guess,sbegin,send,theta_guess,self._params['zeta'],self.dzeta,self.niter,tol)
                else:
                    result = self._newton_method_2(pp,qq,s_guess,sbegin,send,theta_guess,self._params['zeta'],self.dzeta,self.niter,tol)
            except:
                result = None
            
            if result is not None: # if it is successful:
                break
            else: # not successful, change to a random initial condition
                print('Search failed: starting from a random initial guesss!')
                random_guess = np.random.rand(2)
                s_guess = sbegin + (send - sbegin) * random_guess[0]
                if self._params['theta'] is None:
                    theta_guess = random_guess[1] * 2 * np.pi

        # now we go and get all the fixed points by iterating the map
        if result is not None:
            t = self.zeta[0]
            dt = 2*np.pi / self.Nfp

            self.s[0] = result[0]
            self.theta[0] = result[1]
            self.zeta[0] = self._params['zeta']

            ic = np.array([result[0], result[1], 1.0, 0.0, 0.0, 1.0], dtype=np.float64)
            self._integrator.set_initial_value(t,ic)

            # integrate to get a series of fixed points
            for jj in range(1,qq+1):

                # run the integrator
                st = self._integrator.integrate(t + dt)
                
                # extract the result to s theta zeta
                self.s[jj] = st[0]
                self.theta[jj] = st[1]
                self.zeta[jj] = t+dt

                # advance in time
                t = t + dt

            # convert coordinates
            for jj in range(0,qq+1):
                stz = np.array([self.s[jj],self.theta[jj],self.zeta[jj]],dtype=np.float64)
                xyz = self._problem.convert_coords(stz)
                self.x[jj] = xyz[0]
                self.y[jj] = xyz[1]
                self.z[jj] = xyz[2]

            rdata = FixedPoint.OutputData()
            rdata.x = self.x.copy()
            rdata.y = self.y.copy()
            rdata.z = self.z.copy()
            rdata.s = self.s.copy()
            rdata.theta = self.theta.copy()
            rdata.zeta = self.zeta.copy()

            # the jacobian
            rdata.jacobian = np.array([[st[2],st[4]],[st[3],st[5]]], dtype=np.float64)
            # the trace
            rdata.trace = np.trace(rdata.jacobian)

            # set the successful flag
            self.successful = True

        else:
            rdata = None
            print('Fixed point search unsuccessful for pp/qq=',pp,'/',qq)

        return rdata

    def plot(self, plottype=None, xlabel=None, ylabel=None, xlim=None, ylim=None, **kwargs):
        import matplotlib.pyplot as plt
        '''generate the plot for fixed points
        parameters:
            plottype -- which variables to plot: 'RZ' or 'yx', by default using "poincare_plot_type" in problem
            xlabel, ylabel -- what to put for the xlabel and ylabel, by default using "poincare_plot_xlabel" in problem
            xlim, ylim -- the range of plotting, by default plotting the range of all data
            **kwargs -- passed to the plotting routine "plot"
        '''

        if not self.successful:
            raise Exception('A successful call of compute() is needed')

        # default setting
        if plottype is None:
            plottype = self._problem.poincare_plot_type
        if xlabel is None:
            xlabel = self._problem.poincare_plot_xlabel
        if ylabel is None:
            ylabel = self._problem.poincare_plot_ylabel

        if plottype=='RZ':
            xdata = self.x
            ydata = self.z
        elif plottype=='yx':
            xdata = self.y
            ydata = self.x
        else:
            raise Exception('Choose the correct type for plottype')

        if plt.get_fignums():
            fig = plt.gcf()
            ax = plt.gca()
            newfig = False
        else:
            fig, ax = plt.subplots()
            newfig = True

        # set default plotting parameters
        # use x
        if kwargs.get('marker') is None:
            kwargs.update({'marker': 'x'})
        # use gray color
        if kwargs.get('c') is None:
             kwargs.update({'c': 'black'})
            
        xs = ax.plot(xdata, ydata, linestyle="None", **kwargs)

        if not newfig:
            if plottype=='RZ':
                plt.axis('equal')
            if plottype=='yx':
                pass

            plt.xlabel(xlabel, fontsize=20)
            plt.ylabel(ylabel, fontsize=20)
            plt.xticks(fontsize=16)
            plt.yticks(fontsize=16)

            if xlim is not None:
                plt.xlim(xlim)
            if ylim is not None:
                plt.ylim(ylim)

    def _newton_method_1(self, pp, qq, s_guess, sbegin, send, theta, zeta, dzeta, niter, tol):
        '''driver to run Newton's method for one variable s
            pp,qq -- integers, the numerator and denominator of the rotation number
            s_guess -- the guess of s
            sbegin -- the allowed minimum s
            send -- the allowed maximum s
            theta -- the theta value (fixed)
            zeta -- the toroidal plain to investigate
            dzeta -- period in zeta
            niter -- the maximum number of iterations
            tol -- the tolerance of finding a fixed point
        '''

        s = s_guess

        # set up the initial condition
        ic = np.array([s, theta, 1.0, 0.0, 0.0, 1.0], dtype=np.float64)
        self.history.append(ic[0:1].copy())

        t0 = zeta
        dt = dzeta

        succeeded = False

        for ii in range(niter):
            t = t0
            self._integrator.set_initial_value(t0, ic)

            for jj in range(qq):
                output = self._integrator.integrate(t + dt)
                t = t+dt

            dtheta = output[1] - theta - 2 * np.pi * pp
            jacobian = output[3]

            # if the resolution is good enough
            if (abs(dtheta) < tol):
                succeeded = True
                break
            s_new = s - dtheta / jacobian
            s = s_new
           
            if s > send or s < sbegin: # search failed, return None
                return None

            ic = np.array([s, theta, 1.0, 0.0, 0.0, 1.0], dtype=np.float64)
            self.history.append(ic[0:1].copy())

        if succeeded:
            return np.array([s,theta,zeta], dtype=np.float64)
        else:
            return None

    def _newton_method_2(self, pp, qq, s_guess, sbegin, send, theta_guess, zeta, dzeta, niter, tol):
        '''driver to run Newton's method for two variable (s,theta)
            pp,qq -- integers, the numerator and denominator of the rotation number
            s_guess -- the guess of s
            sbegin -- the allowed minimum s
            send -- the allowed maximum s
            theta_guess -- the guess of theta
            zeta -- the toroidal plain to investigate
            dzeta -- period in zeta
            niter -- the maximum number of iterations
            tol -- the tolerance of finding a fixed point
        '''

        self.successful = False

        s = s_guess
        theta = theta_guess

        # set up the initial condition
        ic = np.array([s, theta, 1.0, 0.0, 0.0, 1.0], dtype=np.float64)
        self.history.append(ic[0:1].copy())

        t0 = zeta
        dt = dzeta

        succeeded = False

        st = np.array([s,theta],dtype=np.float64)

        for ii in range(niter):
            t = t0
            self._integrator.set_initial_value(t0, ic)
            for jj in range(qq):
                output = self._integrator.integrate(t + dt)
                t = t + dt

            dtheta = output[1] - theta - 2 * np.pi * pp
            ds = output[0] - s
            dst = np.array([ds,dtheta],dtype=np.float64)
            jacobian = np.array([[output[2],output[4]],[output[3],output[5]]], dtype=np.float64)

            # if the resolution is good enough
            if (np.sqrt(dtheta**2 + ds**2) < tol):
                succeeded = True
                break

            # Newton's step
            st_new = st - np.matmul(np.linalg.inv(jacobian - np.eye(2)), dst)
            s = st_new[0]
            theta = st_new[1]
            st = st_new

            if s > send or s < sbegin: # search failed, return None
                return None

            ic = np.array([s, theta, 1.0, 0.0, 0.0, 1.0], dtype=np.float64)
            self.history.append(ic[0:1].copy())

        if succeeded:
            self.successful = True
            return np.array([s,theta,zeta], dtype=np.float64)
        else:
            return None
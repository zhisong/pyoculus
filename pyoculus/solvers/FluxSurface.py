########################################
# FluxSurface.py: class for finding flux surfaces
# written by @zhisong (zhisong.qu@anu.edu.au)
#

from .BaseSolver import BaseSolver
from .FixedPoint import FixedPoint
import pyoculus.irrationals as ir
import numpy as np

class FluxSurface(BaseSolver):
    """
    Class that used to set up the flux surface finder.

    Call signature:
        my_fluxsurface = FluxSurface(problem, params, integrator, integrator_params) 

    Contains:
        compute -- find the fluxsurface
        plot -- plot the fluxsurface
    """
    def __init__(self, problem,  params=dict(), integrator=None, integrator_params=dict()): 
        '''Set up the class of the flux surface point finder.
        parameters:
            system -- BaseProblem class, the problem to solve
            integrator -- the integrator to use, by default using RKIntegrator
            params -- dict, the parameters used in the ODE solver
            integrator_params -- dict, the parmaters passed to the integrator

            These parameters will be passed to the fixed point finder
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
            raise ValueError('We only support located fixed points for a fixed theta at the moment')

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
        
    def compute(self, iota, n_expand=10, nstart=5, sbegin=-1.0, send=1.0, sguess=0.0, fixed_point_left=None, fixed_point_right=None, tol=None):
        '''Look for the flux surface with a irrational rotation number
        parametes:
            iota -- the irrational! rotation number of the flux surface
            fixed_point_left -- a sucessfully found FixPoint to mark the left bound of the flux surface, 
                                its rotation number needs to be in the convergent sequence of iota
            fixed_point_right -- a sucessfully found FixPoint to mark the right bound of the flux surface,
                                its rotation number needs to be in the convergent sequence of iota and next to fixed_point_left
            n_expand=10 -- the number of terms in the continued fraction expansion of iota, used to approximate the flux surface

        Returns:
            fdata -- a class that contains the results
            fdata.MackayResidue -- the Mackay Residue of the fixed points
            fdata.fixed_points -- all the fixed point located
            fdata.rmnc, rmns, zmnc, zmns -- the Fourier harmonics
        '''

        # check if the user specified fixed points are legal

        # iota will be divided by Nfp
        iota = iota / self.Nfp

        # continued fraction expansion of the input irrational
        ai = ir.expandcf(iota,n_expand+1)

        fpleft = None
        fpright = None

        # determine how the input fixed points fit into the sequence of expansion
        for ii in range(n_expand-1):
            ppqq1 = ir.fromcf(ai[0:ii+1])
            pp1 = ppqq1[0]
            qq1 = ppqq1[1]

            ppqq2 = ir.fromcf(ai[0:ii+2])
            pp2 = ppqq2[0]
            qq2 = ppqq2[1]

            # put the lower order fixed point as fpleft and higher order fpright
            if pp1==fixed_point_left.pp and qq1==fixed_point_left.qq and pp2==fixed_point_right.pp and qq2==fixed_point_right.qq:
                fpleft = fixed_point_left
                fpright = fixed_point_right
                nstart = ii
            elif pp2==fixed_point_left.pp and qq2==fixed_point_left.qq and pp1==fixed_point_right.pp and qq1==fixed_point_right.qq:
                fpleft = fixed_point_right
                fpright = fixed_point_left
                nstart = ii

        if fpleft is None or fpright is None:
            raise ValueError('The input fixed points are illegal')

        fixedpoints = [fpleft, fpright]
        self.nstart = nstart

        for ii in range(nstart+2,n_expand):

            ppqq = ir.fromcf(ai[:ii+1])
            pp = ppqq[0]
            qq = ppqq[1]
            iotatarget = float(pp)/float(qq)

            nextfixedpoint = FixedPoint(self._problem, params=self._params, integrator=self._integrator_type, integrator_params=self._integrator_params)

            sleft = fixedpoints[-2].s[0]
            iotaleft = float(fixedpoints[-2].pp)/float(fixedpoints[-2].qq)
            sright = fixedpoints[-1].s[0]
            iotaright = float(fixedpoints[-1].pp)/float(fixedpoints[-1].qq)

            # interpolate between sleft and sright to get the next guess of s
            sguess = sleft + (sright - sleft) / (iotaright - iotaleft) * (iotatarget - iotaleft)

            # the lower and upper bound of s range
            ssmall = np.min([sleft, sright])
            sbig = np.max([sleft, sright])

            fp = nextfixedpoint.compute(sguess,pp,qq,sbegin=ssmall,send=sbig,tol=tol)

            if not nextfixedpoint.successful:
                raise Exception('Fixed point not found')

            fixedpoints.append(nextfixedpoint)
        
        # save the fixed points found
        self.fixedpoints = fixedpoints

        # assemble the output data
        fdata = FluxSurface.OutputData()
        fdata.fixedpoints = fixedpoints

        # put the flag as successful
        self.successful = True

        return fdata

    def plot(self, plottype=None, xlabel=None, ylabel=None, xlim=None, ylim=None, **kwargs):
        import matplotlib.pyplot as plt
        '''generate the plot for flux surface
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
            xdata = self.fixedpoints[-1].x
            ydata = self.fixedpoints[-1].z
        elif plottype=='yx':
            xdata = self.fixedpoints[-1].y
            ydata = self.fixedpoints[-1].x
        else:
            raise ValueError('Choose the correct type for plottype')

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

    def plot_residue(self):
        import matplotlib.pyplot as plt
        '''generate the plot for residue
        '''

        gamma = ((np.sqrt(5) + 1) / 2)

        xlist_greene = np.arange(1, len(self.fixedpoints)+1)
        greenes_list = np.zeros(len(self.fixedpoints), dtype=np.float64)

        for ii, fp in enumerate(self.fixedpoints):
            greenes_list[ii] = fp.GreenesResidue

        xlist_Mackay = np.arange(2, len(self.fixedpoints)+1)
        Mackay_list = np.zeros(len(self.fixedpoints)-1, dtype=np.float64)

        for ii in range(len(self.fixedpoints)-1):
            Mackay_list[ii] = (self.fixedpoints[ii].GreenesResidue + gamma*self.fixedpoints[ii+1].GreenesResidue) / (1.0+gamma)

        fig, ax = plt.subplots()

        geplot = ax.plot(xlist_greene, greenes_list, label='Greene')
        mcplot = ax.plot(xlist_Mackay, Mackay_list, label='Mackay')

        ax.legend()

        plt.xlabel('Order of fixed point', fontsize=20)
        plt.ylabel('Residue', fontsize=20)
        plt.xticks(fontsize=16)
        plt.yticks(fontsize=16)

## @file qfm.py: class for generating the (weighted) Quadratic Flux Minimising (QFM) surfaces
#  @brief class for generating the QFMs
#  @author Zhisong Qu (zhisong.qu@anu.edu.au)
#

from .base_solver import BaseSolver
import numpy as np

class QFM(BaseSolver):
    def __init__(
        self, problem, params=dict()
    ):
        """! Set up the class of the fixed point finder
        @param problem must inherit pyoculus.problems.BaseProblem, the problem to solve
        @param params dict, the parameters for the solver

        <code> params['ntheta']=100 </code> -- the number of theta points for theta integration
        <code> params['nzeta']=100 </code> -- the number of zeta points for integration
        """

        if "ntheta" not in params.keys():
            params["ntheta"] = 100

        if "nzeta" not in params.keys():
            params["nzeta"] = 100
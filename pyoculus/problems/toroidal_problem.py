## @file toroidal_problem.py
#  @brief containing a problem class with two cyclical coordinates for pyoculus ODE solver
#  @author Zhisong Qu (zhisong.qu@anu.edu.au)
#

from .base_problem import BaseProblem

class ToroidalProblem(BaseProblem):

    def __init__(self):
        """! Set up the problem with two cyclical coordinates, e.g. \f[ (s, \theta, \zeta)  \f]
        """
        self.problem_size = 2
        self.poincare_plot_type = "yx"
        self.poincare_plot_xlabel = "x"
        self.poincare_plot_ylabel = "y"

    def f(self, zeta, st, *args):
        """! Returns ODE RHS
        @param zeta cylindrical angle in ODE
        @param st \f$(s, \theta)\f$ in ODE
        @param *args extra parameters for the ODE
        @returns the RHS of the ODE
        """
        raise NotImplementedError("A problem class should implement member function f")

    def f_tangent(self, t, st, *args):
        """! Returns ODE RHS, with tangent
        @param zeta cylindrical angle in ODE
        @param st \f$(s, \theta, ds_1, d\theta_1, ds_2, d\theta_2)\f$ in ODE
        @param *args extra parameters for the ODE
        @returns the RHS of the ODE, with tangent
        """
        raise NotImplementedError(
            "A problem class should implement member function f_tangent"
        )


    def convert_coords(self, coord1):
        """! Converts coordinates, identity here
        @param coords1 the coordinates to convert
        @returns the converted coordinates
        """
        # by default there is no conversion
        return coord1
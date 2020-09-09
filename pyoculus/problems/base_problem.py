## @file base_problem.py
#  @brief containing base problem class for pyoculus ODE solver
#  @author Zhisong Qu (zhisong.qu@anu.edu.au)
#

## Abstract class that used to drive calculation.
#
# This is an abstract class, should never be used as an instance.
#
# All problems should be a derived type of BaseProblem.
# They should contain the ODEs to solve, as specified by the Hamiltonian system of that problem.
#
# ## Hamilton's equation
# ### 1.5D
# If the Hamiltonian system is 1.5 dimensional, the Hamiltonian takes the form \f$H=H(p,q,t)\f$.
# The Hamilton's equations are given by
#
# \f[ \frac{dq}{dt} = \frac{\partial H}{\partial p}, \quad \frac{dp}{dt} = -\frac{\partial H}{\partial q}. \quad\f]
#
# ### 2D
# If the Hamiltonian system is two-dimensional, the Hamiltonian takes the form \f$H=H(p_1,q_1,p_2,q_2)\f$.
# The Hamilton's equations are given by
#
# \f[ \frac{dq_i}{dt} = \frac{\partial H}{\partial p_i}, \quad \frac{dp_i}{dt} = -\frac{\partial H}{\partial q_i}, \quad\f]
#
# in which \f$i=1,2\f$.
#
# The Hamiltonian \f$H\f$ itself is a constant of motion, i.e. \f$dH/dt=0\f$.
#
# Given \f$ p_1, q_1, q_2\f$, we can invert \f$H=H(p_1,q_1,p_2,q_2)=\mbox{const} \f$ to get \f$p_2 = p_2(p_1, q_1, q_2)\f$.
#
# If further \f$q_2 \f$ is monotonically increasing with \f$t\f$, one can replace \f$t\f$ by \f$q_2\f$ and get
#
# \f[ \frac{dq_1}{dq_2} = \left.\frac{\frac{\partial H}{\partial p_1}}{\frac{\partial H}{\partial p_2}} \right|_{p_2(p_1, q_1, q_2)}, \f]
# \f[ \frac{dp_1}{dq_2} = \left.\frac{-\frac{\partial H}{\partial q_1} }{\frac{\partial H}{\partial p_2}}\right|_{p_2(p_1, q_1, q_2)}, \f]
# after computing which the value of \f$p_2\f$ is substituted into these equations.
#
# ## General ODEs with two variables
# The Hamilton's equations, after converted can be combined and written as
#
# \f[ \frac{d y_1}{dt} = F_1(t, y_1, y_2),\f]
# \f[ \frac{d y_2}{dt} = F_2(t, y_1, y_2).\f]
#
# The tangent of the above ODEs are given by
# \f[ \frac{d}{dt} \begin{bmatrix}\Delta y_1\\ \Delta y_2 \end{bmatrix} = 
#     \begin{bmatrix}   
#         \partial_{y_1}F_1 & \partial_{y_2}F_1  \\
#         \partial_{y_1}F_2 & \partial_{y_2}F_2
#     \end{bmatrix}
#     \cdot
#     \begin{bmatrix}\Delta y_1\\ \Delta y_2 \end{bmatrix},
# \f]
# where the partial derivatives are evaluated at \f$y_1(t), y_2(t)\f$ as the solution of the original ODEs. 
#
# ## Implementing a Problem class
# All problem classes should inherit the BaseProblem class.
#
#     class SomeProblem(BaseProblem):
#         def __init__(self, params):
#             """some codes to initialize the problem"""
#         def f(self, t, y, args=None):
#             """some codes compute the RHS (F) of the ODEs given t and y"""
#         def f_tangent(self, t, y, args=None):
#             """some codes compute the ODEs and the tangent the ODEs given t and y, y contains two tangent vectors"""
#         def convert_coords(self, coord1):
#             """some codes compute the coordinate transformation"""
#
# The member function `f` should be implemented in any case. The `tangent` and `convert_coords` functions are optional for some solvers.
#
class BaseProblem:

    # this is the size of the ODE system

    def __init__(self):
        self.problem_size = 2
        self.poincare_plot_type = "yx"
        self.poincare_plot_xlabel = "y"
        self.poincare_plot_ylabel = "x"

    def f(self, t, y, args=None):
        """! Returns ODE RHS
        @param t time in ODE
        @param y variables in ODE
        @param arg1 parameter for the ODE
        @returns the RHS of the ODE
        """
        raise NotImplementedError("A problem class should implement member function f")

    def f_tangent(self, t, y, args=None):
        """! Returns ODE RHS, with tangent
        @param t time in ODE
        @param y \f$\bf y\f$ variables in ODE and \f$\Delta \mathbf{y}_1\f$, \f$\Delta \mathbf{y}_2\f$
        @param arg1 parameter for the ODE
        @returns the RHS of the ODE, with tangent
        """
        raise NotImplementedError(
            "A problem class should implement member function f_tangent"
        )

    def convert_coords(self, coord1):
        """! Converts coordinates (for example \f$s,\theta,\zeta\f$ to \f$R,Z,\varphi\f$)
        @param coords1 the coordinates to convert
        @returns the converted coordinates
        """
        # by default there is no conversion
        return coord1

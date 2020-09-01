########################################
# BaseProblem.py: base system class for oculus ODE solver
# written by @zhisong (zhisong.qu@anu.edu.au)
#

class BaseProblem:

    """
    Abstract class that used to drive calculation, used in ODE solver.
    Call signature:
        Myproblem = BaseProblem() 

    Contains:
        f - function to compute the RHS of the ODE
        f_tangent - function to compute the RHS of the ODE, with tangent
        coords_convert - function that converts curvilinear coordinates to real coordinates
    """

    # this is the size of the ODE system
    problem_size = 2

    # specifies the way to present the phase space
    poincare_plot_type = 'yx'
    poincare_plot_xlabel = 'y'
    poincare_plot_ylabel = 'x'

    def __init__(self):
        pass

    def f(self, t, y, args=None):
        '''Returns ODE RHS 
        parameters:
            t -- time in ODE
            y -- y(t) in ODE
            arg1 -- parameter for the ODE

        return:
            the RHS of the ODE
        '''
        raise NotImplementedError('A system class should implement member function f')
        return y
    
    
    def f_tangent(self, t, y, args=None):
        '''Returns ODE RHS, with tangent
        parameters:
            t -- time in ODE
            y -- y(t) in ODE, with tangent
            arg1 -- parameter for the ODE

        return:
            the RHS of the ODE, with tangent
        '''
        raise NotImplementedError('A system class should implement member function f_tangent')
        return y

    
    def convert_coords(self, coord1):
        '''Converts coordinates (for example s,theta,zeta to R,Z,phi)
        parameters:
            coords1 - the coordinates to convert

        return:
            the converted coordinates
        '''
        raise NotImplementedError('A system class should implement member function convert_coords')
        return coord1
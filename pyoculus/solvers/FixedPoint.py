########################################
# FixPoint.py: class for finding fixed points
# written by @zhisong (zhisong.qu@anu.edu.au)
#

class FixedPoint:
    """
    Class that used to setup the fixed point finder.
    Call signature:
        MyFixPoint = FixedPoint(integrator, ) 

    Contains:
        bfield - Python wrapper for getting the magnetic field line ode
        bfield_tangent - Python wrapper for getting the magnetic field line ode, with tangent
        pjhfield - Python wrapper for getting the Pressure Jump Hamiltonian (PJH) ode
        pjhfield_tangent - Python wrapper for getting the Pressure Jump Hamiltonian (PJH) ode, with tangent

        get_xyz - Python wrapper for getting the xyz coordinates
    """
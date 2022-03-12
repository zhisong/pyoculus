## @file bfield_problem.py
#  @brief containing a problem class with magnetic fields
#  @author Zhisong Qu (zhisong.qu@anu.edu.au)
#

class BfieldProblem():

    def __init__(self):
        """! Set up the problem
        """
        pass


    def B(self, coords, args=None):
        """! Returns magnetic fields
        @param coordinates
        @param arg1 parameter
        @returns the contravariant magnetic fields
        """
        raise NotImplementedError("A problem class should implement member function B")

    def dBdX(self, coords, args=None):
        """! Returns magnetic fields
        @param coordinates
        @param arg1 parameter
        @returns B, dBdX, the contravariant magnetic fields, the derivatives of them
        """
        raise NotImplementedError(
            "A problem class should implement member function dBdX"
        )
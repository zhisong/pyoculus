## @file bfield_problem.py
#  @brief containing a problem class with magnetic fields
#  @author Zhisong Qu (zhisong.qu@anu.edu.au)
#

class BfieldProblem():

    def __init__(self):
        """! Set up the problem
        """
        ## if the output magnetic field contains the jacobian factor or not
        self.has_jacobian = False


    def B(self, coords, *args):
        """! Returns magnetic fields
        @param coordinates
        @param *args extra parameters
        @returns the contravariant magnetic fields
        """
        raise NotImplementedError("A problem class should implement member function B")

    def dBdX(self, coords, *args):
        """! Returns magnetic fields
        @param coordinates
        @param *args extra parameters
        @returns B, dBdX, the contravariant magnetic fields, the derivatives of them
        """
        raise NotImplementedError(
            "A problem class should implement member function dBdX"
        )

    def B_many(self, x1arr, x2arr, x3arr, input1D=True, *args):
        """! Returns magnetic fields, with multipy coordinate inputs
        @param x1arr the first coordinates. Should have the same length as the other two if input1D=True.
        @param x2arr the second coordinates. Should have the same length as the other two if input1D=True.
        @param x3arr the third coordinates. Should have the same length as the other two if input1D=True.
        @param input1D if False, create a meshgrid with sarr, tarr and zarr, if True, treat them as a list of points
        @param *args extra parameters
        @returns the contravariant magnetic fields
        """
        raise NotImplementedError("B_many is not implemented")

    def dBdX_many(self, x1arr, x2arr, x3arr, input1D=True, *args):
        """! Returns magnetic fields
        @param x1arr the first coordinates. Should have the same length as the other two if input1D=True.
        @param x2arr the second coordinates. Should have the same length as the other two if input1D=True.
        @param x3arr the third coordinates. Should have the same length as the other two if input1D=True.
        @param input1D if False, create a meshgrid with sarr, tarr and zarr, if True, treat them as a list of points
        @param *args extra parameters
        @returns B, dBdX, the contravariant magnetic fields, the derivatives of them
        """
        raise NotImplementedError(
            "dBdX_many is not implemented"
        )
## @file interpolate_coordinates.py: class for keeping interpolating coordinates between a few surfaces
#  @brief class for interpolating coordinates
#  @author Zhisong Qu (zhisong.qu@anu.edu.au)
#

import numpy as np
from scipy.special.orthogonal import jacobi

nax = np.newaxis


class SurfacesToroidal:
    """! Toroidal surfaces object
    Define surfaces \f$s(\vartheta, \zeta) \f$ and the angle transformations \f$ \theta (\vartheta, \zeta) \f$.
    \f$s, \vartheta\f$ are written in terms of Fourier harmonics
    \f[
        s(\vartheta, \zeta) = \sum_{m=0}^{\text{mpol}}\sum_{n=-\text{ntor}}^{\text{ntor}}
        s_{c,m,n} \cos(m \vartheta - n \text{N} \zeta) + s_{s,m,n} \sin(m \vartheta - n \text{N} \zeta),
    \f]
    where \f$N = \text{Nfp}\f$.
    \f[
        \theta(\vartheta, \zeta) = \vartheta + \sum_{m=0}^{\text{mpol}}\sum_{n=-\text{ntor}}^{\text{ntor}}
        \theta_{c,m,n} \cos(m \vartheta - n \text{N} \zeta) + \theta_{s,m,n} \sin(m \vartheta - n \text{N} \zeta).
    \f]
    If stellarator symmetry is assumed, then \f$s\f$ will only have cosine components and \f$\theta \f$ will only have sine components.

    Each surface is assigned a new radial label \f$\rho_i\f$.
    Now a new set of coordinates \f$(\rho, \vartheta, \zeta)\f$ will be constructed
    given the relationship \f$s(\rho, \vartheta, \zeta), \theta(\rho, \vartheta, \zeta)\f$.
    If \f$\rho\f$ is not on one of the known surfaces, an interpolation will be performed on all the Fourier harmonics radially using the method of the user's choice.
    """

    ## Output class
    class CoordsOutput:
        def __init__(self):
            pass

    def __init__(self, nsurfaces=2, mpol=10, ntor=10, Nfp=1, stellar_sym=True):
        """! Initialise the toroidal surface object
        @param mpol the poloidal mode number
        @param ntor the toroidal mode number
        @param Nfp the toroidal periodicity
        @param stellar_sym assuming stellarator symmetry or not
        """
        ## The poloidal mode number
        self.mpol = mpol
        ## The toroidal mode number
        self.ntor = ntor
        ## The field period
        self.Nfp = Nfp
        ## Stellarator symmetry
        self.sym = stellar_sym
        ## Number of surfaces
        self.nsurfaces = nsurfaces

        ## The new s coordinate for each surface
        self.ssurfs = np.linspace(0, 1, nsurfaces)

        ## The cosine coefficients of \f$s\f$, dimension (#interfaces, mpol, 2*ntor+1)
        self.scn = np.zeros([nsurfaces, mpol + 1, 2 * ntor + 1])
        ## The sine coefficients of \f$\vartheta\f$, dimension (#interfaces, mpol, 2*ntor+1)
        self.tsn = self.scn.copy()
        self.scn[0, 0, 0] = 0.0
        self.scn[-1, 0, 0] = 1.0

        if not self.sym:
            ## The sine coefficients of \f$s\f$, dimension (#interfaces, mpol, 2*ntor+1)
            self.ssn = self.tsn.copy()
            ## The cosine coefficients of \f$\vartheta\f$, dimension (#interfaces, mpol, 2*ntor+1)
            self.tcn = self.tsn.copy()

    def add_surface(self, surface, location=0):
        pass

    def get_coords(self, sarr=[0], tarr=[0], zarr=[0], derivative=0, input1D=True):
        """! Compute the coordinates and their derivatives given \f$\rho, \vartheta, \zeta\f$
        @param sarr \f$\rho\f$, an array
        @param tarr \f$\vartheta\f$, an array
        @param zarr \f$\zeta\f$, an array
        @param derivative  0:only return \f$s, \theta\f$, 1: also return the derivatives, 2: also return the second derivatives
        @param input1D if False, create a meshgrid with sarr, tarr and zarr, if True, treat them as a list of points
        @returns coords  coords.s, coords.t contains the coordinates \f$s,\theta\f$. coords.ds and dt contains the first derivatives. coords.dds and ddt contains second derivates. coords.jacobi contains the jacobi matrix \f$J = \partial(s,\theta,\zeta)/\partial(\rho,\vartheta,\zeta)\f$, coords.djacobi contains the derivative.
        """

        coords = self.CoordsOutput()

        sarr = np.atleast_1d(sarr)
        tarr = np.atleast_1d(tarr)
        zarr = np.atleast_1d(zarr)

        mtheta = self._mlist[:, nax] * tarr[nax, :]
        nzeta = self._nlist[:, nax] * zarr[nax, :]
        if input1D:
            alpha = mtheta[:, nax, :] - nzeta[nax, :, :]
            t0arr = tarr
            z0arr = zarr
        else:
            alpha = mtheta[:, nax, nax, :, nax] - nzeta[nax, :, nax, nax, :]
            t0arr = tarr[nax, :, nax]
            z0arr = zarr[nax, nax, :]

        cosalpha = np.cos(alpha)
        sinalpha = np.sin(alpha)
        scns = np.transpose(self._scn_int(sarr), [1, 2, 0])
        tsns = np.transpose(self._tsn_int(sarr), [1, 2, 0])
        if not self.sym:
            ssns = np.transpose(self._ssn_int(sarr), [1, 2, 0])
            tcns = np.transpose(self._tcn_int(sarr), [1, 2, 0])
        if not input1D:
            scns = scns[:, :, :, nax, nax]
            tsns = tsns[:, :, :, nax, nax]
            if not self.sym:
                ssns = ssns[:, :, :, nax, nax]
                tcns = tcns[:, :, :, nax, nax]

        s = np.sum(scns * cosalpha, axis=(0, 1))
        t = t0arr + np.sum(tsns * sinalpha, axis=(0, 1))

        if not self.sym:
            s += np.sum(ssns * sinalpha, axis=(0, 1))
            t += np.sum(tcns * cosalpha, axis=(0, 1))

        coords.s = s
        coords.t = t
        coords.z = np.broadcast_to(z0arr, coords.s.shape)

        if derivative > 0:
            dscns = np.transpose(self._scn_d1(sarr), [1, 2, 0])
            dtsns = np.transpose(self._tsn_d1(sarr), [1, 2, 0])
            if not self.sym:
                dssns = np.transpose(self._ssn_d1(sarr), [1, 2, 0])
                dtcns = np.transpose(self._tcn_d1(sarr), [1, 2, 0])
            if input1D:
                mlist = self._mlist[:, nax, nax]
                nlist = self._nlist[nax, :, nax]
            else:
                mlist = self._mlist[:, nax, nax, nax, nax]
                nlist = self._nlist[nax, :, nax, nax, nax]
                dscns = dscns[:, :, :, nax, nax]
                dtsns = dtsns[:, :, :, nax, nax]
                if not self.sym:
                    dssns = dssns[:, :, :, nax, nax]
                    dtcns = dtcns[:, :, :, nax, nax]

            msinalpha = mlist * sinalpha
            mcosalpha = mlist * cosalpha
            nsinalpha = nlist * sinalpha
            ncosalpha = nlist * cosalpha

            dsdr = np.sum(dscns * cosalpha, axis=(0, 1))
            dtdr = np.sum(dtsns * sinalpha, axis=(0, 1))
            dsdt = -np.sum(scns * msinalpha, axis=(0, 1))
            dtdt = 1 + np.sum(tsns * mcosalpha, axis=(0, 1))
            dsdz = np.sum(scns * nsinalpha, axis=(0, 1))
            dtdz = -np.sum(tsns * ncosalpha, axis=(0, 1))

            if not self.sym:
                dsdr += np.sum(dssns * sinalpha, axis=(0, 1))
                dtdr += np.sum(dtcns * cosalpha, axis=(0, 1))
                dsdt += np.sum(ssns * mcosalpha, axis=(0, 1))
                dtdt += -np.sum(tcns * msinalpha, axis=(0, 1))
                dsdz += -np.sum(ssns * ncosalpha, axis=(0, 1))
                dtdz += np.sum(tcns * nsinalpha, axis=(0, 1))

            zeros = np.zeros_like(coords.s)
            ones = np.ones_like(coords.t)

            coords.ds = np.stack([dsdr, dsdt, dsdz], -1)
            coords.dt = np.stack([dtdr, dtdt, dtdz], -1)
            coords.dz = np.stack([zeros, zeros, ones], -1)
            coords.jacobi = np.stack([coords.ds, coords.dt, coords.dz], -2)
            coords.jacobian = dsdr * dtdt - dtdr * dsdt

        if derivative > 1:
            ddscns = np.transpose(self._scn_d2(sarr), [1, 2, 0])
            ddtsns = np.transpose(self._tsn_d2(sarr), [1, 2, 0])
            if not self.sym:
                ddssns = np.transpose(self._ssn_d2(sarr), [1, 2, 0])
                ddtcns = np.transpose(self._tcn_d2(sarr), [1, 2, 0])
            if not input1D:
                ddscns = ddscns[:, :, :, nax, nax]
                ddtsns = ddtsns[:, :, :, nax, nax]
                if not self.sym:
                    ddssns = dssns[:, :, :, nax, nax]
                    ddtcns = dtcns[:, :, :, nax, nax]

            mmcosalpha = mlist ** 2 * cosalpha
            mmsinalpha = mlist ** 2 * sinalpha
            nncosalpha = nlist ** 2 * cosalpha
            nnsinalpha = nlist ** 2 * sinalpha
            mncosalpha = mlist * nlist * cosalpha
            mnsinalpha = mlist * nlist * sinalpha

            d2sdrdr = np.sum(ddscns * cosalpha, axis=(0, 1))
            d2tdrdr = np.sum(ddtsns * sinalpha, axis=(0, 1))
            d2sdrdt = -np.sum(dscns * msinalpha, axis=(0, 1))
            d2tdrdt = np.sum(dtsns * mcosalpha, axis=(0, 1))
            d2sdrdz = np.sum(dscns * nsinalpha, axis=(0, 1))
            d2tdrdz = -np.sum(dtsns * ncosalpha, axis=(0, 1))
            d2sdtdt = -np.sum(scns * mmcosalpha, axis=(0, 1))
            d2tdtdt = -np.sum(tsns * mmsinalpha, axis=(0, 1))
            d2sdzdz = -np.sum(scns * nncosalpha, axis=(0, 1))
            d2tdzdz = -np.sum(tsns * nnsinalpha, axis=(0, 1))
            d2sdtdz = np.sum(scns * mncosalpha, axis=(0, 1))
            d2tdtdz = np.sum(tsns * mnsinalpha, axis=(0, 1))

            if not self.sym:
                d2sdrdr += np.sum(ddscns * sinalpha, axis=(0, 1))
                d2tdrdr += np.sum(ddtsns * cosalpha, axis=(0, 1))
                d2sdrdt += np.sum(dssns * mcosalpha, axis=(0, 1))
                d2tdrdt += -np.sum(dtcns * msinalpha, axis=(0, 1))
                d2sdrdz += -np.sum(dssns * ncosalpha, axis=(0, 1))
                d2tdrdz += np.sum(dtcns * nsinalpha, axis=(0, 1))
                d2sdtdt += -np.sum(ssns * mmsinalpha, axis=(0, 1))
                d2tdtdt += -np.sum(tcns * mmcosalpha, axis=(0, 1))
                d2sdzdz += -np.sum(ssns * nnsinalpha, axis=(0, 1))
                d2tdzdz += -np.sum(tcns * nncosalpha, axis=(0, 1))
                d2sdtdz += np.sum(ssns * mnsinalpha, axis=(0, 1))
                d2tdtdz += np.sum(tcns * mncosalpha, axis=(0, 1))

            coords.dds = np.moveaxis(
                np.array(
                    [
                        [d2sdrdr, d2sdrdt, d2sdrdz],
                        [d2sdrdt, d2sdtdt, d2sdtdz],
                        [d2sdrdz, d2sdtdz, d2sdzdz],
                    ]
                ),
                (0, 1),
                (-2, -1),
            )
            coords.ddt = np.moveaxis(
                np.array(
                    [
                        [d2tdrdr, d2tdrdt, d2tdrdz],
                        [d2tdrdt, d2tdtdt, d2tdtdz],
                        [d2tdrdz, d2tdtdz, d2tdzdz],
                    ]
                ),
                (0, 1),
                (-2, -1),
            )

            coords.ddz = np.zeros_like(coords.dds)
            coords.djacobi = np.stack([coords.dds, coords.ddt, coords.ddz], -3)
            coords.djacobian = (
                coords.dds[..., 0] * dtdt[..., nax]
                + dsdr[..., nax] * coords.ddt[..., 1]
                - coords.ddt[..., 0] * dsdt[..., nax]
                - dtdr * coords.dds[..., 1]
            )

        return coords

    # def jacobian_transform(self, J, coords: CoordsOutput, derivative=False, dJ=None):
    #     """! Compute the coordinate transformation for the Jacobian
    #     @param the Jacobian \f$|J|\f$ of the old coordinate \f$(s,\theta,\zeta)\f$ to be transformed
    #     @param coords  the output of get_coords, should match the dimension of J
    #     @param derivative  if the derivatives are needed or not. If True, dJ is required.
    #     @param dJ  the derivative of \f$|J|\f$ wrt \f$(s,\theta,\zeta)\f$, have the dimension (3 derivatives,...)
    #     @returns Joutput, dJoutput the transformed Jacobian and the derivatives wrt \f$(\rho,\vartheta,\zeta)\f$
    #     """
    #     Joutput = J * coords.jacobian

        return Joutput

    def metric_transform(self, g, coords: CoordsOutput, derivative=False, dg=None):
        """! Compute the coordinate transformation for a metric tensor \f$g_{ij}\f$ (lower!)
        @param g  the metric \f$g_{ij}\f$ of the old coordinate \f$(s,\theta,\zeta)\f$ to be transformed, should have the dimension (3,3,...)
        @param coords  the output of get_coords, should match the dimension of g
        @param derivative  if the derivatives are needed or not. If True, dg is required.
        @param dg  the derivative of g wrt \f$(s,\theta,\zeta)\f$, have the dimension (3 derivatives, 3,3,...)
        @returns goutput, [dgoutput], the transformed metric and the derivative wrt \f$(\rho,\vartheta,\zeta)\f$

        Let the input metric for \f$(s, \theta, \zeta)\f$  being \f$\mathbf{g}\f$,
        the metric for \f$(\rho, \vartheta, \zeta)\f$ is given by
        \f[
            \mathbf{G} = \mathbf{J}^T \mathbf{g} \mathbf{J},
        \f]
        where
        \f[
            \mathbf{J} = \frac{\partial(s,\theta,\zeta)}{\partial(\rho,\vartheta,\zeta)}
        \f]
        is the Jacobi matrix.
        """

        goutput = np.moveaxis(coords.jacobi, -1, -2) @ g @ coords.jacobi

        return goutput

    def contra_vector_transform(
        self, v, coords: CoordsOutput, has_jacobian=False, derivative=False, dv=None
    ):
        """! Compute the coordinate transformation for a contravariant vector \f$v^i\f$ (upper!)
        @param v  the vector \f$v^i\f$ to be transformed, should have the dimension (3,...)
        @param coords  the output of get_coords, should match the dimension of v
        @param has_jacobian  if the given vector already contains the jacobian
        @param derivative  if the derivatives are needed or not. If True, dv is required.
        @param dv  the derivative of v wrt \f$(s,\theta,\zeta)\f$, have the dimension (3 derivatives, 3 component,...)
        @returns voutput, dvoutput, the transformed metric and the derivative wrt \f$(\rho,\vartheta,\zeta)\f$

        The vector transformation is given by
        \f[
            \left(\begin{array}{c}
            v^{\rho} \\ v^{\vartheta} \\ v^{\zeta}
            \end{array}\right)
            = \mathbf{J}^{-1}
            \left(\begin{array}{c}
            v^{{s}} \\ v^{\theta} \\ v^{\zeta}
            \end{array}\right),
        \f]
        where
        \f[
            \mathbf{J} = \frac{\partial(s,\theta,\zeta)}{\partial(\rho,\vartheta,\zeta)}
        \f]
        is the Jacobi matrix. For computing the result we use matrix solve instead of matrix inversion.

        When the derivative is needed, we note that
        \f[
            \frac{d}{d\rho}\mathbf{J}
            \left(\begin{array}{c}
            v^{\rho} \\ v^{\vartheta} \\ v^{\zeta}
            \end{array}\right)
            +\mathbf{J}
            \frac{d}{d\rho}\left(\begin{array}{c}
            v^{\rho} \\ v^{\vartheta} \\ v^{\zeta}
            \end{array}\right)
            =
            \frac{d}{d\rho}\left(\begin{array}{c}
            v^{{s}} \\ v^{\theta} \\ v^{\zeta}
            \end{array}\right),
        \f]
        therefore
        \f[
            \frac{d}{d\rho}\left(\begin{array}{c}
            v^{\rho} \\ v^{\vartheta} \\ v^{\zeta}
            \end{array}\right)
            =
            \mathbf{J}^{-1}
            \left[\frac{d}{d\rho}\left(\begin{array}{c}
            v^{{s}} \\ v^{\theta} \\ v^{\zeta}
            \end{array}\right)
            -\frac{d}{d\rho}\mathbf{J}
            \left(\begin{array}{c}
            v^{\rho} \\ v^{\vartheta} \\ v^{\zeta}
            \end{array}\right)\right].
        \f]
        Similarly we can compute the \f$\vartheta\f$ and \f$\zeta\f$ derivatives.

        """
        voutput = np.linalg.solve(coords.jacobi, v)

        if has_jacobian:
            voutput = voutput * coords.jacobian[...,nax]

        return voutput

    def construct_interpolant(self, method="cubic_spline", **kwargs):
        """! Construct the interpolant between surfaces for the given method
        @param method  The method being used, can be one of ['cubic_spline', 'cubic_hermite', 'pchip']. The corresponding scipy class is being used.
        @param **kwargs Other parameters passing onto the interpolant constructor
        """

        self._mlist = np.arange(0, self.mpol + 1)
        self._nlist = (
            np.concatenate([np.arange(0, self.ntor + 1), np.arange(-self.ntor, 0)])
            * self.Nfp
        )

        if method == "cubic_spline":
            from scipy.interpolate import CubicSpline

            ## The radial interpolant for s cosine components
            self._scn_int = CubicSpline(self.ssurfs, self.scn, axis=0, **kwargs)
            ## The radial interpolant for theta sine components
            self._tsn_int = CubicSpline(self.ssurfs, self.tsn, axis=0, **kwargs)
            if not self.sym:
                ## The radial interpolant for s sine components
                self._ssn_int = CubicSpline(self.ssurfs, self.ssn, axis=0, **kwargs)
                ## The radial interpolant for theta cosine components
                self._tcn_int = CubicSpline(self.ssurfs, self.tcn, axis=0, **kwargs)
        elif method == "cubic_hermite":  ## this is used in the original Hudson version
            from scipy.interpolate import CubicHermiteSpline

            dydx = np.zeros_like(self.scn)
            self._scn_int = CubicHermiteSpline(
                self.ssurfs, self.scn, dydx, axis=0, **kwargs
            )
            self._tsn_int = CubicHermiteSpline(
                self.ssurfs, self.tsn, dydx, axis=0, **kwargs
            )
            if not self.sym:
                self._ssn_int = CubicHermiteSpline(
                    self.ssurfs, self.ssn, dydx, axis=0, **kwargs
                )
                self._tcn_int = CubicHermiteSpline(
                    self.ssurfs, self.tcn, dydx, axis=0, **kwargs
                )
        elif method == "pchip":
            from scipy.interpolate import PchipInterpolator

            self._scn_int = PchipInterpolator(self.ssurfs, self.scn, axis=0, **kwargs)
            self._tsn_int = PchipInterpolator(self.ssurfs, self.tsn, axis=0, **kwargs)
            if not self.sym:
                self._ssn_int = PchipInterpolator(
                    self.ssurfs, self.ssn, axis=0, **kwargs
                )
                self._tcn_int = PchipInterpolator(
                    self.ssurfs, self.tcn, axis=0, **kwargs
                )
        else:
            raise ValueError("Unknown choice of interpolant")

        # construct the derivative and second derivative
        ## The radial interpolant for s cosine components, first derivative
        self._scn_d1 = self._scn_int.derivative(1)
        ## The radial interpolant for s cosine components, second derivative
        self._scn_d2 = self._scn_int.derivative(2)
        ## The radial interpolant for theta sine components, first derivative
        self._tsn_d1 = self._tsn_int.derivative(1)
        ## The radial interpolant for theta sine components, second derivative
        self._tsn_d2 = self._tsn_int.derivative(2)
        if not self.sym:
            ## The radial interpolant for s sine components,first derivative
            self._ssn_d1 = self._ssn_int.derivative(1)
            ## The radial interpolant for s sine components, second derivative
            self._ssn_d2 = self._ssn_int.derivative(2)
            ## The radial interpolant for theta cosine components, first derivative
            self._tcn_d1 = self._tcn_int.derivative(1)
            ## The radial interpolant for theta cosine components, second derivative
            self._tcn_d2 = self._tcn_int.derivative(2)

    def read_surfaces_from_file(self, filename, Nfp=None):
        """! Read the surfaces into a numpy array on disk
        @param filename  the filename to load from
        @param Nfp  the toroidal periodicity, if known.
        """
        surfaces = np.fromfile(filename)
        dims = surfaces.shape

        if dims[0] == 2:
            self.sym = True
            self.scn = surfaces[0, :]
            self.tsn = surfaces[1, :]
        elif dims[0] == 4:
            self.sym = False
            self.scn = surfaces[0, :]
            self.tsn = surfaces[1, :]
            self.ssn = surfaces[2, :]
            self.tcn = surfaces[3, :]
        else:
            raise ValueError(
                "The first dimension of the array should either be 2 (stellarator symmetry) or 4 (non-stellarator symmetry)"
            )

        # unfold the surface labels
        self.ssurfs = self.tsn[:, 0, 0]
        self.tsn[:, 0, 0] = 0

        self.nsurfaces = dims[1]
        self.mpol = dims[2] - 1
        self.ntor = (dims[3] - 1) // 2

        if not Nfp is None:
            self.Nfp = Nfp

    def write_surfaces_to_file(self, filename):
        """! Save the surfaces into a numpy array on disk
        @param filename  the filename to save to
        """

        if self.sym:
            surfaces = np.concatenate([self.scn[nax, :], self.tsn[nax, :]])
        else:
            surfaces = np.concatenate(
                [self.scn[nax, :], self.tsn[nax, :], self.ssn[nax, :], self.tcn[nax, :]]
            )
        # fold in the s labels, because the m=0,n=0 sine component is not used
        surfaces[1, :, 0, 0] = self.ssurfs
        surfaces.tofile(filename)

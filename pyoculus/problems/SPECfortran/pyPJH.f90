!> @file pyPJH.f90
!> @brief Fortran module for SPEC Pressure Jump Hamiltonian problems
!> @author Zhisong Qu (zhisong.qu@anu.edu.au)
!> @author Stuart Hudson (shudson@pppl.gov)

!> Fortran module for SPEC Pressure Jump Hamiltonian problems
!>
!> Interfaced to Python, see \ref pyoculus.problems.SPECPJH.SPECPJH
!>
!>
!> @todo non-stellarator symmetric cases
!>
!> ## Pressure Jump Hamiltonian (PJH) Fortran module
!>
!> Please note that in the Python module \ref pyoculus.problems.SPECPJH.SPECPJH, \f$\zeta\f$ is equivalent to \f$\varphi\f$
!>
!> ### Usage in Python:
!> To use the pjh module, one needs to initialize the fortran module first in Python (using class pyoculus.problems.SPECProblem or classes derived from it), 
!> then initialize pjh in Python by calling its Python wrap in `pyoculus.problems.SPECfortran.fortran_module`
!>
!>     init_pjh(dp, inside_or_outside, plus_or_minus)
!>
!> @param dp \f$\delta p\f$, the pressure jump
!> @param inside_or_outside for the specified volume, we compute things on the inner interface or outer
!> @param plus_or_minus whether to take the plus or minus sign in computing \f$p_\varphi\f$
!>
!> After initialization, one can call
!>
!>     rhs = get_pjhfield(phi, ptt)
!>
!> or for rhs and tangent,
!> 
!>     rhs = get_pjhfield_tangent(phi, ptt)
!>
!> @param phi the \f$\varphi\f$ angle
!> @param ptt the state vector (\f$p_\theta\f$, \f$\theta\f$), or (\f$p_\theta\f$, \f$\theta\f$, dp1, dt1, dp2, dt2)
!> @returns array of size 2 (for get_pjhfield) or array of size 6 (for get_pjhfield_tangent)
!>
!> To estimate the initial value for \f$p_\theta\f$, sometimes we need \f$B_\theta\f$ on the known side of the interface.
!> This can be obtained by calling
!>
!>     (Btheta, Bphi) = get_covariant_field(theta, phi)
!>
!> @param theta \f$\theta\f$
!> @param phi \f$\varphi\f$
!> @returns the covariant component on the known side of the interface, \f$(B_{1,\theta}, B_{1,\varphi})\f$
!>
MODULE SPECpjh

  USE SPECtypedefns, only : REAL_KIND

  PRIVATE

  INTEGER                          :: innout       !< specifying field on which side (-1) inner side, (+1) outer side
  INTEGER                          :: ioi          !< specifying which interface
  INTEGER                          :: plusminus    !< +1 or -1 sign taken for calculating \f$p_\varphi$\
  REAL(KIND=REAL_KIND)             :: delta_p      !< the jump of pressure

  REAL(KIND=REAL_KIND),ALLOCATABLE :: gBtmne(:)    !< the contravariant components of magnetic field on the interface
  REAL(KIND=REAL_KIND),ALLOCATABLE :: gBzmne(:)    !< the contravariant components of magnetic field on the interface
  REAL(KIND=REAL_KIND),ALLOCATABLE :: gBtmno(:)    !< the contravariant components of magnetic field on the interface
  REAL(KIND=REAL_KIND),ALLOCATABLE :: gBzmno(:)    !< the contravariant components of magnetic field on the interface

  PUBLIC :: get_pjhfield, get_pjhfield_tangent, init_PJH, destroy_PJH, get_b2_interface, get_covariant_field

  CONTAINS

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!> Compute the ODEs for PJH system, \f$dp_\theta/d\varphi\f$, \f$d\theta/d\varphi\f$
  SUBROUTINE get_pjhfield( phi , ptt , dptt )
!f2py threadsafe

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

    USE SPECconstants
    USE SPECvariables
    USE SPECcoords, ONLY: get_metric_interface

    IMPLICIT NONE

    REAL(KIND=REAL_KIND),INTENT(IN) :: phi      !< the \f$\varphi\f$ angle
    REAL(KIND=REAL_KIND),INTENT(IN) :: ptt(2)   !< (\f$p_\theta\f$, \f$\theta\f$)
    REAL(KIND=REAL_KIND),INTENT(OUT) :: dptt(2) !< (\f$dp_\theta/d\varphi\f$, \f$d\theta/d\varphi\f$)

    INTEGER :: imn 

    REAL(KIND=REAL_KIND) :: cosarg(mn), sinarg(mn),arg,pt,pp,theta
    REAL(KIND=REAL_KIND) :: dBB1(0:2),dBB2(0:3),a,b,c,discrim !< test variable

    REAL(KIND=REAL_KIND) :: dR(0:3,0:3,0:3),dZ(0:3,0:3,0:3),gl(3,3,0:2),sg(0:2),dG(0:2),dA(0:1),db1(0:2),dabc(3,0:2),dP(0:2),db2(0:2)
    REAL(KIND=REAL_KIND) :: ftftptpt,ftfpptpp,fpfppppp
    REAL(KIND=REAL_KIND) :: dtdot(0:2),dpdot(0:2),dsdot(0:2)
  
    pt = ptt(1) ! shorthand; "momentum" coordinate;

    theta = ptt(2) ! theta
  
!> - We need the metrics \f$g_{\theta \theta}, g_{\theta \varphi}, g_{\varphi \varphi}\f$ and the Jacobian \f$J\f$, as well as their first \f$\theta\f$ derivatives
!! by calling the subroutine coords::get_metric_interface
    dR=zero ; dZ=zero; gl=zero; sg=zero
    call get_metric_interface(ioi, theta, phi, dR, dZ, gl, sg, 1)
  
!> - Compute \f$ G \f$ and its derivatives.
    dG(0)= gl(2,2,0)*gl(3,3,0)                                                                           -         gl(2,3,0)**2
    dG(1)= gl(2,2,1)*gl(3,3,0)                           + gl(2,2,0)*gl(3,3,1)                           - two *   gl(2,3,0)*gl(2,3,1)
!    dG(2)= gl(2,2,2,2)*gl(3,3,0,0) + gl(2,2,2,0)*gl(3,3,2,0) + gl(2,2,2,0)*gl(3,3,2,0) + gl(2,2,0,0)*gl(3,3,2,2) - two * ( gl(2,3,2,0)*gl(2,3,2,0) + gl(2,3,0,0)*gl(2,3,2,2) )

!> - Compute the squared magnetic field strength \f$B_1^2\f$ on the other side of the interface by calling the subroutine specpjh::get_b2_interface
    call get_b2_interface(dBB1, ioi, theta, phi, gl, sg, 1)

!> - Construct \f$a, b, c\f$ for the quadratic equation for solving \f$p_\varphi\f$, then solve it by applying plus_or_minus sign we specified earlier (stored in dP)
    ! this is a, b, c for computing p_phi

    dabc(1,0) = gl(2,2,0)
    dabc(2,0) = -two * gl(2,3,0) * pt
    dabc(3,0) = gl(3,3,0) * pt**2 - (dBB1(0) + two * delta_p) * dG(0)

    discrim=dabc(2,0)**2-four*dabc(1,0)*dabc(3,0)
   
    if( discrim.lt.zero ) stop "ph00aa :          : WF=? ; discrim.lt.zero ;"
 
    dP(0)  = ( -dabc(2,0) + plusminus *                                                                                   sqrt (discrim)                       ) / (two*dabc(1,0))

!> - Computes \f$b_2\f$
    db2(0)  = gl(3,3,0)*   pt**2 - two*gl(2,3,0)* pt*dP(0)            + gl(2,2,0)*    dP(0)**2

!> - Construct the denominator of \f$d\theta/d\varphi\f$

    dtdot(0) = gl(3,3,0)*two*pt - two*gl(2,3,0)*dP(0)
    dpdot(0) = -two*gl(2,3,0)*pt + gl(2,2,0)*two*dP(0)

!> - Construct the numerators

    ftftptpt= -pt**2 ; ftfpptpp= -pt*dP(0) ; fpfppppp=-dP(0)**2

    dsdot(0) = gl(3,3,1)* ftftptpt                 - two*gl(2,3,1)* ftfpptpp                                          + gl(2,2,1)*fpfppppp                              &
             + dBB1(1)*dG(0) + db2(0)*dG(1)/dG(0)
 
!> - Get right hand side

    dptt(1:2) = (/ dsdot(0) , dtdot(0) /) / dpdot(0)
 
!> .
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  END SUBROUTINE get_pjhfield


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!> Compute the ODEs for PJH system, \f$dp_\theta/d\varphi\f$, \f$d\theta/d\varphi\f$, with tangent
  SUBROUTINE get_pjhfield_tangent( phi , ptt , dptt )
!f2py threadsafe

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

    USE SPECconstants
    USE SPECvariables
    USE SPECcoords, ONLY: get_metric_interface

    IMPLICIT NONE

    REAL(KIND=REAL_KIND),INTENT(IN) :: phi      !< the \f$\varphi\f$ angle
    REAL(KIND=REAL_KIND),INTENT(IN) :: ptt(6)   !< (\f$p_\theta\f$, \f$\theta\f$, \f$\Delta p_{\theta,1}\f$, \f$\Delta \theta_1\f$, \f$\Delta p_{\theta,2}\f$, \f$\Delta \theta_2\f$)
    REAL(KIND=REAL_KIND),INTENT(OUT) :: dptt(6) !< (\f$dp_\theta/d\varphi\f$, \f$d\theta/d\varphi\f$)

    INTEGER :: imn 

    REAL(KIND=REAL_KIND) :: cosarg(mn), sinarg(mn),arg,pt,pp,theta
    REAL(KIND=REAL_KIND) :: dBB1(0:2),dBB2(0:3),a,b,c,discrim !< test variable

    REAL(KIND=REAL_KIND) :: dR(0:3,0:3,0:3),dZ(0:3,0:3,0:3),gl(3,3,0:2),sg(0:2),dG(0:2),dA(0:1),db1(0:2),dabc(3,0:2),dP(0:2),db2(0:2)
    REAL(KIND=REAL_KIND) :: ftftptpt,ftfpptpp,fpfppppp
    REAL(KIND=REAL_KIND) :: dtdot(0:2),dpdot(0:2),dsdot(0:2),TM(2,2)
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

    pt = ptt(1) ! shorthand; "momentum" coordinate;

    theta = ptt(2) ! theta
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

!> - We need the metrics \f$g_{\theta \theta}, g_{\theta \varphi}, g_{\varphi \varphi}\f$ and the Jacobian \f$J\f$, as well as their first \f$\theta\f$ derivatives
!! by calling the subroutine coords::get_metric_interface
    dR=zero ; dZ=zero; gl=zero; sg=zero
    call get_metric_interface(ioi, theta, phi, dR, dZ, gl, sg, 2)
  
!> - Compute \f$ G \f$ and its derivatives.
    dG(0)= gl(2,2,0)*gl(3,3,0)                                                                           -         gl(2,3,0)**2
    dG(1)= gl(2,2,1)*gl(3,3,0)                           + gl(2,2,0)*gl(3,3,1)                           - two *   gl(2,3,0)*gl(2,3,1)
    dG(2)= gl(2,2,2)*gl(3,3,0) + gl(2,2,1)*gl(3,3,1)     + gl(2,2,1)*gl(3,3,1) + gl(2,2,0)*gl(3,3,2)     - two * ( gl(2,3,1)*gl(2,3,1) + gl(2,3,0)*gl(2,3,2) )

!> - Compute the squared magnetic field strength \f$B_1^2\f$ on the other side of the interface by calling the subroutine specpjh::get_b2_interface
    call get_b2_interface(dBB1, ioi, theta, phi, gl, sg, 2)

!> - Construct \f$a, b, c\f$ for the quadratic equation for solving \f$p_\varphi\f$, then solve it by applying plus_or_minus sign we specified earlier (stored in dP)
    ! this is a, b, c for computing p_phi

    dabc(1,0) = gl(2,2,0)
    dabc(1,1) = gl(2,2,1)
 
    dabc(2,0) = -two * gl(2,3,0) * pt
    dabc(2,1) = -two * gl(2,3,1) * pt
    dabc(2,2) = -two * gl(2,3,0)
 
    dabc(3,0) = gl(3,3,0) * pt**2 - (dBB1(0) + two * delta_p) * dG(0)
    dabc(3,1) = gl(3,3,1) * pt**2 - (dBB1(0) + two * delta_p) * dG(1) - dBB1(1) * dG(0)
    dabc(3,2) = gl(3,3,2) * pt*two
   
    discrim=dabc(2,0)**2-four*dabc(1,0)*dabc(3,0)
   
    if( discrim.lt.zero ) stop "ph00aa :          : WF=? ; discrim.lt.zero ;"
 
    dP(0)  = ( -dabc(2,0) + plusminus *                                                                                   sqrt (discrim)                       ) / (two*dabc(1,0))
    dP(1)  = ( -dabc(2,1) + plusminus * half * (two*dabc(2,0)*dabc(2,1)-four*(dabc(1,1)*dabc(3,0)+dabc(1,0)*dabc(3,1))) / sqrt (discrim) - dP(0)*two*dabc(1,1) ) / (two*dabc(1,0))
    dP(2)  = ( -dabc(2,2) + plusminus * half * (two*dabc(2,0)*dabc(2,2)-four*(                    dabc(1,0)*dabc(3,2))) / sqrt (discrim)                       ) / (two*dabc(1,0))
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

!> - Computes \f$b_2\f$
    db2(0)  = gl(3,3,0)*   pt**2 - two*gl(2,3,0)* pt*dP(0)            + gl(2,2,0)*    dP(0)**2
    db2(1)  = gl(3,3,1)*   pt**2 - two*gl(2,3,1)* pt*dP(0)            + gl(2,2,1)*    dP(0)**2 & 
                                 - two*gl(2,3,0)* pt*dP(1)            + gl(2,2,0)*two*dP(0)*dP(1)
    db2(2)  = gl(3,3,0)*two*pt   - two*gl(2,3,0)*(   dP(0)+pt*dP(2) ) + gl(2,2,0)*two*dP(0)*dP(2)

!> - Construct the denominator of \f$d\theta/d\varphi\f$

    dtdot(0) = gl(3,3,0)*two*pt - two*gl(2,3,0)*dP(0)
    dtdot(1) = gl(3,3,1)*two*pt - two*gl(2,3,1)*dP(0) - two*gl(2,3,0)*dP(1)
    dtdot(2) = gl(3,3,0)*two                          - two*gl(2,3,0)*dP(2)
 
    dpdot(0) = -two*gl(2,3,0)*pt + gl(2,2,0)*two*dP(0)
    dpdot(1) = -two*gl(2,3,1)*pt + gl(2,2,1)*two*dP(0) + gl(2,2,0)*two*dP(1)
    dpdot(2) = -two*gl(2,3,0)                          + gl(2,2,0)*two*dP(2)

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

!> - Construct the numerators

    ftftptpt= -pt**2 ; ftfpptpp= -pt*dP(0) ; fpfppppp=-dP(0)**2

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

    dsdot(0) = gl(3,3,1)* ftftptpt                 - two*gl(2,3,1)* ftfpptpp                                          + gl(2,2,1)*fpfppppp                              &
             + dBB1(1)*dG(0) + db2(0)*dG(1)/dG(0)
 
!    dsdot(1) = gl(3,3,2,2)* ftftptpt                 - two*gl(2,3,2,2)* ftfpptpp                                          + gl(2,2,2,2)*fpfppppp                              &
!             + gl(3,3,2,0)*(two*df(2,0,0)*df(2,2,0)) - two*gl(2,3,2,0)*(df(2,2,0)*df(3,0,0)+df(2,0,0)*df(3,2,0)-pt*dP(1)) + gl(2,2,2,0)*two*(df(3,0,0)*df(3,2,0)-dP(0)*dP(1)) & 
!             + dA(1) + (db2(1)-db1(1))*dG(1)/dG(0) + (db2(0)-db1(0))*(dG(2)-dG(1)**2/dG(0))/dG(0)
 
!    dsdot(2) = gl(3,3,2,0)*(-two*pt)                 - two*gl(2,3,2,0)*(-dP(0)-pt*dP(2))                                  + gl(2,2,2,0)*(-two*dP(0)*dP(2))                    &
!                     +  db2(2)        *dG(1)/dG(0) 

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

!> - Get right hand side

    dptt(1:2) = (/ dsdot(0) , dtdot(0) /) / dpdot(0)
 
    TM(1,1) = ( dsdot(2) - dptt(1)*dpdot(2) ) / dpdot(0) ! d radial / d radial ;
    TM(1,2) = ( dsdot(1) - dptt(1)*dpdot(1) ) / dpdot(0) ! d radial / d angle  ;
    TM(2,1) = ( dtdot(2) - dptt(2)*dpdot(2) ) / dpdot(0) ! d angle  / d radial ;
    TM(2,2) = ( dtdot(1) - dptt(2)*dpdot(1) ) / dpdot(0) ! d angle  / d angle  ;
 
    dptt(3)=TM(1,1)*ptt(3) + TM(1,2)*ptt(5) ! standard format;
    dptt(4)=TM(1,1)*ptt(4) + TM(1,2)*ptt(6)
    dptt(5)=TM(2,1)*ptt(3) + TM(2,2)*ptt(5)
    dptt(6)=TM(2,1)*ptt(4) + TM(2,2)*ptt(6)
!> .
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  END SUBROUTINE get_pjhfield_tangent

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!> Get the covariant component of known side of magnetic field \f$\bf B_1\f$, i.e. \f$B_{1,\theta}, B_{1,\varphi}\f$.
!> Can be used to estimate the initial condition for \f$p_\theta\f$. 
  SUBROUTINE get_covariant_field(theta, phi, Bco)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

    USE SPECconstants
    USE SPECvariables
    USE SPECcoords, ONLY: get_metric_interface

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

    IMPLICIT NONE

    REAL(KIND=REAL_KIND), INTENT(IN) :: theta  !< the \f$\theta\f$ angle
    REAL(KIND=REAL_KIND), INTENT(IN) :: phi    !< the \f$\phi\f$ angle
    REAL(KIND=REAL_KIND), INTENT(OUT):: Bco(2) !< the covariant components (\f$B_{1,\theta}, B_{1,\varphi}\f$)

    REAL(KIND=REAL_KIND) :: cosarg(mn), sinarg(mn),pt,pp
    REAL(KIND=REAL_KIND) :: dR(0:3,0:3,0:3),dZ(0:3,0:3,0:3),gl(3,3,0:2),sg(0:2),gB(2:3,0:2)

    call get_metric_interface(ioi, theta, phi, dR, dZ, gl, sg, 1)

    cosarg = COS(im*theta-in1*phi)
    sinarg = SIN(im*theta-in1*phi)

    gB(2,0) = SUM(gBtmne * cosarg)
    gB(3,0) = SUM(gBzmne * cosarg)

    Bco(1) = gB(2,0) * gl(2,2,0) + gB(3,0) * gl(2,3,0)
    Bco(2) = gB(2,0) * gl(2,3,0) + gB(3,0) * gl(3,3,0)

    Bco = Bco / sg(0)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  END SUBROUTINE get_covariant_field

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!> An internal subroutine that computes the squared magnetic field strength on the interface \f$B_1^2\f$ and its first and second \f$\theta\f$ derivatives.
!> This is a private subroutine that is not interfaced with Python.
  SUBROUTINE get_b2_interface(db2, ioi, theta, phi, guvij, sg, ideriv)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

    USE SPECtypedefns, ONLY : REAL_KIND
    USE SPECconstants
    USE SPECvariables

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

    IMPLICIT NONE
    REAL(KIND=REAL_KIND),    INTENT(IN) :: theta, phi ! the three coordinates (s,theta,varphi)
    INTEGER             ,    INTENT(IN) :: ioi ! on which interface
    INTEGER             ,    INTENT(IN) :: ideriv ! the level of derivatives. 0 for none, 1 for first, 2 for second

    REAL(KIND=REAL_KIND),    INTENT(IN):: guvij(3,3,0:2), sg(0:2)

    REAL(KIND=REAL_KIND),    INTENT(OUT):: db2(0:2)

    REAL(KIND=REAL_KIND) :: lss, rsign
    REAL(KIND=REAL_KIND) :: cosarg(mn), sinarg(mn), gB(2:3,0:2)

    INTEGER :: ii, jj

    db2 = zero

    cosarg = COS(im*theta-in1*phi)
    sinarg = SIN(im*theta-in1*phi)

    gB(2,0) = SUM(gBtmne * cosarg)
    gB(3,0) = SUM(gBzmne * cosarg)

    DO ii = 2, 3
      DO jj = 2, 3
        db2(0) = db2(0) + guvij(ii,jj,0) * gB(ii,0) * gB(jj,0)
      ENDDO
    ENDDO
    db2(0) = db2(0) / sg(0)**2

    IF (ideriv .GE. 1) THEN
      gB(2,1) = SUM(gBtmne * (-im) * sinarg)
      gB(3,1) = SUM(gBzmne * (-im) * sinarg)

      DO ii = 2, 3
        DO jj = 2, 3
          db2(1) = db2(1) + guvij(ii,jj,1) * gB(ii,0) * gB(jj,0) + guvij(ii,jj,0) * (gB(ii,1) * gB(jj,0) + gB(ii,0) * gB(jj,1))
        ENDDO
      ENDDO
      db2(1) = db2(1) / sg(0)**2 - 2 * db2(0) * sg(1) / sg(0)
    ENDIF

    IF (ideriv .GE. 2) THEN

      gB(2,2) = SUM(gBtmne * (-im) * (im) * cosarg)
      gB(3,2) = SUM(gBzmne * (-im) * (im) * cosarg)

      DO ii = 2, 3
        DO jj = 2, 3
          db2(2) = db2(2) + guvij(ii,jj,2) * gB(ii,0) * gB(jj,0) + guvij(ii,jj,0) * (gB(ii,2) * gB(jj,0) + gB(ii,0) * gB(jj,2) + two * gB(ii,1) * gB(jj,1)) &
                      +two* guvij(ii,jj,1) * (gB(ii,1) * gB(jj,0) + gB(ii,0) * gB(jj,1)) 
        ENDDO
      ENDDO
      db2(2) = db2(2) / sg(0)**2 - 2 * db2(0) * sg(1) / sg(0)

    ENDIF


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  END SUBROUTINE get_b2_interface

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!> Initialize the PJH system.
!> Whenever \f$\delta p\f$ is modified, you will need to call this subroutine (in Python interface) again.
  SUBROUTINE init_pjh(dp, inside_or_outside, plus_or_minus)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
    USE SPECconstants
    USE SPECvariables
    USE SPECbasefn

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

    IMPLICIT NONE

    REAL(KIND=REAL_KIND), INTENT(IN) :: dp   !< the pressure jump \f$\delta p\f$
    INTEGER, INTENT(IN) :: inside_or_outside !< 0 for inside interface of the selected volume ivol (set in SPECvariables::ivol), +1 for outside
    INTEGER, INTENT(IN) :: plus_or_minus     !< set the sign (-1 or +1) for computing \f$p_\varphi\f$ using the quadratic root equation.

    INTEGER :: ii, jj, ll, mi, ni

    REAL(KIND=REAL_KIND) :: sbar, lss
    REAL(KIND=REAL_KIND) :: cheby(0:Lrad,0:1), zernike(0:Lrad,0:Mpol,0:1)
    
    REAL(KIND=REAL_KIND) :: TT(0:Lrad,0:1)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

    plusminus = plus_or_minus
    delta_p = dp
    ioi = ivol + inside_or_outside       ! the index of the interface

    IF (ALLOCATED(gBtmne)) THEN
      DEALLOCATE(gBtmne)
    ENDIF
    IF (ALLOCATED(gBzmne)) THEN
      DEALLOCATE(gBzmne)
    ENDIF
    ALLOCATE(gBtmne(1:mn))
    ALLOCATE(gBzmne(1:mn))

    IF (NOTstellsym) THEN
      IF (ALLOCATED(gBtmno)) THEN
        DEALLOCATE(gBtmno)
      ENDIF
      IF (ALLOCATED(gBzmno)) THEN
        DEALLOCATE(gBzmno)
      ENDIF
      ALLOCATE(gBtmno(1:mn))
      ALLOCATE(gBzmno(1:mn))
    ENDIF

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!> The contravariant component of \f$\bf B_1\f$ is computed on the interface.
!> The Fourier components of \f$ J B^{\theta}_1\f$ and \f$ J B^{\varphi}_1\f$ are stored for further use.

    IF (Lcoordinatesingularity) THEN
      sbar = FLOAT(inside_or_outside)
      CALL get_zernike(sbar, Lrad, Mpol, zernike(:,:,0:1))
    ELSE
      lss = FLOAT(inside_or_outside) * two - one
      CALL get_cheby(lss, Lrad, cheby(0:Lrad,0:1))
    END IF


    DO ii = 1, mn ; mi = im(ii) ; ni = in1(ii) 


      IF( Lcoordinatesingularity ) THEN ! regularization factor depends on mi; 17 Dec 15;

        DO ll = 0, Lrad ; TT(ll,0:1) = (/        zernike(ll,mi,0),        zernike(ll,mi,1)*half                      /)
        ENDDO

      ELSE

        DO ll = 0, Lrad ; TT(ll,0:1) = (/             cheby(ll,0),             cheby(ll,1)                           /)
        ENDDO

      ENDIF ! end of IF( Lcoordinatesingularity ) ; 16 Jan 15;
  
      gBtmne(ii) = SUM( -Aze(1:Lrad+1,ii) * TT(0:Lrad, 1) )
      gBzmne(ii) = SUM( +Ate(1:Lrad+1,ii) * TT(0:Lrad, 1) )

      IF (NOTstellsym) THEN
        gBtmno(ii) = SUM( -Azo(1:Lrad+1,ii) * TT(0:Lrad, 1) )
        gBzmno(ii) = SUM( +Ato(1:Lrad+1,ii) * TT(0:Lrad, 1) )
      END IF

    ENDDO

  END SUBROUTINE init_pjh

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!> This subroutine frees the memory space allocated by specpjh::init_pjh.
  SUBROUTINE destroy_pjh()

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

    USE SPECvariables

    IMPLICIT NONE

    IF (ALLOCATED(gBtmne)) THEN
      DEALLOCATE(gBtmne)
    ENDIF
    IF (ALLOCATED(gBzmne)) THEN
      DEALLOCATE(gBzmne)
    ENDIF

    IF (NOTstellsym) THEN
      IF (ALLOCATED(gBtmno)) THEN
        DEALLOCATE(gBtmno)
      ENDIF
      IF (ALLOCATED(gBzmno)) THEN
        DEALLOCATE(gBzmno)
      ENDIF
    ENDIF

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  END SUBROUTINE destroy_pjh

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

END MODULE SPECpjh

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
! Original ph00aa.h subroutine
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \item Pressure-jump Hamiltonian o.d.e.'s, and tangent map,

!latex \subsubheading{Pressure jump Hamiltonian}
!latex \item The pressure-jump Hamiltonian is derived from the force-balance condition at the ideal interfaces.
!latex Let $p_1$ and ${\bf B}_1$ be the pressure and field immediately inside the interface, and $p_1$ and ${\bf B}_1$ be the pressure and field immediately outside, 
!latex then the force balance condition requires that
!latex \be H \equiv 2 \, \delta p =  2 ( p_1 - p_2 ) = B_2^2 - B_1^2
!latex \ee
!late be a constant.
!latex For Beltrami fields, which satisfy $\nabla \times {\bf B}=\mu {\bf B}$, the magnitude of the field, $B$, on the interface (where we assume that $B^s=0$) may be written
!latex \be 
!latex B^2 = \frac{g_{\varphi\varphi} f_\theta f_\theta - 2 g_{\theta\varphi}f_\theta f_\varphi + g_{\theta\theta} f_\varphi f_\varphi}{g_{\theta\theta}g_{\varphi\varphi}-g_{\theta\varphi}g_{\theta\varphi}}
!latex \ee
!latex where $f$ is a surface potential and $g_{\theta\theta}$, $g_{\theta\varphi}$ and $g_{\varphi\varphi}$ are metric elements local to the interface.
!latex \item Assuming that the field $B_1$ is known on the `inside' of the interface, ie. $B_{1\theta}=f_\theta$, $B_{1\varphi}=f_\varphi$ and $f$ is known, 
!latex it is required to determine the tangential field, $p_\theta = B_\theta$ and $p_\varphi = B_\varphi$, on the `outside' of the interface.
!latex \item The o.d.e.'s are given by Hamilton's equations \be
!latex \dot \theta   =  \frac{\partial H}{\partial p_\theta}\Big|_{\theta,\varphi,p_\varphi}, \;\;
!latex \dot p_\theta = -\frac{\partial H}{\partial \theta}\Big|_{p_\theta,\varphi,p_\varphi}, \;\;
!latex \dot \varphi     =  \frac{\partial H}{\partial p_\varphi}\Big|_{\theta,p_\theta,\varphi}, \;\;
!latex \dot p_\varphi   = -\frac{\partial H}{\partial \varphi}\Big|_{\theta,p_\theta,p_\varphi}, \ee
!latex where the `dot' denotes derivative with respect to `time'.
!latex \item This is reduced to a $1\frac{1}{2}$ dimensional system by using $\varphi$ as the time-like integration parameter, and replacing the equation for $\dot p_\varphi$ with 
!latex \be p_\varphi= P(\theta,p_\theta,\varphi; \delta p) = \frac{-b\pm\sqrt{b^2-4ac}}{2a} \label{eq:pphi} \ee
!latex where \mbox{$a=g_{\theta\theta}$}, \mbox{$b=-2 g_{\theta\varphi}p_\theta$} 
!latex and \mbox{$c=g_{\varphi\varphi} p_{\theta}^2 - b_1 - 2 \, \delta p \, G$} (see below for definition of $b_1$ and $G$).
!latex \item The o.d.e.'s that then need to be followed are (see below for definition of $A$ and $b_2$) \be
!latex \frac{d   \theta}{d\varphi}&=& \frac{g_{\varphi\varphi} p_{\theta} - g_{\theta\varphi} p_{\varphi}}{-g_{\theta\varphi}p_{\theta}+g_{\theta\theta}p_{\varphi}},\\
!latex \frac{d p_\theta}{d\varphi}&=& \frac{g_{\varphi\varphi,\theta} (f_{\theta}^2-p_{\theta}^2)-2g_{\theta\varphi,\theta}(f_{\theta}f_{\varphi}-p_{\theta}p_{\varphi})
!latex +g_{\theta\theta,\theta}(f_{\varphi}^2-p_{\varphi}^2)+A+(b_2-b_1)G_{,\theta} / G }
!latex {-2g_{\theta\varphi}p_\theta+g_{\theta\theta}2p_{\varphi}}. \ee

!latex \item Note that $d\theta / d \varphi = B^\theta / B^\varphi $; there is a fundamental relation between the pressure-jump Hamiltonian and the field-line Hamiltonian. 
!latex (Furthermore, in many cases the surface will be given in straight field line coordinates, so $d \theta / d\varphi = const.$.)

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

! subroutine ph00aa( phi , ptt , dptt ) ! fixed format;

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
!   use SPECconstants
!   use fileunits, only : ounit
!   use inputlist, only : Wph00aa, pressure, pscale
!   use allglobal, only : myid, ncpu, ivol, cpus, mn, im, in, iRbc, iZbs

!   use pjhamilto, only : Itor,Gpol,spmn,io

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

!   LOCALS

!   REAL,intent(in) :: phi,ptt(6)
!   REAL,intent(out) :: dptt(6)

!   INTEGER,parameter :: plusminus=1

!   REAL :: delta_p

!   INTEGER :: imn,ioi

!   REAL :: cs(2),arg,pt,pp
!   REAL :: dBB1(0:2),dBB2(0:3),a,b,c,discrim

!   REAL :: dR(0:3,0:3,0:3),dZ(0:3,0:3,0:3),df(0:3,0:3,0:3),gl(3,3,0:3,0:3),dG(0:2),dA(0:1),db1(0:3),dabc(3,0:2),dP(0:2),db2(0:2)
!   REAL :: ftftptpt,ftfpptpp,fpfppppp
!   REAL :: dtdot(0:2),dpdot(0:2),dsdot(0:2),TM(2,2)
  
!   BEGIN(ph00aa)
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
! !latex \item The conserved `energy' $\delta p$ is given by the jump in pressure across the interface.

!   if( io.eq.0 ) delta_p = ( pressure(ivol) - pressure(ivol-1) ) * pscale
!   if( io.eq.1 ) delta_p = ( pressure(ivol) - pressure(ivol+1) ) * pscale
  
!   ioi = ivol+io-1

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
!   pt = ptt(1) ! shorthand; "momentum" coordinate;
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
!   dR=zero ; dZ=zero ; df=zero ; df(2,0,0)=Itor(ivol,io) ; df(3,0,0)=Gpol(ivol,io)
  
!   do imn = 1,mn ; arg = im(imn)*ptt(2)-in(imn)*phi ; cs = (/cos(arg),sin(arg)/)
   
!    dR(0,0,0) = dR(0,0,0) + iRbc(ioi,imn)*cs(1)
   
!    dR(2,0,0) = dR(2,0,0) + iRbc(ioi,imn)*cs(2)*(-im(imn))
!    dR(3,0,0) = dR(3,0,0) + iRbc(ioi,imn)*cs(2)*(+in(imn))
   
!    dR(2,2,0) = dR(2,2,0) + iRbc(ioi,imn)*cs(1)*(-im(imn))*(+im(imn))
!    dR(3,2,0) = dR(3,2,0) + iRbc(ioi,imn)*cs(1)*(+in(imn))*(+im(imn))
   
!    dR(2,2,2) = dR(2,2,2) + iRbc(ioi,imn)*cs(2)*(-im(imn))*(+im(imn))*(-im(imn))
!    dR(3,2,2) = dR(3,2,2) + iRbc(ioi,imn)*cs(2)*(+in(imn))*(+im(imn))*(-im(imn))
   
! !  dZ(0,0,0) = dZ(0,0,0) + iZbs(ioi,imn)*cs(2)
   
!    dZ(2,0,0) = dZ(2,0,0) + iZbs(ioi,imn)*cs(1)*(+im(imn))
!    dZ(3,0,0) = dZ(3,0,0) + iZbs(ioi,imn)*cs(1)*(-in(imn))
   
!    dZ(2,2,0) = dZ(2,2,0) + iZbs(ioi,imn)*cs(2)*(+im(imn))*(-im(imn))
!    dZ(3,2,0) = dZ(3,2,0) + iZbs(ioi,imn)*cs(2)*(-in(imn))*(-im(imn))
   
!    dZ(2,2,2) = dZ(2,2,2) + iZbs(ioi,imn)*cs(1)*(+im(imn))*(-im(imn))*(+im(imn))
!    dZ(3,2,2) = dZ(3,2,2) + iZbs(ioi,imn)*cs(1)*(-in(imn))*(-im(imn))*(+im(imn))
   
! !df(0,0,0) = df(0,0,0)+spmn(ivol,io,imn)*cs(2)
   
!    df(2,0,0) = df(2,0,0)+spmn(ivol,io,imn)*cs(1)*(+im(imn))
!    df(3,0,0) = df(3,0,0)+spmn(ivol,io,imn)*cs(1)*(-in(imn))

!    df(2,2,0) = df(2,2,0)+spmn(ivol,io,imn)*cs(2)*(+im(imn))*(-im(imn))
!    df(3,2,0) = df(3,2,0)+spmn(ivol,io,imn)*cs(2)*(-in(imn))*(-im(imn))
   
!    df(2,2,2) = df(2,2,2)+spmn(ivol,io,imn)*cs(1)*(+im(imn))*(-im(imn))*(+im(imn))
!    df(3,2,2) = df(3,2,2)+spmn(ivol,io,imn)*cs(1)*(-in(imn))*(-im(imn))*(+im(imn))
   
!   enddo
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
 
!   gl(2,2,0,0)=         dR(2,0,0)*dR(2,0,0)                       + dZ(2,0,0)*dZ(2,0,0)
!   gl(2,2,2,0)= two * ( dR(2,0,0)*dR(2,2,0)                       + dZ(2,0,0)*dZ(2,2,0))
!   gl(2,2,2,2)= two * ( dR(2,2,0)*dR(2,2,0) + dR(2,0,0)*dR(2,2,2) + dZ(2,2,0)*dZ(2,2,0) + dZ(2,0,0)*dZ(2,2,2))

!   gl(2,3,0,0)= dR(2,0,0)*dR(3,0,0)                                                                   + dZ(2,0,0)*dZ(3,0,0)
!   gl(2,3,2,0)= dR(2,2,0)*dR(3,0,0)                       + dR(2,0,0)*dR(3,2,0)                       + dZ(2,2,0)*dZ(3,0,0)                       + dZ(2,0,0)*dZ(3,2,0)
!   gl(2,3,2,2)= dR(2,2,2)*dR(3,0,0) + dR(2,2,0)*dR(3,2,0) + dR(2,2,0)*dR(3,2,0) + dR(2,0,0)*dR(3,2,2) &
!              + dZ(2,2,2)*dZ(3,0,0) + dZ(2,2,0)*dZ(3,2,0) + dZ(2,2,0)*dZ(3,2,0) + dZ(2,0,0)*dZ(3,2,2)

!   gl(3,3,0,0)=         dR(3,0,0)*dR(3,0,0)                       + dR(0,0,0)*dR(0,0,0)                       + dZ(3,0,0)*dZ(0,0,0)
!   gl(3,3,2,0)= two * ( dR(3,0,0)*dR(3,2,0)                       + dR(0,0,0)*dR(2,0,0)                       + dZ(3,0,0)*dZ(3,2,0)                       )
!   gl(3,3,2,2)= two * ( dR(3,2,0)*dR(3,2,0) + dR(3,0,0)*dR(3,2,2) + dR(2,0,0)*dR(2,0,0) + dR(0,0,0)*dR(2,2,0) + dZ(3,2,0)*dZ(3,2,0) + dZ(3,0,0)*dZ(3,2,2) )
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

! !latex \item Some description of the internal notation follows :
! !latex \be G = g_{\theta\theta} g_{\varphi\varphi} - g_{\theta\varphi} g_{\theta\varphi}, \ee

!    dG(0)= gl(2,2,0,0)*gl(3,3,0,0)                                                                               -         gl(2,3,0,0)**2
!    dG(1)= gl(2,2,2,0)*gl(3,3,0,0)                           + gl(2,2,0,0)*gl(3,3,2,0)                           - two *   gl(2,3,0,0)*gl(2,3,2,0)
!    dG(2)= gl(2,2,2,2)*gl(3,3,0,0) + gl(2,2,2,0)*gl(3,3,2,0) + gl(2,2,2,0)*gl(3,3,2,0) + gl(2,2,0,0)*gl(3,3,2,2) - two * ( gl(2,3,2,0)*gl(2,3,2,0) + gl(2,3,0,0)*gl(2,3,2,2) )

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

! !latex \be A = g_{\varphi\varphi} 2 f_{\theta} f_{\theta\theta} - 2 g_{\theta\varphi} ( f_{\theta\theta} f_{\varphi} + f_{\theta} f_{\varphi\theta} ) 
! !latex + g_{\theta\theta} 2 f_{\varphi} f_{\varphi\theta}, \ee

!    dA(0) = gl(3,3,0,0)*two*df(2,0,0)*df(2,2,0)&
!          - two*gl(2,3,0,0)*( df(2,2,0)*df(3,0,0)+df(2,0,0)*df(3,2,0) ) &
!          + gl(2,2,0,0)*two*df(3,0,0)*df(3,2,0)
 
!    dA(1) = gl(3,3,2,0)*two*df(2,0,0)*df(2,2,0) + gl(3,3,0,0)*two*( df(2,2,0)*df(2,2,0) + df(2,0,0)*df(2,2,2) ) & 
!          - two* ( gl(2,3,2,0)*( df(2,2,0)*df(3,0,0)+df(2,0,0)*df(3,2,0) ) + gl(2,3,0,0)*( df(2,2,2)*df(3,0,0)+two*df(2,2,0)*df(3,2,0)+df(2,0,0)*df(3,2,2) ) ) &
!          + gl(2,2,2,0)*two*df(3,0,0)*df(3,2,0) + gl(2,2,0,0)*two*( df(3,2,0)*df(3,2,0) + df(3,0,0)*df(3,2,2) )

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

! !latex \be b_1 = g_{\varphi\varphi} f_{\theta}^2 - 2 g_{\theta\varphi} f_{\theta} f_{\varphi} + g_{\theta\theta} f_{\varphi}^2, \ee

!    db1(0) = gl(3,3,0,0)*df(2,0,0)**2 - two*gl(2,3,0,0)*df(2,0,0)*df(3,0,0) + gl(2,2,0,0)*df(3,0,0)**2
!    db1(1) = gl(3,3,2,0)*df(2,0,0)**2 - two*gl(2,3,2,0)*df(2,0,0)*df(3,0,0) + gl(2,2,2,0)*df(3,0,0)**2 + dA(0)

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

!    dabc(1,0) = gl(2,2,0,0)
!    dabc(1,1) = gl(2,2,2,0)
 
!    dabc(2,0) = -two * gl(2,3,0,0) * pt
!    dabc(2,1) = -two * gl(2,3,2,0) * pt
!    dabc(2,2) = -two * gl(2,3,0,0)
 
!    dabc(3,0) = gl(3,3,0,0) * pt**2 - db1(0) - two * delta_p * dG(0)
!    dabc(3,1) = gl(3,3,2,0) * pt**2 - db1(1) - two * delta_p * dG(1)
!    dabc(3,2) = gl(3,3,0,0) * pt*two
   
!    discrim=dabc(2,0)**2-four*dabc(1,0)*dabc(3,0)
   
!    if( discrim.lt.zero ) stop "ph00aa :          : WF=? ; discrim.lt.zero ;"
 
!    dP(0)  = ( -dabc(2,0) + plusminus *                                                                                   sqrt (discrim)                       ) / (two*dabc(1,0))
!    dP(1)  = ( -dabc(2,1) + plusminus * half * (two*dabc(2,0)*dabc(2,1)-four*(dabc(1,1)*dabc(3,0)+dabc(1,0)*dabc(3,1))) / sqrt (discrim) - dP(0)*two*dabc(1,1) ) / (two*dabc(1,0))
!    dP(2)  = ( -dabc(2,2) + plusminus * half * (two*dabc(2,0)*dabc(2,2)-four*(                    dabc(1,0)*dabc(3,2))) / sqrt (discrim)                       ) / (two*dabc(1,0))
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

!   !pphi=dP(0) ! returned through allglobal;

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

! !latex \be b_2 = g_{\varphi\varphi} p_{\theta}^2 - 2 g_{\theta\varphi} p_{\theta} p_{\varphi} + g_{\theta\theta} p_{\varphi}^2, \ee

!    db2(0) = gl(3,3,0,0)*   pt**2 - two*gl(2,3,0,0)* pt*dP(0)            + gl(2,2,0,0)*    dP(0)**2
!    db2(1) = gl(3,3,2,0)*   pt**2 - two*gl(2,3,2,0)* pt*dP(0)            + gl(2,2,2,0)*    dP(0)**2 & 
!                                  - two*gl(2,3,0,0)* pt*dP(1)            + gl(2,2,0,0)*two*dP(0)*dP(1)
!    db2(2) = gl(3,3,0,0)*two*pt   - two*gl(2,3,0,0)*(   dP(0)+pt*dP(2) ) + gl(2,2,0,0)*two*dP(0)*dP(2)

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

!    dtdot(0) = gl(3,3,0,0)*two*pt - two*gl(2,3,0,0)*dP(0)
!    dtdot(1) = gl(3,3,2,0)*two*pt - two*gl(2,3,2,0)*dP(0) - two*gl(2,3,0,0)*dP(1)
!    dtdot(2) = gl(3,3,0,0)*two                            - two*gl(2,3,0,0)*dP(2)
 
!    dpdot(0) = -two*gl(2,3,0,0)*pt + gl(2,2,0,0)*two*dP(0)
!    dpdot(1) = -two*gl(2,3,2,0)*pt + gl(2,2,2,0)*two*dP(0) + gl(2,2,0,0)*two*dP(1)
!    dpdot(2) = -two*gl(2,3,0,0)                            + gl(2,2,0,0)*two*dP(2)

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
 
!    ftftptpt=df(2,0,0)**2-pt**2 ; ftfpptpp= df(2,0,0)*df(3,0,0)-pt*dP(0) ; fpfppppp=df(3,0,0)**2-dP(0)**2

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

!    dsdot(0) = gl(3,3,2,0)* ftftptpt                 - two*gl(2,3,2,0)* ftfpptpp                                          + gl(2,2,2,0)*fpfppppp                              &
!             + dA(0) + (db2(0)-db1(0))*dG(1)/dG(0)
 
!    dsdot(1) = gl(3,3,2,2)* ftftptpt                 - two*gl(2,3,2,2)* ftfpptpp                                          + gl(2,2,2,2)*fpfppppp                              &
!             + gl(3,3,2,0)*(two*df(2,0,0)*df(2,2,0)) - two*gl(2,3,2,0)*(df(2,2,0)*df(3,0,0)+df(2,0,0)*df(3,2,0)-pt*dP(1)) + gl(2,2,2,0)*two*(df(3,0,0)*df(3,2,0)-dP(0)*dP(1)) & 
!             + dA(1) + (db2(1)-db1(1))*dG(1)/dG(0) + (db2(0)-db1(0))*(dG(2)-dG(1)**2/dG(0))/dG(0)
 
!    dsdot(2) = gl(3,3,2,0)*(-two*pt)                 - two*gl(2,3,2,0)*(-dP(0)-pt*dP(2))                                  + gl(2,2,2,0)*(-two*dP(0)*dP(2))                    &
!                     +  db2(2)        *dG(1)/dG(0) 

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

!    dptt(1:2) = (/ dsdot(0) , dtdot(0) /) / dpdot(0)
 
!    TM(1,1) = ( dsdot(2) - dptt(1)*dpdot(2) ) / dpdot(0) ! d radial / d radial ;
!    TM(1,2) = ( dsdot(1) - dptt(1)*dpdot(1) ) / dpdot(0) ! d radial / d angle  ;
!    TM(2,1) = ( dtdot(2) - dptt(2)*dpdot(2) ) / dpdot(0) ! d angle  / d radial ;
!    TM(2,2) = ( dtdot(1) - dptt(2)*dpdot(1) ) / dpdot(0) ! d angle  / d angle  ;
 
!    dptt(3)=TM(1,1)*ptt(3) + TM(1,2)*ptt(5) ! standard format;
!    dptt(4)=TM(1,1)*ptt(4) + TM(1,2)*ptt(6)
!    dptt(5)=TM(2,1)*ptt(3) + TM(2,2)*ptt(5)
!    dptt(6)=TM(2,1)*ptt(4) + TM(2,2)*ptt(6)

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
!   return
  
! end subroutine ph00aa

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

! !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! subroutine ph00aa_out(zeta,ptt)

!   use inputlist, only : 
!   use fileunits, only : ounit
!   use allglobal, only : pi2nfp

!   LOCALS

!   REAL, intent(inout) :: zeta
!   REAL, intent(in) :: ptt(6)
  
!   zeta = zeta + pi2nfp ! next intermediate output location;

!   return

! end subroutine ph00aa_out

! !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! REAL function ph00aa_end(zeta,ptt)

!   use allglobal, only : pi2nfp

!   LOCALS

!   REAL, intent(inout) :: zeta
!   REAL, intent(in) :: ptt(6)
  
!   ph00aa_end = zeta - pi2nfp ! integration termination;

!   return

! end FUNCTION ph00aa_end

! !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!


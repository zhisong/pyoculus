!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
! Fortran module for SPEC problems
! outputs ODEs for PJH
! written by Zhisong Qu (zhisong.qu@anu.edu.au)
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
! Usage:
!     To use the pjh module, one needs to initialize the fortran module first in Python, then initialize pjh in Python by
!       init_pjh(dp, inside_or_outside, plus_or_minus)
!     Parameters:
!       dp - delta p, the pressure jump
!       inside_or_outside - for the specified volume, we compute things on the inner interface or outer
!       plus_or_minus - whether to take the plus or minus sign in computing p_zeta
!
!     After initialization, one can call
!       rhs = get_pjhfield(phi , ptt)
!     or for rhs and tangent, 
!       rhs = get_pjhfield_tangent(phi , ptt)
!     Parameters:
!       phi - the zeta angle
!       ptt - the initial condition (theta, p_theta)
!     Returns:
!       rhs - array of size 2 (for get_pjhfield) or array of size 6 (for get_pjhfield_tangent)
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

MODULE pjh

  USE typedefns, only : REAL_KIND

  PRIVATE

  INTEGER                          :: innout       ! specifying field on which side (-1) inner side, (+1) outer side
  INTEGER                          :: ioi          ! specifying which interface
  INTEGER                          :: plusminus    ! +1 or -1 sign taken for calculating p_zeta
  REAL(KIND=REAL_KIND)             :: delta_p      ! the jump of pressure

  REAL(KIND=REAL_KIND),ALLOCATABLE :: gBtmne(:), gBzmne(:)  ! the contravariant components of magnetic field on the interface
  REAL(KIND=REAL_KIND),ALLOCATABLE :: gBtmno(:), gBzmno(:)  ! the contravariant components of magnetic field on the interface

  PUBLIC :: get_pjhfield, init_PJH, destroy_PJH, get_b2_interface

  CONTAINS

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  SUBROUTINE get_pjhfield( phi , ptt , dptt )
!f2py threadsafe

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

    USE constants
    USE variables
    USE coords, ONLY: get_metric_interface

    IMPLICIT NONE

    REAL(KIND=REAL_KIND),INTENT(IN) :: phi,ptt(2)
    REAL(KIND=REAL_KIND),INTENT(OUT) :: dptt(2)

    REAL :: delta_p

    INTEGER :: imn

    REAL(KIND=REAL_KIND) :: cosarg(mn), sinarg(mn),arg,pt,pp,theta
    REAL(KIND=REAL_KIND) :: dBB1(0:2),dBB2(0:3),a,b,c,discrim

    REAL(KIND=REAL_KIND) :: dR(0:3,0:3,0:3),dZ(0:3,0:3,0:3),df(0:3,0:3,0:3),gl(3,3,0:2),sg(0:2),dG(0:2),dA(0:1),db1(0:2),dabc(3,0:2),dP(0:2),db2(0:2)
    REAL(KIND=REAL_KIND) :: ftftptpt,ftfpptpp,fpfppppp
    REAL(KIND=REAL_KIND) :: dtdot(0:2),dpdot(0:2),dsdot(0:2),TM(2,2)
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

    pt = ptt(1) ! shorthand; "momentum" coordinate;

    theta = ptt(2) ! theta
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
    dR=zero ; dZ=zero; gl=zero; sg=zero
    call get_metric_interface(ioi, theta, phi, dR, dZ, gl, sg, 1)
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

! !latex \item Some description of the internal notation follows :
! !latex \be G = g_{\theta\theta} g_{\phi\phi} - g_{\theta\phi} g_{\theta\phi}, \ee

    dG(0)= gl(2,2,0)*gl(3,3,0)                                                                           -         gl(2,3,0)**2
    dG(1)= gl(2,2,1)*gl(3,3,0)                           + gl(2,2,0)*gl(3,3,1)                           - two *   gl(2,3,0)*gl(2,3,1)
!    dG(2)= gl(2,2,2,2)*gl(3,3,0,0) + gl(2,2,2,0)*gl(3,3,2,0) + gl(2,2,2,0)*gl(3,3,2,0) + gl(2,2,0,0)*gl(3,3,2,2) - two * ( gl(2,3,2,0)*gl(2,3,2,0) + gl(2,3,0,0)*gl(2,3,2,2) )

    call get_b2_interface(dBB1, ioi, theta, phi, gl, sg, 1)

    ! this is a, b, c for computing p_phi

    dabc(1,0) = gl(2,2,0)
    !dabc(1,1) = gl(2,2,1)
 
    dabc(2,0) = -two * gl(2,3,0) * pt
    !dabc(2,1) = -two * gl(2,3,1) * pt
    !dabc(2,2) = -two * gl(2,3,0)
 
    dabc(3,0) = gl(3,3,0) * pt**2 - (dBB1(0) + two * delta_p) * dG(0)
    !dabc(3,1) = gl(3,3,1) * pt**2 - (db1(0) + two * delta_p) * dG(1) - db1(1) * dG(0)
    !dabc(3,2) = gl(3,3,2) * pt*two
   
    discrim=dabc(2,0)**2-four*dabc(1,0)*dabc(3,0)
   
    if( discrim.lt.zero ) stop "ph00aa :          : WF=? ; discrim.lt.zero ;"
 
    dP(0)  = ( -dabc(2,0) + plusminus *                                                                                   sqrt (discrim)                       ) / (two*dabc(1,0))
!    dP(1)  = ( -dabc(2,1) + plusminus * half * (two*dabc(2,0)*dabc(2,1)-four*(dabc(1,1)*dabc(3,0)+dabc(1,0)*dabc(3,1))) / sqrt (discrim) - dP(0)*two*dabc(1,1) ) / (two*dabc(1,0))
!    dP(2)  = ( -dabc(2,2) + plusminus * half * (two*dabc(2,0)*dabc(2,2)-four*(                    dabc(1,0)*dabc(3,2))) / sqrt (discrim)                       ) / (two*dabc(1,0))
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

! !latex \be b_2 = g_{\phi\phi} p_{\theta}^2 - 2 g_{\theta\phi} p_{\theta} p_{\phi} + g_{\theta\theta} p_{\phi}^2, \ee

    db2(0)  = gl(3,3,0)*   pt**2 - two*gl(2,3,0)* pt*dP(0)            + gl(2,2,0)*    dP(0)**2
!    db2(1) = gl(3,3,2,0)*   pt**2 - two*gl(2,3,2,0)* pt*dP(0)            + gl(2,2,2,0)*    dP(0)**2 & 
!                                  - two*gl(2,3,0,0)* pt*dP(1)            + gl(2,2,0,0)*two*dP(0)*dP(1)
!    db2(2) = gl(3,3,0,0)*two*pt   - two*gl(2,3,0,0)*(   dP(0)+pt*dP(2) ) + gl(2,2,0,0)*two*dP(0)*dP(2)

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

    dtdot(0) = gl(3,3,0)*two*pt - two*gl(2,3,0)*dP(0)
!    dtdot(1) = gl(3,3,2,0)*two*pt - two*gl(2,3,2,0)*dP(0) - two*gl(2,3,0,0)*dP(1)
!    dtdot(2) = gl(3,3,0,0)*two                            - two*gl(2,3,0,0)*dP(2)
 
    dpdot(0) = -two*gl(2,3,0)*pt + gl(2,2,0)*two*dP(0)
!    dpdot(1) = -two*gl(2,3,2,0)*pt + gl(2,2,2,0)*two*dP(0) + gl(2,2,0,0)*two*dP(1)
!    dpdot(2) = -two*gl(2,3,0,0)                            + gl(2,2,0,0)*two*dP(2)

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
 
    ftftptpt= -pt**2 ; ftfpptpp= -pt*dP(0) ; fpfppppp=-dP(0)**2

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

    dsdot(0) = gl(3,3,1)* ftftptpt                 - two*gl(2,3,1)* ftfpptpp                                          + gl(2,2,1)*fpfppppp                              &
             + dBB1(1)*G(0) + db2(0)*dG(1)/dG(0)
 
!    dsdot(1) = gl(3,3,2,2)* ftftptpt                 - two*gl(2,3,2,2)* ftfpptpp                                          + gl(2,2,2,2)*fpfppppp                              &
!             + gl(3,3,2,0)*(two*df(2,0,0)*df(2,2,0)) - two*gl(2,3,2,0)*(df(2,2,0)*df(3,0,0)+df(2,0,0)*df(3,2,0)-pt*dP(1)) + gl(2,2,2,0)*two*(df(3,0,0)*df(3,2,0)-dP(0)*dP(1)) & 
!             + dA(1) + (db2(1)-db1(1))*dG(1)/dG(0) + (db2(0)-db1(0))*(dG(2)-dG(1)**2/dG(0))/dG(0)
 
!    dsdot(2) = gl(3,3,2,0)*(-two*pt)                 - two*gl(2,3,2,0)*(-dP(0)-pt*dP(2))                                  + gl(2,2,2,0)*(-two*dP(0)*dP(2))                    &
!                     +  db2(2)        *dG(1)/dG(0) 

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

    dptt(1:2) = (/ dsdot(0) , dtdot(0) /) / dpdot(0)
 
!    TM(1,1) = ( dsdot(2) - dptt(1)*dpdot(2) ) / dpdot(0) ! d radial / d radial ;
!    TM(1,2) = ( dsdot(1) - dptt(1)*dpdot(1) ) / dpdot(0) ! d radial / d angle  ;
!    TM(2,1) = ( dtdot(2) - dptt(2)*dpdot(2) ) / dpdot(0) ! d angle  / d radial ;
!    TM(2,2) = ( dtdot(1) - dptt(2)*dpdot(1) ) / dpdot(0) ! d angle  / d angle  ;
 
!    dptt(3)=TM(1,1)*ptt(3) + TM(1,2)*ptt(5) ! standard format;
!    dptt(4)=TM(1,1)*ptt(4) + TM(1,2)*ptt(6)
!    dptt(5)=TM(2,1)*ptt(3) + TM(2,2)*ptt(5)
!    dptt(6)=TM(2,1)*ptt(4) + TM(2,2)*ptt(6)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  END SUBROUTINE get_pjhfield

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  SUBROUTINE get_b2_interface(db2, ioi, theta, zeta, guvij, sg, ideriv)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

    USE typedefns, ONLY : REAL_KIND
    USE constants
    USE variables

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

    IMPLICIT NONE
    REAL(KIND=REAL_KIND),    INTENT(IN) :: theta, zeta ! the three coordinates (s,theta,zeta)
    INTEGER             ,    INTENT(IN) :: ioi ! on which interface
    INTEGER             ,    INTENT(IN) :: ideriv ! the level of derivatives. 0 for none, 1 for first, 2 for second

    REAL(KIND=REAL_KIND),    INTENT(IN):: guvij(3,3,0:2), sg(0:2)

    REAL(KIND=REAL_KIND),    INTENT(OUT):: db2(0:2)

    REAL(KIND=REAL_KIND) :: lss, rsign
    REAL(KIND=REAL_KIND) :: cosarg(mn), sinarg(mn), gB(2:3,0:2)

    INTEGER :: ii, jj

    db2 = zero

    cosarg = COS(im*theta-in1*zeta)
    sinarg = SIN(im*theta-in1*zeta)

    gB(2,0) = SUM(gBtmne * cosarg)
    gB(3,0) = SUM(gBzmne * cosarg)

    DO ii = 2, 3
      DO jj = 2, 3
        db2(0) = db2(0) + guvij(ii,jj,0) * gB(ii,0) * gB(jj,0)
      ENDDO
    ENDDO
    db2(0) = db2(0) / sg(0)**2

    IF (ideriv .ge. 1) THEN
      gB(2,1) = SUM(gBtmne * (-im) * sinarg)
      gB(3,1) = SUM(gBzmne * (-im) * sinarg)

      DO ii = 2, 3
        DO jj = 2, 3
          db2(1) = db2(1) + guvij(ii,jj,1) * gB(ii,0) * gB(jj,0) + guvij(ii,jj,0) * (gB(ii,1) * gB(jj,0) + gB(ii,0) * gB(jj,1))
        ENDDO
      ENDDO
      db2(1) = db2(1) / sg(0)**2 - 2 * db2(0) * sg(1) / sg(0)
    ENDIF

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  END SUBROUTINE get_b2_interface

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  SUBROUTINE init_pjh(dp, inside_or_outside, plus_or_minus)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
    USE constants
    USE variables
    USE basefn

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: inside_or_outside ! take the inner surface or outer surface of a volume
    INTEGER, INTENT(IN) :: plus_or_minus     ! choose plus sign or minus sign for computing p_zeta
    REAL(KIND=REAL_KIND), INTENT(IN) :: dp   ! delta_p, the pressure jump

    INTEGER :: ii, jj, ll, mi, ni

    REAL(KIND=REAL_KIND) :: sbar, lss
    REAL(KIND=REAL_KIND) :: cheby(0:Lrad,0:1), zernike(0:Lrad,0:Mpol,0:1)
    
    REAL(KIND=REAL_KIND) :: TT(0:Lrad,0:1)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

    plusminus = plus_or_minus
    delta_p = dp
    ioi = ivol + inside_or_outside - 1       ! the index of the interface

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

  SUBROUTINE destroy_pjh()

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

    USE variables

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

END MODULE pjh

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
!latex B^2 = \frac{g_{\phi\phi} f_\theta f_\theta - 2 g_{\theta\phi}f_\theta f_\phi + g_{\theta\theta} f_\phi f_\phi}{g_{\theta\theta}g_{\phi\phi}-g_{\theta\phi}g_{\theta\phi}}
!latex \ee
!latex where $f$ is a surface potential and $g_{\theta\theta}$, $g_{\theta\phi}$ and $g_{\phi\phi}$ are metric elements local to the interface.
!latex \item Assuming that the field $B_1$ is known on the `inside' of the interface, ie. $B_{1\theta}=f_\theta$, $B_{1\phi}=f_\phi$ and $f$ is known, 
!latex it is required to determine the tangential field, $p_\theta = B_\theta$ and $p_\phi = B_\phi$, on the `outside' of the interface.
!latex \item The o.d.e.'s are given by Hamilton's equations \be
!latex \dot \theta   =  \frac{\partial H}{\partial p_\theta}\Big|_{\theta,\phi,p_\phi}, \;\;
!latex \dot p_\theta = -\frac{\partial H}{\partial \theta}\Big|_{p_\theta,\phi,p_\phi}, \;\;
!latex \dot \phi     =  \frac{\partial H}{\partial p_\phi}\Big|_{\theta,p_\theta,\phi}, \;\;
!latex \dot p_\phi   = -\frac{\partial H}{\partial \phi}\Big|_{\theta,p_\theta,p_\phi}, \ee
!latex where the `dot' denotes derivative with respect to `time'.
!latex \item This is reduced to a $1\frac{1}{2}$ dimensional system by using $\phi$ as the time-like integration parameter, and replacing the equation for $\dot p_\phi$ with 
!latex \be p_\phi= P(\theta,p_\theta,\phi; \delta p) = \frac{-b\pm\sqrt{b^2-4ac}}{2a} \label{eq:pphi} \ee
!latex where \mbox{$a=g_{\theta\theta}$}, \mbox{$b=-2 g_{\theta\phi}p_\theta$} 
!latex and \mbox{$c=g_{\phi\phi} p_{\theta}^2 - b_1 - 2 \, \delta p \, G$} (see below for definition of $b_1$ and $G$).
!latex \item The o.d.e.'s that then need to be followed are (see below for definition of $A$ and $b_2$) \be
!latex \frac{d   \theta}{d\phi}&=& \frac{g_{\phi\phi} p_{\theta} - g_{\theta\phi} p_{\phi}}{-g_{\theta\phi}p_{\theta}+g_{\theta\theta}p_{\phi}},\\
!latex \frac{d p_\theta}{d\phi}&=& \frac{g_{\phi\phi,\theta} (f_{\theta}^2-p_{\theta}^2)-2g_{\theta\phi,\theta}(f_{\theta}f_{\phi}-p_{\theta}p_{\phi})
!latex +g_{\theta\theta,\theta}(f_{\phi}^2-p_{\phi}^2)+A+(b_2-b_1)G_{,\theta} / G }
!latex {-2g_{\theta\phi}p_\theta+g_{\theta\theta}2p_{\phi}}. \ee

!latex \item Note that $d\theta / d \phi = B^\theta / B^\phi $; there is a fundamental relation between the pressure-jump Hamiltonian and the field-line Hamiltonian. 
!latex (Furthermore, in many cases the surface will be given in straight field line coordinates, so $d \theta / d\phi = const.$.)

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

! subroutine ph00aa( phi , ptt , dptt ) ! fixed format;

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
!   use constants
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
! !latex \be G = g_{\theta\theta} g_{\phi\phi} - g_{\theta\phi} g_{\theta\phi}, \ee

!    dG(0)= gl(2,2,0,0)*gl(3,3,0,0)                                                                               -         gl(2,3,0,0)**2
!    dG(1)= gl(2,2,2,0)*gl(3,3,0,0)                           + gl(2,2,0,0)*gl(3,3,2,0)                           - two *   gl(2,3,0,0)*gl(2,3,2,0)
!    dG(2)= gl(2,2,2,2)*gl(3,3,0,0) + gl(2,2,2,0)*gl(3,3,2,0) + gl(2,2,2,0)*gl(3,3,2,0) + gl(2,2,0,0)*gl(3,3,2,2) - two * ( gl(2,3,2,0)*gl(2,3,2,0) + gl(2,3,0,0)*gl(2,3,2,2) )

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

! !latex \be A = g_{\phi\phi} 2 f_{\theta} f_{\theta\theta} - 2 g_{\theta\phi} ( f_{\theta\theta} f_{\phi} + f_{\theta} f_{\phi\theta} ) 
! !latex + g_{\theta\theta} 2 f_{\phi} f_{\phi\theta}, \ee

!    dA(0) = gl(3,3,0,0)*two*df(2,0,0)*df(2,2,0)&
!          - two*gl(2,3,0,0)*( df(2,2,0)*df(3,0,0)+df(2,0,0)*df(3,2,0) ) &
!          + gl(2,2,0,0)*two*df(3,0,0)*df(3,2,0)
 
!    dA(1) = gl(3,3,2,0)*two*df(2,0,0)*df(2,2,0) + gl(3,3,0,0)*two*( df(2,2,0)*df(2,2,0) + df(2,0,0)*df(2,2,2) ) & 
!          - two* ( gl(2,3,2,0)*( df(2,2,0)*df(3,0,0)+df(2,0,0)*df(3,2,0) ) + gl(2,3,0,0)*( df(2,2,2)*df(3,0,0)+two*df(2,2,0)*df(3,2,0)+df(2,0,0)*df(3,2,2) ) ) &
!          + gl(2,2,2,0)*two*df(3,0,0)*df(3,2,0) + gl(2,2,0,0)*two*( df(3,2,0)*df(3,2,0) + df(3,0,0)*df(3,2,2) )

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

! !latex \be b_1 = g_{\phi\phi} f_{\theta}^2 - 2 g_{\theta\phi} f_{\theta} f_{\phi} + g_{\theta\theta} f_{\phi}^2, \ee

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

! !latex \be b_2 = g_{\phi\phi} p_{\theta}^2 - 2 g_{\theta\phi} p_{\theta} p_{\phi} + g_{\theta\theta} p_{\phi}^2, \ee

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


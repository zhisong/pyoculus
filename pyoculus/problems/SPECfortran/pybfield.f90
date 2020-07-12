!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
! Fortran module for SPEC problems
! outputs contravariant Bfield given (s, theta zeta) coordinates
! adapted bfield.f90 in SPEC main code
! written by Zhisong Qu (zhisong.qu@anu.edu.au)
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

MODULE bfield

  CONTAINS

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  SUBROUTINE get_bfield( zeta, st, Bst )
!f2py threadsafe
    
  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
    
    USE typedefns, ONLY : REAL_KIND

    USE constants, ONLY : zero, one, half, two
    
    USE variables

    USE basefn
    
    IMPLICIT NONE
    
    REAL(KIND=REAL_KIND), INTENT(IN)   :: zeta,  st(2)
    REAL(KIND=REAL_KIND), INTENT(OUT)  ::       Bst(2)
    
    INTEGER              :: lvol, ii, ll, mi, ni
    REAL(KIND=REAL_KIND) :: teta, lss, sbar, arg, carg, sarg, dBu(1:3)
    REAL(KIND=REAL_KIND) :: cheby(0:Lrad,0:1), zernike(0:Lrad,0:Mpol,0:1)
    
    REAL               :: TT(0:Lrad,0:1) ! this is almost identical to cheby; 17 Dec 15;

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

    lvol = ivol ;  ! short hand

    Bst(1:2) = (/ zero , zero /) ! set default intent out; this should cause a compilation error IF Node.ne.2;
  
  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
    
    lss = st(1) ; teta = st(2) 
      
    IF( Lcoordinatesingularity ) sbar = MAX( ( lss + one ) * half, zero )

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
    
    IF (Lcoordinatesingularity) THEN
      CALL get_zernike(sbar, Lrad, Mpol, zernike(:,:,0:1))
    ELSE
      CALL get_cheby(lss, Lrad, cheby(0:Lrad,0:1))
    end IF

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

    dBu(1:3) = zero ! initialize summation;

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
    
    DO ii = 1, mn ; mi = im(ii) ; ni = in1(ii) ; arg = mi * teta - ni * zeta ; carg = COS(arg) ; sarg = SIN(arg) ! shorthand; 20 Apr 13;
    
  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
    
      IF( Lcoordinatesingularity ) THEN ! regularization factor depends on mi; 17 Dec 15;

        DO ll = 0, Lrad ; TT(ll,0:1) = (/        zernike(ll,mi,0),        zernike(ll,mi,1)*half                      /)
        ENDDO

      ELSE

        DO ll = 0, Lrad ; TT(ll,0:1) = (/             cheby(ll,0),             cheby(ll,1)                           /)
        ENDDO

      ENDIF ! end of IF( Lcoordinatesingularity ) ; 16 Jan 15;

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

    ! DO ll = 0, Lrad ! loop over Chebyshev summation; 20 Feb 13;
    !   ;dBu(1) = dBu(1) + ( - mi * Aze(ll+1,ii) - ni * Ate(ll+1,ii) ) * TT(ll,0) * sarg
    !   ;dBu(2) = dBu(2) + (                                  -      Aze(ll+1,ii) ) * TT(ll,1) * carg
    !   ;dBu(3) = dBu(3) + (        Ate(ll+1,ii)                                  ) * TT(ll,1) * carg
    !   IF( NOTstellsym ) THEN ! include non-symmetric harmonics; 28 Jan 13;
    !   dBu(1) = dBu(1) + ( + mi * Azo(ll+1,ii) + ni * Ato(ll+1,ii) ) * TT(ll,0) * carg
    !   dBu(2) = dBu(2) + (                                  -      Azo(ll+1,ii) ) * TT(ll,1) * sarg
    !   dBu(3) = dBu(3) + (        Ato(ll+1,ii)                                  ) * TT(ll,1) * sarg
    !   ENDIF
    ! ENDDO ! end of DO ll; 10 Dec 15;

      dBu(1) = dBu(1) + SUM(( - mi * Aze(1:Lrad+1,ii)- ni * Ate(1:Lrad+1,ii) ) * TT(0:Lrad, 0)) * sarg
      dBu(2) = dBu(2) + SUM((                        -      Aze(1:Lrad+1,ii) ) * TT(0:Lrad, 1)) * carg
      dBu(3) = dBu(3) + SUM((        Ate(1:Lrad+1,ii)                        ) * TT(0:Lrad, 1)) * carg

      IF( NOTstellsym ) THEN
        dBu(1) = dBu(1) + SUM(( + mi * Azo(1:Lrad+1,ii) + ni * Ato(1:Lrad+1,ii) ) * TT(0:Lrad,0)) * carg
        dBu(2) = dBu(2) + SUM((                         -      Azo(1:Lrad+1,ii) ) * TT(0:Lrad,1)) * sarg
        dBu(3) = dBu(3) + SUM((        Ato(1:Lrad+1,ii)                         ) * TT(0:Lrad,1)) * sarg
      ENDIF
  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
    
    ENDDO ! end of DO ii = 1, mn;
    
  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

    Bst(1:2) = dBu(1:2) / dBu(3) ! normalize field line equations to toroidal field; 20 Apr 13;
    
  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  END SUBROUTINE get_bfield

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

END MODULE bfield

! subroutine bfield_tangent( zeta, st, Bst ) 
!f2py threadsafe 
! !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!   use constants, only : zero, one, half, two
  
!   use numerical, only : vsmall, small
  
!   use fileunits, only : ounit
  
!   use inputlist, only : Wmacros, Wbfield, Lrad, Mpol
  
!   use cputiming, only : Tbfield
  
!   use allglobal, only : myid, ncpu, cpus, mn, im, in, halfmm, regumm, &
!                         ivol, gBzeta, Ate, Aze, Ato, Azo, &
!                         NOTstellsym, &
!                         Lcoordinatesingularity, Mvol, &
!                         Node ! 17 Dec 15;
  
! !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!   LOCALS
  
!   REAL, intent(in)   :: zeta,  st(1:6)
!   REAL, intent(out)  ::       Bst(1:6)
  
!   INTEGER            :: lvol, ii, ll, mi, ni, ideriv
!   REAL               :: teta, lss, sbar, sbarhm(0:1), arg, carg, sarg, dBu(1:3,0:2)
!   REAL               :: cheby(0:Lrad(ivol),0:2), zernike(0:Lrad(1),0:Mpol,0:2)

!   REAL               :: M(2,2), deltax(2,2)
  
!   REAL               :: TT(0:Lrad(ivol),0:2) ! this is almost identical to cheby; 17 Dec 15;

!   BEGIN(bfield)

! !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! #IFdef DEBUG
!   FATAL( bfield, ivol.lt.1 .or. ivol.gt.Mvol, invalid ivol )
! #ENDIF

! !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!   lvol = ivol ; ideriv = 0 ! the argument list of bfield is fixed by NAG requirements, but volume index is required below;

!   Bst = zero ! set default intent out; this should cause a compilation error IF Node.ne.2;
 
! !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!   lss = st(1) ; teta = st(2) ;
  
!   ! the perturbation
!   deltax(1:2,1) = st(3:4);
!   deltax(1:2,2) = st(5:6);
  
!   IF( abs(lss).gt.one ) goto 9999 ! out of DOmain;
    
!   IF( Lcoordinatesingularity ) sbar = max( ( lss + one ) * half, small )

! !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!   IF (Lcoordinatesingularity) THEN
!     call get_zernike_d2(sbar, Lrad(lvol), Mpol, zernike(:,:,0:2))
!   ELSE
!     call get_cheby_d2(lss, Lrad(lvol), cheby(0:Lrad(lvol),0:2))
!   end IF

! !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!   dBu(1:3,0:2) = zero ! initialize summation;
  
! !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!   DO ii = 1, mn ; mi = im(ii) ; ni = in(ii) ; arg = mi * teta - ni * zeta ; carg = cos(arg) ; sarg = sin(arg) ! shorthand; 20 Apr 13;
   
! !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
!    IF( Lcoordinatesingularity ) THEN ! regularization factor depends on mi; 17 Dec 15;

!     FATAL( bfield, abs(sbar).lt.vsmall, need to avoid divide-by-zero )

!     DO ll = 0, Lrad(lvol) ; TT(ll,0:2) = (/        zernike(ll,mi,0),        zernike(ll,mi,1)*half , zernike(ll,mi,2)*half*half /)
!     ENDDO

!    ELSE

!     DO ll = 0, Lrad(lvol) ; TT(ll,0:2) = (/             cheby(ll,0),             cheby(ll,1)       ,cheby(ll,2)               /)
!     ENDDO

!    ENDIF ! end of IF( Lcoordinatesingularity ) ; 16 Jan 15;

! !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
!    DO ll = 0, Lrad(lvol) ! loop over Chebyshev summation; 20 Feb 13;
!     ! no derivative
!     ;dBu(1,0) = dBu(1,0) + ( - mi * Aze(ll,ii) - ni * Ate(ll,ii) ) * TT(ll,0) * sarg
!     ;dBu(2,0) = dBu(2,0) + (                                  -      Aze(ll,ii) ) * TT(ll,1) * carg
!     ;dBu(3,0) = dBu(3,0) + (        Ate(ll,ii)                                  ) * TT(ll,1) * carg
!     ! ds
!     ;dBu(1,1) = dBu(1,1) + ( - mi * Aze(ll,ii) - ni * Ate(ll,ii) ) * TT(ll,1) * sarg
!     ;dBu(2,1) = dBu(2,1) + (                                  -      Aze(ll,ii) ) * TT(ll,2) * carg
!     ;dBu(3,1) = dBu(3,1) + (        Ate(ll,ii)                                  ) * TT(ll,2) * carg
!     ! dtheta
!     ;dBu(1,2) = dBu(1,2) + mi * ( - mi * Aze(ll,ii) - ni * Ate(ll,ii) ) * TT(ll,0) * carg
!     ;dBu(2,2) = dBu(2,2) - mi * (                                  -      Aze(ll,ii) ) * TT(ll,1) * sarg
!     ;dBu(3,2) = dBu(3,2) - mi * (        Ate(ll,ii)                                  ) * TT(ll,1) * sarg
!     IF( NOTstellsym ) THEN ! include non-symmetric harmonics; 28 Jan 13;
!      ! no derivative
!      dBu(1,0) = dBu(1,0) + ( + mi * Azo(ll,ii) + ni * Ato(ll,ii) ) * TT(ll,0) * carg
!      dBu(2,0) = dBu(2,0) + (                                  -      Azo(ll,ii) ) * TT(ll,1) * sarg
!      dBu(3,0) = dBu(3,0) + (        Ato(ll,ii)                                  ) * TT(ll,1) * sarg
!      ! ds
!      dBu(1,1) = dBu(1,1) + ( + mi * Azo(ll,ii) + ni * Ato(ll,ii) ) * TT(ll,1) * carg
!      dBu(2,1) = dBu(2,1) + (                                  -      Azo(ll,ii) ) * TT(ll,2) * sarg
!      dBu(3,1) = dBu(3,1) + (        Ato(ll,ii)                                  ) * TT(ll,2) * sarg
!      ! dtheta
!      dBu(1,2) = dBu(1,2) - mi * ( + mi * Azo(ll,ii) + ni * Ato(ll,ii) ) * TT(ll,0) * sarg
!      dBu(2,2) = dBu(2,2) + mi * (                                  -      Azo(ll,ii) ) * TT(ll,1) * carg
!      dBu(3,2) = dBu(3,2) + mi * (        Ato(ll,ii)                                  ) * TT(ll,1) * carg
!     ENDIF
!    ENDDO ! end of DO ll; 10 Dec 15;
   
! !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
!   ENDDO ! end of DO ii = 1, mn;
  
! !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!   gBzeta = dBu(3,0) ! gBzeta is returned through global; 20 Apr 13;

! !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!   IF( abs(gBzeta).lt.vsmall ) THEN

!    cput = GETTIME

!    write(ounit,'("bfield : ",f10.2," : lvol=",i3," ; zeta="es23.15" ; (s,t)=("es23.15" ,"es23.15" ) ; B^z="es23.15" ;")') &
!                              cput-cpus, lvol,        zeta,             st(1:2),                       gBzeta

!    FATAL( bfield, abs(gBzeta).lt.vsmall, field is not toroidal )

!   ENDIF
  
! !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!   Bst(1:2) = dBu(1:2,0) / gBzeta ! normalize field line equations to toroidal field; 20 Apr 13;


!   ! assemble the tangent matrix
!   M(1,1) = (dBu(1,1) * dBu(3,0) - dBu(3,1) * dBu(1,0)) / gBzeta**2   ! d(Bs/Bz)/ds
!   M(1,2) = (dBu(1,2) * dBu(3,0) - dBu(3,2) * dBu(1,0)) / gBzeta**2   ! d(Bs/Bz)/dt
!   M(2,1) = (dBu(2,1) * dBu(3,0) - dBu(3,1) * dBu(2,0)) / gBzeta**2   ! d(Bt/Bz)/ds
!   M(2,2) = (dBu(2,2) * dBu(3,0) - dBu(3,2) * dBu(2,0)) / gBzeta**2   ! d(Bt/Bz)/dt

!   deltax = MATMUL(M, deltax)

!   Bst(3:4) = deltax(1:2,1)
!   Bst(5:6) = deltax(1:2,2)

! !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!   RETURN(bfield)
  
! !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! end subroutine bfield_tangent

! !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

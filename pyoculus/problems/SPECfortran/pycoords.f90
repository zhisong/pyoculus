!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
! Fortran module for SPEC problems
! outputs the Cartesian coordinates given (s, theta zeta) coordinates
! output metrics
! written by Zhisong Qu (zhisong.qu@anu.edu.au)
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

MODULE coords

  CONTAINS

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  SUBROUTINE get_xyz(stz , RpZ )
!f2py threadsafe  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
    use typedefns, only : REAL_KIND
    use constants, only : zero, one, half
    use variables, only : Igeometry, Ntor, mn, im, in1, &
                          iRbc, iZbs, iRbs, iZbc, &
                          Lcoordinatesingularity, &
                          NOTstellsym, &
                          Mvol, ivol
    
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
    IMPLICIT NONE  
  
    
    REAL(KIND=REAL_KIND),    INTENT(IN)  :: stz(1:3)
    REAL(KIND=REAL_KIND),    INTENT(OUT) :: RpZ(1:3)
  
    INTEGER              :: ii, mi, ni, lvol
    REAL(KIND=REAL_KIND) :: Remn, Zomn, Romn, Zemn, RR, phi, ZZ, arg, carg, sarg, lss, alss, blss, sbar, sbarhim, fj
  
    lvol = ivol

    RpZ(1:3) = zero ; RR = zero ; phi = stz(3) ; ZZ = zero ! initialize intent(out), summations ; 17 Dec 15;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
    lss = stz(1) ; alss = half - lss * half ; blss = half + lss * half ! shorthand; 17 Dec 15;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

    DO ii = 1, mn ; mi = im(ii) ; ni = in1(ii) ; arg = mi * stz(2) - ni * phi ; carg = COS(arg) ; sarg = SIN(arg) ! Fourier sum; 17 Dec 15;
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
      Remn = zero ; Zomn = zero
      Romn = zero ; Zemn = zero
   
      IF( Lcoordinatesingularity ) THEN
    
        sbar = ( lss + one ) * half
    
        IF( mi.EQ.0 ) THEN
          IF( Igeometry.EQ.2 ) THEN ; fj = sbar
          ELSE                      ; fj = sbar**2
          ENDIF
        ELSE
          IF( Igeometry.EQ.2 ) THEN ; fj = sbar**im(ii)
          ELSE                      ; fj = sbar**(im(ii)+1)
          ENDIF                        ; 
        ENDIF
    
        Remn = iRbc(ii,1) + ( iRbc(ii,2) - iRbc(ii,1) ) * fj
        IF( NOTstellsym ) THEN
          Romn = iRbs(ii,1) + ( iRbs(ii,2) - iRbs(ii,1) ) * fj
        ENDIF
        IF( Igeometry.eq.3 ) THEN ! recall that for cylindrical geometry there is no need for Z; 20 Apr 13;
          Zomn = iZbs(ii,1) + ( iZbs(ii,2) - iZbs(ii,1) ) * fj
          IF( NOTstellsym ) THEN
            Zemn = iZbc(ii,1) + ( iZbc(ii,2) - iZbc(ii,1) ) * fj
          ENDIF
        ENDIF

      ELSE ! matches IF( Lcoordinatesingularity) ; 20 Apr 13;
    
        Remn =   alss * iRbc(ii,lvol) + blss * iRbc(ii,lvol+1)
        IF( NOTstellsym ) THEN
          Romn =   alss * iRbs(ii,lvol) + blss * iRbs(ii,lvol+1)
        ELSE
          Romn = zero
        ENDIF ! end of IF( NOTstellsym ) ; 22 Apr 13;
        IF( Igeometry.eq.3 ) THEN
          Zomn =   alss * iZbs(ii,lvol) + blss * iZbs(ii,lvol+1)
          IF( NOTstellsym ) THEN
            Zemn =   alss * iZbc(ii,lvol) + blss * iZbc(ii,lvol+1)
          ELSE
            Zemn = zero
          ENDIF ! end of IF( NOTstellsym ) ; 22 Apr 13; 
        ELSE    
          Zomn = zero
          Zemn = zero
        ENDIF ! end of IF( Igeometry.eq.3 ) ; 22 Apr 13;
    
      ENDIF ! end of IF( Lcoordinatesingularity ) ; 20 Feb 13;
   
      RR = RR + Remn * carg + Romn * sarg
      IF( Igeometry.eq.3 ) THEN
        ZZ = ZZ + Zemn * carg + Zomn * sarg
      ENDIF
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
    ENDDO ! end of do ii = 1, mn ; Fourier summation loop;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

    RpZ(1:3) = (/ RR, phi, ZZ /) ! return coordinate functions;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  END SUBROUTINE get_xyz

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

END MODULE coords
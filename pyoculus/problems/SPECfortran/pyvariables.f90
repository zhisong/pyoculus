!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
! Fortran MODULE for SPEC problems
! a MODULE that keeps all the global variables for other interfaced fortran SUBROUTINE
! written by Zhisong Qu (zhisong.qu@anu.edu.au)
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

MODULE typedefns
  INTEGER, PARAMETER :: REAL_KIND = 8
END MODULE typedefns

MODULE constants

  USE typedefns, only : REAL_KIND

  IMPLICIT NONE

  REAL(KIND=REAL_KIND), PARAMETER :: zero       =    0.0
  REAL(KIND=REAL_KIND), PARAMETER :: one        =    1.0
  REAL(KIND=REAL_KIND), PARAMETER :: two        =    2.0
  REAL(KIND=REAL_KIND), PARAMETER :: three      =    3.0
  REAL(KIND=REAL_KIND), PARAMETER :: four       =    4.0
  REAL(KIND=REAL_KIND), PARAMETER :: five       =    5.0
  REAL(KIND=REAL_KIND), PARAMETER :: six        =    6.0
  REAL(KIND=REAL_KIND), PARAMETER :: seven      =    7.0
  REAL(KIND=REAL_KIND), PARAMETER :: eight      =    8.0
  REAL(KIND=REAL_KIND), PARAMETER :: nine       =    9.0
  REAL(KIND=REAL_KIND), PARAMETER :: ten        =   10.0

  REAL(KIND=REAL_KIND), PARAMETER :: eleven     =   11.0
  REAL(KIND=REAL_KIND), PARAMETER :: twelve     =   12.0

  REAL(KIND=REAL_KIND), PARAMETER :: hundred    =  100.0
  REAL(KIND=REAL_KIND), PARAMETER :: thousand   = 1000.0

  REAL(KIND=REAL_KIND), PARAMETER :: half       =   one / two
  REAL(KIND=REAL_KIND), PARAMETER :: third      =   one / three 
  REAL(KIND=REAL_KIND), PARAMETER :: quart      =   one / four
  REAL(KIND=REAL_KIND), PARAMETER :: fifth      =   one / five
  REAL(KIND=REAL_KIND), PARAMETER :: sixth      =   one / six

  REAL(KIND=REAL_KIND), PARAMETER :: pi2        =   6.28318530717958623
  REAL(KIND=REAL_KIND), PARAMETER :: pi         =   pi2 / two
  REAL(KIND=REAL_KIND), PARAMETER :: mu0        =   2.0E-07 * pi2
  REAL(KIND=REAL_KIND), PARAMETER :: goldenmean =   1.618033988749895 ! golden mean = ( one + sqrt(five) ) / two ;

  REAL(KIND=REAL_KIND), PARAMETER :: version    =   2.20

END MODULE constants

MODULE variables

  USE typedefns, only : REAL_KIND

  ! some general variables that exist in SPEC input and output
  INTEGER                         :: Igeometry
  INTEGER                         :: Mvol
  INTEGER                         :: Ntor
  INTEGER                         :: Mpol
  INTEGER                         :: Nfp
  INTEGER                         :: mn
  LOGICAL                         :: NOTstellsym
  INTEGER, allocatable            :: im(:), in1(:)            ! Fourier modes;
  REAL(KIND=REAL_KIND),allocatable:: iRbc(:,:) , iZbs(:,:)   ! interface surface geometry;     stellarator symmetric;
  REAL(KIND=REAL_KIND),allocatable:: iRbs(:,:) , iZbc(:,:)   ! interface surface geometry; non-stellarator symmetric;

  ! variables for the volume we are interested in
  INTEGER                         :: ivol       ! which volume we are stuyding
  INTEGER                         :: Lrad       ! what is the radial resolution for that volume
  LOGICAL                         :: Lcoordinatesingularity
  REAL(KIND=REAL_KIND),allocatable:: Ate(:,:), Aze(:,:), Ato(:,:), Azo(:,:)

CONTAINS

END MODULE
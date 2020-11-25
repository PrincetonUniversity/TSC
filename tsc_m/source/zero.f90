      subroutine zero
!......5.80 zero
!***********************************************************
!                                                          *
!.....zero out arrays                                      *
!                                                          *
!***********************************************************
      USE CLINAM
      USE SCRATCH
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!l    zero out equlibrium quantities:
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER i,j
!============
      do 100 i=1,nxp
      xary(i)    = 0._R8
      do 100 j=1,nzp
      psi(i,j)   = 0._R8
      vecx(i,j)  = 0._R8
      g(i,j)     = 0._R8
      q(i,j)     = 0._R8
      psio(i,j)  = 0._R8
  100 continue
      do 101 j=1,nzp
      zary(j)    = 0._R8
  101 continue
!
!l    zero out quantities used in time-dependent part:
      do 200 i=1,nxp
      do 200 j=1,nzp
      abig(i,j)  = 0._R8
      vecz(i,j)  = 0._R8
      w(i,j)     = 0._R8
      u(i,j)     = 0._R8
      wo(i,j)    = 0._R8
      go(i,j)    = 0._R8
      uo(i,j)    = 0._R8
      ajphi(i,j) = 0._R8
      b(i,j)     = 0._R8
      bo(i,j)    = 0._R8
      omeg(i,j)  = 0._R8
      r(i,j)     = 0._R8
      ro(i,j)    = 0._R8
      qo(i,j)    = 0._R8
  200 continue
!
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

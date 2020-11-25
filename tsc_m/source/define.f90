      subroutine define
!......5.20 define
!***********************************************************
!                                                          *
!.....define constants and initialize counters             *
!                                                          *
!***********************************************************
      USE CLINAM
      USE POLCUR
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER nthvar
!============
      pi      = 3.1415926535897_R8
      tpi     = 2._R8*pi
      sqrt2   = sqrt(2._R8)
!
      kcycle  = -1
      time    = tpro(istart)
      times = time*udst
      ineg    = 0
      iprnt   = nskipr
      iplt    = nskipl
      dt      = 0
      nthvar  = pnsave
      nskip2  = ncycle/nthvar + 1
      iplt2   = nskip2
      iskipsf = nskipsfi
      itimets = 10
      nqflag  = 40
      itpest = 1
      mframe = 0
      nc0 = 1
      fbchi = 1._R8
      polcurchi = 0._R8
!
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

      function log10a(f1)
!
!c purpose: get an integer approximation of an antilogarithm
!
!c program statements:
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      REAL log10a
!============
! idecl:  explicitize implicit REAL declarations:
      REAL   f1,g1
!============
      g1 = f1
      if (f1 .lt. 0.0 ) g1 = -f1
      if (f1 .eq. 0.0 ) then
        print *,'log10awarning: trying to take antilog'
        print *,'of 0.0 - will return 0'
        log10a= 0
        return
      endif
!
! get antilog
!
      g1 = alog10(g1)
      log10a= 0
      if (g1 .lt. 0.0 ) log10a= g1-0.9999999999 
      if (g1 .ge. 0.0 ) log10a= g1+0.0000000001 
!
      return
      END
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

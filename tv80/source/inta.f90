      function inta (f1)
!
!c purpose: to get an integer approximation of a floater
!
!c program statements:
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER inta
!============
! idecl:  explicitize implicit REAL declarations:
      REAL   f1
!============
      inta = int(f1+0.9999999999 )
      if (f1 .lt. 0.0 ) inta = int(f1-0.9999999999 )
      return
 
      END
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

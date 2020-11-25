      function icomp (f1, f2)
!
!c purpose: to determine the numeric relationship
!           between two floaters
!
!c program statements:
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER icomp
!============
! idecl:  explicitize implicit REAL declarations:
      REAL   f1,f2,g1,g2
!============
      icomp = 0
      if (f1 .eq. f2) return
      icomp = 1
      if (f1 .lt. f2) icomp = -1
!
! see if floaters are equivalent within a reasonable value
      g1 = 1.00001 
      if (icomp .eq. 1 .and. 0.0 .le. f1                                 &  
     & .or. icomp .eq. -1 .and. f1 .lt. 0.0 )                            &  
     & g1 = 0.99999 
      g2 = g1*f1
      if (icomp .eq. 1 .and. g2 .lt. f2                                  &  
     & .or. icomp .eq. -1 .and. f2 .le. g2)                              &  
     & icomp = 0
!
      return
      END
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

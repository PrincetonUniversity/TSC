        subroutine makidmat (matrix)
 
!**********************************************************************
!
!  makidmat  -  Make an identity matrix
!
!  synopsis     call makidmat (matrix)
!               real matrix(3,3)        Matrix to change into identity
!
!**********************************************************************
 
 
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER i,j
!============
        REAL   matrix(3,3)
 
        do 10 i = 1,3
          do 20 j = 1,3
            if (i .eq. j) then
              matrix(i,j) = 1.0 
            else
              matrix(i,j) = 0.0 
            endif
20        continue
10      continue
 
        return
        end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

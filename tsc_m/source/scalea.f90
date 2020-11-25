      subroutine scalea(arrayv,arrays,n,scaledmax,arraymax,amult)
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER n,j
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 arrays,scaledmax,arraymax,amult,arrayv,amax
!============
      dimension arrayv(*),arrays(*)
!============      
!
!....scales all elements in arrayv to maximum value of scaledmax
!    returns maximum of original arrayv as arraymax*scaledmax
!
      amax = -1.E60_R8
      do 10 j=1,n
      amax = max(amax,arrayv(j))
   10 continue
      arraymax = amax/scaledmax
      if(arraymax.le.0) go to 30
      do 20 j=1,n
      arrays(j) = arrayv(j)/arraymax
   20 continue
      arraymax = arraymax*amult
      return
   30 continue
      do 40 j=1,n
   40 arrays(j) = 0._R8
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

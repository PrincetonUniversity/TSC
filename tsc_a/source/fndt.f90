      subroutine fndt(x,max,i,t)
!.....
!.....     Routine to find interval of piecewise linear data x(1:max)
!.....     within which t falls.
!.....
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER i, max
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 t,x
!============
      dimension x(max)
!.....
      i=0
      if(t.lt.x(1)) return
      i=max
      if(t.ge.x(max)) return
      do 10 i=1,max-1
        if(t.ge.x(i).and.t.lt.x(i+1)) return
   10 continue
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

      subroutine dif1(x,y,yp,n)
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER n,i
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 y,yp,x
!============
      dimension x(*),y(*),yp(*)
!============      
      do 20 i=2,n-1
        yp(i)=( y(i+1)-y(i-1) )/( x(i+1)-x(i-1) )
20    continue
      yp(1)=( y(2)-y(1) )/( x(2)-x(1) )
      yp(n)=( y(n)-y(n-1) )/( x(n)-x(n-1) )
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

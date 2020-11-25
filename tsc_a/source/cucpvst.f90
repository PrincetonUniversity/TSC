      subroutine cucpvst(tcu,cpcu)
!.....
!.....     tcu  - temperature (K)
!.....     cpcu - Specific heat of copper (J/kg/K)
!.....     Curve fit data for copper from J.H.Schultz - 6/85
!.....
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER i
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 cpcu,tcu,x,a,b
!============
      dimension x(14),a(14),b(14)
      data x/ 70._R8, 80._R8, 90._R8,100._R8,120._R8,140._R8,160._R8,    &  
     & 180._R8,200._R8,220._R8,240._R8,                                  &  
     &       260._R8,280._R8,300._R8/
      data a/173._R8,205._R8,232._R8,254._R8,288._R8,313._R8,332._R8,    &  
     & 346._R8,356._R8,364._R8,371._R8,                                  &  
     &       376._R8,381._R8,386._R8/
      data b/ 3.2_R8, 2.7_R8, 2.2_R8, 1.7_R8,1.25_R8,0.95_R8,0.70_R8,    &  
     & 0.50_R8,0.40_R8,0.35_R8,0.25_R8,                                  &  
     &       0.25_R8,0.25_R8,0.1046_R8/
!.....
      if(tcu.le.70.0_R8) then
        cpcu=173.0_R8*(tcu/70.0_R8)**3
      else
        call fndt(x,14,i,tcu)
        cpcu=a(i)+b(i)*(tcu-x(i))
      end if
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

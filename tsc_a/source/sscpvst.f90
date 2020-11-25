      subroutine sscpvst(tss,cpss)
!.....
!.....     tss  - temperature (K)
!.....     cpss - Specific heat (J/kg/K)
!.....     Curve fit data for INCONEL X-750 Handbook on Materials
!.....     for Superconducting Machinery
!.....
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER i,i1
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 cpss,tss,x,sscp,dcpdt
!============
      dimension x(8),sscp(8)
      data x   / 50._R8,100._R8,150._R8,200._R8,250._R8,300._R8,350._R8,  &  
     & 400._R8/
      data sscp/ 91._R8,230._R8,311._R8,390._R8,422._R8,459._R8,481._R8,  &  
     & 500._R8/
!.....
      call fndt(x,8,i,tss)
      i1=i+1
      if(i.eq.0) then
        i1=2
        i=1
      end if
      if(i.eq.8) then
        i1=8
        i=7
      end if
      dcpdt=(sscp(i1)-sscp(i))/(x(i1)-x(i))
      cpss=dcpdt*(tss-x(i))+sscp(i)
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

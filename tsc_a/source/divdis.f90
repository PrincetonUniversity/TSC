      subroutine divdis(n,ai,xc,zc,dis)
!
!.....called by subroutine divplat in strike point search algorithm
      USE CLINAM
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER n,i
      INTEGER int
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 ai,xc,zc,dis,dxseg,dzseg,dis1
!============
      int = ai
      if(int.gt.nseg(n)) int = nseg(n)
      if(int.lt.1    ) int = 1
      dxseg = xsega(n,int+1) - xsega(n,int)
      dzseg = zsega(n,int+1) - zsega(n,int)
      dis1 = sqrt(dxseg**2 + dzseg**2)
      xc = xsega(n,int) + (ai-int)*dxseg
      zc = zsega(n,int) + (ai-int)*dzseg
      dis = sqrt((xc-xsega(n,int))**2                                    &  
     &         + (zc-zsega(n,int))**2)
      if(int.eq.1) return
!
      do 100 i=1,int-1
      dis = dis + sqrt((xsega(n,i+1)-xsega(n,i))**2                      &  
     &               + (zsega(n,i+1)-zsega(n,i))**2)
  100 continue
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

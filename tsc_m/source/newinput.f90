!#include "f77_dcomplx.h"
      subroutine newinput
!}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}   ROS   }}}}}}}}}}}}}}}}}}}}}}
!
!.....writes new type 10 cards for restart input file
!
      USE CLINAM
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER ios,m,ngr,iabs,l
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 pprom,agroupm,cwicsm,am
      REAL*8 AREAL
!============
      if ( numargs .lt. 1 ) then
         filename = 'newtype10'
      else 
         filename = 'newtype10' // '.' // trim(suffix)
      end if
      open(50,file=trim(filename),status='unknown',iostat=ios)
      write(50,1000)
 1000 format("c.....replace existing type 10 cards with the following")
!
      do 100 m=1,nwire
      ngr = iabs(igroupw(m))
      pprom = 0
      do 50 l=1,ntpts
      pprom = pprom + fact(l)*gcur(l,ngr)*aturnsw(m)
   50 continue
      pprom = pprom +  gcurfb(ngr)*udsi*aturnsw(m)
      agroupm = AREAL(igroupw(m))
      cwicsm  = (ccoil(ncoil-nwire+m)*udsi- pprom)*1.E-3_R8
      am = m
      write(50,1001) am,xwire(m),zwire(m),agroupm,aturnsw(m),            &  
     &               rswires(m),cwicsm
  100 continue
 1001 format("10",8x,5f10.3,f10.7,f10.5)
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
! 15Apr2005 fgtok

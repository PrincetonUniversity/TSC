      subroutine jeton
!
!.....define jet source terms
      USE CLINAM
      USE SAPROP
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER j
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 acon,sumane,sumpro,dtjet,delne
!============
!     dimension xs(ppsi), apf(ppsi)
      data acon  / 0.0570814_R8/
!============      
      INTEGER :: istat = 0 
      REAL*8, ALLOCATABLE, DIMENSION(:) :: xs
      REAL*8, ALLOCATABLE, DIMENSION(:) :: apf
!============      
      IF(.not.ALLOCATED(xs)) ALLOCATE( xs(ppsi), STAT=istat)
      IF(.not.ALLOCATED(apf)) ALLOCATE( apf(ppsi), STAT=istat)
!============      
      if (istat .ne. 0) stop 'Allocation Error : jeton  ' 
!============      
!
! * * define minor radii
      do 10 j = 3,npsit
      xs(j)= sqrt((j-1.5_R8)/(npsit-1))
   10 continue
!
      do 20 j = 2,npsit
      if(j.eq.2) then
      apf(j) = 0._R8
          else
      apf(j) = (acon*(1._R8-xs(j)**2)/xs(j)**(3._R8/7._R8)) * ( 98._R8**  &  
     & (7._R8/4._R8)                                                     &  
     & - (98._R8-7._R8*xs(j)**(11._R8/7._R8)*(25._R8-11._R8*xs(j)**2))**  &  
     & (7._R8/4._R8))**(3._R8/7._R8)
          endif
   20 continue
!
      sumane  = 0._R8
      sumpro  = 0._R8
      do 30 j = 2,npsit
      sumane  = sumane + vp(j)*ane(j)
      sumpro  = sumpro + apf(j)*vp(j)
   30 continue
      dtjet   = acoef(818)-acoef(817)
      delne   = acoef(816)*(1._R8-acoef(819))                            &  
     &            * (sumane/sumpro)*(usdd/(dtjet*usdt))
      sjane0  = sumane
!
      do 40 j = 2,npsit
      sravejet(j) = delne*apf(j)
   40 continue
      sravejet(1) = sravejet(2)
!
      do 50 j = npsit+1, npsi
      sravejet(j) = 0._R8
   50 continue
!
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

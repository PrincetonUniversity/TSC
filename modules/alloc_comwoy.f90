      SUBROUTINE alloc_comwoy  
      
      USE COMWOY
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
                                                                           
      ALLOCATE ( ipow(pngroup), STAT=istat)
        ipow = 0
      ALLOCATE ( ipowm1(pngroup), STAT=istat)
        ipowm1 = 0
      ALLOCATE ( ipowm2(pngroup), STAT=istat)
        ipowm2 = 0
      ALLOCATE ( ipowp(pngroup), STAT=istat)
        ipowp = 0
      ALLOCATE ( flxp(pngroup), STAT=istat)
        flxp = 0
      ALLOCATE ( vegint(pngroup), STAT=istat)
        vegint = 0

      if (istat .ne. 0) then 
         print *, 'Allocation Error : alloc_comwoy   ' 
         stop
      end if
            
      return
      END SUBROUTINE alloc_comwoy  

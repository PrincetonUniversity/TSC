      SUBROUTINE alloc_fvv1  
      
      USE FVV1
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
        ALLOCATE ( bpolr(pnwire), STAT=istat)
        bpolr = 0 
        ALLOCATE ( bpolz(pnwire), STAT=istat)
        bpolz = 0 
        ALLOCATE ( fwirer(pnwire), STAT=istat)
        fwirer = 0 
        ALLOCATE ( fwirez(pnwire), STAT=istat)
        fwirez = 0 
        ALLOCATE ( vvdum(16), STAT=istat)
        vvdum = 0                                                          
      if (istat .ne. 0) then 
         print *, 'Allocation Error : alloc_fvv1   ' 
         stop
      end if
            
      return
      END SUBROUTINE alloc_fvv1  

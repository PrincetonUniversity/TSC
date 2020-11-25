      SUBROUTINE alloc_comsta  
      
      USE COMSTA
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
        ALLOCATE ( tpl(2*pnsave), STAT=istat)
        tpl = 0 
        ALLOCATE ( fstabpl(2*pnsave), STAT=istat)
        fstabpl = 0 
        ALLOCATE ( fdestpl(2*pnsave), STAT=istat)
        fdestpl = 0 
        ALLOCATE ( fsumpl(2*pnsave), STAT=istat)
        fsumpl = 0 
        ALLOCATE ( taupl(2*pnsave), STAT=istat)
        taupl = 0 
        ALLOCATE ( cpasopl(2*pnsave), STAT=istat)
        cpasopl = 0 
        ALLOCATE ( cpasupl(2*pnsave), STAT=istat)
        cpasupl = 0                                                        
      if (istat .ne. 0) then 
         print *, 'Allocation Error : alloc_comsta   ' 
         stop
      end if
            
      return
      END SUBROUTINE alloc_comsta  

      SUBROUTINE alloc_scadvan  
      
      USE SCADVAN
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
        ALLOCATE ( iforce(pnx,pnz), STAT=istat)
        iforce = 0 
        ALLOCATE ( dxa(pnx,pnz), STAT=istat)
        dxa = 0 
        ALLOCATE ( dxc(pnx,pnz), STAT=istat)
        dxc = 0 
        ALLOCATE ( dzb(pnx,pnz), STAT=istat)
        dzb = 0 
        ALLOCATE ( dzd(pnx,pnz), STAT=istat)
        dzd = 0 
        ALLOCATE ( face(pnx,pnz), STAT=istat)
        face = 0                                                           
      if (istat .ne. 0) then 
         print *, 'Allocation Error : alloc_scadvan   ' 
         stop
      end if
            
      return
      END SUBROUTINE alloc_scadvan  

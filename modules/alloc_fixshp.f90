      SUBROUTINE alloc_fixshp  
      
      USE FIXSHP
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
        ALLOCATE ( amatsc(pngroup,pngroup), STAT=istat)
        amatsc = 0 
        ALLOCATE ( uuscfb(pngroup,pngroup), STAT=istat)
        uuscfb = 0 
        ALLOCATE ( vvscfb(pngroup,pngroup), STAT=istat)
        vvscfb = 0 
        ALLOCATE ( bvecsc(pngroup), STAT=istat)
        bvecsc = 0 
        ALLOCATE ( gsumg(pngroup), STAT=istat)
        gsumg = 0 
        ALLOCATE ( grsumg(pngroup), STAT=istat)
        grsumg = 0 
        ALLOCATE ( gzsumg(pngroup), STAT=istat)
        gzsumg = 0 
        ALLOCATE ( grrsumg(pngroup), STAT=istat)
        grrsumg = 0 
        ALLOCATE ( grho(pngroup), STAT=istat)
        grho = 0 
        ALLOCATE ( grho2(pngroup), STAT=istat)
        grho2 = 0 
        ALLOCATE ( sigma(pngroup), STAT=istat)
        sigma = 0 
        ALLOCATE ( xvecsc(pngroup), STAT=istat)
        xvecsc = 0 
        ALLOCATE ( xvecsco(pngroup), STAT=istat)
        xvecsco = 0 
        ALLOCATE ( bvecsco(pngroup), STAT=istat)
        bvecsco = 0 
        ALLOCATE ( rv1scfb(pngroup), STAT=istat)
        rv1scfb = 0                                                        
      if (istat .ne. 0) then 
         print *, 'Allocation Error : alloc_fixshp   ' 
         stop
      end if
            
      return
      END SUBROUTINE alloc_fixshp  

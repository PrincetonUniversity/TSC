      SUBROUTINE alloc_scr3  
      
      USE SCR3
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
        ALLOCATE ( amat(pfour,pfour,ppsi), STAT=istat)
        amat = 0 
        ALLOCATE ( bmat(pfour,pfour,ppsi), STAT=istat)
        bmat = 0 
        ALLOCATE ( cmat(pfour,pfour,ppsi), STAT=istat)
        cmat = 0 
        ALLOCATE ( emat(pfour,pfour,ppsi), STAT=istat)
        emat = 0 
        ALLOCATE ( fvec(pfour,ppsi), STAT=istat)
        fvec = 0 
        ALLOCATE ( dvec(pfour,ppsi), STAT=istat)
        dvec = 0 
        ALLOCATE ( pvec(pfour,ppsi), STAT=istat)
        pvec = 0 
        ALLOCATE ( tmp1(pfour,pfour), STAT=istat)
        tmp1 = 0 
        ALLOCATE ( tmp4(pfour,pfour+1), STAT=istat)
        tmp4 = 0 
        ALLOCATE ( tmp5(pfour,pfour+1), STAT=istat)
        tmp5 = 0 
        ALLOCATE ( tmp3(pfour), STAT=istat)
        tmp3 = 0 
        ALLOCATE ( xs6(ppsi), STAT=istat)
        xs6 = 0 
        ALLOCATE ( as6(ppsi,pneq), STAT=istat)
        as6 = 0 
        ALLOCATE ( bs6(ppsi,pneq), STAT=istat)
        bs6 = 0 
        ALLOCATE ( cs6(ppsi,pneq), STAT=istat)
        cs6 = 0 
        ALLOCATE ( ds6(ppsi,pneq), STAT=istat)
        ds6 = 0                                                            
      if (istat .ne. 0) then 
         print *, 'Allocation Error : alloc_scr3   ' 
         stop
      end if
            
      return
      END SUBROUTINE alloc_scr3  

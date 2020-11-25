      SUBROUTINE alloc_iter1  
      
      USE ITER1
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
        ALLOCATE ( rz(pncoil,20), STAT=istat)
        rz = 0 
        ALLOCATE ( zz(pncoil,20), STAT=istat)
        zz = 0 
        ALLOCATE ( bpolav(pncoil), STAT=istat)
        bpolav = 0 
        ALLOCATE ( bpolavr(pncoil), STAT=istat)
        bpolavr = 0 
        ALLOCATE ( bpolavz(pncoil), STAT=istat)
        bpolavz = 0 
        ALLOCATE ( bpolmax(pncoil), STAT=istat)
        bpolmax = 0 
        ALLOCATE ( bpolr(5,pncoil,pncoil), STAT=istat)
        bpolr = 0 
        ALLOCATE ( bpolz(5,pncoil,pncoil), STAT=istat)
        bpolz = 0 
        ALLOCATE ( sumr(5), STAT=istat)
        sumr = 0 
        ALLOCATE ( sumz(5), STAT=istat)
        sumz = 0 
        ALLOCATE ( ansr(5), STAT=istat)
        ansr = 0 
        ALLOCATE ( ansz(5), STAT=istat)
        ansz = 0 
        ALLOCATE ( rfz(pncoil,5), STAT=istat)
        rfz = 0 
        ALLOCATE ( zfz(pncoil,5), STAT=istat)
        zfz = 0                                                            
      if (istat .ne. 0) then 
         print *, 'Allocation Error : alloc_iter1   ' 
         stop
      end if
            
      return
      END SUBROUTINE alloc_iter1  

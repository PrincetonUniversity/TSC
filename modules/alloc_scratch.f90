      SUBROUTINE alloc_scratch  
      
      USE SCRATCH
      
      IMPLICIT NONE
      INTEGER :: istat = 0 
      
        ALLOCATE ( bigcom(1), STAT=istat)
        bigcom = 0 
        ALLOCATE ( vecx(penx,penz), STAT=istat)
        vecx = 0 
        ALLOCATE ( vecz(penx,penz), STAT=istat)
        vecz = 0 
        ALLOCATE ( bo2(penx,penz), STAT=istat)
        bo2 = 0                                                            
      if (istat .ne. 0) then 
         print *, 'Allocation Error : alloc_scratch   ' 
         stop
      end if
            
      return
      END SUBROUTINE alloc_scratch  

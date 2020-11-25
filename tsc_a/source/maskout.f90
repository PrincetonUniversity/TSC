      subroutine  maskout
!***************************************************************************
      USE CLINAM
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER j,i
!============
        write (nout, 520)   (j, j=1,9)
  520   format(/' masking array :  iexv(i,j) = 0  says plasma can        &  
     & exist'                                                            &  
     &  /10x, 10i10)
        write (nout, 522)
  522   format (   '  j  ', 6x, 11('        0 ') )
        write (nout, 540)
        do 600 j=nzp,2,-1
        write (nout, 540)   j, zary(j), (iexv(i,j), i=2,nxp)
  540   format (i4, f6.2, 1x, 120i1)
  600   continue
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

      subroutine parmchk1
!......1.10 parmchk1
!
!
!.....writes values of all parameters as part of initialization
!.....then returns to main program
!
      USE CLINAM
      USE NONCOR
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER i
!============
      nnp(1) = iformat
      nnp(2) = pdelay
      nnp(3) = pnx
      nnp(4) = pnz
      nnp(5) = pncoil
      nnp(6) = pobs
      nnp(7) = pngroup
      nnp(8) = pnlim
      nnp(9) = ptpts
      nnp(10) = pglob
      nnp(11) = pnsave
      nnp(12) = pneq
      nnp(13) = pnthe
      nnp(14) = ppsi
      nnp(15) = pkw
      nnp(16) = pfour
      nnp(17) = plw
      nnp(18) = pimp
      nnp(19) = pne
      nnp(20) = pte
      nnp(21) = pnode
      nnp(22) = pnplat
      nnp(23) = pnseg
      nnp(24) = pnsep
      nnp(25) = pnfeed
      nnp(26) = pngrps
      nnp(27) = pnwire
!
      write(nout,1001)
 1001 format(5x," iformat pdelay pnx pnz coil pobs pngr nlim tpts",      &  
     &" glob nsav pneq nthe ppsi  pkw four  plw pimp  pne  pte ",/,12x,  &  
     &" node plat  seg  nsep numfb grps ")
      write(nout,1000) (nnp(i),i=1,27)
 1000 format(" parameters:",20i5,/,12x,7i5)
!
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

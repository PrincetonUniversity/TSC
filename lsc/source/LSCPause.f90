       SUBROUTINE LSCPause
      USE Doflags
      USE params
      USE tscunits
       IMPLICIT NONE
       CHARACTER*1 ch
       logical ::  DoPaus = .FALSE.

       if (DoPaus) then
         write(nTSCscrn,'('' ... hit <cr> to continue ...'')')
         read (5,'(a1)') ch
       endif
       return
       END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"

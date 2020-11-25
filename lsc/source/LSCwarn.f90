!     ------------------------------------------------------------------
      SUBROUTINE LSCwarn ( ErrMsg )
      USE tscunits
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!============
! idecl:  explicitize implicit INTEGER declarations:
!============
      CHARACTER*(*) ErrMsg
      INTEGER nWarning, maxWarn, i, ilen
      SAVE
      DATA    nWarning, maxWarn                                          &  
     &            / 0 ,     1000 /
      ilen = len(ErrMsg)
      nWarning = nWarning + 1
      if (nWarning .gt. maxWarn ) then
          write(nTSCscrn,'('' LSCquit: '',70a1 )')                       &  
     &                       (ErrMsg(i:i), i=1,ilen)
          call LSCstop ( ' too many warnings issued! ' )
          nWarning = 0
      else if (iRayTrsi .ne. 0) then
          write(nTSCscrn,'('' LSCwarn: '',70a1 )')                       &  
     &                       (ErrMsg(i:i), i=1,ilen)
      else
!     .                                 Do not issue warnings if only re-doing
!     .                                 the J of E calculation:  iRayTrsi=0
          return
      endif
!
!
      return
!     ------------------------------------------------------------------
      ENTRY LSCclear
      nWarning = 0
      iEndRy = 0
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"

!
      SUBROUTINE GrafParm
      USE Doflags
      USE DqlBins
      USE FeBins
      USE Jrf
      USE params
      USE PlPr
      USE Ramppwr
      USE RayBins
      USE RayWrk
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      call GrafFixP('BEGIN PLOT', 0    )
      call GrafFltP(' fghz     ', fghz)
      call GrafFltP(' HstpLH   ', HstpLH)
      call GrafFixP(' nstep    ', nstep)
      call GrafFixP(' nrays    ', nrays)
      call GrafFixP(' ntors    ', ntors)
      call GrafFixP(' npols    ', npols)
      call GrafFixP(' nGrps    ', nGrps)
      call GrafFixP(' npsi     ', npsi)
      call GrafFixP(' nzones   ', nzones)
      call GrafFixP(' nv       ', nv)
      call GrafFixP(' nsmoo    ', nsmoo)
      call GrafFixP(' nsmw     ', nsmw)
      call GrafFixP(' nRampUp  ', nRampUp)
      call GrafFixP(' nFlat    ', nFlat)
      call GrafFixP(' idiag(1) ', idiag(1))
      call GrafFixP(' idiag(2) ', idiag(2))
      call GrafFltP(' WeghtItr ', WeghtItr)
      call GrafFltP(' ScatKdeg ', ScatKdeg)
      call GrafFltP(' DiffuJrf ', DiffuJrf)
      call GrafFltP(' PrfSpred ', PrfSpred)
      call GrafFixP(' DoTRAN   ', DoTRAN)
      call GrafFixP(' DoXcam   ', DoXcam)
      call GrafFixP(' Do1Rpr   ', Do1Rpr)
      call GrafFixP(' Do0Edc   ', Do0Edc)
!
      if ( DoBram .eq. 0 ) then
      call GrafFixP(' DoBram   ', DoBram)
      call GrafFltP(' nparmax  ', nparmax)
      call GrafFltP(' nparmin  ', nparmin)
      call GrafFltP(' npolmax  ', npolmax)
      call GrafFltP(' npolmin  ', npolmin)
      call GrafFixP(' nGrps    ', nGrps  )
      call GrafFltP(' centers-1', centers(1))
      call GrafFltP(' widths-1 ', widths (1))
      call GrafFltP(' powers-1 ', powers (1))
      if ( nGrps .ge. 2 ) then
      call GrafFltP(' centers-2', centers(2))
      call GrafFltP(' widths-2 ', widths (2))
      call GrafFltP(' powers-2 ', powers (2))
      endif
      if ( nGrps .ge. 3 ) then
      call GrafFltP(' centers-3', centers(3))
      call GrafFltP(' widths-3 ', widths (3))
      call GrafFltP(' powers-3 ', powers (3))
      endif
      endif
!
      if ( DoBram .eq. 1 ) then
      call GrafFixP(' DoBram   ', DoBram)
      call GrafFixP(' nslices  ', nslices)
      call GrafFixP(' nGrps    ', nGrps  )
      call GrafChrS(' couplers1', couplers(1))
      call GrafFltP(' paseDeg-1', phaseDeg(1))
      call GrafFltP(' powers-1 ', powers (1))
      if ( nGrps .ge. 2 ) then
      call GrafChrS(' couplers2', couplers(2))
      call GrafFltP(' paseDeg-2', phaseDeg(2))
      call GrafFltP(' powers-2 ', powers (2))
      endif
      if ( nGrps .ge. 3 ) then
      call GrafChrS(' couplers3', couplers(3))
      call GrafFltP(' paseDeg-3', phaseDeg(3))
      call GrafFltP(' powers-3 ', powers (3))
      endif
      endif
!
      call GrafFixP(' TurnNegs ', TurnNegs)
      if (TailTeps*TailPeps*TailNeps .gt. 0.00_R8) then
      call GrafFltP(' TailTeps ', TailTeps)
      call GrafFltP(' TailPeps ', TailPeps)
      call GrafFltP(' TailNeps ', TailNeps)
      endif
      call GrafFltP(' Vmin     ', Vmin)
      call GrafFltP(' Vmax     ', Vmax)
      call GrafFixP(' nfreq    ', nfreq)
      call GrafFixP(' lfast    ', lfast)
!     call GrafFltP(' thet0    ', thet0)
      call GrafFltP(' thet0unus', thet0)
!
      call GrafFixP(' PPSI     ', PPSI)
      call GrafFixP(' PNX      ', PNX )
      call GrafFixP(' PNZ      ', PNZ)
      call GrafFixP(' PIMP     ', PIMP )
      call GrafFixP(' PWORDS   ', PWORDS)
      call GrafFixP(' NPSIDIM  ', NPSIDIM)
      call GrafFixP(' NRAYDIM  ', NRAYDIM)
      call GrafFixP(' NTORDIM  ', NTORDIM)
      call GrafFixP(' NPOLDIM  ', NPOLDIM)
      call GrafFixP(' NVELDIM  ', NVELDIM)
      call GrafFixP(' NZONDIM  ', NZONDIM)
      call GrafFixP(' NPLTDIM  ', NPLTDIM)
      call GrafFixP(' NGRPDIM  ', NGRPDIM)
 
      call MkGrfLst(' List of LSC inputs, params   ')
      call GrafFixP('END PLOT  ', 0   )
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"

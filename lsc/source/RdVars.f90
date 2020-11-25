!                                                                      |
!                                                                      |
!     block data program                    ---------------------------|
 
!     RdVars    Read Variables, sets constants ------------------------|
!                                                                      |
!                                                                      |
      SUBROUTINE RdVars
      USE params
      USE CGSetc
      USE Inpval
      USE MKSetc
      USE PIetc
      USE ProfBody
      USE tscunits
      USE xparams
      USE Xray
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER iabs
!============
      CHARACTER*12 FileNa1
      INTEGER ifirst,i
      REAL*8    tmpdum
      DATA    ifirst /1/
      DATA    FileNa1 /'input.lhh'/
!
!     NAMELIST /inpvalue/  goes here from the include file:
!     NAMELIST /inpexprt/  goes here from the include file:
!     .                                 Reset nrays,ntors,npols so that
!     .                                 the error checker can detect whether
!     .                                 ntors&npols are basic inputs or
!     .                                 nrays is basic input.  This for back-
!     .                                 compatibility with input files.
      if (DoTRAN .ne. TRUE) then
         nrays = -1
         ntors = -1
         npols = -1
         TailPeps = 0.00_R8
         TailTeps = 0.00_R8
         TailNeps = 0.00_R8
      endif
!     .
      if (DoTRAN .ne. TRUE) then
         open (nTSCunus,status='old',file = FileNa1, err=1300 )
         read (nTSCunus, inpvalue)
         read (nTSCunus, inpexprt)
         close (nTSCunus)
      endif
 
 12   begin = 0._R8
 
      woc2 = (TWOPI * fghz / CLIGHT)**2 * 1.E2_R8
      woc4 =  woc2**2
      omega=  TWOPI*fghz * 1.E9_R8
 15   continue
!
      if (DiffuJrf .lt. 0.0000_R8) then
              DiffuJrf = 0.00_R8
              call LSCwarn(' DiffuJrf set to 0')
              endif
      if (PrfSpred .lt. 0.0000_R8.or.                                    &  
     &   (DiffuJrf .eq. 0.0000_R8.and.                                   &  
     &    PrfSpred .gt. 0.0000_R8))then
              PrfSpred = 0.00_R8
              call LSCwarn(' PrfSpred set to 0')
              endif
      if (PrfSpred .gt. 1.0000_R8) then
              PrfSpred = 1.00_R8
              call LSCwarn(' PrfSpred set to 1')
              endif
      if (nparmin .ge. nparmax .and. DoBram .eq. 0) then
              tmpdum = nparmax
              nparmax = nparmin+0.001_R8
              nparmin = tmpdum -0.001_R8
              call LSCwarn(' nparmin/max reversed')
              endif
!     Check dimensions now
      if ( npsi   .gt. NPSIDIM   ) then
              npsi   = NPSIDIM
              call LSCwarn ( ' npsi set to NSPIDIM ')
              endif
!
!     .                                 The if-else-if sequence below
!     .                                 is to maintain compatibility w.
!     .                                 old usage npols=1,ntors=nrays,
!     .                                 ntors not explicitly used;
!     .                                 IF ntors&nrays not given then
!     .                                 use NTORDIM
 
      if ( ntors  .gt. NTORDIM   ) then
              ntors  = NTORDIM
              call LSCwarn ( ' ntors set to NTORDIM ')
              endif
      if ( npols  .gt. NPOLDIM   ) then
              npols  = NPOLDIM
              call LSCwarn ( ' npols set to NPOLDIM ')
              endif
      if ( nrays  .gt. NRAYDIM   ) then
              nrays  = NRAYDIM
              call LSCwarn ( ' nrays set to NRAYDIM ')
              endif
      if ( npols  .lt. 1         ) then
              npols  = 1
              endif
      if ( nrays .eq. -1 .and. ntors .eq. -1 ) then
                  npols = 1
                  ntors = NTORDIM
                  nrays = NTORDIM
!     .                                 IF ntors given calculate nrays
        else if ( nrays  .eq. -1 ) then
                  nrays  = ntors*npols
!     .                                 IF nrays given calculate ntors
        else if ( ntors .eq. -1)    then
                  ntors  = nrays/npols
                  if ( ntors  .gt. NTORDIM   ) then
                     ntors  = NTORDIM
                     nrays  = ntors*npols
                     call LSCwarn ( ' ntors set to NTORDIM ')
                  endif
      endif
!
      if (nrays .ne. ntors*npols) then
                  nrays = ntors*npols
                  call LSCwarn ( ' nrays set to ntors*npols ')
              endif
!
      if ( nv     .gt. NVELDIM     ) then
              nv     = NVELDIM
              call LSCwarn ( ' nv set to NVELDIM ')
              endif
      if ( nzones .gt. NZONDIM   ) then
              nzones = NZONDIM
              call LSCwarn ( ' nzones set to NZONDIM ')
              endif
      if ( nRampUp  .gt. NRAMPDIM ) then
              nRampUp = NRAMPDIM
              call LSCwarn ( ' nRampUp set to NRAMPDIM ')
              endif
      if ( nsmoo  .gt. nv/3   .or.                                       &  
     &     nsmoo  .gt. nzones/3    ) then
                call LSCwarn (' nsmoo set to min(nv,nzones)/3 ')
                nsmoo = min(nv, nzones) / 3
                endif
!                                       The relative size of these dimensions
!                                       is implicitly assumed in work arrays.
!     if ( NVELDIM .lt. (NPSIDIM + NZONDIM)) then           ! NWKVDIM avoids
!               call LSCstop (' NVELDIM < NPSIDIM+NZONDIM' )! this clumsiness
!               endif                             ! =NVELDIM+NPSIDIM+NZONDIM
!
      if(mod(nsmoo, 2) .eq. 0) then
                call LSCwarn (' nsmoo MUST BE ODD ')
                nsmoo = iabs(nsmoo - 1)
                endif
      if ( nsmw .lt. nsmoo/8 ) then
                call LSCwarn (' nsmoo-width seems too SMALL ')
                endif
      if ( nsmw .ge. nsmoo   ) then
                call LSCwarn (' nsmoo-width seems too LARGE ')
                endif
      if ( WeghtItr .gt. 1.0_R8.or. WeghtItr .lt. 0.0_R8) then
                WeghtItr = 0.5_R8
                call LSCwarn (' WeghtItr set to 0.50 ')
                endif
      if ( TailTeps .gt. 0.00_R8.and. TailNeps .gt. 0.00_R8) then
        TailPeps = TailNeps/TailTeps
      else if (TailTeps .gt. 0.00_R8.and. TailPeps .gt. 0.00_R8) then
        TailNeps = TailPeps*TailTeps
      else if (TailNeps .gt. 0.00_R8.and. TailPeps .gt. 0.00_R8) then
        TailTeps = TailNeps/TailPeps
      else
        TailTeps = 0.00_R8
        TailNeps = 0.00_R8
        TailPeps = 0.00_R8
      endif
      if (TailTeps .gt. 0.3_R8.or.                                       &  
     &    TailPeps .gt. 0.3_R8.or.                                       &  
     &    TailNeps .gt. 0.1_R8) then
        call LSCwarn(' Too much fast electron tail')
        TailTeps = 0.00_R8
        TailNeps = 0.00_R8
        TailPeps = 0.00_R8
      endif
      if ( nGrps .gt. NGRPDIM .or. nGrps .le. 0 ) then
                nGrps = 1
                call LSCwarn (' nGrps being set to 1 ')
                endif
      do i=1,nGrps
         if(powers(i) .le. 0.00_R8) then
            powers(i) = 0.001_R8
            call LSCwarn( ' powers(i) being set to 0.001 ')
         endif
      enddo
 
      if (DoBram .eq. 1 ) then
        do i=1,nGrps
 
          if (                                                           &  
     & couplers(i).ne.'PBXMFAST'.and.couplers(i).ne.'PBXMSLOW'.and.      &  
     & couplers(i).ne.'SLOWSLOW'.and.couplers(i).ne.'TORSUPRA'.and.      &  
     & couplers(i).ne.'TOKDEVAR'.and.couplers(i).ne.'TFTRLHCD'.and.      &  
     & couplers(i).ne.'JET_LHCD'.and.couplers(i).ne.'USRSPEC1'.and.      &  
     & couplers(i).ne.'USRSPEC2'.and.couplers(i).ne.'USRSPEC3')then
              couplers(i)    = 'PBXMSLOW'
              call LSCwarn (' couplers(i) being set to PBXMSLOW ')
          endif
 
        enddo
      endif
 
      if(ifirst .eq. 1) then
         ifirst= 0
         PIO4  = atan(1._R8)
         PI    = 4._R8* PIO4
         TWOPI = 2._R8* PI
         RTPI  = sqrt(PI)
         dompesqdn = 4._R8* PI * (ECHARG/1.E-10_R8)**2 /                 &  
     &                  (EMASS/1.E-28_R8) *                              &  
     &                                                        1.E8_R8
         RESTenergy = (EMASS/1.E-28_R8)*(CLIcgs/1.E10_R8)**2 *1.E-8_R8
         ERGperEV = 1.6E-12_R8
         keVperERG = 1._R8/ (1000._R8* ERGperEV)
         if (DoTRAN .eq. 0 ) then
           call GrafParm
         endif
      endif
 
 
      return
 1300 call LSCstop (' LSC cant find input.lhh ')
      return
      END
!     -                                                                |
!     RdVars    ends                        ---------------------------|
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"

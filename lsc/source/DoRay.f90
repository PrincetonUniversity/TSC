 
!                                                                      |
!     Finished Initializing ray zones                ------------------|
!
!                                                                      |
!                                                                      |
!     Raystore.F                        -------------------------------|
!     Raytrace.F  begins                    ---------------------------|
!
!     DoRay     begins                      ---------------------------|
!                                                                      |
!                                                                      |
      SUBROUTINE DoRay(iDoRay)
      USE dielec
      USE Doflags
      USE params
      USE PlPr
      USE RayBins
      USE RayWrk
      USE tscunits
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      EXTERNAL RayIni, PredcLSC, EZfini
      INTEGER iry, ity, ipy, RayIniErr, iFillPsi, iDoRay, iErrCount
      iFillPsi =1
!     iDoRay : if it has a value 1 thru nrays==npols*ntors
!              then trace just the ray with that index
!              if it is any other value, then trace all rays
!
      iErrCount = 1
      iry = 0
      do 70 ity =1,ntors
      do 69 ipy =1,npols
        iry=iry+1
        iray = iry
        if ( iDoRay .eq. iray .or.                                       &  
     &       iDoRay .lt. 1    .or.                                       &  
     &       iDoRay .gt. nrays     ) then
          mstpLH = nfreq
          call RyZnInit    !! 9sep 93
          iplt = 0
          lnewray = iray
          enpar = ntor (ity)
          enpol = npol (ipy)
          izone = 1
          call RayIni  (RayIniErr)
          if (RayIniErr .ge. 1) then
             iErrCount = iErrCount+1
             go to 69
          endif
          call PredcLSC
!         call LSCPause
          if(PlFlg(RAYPL) .eq. TRUE)then
             call PlRay(iFillPsi)
                iFillPsi=0
                if (iError .gt. 0)  return
          endif
        endif
 69   continue
 70   continue
      if (iErrCount .ge. ntors*npols) then
        iError = 1
        call LSCtrace('RayIni')
      endif
      return
      END
!                                                                      |
!                                                                      |
!     DoRay     ends                        ---------------------------|
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"

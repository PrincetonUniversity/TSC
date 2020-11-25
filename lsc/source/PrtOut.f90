 
 
 
!     Raystore.F                        -------------------------------|
!                                                                      |
!                                                                      |
!
!     PrtOut                                ---------------------------|
!                                                                      |
!                                                                      |
      SUBROUTINE PrtOut(iprint)
      USE dielec
      USE Doflags
      USE params
      USE PIetc
      USE PlPr
      USE RayBins
      USE RayWrk
      USE tscunits
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
            EXTERNAL plasma2d, DispRela, E2byPr
      INTEGER icount, iprint
      REAL*8    NperFs, NperSl
!     If iprint = 0, print only if the count is reached
!                 1, print at this point regardless of the count
!                                       If starting a new ray, set
!                                       the print counter for headings
!                                       icount
!                                       to 0.
      if ( lnewray .gt. 0 ) icount = 0
 
        call E2byPr
 
      mstpLH = mstpLH + 1
      if ( iprint .eq. 1 ) go to 15
      if ( mstpLH .lt. nfreq ) return
      mstpLH = 0
 
15    continue
      iplt = iplt + 1
      if ( iplt .gt.IDM ) iplt = IDM
      yrr(iplt)   = y(1)
      yzz(iplt)   = y(2)
      ypp(iplt)   = y(3)
      yps(iplt)   = rtPsRy(izone)
      ypr(iplt)   = NperRy(izone)
      ypl(iplt)   = npar  (izone,iray)
      yed(iplt)   = power (izone,iray)
 
      if ( mod(icount,10) .eq. 0 .and. PrFlg(RAYWR) .ge. TRUE ) then
            write(nLSCcomm, 610) iray
      endif
      icount = icount + 1
!
!
      if (PrFlg(RAYWR) .ge. TRUE )                                       &  
     &                    write (nLSCcomm, 611)                          &  
     & izone,izind(izone,iray), yps(iplt),                               &  
     & npar(izone,iray), NperRy(izone),                                  &  
     & ezsq(izone,iray),                                                 &  
     & dtdV, DetrRy(izone), DistRy(izone)
 
!izone izind  rpsi  npar  nper     Ez2    Err2    Ert2    dtdV     det path ln
!  199   101 +0.99 -10.1 +100. +6.e-09 -1.e-00 -1.e-00 +1.e-00 -1.e-05  26.001
 
!izone izind  rpsi  npar  nper     Ez2    dtdV     det path ln
!  199   101 +0.99 -10.1 +100. +6.e-09 +1.e-00 -1.e-05  26.001
 
 
 610  format     (' ray',i2,/,                                           &  
     &' izone izind  rpsi  npar  nper',                                  &  
     &'     Ez2    dtdV     det path ln' )
 
 611  format     (                                                       &  
     &1x,i5,                                                             &  
     &1x,i5,                                                             &  
     &1x,f5.2,                                                           &  
     &1x,f5.1,                                                           &  
     &1x,f5.0,                                                           &  
     &1x,1pe7.0,                                                         &  
     &1x,1pe7.0,                                                         &  
     &1x,1pe7.0,                                                         &  
     &   0pf8.3   )
 
!     diagnostic prints
!
!izone izind        epsZ        epsL      NperRy      NperFs    NperSl     dtdV
!  199   101 +1.2345e-05 +1.2345e-12 +1.2345e-12 +1.2345e-11 +1.23e-11 1.23e-00
 1311 format( 1x,i5,i5,5(1pe11.4,1x),e10.2)
 1312 format('^',i5,i5,5(1pe11.4,1x),e10.2)
 1310 format(' ray',i2,/,                                                &  
     &' izone izind      epsZ        epsL      NperRy',                  &  
     &'      NperFs      NperSl       dtdV'  )
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"

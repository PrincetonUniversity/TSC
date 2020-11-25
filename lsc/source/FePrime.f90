!     -----------------------------------------------------------------
      SUBROUTINE FePrime
      USE DqlBins
      USE FeBins
      USE params
      USE ProfBody
      USE RayBins
      USE WkAry
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER ips, iv
      INTEGER iNotSmoo, iYesSmoo
      DATA    iNotSmoo, iYesSmoo / 1 , 2 /
      EXTERNAL BrodCast
!
      call FePru
!     compute unsmoothed derivative
!
 
!
      do 10 ips = 1, npsi
        do 5 iv = 1, nv
          dfdv(iv,ips,iYesSmoo) = vpar(iv)*dfdv(iv,ips,iNotSmoo)
  5     continue
 10   continue
 
      call Smooth(dfdv(1,1,iYesSmoo),                                    &  
     &             nv, NVELDIM, npsi, nsmoo, qlsm)
 
!     convolve with smoothing function
!
!     NOTE: For energy conservation, the smoothing done here must
!           be identical to that done in constructing Dql.
!           Conceptually equivalent to inclusion of resonance broadening
!           in wave-particle interaction.  Loops added at 5,14,15,16
!           make smoothing over v df/dv, rather than over df/dv alone.
!
      do 20 ips = 1, npsi
         do 14 iv = 1, ivZero-1
           dfdv(iv,ips,iYesSmoo) = dfdv(iv,ips,iYesSmoo)/vpar(iv)
 14      continue
 15      dfdv(ivZero,ips,2) = 0.00_R8
         do 16 iv = ivZero+1 , nv
           dfdv(iv,ips,iYesSmoo) = dfdv(iv,ips,iYesSmoo)/vpar(iv)
 16      continue
 20   continue
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"

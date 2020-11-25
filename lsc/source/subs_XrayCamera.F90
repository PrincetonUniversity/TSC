!#include "f77_dcomplx.h"
#define COMPILE_X_RAY_CAMERA
#ifdef COMPILE_X_RAY_CAMERA
      SUBROUTINE XrayCam2
      USE params
      USE xparams
      USE RayBins
      USE FeBins
      USE ProfBody
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      call xRayInp
      call xBremSetup(ZbrAry(1))
      call xRaySetup
      call xIntens(NVELDIM,NPSIDIM,iITR,fe, nv,                          &  
     &             npsi, NeAry, ZbrAry, Vpar)
      call xEmit
!
      return
      END
!
!     ------------------------------------------------------------------
      SUBROUTINE xRayInp
      USE params
      USE TSCgrap
      USE tscunits
      USE emparams
      USE xparams
      USE Xray
      USE emitter
      USE camera
      USE numerics
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!     NameList begins                       ---------------------------|
!                                                                      |
!                                                                      |
      NAMELIST /inpxry/                                                  &  
     &        dFoilTCM, FoilCode, iAbsXray,                              &  
     &        nEbins, nMUbins, E_max, E_min, nr_source, nz_source,       &  
     &        Z_bound_min, Z_bound_max, R_bound_max, R_bound_min,        &  
     &        n_pixel_x, n_pixel_y,                                      &  
     &        E_ph_min, dE_ph,                                           &  
     &        Rpinhole, Zpinhole, phi_pinhole, pinhole_size,             &  
     &        focal_length, screen_d, R_tangent, z_tangent,              &  
     &        dlxray, npoints_max,                                       &  
     &        PusherMajor, PusherMinor
!                                                                      |
!                                                                      |
!     NameList ends                         ---------------------------|
!
!      include "emitter.inc"
!
      Z_plasm_max = max(ZlcfsMax,-ZlcfsMin)
      Z_plasm_min = -Z_plasm_max
      R_plasm_max = RlcfsMax
      R_plasm_min = RlcfsMin
 
      Z_bound_max = Z_plasm_max*1.005_R8
      Z_bound_min =-Z_plasm_max*1.005_R8
      R_bound_max = RlcfsMax*1.005_R8
      R_bound_min = min(RlcfsMin/1.005_R8,R_bound_min)
!
!      include "camera.inc"
!
!
!      include "numerics.inc"
      open (nTSCunus,status='old',file = 'input.xry', err=1300 )
      read (nTSCunus, inpxry)
      close (nTSCunus)
!
      goto 1301
 1300 call LSCwarn (' cant find input.xry; use defaults! ')
 1301 continue
      if (n_pixel_x .gt. NPIXDIM) then
        call LSCwarn('n_pixel_x reset to NPIXDIM')
        n_pixel_x = NPIXDIM
        endif
 
      if (n_pixel_y .gt. NPIXDIM) then
        call LSCwarn('n_pixel_y reset to NPIXDIM')
        n_pixel_y = NPIXDIM
        endif
 
      if (nr_source .gt. NRDIM) then
        call LSCwarn('nr_source reset to NRDIM')
        nr_source = NRDIM
        endif
 
      if (nz_source .gt. NZDIM) then
        call LSCwarn('nz_source reset to NZDIM')
        nz_source = NZDIM
        endif
 
      if (nMUbins .gt. NMUDIM) then
        call LSCwarn('nMUbins reset to NMUDIM')
        nMUbins = NMUDIM
        endif
 
      if (2*(nMUbins/2) .eq. nMUbins) then
        call LSCwarn('nMUbins must be odd; subtracting 1')
        nMUbins = nMUbins-1
        endif
 
      if (nEbins .gt. NENDIM) then
        call LSCwarn('nEbins reset to NENDIM')
        nEbins = NENDIM
        endif
 
      if (E_max .gt. 0.255_R8) then
        call LSCwarn('E_max reset to 0.255 MeV')
        E_max = 0.255_R8
        endif
!
      pix_size_x = screen_d/n_pixel_x
      pix_size_y = screen_d/n_pixel_y
!
      call GrafFixP('BEGIN PLOT', 0    )
      call GrafFixP(' NRDIM    ', NRDIM)
      call GrafFixP(' NZDIM    ', NZDIM)
      call GrafFixP(' NRZMXDIM ', NRZMXDIM)
      call GrafFixP(' NPIXDIM  ', NPIXDIM)
      call GrafFixP(' NRCHORDIM', NCHORDIM)
      call GrafFixP(' MAXPOINTS', MAXPOINTS)
      call GrafFixP(' SPACEDIM ', SPACEDIM)
      call GrafFltP(' ACC_ERROR', ACCEPTABLE_ERROR)
      call GrafFixP(' NMUDIM   ', NMUDIM)
      call GrafFixP(' NENDIM   ', NENDIM)
      call GrafFixP(' iAbsXray ', iAbsXray)
      call GrafChrS(' FoilCode ', FoilCode)
      call GrafFltP(' dFoilTCM ', dFoilTCM)
      call GrafFixP(' nMUbins  ', nMUbins)
      call GrafFixP(' nEbins   ', nEbins)
      call GrafFltP(' mu_min   ', mu_min)
      call GrafFltP(' mu_max   ', mu_max)
      call GrafFltP(' E_min    ', E_min)
      call GrafFltP(' E_max    ', E_max)
      call GrafFltP(' E_ph_min ', E_ph_min)
      call GrafFltP(' dE_ph    ', dE_ph)
      call GrafFixP(' nr_source', nr_source)
      call GrafFixP(' nz_source', nz_source)
      call GrafFltP(' mu_0     ', mu_0)
      call GrafFltP(' mu_width ', mu_width)
      call GrafFltP(' ZboundMin', Z_bound_min)
      call GrafFltP(' ZboundMax', Z_bound_max)
      call GrafFltP(' RboundMin', R_bound_min)
      call GrafFltP(' RboundMax', R_bound_max)
      call GrafFltP(' ZplasmMin', Z_plasm_min)
      call GrafFltP(' ZplasmMax', Z_plasm_max)
      call GrafFltP(' RplasmMin', R_plasm_min)
      call GrafFltP(' RplasmMax', R_plasm_max)
      call GrafFltP(' ZlcfsMin ', ZlcfsMin)
      call GrafFltP(' ZlcfsMax ', ZlcfsMax)
      call GrafFltP(' RlcfsMin ', RlcfsMin)
      call GrafFltP(' RlcfsMax ', RlcfsMax)
      call GrafFltP(' PusherMaj', PusherMajor)
      call GrafFltP(' PusherMin', PusherMinor)
      call GrafFixP(' n_pixel_x', n_pixel_x)
      call GrafFixP(' n_pixel_y', n_pixel_y)
      call GrafFltP(' focal_len', focal_length)
      call GrafFltP(' screen_d ', screen_d)
      call GrafFltP(' pinhole_s', pinhole_size)
      call GrafFltP(' R_tangent', R_tangent)
      call GrafFltP(' Z_tangent', Z_tangent)
      call GrafFltP(' RcntChrd ', RcntChrd)
      call GrafFltP(' Rpinhole ', Rpinhole)
      call GrafFltP(' Zpinhole ', Zpinhole)
      call GrafFltP(' pix_sizex', pix_size_x)
      call GrafFltP(' pix_sizey', pix_size_y)
      call GrafFixP(' npointsMx', npoints_max)
      call GrafFltP(' dlxray   ', dlxray)
!     call GrafFixP(' .', .)
!     call GrafFltP(' .', .)
!     call GrafChrS(' .', .)
      call MkGrfLst(' List X ray Camera parameters ')
      call GrafFixP('END PLOT  ', 0   )
!
      return
      END
!     ..................................................................
      BLOCK DATA xRayBk
      USE params
      USE xparams
      USE emparams
      USE Xray
 
      USE emitter
 
      USE camera
 
      USE numerics
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
 
!     DATA mu_min, mu_max / -1._R8, 1._R8/
!     DATA nMUbins / NMUDIM /
!     DATA nEbins  / NENDIM /
!     DATA E_ph_min / 1.E-6_R8/
!     DATA dE_ph / .01_R8/
!     DATA E_min, E_max / 0.01_R8, 0.2_R8/
!     DATA Z_bound_min, Z_bound_max /                                    &  
!    &           -1.00_R8,        1.00_R8/
!     DATA R_bound_min, R_bound_max /                                    &  
!    &            1.10_R8,        2.00_R8/
!     DATA mu_0 / 1._R8/ mu_width / .2_R8/
!     DATA nr_source / NRDIM / nz_source / NZDIM /
!
!     DATA Rpinhole / 2.661_R8/
!     DATA Zpinhole /   0._R8/
!     DATA phi_pinhole / 0._R8/
!     DATA pinhole_size / 0.00625_R8/
!     DATA focal_length, screen_d  /                                     &  
!    &          0.4905_R8,  0.215_R8/
!     DATA n_pixel_x, n_pixel_y  / NPIXDIM, NPIXDIM /
!     DATA R_tangent / 1.524_R8/
!     DATA RcntChrd  / 1.65_R8/
!     DATA z_tangent / 0._R8/
!     DATA pix_size_x, pix_size_y / 0.01_R8, 0.01_R8/
!
!     DATA npoints_max / MAXPOINTS /
!     DATA eps / ACCEPTABLE_ERROR /
!     DATA dlxray / 0.005_R8/
!     DATA  PusherMajor, PusherMinor / 1.152_R8, 0.190_R8/
!     DATA FoilCode / 'AG' /
!     DATA dFoilTCM / 0.000_R8/
!     DATA iAbsXray / 1 /
!
      END
!     .         --------------------------------------------------------
      SUBROUTINE xRaySetup
      USE params
      USE FeBins
      USE xparams
      USE Xray
      USE emparams
      USE emitter
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      REAL*8 AREAL
!     generate grids, cross-section database
!     and adjust the sigtot by the transmission and illumination factor
!
      call ugrid(muvec, nMUbins, mu_min, mu_max)
      dmu_inv = AREAL(nMUbins - 1) / (mu_max - mu_min)
!
      call ugrid(E_incvec, nEbins, E_min, E_max)
      call xSetupEindex
!     set up indices for table lookup into sigma
      call xGenEofv(nv, Eofv, dEofv, Vpar)
      if (iAbsXray .gt. 0) then
        call GetXmnFac(nEbins, E_incvec, XmnFac, dFoilTCM, FoilCode)
!        call lscpause
!        call vwrite(nEbins, E_incvec, 'E_incvec')
!        call vwrite(nEbins, XmnFac, 'XmnFac')
!        call lscpause
      endif
!     compute vector XmnFac of absorption coefficients
      call xSigma
 
      return
      END
!     .         --------------------------------------------------------
      SUBROUTINE xGenEofv(n, E, dE, v)
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER n, i
      REAL*8   signorm, Mevpermcsq
      COMMON / sigbrcom0 / signorm, Mevpermcsq
      REAL*8   E(n), dE(n), v(n)
      do i = 1, n
         E(i) = 0.5_R8* Mevpermcsq * v(i) * v(i)
      enddo
      do i = 2, n - 1
         dE(i) = 0.5_R8* abs(E(i + 1) - E(i - 1))
!     accounts for negative velocity
      enddo
      dE(1) = 0.5_R8* abs(E(2) - E(1))
      dE(n) = 0.5_R8* abs(E(n) - E(n - 1))
      return
      END
!     .         --------------------------------------------------------
      SUBROUTINE xSetupEindex
      USE params
      USE xparams
      USE FeBins
      USE Xray
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER i, j
      REAL*8   Einc, MevPerMcsq
      PARAMETER ( MevPerMcsq = 0.511_R8)
 
      do i = 1, nv
         Einc = 0.5_R8* Mevpermcsq * Vpar(i) * Vpar(i)
         j = 1
         do while((E_incvec(j) .le. Einc) .and. (j .lt. nEbins))
            j = j + 1
         enddo
         j = j - 1
         Eint(i) = j
         if(j .ge. 1 .and. j .lt. nEbins) then
           Efrac(i) = (Einc - E_incvec(j)) / (E_incvec(j + 1) -          &  
     &        E_incvec(j))
         else
           Efrac(i) = 0.00_R8
         endif
      enddo
      return
      END
!     .         --------------------------------------------------------
      REAL*8 FUNCTION Xraytransmission(EmeV)
      USE params
      USE xparams
      USE Xray
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER ie
      REAL*8 frac, EmeV
      if(iAbsXray .gt. 0)then
         call Lookupindex(EmeV, E_incvec, nEbins, ie, frac)
         Xraytransmission = XmnFac(ie) * (1._R8- frac) +                 &  
     &        XmnFac(ie + 1) * frac
      else
         Xraytransmission = 1._R8
      endif
      return
      END
!     .         --------------------------------------------------------
      SUBROUTINE Lookupindex(c, vec, n, i, frac)
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER n, i
      REAL*8                                                             &  
     &        c, vec(n), frac, ri
      REAL*8    RE41, AREAL
      ri = (c - vec(1)) / (vec(2) - vec(1))
!     i = ifix(ri)
      RE41 = ri
      i = int(RE41)
      frac = ri - AREAL(i)
      if(i .lt. 1)then
         i = 1
         frac = 0._R8
      else if(i .ge. n)then
         i = n - 1
         frac = 1._R8
      endif
      return
      END
!     .         --------------------------------------------------------
      SUBROUTINE xSigma
      USE params
      USE xparams
      USE Xray
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      REAL*8    sigg, E_photon, Xraytransmission
      INTEGER  jj, kk
      do 10 kk = 1, nEbins
 
         do 20 jj = 1, nMUbins
            sigtot(jj,kk) = 0._R8
            E_photon = E_ph_min + 0.5_R8* dE_ph
            do while(E_photon .lt. E_incvec(kk))
               call xSigGen(E_incvec(kk), E_photon, muvec(jj),           &  
     &              sigg)
               sigtot(jj,kk) = sigtot(jj,kk) + sigg * E_photon * dE_ph   &  
     &              * Xraytransmission(E_photon       )
!     weighted to yield rate of energy emission.
!
               E_photon = E_photon + dE_ph
            enddo
 20      continue
 10   continue
 
!c    Ernies contour plot of sigtot contours 28Mar93 begin
!     call X_pl_intensity_cntrs(sigtot, NMUDIM, nMUbins, NEbins)
!     call EZfini(0,0)
!c    Ernies contour plot of sigtot contours 28Mar93 begin
 
      return
      END
 
!     +---------------------------------------------------------------+
 
      SUBROUTINE xIntens(NVDIM,NPDIM,iIterate,FeOfv, nv, np,             &  
     &     NeAry, ZbrAry, Vpar)
!     intensity subroutine calculates the value of the photon intensity
!     over a range of scattering angles mu and for each psi flux surface
!     Actually, the values given are I(mu,psi)/Clight.
 
      USE params
      USE xparams
      USE Xray
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER NPDIM, np, NVDIM, nv, iv, j, ip, iIterate
      INTEGER iv0
      INTEGER ith
      REAL*8    FeOfv(NVDIM, NPDIM, 2), sig_th(NMUDIM), Ifactor
      REAL*8    NeAry(NPDIM), ZbrAry(NPDIM), Vpar(NVDIM)
!     PARAMETER (MevPerMcsq = 0.511) delete this; add lines below
      REAL*8   signorm, Mevpermcsq
      COMMON / sigbrcom0 / signorm, Mevpermcsq
      REAL*8                                                             &  
     &       ONE
      DATA   ONE/                                                        &  
     &       1.0_R8/
      iv0 = (nv + 1) / 2
!     velocity index for which Vpar(iv0) = 0
!     mu is measured relative to positive Vpar axis.
!     sigtot   INPUT  differential cross section times photon energy
!                     times transmission factor (.le. 1.0)
!                     w.r.t. solid angle(str) in
!                     MeV - cm^2 per atom per incident electron
!     sig_th   OUTPUT differential cross section times photon energy
!                     times transmission factor (.le. 1.0)
!                     w.r.t. solid angle(str) interpolated velocity in
!                     MeV - cm^2 per atom per incident electron
!     inten    OUTPUT joules/sec per cm^3 per steradian
 
      call BrodCast(NMUDIM * NPDIM, inten, 0._R8)
      do 100 iv = 1, nv
!.
        do 10 ith = 1, nMUbins
          if(Eint(iv) .ge. 1 .and. Eint(iv) .lt. nEbins) then
            sig_th(ith) = sigtot(ith, Eint(iv)) * (1._R8- Efrac(iv)) +   &  
     &                    Efrac(iv) *sigtot(ith, Eint(iv) + 1)
          else
            sig_th(ith) = 0.00_R8
          endif
 10     continue
!.
        if(iv .le. iv0) then
          call xVecRefl(nMUbins, sig_th)
          Ifactor = 0.5_R8* abs((Vpar(iv + 1) + Vpar(iv)) *              &  
     &         (Vpar(iv + 1) - Vpar(iv)))
       else
          Ifactor = 0.5_R8* abs((Vpar(iv) + Vpar(iv - 1)) *              &  
     &         (Vpar(iv) - Vpar(iv - 1)))
        endif
!     note: dEofv ALWAYS .gt. 0, and Vpar(iv) increases from (-1., 1.),
!
        do 50 ip = 1, np
          do 40 j = 1, nMUbins
            inten(j, ip) = inten(j, ip) +                                &  
     &                     FeOfv(iv, ip, iIterate) * sig_th(j) *         &  
     &                     Ifactor
 40      continue
 50     continue
 100  continue
!
      do 150 ip = 1, np
        do 140 j = 1, nMUbins
          inten(j, ip) = inten(j, ip) * NeAry(ip)                        &  
     &                   / max(ONE ,ZbrAry(ip)) * 4.8E-3_R8
!    ^                   / max(1.00,ZbrAry(ip)) * 4.8e-3
 
 140    continue
 150  continue
      return
      END
 
!     +---------------------------------------------------------------+
 
      SUBROUTINE xVecRefl(n, vec)
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER i, n, i_comp
      REAL*8    vec(n), vec_tmp
      do i = 1, (n + 1) / 2
         i_comp = n + 1 - i
         vec_tmp = vec(i_comp)
         vec(i_comp) = vec(i)
         vec(i) = vec_tmp
      enddo
      return
      END
 
!     +---------------------------------------------------------------+
 
      REAL*8   FUNCTION delta_em(rr, zz, mu)
      USE params
      USE RayWrk
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER ip, imu
      REAL*8    rr, zz, mu, ip_frac, imu_frac
      call xGetPsi(rr, zz, mu, ip, ip_frac, imu, imu_frac)
      call xPsiIntp(ip, ip_frac, imu, imu_frac, delta_em)
      delta_em = delta_em * Rmaj / rr
!     to account for convergence of flux tubes
      return
      END
!     .         --------------------------------------------------------
      SUBROUTINE xGetPsi(r, z, mu, ip, ip_frac, imu, imu_frac)
!     find psi of r,z
      USE params
      USE RayBins
      USE ProfBody
      USE TSCgrap
      USE emparams
      USE xparams
      USE Xray
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER ip, imu, jlo, jhi
      REAL*8    mu, ip_frac, imu_frac, rmu
      REAL*8    r,z,psi,Br,Bz,RBphi,omc,Tee,pe2,pi2,aio,ael
      REAL*8 AREAL
      DATA    jlo / -1 /
      call plasma2d(r,z,psi,Br,Bz,RBphi,omc,Tee,pe2,pi2,aio,ael)
      if (psi .ge. PsiAry(npsi) .or. pe2 .le. 0.00_R8.or.                &  
     &      r .ge. RlcfsMax     .or.   r .le. RlcfsMin .or.              &  
     &      z .le. ZlcfsMin     .or.   z .ge. ZlcfsMax ) then
        ip      = 0
        ip_frac = 0.00_R8
        imu     = 1
        imu_frac= 0.00_R8
        return
      endif
      call huntnr(PsiAry, npsi, psi, jlo)
      ip = jlo
      jhi= jlo+1
      if (jhi .gt. npsi) then
        ip_frac = 0.00_R8
      else
        ip_frac = (psi - PsiAry(jlo))/(PsiAry(jlo+1) - PsiAry(jlo))
      endif
      rmu = (mu - mu_min) * dmu_inv + 1.00_R8
      imu = int(rmu)
      imu_frac = rmu - AREAL(imu)
      return
      END
 
      SUBROUTINE xPsiIntp(ip, ip_frac, imu, imu_frac, delta_em)
      USE params
      USE RayBins
      USE xparams
      USE Xray
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER ip, imu, ip_plus, imu_plus
      REAL*8    ip_frac, imu_frac, delta_em
!     linear intepolation
      imu_plus = min(imu + 1, nMUbins)
      ip_plus = min(ip + 1, npsi)
      if (ip .le. 0 ) then
        delta_em = 0.00_R8
      else
        delta_em = (1._R8- ip_frac - imu_frac) * inten(imu, ip)          &  
     &     + ip_frac * inten(imu, ip_plus)                               &  
     &     + imu_frac *inten(imu_plus, ip)
      endif
      return
      END
!
!     emit.F          --------------------------------------------------
!
!     produce bremstrahlung emission pattern from prescribed source
!     as observed by X-ray pinhole camera
 
      SUBROUTINE xEmit
!NCAR      call opngks
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
           call EZinit
      call X_setup_camera
      call X_compute_pattern
      call X_diagnose_trajectories
      call X_plot_pixel_data
      call X_plot_source_profile
           call MkGrfLst('X ray contour plots')
           call EZfini(0,0)
!!SLICES HERE?
           call EZinit
      call xSlices
      call eSlices
!     call XmnGraf
 
!     call X_wrres
!NCAR      call clsgks
           call EZfini(0,0)
 
      return
      END
!     .         --------------------------------------------------------
      SUBROUTINE X_rzphi_to_xyz(xyz, rzphi)
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      REAL*8   xyz(*), rzphi(*)
!     transform from cylindrical to cartesian coordinates
      xyz(3) = rzphi(2)
      xyz(1) = rzphi(1) * cos(rzphi(3))
      xyz(2) = rzphi(1) * sin(rzphi(3))
      return
      END
!     .         --------------------------------------------------------
      SUBROUTINE X_xyz_norm(norm, xyz)
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      REAL*8   norm, xyz(*)
      norm = sqrt(xyz(1) * xyz(1) + xyz(2) * xyz(2) + xyz(3) *           &  
     &     xyz(3))
      return
      END
!     .         --------------------------------------------------------
      SUBROUTINE X_rzphi_norm(norm, rzphi)
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      REAL*8   norm, rzphi(*)
      norm = sqrt(rzphi(1) * rzphi(1) + rzphi(2) * rzphi(2))
      return
      END
!     .         --------------------------------------------------------
      SUBROUTINE X_det_33(det, mat)
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      REAL*8   det, mat(3, *)
      det = mat(1, 1) * mat(2, 2) * mat(3, 3) +                          &  
     &     mat(1, 2) * mat(2, 3) * mat(3, 1) + mat(1, 3) *               &  
     &     mat(2, 1) * mat(3, 2) - mat(1, 1) * mat(3, 2) *               &  
     &     mat(2, 3) - mat(2, 1) * mat(1, 2) * mat(3, 3) -               &  
     &     mat(3, 1) * mat(2, 2) * mat(1, 3)
      return
      END
!     .         --------------------------------------------------------
      SUBROUTINE X_xyz_to_rzphi(rzphi, xyz)
      USE numerics
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      REAL*8   xyz(*), rzphi(*)
!     transform from cartesian to cylindrical coordinates
      rzphi(2) = xyz(3)
      rzphi(1) = sqrt(xyz(1) * xyz(1) + xyz(2) * xyz(2))
      if(rzphi(1) .eq. 0._R8)then
         rzphi(3) = 0._R8
      else if(xyz(1) .gt. 0._R8)then
         rzphi(3) = asin(xyz(2) / rzphi(1))
      else
         rzphi(3) = pi - asin(xyz(2) / rzphi(1))
      endif
      return
      END
!     .         --------------------------------------------------------
      SUBROUTINE X_setup_camera
      USE emparams
      USE emitter
      USE camera
      USE numerics
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
 
      INTEGER ix, iy, i, j
      REAL*8   chi, mag_xy, del_R, del_z, rpy_sq(NPIXDIM),               &  
     &     FocLen2, axis_orientation(3, 3),                              &  
     &     chord_orientation(3), x_pix_min, x_pix_max, y_pix_min,        &  
     &     y_pix_max, rpx_sq, rp_sq, phi_pix, mu_pix, xyz_pinhole(3),    &  
     &     distance_pix_to_pinhole, det, norm, y_hat_xyz_pinhole(3)
      REAL*8    rt_rp_sq, asinarg, SMALL, AREAL
      DATA    SMALL /1.0E-10_R8/
 
      PI = 4._R8* atan(1._R8)
      R_bound_min_sq = R_bound_min * R_bound_min
      R_bound_max_sq = R_bound_max * R_bound_max
 
      call X_rzphi_to_xyz(xyz_pinhole, pinhole_loc)
!
!     begin by:
!     determining orientation of camera symmetry axis relative to
!     orientation of unit vector from pinhole (R, z, phi) to
!     tangency point prescribed by R_tangent, z_tangent
!
!     set up camera axis orientation -- spherical angles
      call X_normvec_xyz(3, y_hat_xyz_pinhole, xyz_pinhole)
      call X_xyz_norm(norm, y_hat_xyz_pinhole)
      chi = phi_pinhole + acos(R_tangent / Rpinhole)
      del_R = Rpinhole * cos(pi / 2._R8- (chi - phi_pinhole))
      del_z = (z_tangent - Zpinhole)
      mag_xy = sqrt(del_R * del_R + del_z * del_z)
      phi_axis = pi / 2._R8+ chi
      mu_axis = acos(del_z / mag_xy)
      axis_orientation(1, 3) = sin(mu_axis) * cos(phi_axis)
      axis_orientation(2, 3) = sin(mu_axis) * sin(phi_axis)
      axis_orientation(3, 3) = cos(mu_axis)
      axis_orientation(1, 2) = sin(phi_axis)
      axis_orientation(2, 2) = - cos(phi_axis)
      axis_orientation(3, 2) = 0._R8
      axis_orientation(1, 1) = - cos(phi_axis) * cos(mu_axis)
      axis_orientation(2, 1) = - sin(phi_axis) * cos(mu_axis)
      axis_orientation(3, 1) = sin(mu_axis)
      call X_det_33(det, axis_orientation)
!     gives the (x, y, z) projections of the three camera axes
!
!     next:
!     set up to compute orientation (cartesian coordinates) of chords
!     originating from each element
!
!     first step:
!     determine angular offset, relative to camera axis, of chords
!     eminating from each CCD element
!     orientation of axes is that positive x axis lies in vertical plane
!     formed by camera symmetry axis and tokamak centerline (z axis)
!     in camera coordinates then:
!     symmetry axis is polar axis, x axis is origin for phi and
!     y = z cross vx
!     the pixels are assumed located on the x-y plane, a distance
!     focal_length from the pinhole, with dimensions
!     pix_sixe_x and pix_size_y
      FocLen2 = focal_length * focal_length
      x_pix_min = - 0.5_R8* AREAL(n_pixel_x) * pix_size_x
      x_pix_max = - x_pix_min
      y_pix_min = - 0.5_R8* AREAL(n_pixel_y) * pix_size_y
      y_pix_max = - y_pix_min
      call ugrid(x_pix_vec, n_pixel_x, x_pix_min, x_pix_max)
      call ugrid(y_pix_vec, n_pixel_y, y_pix_min, y_pix_max)
      do iy = 1, n_pixel_y
         rpy_sq(iy) = y_pix_vec(iy) * y_pix_vec(iy)
      enddo
      do ix = 1, n_pixel_x
         rpx_sq = x_pix_vec(ix) * x_pix_vec(ix)
         do iy = 1, n_pixel_y
            rp_sq = rpx_sq + rpy_sq(iy)
            rt_rp_sq = sqrt(rp_sq) * ( 1.0_R8+ SMALL )
!           mu_pix = asin(sqrt(rp_sq / (rp_sq + FocLen2)))
            asinarg = sqrt( rp_sq / (rp_sq + FocLen2) )
            mu_pix = asin(asinarg)
            if(rp_sq .eq. 0._R8)then
               phi_pix = 0._R8
            else if(x_pix_vec(ix) .ge. 0._R8)then
               asinarg = y_pix_vec(iy) / rt_rp_sq
!              call asinDebg(asinarg,ix,iy,n_pixel_x,n_pixel_y,2)
!              phi_pix = asin(y_pix_vec(iy) / sqrt(rp_sq))
               phi_pix = asin(asinarg)
            else
               asinarg = y_pix_vec(iy) / rt_rp_sq
!              call asinDebg(asinarg,ix,iy,n_pixel_x,n_pixel_y,3)
!              phi_pix = pi - asin(y_pix_vec(iy) / sqrt(rp_sq))
               phi_pix = pi - asin(asinarg)
            endif
!     generate the projection of chord orientation onto camera
!     axes
            phi_pix = pi + phi_pix
!     chord originates at pixel plane
            chord_orientation(1) = sin(mu_pix) * cos(phi_pix)
            chord_orientation(2) = sin(mu_pix) * sin(phi_pix)
            chord_orientation(3) = cos(mu_pix)
            call X_xyz_norm(norm, chord_orientation)
            distance_pix_to_pinhole = sqrt(rp_sq + FocLen2)
            do i = 1, 3
               y_hat_camera(i, ix, iy) = 0._R8
               do j = 1, 3
                  y_hat_camera(i, ix, iy) = y_hat_camera(i, ix, iy) +    &  
     &                 axis_orientation(i, j) * chord_orientation(j)
!     in cartesian coordinates
               enddo
            enddo
            call X_xyz_norm(norm, y_hat_camera(1, ix, iy))
            do i = 1, 3
               chord_origin(i, ix, iy) = xyz_pinhole(i) -                &  
     &              distance_pix_to_pinhole * y_hat_camera(i, ix, iy)
!     in cartesian coordinates
            enddo
         enddo
      enddo
      return
      END
!     .         --------------------------------------------------------
!      SUBROUTINE asinDebg(asinarg,ix,iy,n_pixel_x,n_pixel_y,icall)
!      if (asinarg .gt. 1.0 ) then
!        write(6,'('' asinarg,ix,iy,npixx,npixy,icall:''
!     ^  e20.13,5i4)')asinarg,ix,iy,n_pixel_x,n_pixel_y,icall
!        asinarg=+1.00
!      endif
!      if (asinarg .lt.-1.0 ) then
!        write(6,'('' asinarg,ix,iy,npixx,npixy,icall:''
!     ^  e20.13,5i4)')asinarg,ix,iy,n_pixel_x,n_pixel_y,icall
!        asinarg=-1.00
!      endif
!      return
!      END
!     .         --------------------------------------------------------
      SUBROUTINE X_normvec_xyz(n, y_hat, y)
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER n
      REAL*8   norm, y_hat(n), y(n), ynorminv
      call X_xyz_norm(norm, y)
      ynorminv = 1._R8/ norm
      call X_vsmult(n, y_hat, y, ynorminv)
      return
      END
!     .         --------------------------------------------------------
      SUBROUTINE X_vsmult(n, vout, vin, c)
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER n, i
      REAL*8   vout(n), vin(n), c
      do i = 1, n
         vout(i) = c * vin(i)
      enddo
      return
      END
!     .         --------------------------------------------------------
      SUBROUTINE X_diagnose_trajectories
      USE emparams
      USE camera
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER ix, iy
      INTEGER pix_count(NPIXDIM)
      INTEGER iDBG
      DATA iDBG/0/
      if (iDBG .eq. 0 ) return
      do ix = 1, NPIXDIM
         pix_count(ix) = ix
      enddo
      write(4, 101)(pix_count(iy), pix_count(iy), iy = 1, n_pixel_y)
 101  format(10(i5, 1x, i5, 2x))
      write(4, 102)
 102  format(/)
      do ix = 1, n_pixel_x
            write(4, 100)(r_tangent_actual(ix, iy),                      &  
     &           z_tangent_actual(ix, iy), iy = 1, n_pixel_y)
 100        format(10(f5.1, 1x, f5.1, 2x))
      enddo
      return
      END
!     .         --------------------------------------------------------
      SUBROUTINE X_compute_pattern
      USE emparams
      USE camera
      USE numerics
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER ix, iy
      REAL*8   y_beg_rzphi(3), y_hat_xyz(3), y_end(3), SolAng
      solAng = 3.14159_R8*(pinhole_size/focal_length)**2 * 0.25E+2_R8
      do ix = 1, n_pixel_x
         do iy = 1, n_pixel_y
            call X_init_chord(ix, iy, y_beg_rzphi, y_hat_xyz)
            call X_c_phot(y_beg_rzphi, y_hat_xyz, dlxray, y_end, 3,      &  
     &           photon_count(ix, iy), npoints(ix, iy),                  &  
     &           chord(1, 1, ix, iy), r_tangent_actual(ix, iy),          &  
     &           z_tangent_actual(ix, iy))
            photon_count(ix, iy) =  photon_count(ix, iy)*SolAng
         enddo
      enddo
!      call LSCpause
!      do ix = 1, n_pixel_x
!         do iy = 1, n_pixel_y
!
!        if (pc(ix,iy) .gt. 1.e+20 .or.
!     ^      pc(ix,iy) .lt.-1.e+20  ) then
!            write(6,'('' pc is bad at i,j;'', e10.3, i4, i4)')
!     ^      pc(ix,iy), ix, iy
!            pc(ix,iy) = 0.
!         endif
!         enddo
!      enddo
!
!      call LSCpause
 
      return
      END
!     .         --------------------------------------------------------
      SUBROUTINE X_init_chord(ix, iy, y_beg_rzphi, y_hat_xyz)
      USE emparams
      USE camera
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER ix, iy
      REAL*8    y_beg_rzphi(*), y_hat_xyz(*)
      call X_xyz_to_rzphi(y_beg_rzphi, chord_origin(1, ix, iy))
      call xVecCopy(3, y_hat_xyz, y_hat_camera(1, ix, iy))
      return
      END
!     .         --------------------------------------------------------
      INTEGER FUNCTION inbound(rsq, rr, zz)
      USE emparams
      USE emitter
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      REAL*8   rsq, rr, zz
      if((rsq .ge. R_bound_min_sq) .and. (rsq .le. R_bound_max_sq)       &  
     &     .and.                                                         &  
     &     (zz .ge. Z_bound_min) .and. (zz .le. Z_bound_max)             &  
     &     .and.                                                         &  
     &     ( (rr-PusherMajor)**2 + zz**2 .ge. PusherMinor**2) ) then
         inbound = 1
      else
         inbound = 0
      endif
      return
      END
!     .         --------------------------------------------------------
      SUBROUTINE X_plot_source_profile
      USE emparams
      USE emitter
      USE camera
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER iz, ir
      REAL*8   delta_em, rr, zz
      call ugrid(r_source, nr_source, R_plasm_min, R_plasm_max)
      call ugrid(z_source, nz_source, Z_plasm_min, Z_plasm_max)
      do iz = 1, nz_source
         zz = z_source(iz)
         do ir = 1, nr_source
            rr = r_source(ir)
            source_profile(ir, iz) = delta_em(rr, zz, mu_0)
         enddo
      enddo
      call X_pl_emisivity_cntrs(source_profile, NRDIM, nr_source,        &  
     &     nz_source,r_source,z_source)
      return
      END
!     .         --------------------------------------------------------
      SUBROUTINE X_plot_pixel_data
      USE emparams
      USE emitter
      USE camera
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER nr, DOCHORDPLOT
      REAL*8    Radii(2), Z_height
      DATA DOCHORDPLOT / 0 /
!     plot chords in (x, y) and (r, z) views
      nr = 2
      Radii(1) = R_bound_min
      Radii(2) = R_bound_max
      Z_height = Z_bound_max
      if(DOCHORDPLOT .eq. 1 ) then
      call X_plan_view(nr, Radii, npoints, n_pixel_x, n_pixel_y, chord,  &  
     &     MAXPOINTS, NPIXDIM)
      call X_elevation_view(nr, Radii, - Z_height, Z_height, npoints,    &  
     &     n_pixel_x, n_pixel_y, chord, MAXPOINTS, NPIXDIM)
      endif
      call X_pl_intensity_cntrs(photon_count, NPIXDIM, n_pixel_x,        &  
     &     n_pixel_y)
      return
      END
!     .         --------------------------------------------------------
      SUBROUTINE xSlices
      USE emparams
      USE emitter
      USE camera
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!     DATA RcntChrd  / 1.65 / ! the nominal radius of the usual central
!                             ! chord plotted in the time development of the
!                             ! xray signal
      CHARACTER*70 MyString
      INTEGER i,j,icenter,jcenter, ioffset, joffset
      REAL*8    Tan2Pin, SHAVE
      REAL*8    slice(NPIXDIM,2),                                        &  
     &        signal(NPIXDIM, NCHORDIM), sigMax, xmaxy
      DATA    SHAVE / 0.985_R8/
      REAL*8    ZERO, ONE, AREAL
      DATA    ZERO, ONE /                                                &  
     &        0.0_R8, 1.0_R8/
 
      Tan2Pin = sqrt ( Rpinhole**2 - R_tangent**2 )
!23456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789
      joffset = n_pixel_y *                                              &  
     &  abs(RcntChrd-R_tangent) / Tan2Pin * focal_length / (screen_d)
      ioffset = n_pixel_x *                                              &  
     &  abs(RcntChrd-R_tangent) / Tan2Pin * focal_length / (screen_d)
      icenter = (n_pixel_x + 1)/2
      jcenter = (n_pixel_y + 1)/2
!
!
      sigMax=0.00_R8
      do 10 i=1, n_pixel_x
        slice(i,1) = AREAL(i)
        signal(i,1) = photon_count(i,(jcenter))
        signal(i,2) = photon_count(i,(jcenter-joffset))
        signal(i,3) = photon_count(i,(jcenter-joffset/2))
        signal(i,4) = photon_count(i,(jcenter+joffset))
        signal(i,5) = photon_count(i,(jcenter+joffset/2))
        sigMax= max(sigMax, signal(i,1) )
        sigMax= max(sigMax, signal(i,2) )
        sigMax= max(sigMax, signal(i,3) )
        sigMax= max(sigMax, signal(i,4) )
        sigMax= max(sigMax, signal(i,5) )
 10   continue
      do 11 j=1, n_pixel_y
        slice(j,2) = AREAL(j)
        signal(j,6) = photon_count(icenter,j)
        sigMax= max(sigMax, signal(j,6))
 11   continue
      call XRopen(8,'XRinput')
      call XRdscw(slice(1,1), signal(1,2), n_pixel_x,                    &  
     &            'Vertical x ray slice')
      call XRclos(8)
 
      do 15 i=1, n_pixel_x
        signal(i,1) = signal(i,1)/(sigMax+1.E-30_R8)*SHAVE
        signal(i,2) = signal(i,2)/(sigMax+1.E-30_R8)*SHAVE
        signal(i,3) = signal(i,3)/(sigMax+1.E-30_R8)*SHAVE
        signal(i,4) = signal(i,4)/(sigMax+1.E-30_R8)*SHAVE
        signal(i,5) = signal(i,5)/(sigMax+1.E-30_R8)*SHAVE
 15   continue
      do 16 j=1, n_pixel_y
        signal(j,6) = signal(j,6)/(sigMax+1.E-30_R8)*SHAVE
 16   continue
!
!     ULH corner
      call EZrnd2 (AREAL(n_pixel_x),xmaxy,i,j)
      call EZsets(100,500,450,700,                                       &  
     & ZERO,xmaxy,ZERO, ONE, 1)
!    ^ 0.,xmaxy,0.00, 1.00, 1)
      call EZaxes(1,2,1,5)
      call EZcurv(slice(1,1),signal(1,1),n_pixel_x)
      write(MyString,'(''Vert slice /'',1pe9.2,                          &  
     &    '' at tangent$'')') sigMax
      call EZwrit(100,725,MyString,0,0)
!     URH corner
      call EZsets(600,1000,450,700,                                      &  
     & ZERO,xmaxy,ZERO, ONE, 1)
!    ^ 0.,xmaxy,0.00, 1.00, 1)
      call EZaxes(1,2,1,5)
      call EZcurv(slice(1,1),signal(1,2),n_pixel_x)
      call EZcros(slice(1,1),signal(1,3),n_pixel_x)
 
      write(MyString, '(''at chord'',i3,'' and'',i3,'' (dots)$'')')      &  
     &                   (jcenter-joffset), (jcenter-joffset/2)
      call EZwrit(600,725,MyString,0,0)
!
!     LLH corner
      call EZsets(100,500,100,350,                                       &  
     & ZERO,xmaxy,ZERO, ONE, 1)
!    ^ 0.,xmaxy,0.00, 1.00, 1)
      call EZaxes(1,2,1,5)
      call EZcurv(slice(1,1),signal(1,4),n_pixel_x)
      call EZcros(slice(1,1),signal(1,5),n_pixel_x)
      write(MyString, '(''at chord'',i3,'' and'',i3,'' (dots)$'')')      &  
     &                   (jcenter+joffset), (jcenter+joffset/2)
      call EZwrit(100,375,MyString,0,0)
!
!     LRH corner
      call EZrnd2 (AREAL(n_pixel_y),xmaxy,i,j)
      call EZsets(600,1000,100,350,                                      &  
     & ZERO,xmaxy,ZERO, ONE, 1)
!    ^ 0.,xmaxy,0.00, 1.00, 1)
      call EZaxes(1,2,1,5)
      call EZcurv(slice(1,2),signal(1,6),n_pixel_y)
      write(MyString,'(''Horiz slice '',''at mid plane$'')')
      call EZwrit(600,375,MyString,0,0)
!
           call MkGrfLst('X ray slices; 4 graphs')
           call EZfini(0,0)
!
!
      return
      END
!     .         --------------------------------------------------------
      SUBROUTINE eSlices
      USE emparams
      USE emitter
!      USE camera
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      CHARACTER*70 MyString
      INTEGER i,j,icenter,jcenter, ioffset, joffset
      REAL*8    SHAVE
      REAL*8    signal(NRZMXDIM, NCHORDIM), sigMax
      REAL*8    rminy, rmaxy, rdumy, zmaxy
      DATA    SHAVE / 0.985_R8/
      REAL*8    ZERO, ONE
      DATA    ZERO, ONE /                                                &  
     &        0.0_R8, 1.0_R8/
 
      ioffset = nr_source / 4
      joffset = nz_source / 4
      icenter = (nr_source + 1)/2
      jcenter = (nz_source + 1)/2
!
!
      sigMax=0.00_R8
      do 10 j=1, nz_source
        signal(j,1) = source_profile((icenter          ),j)
        signal(j,2) = source_profile((icenter+ioffset  ),j)
        signal(j,3) = source_profile((icenter+ioffset/2),j)
        signal(j,4) = source_profile((icenter-ioffset  ),j)
        signal(j,5) = source_profile((icenter-ioffset/2),j)
        sigMax= max(sigMax, signal(j,1) )
        sigMax= max(sigMax, signal(j,2) )
        sigMax= max(sigMax, signal(j,3) )
        sigMax= max(sigMax, signal(j,4) )
        sigMax= max(sigMax, signal(j,5) )
 10   continue
      do 11 i=1, nr_source
        signal(i,6) = source_profile(i,jcenter)
        sigMax= max(sigMax, signal(i,6))
 11   continue
 
      do 15 j=1, nz_source
        signal(j,1) = signal(j,1)/(sigMax+1.E-30_R8)*SHAVE
        signal(j,2) = signal(j,2)/(sigMax+1.E-30_R8)*SHAVE
        signal(j,3) = signal(j,3)/(sigMax+1.E-30_R8)*SHAVE
        signal(j,4) = signal(j,4)/(sigMax+1.E-30_R8)*SHAVE
        signal(j,5) = signal(j,5)/(sigMax+1.E-30_R8)*SHAVE
 15   continue
      do 16 i=1, nr_source
        signal(i,6) = signal(i,6)/(sigMax+1.E-30_R8)*SHAVE
 16   continue
!
!     ULH corner
      call EZrnd2( z_source(nz_source),zmaxy,i,j )
      call EZsets(100,500,450,700,                                       &  
     & -zmaxy,zmaxy,ZERO, ONE, 1)
!    ^ -zmaxy,zmaxy,0.00, 1.00, 1)
      call EZaxes(i,j,1,5)
      call EZcurv(z_source,signal(1,1),nz_source)
      write(MyString,'(''Vert slice /'',1pe9.2,                          &  
     &    '' at center$'')') sigMax
      call EZwrit(100,725,MyString,0,0)
!
!     URH corner
      call EZsets(600,1000,450,700,                                      &  
     & -zmaxy,zmaxy,ZERO,  ONE, 1)
!    ^ -zmaxy,zmaxy,0.00, 1.00, 1)
      call EZaxes(i,j,1,5)
      call EZcurv(z_source,signal(1,2),nz_source)
      call EZcros(z_source,signal(1,3),nz_source)
 
      write(MyString, '(''at r_loc'',i3,'' and'',i3,'' (dots)$'')')      &  
     &                   (icenter+ioffset), (icenter+ioffset/2)
      call EZwrit(600,725,MyString,0,0)
!
!     LLH corner
      call EZsets(100,500,100,350,                                       &  
     & -zmaxy,zmaxy,ZERO,  ONE, 1)
!    ^ -zmaxy,zmaxy,0.00, 1.00, 1)
      call EZaxes(i,j,1,5)
      call EZcurv(z_source,signal(1,4),nz_source)
      call EZcros(z_source,signal(1,5),nz_source)
      write(MyString, '(''at r_loc'',i3,'' and'',i3,'' (dots)$'')')      &  
     &                   (icenter-ioffset), (icenter-ioffset/2)
      call EZwrit(100,375,MyString,0,0)
!
!     LRH corner
      call EZrnd2 (r_source(nr_source), rmaxy, i,j)
      call EZrnd2 (rmaxy-r_source(1),rdumy, i,j)
      rminy = rmaxy - rdumy
      call EZsets(600,1000,100,350,                                      &  
     & rminy,rmaxy,ZERO,  ONE, 1)
!    ^ rminy,rmaxy,0.00, 1.00, 1)
      call EZaxes(i,j,1,5)
      call EZcurv(r_source,signal(1,6),nr_source)
      write(MyString,'(''Horiz slice '',''at mid plane$'')')
      call EZwrit(600,375,MyString,0,0)
!
           call MkGrfLst('Emissivity slices; 4 graphs')
           call EZfini(0,0)
!
      return
      END
!     .         --------------------------------------------------------
!      SUBROUTINE XmnGraf
!      include "xparams.inc"
!      include "Xray.inc"
 
!...
!      return
!      END
!     .         --------------------------------------------------------
      SUBROUTINE X_wrres
      USE params
      USE emparams
      USE xparams
      USE Xray
      USE emitter
      USE camera
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!     write results to disk file
 
      INTEGER ix, iy
      INTEGER iDBG
      DATA iDBG/1/
      if (iDBG .eq. 0 ) return
!.
      write(4, 102) 'photon_c'
 102  format(/,a10)
      do ix = 1, NPIXDIM
      write(4, 101)(photon_count(ix,iy), iy = 1, NPIXDIM)
 101  format(8(1pe10.1))
      enddo
      write(4, 102)
!.
      write(4, 102) 'source_p'
      do ix = 1, NRDIM
      write(4, 101)(source_profile(ix,iy), iy = 1, NZDIM)
      enddo
      write(4, 102)
!.
      write(4, 102) 'sigtot'
      do ix = 1, NMUDIM
      write(4, 101)(sigtot(ix,iy), iy = 1, NENDIM)
      enddo
      write(4, 102)
!.
      write(4, 102) 'inten'
      do ix = 1, NMUDIM
      write(4, 101)(inten(ix,iy), iy = 1, NPSIDIM)
      enddo
      write(4, 102)
!.
 
      return
      END
!     .         --------------------------------------------------------
      SUBROUTINE X_c_phot(y_beg, y_hat, dlchord, y_end, n_space, pc,     &  
     &     npoints, chord, r_tangent_actual, z_tangent_actual)
      USE emparams
      USE numerics
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER n_space, inbound, i, count, state,                         &  
     &     current_status, previous_status
!     INTEGER*4 npoints
      INTEGER   npoints
      REAL*8                                                             &  
     &     y_beg(n_space), dlchord, y(SPACEDIM), y_hat(SPACEDIM),        &  
     &     delta_em, y_end(n_space), mu_yb, dy_xyz(3),                   &  
     &     r_tangent_actual, z_tangent_actual, radius_sq,                &  
     &     radius_sq_min, rr, zz
!     REAL*4 chord(3, *), pc
      REAL*8   chord(3, *), pc
      call X_rzphi_to_xyz(y, y_beg)
      do i = 1, 3
         dy_xyz(i) = y_hat(i) * dlchord
      enddo
      pc = 0._R8
      count = 0
      state = 1
      radius_sq = y(1) * y(1) + y(2) * y(2)
      rr = sqrt (radius_sq)
      radius_sq_min = radius_sq
      previous_status = inbound(radius_sq, rr , y(3))
      current_status = previous_status
      do while((state .lt. 3) .and. (count .lt. npoints_max))
         do i = 1, 3
            chord(i, count + 1) = y(i)
            y(i) = y(i) + dy_xyz(i)
         enddo
         radius_sq = y(1) * y(1) + y(2) * y(2)
         rr = sqrt (radius_sq)
         if(radius_sq_min .gt. radius_sq)then
            radius_sq_min = radius_sq
            z_tangent_actual = y(3)
         endif
         count = count + 1
         if(state .eq. 2)then
!           rr = sqrt(radius_sq) ! we already did this
            zz = y(3)
            pc = pc + delta_em(rr, zz, mu_yb(y, y_hat)) * dlchord
         endif
         previous_status = current_status
         current_status = inbound(radius_sq, rr, y(3))
         if(current_status .ne. previous_status)  then
            state = state + 1
         endif
      enddo
      r_tangent_actual = sqrt(radius_sq_min)
      npoints = count
      call X_xyz_to_rzphi(y_end, y)
      return
      END
!     .         --------------------------------------------------------
      SUBROUTINE X_advance_y(y, dy_xyz)
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER i
      REAL*8   y(*), dy_xyz(*), y_old(3)
      call X_rzphi_to_xyz(y_old, y)
!     call X_project_dy(y, dy_xyz, dy_rzphi)
      do i = 1, 3
         y(i) = y_old(i) + dy_xyz(i)
      enddo
      return
      END
!     .         --------------------------------------------------------
      SUBROUTINE X_project_dy(y_rzphi, dy_xyz, dy_rzphi)
      USE emparams
      USE numerics
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!     given y (R, z, phi) and dy (dx, dy, dz), compute (dR, dz, dphi)
      INTEGER i, j
      REAL*8   y_rzphi(*), dy_xyz(*), dy_rzphi(*), T(3, 3), dRdx, dRdy,  &  
     &     dRdz, dzdx, dzdy, dzdz, dphidx, dphidy, dphidz,               &  
     &     xx, yy, Rsq_inv
      EQUIVALENCE                                                        &  
     &     (T(1, 1), dRdx), (T(1, 2), dRdy), (T(1, 3), dRdz),            &  
     &     (T(2, 1), dzdx), (T(2, 2), dzdy), (T(2, 3), dzdz),            &  
     &     (T(3, 1), dphidx), (T(3, 2), dphidy), (T(3, 3), dphidz)
      DATA dRdz / 0._R8/ dzdx, dzdy, dzdz / 0, 0, 1._R8/ dphidz / 0._R8/
      xx = y_rzphi(1) * cos(y_rzphi(3))
      yy = y_rzphi(1) * sin(y_rzphi(3))
      dRdx = xx / y_rzphi(1)
      dRdy = yy / y_rzphi(1)
      Rsq_inv = 1._R8/ (y_rzphi(1) * y_rzphi(1))
      dphidx = - yy * Rsq_inv
      dphidy = xx * Rsq_inv
      do i = 1, 3
         dy_rzphi(i) = 0._R8
         do j = 1, 3
            dy_rzphi(i) = dy_rzphi(i) + T(i, j) * dy_xyz(j)
         enddo
      enddo
      return
      END
!     .         --------------------------------------------------------
      REAL*8   FUNCTION mu_yb(y, y_hat)
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER i
      REAL*8   y(*), y_hat(*), inner_prod, B_vector(3), rr
!     B(1) = B_x, B(2) = B_y, B(3) = B_z
!     construct cos angle between magnetic field (assumed in phi direction)
!     and y_hat at point y (in r, Z, phi, coordinates) and
!     unit vector y_hat (in cartesian x, y, z, coordinates)
      rr = sqrt(y(1) * y(1) + y(2) * y(2))
      B_vector(1) = y(2) / rr
      B_vector(2) = - y(1) / rr
      B_vector(3) = 0._R8
      inner_prod = 0._R8
      do i = 1, 3
         inner_prod = inner_prod + B_vector(i) * y_hat(i)
      enddo
      mu_yb = inner_prod
      return
      END
!     .         --------------------------------------------------------
      SUBROUTINE xVecCopy(n, tovec, frvec)
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER n, i
      REAL*8    tovec(n), frvec(n)
      do i = 1, n
         tovec(i) = frvec(i)
      enddo
      return
      END
!     .         --------------------------------------------------------
!
!     graphs    --------------------------------------------------------
      SUBROUTINE X_plan_view(nr, Radii, npoints,nx,ny, chord4, n1,n34)
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER n34, nx, ny, nr, n1
!     INTEGER*4 npoints(n34, n34)
      INTEGER   npoints(n34, n34)
      REAL*8    R_min, R_max, Radii(nr)
!     REAL*4  chord4(3, n1, n34, n34)
      REAL*8    chord4(3, n1, n34, n34)
      call VecMnMx( Radii, nr, R_min, R_max)
      call xPlanSet (R_min, R_max)
      call PlPlBnds(nr, Radii)
      call PlPlChrds (npoints, nx, ny, chord4, n1, n34)
      return
      END
!     .         --------------------------------------------------------
      SUBROUTINE xPlanSet (R_min, R_max)
!     INTEGER*4 I1, I2
!     DATA      I1, I2 / 1, 2 /
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      REAL*8   R_min, R_max
!     REAL*4 xdum4(2), ydum4(2), twormax
      REAL*8   twormax
 
 
      twormax = R_max + R_max
!NCAR      call agsetf('X/MINIMUM.', - twormax)
!NCAR      call agsetf('X/MAXIMUM.', twormax)
!NCAR      call agsetf('Y/MINIMUM.', - twormax)
!NCAR      call agsetf('Y/MAXIMUM.', twormax)
!NCAR      call agstup(xdum4, 1, 1, 2, 1, ydum4, 1, 1, 2, 1)
!NCAR      call agback
           call EZsets (50,550,  50,550,                                 &  
     &        -twormax, twormax, -twormax, twormax, 1)
      return
      END
!     .         --------------------------------------------------------
      SUBROUTINE PlPlBnds(nr, Radii)
      USE xgraphs
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!     INTEGER*4 nthet4, I1
      INTEGER   nthet4, I1
      INTEGER   i, nr, ir
      REAL*8      Radii(nr)
!     REAL*4  xva4(NWKDIX), yva4(NWKDIX),
!    ^        costhet4(NWKDIX), sinthet4(NWKDIX)
      REAL*8    xva4(NWKDIX), yva4(NWKDIX),                              &  
     &        costhet4(NWKDIX), sinthet4(NWKDIX)
!     REAL*4    dthet4, pi, rim1, thet4
      REAL*8      dthet4, pi, rim1, thet4
!     EQUIVALENCE (xva4(1), wkarx(1, 1)), (yva4(1), wkarx(1, 2)),        &  
!    &     (costhet4(1), wkarx(1, 3)), (sinthet4(1), wkarx(1, 4))
      REAL*8 AREAL
      DATA I1/ 1 /
      nthet4 = NWKDIx
      pi = 4._R8* atan(1._R8)
      dthet4 = (2._R8* pi) / AREAL(nthet4 - 1)
      do i = 1, nthet4
         rim1 = AREAL(i - 1)
         thet4 = dthet4 * rim1
         costhet4(i) = cos(thet4)
         sinthet4(i) = sin(thet4)
         wkarx(i,3) = costhet4(i)
         wkarx(i,4) = sinthet4(i)
      enddo
      do ir = 1, nr
         do i = 1, nthet4
            xva4(i) = Radii(ir) * costhet4(i)
            yva4(i) = Radii(ir) * sinthet4(i)
            wkarx(i,1) = xva4(i)
            wkarx(i,2) = yva4(i)
         enddo
!NCAR         call agcurv(xva4, 1, yva4, 1, nthet4, 1)
              call EZcurv(xva4, yva4, nthet4)
      enddo
      return
      END
!     .         --------------------------------------------------------
      SUBROUTINE PlPlChrds(npoints, nx, ny, chord4, n1, n34)
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER   n1, n34
      INTEGER   npoints(n34, n34), nx, ny, ix, iy, np
      REAL*8      chord4(3, n1, n34, n34)
 
      do ix = 1, nx, 2
        do iy = 1, ny, 2
!ncar     call agcurv(chord4(1, 1, ix, iy), 3, chord4(2, 1, ix, iy),3,
!ncar^    npoints(ix, iy), 1)
!         call EZcurvmd(chord4(1, 1, ix, iy), 3, chord4(2, 1, ix, iy),3,
!    ^    npoints(ix, iy))
 
          np = npoints(ix, iy) - 1
 
          call EZdra(chord4(1, 1, ix, iy), chord4(2, 1, ix, iy), 0 )
          call EZdra(chord4(1, 1, ix, iy), chord4(2, 1, ix, iy), 1 )
          call EZdra(chord4(1,np, ix, iy), chord4(2,np, ix, iy), 1 )
 
        enddo
      enddo
      return
      END
!     .         --------------------------------------------------------
      SUBROUTINE X_elevation_view(nr, Radii, Z_min, Z_max, npoints, nx,  &  
     &     ny, chord4, n1, n34)
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER   nr, n34, nx, ny, n1
      INTEGER*4 npoints(n34, n34)
      REAL*8   R_min, R_max, Radii(nr), Z_min, Z_max
!     REAL*4 chord4(3, n1, n34, n34)
      REAL*8   chord4
      REAL*8    ZERO
      DATA    ZERO /                                                     &  
     &        0.0_R8/
      call VecMnMx ( Radii, nr,  R_min, R_max)
!     call xElevSet (0., R_max, Z_min, Z_max)
      call xElevSet (ZERO, R_max, Z_min, Z_max)
      call PlElBnds (R_min, R_max, Z_min, Z_max)
      call PlElChrds (npoints, nx, ny, chord4, n1, n34)
           call MkGrfLst('Xray plan, elevation chords')
           call EZfini(0,0)
      return
      END
!     .         --------------------------------------------------------
      SUBROUTINE PlElBnds (R_min, R_max, Z_min, Z_max)
!     INTEGER*4 i, I1, I2, I3, I4, I5
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER   i, I1, I2, I3, I4, I5
      REAL*8    R_min, R_max, Z_min, Z_max, ris
!     REAL*4  xva4 (5, 2), yva4 (5, 2)
      REAL*8    xva4 (5, 2), yva4 (5, 2)
      DATA I1, I2, I3, I4, I5 / 1, 2, 3, 4, 5 /
      REAL*8    ZERO
      DATA    ZERO /                                                     &  
     &        0.0_R8/
 
!     i = 1: right half
!     i = 2: left half
      do i = I1, I2
         if(i .eq. I1) then
            ris = 1._R8
         else
            ris = -1._R8
         endif
         xva4(I1, i) = ris * R_min
         yva4(I1, i) = Z_min
         xva4(I2, i) = ris * R_max
         yva4(I2, i) = Z_min
         xva4(I3, i) = ris * R_max
         yva4(I3, i) = Z_max
         xva4(I4, i) = ris * R_min
         yva4(I4, i) = Z_max
         xva4(I5, i) = xva4(I1, i)
         yva4(I5, i) = yva4(I1, i)
 
!ncar      call agcurv  (xva4(I1, i),I1, yva4(I1, i),I1,I5,I1)
           call EZcurvmd(xva4(I1, i),I1, yva4(I1, i),I1,I5)
 
      enddo
!          call EZdra(0.00, Z_min    , 0)
!          call EZdra(0.00, Z_min    , 1)
!          call EZdra(0.00, Z_min/10., 1)
!          call EZdra(0.00, 0.00     , 0)
!          call EZdra(0.00, 0.00     , 1)
!          call EZdra(0.00, Z_max/10., 0)
!          call EZdra(0.00, Z_max/10., 1)
!          call EZdra(0.00, Z_max    , 1)
 
           call EZdra(ZERO, Z_min    , 0)
           call EZdra(ZERO, Z_min    , 1)
           call EZdra(ZERO, Z_min/10._R8, 1)
           call EZdra(ZERO, ZERO     , 0)
           call EZdra(ZERO, ZERO     , 1)
           call EZdra(ZERO, Z_max/10._R8, 0)
           call EZdra(ZERO, Z_max/10._R8, 1)
           call EZdra(ZERO, Z_max    , 1)
 
      return
      END
!     .         --------------------------------------------------------
      SUBROUTINE xElevSet (R_min_pl, R_max, Z_min, Z_max)
!     INTEGER*4 I1, I2
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER   I1, I2
      REAL*8   R_min_pl, R_max, Z_min, Z_max, Z_mn, Z_mx
!     REAL*4 xdum4(2), ydum4(2), twormax, twozmin, twozmax
!     REAL*8   xdum4(2), ydum4(2)
      REAL*8   twormax, twozmin, twozmax
!     REAL*4 R_min_pl4
      REAL*8   R_min_pl4
      DATA I1, I2 / 1,2 /
 
      Z_mn = Z_min - 0.5_R8* (Z_max - Z_min)
      Z_mx = Z_max + 0.5_R8* (Z_max - Z_min)
 
      R_min_pl4 = R_min_pl
 
      twormax = 2._R8* R_max
      twozmin = 2._R8* Z_mn
      twozmax = 2._R8* Z_mx
 
 
!NCAR      call agsetf('X/MINIMUM.', R_min_pl4)
!NCAR      call agsetf('X/MAXIMUM.', twormax)
!NCAR      call agsetf('Y/MINIMUM.', twozmin)
!NCAR      call agsetf('Y/MAXIMUM.', twozmax)
!ncar      call agstup(xdum4,I1,I1,I2,I1, ydum4,I1,I1,I2,I1)
!NCAR      call agback
           call EZsets (700,1000, 150,450,                               &  
     &        R_min_pl4, twormax, twozmin, twozmax, 1)
      return
      END
!     .         --------------------------------------------------------
      SUBROUTINE PlElChrds (npoints, nx, ny, chord4, n1, n34)
      USE xgraphs
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER   n34, nx, ny, ix, iy, n1, ip, np
!     INTEGER*4 npoints(n34, n34), I1
      INTEGER   npoints(n34, n34), I1
!     REAL*4 chord4(3, n1, n34, n34)
      REAL*8   chord4(3, n1, n34, n34)
!     REAL*4 rr(1), zz(1)
      REAL*8   rr(NWKDIx), zz(NWKDIx)
!     EQUIVALENCE (rr(1), wkarx(1, 1)), (zz(1), wkarx(1, 2))
      DATA I1 / 1 /
      do ix = 1, nx, 4
         do iy = 1, ny, 4
            do ip = 1, npoints(ix, iy)
               call xyz2rz4(chord4(1, ip, ix, iy),I1, rr(ip), zz(ip))
               wkarx(ip,1) = rr(ip)
               wkarx(ip,2) = zz(ip)
            enddo
!ncar            call agcurv  (rr,I1, zz,I1, npoints(ix, iy),I1)
!                call EZcurvmd(rr,I1, zz,I1, npoints(ix, iy))
 
            np = npoints(ix, iy) - 1
 
            call EZdra(rr( 1), zz( 1), 0 )
            do ip=1,np,30
              call EZdra(rr(ip), zz(ip), 1 )
            enddo
         enddo
      enddo
      return
      END
!     .         --------------------------------------------------------
      SUBROUTINE xyz2rz4(xyz, n, rr, zz)
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER n, i
      REAL*8    xyz(3, n)
!     REAL*4  rr(n), zz(n)
      REAL*8    rr(n), zz(n)
      do i = 1, n
         rr(i) = sqrt(xyz(1, i) * xyz(1, i) + xyz(2, i) * xyz(2, i))
         zz(i) = xyz(3, i)
      enddo
      return
      END
!     .         --------------------------------------------------------
      SUBROUTINE X_pl_intensity_cntrs(pc, nxdim, nx, ny)
      USE emparams
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!============
! idecl:  explicitize implicit REAL declarations:
!     REAL*8 x_pl_emisivity_cntrs
!============
      CHARACTER*5 whichside
      CHARACTER*30 ChStr
      INTEGER   nxdim, nx, ny, iside, ix1(2), ix2(2), iy1, iy2
!     REAL*4 pc (nxdim, ny)
      REAL*8   pc (nxdim, ny)
!.
      REAL*8    pcmin, pcmax, pcave
      REAL*8    xminy, xmaxy, xdumy, yminy, ymaxy
      REAL*8    xAry(NPIXDIM), yAry(NPIXDIM), clevelin(20)
      REAL*8    xGiv(*)      , yGiv(*)
      REAL*8    AREAL
      INTEGER i,j, kclev1, kclev2
      INTEGER ixmax,ixmin,ixstep
      INTEGER jymax,jymin,jystep
      DATA    ix1(1), ix2(1), ix1(2), ix2(2), iy1, iy2 /                 &  
     &           100,  500,    600,    1000 , 200, 600 /
      whichside = 'LEFT'
      iside = 1
      do 10 i=1,nx
        xAry(i) = AREAL(i)
 10   continue
      do 11 j=1,ny
        yAry(j) = AREAL(j)
 11   continue
      call EZrnd2 (xAry(nx)  , xmaxy, i,j)
      call EZrnd2 (yAry(ny)  , ymaxy, i,j)
      xminy = 0.00_R8
      yminy = 0.00_R8
 
      kclev1 = -5
      kclev2 = 2
 
!
      goto 15
      ENTRY      X_pl_emisivity_cntrs(pc, nxdim, nx, ny,                 &  
     &           xGiv, yGiv)
      whichside = 'RIGHT'
      iside = 2
 
!  round up from max found--> nice value; #major divns; #minor divns
      call EZrnd2 (xGiv(nx)  , xmaxy, i,j)
      call EZrnd2 (xmaxy-xGiv(1), xdumy, i,j)
      call EZrnd2 (yGiv(ny)  , ymaxy, i,j)
      if (2._R8*ymaxy .gt. xdumy) then
        xminy = xmaxy - 2.0_R8*ymaxy
        yminy = - ymaxy
      else
        xminy = xmaxy - xdumy
        ymaxy = xdumy / 2._R8
        yminy =-ymaxy
      endif
 
      kclev1 = -5
      kclev2 = 2
!
!     ...................
!
 15   continue
 
      pcmin = pc(1,1)
      pcmax = pc(1,1)
      pcave = 0.00_R8
      do 20 i=1,nx
      do 19 j=1,ny
        pcmin = min(pcmin, pc(i,j))
        pcmax = max(pcmax, pc(i,j))
        pcave = pcave +      pc(i,j)
 19   continue
 20   continue
      pcave = pcave / AREAL(nx *  ny)
 
!
!     NCAR Code STARTS:
!     Force PLOTCHAR to use characters of the lowest quality.
!NCAR      call pcseti ('QU - QUALITY FLAG',2)
!     initialize the drawing of the contour plot.
!NCAR      call cprect (pc, nxdim, nx, ny, wkarx, NWKDIX * 4,
!NCAR     ^     iwkarx, NWKDIX)
!     Draw the default background.
!NCAR      call cpback (pc, wkarx, iwkarx)
!     Draw contour lines and labels.
!NCAR      call cpcldr (pc, wkarx, iwkarx)
!     Add the informational label and the high/low labels.
!NCAR      call cplbdr (pc, wkarx, iwkarx)
!     Compute and print statistics for the plot, label it, and put a
!     boundary line around the edge of the plotter frame.
!     cALL CAPSAP ('counts', TIME, IAMA, 0)
!     cALL LABTOP ('EXAMPLE 1-1',.017)
!     cALL BNDARY
!     NCAR Code ENDS:
!
!     ------------------------------------------------------------------
!     kclev1: >0 --> clevelin contains levels to be used;
!                    dots for index less than kclev2;
!                    solid for index greater or equal kclev2
!     kclev1: =0 --> first contour at clevelin(1) with
!                    next one up by clevelin(2), and so on and so on
!     kclev1: <0 --> rcontr to choose -kclev1 equally spaced values between
!                    clevelin(1) and clevelin(2)
!     clevelin:      array of contour levels; this is output if kclev1<0
!     kclev2:        separates dots from solid lines
!     ------------------------------------------------------------------
!      previous to Mar 93 usage:
!      kclev1 =-10
!      kclev2 = 2
!      pcmin = pcmin + (pcmax-pcmin) / 10.
!      clevelin(1) = pcmin
!      clevelin(2) = pcmax-pcmin
 
 
      pcmin = pcmax / 100._R8
      clevelin(1) = pcmin
      clevelin(2) = pcmax*0.91_R8
 
      ixmin = 1
      ixmax = nx
      jymin = 1
      jymax = ny
      ixstep= 1
      jystep= 1
 
      if (iside .eq. 1) then
      call EZrcon(ix1(iside),ix2(iside),     iy1,        iy2,            &  
     &                 xminy,     xmaxy,   yminy,      ymaxy,            &  
     &  kclev1,clevelin,kclev2,                                          &  
     &  pc,                                                              &  
     &  nxdim,                                                           &  
     &  xAry, ixmin, ixmax, ixstep,                                      &  
     &  yAry, jymin, jymax, jystep )
      else if (iside .eq. 2) then
      call EZrcon(ix1(iside),ix2(iside),     iy1,        iy2,            &  
     &                 xminy,     xmaxy,   yminy,      ymaxy,            &  
     &  kclev1,clevelin,kclev2,                                          &  
     &  pc,                                                              &  
     &  nxdim,                                                           &  
     &  xGiv, ixmin, ixmax, ixstep,                                      &  
     &  yGiv, jymin, jymax, jystep )
      else
         continue
      endif
!
      if (iside .eq. 1)                                                  &  
     &      write(ChStr,'(''min: '',1pe9.2,'' w/cm^2 ...$'')')pcmin
      if (iside .eq. 2)                                                  &  
     &      write(ChStr,'(''min: '',1pe9.2,'' w/cm^3 ...$'')')pcmin
      call EZwrit(ix1(iside), 750, ChStr, 0,0)
!
      write(ChStr,'(''max: '',1pe9.2,''  $'')')pcmax
      call EZwrit(ix1(iside), 725, ChStr, 0,0)
!
      write(ChStr,'(''ave: '', 1pe9.2,'' $'')')pcave
      call EZwrit(ix1(iside), 700, ChStr, 0,0)
 
      if (iside .eq. 1) then
      write(ChStr,'(''Xmx/mn:'',2(1x,1pe9.2),''$'')')                    &  
     &               xAry(ixmax),xary(1)
      call EZwrit(ix1(iside), 675, ChStr, 0,0)
      write(ChStr,'(''Ymx/mn:'',2(1x,1pe9.2),''$'')')                    &  
     &               yAry(jymax),yary(1)
      call EZwrit(ix1(iside), 650, ChStr, 0,0)
      else if (iside .eq. 2 ) then
      write(ChStr,'(''Xmx/mn:'',2(1x,1pe9.2),''$'')')                    &  
     &               xGiv(ixmax),xGiv(1)
      call EZwrit(ix1(iside), 675, ChStr, 0,0)
      write(ChStr,'(''Ymx/mn:'',2(1x,1pe9.2),''$'')')                    &  
     &               yGiv(jymax),yGiv(1)
      call EZwrit(ix1(iside), 650, ChStr, 0,0)
 
      endif
 
!
      return
      END
!
!     +---------------------------------------------------------------+
!
      SUBROUTINE xSigGen(Emev, Ephoton, mu, sigma)
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      REAL*8   EMev, Ephoton, mu, sigma, T0, k
      REAL*8   signorm, Mevpermcsq
      COMMON / sigbrcom0 / signorm, Mevpermcsq
!
!     converts to / from Mev to Koch, Motz units
!     input energies in Mev, angle in radians
!     return sigma in cm**2 / sR / MeV
!
      T0 = Emev / Mevpermcsq
      k = Ephoton / Mevpermcsq
      call xSigBrem(T0, k, mu, sigma)
      sigma = sigma / Mevpermcsq
      return
      END
 
!     +---------------------------------------------------------------+
 
      SUBROUTINE xBremSetup(Z)
 
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      REAL*8   Z
      REAL*8   pie_local, r0, ECHARGE, EMASS, CLIGHT, HBAR, alpha
      PARAMETER(ECHARGE = 4.8E-10_R8, EMASS = 9.1E-28_R8,                &  
     &     CLIGHT = 3.E10_R8, HBAR = 1.05E-27_R8)
      REAL*8   signorm, Mevpermcsq
      COMMON / sigbrcom0 / signorm, Mevpermcsq
!
      MevPerMcsq = .511_R8
      pie_local = 4._R8* atan(1._R8)
      alpha = ECHARGE * ECHARGE / (HBAR * CLIGHT)
      r0 = ECHARGE * ECHARGE / (EMASS * CLIGHT * CLIGHT)
      signorm = alpha * r0 * r0 * Z * Z / ( 8._R8* pie_local )
      return
      END
 
!     +-----------------------------------------------------------------+
 
      SUBROUTINE xSigBrem(T0, k, mu, sigma)
!     Reference:
!     H. W. Koch and J. W. Motz,
!     ``Bremsstrahlung Cross-Section Formulas and Related Data,''
!     \RMP{31}{920}{59}.
!     Formula 2BN on p 924
!     '2' denotes differential in photon energy and angle
!     'B' denotes Born approximation
!     'N' denotes no screening
!     T0       INPUT  initial k-e of electron in a collision; m_0 c^2 units
!     k        INPUT  energy of emitted photon;               m_0 c^2 units
!     mu       INPUT  cos theta_0; angle between initial e-momentum and emitted
!                     photon
!     sigma    OUTPUT differential cross section
!                     w.r.t. photon energy(k) & solid angle(\Omega_k) in
!                     cm^2 per atom per incident electron
!     signorm  Z^2 r_0^ /(8 \pi 137); units of cm^2
 
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      REAL*8   sigma, k, mu, T0
      REAL*8   ksq, L, Delta0, p0, p0sq, eps, epsQ, Q, p, E, E0, T,      &  
     &     Qsq, beta, beta0, Cth0, Sth0, EE0, term(13), E0sq,            &  
     &     Esq, Delta0sq, Delta0_4, p0sqDsq, p0sqD4, Sth0sq
      REAL*8   signorm, Mevpermcsq
      COMMON / sigbrcom0 / signorm, Mevpermcsq
!
!     computation of much used factors
      ksq = k * k
      p0 = sqrt(T0 * (T0 + 2._R8))
      p0sq = p0 * p0
      E0 = 1._R8+ T0
      E0sq = E0 * E0
      E = E0 - k
      Esq = E * E
      EE0 = E * E0
      T = E - 1._R8
      p =  sqrt(T * (T + 2._R8))
      beta0 = sqrt(1._R8- 1._R8/ (E0 * E0))
      beta = sqrt(1._R8- 1._R8/ (E * E))
      Cth0 = mu
      Qsq = p0sq + ksq - 2._R8* p0 * k * Cth0
      Q = sqrt(Qsq)
      Sth0 = sqrt(1._R8- mu * mu)
      Sth0sq = Sth0 * Sth0
      Delta0 = E0 - p0 * cth0
      Delta0sq = Delta0 * Delta0
      Delta0_4 = Delta0sq * Delta0sq
      p0sqDsq = p0sq * Delta0sq
      p0sqD4 = p0sq * Delta0_4
      eps = log((E + p) / (E - p))
      epsQ = log((Q + p) / (Q - p))
      L = log((EE0 - 1._R8+ p * p0) / (EE0 - 1._R8- p * p0))
      term(1) = 8._R8* Sth0sq * ( 2._R8* E0sq + 1._R8) / p0sqD4
      term(2) = - 2._R8* (5._R8* E0sq + 2._R8* EE0 + 3._R8) / p0sqDsq
      term(3) = - 2._R8* (p0sq - ksq) / (Qsq * Delta0sq)
      term(4) = 4._R8* E / (p0sq * Delta0)
      term(5) = 4._R8* E0 * Sth0sq * (3._R8* k - p0sq * E) / p0sqD4
      term(6) = 4._R8* E0sq * (E0sq + Esq) / p0sqDsq
      term(7) = 2._R8* (1._R8- (7._R8* E0sq - 3._R8* EE0 + Esq))         &  
     &          / p0sqDsq
      term(8) = 2._R8* k * (E0sq + EE0 - 1._R8) / (p0sq * Delta0)
      term(9) = - 4._R8* eps / (p * Delta0)
      term(10) = epsQ / (p * Q)
      term(11) = 4._R8/ Delta0sq - 6._R8* k / Delta0 -                   &  
     &      2._R8* k * (p0sq - ksq) / (Qsq * Delta0)
      term(12) = L / (p * p0)
      term(13) = signorm * p / (k * p0)
      sigma = term(13) * (                                               &  
     &     term(1) + term(2) + term(3) + term(4)                         &  
     &     + term(12) * ( term(5) + term(6) + term(7) + term(8) )        &  
     &     + term(9) + term(10) * term(11)                               &  
     &     )
      return
      END
!
!     graphs          --------------------------------------------------
!     End Xray Camera --------------------------------------------------
!
!
!     . ----------------------------------------------------------------
!     . ----------------------------------------------------------------
!     . Absorbers and Scintilator Treated Here -------------------------
!     . ----------------------------------------------------------------
!     . ----------------------------------------------------------------
!     .
!     Copyright D. W. Ignat, S. von Goeler, J. E. Stevens,
!     P. G. Roney, E. J. Valeo
!     Princeton University, Plasma Physics Laboratory, 1992, 1993
!
!     ------------------------------------------------------------------
!
      SUBROUTINE GetXmnFac(nEs, EbyMeV, XmnFac, dcu, FC )
      USE xparams
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER nEs
      REAL*8    EbyMeV(nEs), XmnFac(nEs)
      REAL*8 DCU
      CHARACTER*2 FC
      CHARACTER*32 ErrorTxt
!                       ! DCU: thickness of foil (cm)
!                       ! FC:  FoilCode: CU,AG,TA,MO, or 00 for no foil
      INTEGER  nds, ndsDIM
      PARAMETER ( NDS = 12, ndsDIM = 12)
      REAL*8    dcm(nds)
      INTEGER i
      INTEGER iVACUUM,  iCARBON,iALUMINM,   iIRON, iNICKEL,              &  
     &                  iCOPPER, iIODINE, iCESIUM, iCSI,                 &  
     &                  iTANTALUM, iMOLYBDENUM, iSILVER
      PARAMETER (iVACUUM = 1, iCARBON = 2, iALUMINM = 3)
      PARAMETER (iIRON = 4, iNICKEL = 5)
      PARAMETER (iCOPPER = 6, iIODINE = 7 , iCESIUM = 8, iCSI = 9)
      PARAMETER (iTANTALUM = 10, iMOLYBDENUM = 11, iSILVER = 12)
      INTEGER nSigInt
      PARAMETER (nSigInt = NENDIM)
      REAL*8    Ekv(nSigInt)
      REAL*8    rhoCsCSI, rhoICsI, rhoCarbon
      DATA    rhoCsCSI, rhoICsI, rhoCarbon /                             &  
     &           2.307_R8,   2.203_R8,       1.0_R8/
!                                            !must be looked up
!     k Edges at  33.169,       35.985,  kV
!     These Data Give the Shield Always in Place in the Actual Experiment:
!     DATA    dcm(iVACUUM), dcm(iCARBON), dcm(iALUMINM) /
!    ^                0.0 ,         0.2 ,     0.37489   /
!     DATA    dcm(iIRON)  , dcm(iNICKEL), dcm(iCOPPER)  /
!    ^            0.0031  ,       0.0137,         0.0   /
!     DATA    dcm(iIODINE), dcm(iCESIUM), dcm(iCSI)     /
!    ^              0.0   ,         0.0 ,     0.0297    /
!     DATA  dcm(iTANTALUM),dcm(iMOLYBDENUM),dcm(iSILVER)/
!    ^              0.0   ,         0.0 ,         0.0   /
!
!     If we wanted to have phony foils to look at Photon Temperature
!     at low energy then SvG recomments the following:
!        10 kv 0.014 cm Al; 20 kv 0.109 cm Al;
!        30 kv 0.335 cm Al; 40 kv 0.650 cm Al.
!
 
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 
!     put in thicknesses of magnetic screen over the xray tube --
!     identical to data statements above
      dcm(iVACUUM)   = 0.0_R8
      dcm(iCARBON)   = 0.2_R8
      dcm(iALUMINM)  = 0.37489_R8
      dcm(iIRON)     = 0.0031_R8
      dcm(iNICKEL)   = 0.0137_R8
      dcm(iMOLYBDENUM)= 0.0006_R8
!
!     put in active material on the screen
      dcm(iCSI)      = 0.0297_R8
!
!     zero out thicknesses for optional materials -
      dcm(iCOPPER)    = 0.0_R8
      dcm(iIODINE)    = 0.0_R8
      dcm(iCESIUM)    = 0.0_R8
      dcm(iTANTALUM)  = 0.0_R8
      dcm(iSILVER)    = 0.0_R8
!     .                                 determine subscript into dcm
!     .                                 corresponding to foil code, and
!     .                                 put foil thickness into place;
!     .
!     .                                 If PhonyFoil, ie  FC='PF', then
!     .                                 zero out the default thickness for
!     .                                 everything, and put in an Aluminum foil
!     .                                 of the given thickness (SvGoeler Dec93)
      if (FC.eq.'CU') then
        dcm(iCOPPER) = dcu
      else if (FC.eq.'TA') then
        dcm(iTANTALUM) = dcu
      else if (FC.eq.'MO') then
        dcm(iMOLYBDENUM) = dcu + dcm(iMOLYBDENUM)
      else if (FC.eq.'AG') then
        dcm(iSILVER) = dcu
!     .                                Phony Foil Here
      else if (FC.eq.'PF') then
        do i=1,NDS
           dcm(i)=0.00_R8
        enddo
!     .                                put back in thickness of the
!     .                                aluminum phony foil 'PF'
           dcm(iALUMINM) = dcu
!     .                                put back in active screen material
           dcm(iCSI)      = 0.0297_R8
!     .
      else if (FC.eq.'00') then
        continue
      else
        write(ErrorTxt,                                                  &  
     &      '(a2,'' !! Bad foil code in GetXmnFac'')') FC
        call LSCwarn (ErrorTxt)
      endif
 
 
!     .                                 Convert units of mc^2? or MeV to keV
      do 10 i = 1,nEs
        Ekv(i)  =  EbyMeV(i)*1000._R8
 10   continue
 
        call GetXmiss (nEs, Ekv, XmnFac,  nds, dcm )
 
!
!DBG        write(4,'(('' Ekv, Xmn: '', 2(1pe10.3,1x),/))')
!DBG     ^  (Ekv(i), XmnFac(i), i=1,nEs )
      return
      END
!
!     ------------------------------------------------------------------
!
      SUBROUTINE Absorb12(Energy, iNAME,                                 &  
     &                   BarnPerA, MuByRho, nSigma, iError, ErrorTxt )
!     Reference:
!     X ray crosss sections based on W. H. McMaster, el al, UCRL 50174
!     Sec. II rev. 1
!     May 1969
!     \bibitem{mcmaster}
!     W. H. McMaster, \etal , UCRL 50174,
!     Sec. II rev. 1,  May 1969.
 
!     Copyright D. W. Ignat, S. von Goeler, J. E. Stevens,
!     P. G. Roney, E. J. Valeo
!     Princeton University, Plasma Physics Laboratory, 1992, 1993
!
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER iNAME, iError
      CHARACTER*32  BLANK
      CHARACTER*(*) ErrorTxt
      REAL*8    Energy, BarnPerA, MuByRho, nSigma
      INTEGER iVACUUM,  iCARBON,iALUMINM,   iIRON, iNICKEL,              &  
     &                  iCOPPER, iIODINE, iCESIUM, iCSI,                 &  
     &                  iTANTALUM, iMOLYBDENUM, iSILVER
      PARAMETER (iVACUUM = 1, iCARBON = 2, iALUMINM = 3)
      PARAMETER (iIRON = 4, iNICKEL = 5)
      PARAMETER (iCOPPER = 6, iIODINE = 7 , iCESIUM = 8, iCSI = 9)
      PARAMETER (iTANTALUM = 10, iMOLYBDENUM = 11, iSILVER = 12)
      INTEGER iLshell, iKshell,  iCoher ,  iInCoh
      INTEGER nLshells, nLatoms
      PARAMETER (iLshell = 1, iKshell = 2,  iCoher = 3, iInCoh = 4)
      PARAMETER (nLshells = 4, nLatoms = 12)
!     nLshells is a misnomer: 1 refers to L interactions below k edge
!     .                       2 refers to K interactions above k edge
!     .                       3 refers to coherent scattering
!     .                       4 refers to incoherent scattering
!
!     If M shell interactions would be considered, then
!     nMshells would be 5: M, L, K, coherent, incoherent
!     Note, we are ignoring M shells in Cs and I
!
!
      REAL*8    kEdge(nLatoms), grPcc(nLatoms), AtmWt(nLatoms)
      REAL*8    Avogadro
 
      INTEGER iShell, i,j
      REAL*8    lnEnergy, lnSigma, epower(0:3)
 
      REAL*8    aL(4,nLshells,nLatoms)
      DATA    BLANK /'                                '/
!.    .              0123456789 123456789 123456789 12
      DATA    Avogadro / 6.02252E-01_R8/
!     DATA    Avogadro / 6.02252e-01 /  =  6.02 10^23 x 10^-24 for barns
      DATA     kEdge(iALUMINM), kEdge(iIRON), kEdge(iCOPPER)/            &  
     &                   1.560_R8,        7.112_R8,          8.979_R8/
      DATA     kEdge(iIODINE ), kEdge(iCESIUM),kEdge(iNICKEL)/           &  
     &                  33.169_R8,       35.985_R8,          8.333_R8/
      DATA kEdge(iTANTALUM),kEdge(iMOLYBDENUM),kEdge(iSILVER)/           &  
     &                  67.414_R8,       19.999_R8,         25.514_R8/
 
      DATA     grPcc(iALUMINM), grPcc(iIRON), grPcc(iCOPPER)/            &  
     &                 2.702_R8,       7.860_R8,         8.940_R8/
      DATA     grPcc(iIODINE ), grPcc(iCESIUM),grPcc(iNICKEL)/           &  
     &                 4.94_R8,      1.873_R8,        8.90_R8/
      DATA grPcc(iTANTALUM),grPcc(iMOLYBDENUM),grPcc(iSILVER) /          &  
     &                16.600_R8,       10.220_R8,          10.500_R8/
 
      DATA     AtmWt(iALUMINM), AtmWt(iIRON), AtmWt(iCOPPER)/            &  
     &                26.970_R8,      55.850_R8,        63.540_R8/
      DATA     AtmWt(iIODINE ), AtmWt(iCESIUM),AtmWt(iNICKEL)/           &  
     &                126.90_R8,      132.910_R8,        58.690_R8/
      DATA AtmWt(iTANTALUM),AtmWt(iMOLYBDENUM),AtmWt(iSILVER)/           &  
     &                180.950_R8,       95.950_R8,       107.880_R8/
 
      DATA     kEdge(iCARBON), grPcc(iCARBON) , AtmWt(iCARBON)  /        &  
     &                  0.00_R8,         1.580_R8,      12.010_R8/
 
      DATA   ((aL(i,j      ,iVACUUM),i = 1,4),j = 1,nLshells)  /         &  
     &        16*0.0_R8/
 
      DATA   ((aL(i,j,iCARBON ),i = 1,4),j = 1,nLshells)          /      &  
     &         0.0_R8,        0.0_R8,        0.0_R8,        0.0_R8,      &  
     &         1.06879E1_R8, -2.71400_R8,   -2.00530E-1_R8,              &  
     &         2.07248E-2_R8,                                            &  
     &         3.10861_R8,   -2.60580E-1_R8,-2.71974E-1_R8,              &  
     &         1.35181E-2_R8,                                            &  
     &        -9.82878E-1_R8, 1.46693_R8,   -2.93743E-1_R8,              &  
     &         1.56005E-2_R8/
 
      DATA   ((aL(i,j,iALUMINM),i = 1,4),j = 1,nLshells)          /      &  
     &         1.08711E1_R8, -2.77860_R8,    1.75853E-1_R8, 0.0_R8,      &  
     &         1.31738E1_R8, -2.18203_R8,   -2.5896E-1_R8,               &  
     &         2.2284E-2_R8,                                             &  
     &         4.51995_R8,  1.40549E-1_R8,-3.5244E-1_R8,  1.9369E-2_R8,  &  
     &        -4.39322E-1_R8, 1.30867_R8,  -2.22648E-1_R8,               &  
     &         7.54210E-3_R8/
 
      DATA   ((aL(i,j,iIRON   ),i = 1,4),j = 1,nLshells)          /      &  
     &         1.36696E1_R8, -2.39195_R8,   -1.37648E-1_R8, 0._R8,       &  
     &         1.43456E1_R8, -1.23491_R8,   -4.18785E-1_R8,              &  
     &         3.21662E-2_R8,                                            &  
     &         5.93292_R8,    2.25038E-1_R8,-3.61748E-1_R8,              &  
     &         1.93024E-2_R8,                                            &  
     &        -3.42379E-1_R8,1.57245_R8,   -2.53198E-1_R8,               &  
     &         9.85822E-3_R8/
 
      DATA   ((aL(i,j,iNICKEL ),i = 1,4),j = 1,nLshells)          /      &  
     &         1.39848E1_R8, -2.48080_R8,   -8.88115E-2_R8, 0.0_R8,      &  
     &         1.42388E1_R8, -9.67736E-1_R8,-4.78070E-1_R8,              &  
     &         3.66138E-2_R8,                                            &  
     &         6.09204_R8, 2.52277E-1_R8, -3.66568E-1_R8,                &  
     &         1.96586E-2_R8,                                            &  
     &        -5.0436E-1_R8, 1.70040_R8,    -2.76443E-1_R8,              &  
     &         1.12628E-2_R8/
 
      DATA   ((aL(i,j,iCOPPER ),i = 1,4),j = 1,nLshells)           /     &  
     &        14.2439_R8,    -2.58677_R8,   -6.67398E-2_R8, 0._R8,       &  
     &        14.5808_R8,    -1.18375_R8,   -0.41385_R8,                 &  
     &         3.12088E-2_R8,                                            &  
     &         6.17739_R8,    2.73123E-1_R8,-3.7236E-1_R8,               &  
     &         2.01638E-2_R8,                                            &  
     &        -5.7021E-1_R8,  1.75042_R8,   -2.84555E-1_R8,              &  
     &         1.1693E-2_R8/
 
      DATA   ((aL(i,j,iIODINE ),i = 1,4),j = 1,nLshells)           /     &  
     &         1.64086E1_R8, -2.48214_R8,   -5.07179E-2_R8, 0.0_R8,      &  
     &         1.21075E1_R8,  1.43635_R8,   -8.82038E-1_R8,              &  
     &         6.03575E-2_R8,                                            &  
     &         7.27415_R8,    3.77223E-1_R8,-3.69728E-1_R8,              &  
     &         1.86280E-2_R8,                                            &  
     &        -4.04420E-2_R8, 1.65596_R8,   -2.510676_R8,                &  
     &         9.04874E-3_R8/
 
      DATA   ((aL(i,j,iCESIUM ),i = 1,4),j = 1,nLshells)           /     &  
     &         1.65418E1_R8, -2.46363_R8,   -5.42849E-2_R8, 0.0_R8,      &  
     &         1.13757E1_R8,  1.94161_R8,   -9.83232E-1_R8,              &  
     &         6.71986E-2_R8,                                            &  
     &         7.33490_R8,    3.76825E-1_R8,-3.65713E-1_R8,              &  
     &         1.81843E-2_R8,                                            &  
     &         1.84861E-1_R8, 1.50030_R8,   -2.13333E-1_R8,              &  
     &         6.24264E-3_R8/
 
      DATA  ((aL(i,j,iTANTALUM),i = 1,4),j = 1,nLshells)           /     &  
     &         1.72410E1_R8, -2.30313_R8,   -5.91006E-2_R8, 0.0_R8,      &  
     &         8.65271_R8,    3.73117_R8,   -1.26359_R8,                 &  
     &         8.23539E-2_R8,                                            &  
     &         7.94534_R8,    3.87299E-1_R8,-3.47926E-1_R8,              &  
     &         1.63299E-2_R8,                                            &  
     &         1.96871E-1_R8, 1.50623_R8,   -1.91396E-1_R8,              &  
     &         3.70889E-3_R8/
 
      DATA  ((aL(i,j,iSILVER),i = 1,4),j = 1,nLshells)             /     &  
     &         1.56869E+1_R8,-2.22636_R8,   -1.12223E-1_R8, 0.0_R8,      &  
     &         1.33926E+1_R8, 4.41380E-1_R8,-6.93711E-1_R8,              &  
     &         4.82085E-2_R8,                                            &  
     &         7.06446_R8,    3.63456E-1_R8,-3.73597E-1_R8,              &  
     &         1.92478E-2_R8,                                            &  
     &        -1.66475E-1_R8, 1.65794_R8,   -2.48740E-1_R8,              &  
     &         8.66218E-3_R8/
 
      DATA  ((aL(i,j,iMOLYBDENUM),i = 1,4),j = 1,nLshells)         /     &  
     &         1.53494E+1_R8,-2.26646_R8,   -1.16881E-1_R8, 0.0_R8,      &  
     &         1.39853E+1_R8,-1.17426E-1_R8,-5.91094E-1_R8,              &  
     &         4.17843E-2_R8,                                            &  
     &         6.84600_R8,    3.02797E-1_R8,-3.51131E-1_R8,              &  
     &         1.74403E-2_R8,                                            &  
     &        -5.62860E-2_R8, 1.55778_R8,   -2.33341E-1_R8,              &  
     &         7.85506E-3_R8/
!
!     ............................................. Execution begins here
!
      iError  =  0
      ErrorTxt= BLANK
 
      if (iNAME .eq. 1) then
!     .                                 The practice of including 'data' for
!     .                                 vacuum came from Jim Stevens, and we
!     .                                 are not sure why.  Here it is bypassed
        BarnPerA  =  0._R8
        MuByRho   =  0._R8
        nSigma    =  0._R8
        return
      endif
 
      if (Energy .lt. 1._R8) then
        write(ErrorTxt,'('' invalid low energy: '',1pe10.3)')  Energy
        BarnPerA  =  0._R8
        MuByRho   =  0._R8
        nSigma    =  0._R8
        iError  =  1
        return
      endif
 
      if (Energy .gt. 1000._R8) then
        write(ErrorTxt,'('' invalid high energy: '',1pe10.3)') Energy
        BarnPerA  =  0._R8
        MuByRho   =  0._R8
        nSigma    =  0._R8
        iError  =  2
        return
      endif
 
      if (iNAME .gt. 12) then
        write(ErrorTxt,'('' invalid atom: '', i4)')           iNAME
        BarnPerA  =  0._R8
        MuByRho   =  0._R8
        nSigma    =  0._R8
        iError  =  3
        return
      endif
 
 
!     !                                 Major branch based on number of shells
!     !                                 (branch not used; M shell ignored now)
 
      IF (iNAME .le. nLatoms) then
!     .                                 Use this code for L-shell atoms
        if (Energy .lt. kEdge(iNAME)) then
          iShell  =  1
        else
          iShell  =  2
        endif
!
        lnEnergy  =  log(Energy)
        epower(0)  =  1._R8
        do 5 i = 1,3
          epower(i) = epower(i-1)*lnEnergy
 5      continue
 
          lnSigma  =  0._R8
        do 10 i = 1,4
          lnSigma  =  lnSigma + aL(i,iShell,iNAME)*epower(i-1)
  10    continue
 
        BarnPerA  =  exp(lnSigma)
 
        do 40 iShell = nLshells-1, nLshells
            lnSigma  =  0._R8
          do 30 i = 1,4
            lnSigma  =  lnSigma + aL(i,iShell,iNAME)*epower(i-1)
  30      continue
        BarnPerA  =  BarnPerA + exp(lnSigma)
  40    continue
 
!     .                                 Use this code for M-shell atoms
!     .                                 if and when this is added
      ELSE
        write(ErrorTxt,'('' invalid atom: '', i4)') iNAME
        BarnPerA  =  0._R8
        MuByRho   =  0._R8
        nSigma    =  0._R8
        iError  =  4
        return
      ENDIF
!     .                                 Finally, get N rho / AtWt * Sigma
 
 
 
      MuByRho   =  BarnPerA * Avogadro / AtmWt(iNAME)
      nSigma    =  MuByRho  * grPcc(iNAME)
 
 
      return
      END
!
!     ------------------------------------------------------------------
!
      SUBROUTINE GetXmiss (nSigInt, Eaxis, Xmiss,                        &  
     &                     nds, dcm )
!     Copyright D. W. Ignat, S. von Goeler, J. E. Stevens,
!     P. G. Roney, E. J. Valeo
!     Princeton University, Plasma Physics Laboratory, 1992, 1993
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER nSigInt
      INTEGER nds
      CHARACTER*32 ErrorTxt
      REAL*8    dcm(nds)
      REAL*8    Eaxis(nSigInt), Xmiss(nSigInt)
      INTEGER iVACUUM,  iCARBON,iALUMINM,   iIRON, iNICKEL,              &  
     &                  iCOPPER, iIODINE, iCESIUM, iCSI,                 &  
     &                  iTANTALUM, iMOLYBDENUM, iSILVER
      PARAMETER (iVACUUM = 1, iCARBON = 2, iALUMINM = 3)
      PARAMETER (iIRON = 4, iNICKEL = 5)
      PARAMETER (iCOPPER = 6, iIODINE = 7 , iCESIUM = 8, iCSI = 9)
      PARAMETER (iTANTALUM = 10, iMOLYBDENUM = 11, iSILVER = 12)
 
      INTEGER i, j, iError
      REAL*8    BarnPerA, MuByRho, nSigma
      REAL*8    MuRhoD, MuRhoCsI
      REAL*8    rhoCsCSI, rhoICsI, rhoCarbon
      DATA    rhoCsCSI, rhoICsI, rhoCarbon /                             &  
     &           2.307_R8,   2.203_R8,       1.0_R8/
!                                            !must be looked up
!     k Edges at  33.169,       35.985,  kV
!
      do 70 i = 1,nSigInt
!
      MuRhoD    =  0.00_R8
      MuRhoCsI  =  0.00_R8
 
      if (Eaxis(i) .gt. 1000._R8) then
        Xmiss(i) = 1.00_R8
      else if(Eaxis(i) .lt. 1.00_R8) then
        Xmiss(i) = 0.00_R8
      else
!
        do 50 j = 1,nds
          if (dcm(j) .gt. 0._R8.and. j .ne. iCsI) then
            call Absorb12(Eaxis(i), j,                                   &  
     &                   BarnPerA, MuByRho, nSigma, iError, ErrorTxt)
!
            if (j .eq. iCARBON) then
              MuRhoD  =  MuRhoD + MuByRho*rhoCarbon*dcm(j)
            else
              MuRhoD  =  MuRhoD + nSigma*dcm(j)
            endif
            if (iError .ge. 1) call LSCwarn(ErrorTxt)
          endif
!
          if (dcm(j) .gt. 0._R8.and. j .eq. iCSI) then
            call Absorb12(Eaxis(i), iIODINE,                             &  
     &                   BarnPerA, MuByRho, nSigma, iError, ErrorTxt)
!
            MuRhoCsI  =  MuRhoCsI + MuByRho*rhoICsI*dcm(j)
!
            call Absorb12(Eaxis(i), iCESIUM,                             &  
     &                   BarnPerA, MuByRho, nSigma, iError, ErrorTxt)
!
            MuRhoCsI  =  MuRhoCsI + MuByRho*rhoCsCsI*dcm(j)
!
            if (iError .ge. 1) call LSCwarn(ErrorTxt)
          endif
!
 50     continue
!
!
        if (MuRhoCsI .eq. 0._R8) then
          call LSCwarn(' no phosphor screen found; 1e-5 assumed')
          MuRhoCsI  =  1.E-05_R8
        endif
 
        Xmiss(i)  =  exp ( - MuRhoD )*(1._R8- exp(-MuRhoCsI) )
      endif
 70   continue
!
      return
      END
 
!     ------------------------------------------------------------------
!     . ----------------------------------------------------------------
!     . Absorbers and Scintilators End Here ----------------------------
!     . ----------------------------------------------------------------
 
 
 
!                                                                      |
!     -----------------------------------------------------------------|
!                                                                      |
      SUBROUTINE XRopen(GivnUnit, FileFrstNa)
!============
! idecl:  explicitize implicit REAL declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!     REAL*8 xrclos
!============
      INTEGER GivnUnit, SomeUnit, linepointer
      CHARACTER*1 PERIOD, SPACE
      CHARACTER*8  FileFrstNa, FileScndNa
      CHARACTER*8  mmddyy, hhmmss
      CHARACTER*17 FileCompNa
      CHARACTER*75 SomeString
 
      COMMON / XRrwcom / linepointer,SomeUnit
      DATA        PERIOD, SPACE / '.' , ' ' /
 
      CHARACTER*8  zzdate
      CHARACTER*10 zztime
      CHARACTER*5  zzzone
!
!     zzdate  ccyymmdd    cc=century yy=year mm=month dd=day
!     zztime  hhmmss.ttt  hh=hour mm=minute ss=second ttt=millisecond
!     zzzone  Shhmm       S=sign (-=west, +=east) hh=hour mm=minute
!
      INTEGER ivals(8)
!                             EXAMPLE
!     ivals ( 1 ) year        "2000"
!     ivals ( 2 ) month       "10"
!     ivals ( 3 ) day         "12"
!     ivals ( 4 ) t - UTC     "-0400"
!     ivals ( 5 ) hour        "8"
!     ivals ( 6 ) minute      "21"
!     ivals ( 7 ) second      "23"
!     ivals ( 8 ) millisecond "621"
!
      call date_and_time(zzdate,zztime,zzzone,ivals)
      write(mmddyy,'(a)') zzdate(5:8)//zzdate(3:4)
      write(hhmmss,'(a)') zztime(1:6)
 
      SomeUnit = GivnUnit
 
      if (FileFrstNa .eq. 'XRinput') then
        open (SomeUnit,file='XRinput',                                   &  
     &        status='unknown')
        linepointer=0
      else
        FileScndNa = hhmmss
        open (SomeUnit,file=FileFrstNa//PERIOD//FileScndNa,              &  
     &        status='unknown',err=1300 )
        linepointer=0
        SomeString=FileFrstNa//SPACE//mmddyy//SPACE//hhmmss
        write(SomeUnit,'(a75)') SomeString
        linepointer=1
      endif
 1300 continue
!1300 continue need something real here
      return
 
      ENTRY      XRclos(GivnUnit)
      close(GivnUnit)
      return
      END
 
      SUBROUTINE XRdscw(xarray, yarray, npoints, SomeString)
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER SomeUnit, linepointer
      INTEGER i,npoints
      REAL*8    xarray(npoints), yarray(npoints)
      CHARACTER*(*) SomeString
      COMMON / XRrwcom / linepointer,SomeUnit
      write(SomeUnit,'(a75)') SomeString
        linepointer=linepointer+1
      do 10 i=1,npoints
        write(SomeUnit,'(1x,1pe10.3,1x,1pe10.3)')xarray(i),yarray(i)
        linepointer=linepointer+1
 10   continue
      return
      END
!
!      SUBROUTINE XRdscr(xarray, yarray, npoints, SomeString)
!      IMPLICIT NONE
!      INTEGER SomeUnit, linepointer
!      INTEGER i,npoints
!      REAL*8    xarray(npoints), yarray(npoints)
!      CHARACTER*75 SomeString
!      COMMON / XRrwcom / linepointer,SomeUnit
!
!      if(linepointer .eq. 0) then
!        read(SomeUnit,'(a75)') SomeString
!        linepointer=linepointer+1
!      endif
!
!      read(SomeUnit,'(a75)') SomeString
!        linepointer=linepointer+1
!      do 10 i=1,npoints
!        read (SomeUnit,'(1x,1pe10.3,1x,1pe10.3)')xarray(i),yarray(i)
!        linepointer=linepointer+1
! 10   continue
!
!      return
!      END
 
#else
      SUBROUTINE XrayCam2
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      return
      END
      SUBROUTINE XRopen(i,a)
!============
! idecl:  explicitize implicit REAL declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!============
      INTEGER i
      CHARACTER*24 a
      return
      ENTRY XRclos(i)
      return
      END
      SUBROUTINE XRdscw(x, y, n, a)
      USE camera
      USE emitter
      USE emparams
      USE FeBins
      USE numerics
      USE params
      USE ProfBody
      USE RayBins
      USE RayWrk
      USE TSCgrap
      USE tscunits
      USE xgraphs
      USE xparams
      USE Xray
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER n
      REAL*8    x(n), y(n)
      CHARACTER*(*) a
      return
      END
!      SUBROUTINE XRdscr
!      return
!      END
#endif
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
! 24May2005 fgtok

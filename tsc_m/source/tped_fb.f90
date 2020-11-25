      subroutine tped_fb(a3208,tped_corr)
!
      USE CLINAM
      USE SAPROP
!     USE COMWOY
!     USE SVDCOM
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER l,ll,llsav,li,ii,ngr,iabs,ngrvwl,ngrvcl,ig,iigr,igr,i
      INTEGER icnt8,icnt9,iw,iswtch,n,nl,ie,np,iem10,iem20,iii
      INTEGER index, lsav
!============
!============      
!
      REAL*8, DIMENSION(20) :: eped_cur,eped_den,eped_tem
      REAL*8 ne_ped,te_ped,te_ped1,te_ped2,t_eped,cur_t0,tped_corr,AREAL
      REAL*8 eped_den_min,eped_den_max,eped_cur_min,eped_cur_max,tint,a3208
!     REAL*8 coiturn,denom,denom1,term1,term1d,denom2,term2
!     REAL*8 term2d,denom3,term3,term3d,sumturn,sumcur,fbchix
      REAL*8 rlin,rdis,pval,rval,rpval,zden,zr,za,zip,zk,zd
      REAL*8 zp,zb,zq,ztohmic,ztau,pped,zlh,znp,zbe,zli
!     REAL*8 dtramp,aplfb,fluxdif,psinn,pinterp,psinp,ellip,delta
!     REAL*8 x1,x2,z1,z2,psi1,psi2,psi3,psi4,psisn1,alphac,term
!     REAL*8 fluxdif1,formf,fac,facrxw,fluxoff,termd,atuv,ppcur
!     REAL*8 aplav,plcur,rebmax,rebmin

      REAL*8, DIMENSION(5) :: te_dum
      INTEGER j,j1,j2,npsihm,jd1,jd2,jc1,jc2,flag_tped


!-------------------------------------------
!	interpolate 
! needed: total current at time t
!         pedestal density
!         a3003
!	  pedestal location (for pedestal density)
     
      npsihm = int(pwidthc * AREAL(npsit))
      ne_ped = ane(npsihm)
      cur_t0 = apld0*udsi
      te_ped = 0.5*(te(npsihm)+ti(npsihm))
      zip = tcurdtp*tpi*udsi*1.E-6_R8
      zr = rmajor
      za = rminor
      zk = shape5
      zd = shape3
      zb = gzero/zr
      zip = tcurdtp*tpi*udsi*1.E-6_R8
      flag_tped = int(a3208)
      SELECT CASE (flag_tped)
      CASE (1)
!     P. Snyder EPED1 model
      call eped1_fb(t_eped)      
      CASE (2)
!     Sugihara scaling. Ref: Plasma Phys. COntrol. Fusion 45 (2003) L55-L62
        znp = 1.E-20_R8*ane(npsihm)
        rlin = 0
        rdis = 0._R8
        do i=iminn,imaxx
          if(iexv(i,nh).eq.1 .or. iexs(i,nh).eq.1) cycle
          pval = psi(i,nh)
          if(pval.gt.psilim) cycle
          call reval(pval,idens,isurf,rval,rpval,i,nh)
          rlin = rlin + rval*deex*udsd
          rdis = rdis + deex
        enddo
        if(rdis.le.deex) rdis = deex
        zden = (rlin/rdis)/udsd
        zp = (pohmic+paux+prad+palpha)*1.E-6_R8
!       zlh= 1.38*zden**0.77 * zb**0.92 * zr**1.23 * za**0.76 
        zlh= 2.84/2.5*zden**0.58 * zb**0.82 * zr * za**0.81 

        pped = 2.41E+03 * 2.5**0.33 * znp**(-0.33) *                    &
     &                 zr**1.33 * za**(-4) *                            &
     &                 zip**2 * (zr/za)**(-2.94) *                      &
     &                 zk**3.62 * (0.5*(1.0+zk**2))**(-2.33) *          &
     &                 (1.0+zd)**3.2 *                                  &
     &                 (zp/zlh)**0.06 
        t_eped = 0.5*pped/1.6e-19/ne_ped
        IF (zip .lt. 7.0 .AND. t_eped .gt. 4.0E+03) t_eped = 4.0E+03
      CASE (3)
!     Maget scaling. Ref: Nucl. Fusion 53 (2013) 093011
        zbe = betapol*2.0*zb/zip
        pped = 117.1*zip**1.59*zb**0.66*(2*ali2)**0.57*zbe**0.02*acoef(3209)**1.06
        t_eped = 0.5*pped/1.6e-19/ne_ped
        IF (zip .lt. 7.0 .AND. t_eped .gt. 4.0E+03) t_eped = 4.0E+03
      CASE (10:)
        t_eped = a3208
      END SELECT     
      tped_corr = te_ped/t_eped

      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

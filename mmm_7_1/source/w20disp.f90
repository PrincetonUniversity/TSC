subroutine w20disp( zamr, zami, zbmr, zbmi, wz, G_ave, kps )
  use w20data
  implicit none

  !Equation matrices
  real*8, dimension(neq,neq) :: zamr, zami, zbmr, zbmi

  !Input frequency
  complex*16, intent(in) :: wz

  REAL*8 RAV, G_ave
  REAL*8 H1,XH,R,HQR,HQI,WM

  COMPLEX*16 ALPC,ALPHA,HQ,IU,H2
  COMPLEX*16 A_lpc, A_lpk

  REAL*8 wr,wi,H

  real*8 :: alp, alf, kps
!  REAL*8 alp,kps,alf

  real*8 zdenom

  integer :: j1, j2

!
!    variables i=1,6: e Phi/Te, dTi/Ti, dni/ni, dTe/Te, dnq/nq, dTq/Tq
!    variables j=1,6 same as for i
!

  IU=(0.,1.)
!
!..print header
!
  if ( lprint .gt. 2 ) then
     write (6,*)
     write (6,*) &
          'Weiland-Nordman eigenvalue equations, subroutine etawn6'
     write (6,*) '(all frequencies normalized by omega_{De})'
     write (6,*) '(all diffusivities normalized by ' &
          ,'omega_{De} / k_y^2'
  endif


!     zflz   = aimp * zflh / zimp**2
!
!..diagnostic output
!

  if ( lprint .gt. 6 ) then
     write (6,*)
     write (6,*) '--------------------------------------'
     write (6,*)
     write (6,*)
     write (6,*) gne,' = gne'
     write (6,*) gni,' = gni'
     write (6,*) gnz,' = gnz'
     write (6,*) gte,' = gte'
     write (6,*) gti,' = gti'
     write (6,*) gtz,' = gtz'
     write (6,*)
     write (6,*) vef,' = vef'
     write (6,*) betae, '=betae'
     write (6,*) zimp,' = zimp'
     write (6,*) aimp,' = aimp'
     write (6,*) 
     write (6,*) zflh,' = zflh'
     write (6,*) zflz,' = zflz'
     write (6,*)
     write (6,*) bt,' =bt'
     write (6,*) bt1,'  =bt1'
     write (6,*) kappa,'  =kappa'
     write (6,*) wz,'   =wz'
  endif

  !Reset system of equations to zero
  zamr(1:neq,1:neq) = 0.0
  zami(1:neq,1:neq) = 0.0
  zbmr(1:neq,1:neq) = 0.0
  zbmi(1:neq,1:neq) = 0.0

!
!..Nine  equations with impurities, trapped electrons, parallel ion
!  motion, collisions,  FLR , finite beta and parallel motion of impurities
!
!  Equations for e phi/Te, Ti, ni, Tet, n_z, T_z, F, Av, K
!  Here F = omega*gamma (collisions), K = omega*Av and Av is the parallel
!   magnetic vector potential.
!

      H1=4.*q*q*zflh

      H=0.5*ABS(shat)/q
      H2=IU*H

      A_lpk=0.5*shat*SQRT(H1*XT*zflh*(1.0+TAUH*(gni+gti)/(2*WZ)) )
!      A_lpk=0.5*shat*SQRT(H1*XT*zflh*(1.0+TAUH*(gni+gti))/(2*WZ) )

      IF ( DREAL(A_lpk) .LT. 0) then
         A_lpk=-A_lpk
      END IF

      if ( abs(dreal(A_lpk)).lt.0.01 ) then
         A_lpk = A_lpk - dreal(A_lpk) + 0.01
      end if

      ALPC=-IU*A_lpk
      ALPHA=-IU*ABS(shat)*q*zflh
      XH=ABS(ALPHA/ALPC)
      ALPC=XH*ALPC
      R=2.*ABS(DREAL(WZ*ALPC))
      IF(R.LT.0.001) R=0.001    !! NEW 01.03.8
      HQ=2.*ALPC/H1

      !Geometric average multiplying \omega_{De}
      G_ave=(1.0+0.5*shear/R)*EXP(-0.25/R)-0.5*alpha_MHD*(1.-EXP(-1/R))
      G_ave=max(G_ave,1.0E-2)
!      G_ave = max(G_ave,1.0E-4)
! Temporary investigation
!      G_ave = 1.0

      alp=max(0.1,0.5*R)

      zdenom = (2.D0*zflh*q*q*betae*(1.D0 - fte))
      alf = (alp+1.0E-6) / (zdenom+1.0E-6)
      kps = 0.5*SQRT(alp/zflh)/q
      kpc = 1.D0
      RAV=1.+0.25*shat**2.0/alp

!
!  *********
      HQR=DREAL(HQ)
      HQI=DIMAG(HQ)
!  *********


      if ( lprint > 6 ) then
         WRITE(6,10002) alp,shat**2,RAV
10002    FORMAT(2X,'alp=',ES15.6,'  SH2=',ES15.6,' RAV=',ES15.6)
         WRITE(6,10003) XH,G_ave,alf
10003    FORMAT(2X,'XH=',ES15.6,' G_ave=',ES15.6,' alf=',ES15.6)
         write(6,10005) alpha_MHD
10005    FORMAT(2X,'MHD Alpha=',ES15.6)
         write(6,10006) DREAL(A_lpk),R
10006    FORMAT(2X,' A_lpk=',ES15.6,' R='ES15.6)
         write(6,*)"HQR=",hqr,"   HQI=",hqi
      end if

!
!--- WE NOW DEFINE THE MATRIX ELEMENTS --------------
! hydrogen density
!
!
!...Equations rewritten in terms of the profile gradients
!
        zamr(1,1) = - G_ave+HQR + 0.5 * ( gni - zflh * ztauh * ( gni + gti ) )
        zami(1,1) = HQI
        zamr(1,2) = (HQR-G_ave)*ztauh
        zami(1,2) = ztauh*HQI
        zamr(1,3) = (HQR-G_ave)*ztauh
        zami(1,3) = ztauh*HQI
        zamr(1,8) = -em*ztauh*HQR*0.5*(gni+gti)/kpc
        zami(1,8) = -em*ztauh*HQI*0.5*(gni+gti)/kpc
        zamr(1,9) = -em*HQR/kpc
        zami(1,9) = -em*HQI/kpc
!
        zbmr(1,1) = zflh
        zbmr(1,3) = 1.
!
!  hydrogen energy
!
        zamr(2,1) = 0.5*(gti - tvr * gni)
        zamr(2,2) = - ztauh * ftr
!
        zbmr(2,2) = 1.
        zbmr(2,3) = - tvr
!
!  total electron density expressed in ion density and imp density
!
        zamr(3,1) = -1.0 + 0.5*fte*gne
       zami(3,1) = vef*(1.-fte)
       zamr(3,3) = 1. - znz_ne -zns_ne
       zami(3,3) = -vef*(1. - znz_ne - zns_ne)
       zamr(3,4) = fte
       zamr(3,5) = znz_ne
       zami(3,5) = -vef*znz_ne
       zami(3,7) = vef*fte
       zamr(3,8) = -em*0.5*gne*(1. - fte)/kpc
       zami(3,8) =  em*0.5*gne*(1. - fte)*vef/kpc
       zamr(3,9) =  em* (1. - fte)* (1.+0.5*gne) / kpc
       zami(3,9) = -em* (1. - fte)* vef / kpc

      zbmr(3,1) = fte - 1.
      zbmr(3,3) = 1. - znz_ne - zns_ne
      zbmr(3,5) = znz_ne
      zbmr(3,9) = em*(1. - fte)/kpc
!
!  trapped electron energy
!
      zamr(4,1) = fte * 0.5 * ( gte - tvr*gne )
      zami(4,1) = vef*tvr*(bt-2.5*(1.-fte))
      zami(4,3) = -vef*tvr*bt1*(1.-znz_ne -zns_ne)
      zamr(4,4) = fte*ftr
      zami(4,5) = -vef*tvr*bt1*znz_ne
      zami(4,7) = -ftr*vef*fte
!
      zbmr(4,1) = (1. - fte)*tvr
      zbmr(4,3) = -(1. - znz_ne -zns_ne)*tvr
      zbmr(4,4) = fte
      zbmr(4,5) = -znz_ne*tvr
!
!  impurity density
!
      zamr(5,1) = - G_ave + zimp*HQR/aimp        &
           +0.5*( gnz - zflz*ztauz*(gnz+gtz) )
      zami(5,1) = zimp*HQI/aimp

      zamr(5,5) = (HQR*zimp/aimp-G_ave)*ztauz
      zami(5,5) = zimp*ztauz*HQI/aimp

      zamr(5,6) = (HQR*zimp/aimp-G_ave)*ztauz
      zami(5,6) = zimp*ztauz*HQI/aimp

      zamr(5,8) = -em*HQR*zimp*ztauz*0.5*(gnz+gtz)/(kpc*aimp)
      zami(5,8) = -em*HQI*zimp*ztauz*0.5*(gnz+gtz)/(kpc*aimp)

      zamr(5,9) = -em*HQR*zimp/(kpc*aimp)
      zami(5,9) = -em*HQI*zimp/(kpc*aimp)
 
      zbmr(5,1) = zflz
      zbmr(5,5) = 1.
!
!  impurity energy
!
      zamr(6,1) = 0.5*(gtz - tvr * gnz)
      zamr(6,6) = -ztauz*ftr
!
      zbmr(6,5) = -tvr
      zbmr(6,6) = 1.
!
!  variable F
!
      zamr(7,1) = 0.5*gte - 1.0
      zami(7,1) = vef
      ZAMR(7,7) = 1.
      zami(7,7) = -vef
!
      zbmr(7,1) = -1.
      zbmr(7,7) = 1.
!
!
!  electromagnetic parallel vectorpotential Av = e A_par/Te
!
      zamr(8,1) = em1*kpc*(0.5*gne+HQR*(fft+zimp*fzft/aimp))
      zami(8,1) = em1*HQI*(fft+zimp*fzft/aimp)*kpc

      zamr(8,2) = em1*HQR*ztauh*fft*kpc
      zami(8,2) = em1*HQI*ztauh*fft*kpc

      zamr(8,3) = em1*HQR*ztauh*fft*kpc
      zami(8,3) = em1*HQI*ztauh*fft*kpc

      zamr(8,5) = em1*HQR*zimp*ztauz*fzft*kpc/aimp
      zami(8,5) = em1*HQI*zimp*ztauz*fzft*kpc/aimp

      zamr(8,6) = em1*HQR*zimp*ztauz*fzft*kpc/aimp
      zami(8,6) = em1*HQI*zimp*ztauz*fzft*kpc/aimp

      zamr(8,8) = em1*(0.5*(gne + gte) - alf*zflh*RAV ) &
           - em1*HQR*(fft*ztauh*0.5*(gni+gti)           &
           + zimp*fzft*ztauz*0.5*(gtz+gnz)/aimp)


      zami(8,8) = -em1*HQI*(fft*ztauh*0.5*(gni + gti)   &
           +zimp*fzft*ztauz*0.5*(gnz + gtz)/aimp)


      zamr(8,9) = -em1*(0.5*gne+HQR*(fft+zimp*fzft/aimp))
      zami(8,9) = -em1*HQI*(fft+zimp*fzft/aimp)

      zbmr(8,1) = em1*kpc
      zbmr(8,8) = em1
      zbmr(8,9)= -em1
!
!     K = omega*Av
!
      zamr(9,9) = em1
!
      zbmr(9,8) = em1
!
      if ( lprint .gt. 6 ) then
        write (6,*)
        write (6,*) ' zamr(j1,j2)  j2 ->'
        do j1=1,9
          write (6,192) (zamr(j1,j2),j2=1,9)
        enddo

        write (6,*)

        write (6,*) ' zami(j1,j2)  j2 ->'
        do j1=1,9
          write(6,192) (zami(j1,j2),j2=1,9)
        enddo

       write (6,*)
        write (6,*) ' zbmr(j1,j2)  j2->'
        do j1=1,9
          write (6,192) (zbmr(j1,j2),j2=1,9)
        enddo

       write (6,*)
       write(6,*) ' zbmi(j1,j2)  j2->'
       do j1=1,9
         write(6,192) (zbmi(j1,j2),j2=1,9)
       enddo

 192  format (1p10e12.4)
      endif

    End subroutine w20disp

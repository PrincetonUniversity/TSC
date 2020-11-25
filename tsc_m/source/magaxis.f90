      subroutine magaxis
!              from version 9.06
!              mod 28 Oct 93 to patch around bug in newton iter.
!..   find magnetic axis location by newtons iteration
!
      USE CLINAM
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!     -------------------------------------------------------------
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER ndima,kgone,kcycgone,isw,jfault,ntime,minc,mloop
      INTEGER imin2,imax2,jmin2,jmax2,i,j,imagt,jmagt,iguess,jguess
      INTEGER ii,iw,jw,n
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 pfx,pfy,pfxy,pds,pmin,xmaglst,zmaglst,pave,zdpsipm
      REAL*8 zdpsi0m,pmneg,xguess,zguess,xtry,ztry,gradsq,dpsidx
      REAL*8 dpsidz,gsval,psval,psixz,psixx,psizz,denom,delx,delz
      REAL*8 xmag1,zmag1,zdxm,zdzm,zdpsimin
!============
      parameter  (ndima=6)
      dimension  pfx(ndima),pfy(ndima),pfxy(ndima,ndima),pds(ndima)
!     -------------------------------------------------------------
      data  kgone, kcycgone /0, 9999999/
!============      
      isw = 1
      jfault = 0
!
!.....first look for min psi in plasma domain
      pmin = 1.E20_R8
      ntime = 0
      xmaglst = xmag
      zmaglst = zmag
!
     if(acoef(501).eq.0 .and. acoef(504).ne.0) then
        xmag = xmaglst
        zmag = zmaglst
        imag = (xmag-ccon)/deex + 2
        jmag = nh
        psimin = -1.e8
        return
     endif
!
      pmin = 1.E20_R8
      minc = 2
      if(kcycle.le.0) minc = 4
      mloop = 0
  121 continue
!
      imin2 = imag - minc
      imax2 = imag + minc
      jmin2 = jmag - minc
      jmax2 = jmag + minc
      if(imin2 .lt. iminn+2 .or. igone.eq.1) imin2 = iminn+2
      if(imax2 .gt. imaxx-2 .or. igone.eq.1) imax2 = imaxx-2
      if(jmin2 .lt. jminn   .or. igone.eq.1) jmin2 = jminn
      if(jmax2 .gt. jmaxx   .or. igone.eq.1) jmax2 = jmaxx
!
!
      if(isym .eq. 1 .and. zplas .eq. 0._R8.and. kcycle.le.0             &  
     &      .and. irst1.ne.2) jmax2 = 2
!
      do 120 i=imin2,imax2
      do 130 j=jmin2,jmax2
      if(iexv(i,j).eq.1 .or. iexs(i,j).eq.1) go to 130
      if(psi(i,j).gt.pmin) go to 130
      imagt = i
      jmagt = j
      pmin = psi(i,j)
  130 continue
  120 continue
      if(igone .eq. 1 .or. lrswtch.gt.0) then
           psimin = pmin
           imag = imagt
           jmag = jmagt
           xmag = xary(imagt)
           zmag = zary(jmagt)
           return
      endif
!
      pave = .25_R8*(psi(imagt+1,jmagt) + psi(imagt-1,jmagt)             &  
     &          + psi(imagt,jmagt+1) + psi(imagt,jmagt-1))
      if(pmin .le. pave) go to 122
      mloop = mloop + 1
      if(mloop.gt.4) go to 112
      minc = minc*2
      go to 121
!
!
  122 continue
      imag = imagt
      jmag = jmagt
      psimin = pmin
      xmag = xary(imag)
      zmag = zary(jmag)
!
!
!
!.....CHECK FOR MINIMUM SLOPE TO NEIGHBORING POINTS:
         iguess= imag
         jguess= jmag
      if (isym.eq. 0)  go to 200
! ....for symmetric cases:
         zdpsipm= psi(iguess+1,jmag)-psi(iguess-1,jmag)
         zdpsi0m= psi(iguess  ,jmag)-psi(iguess-1,jmag)
      if (zdpsipm.gt. 0.0_R8)  imag= iguess-1
      if (abs(zdpsipm).lt. 0.25_R8*abs(zdpsi0m))  imag= iguess
!---- if (psi(imag+1,jmag).gt.psi(imag-1,jmag))  imag= imag-1
!
!<777 format (1x,4i4,f10.5,2x,6f10.4,2x,2f10.5)
      go to 205
! ....for non-symmetric cases test with corner-points:
  200 continue
!--      iguess= imag
!--      jguess= jmag
         pmneg = min(psi(imag-1,jmag-1),                                 &  
     &                 psi(imag+1,jmag-1),                               &  
     &                 psi(imag-1,jmag+1),                               &  
     &                 psi(imag+1,jmag+1))
      if (psi(iguess+1,jguess+1).gt.pmneg)  go to 201
         imag= iguess
         jmag= jguess
      go to 205
  201 if (psi(iguess-1,jguess+1).gt.pmneg)  go to 202
         imag= iguess-1
         jmag= jguess
      go to 205
  202 if (psi(iguess+1,jguess-1).gt.pmneg)  go to 203
         imag= iguess
         jmag= jguess-1
      go to 205
  203 if (psi(iguess-1,jguess-1).gt.pmneg)  go to 205
         imag= iguess-1
         jmag= jguess-1
!
!.....check for axis out of range
  205 continue
      if(imag.eq.imin2) go to 112
      if(imag.eq.imax2) go to 112
      if(isym.eq.0 .and. jmag.eq.jmax2) go to 112
      if(isym.eq.0 .and. jmag.eq.jmin2) go to 112
!
      if(acoef(895) .eq. 0) then
      do 135 ii=1,nwire
      iw = iwire(ii)
      jw = jwire(ii)
      if(iw.eq.imag .and. jw.eq.jmag) go to 112
  135 continue
      endif
!
!.....begin newtons iteration for mag axis location
      xguess = xary(imag)
      zguess = zary(jmag)
      xtry = xguess
      ztry = zguess
      call grap(3,ztry,xtry,gradsq,dpsidx,dpsidz,gsval,psval,            &  
     &          psixz,psixx,psizz,isw)
      isw = 0
      do 100 n=1,8
      nn = n
      denom = psixx*psizz-psixz**2
      if(denom.eq.0) go to 111
      delx = (-psizz*dpsidx + psixz*dpsidz)/denom
      delz = ( psixz*dpsidx - psixx*dpsidz)/denom
      xtry = xtry + delx
      ztry = ztry + delz
      if(abs(xtry-xguess).gt. 2._R8*deex) go to 111
      if(abs(ztry-zguess).gt. 2._R8*deez) go to 111
      if(abs(delx)+abs(delz) .lt. 1.E-9_R8* deex) go to 101
      call grap(3,ztry,xtry,gradsq,dpsidx,dpsidz,gsval,psval,            &  
     &          psixz,psixx,psizz,isw)
!     if(gradsq .le. 1.e-20) go to 101
  100 continue
      go to 111
  101 continue
      denom = psixx*psizz-psixz**2
      if(denom.lt.0) go to 110
!
!.....found minimum:
      xmag = xtry
      zmag = ztry
      psimin = psval
      if(acoef(97).gt.0.and.ntime.eq.0.and.psimin.gt.psilim) go to 123
      return
!
  123 ntime = 1
      if(kcycle.gt.0 .and. igone.eq.0) then
      igone = 1
      write(nout,1124) kcycle, jfault, imag,jmag
      write(nterm,1124) kcycle, jfault, imag,jmag
 1124 format("  * * * igone = 1 at cycle =",i7," * * *     jfault =",i4,  &  
     &                                                                   &  
     &                  "   imag,jmag =", 2i4)
                  endif
      xmag = xary(imag)
      zmag = zary(jmag)
      if (igone.eq.0)   then
      write( nout,1126) xmaglst,zmaglst,xmag,zmag, psimin, phalo,        &  
     & psilim
      write(nterm,1126) xmaglst,zmaglst,xmag,zmag, psimin, phalo,        &  
     & psilim
            endif
!
      return
!
!.....x-point
  110    jfault = 110
      if(acoef(97).gt.0 .and. ntime.eq.0) go to 123
      if(lrswtch.gt.0 .or. acoef(97).gt.0) return
      write(nout,1212) xtry,ztry,psval,psixz,psixx,psizz
 1212 format(' xpoint,xz,ps,pxz,pxx,pzz=',1p6e11.3)
!
         xmag = xary(iguess)
         zmag = zary(jguess)
         psimin = psi(iguess,jguess)
      return
!
!.....error exit
  112    jfault = 112
      if(acoef(97).gt.0 .and. ntime.eq.0) go to 123
      if(lrswtch.gt.0 .or. acoef(97).gt.0) return
      write(nout,1112) imag,jmag
 1112 format(" error, mag axis out of range or at coil location",        &  
     &       "   imag,jmag =",2i4)
      ineg=18
!
!.....warning exit
  111    jfault = 111
      kgone = kgone + 1
!ccccccc    if(acoef(97).gt.0 .and. ntime.eq.0) go to 123
!ccccccc    if(lrswtch.gt.0 .or. acoef(97).gt.0) return
      if(acoef(97).gt.0 .and. kcycle.eq.kcycgone) go to 123
      if(lrswtch.gt.0) return
      if (kgone.eq.1 .or. mod(kgone,50).eq.0)   then
        if (kgone.eq.1)                                                  &  
     &        write(nout,1213) nn,xtry,xguess,ztry,zguess,denom,delx,    &  
     & delz
 1213 format(" warning error in magaxis ,n=",i2,1p7e12.4)
      write(nout,1125) kgone, kcycle, jfault, imag,jmag
        if (mod(kgone,200).eq.0)                                         &  
     &      write(nterm,1125) kgone, kcycle, jfault, imag,jmag
 1125 format("  * * * kgone =",i4,"  at cycle =",i7," * * *    jfault    &  
     & ="                                                                &  
     &,i4,"   imag,jmag =", 2i4)
                        endif
!
         xmag = xary(iguess)
         zmag = zary(jguess)
         psimin = psi(iguess,jguess)
      call axm2d2 (psi,penx,iguess,jguess, xmag , zmag ,xmag,zmag,       &  
     &                 xmag1,zmag1, pds,                                 &  
     &                 pfx,pfy,pfxy,ndima,deex,deez,isym)
         zdxm= xmag1-xmag
         zdzm= zmag1-zmag
         zdpsimin= psimin-pds(1)
         psimin = pds(1)
         xmag   = xmag1
         zmag   = zmag1
!-----
      if (kgone.ne.1 .and. mod(kgone,50).ne.0)   return
      write( nout,1126) xmaglst,zmaglst,xmag,zmag, psimin, phalo,        &  
     & psilim
      if (mod(kgone,200).eq.0)                                           &  
     &write(nterm,1126) xmaglst,zmaglst,xmag,zmag, psimin, phalo,        &  
     & psilim
 1126 format("  *****  xmaglst  zmaglst     xmag     zmag      psimin    &  
     &   phalo     psilim"/ 7x, 4f9.5, 1x, 1p3e11.4)
!----
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

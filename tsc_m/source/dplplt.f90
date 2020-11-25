      subroutine dplplt
!.....6.904.1 dplplt
!
!.....make a plot of divertor plate region
!
      USE CLINAM
      USE S50COM
      USE SVDCOM
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER nplot,inplt,k1,k2,imin,imax,jmin,jmax,i,j,n,l
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 cval,xmin,xmax,zmin,zmax,dx50,dz50,dum1,dum2
      REAL*8 dum3,dum4,dum5,dum6,dum7,hfluxm,dmax,dis,fnorm
      REAL*8 thseg,xm,zm,al,xp,zp,xco,zco,rc,zc,x5,z5,x6,z6,x3,z3
      REAL*8 x4,z4,psp,rplusa,radius
      REAL*8 pinterp
!============
      dimension cval(7), psp(6)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: cmat,bmat
!============      
      INTEGER :: istat = 0 
!============      
!
      if(idiv.eq.0.or.noplot(88).gt.0._R8) return
      if(iplate.eq.0 .or. nplate.eq.0) return
!
      if(.NOT.ALLOCATED(cmat)) ALLOCATE (cmat(pnx,pnx), STAT=istat)
      if(.NOT.ALLOCATED(bmat)) ALLOCATE (bmat(pnx,pnx), STAT=istat) 
      if ( istat .ne. 0 ) stop "Allocation Error : dplplt"

      nplot=acoef(710)
      do 666 inplt=1,nplot
      xmin=acoef(711+(inplt-1)*4)
      xmax=acoef(712+(inplt-1)*4)
      zmin=acoef(713+(inplt-1)*4)
      zmax=acoef(714+(inplt-1)*4)
      if(xmin.ge.xmax .or. zmin.ge.zmax) go to 666
!
      call maps(xmin,xmax,zmin,zmax,.142_R8,.858_R8,.274_R8,1.0_R8)
!
      k1 = 7
      k2 = 7
      rplusa = rmajor + rminor
      do 901 i=1,6
      radius = rplusa + float(i)*0.005
      psp(i) = pinterp(radius,zmag,imag,jmag)
 901  continue
      cval(7) = psisep
      cval(6) = psp(1)
      cval(5) = psp(2)
      cval(4) = psp(3)
      cval(3) = psp(4)
      cval(2) = psp(5)
      cval(1) = psp(6)
!     cval(3) = psisep
!     cval(2) = psisep + pscrape
!     cval(1)=psisep + pscrape*.02_R8/.006_R8
!
      imin = (xmin-ccon)/deex + 3
      imax = (xmax-ccon)/deex + 2
      jmin = (zmin-zzero)/deez + 3
      jmax = (zmax-zzero)/deez + 2
!
      dx50 = (xary(imax)-xary(imin))/(pnx-1)
      dz50 = (zary(jmax)-zary(jmin))/(pnx-1)
      do 161 i=1,pnx
      r1(i) = xary(imin) + (i-1)*dx50
  161 z1(i) = zary(jmin) + (i-1)*dz50
      do 171 i=1,pnx
      do 171 j=1,pnx
      call grap(1,z1(j),r1(i),dum1,dum2,dum3,dum4,cmat(i,j),             &  
     &                      dum5,dum6,dum7,0)
  171 continue
        if (imovie.eq.7)   call  colora("cyan")
      call rcontr(k1,cval,k2,cmat,pnx,r1,1,pnx,1,                        &  
     &                                z1,1,pnx,1)
!
        if (imovie.eq.7)   call  colora("white")
      do 201 n=1,nplate
      call setcrt(xsega(n,1),zsega(n,1))
      do 8111 l=2,nseg(n)+1
 8111 call vector(xsega(n,l),zsega(n,l))
!
!.....plot distribution of heat flux
        if (imovie.eq.7)   call  colora("red")
      hfluxm = 0._R8
      do 160 l=1,nseg(n)
      hfluxm = max(hfluxm,hplate(n,l))
  160 continue
      if(hfluxm .le. 0) go to 199
      dmax = 0.5_R8*max(dsep(n,1),dsep(n,2))
      if(dmax.gt.0) go to 198
      dmax = 1.E6_R8
      do 197 l=1,nseg(n)+1
      dis = sqrt(zsega(n,l)**2 + (xplas-xsega(n,l))**2)
      dmax = min(dmax,.2_R8*dis)
  197 continue
  198 continue
      if(dmax.le.0) go to 199
      fnorm = dmax/hfluxm
!
      call setcrt(xsega(n,1),zsega(n,1))
      do 170 l=1,nseg(n)
      thseg = atan2(zsega(n,l+1)-zsega(n,l),xsega(n,l+1)-xsega(n,l))
      xm = 0.5_R8*(xsega(n,l+1)+xsega(n,l))
      zm = 0.5_R8*(zsega(n,l+1)+zsega(n,l))
      al = hplate(n,l)*fnorm
      xp = xm + al*sin(thseg)
      zp = zm - al*cos(thseg)
      call vector(xp,zp)
      call vector(xm,zm)
      call setcrt(xp,zp)
  170 continue
  199 continue
        xco = xmin + 0.10_R8*(xmax-xmin)
        zco = zmax - 0.05_R8*(zmax-zmin)
        call setold(xco,zco,1,0,1,0)
        write (s100,1000)   times
      call gtext(s100,80,0)
 1000   format (" time =", f8.3)
      xco = xmin - 0.15_R8*(xmax-xmin)
      zco = zmin
      call setold(xco,zco,1,0,1,1)
      write(s100,1001) n,hfluxm
      call gtext(s100,80,0)
 1001 format(" heat flux on plate",i3,"  max=",f6.2," Mwatt/m**2")
      write(nsc1,1002) n,kcycle
 1002 format(" heat flux, plate",i3," cycle",i7)
!
  201 continue
      call frscj(6)
666   continue
!
        if (imovie.eq.7)   then
!        DEALLOCATE (cmat, bmat)
         return
        end if 
      if(isvd.le.0) then
!        DEALLOCATE (cmat, bmat)
         return
       end if
      if(rcocom(1).le.0) then
!        DEALLOCATE (cmat, bmat)
         return
      end if
!
!.....compare fit with true value
      call maps(xmin,xmax,zmin,zmax,.142_R8,.858_R8,.274_R8,1.0_R8)
      k1 = 1
      k2 = 1
      cval(1) = psisep
      call rcontr(k1,cval,k2,cmat,pnx,r1,1,pnx,1,z1,1,pnx,1)
      do 8112 n=1,nplate
      call setcrt(xsega(n,1),zsega(n,1))
      do 8112 l=2,nseg(n)+1
 8112 call vector(xsega(n,l),zsega(n,l))
      do 172 i=1,pnx
      rc = r1(i)
      do 172 j=1,pnx
      zc = z1(j)
      cmat(i,j) =                                                        &  
!.....flux basis vectors at point i
     & + coef(1)                                                         &  
     & + coef(2)*((rc-rs) + 0.5_R8*rsi*(zc-zs)**2*(1._R8-rsi*(rc-rs)))   &  
     & + coef(3)*((zc-zs))                                               &  
     & + coef(4)*((rc-rs)**2 - (zc-zs)**2 + rsi*(rc-rs)*(zc-zs)**2)      &  
     & + coef(5)*((rc-rs)*(zc-zs)*(1._R8+0.5_R8*rsi*(rc-rs)))            &  
     & + coef(6)*((zc-zs)*((zc-zs)**2 - 3._R8*(rc-rs)**2))               &  
     & + coef(7)*((rc-rs)*((rc-rs)**2 - 3._R8*(zc-zs)**2))
  172 continue
      call setold(rs,zs,1,0,1,1)
      write(s100,2172)
      call gtext(s100,80,0)
 2172 format(1hs)
      k1 = 1
      k2 = 1
      cval(1) = psepcal
      call rcontr(k1,cval,k2,cmat,pnx,r1,1,pnx,1,z1,1,pnx,1)
      call setold(xco,zco,1,0,1,1)
      write(s100,1003) kcycle
      call gtext(s100,80,0)
 1003 format(" actual and calc flux,cycle",i7)
      write(s100,1004) psisep,psepcal
      call gtext(s100,80,0)
 1004 format(" psisep= ",1pe12.4," psepcal= ",1pe12.4)
      write(nsc1,1003) kcycle
      do 82 n=1,nptsp
      if(rcocom(n).lt.xmin .or. rcocom(n).gt.xmax) go to 82
      if(zcocom(n).lt.zmin .or. zcocom(n).gt.zmax) go to 82
      x5 = rcocom(n) + .5_R8*deex
      z5 = zcocom(n) + .5_R8*deez
      x6 = rcocom(n) - .5_R8*deex
      z6 = zcocom(n) - .5_R8*deez
      x3 = rcocom(n) + .5_R8*deex
      z3 = zcocom(n) - .5_R8*deez
      x4 = rcocom(n) - .5_R8*deex
      z4 = zcocom(n) + .5_R8*deez
      call setcrt(x5,z5)
      call vector(x6,z6)
      call setcrt(x3,z3)
      call vector(x4,z4)
   82 continue
      call frscj(6)
!     DEALLOCATE (cmat,bmat)
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

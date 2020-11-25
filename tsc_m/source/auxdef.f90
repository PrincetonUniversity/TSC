      subroutine auxdef
!.....3.71 divplat
!
!......define auxialliary arrays
!
      USE CLINAM
      USE SCR1
      USE SCR15
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER i,j,iajdef,ii,igsdef,nzp1,nxp1
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 rojmin,exp2,aj5,rmid,psimid,rval,rrval,prmin,qmin,rmin
      REAL*8 pval,ppval,xe,xa,xc,v2t,v4t,v6t,v8t,v5t,dels
      REAL*8 ppi,ppj,aje,pvac,qmid
      REAL*8 exp
!============
      if(lrswtch.ge.1) go to 339
      rojmin = fracn0*r0
      exp2 = 3._R8/5._R8
      exp = 5._R8/3._R8
!
!
!.....only evaluate rsurf and qsurf every 10 cycles
      if(kcycle.le.0) ifunct = 10
      if(isurf.eq.0 .and. (idens .eq. 0 .or. idens .eq. 3)) go to 16
      ifunct = ifunct + 1
      if(ifunct.lt.10) go to 16
      ifunct = 0
!
!.....redefine density arrays for isurf.ne. 0 or idens=1,2
      do 15 i=3,nxp
      aj5 = dxdz*xarh(i)
      do 15 j=3,nzp
      rmid = rojmin*ajey(i)
      if(igone.eq.1) go to 14
      if(iexvc(i,j).gt.0 .or. iexs(i,j).eq.1) go to 14
      psimid = .25_R8*(psi(i,j)+psi(i-1,j)+psi(i,j-1)+psi(i-1,j-1))
      if(psimid.gt.psilim) go to 14
      call reval(psimid,idens,isurf,rval,rrval,i,j)
      rmid = rval*ajey(i)
   14 rsurf(i,j) = rmid
      if(kcycle.gt.0) go to 15
      r(i,j) = rmid
      ro(i,j) = rmid
   15 continue
   16 continue
!
!.....define density and pressure roj and pr from r and q arrays
      prmin = tevv*rojmin*udsd/(.5_R8*udsh)
      do 1 i=3,nxp
      qmin = prmin**exp2*ajey(i)
      do 1 j=3,nzp
      rmin = rojmin*ajey(i)
      if(r(i,j).lt.rmin) r(i,j) = rmin
      roj(i,j) = ro(i,j)/ajey(i)
      pro(i,j) = pr(i,j)
      if(igone.eq.1) go to 2
      psimid = .25_R8*(psi(i,j)+psi(i-1,j)+psi(i,j-1)+psi(i-1,j-1))
      if(ipres.eq.1 .or. isurf.eq.1) go to 3
      if(psimid.gt.psilim) go to 2
      if(q(i,j).lt.qmin) go to 2
      if(iexvc(i,j).gt.0) go to 2
!
      pr(i,j) = (q(i,j)/ajey(i))**exp
!
      go to 1
    3 call peval(psimid,1+isurf,pval,ppval,i,j)
      pr(i,j) = pval
      go to 1
    2 q(i,j) = qmin
      pr(i,j) = prmin
    1 continue
!
!
  339 continue
!
!.....   ---> ajdef entry for defining ajphi array only
      iajdef = 0
      go to 430
      entry ajdef
      iajdef = 1
  430 continue
!
!.....define toroidal current array
!
      do 300 i=3,nx
      xe = xary(i)
      xa = xarh(i+1)
      xc = xarh(i)
      v2t = dzisq
      v4t = (xe/xc)*dxisq
      v6t = (xe/xa)*dxisq
      v8t = dzisq
      v5t = -(v2t+v4t+v6t+v8t)
      do 400 j=3-isym,nz
      dels =                  v2t*psi(i,j-1)                             &  
     &     + v4t*psi(i-1,j  )+v5t*psi(i,j  )+v6t*psi(i+1,j  )            &  
     &                       +v8t*psi(i,j+1)
      ajphi(i,j) = dels/xary(i)
  400 continue
  300 continue
!cj dir$ ivdep
       do 311 i=3,nx
      ajphi(i,nzp) = ajphi(i,nz)
  311 continue
      if(isym.ne.0) go to 313
!cj dir$ ivdep
      do 312 i=3,nx
  312 ajphi(i,2)=ajphi(i,3)
  313 continue
!cj dir$ ivdep
      do 321 j=2,nzp
      ajphi(2,j) = ajphi(3,j)
      ajphi(nxp,j) = ajphi(nx,j)
  321 continue
!
!
      if(nwire.le.0) go to 404
      do 405 ii=1,nwire
      i = iwire(ii)
      j = jwire(ii)
!
!.....store current in ccoil and set ajphi=0 for force calculation
      ccoil(ncoil-nwire+ii) = ajphi(i,j)*dxdz
      ajphi(i,j) = 0._R8
  405 continue
  404 continue
      if(iajdef.eq.1) return
!
!.....   ---> gsdef entry for defining gs array only
      igsdef = 0
      go to 440
      entry gsdef
      igsdef = 1
  440 continue
!
!.....define gs array at vertices
!
      do 450 i=3,nx
      do 450 j=3-isym,nz
      gs(i,j) = .25_R8*(g(i,j)*xsqoj(i)                                  &  
     &             + g(i,j+1)*xsqoj(i)                                   &  
     &             + g(i+1,j)*xsqoj(i+1)                                 &  
     &             + g(i+1,j+1)*xsqoj(i+1))
  450 continue
      do 470 j=2,nzp
      gs(nxp,j) = gzero
  470 gs(2,j) = gzero
      do 475 i=2,nxp
      gs(i,1) = gzero
  475 gs(i,nzp) = gzero
      if(isym.eq.0) go to 476
!cj dir$ ivdep
      do 477 i=2,nxp
!..boundary condition for antisymmetric toroidal field
      gs(i,1) = gs(i,3)
      if(jsym.eq.-1) gs(i,1) = -gs(i,3)
477   continue
      go to 478
  476 do 479 i=2,nxp
      gs(i,2) = gzero
!..boundary condition for antisymmetric toroidal field
      if(jsym.eq.-1) gs(i,2) = 0._R8
479   continue
  478 continue
      if(igsdef .eq.1) return
      if(isurf.ne.0 .or. ipres.eq.1) go to 754
!
!.....define bmagy array
!
      nzp1=nzp+1
      do 751 i=3,nx
      do 750 j=2,nzp1
      ppi = .5_R8*(psi(i,j)+psi(i,j-1)-psi(i-1,j)-psi(i-1,j-1))
      ppj = .5_R8*(psi(i,j)+psi(i-1,j)-psi(i,j-1)-psi(i-1,j-1))
      xe = xarh(i)
      aje = xe*dxdz
      bmagy(i,j) = ((deex*ppj )**2                                       &  
     &           +  ( deez*ppi)**2)/aje**2
  750 continue
  751 continue
!cj dir$ ivdep
      do 752 j=2,nzp1
      bmagy(2,j) = bmagy(3,j)
  752 bmagy(nx+1,j) = bmagy(nx,j)
      nxp1 = nxp+1
!cj dir$ ivdep
      do 753 i=2,nxp1
      bmagy(i,2) = bmagy(i,3)
      bmagy(i,nz+1) = bmagy(i,nz)
  753 continue
      return
  754 continue
!
      if(ifunct.ne.0) return
      entry prdef
      pvac = (tevv*udsd*rojmin/(.5_R8*udsh))
      do 5 i=3,nxp
      aj5 = dxdz*xarh(i)
      do 5 j=3,nzp
      qmid = aj5*pvac**(3._R8/5._R8)
      if(igone.eq.1) go to 4
      psimid = .25_R8*(psi(i,j)+psi(i-1,j)+psi(i,j-1)+psi(i-1,j-1))
      if(iexvc(i,j).gt.0 .or. iexs(i,j).eq.1) go to 4
      if(psimid.gt.psilim) go to 4
      call peval(psimid,1+isurf,pval,ppval,i,j)
      qmid = max(pval,pvac)**(3._R8/5._R8)*aj5
    4 qsurf(i,j) = qmid
    5 continue
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

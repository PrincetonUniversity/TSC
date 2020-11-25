      subroutine tstep
!......3.20 tstep
!
!.....determine time step for next cycle
!
      USE CLINAM
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER i,j
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 asize,t2min,vsqm,bsqm,btmax,aje,ai,ppi,aj,ppj,vsq,bsq
      REAL*8 btry,tstry,allam,anum,denom,ratio,factor,ffactry
      REAL*8 dtgrow,dttry,vmax,bmax
!============
      if(ineg.ne.0) return
!.....save previous time step in dtold
!
      dtold = dt
      if(ifrst(1).eq.1 .and. irst1.ne.1)go to 300
      itimets = itimets + 1
      if(itimets .lt. 10) return
      itimets = 0
      dt1 = dtfac*ndiv*min(deex,deez)**2/(etav*8._R8*(1._R8-th))
      dt = (min(dt1,1.2_R8*dtold))
      dt2 = 3._R8*dt1
      if(lrswtch.gt.0) go to 103
!
!.....search for largest allowable time step
!.....limit time step increase to 20 percent per 10 cycles
!
      snum   = 0._R8
      asize  = sqrt(alx*alz)
      t2min  = dtmax**2*3._R8/dtfac**2
      idtmin = 0._R8
      jdtmin = 0._R8
      vsqm   = 0._R8
      bsqm   = 0._R8
      btmax = 0._R8
      do 100 i=3,nx
      aje  =   xary(i)*dxdz
      do 200 j=3-isym,nz
      ai   = (abig(i+1,j)-abig(i-1,j))*.5_R8
      ppi  = ( psi(i+1,j)- psi(i-1,j))*.5_R8
      aj   = (abig(i,j+1)-abig(i,j-1))*.5_R8
      ppj  = ( psi(i,j+1)- psi(i,j-1))*.5_R8
      vsq = ((deex*aj )**2                                               &  
     &    +  ( deez*ai)**2)/aje**2
      bsq = ((deex*ppj )**2                                              &  
     &    +  ( deez*ppi)**2)/aje**2 + 1.E-8_R8
      btry = g(i,j)*xsqoj(i)/xarh(i)
!
      bsqsv(i,j) = bsq
      if(iexvc(i,j).gt.0) go to 200
      if(btry.gt.btmax) btmax = btry
!
!.....look for maximum for plasma points only
      if(vsq.gt.vsqm) vsqm = vsq
      if(bsq.gt.bsqm) bsqm = bsq
      tstry = dxdz/(vsq+2._R8*bsq)
      if(tstry.gt.t2min) go to 200
      t2min = tstry
      idtmin = i
      jdtmin = j
  200 continue
  100 continue
      dt2 = dtfac*sqrt(t2min)
  103 continue
!
      if(btmax.le.0) btmax = 1._R8
      dt3 = dtfac*ndiv*min(deex,deez)/(2.8_R8*btmax)
      if(itevv.eq.0) go to 101
!
!.....for itevv=1 choose tevv so as not to limit time step
      allam = 24._R8-log(1.E-3_R8*sqrt(udsd)/tevv)
      anum = 0.5_R8*1.03E-4_R8*allam*usdr*8._R8*(1._R8-th)*1.5_R8*       &  
     & min(dt2,dt3)
      denom = dtfac*ndiv*min(deex,deez)**2
      if(anum.gt.0 .and. denom.gt.0) tevv=(anum/denom)**(2._R8/3._R8)
      dt1 = 1.5_R8*min(dt2,dt3)
!
  101 continue
      if(amach.le.0 .or. kcycle.le.0) go to 102
      if(iffac.eq.0) go to 102
!
!.....for iffac=1 choose ffac to keep amach.le.acoef(801)
      ratio = acoef(801)/amach - 1._R8
      factor = ratio*ratio*ratio + 1._R8
      if(factor.lt.acoef(802)) factor = acoef(802)
      if(factor.gt.acoef(803)) factor = acoef(803)
      ffactry = factor*ffac
      if(ffactry.gt.acoef(804)) ffactry = acoef(804)
      if(ffactry.lt.acoef(805)) ffactry = acoef(805)
      ffac = ffactry
  102 continue
      dtgrow = 1.2_R8*dt
      dttry = sqrt(1._R8/(1._R8/dt1**2 + 1._R8/dt2**2 + 1._R8/dt3**2))
      dt = min(dttry,dtgrow,dtmax)
      if(lrswtch .ne. 0) return
      vmax = 0._R8
      if(vsqm.gt.0) vmax = sqrt(vsqm)
      bmax = sqrt(bsqm)
      amach = vmax/bmax
!
!.....define reciprocal density so that alfven wave
!.....in vacuum will not limit time step
      do 400 i=3,nx
      do 350 j=3-isym,nz
      rdens(i,j) = 1.0_R8
      if(iexvc(i,j).gt.0) go to 250
      go to 350
  250 if(bsqsv(i,j) .gt. bsqm) rdens(i,j) = bsqm/bsqsv(i,j)
  350 continue
      if(isym.eq.1) rdens(i,1) = rdens(i,3)
  400 continue
!.....now smooth, first store density in bsqsv,
!     then average over nearest neighbors
      do 410 i=3,nx
      do 410 j=3-2*isym,nz
  410 bsqsv(i,j) = rdens(i,j)
      do 420 i=4,nx-1
      do 430 j=4-2*isym,nz-1
  430 rdens(i,j) = .5_R8*bsqsv(i,j) + .125_R8*(bsqsv(i-1,j) + bsqsv(i+1,  &  
     & j)                                                                &  
     &                                 + bsqsv(i,j-1) + bsqsv(i,j+1) )
      if(isym.eq.1) rdens(i,1) = rdens(i,3)
  420 continue
!
!
      return
!
  300 continue
      ifrst(1)=0
      dt = 2._R8*dtmin
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

      subroutine groupcur(grsum,grsum0,gvsum,gvsum0,gcur0ka,gcurfka)
!
!.....calculates total group current and preprogrammed group current
!
      USE CLINAM
      USE SCRATCH
      USE SCR11, ONLY : gcurr
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER iig,ig,l,jg,iw1,iabs,ic1,ic
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 sum0,sumv,sumt,sumve,fac,sumvi
      REAL*8 sum
      REAL*8, DIMENSION(*) :: grsum,grsum0,gvsum,gvsum0
      REAL*8, DIMENSION(*) :: gcur0ka, gcurfka
!============
      do 100 iig=1,ngroupt
      ig=nogroupt(iig)
      gcurr(ig) = 0._R8
      do 90 l=1,ntpts
   90 gcurr(ig) = gcurr(ig) + fact(l)*gcur(l,ig)
      gcur0ka(ig) = gcurr(ig)*.001_R8
      gcurfka(ig) = gcurfb(ig)*udsi*.001_R8
  100 continue
      do 186 iig=1,ngroupt
      ig=nogroupt(iig)
      sum = 0._R8
      sum0 = 0._R8
      sumv = 0._R8
      sumt = 0._R8
      sumve = 0._R8
!
      if(kcycle.le.0 .or. icirc.eq.0) go to 106
      do 105 jg=1,nreg
      if(ig.ne.jnreg(jg)) go to 105
      sumve = veg(jg)*udsv*.001_R8
      if(aturneg(jg).ne.0) sumve = sumve/aturneg(jg)
  105 continue
  106 continue
      if(nwire.eq.0) go to 188
      do 187 iw1=1,nwire
      if(iabs(igroupw(iw1)).ne.ig) go to 187
      ic1 = iw1 + ncoil-nwire
      fac = 1._R8
      if(isym.eq.1 .and. zwire(iw1).ne.0) fac=2._R8
      sum = sum + (ccoil(ic1)*udsi*.001_R8*fac)
      sum0 = sum0 + (aturnsw(iw1)*gcurr(ig)*fac*.001_R8)
      sumv = sumv + (aturnsw(iw1)*resave(iw1)*tpi*udsv*.001_R8)
      sumt = sumt + (aturnsw(iw1))
  187 continue
  188 continue
      sumvi = 0._R8
      if(sumt.eq.0) go to 189
      sumvi = sumv/sumt
  189 continue
      grsum(ig) = sum
      grsum0(ig) = sum0
      gvsum(ig) = sumvi + sumve
      gvsum0(ig) = gvolt0(ig)*.001_R8
  186 continue
!
   88 if(ncoil-nwire.le.0) go to 89
!
      do 286 iig=1,ngroupt
      ig=nogroupt(iig)
      do 287 ic=1,ncoil-nwire
      if(iabs(igroupc(ic)).ne.ig) go to 287
      fac = 1._R8
      if(isym.eq.1 .and. zcoil(ic).ne.0) fac=2._R8
      grsum(ig) = grsum(ig)+ccoil(ic)*udsi*.001_R8*fac
      grsum0(ig) = grsum0(ig) + aturnsc(ic)*gcurr(ig)*fac*.001_R8
      if(icirc.ge.1) go to 287
      gvsum(ig) = gvsum(ig) + ccoil(ic)*rscoil(ic)*fac*udsv*.001_R8
  287 continue
  286 continue
   89 continue
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

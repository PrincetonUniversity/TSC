      subroutine limpsi
!
      USE CLINAM
      USE WALLCL
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!.................................................................
!....calculate psilim
!.................................................................
!....note  iplim = 0  to initialize
!                = 2  top divertor
!                = 3  bot divertor
!                = 4  top corner
!                = 5  bot corner
!                = 6  other limiter
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER i,j,n,numsepo,nsep,nsepo,imin,imax,jmin,jmax,jzz,jz
      INTEGER izz,iz,itside,itshad,isearch,iplus,jplus,iindx,jindx
      INTEGER isw,iabs,izone,jzone,itest,jtest,l,ll,iplimm,nnsave
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 psinn,pinterp,xtry,ztry,gradsq,dpsidx,dpsidz,gsval
      REAL*8 psval,psixz,psixx,psizz,denom,delx,delz,sinth,costh
      REAL*8 xtest,ztest,ptest,term1,term2,dis,sepsqr,sepsqd,sepsqu
      REAL*8 sepsql,psitry
!============
      do 301 i=2,nxp
      do 301 j=2,nzp
  301 iexs(i,j)=0
!
      psilim = 1.E20_R8
      iplim = 0
      if(igone.eq.1 .or. lrswtch .ge. 1) then
      psilim = -1.E20_R8
      do   nn=1,nlim
        psinn = pinterp(xlima(nn),zlima(nn),ilima(nn),jlima(nn))
        psilim = max(psilim,psinn)
      enddo
      if(idiv.eq.0) return
!
      endif
!
!.....bypass separatrix logic for idiv=0
      if(idiv.eq.0) go to 242
!
!
!.....Reset some variables. (Not more than 4 x-points assumed.)
!.....xsep(n) & zsep(n) are coordinates for x-points.
!.....izsep(n) & jzsep(n) are the corresponding grid numbers.
!.....psep(n) are psi-values at x-points.
!.....nsep is number of x-points found so far.
!.....numsep is set to total number of x-points.
!.....Variables ending in "o" are used to store values from last time-step.
      do 210 n=1,nsepmax
        xsep(n) = 0._R8
        zsep(n) = 0._R8
        psep(n) = 1.E20_R8
        izsepo(n) = izsep(n)
        izsep(n) = 0._R8
        jzsepo(n) = jzsep(n)
        jzsep(n) = 0._R8
  210 continue
      numsepo = numsep
      numsep = 0
      nsep=1
      nsepo = 0
!
!
!.....Start of main loop. One x-point is searched for on each pass.
  211 continue
      nsepo = nsepo + 1
!.....Limit search to square specified in input data.
!.....If no symmetry, use mirror of upper border as lower border.
      imin = iminsep
      imax = imaxsep
      jmin = jminsep
      if (isym.ne.1) jmin = 2*nh-jmaxsep
      jmax = jmaxsep
!.....If equilibrium calculation, don't limit x-point search further.
      if(kcycle.le.0) go to 219
!
!.....only search for new xpoints every 100 cycles
      if(nsepo.gt.numsepo .and. mod(kcycle,10) .ne. 0) go to 227
!.....If searching for x-point not found in previous time-step,
!.....we cannot limit search further.
      if(nsepo.gt.numsepo) go to 219
!.....Limit search for x-point to small square around x-point found in
!.....last time-step.
!     imin = izsepo(nsepo) - int(acoef(45))
!     imax = izsepo(nsepo) + int(acoef(45))
!     jmin = jzsepo(nsepo) - int(acoef(45))
!     jmax = jzsepo(nsepo) + int(acoef(45))
      imin = izsepo(nsepo) -     acoef(45) 
      imax = izsepo(nsepo) +     acoef(45) 
      jmin = jzsepo(nsepo) -     acoef(45) 
      jmax = jzsepo(nsepo) +     acoef(45) 
      imin = max(imin,3)
      imax = min(imax,nx)
      jmin = max(jmin,2)
      jmax = min(jmax,nz)
  219 continue
!.....Search in square defined above.
      do 226 jzz=jmin,jmax
        jz = jzz
!.......Don't search between lower boarder given in input, and it's
!.......mirror image.
        if (isym.ne.1.and.jz.lt.jminsep.and.jz.gt.2*nh-jminsep)          &  
     &  go to 226
        do 101 izz=imin,imax
      iz = izz
!
!....limit search to plasma region, avoiding region shaddowed by x-point
!
      itside=1
      itshad=1
      isearch=1
      do 99 iplus=-isearch,isearch
      do 99 jplus=-isearch,isearch
      iindx = iz+iplus
      if(iindx.gt.nxp) iindx=nxp
      if(iindx.lt.2 ) iindx=2
      jindx = jz+jplus
      if(jindx.gt.nz) jindx=nz
      if(jindx.lt.2 ) jindx=2
      itside = itside*iexv(iindx,jindx)
      itshad = itshad*iexs(iindx,jindx)
   99 continue
      if(itshad.gt.0) go to 101
      if(itside.gt.0) go to 101
!
!.........do newtons iteration to locate x-point.
          isw = 1
          xtry = .5_R8*(xary(iz)+xary(iz+1))
          ztry = .5_R8*(zary(jz)+zary(jz+1))
          call grap(3,ztry,xtry,gradsq,dpsidx,dpsidz,gsval,psval,        &  
     &    psixz,psixx,psizz,isw)
          isw = 0
          do 100 n=1,5
            denom = psixx*psizz-psixz**2
            if(denom.eq.0) go to 101
            delx = (-psizz*dpsidx + psixz*dpsidz)/denom
            delz = ( psixz*dpsidx - psixx*dpsidz)/denom
      if(abs(delx).lt.1.E-6_R8*deex .and. abs(delz).lt.1.E-6_R8*deex)    &  
     &   go to 1011
            xtry = xtry + delx
            ztry = ztry + delz
            if(xtry.lt.xary(iz-1) .or. xtry.gt.xary(iz+2)) go to 101
            if(ztry.lt.zary(jz-1) .or. ztry.gt.zary(jz+2)) go to 101
            call grap(3,ztry,xtry,gradsq,dpsidx,dpsidz,gsval,psval,      &  
     &      psixz,psixx,psizz,isw)
  100     continue
 1011 continue
          if(denom.gt.0) go to 101
          if(gradsq.gt.1.E-8_R8) go to 101
!.........Discard x-point if same as found earlier.
          if (nsep.le.1) go to 221
          do 220 n=1,nsep-1
            if (iabs(izsep(n)-iz).lt.3.and.iabs(jzsep(n)-jz).lt.3)       &  
     &      go to 101
  220     continue
!
  221     continue
      izone = (xtry-ccon)/deex+3
      jzone = (ztry-zzero)/deez+3
!
      if(acoef(895).eq.1.0_R8) go to 1222
      iplus = -1
      if(xtry.lt.xmag) iplus = 0
      jplus = -1
      if(ztry.lt.zmag) jplus = 0
      if(iexv(izone+iplus,jzone+jplus).eq.1) go to 3101
      if(xmag.le.0) go to 223
 1222 continue
!.......check if psi decreases towards plasma
        denom = sqrt((ztry-zmag)**2+(xtry-xmag)**2)
        sinth = (ztry-zmag)/denom
        costh = (xtry-xmag)/denom
        xtest = xtry - deex*costh
        ztest = ztry - deex*sinth
        itest = (xtest-ccon)/deex + 2
        jtest = (ztest+alz*(1-isym))/deez + 2
        ptest = pinterp(xtest,ztest,itest,jtest)
        if(ptest.gt.psval) go to 101
  223   continue
!
!.........Save x-point found.
          xsep(nsep) = xtry
          zsep(nsep) = ztry
          psep(nsep) = psval
          izsep(nsep) = iz
          jzsep(nsep) = jz
      if(xmag.le.0) go to 303
      do 300 i=2,nxp
      do 300 j=2,nzp
      term1 = (zmag-zsep(nsep))*(zary(j)-zsep(nsep))
      term2 =-(xmag-xsep(nsep))*(xary(i)-xsep(nsep))
      dis = sqrt((zary(j)-zsep(nsep))**2 + (xary(i)-xsep(nsep))**2)
      if(term1 .lt. term2 .and. dis .gt. deex) iexs(i,j)=1
  300 continue
  303 continue
          nsep=nsep+1
!.........If not more than nsepmax-1 x-points have been found, loop back
!.........and search for next.
          if(nsep.le.nsepmax) go to 211
!.........Ready, found nsepmax x-points.
          go to 227
 3101 continue
  101   continue
  226 continue
!.....Ready found nsep-1 x-points.
      if (nsepo.le.numsepo) go to 211
  227 numsep=nsep-1
!
!
!.....Bypass rest of routine if no x-point.
      psisep = psimin
      if (numsep.eq.0) go to 242
!.....Determine smallest square that encloses all x-points and smallest
!.....psi at x-points.
      psisep=1.E20_R8
      sepsqr=-1.E20_R8
!.....Set psilim to smallest value of psi at x-points
!
      sepsqd=1.E20_R8
      sepsqu=-1.E20_R8
      do 231 n=1,numsep
        sepsql=min(sepsql,xsep(n))
        sepsqr=max(sepsqr,xsep(n))
        sepsqd=min(sepsqd,zsep(n))
        sepsqu=max(sepsqu,zsep(n))
      if(psep(n).gt.psisep) go to 231
        psisep=min(psep(n),psisep)
!....set iplim to minus the limiting x-point
      iplim = -n
  231 continue
!
!
      if(igone.gt.0) return
!
!
      psilim=psirat*(psisep-psimin) + psimin
      if(psimin .eq. 0) psilim = psisep
  242 continue
      if(igone.gt.0) return
!....set up mask array to exclude plasma on far side of separatrix
!
!
!
      if(acoef(811).le.0) go to 245
      if(isurf     .eq.0) go to 245
!
!.....for acoef(811)>0, limit plasma at surface where q=acoef(811)
      do 244 l=2,npsit
      ll = l
      if(qprof2(l).gt.acoef(811)) go to 246
  244 continue
  246 continue
      denom = qprof2(ll)-qprof2(ll-1)
      if(denom.le.0) go to 245
      psitry = xsv2(ll-1) + (acoef(811)-qprof2(ll-1))*                   &  
     &                      (xsv2(ll)-xsv2(ll-1))/denom
      if(psitry.gt.psilim) go to 245
      psilim = psitry
      iplimm = 0
  245 continue
      if(nlim.le.0) go to 241
      plimiter = 1.E60_R8
      do 240 nn=1,nlim
!.......added 1/6/87
      if(iexs(ilima(nn),jlima(nn)) .eq. 1) go to 240
      if(acoef(1) .eq. 1._R8.and. zmag*zlima(nn) .lt. 0) go to 240
        psinn = pinterp(xlima(nn),zlima(nn),ilima(nn),jlima(nn))
        if(psinn.gt.plimiter) go to 240
        nnsave = nn
      if(xmag.le.0.or. igone.ne.0) go to 239
!.......check if psi decreases towards plasma
        denom = sqrt((zlima(nnsave)-zmag)**2+(xlima(nnsave)-xmag)**2)
        sinth = (zlima(nnsave)-zmag)/denom
        costh = (xlima(nnsave)-xmag)/denom
        xtest = xlima(nnsave) - deex*costh
        ztest = zlima(nnsave) - deex*sinth
        itest = (xtest-ccon)/deex + 2
        jtest = (ztest+alz*(1-isym))/deez + 2
        ptest = pinterp(xtest,ztest,itest,jtest)
        if(ptest.gt.psinn) go to 240
  239   plimiter = psinn
        iplimm = nn
  240 continue
      if(plimiter .ge. psilim) go to 241
      psilim = plimiter
      iplim = iplimm
  241 continue
!     write(nterm,1241) iplim,numsep
!1241 format(" iplim,numsep =",2i4)
      if(lrswtch.eq.0 .and. psilim.ge.1.E20_R8) ineg=34
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

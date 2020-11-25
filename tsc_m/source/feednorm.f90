      subroutine feednorm
!......5.91 feednorm
!
! until now  azerv resp. rzerv replaced xmag; now suppressed.           iern
!.....calculates normalized feedback gains
!
      USE CLINAM
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER nfeed,l,nnn,ig,ie,ll,n,np,ii,ngr,iabs,ngrvwl,iepass
      INTEGER ngrvcl
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 term,tcoilii,xwireii,zwireii,an,anp,aindx,xuse,zuse
      REAL*8 ause,gradpsi,psipass,rmajpass,aminpas
!     REAL*8 sum
!============
!     dimension sum(ptpts)
!============      
      INTEGER :: istat = 0 
      REAL*8, ALLOCATABLE, DIMENSION(:) :: sum
!============      
      IF(.not.ALLOCATED(sum)) ALLOCATE( sum(ptpts), STAT=istat)
!============      
      if (istat .ne. 0) stop 'Allocation Error : feednorm  ' 
!============      
!
      do 2000 nfeed=1,3
!
      go to(101,102,103),nfeed
  101 write(nout,1500)
      go to 104
  102 write(nout,1501)
      go to 104
  103 write(nout,1502)
  104 continue
      do 811 l=1,numfb
      do 813 nnn=1,ngroupt
      nn = nogroupt(nnn)
  813 gcurfb(nn) = 0._R8
!
      go to(201,202,203),nfeed
  201 term = fbfac(l)
      go to 204
  202 term = fbfacd(l)
      go to 204
  203 term = fbfaci(l)
  204 continue
      if(term.eq.0) go to 811
!
      if(nrfb(l).le.0) go to 811
      ig = nrfb(l)
      gcurfb(ig) = term*usdi
!
      ie = ipext(l)
      if(ie.eq.2 .or. ie.eq.3) go to 811
      if(ie.ge.15 .and. ie.le.20) go to 811
      if(ie.gt.24)               go to 811
      do 812 ll=1,ntpts
  812 sum(ll) = 0._R8
!
      do 500 ll=1,ntpts
      if(tfbon(l).ge.0 .and. tpro(ll).lt.tfbon(l)) go to 500
      if(tfbof(l).ge.0 .and. tpro(ll).gt.tfbof(l)) go to 500
!
      n = 2*nfeedv(ll,l)-1
      np = n+1
!
      if(ie.eq.1) then
        if(n .le. 0) go to 700
        if(xobs(n) .le.0) go to 700
        if(xobs(np).le.0) go to 700
      endif
!
      do 830 ii=1,nwire
      ngr = iabs(igroupw(ii))
      tcoilii = (gcurfb(ngr))*aturnsw(ii)
!
      ngrvwl = ngrvw1(ii)
      if(ngrvwl.ne.ig) go to 829
      tcoilii = tcoilii+atnvw1(ll,ii)*gcurfb(ngrvwl)
!
  829 continue
      ngrvwl = ngrvw2(ii)
      if(ngrvwl.ne.ig) go to 839
      tcoilii = tcoilii+atnvw2(ll,ii)*gcurfb(ngrvwl)
!
  839 continue
      ngrvwl = ngrvw3(ii)
      if(ngrvwl.ne.ig) go to 849
      tcoilii = tcoilii+atnvw3(ll,ii)*gcurfb(ngrvwl)
!
  849 continue
      ngrvwl = ngrvw4(ii)
      if(ngrvwl.ne.ig) go to 859
      tcoilii = tcoilii+atnvw4(ll,ii)*gcurfb(ngrvwl)
!
  859 continue
!
      if(tcoilii.eq.0) go to 830
      xwireii = xary(iwire(ii))
      zwireii = zary(jwire(ii))
      if(ie.eq.1) then
      call gf(ineg,nmult,xobs(n ),zobs(n ),xwireii,zwireii,an )
      call gf(ineg,nmult,xobs(np),zobs(np),xwireii,zwireii,anp)
      sum(ll) = sum(ll) + tcoilii*(an -anp)/tpi
      endif
!
      if(ie.eq.4) then
      call gf(ineg,nmult,xplas,zplas,xwireii,zwireii,an)
      aindx = xplas*(log(8._R8*xplas/aminor)-1.5_R8)*usdi
      sum(ll) = sum(ll) + tcoilii*an/aindx
      endif
      if(ie.le.4  .or. ie.ge.15) go to 910
      if(ie.le.10 .and. ie.ge. 7)go to 910
!
      xuse = xplas
      zuse = zplas
      ause = aminor
      if(ie.eq.5 .and. xmagz(ll).gt.0) xuse = xmagz(ll)
!@@@a
!     if(rzerv(ll).gt.0) xuse = rzerv(ll)
!     if(azerv(ll).gt.0) xuse = azerv(ll)
!@@@e
      gradpsi = xplas*pcur(ll)*usdi/(tpi*aminor)
      if(ie.eq.5) iepass=1
      if(ie.eq.6) iepass=5
      if(ie.eq.11) iepass=1
      if(ie.eq.12) iepass=4
      if(ie.eq.13) iepass=2
      if(ie.eq.14) iepass=3
      call fluxcal(xwireii,zwireii,iepass,psipass,xuse,ause)
      sum(ll) = sum(ll) - tcoilii*psipass/tpi/gradpsi
  910 continue
!
!
      if(ie.gt.20) then
      iepass = ie-20
      call fluxcal(xwireii,zwireii,iepass,psipass,rmajpass,aminpas)
      sum(ll) = sum(ll) + tcoilii*psipass/tpi
      endif
!
      if(isym.eq.0) go to 830
      if(zwireii.eq.0) go to 830
      if(ie.eq.1) then
      call gf(ineg,nmult,xobs(n ),zobs(n ),xwireii,-zwireii,an )
      call gf(ineg,nmult,xobs(np),zobs(np),xwireii,-zwireii,anp)
      sum(ll) = sum(ll) + tcoilii*(an -anp)/tpi
      endif
!
      if(ie.eq.4) then
      call gf(ineg,nmult,xplas,zplas,xwireii,-zwireii,an)
      aindx = xplas*(log(8._R8*xplas/aminor)-1.5_R8)*usdi
      sum(ll) = sum(ll) + tcoilii*an/aindx
      endif
      if(ie.le.4  .or. ie.ge.15) go to 911
      if(ie.le.10 .and. ie.ge. 7)go to 911
!
      xuse = xplas
      zuse = zplas
      ause = aminor
      if(ie.eq.5 .and. xmagz(ll).gt.0) xuse = xmagz(ll)
!@@@a
!     if(rzerv(ll).gt.0) xuse = rzerv(ll)
!     if(azerv(ll).gt.0) xuse = azerv(ll)
!@@@e
      gradpsi = xplas*pcur(ll)*usdi/(tpi*aminor)
      if(ie.eq.5) iepass=1
      if(ie.eq.6) iepass=5
      if(ie.eq.11) iepass=1
      if(ie.eq.12) iepass=4
      if(ie.eq.13) iepass=2
      if(ie.eq.14) iepass=3
      call fluxcal(xwireii,-zwireii,iepass,psipass,xuse,ause)
      sum(ll) = sum(ll) - tcoilii*psipass/tpi/gradpsi
  911 continue
!
!
      if(ie.gt.20) then
      iepass = ie-20
      call fluxcal(xwireii,-zwireii,iepass,psipass,rmajpass,aminpas)
      sum(ll) = sum(ll) + tcoilii*psipass/tpi
      endif
  830 continue
!
!.....recompute external coils
      if(ncoil.eq.nwire) go to 890
      do 896 ii=1,ncoil-nwire
      ngr = iabs(igroupc(ii))
      tcoilii = aturnsc(ii)*(gcurfb(ngr))
!
      ngrvcl = ngrvc1(ii)
      if(ngrvcl.ne.ig) go to 929
      tcoilii = tcoilii+atnvc1(ll,ii)*gcurfb(ngrvcl)
!
  929 continue
      ngrvcl = ngrvc2(ii)
      if(ngrvcl.ne.ig) go to 939
      tcoilii = tcoilii+atnvc2(ll,ii)*gcurfb(ngrvcl)
!
  939 continue
      ngrvcl = ngrvc3(ii)
      if(ngrvcl.ne.ig) go to 949
      tcoilii = tcoilii+atnvc3(ll,ii)*gcurfb(ngrvcl)
!
  949 continue
      ngrvcl = ngrvc4(ii)
      if(ngrvcl.ne.ig) go to 959
      tcoilii = tcoilii+atnvc4(ll,ii)*gcurfb(ngrvcl)
!
  959 continue
      if(tcoilii.eq.0) go to 896
      if(ie.eq.1) then
      call gf(ineg,nmult,xobs(n ),zobs(n ),xcoil(ii),zcoil(ii),an )
      call gf(ineg,nmult,xobs(np),zobs(np),xcoil(ii),zcoil(ii),anp)
      sum(ll) = sum(ll) + tcoilii*(an -anp)/tpi
      endif
!
      if(ie.eq.4) then
      call gf(ineg,nmult,xplas,zplas,xcoil(ii),zcoil(ii),an)
      aindx = xplas*(log(8._R8*xplas/aminor)-1.5_R8)*usdi
      sum(ll) = sum(ll) + tcoilii*an/aindx
      endif
      if(ie.le.4  .or. ie.ge.15) go to 912
      if(ie.le.10 .and. ie.ge. 7)go to 912
!
      xuse = xplas
      zuse = zplas
      ause = aminor
      if(ie.eq.5 .and. xmagz(ll).gt.0) xuse = xmagz(ll)
!@@@a
!     if(rzerv(ll).gt.0) xuse = rzerv(ll)
!     if(azerv(ll).gt.0) xuse = azerv(ll)
!@@@e
      gradpsi = xplas*pcur(ll)*usdi/(tpi*aminor)
      if(ie.eq.5) iepass=1
      if(ie.eq.6) iepass=5
      if(ie.eq.11) iepass=1
      if(ie.eq.12) iepass=4
      if(ie.eq.13) iepass=2
      if(ie.eq.14) iepass=3
      call fluxcal(xcoil(ii),zcoil(ii),iepass,psipass,xuse,ause)
      sum(ll) = sum(ll) - tcoilii*psipass/tpi/gradpsi
  912 continue
!
!
      if(ie.gt.20) then
      iepass = ie-20
      call fluxcal(xcoil(ii),zcoil(ii),iepass,psipass,rmajpass,aminpas)
      sum(ll) = sum(ll) + tcoilii*psipass/tpi
      endif
!
      if(isym.eq.0) go to 896
      if(zcoil(ii).eq.0) go to 896
      if(ie.eq.1) then
      call gf(ineg,nmult,xobs(n ),zobs(n ),xcoil(ii),-zcoil(ii),an )
      call gf(ineg,nmult,xobs(np),zobs(np),xcoil(ii),-zcoil(ii),anp)
      sum(ll) = sum(ll) + tcoilii*(an -anp)/tpi
      endif
!
      if(ie.eq.4) then
      call gf(ineg,nmult,xplas,zplas,xcoil(ii),-zcoil(ii),an)
      aindx = xplas*(log(8._R8*xplas/aminor)-1.5_R8)*usdi
      sum(ll) = sum(ll) + tcoilii*an/aindx
      endif
      if(ie.le.4  .or. ie.ge.15) go to 913
      if(ie.le.10 .and. ie.ge. 7)go to 913
!
      xuse = xplas
      zuse = zplas
      ause = aminor
      if(ie.eq.5 .and. xmagz(ll).gt.0) xuse = xmagz(ll)
!@@@a
!     if(rzerv(ll).gt.0) xuse = rzerv(ll)
!     if(azerv(ll).gt.0) xuse = azerv(ll)
!@@@e
      gradpsi = xplas*pcur(ll)*usdi/(tpi*aminor)
      if(ie.eq.5) iepass=1
      if(ie.eq.6) iepass=5
      if(ie.eq.11) iepass=1
      if(ie.eq.12) iepass=4
      if(ie.eq.13) iepass=2
      if(ie.eq.14) iepass=3
      call fluxcal(xcoil(ii),-zcoil(ii),iepass,psipass,xuse,ause)
      sum(ll) = sum(ll) - tcoilii*psipass/tpi/gradpsi
  913 continue
!
!
      if(ie.gt.20) then
      iepass = ie-20
      call fluxcal(xcoil(ii),-zcoil(ii),iepass,psipass,rmajpass,aminpas)     
      sum(ll) = sum(ll) + tcoilii*psipass/tpi
      endif
!
  896 continue
  890 continue
!
  500 continue
!
      write(nout,2500) l,(sum(ll),ll=1,ntpts)
  811 continue
 2000 continue
 1500 format(/" proportional gains for fback(dimensionless)ll=1,ntpts")
 1501 format(/" derivative gains for feedback(sec)**(+1)   ll=1,ntpts")
 1502 format(/" integral gains for feedback  (sec)**(-1)   ll=1,ntpts")
 2500 format(2x,i5,10f10.4,10(/,7x,10f10.4))
!
      return
!
  700 continue
      ineg=20
      write(nout,1700) l
 1700 format(" observation point not defined for feedback system",i3)
      return
510      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

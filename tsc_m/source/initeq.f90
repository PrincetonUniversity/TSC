      subroutine initeq
!
!.....solve for initial equilibrium
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      USE CLINAM
      USE WALLCL
      USE CBSVD
      USE SCR11

      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER ii,n,i,j,itmax,iabs,itag,iter,iquit
      INTEGER istrike,itval,k,iold,inew,i1,i2
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 gppmin,psimin0,zcurf0,rstr,zstr,eps1,eps2,denom
      REAL*8 denom2,epsmax,gpt,pr1,pr2
!============
!     common/cbsvd/sigmax,nranksc,ngscfb
      if(irst1.eq.1) go to 9544
!
!.....................................................................
!.....solve elliptic problem for poloidal flx
!.....................................................................
!
!.....define inhomogeneous terms for poloidal flux equation
      do 150 ii=1,nwire
      n = ncoil-nwire+ii
      i = iwire(ii)
      j = jwire(ii)
  150 vd(i,j) = ccoil(n)*xary(i)
 8154 format(1p,10e12.4)
!
!.....initialize GS array
      do 155 i=1,nxp
      do 155 j=1,nzp
  155 gs(i,j) = gzero
!
      xmag = xplas
      zmag = zplas
      imag = (xplas-ccon)/deex + 2
      jmag = (zplas-zzero)/deez + 2.0000001_R8
      if(lrswtch.ge.1.or.igone.ge.1) go to 151
      vd(imag,jmag) = pcur(istart)*usdi*xplas
      if(xzeric.ne.0) then
      call curinit
      xmag = xzeric
      zmag = zzeric
      imag = (xzeric-ccon)/deex + 2
      jmag = (zzeric-zzero)/deez + 2.0000001_R8
      endif
  151 continue
!
      if(iflux.eq.4) call elliptic(4,psizer)
      call elliptic(1,psi)
!..rxw/05/10/87
!===      if(lrswtch.ge.1) go to 16
!.....initialize subrotine appvolt2 if acoef(290)=2:
      if(acoef(290).eq.2) call appvolt2
!...rxw/end
!
!......................................................................
!......start iteration loop to solve free boundary grad-shafranov eqn
!......................................................................
 9544 continue
      if(ifunc.eq.2) then
      write(nterm,1103)
      write(nout, 1102)
      else
      write (nterm, 1003)
      write(nout,1002)
!
 1003 format(" itr       gp1         gp2      psimin        psilim",     &  
     &"     xmag       zmag        eps")
 1002 format("1      gp1          p0         psimin      xmag",          &  
     & "        zmag         eps      psilim      xcurf       zcurf",    &  
     & " ipl" )
!
 1103 format(" itr       gp1          p0      psimin        psilim",     &  
     &"     xmag       zmag        eps")
 1102 format("1      gp1         gp2         psimin      xmag",          &  
     & "        zmag         eps      psilim      xcurf       zcurf",    &  
     & " ipl" )
      endif
      itmax = 200
      if(iabs(neqmax).gt.1._R8) itmax = iabs(neqmax)
      itag=0
      do 15 iter=1,itmax
      call newj2(iter)
      call icalc
      if(ineg.ne.0) return
      if(iflux.eq.4) call elliptic(4,psizer)
      call elliptic(1,psi)
!..rxw/05/10/87
      if(lrswtch.ge.1.or.igone.ge.1) go to 15
!...rxw/end
      call geval(psimin,1,gmin,gpmin,gppmin,imag,jmag)

      if(abs(psimin-psimin0).lt.acoef(23)*acoef(50)*abs(psilim-psimin)   &  
     &   .and.   abs(zcurf-zcurf0).lt.acoef(24)*alz) then
        if(acoef(901).eq.0._R8.or.iter.lt.acoef(907)) goto 160
        write(nterm,1843)
 1843 format('convergence criterion is satisfied, but check svd          &  
     & residual                                                          &  
     &s')
        write(nterm,1844)
 1844 format('type 0 to end, 1 to change input parameters')
        read(nterm,*) iquit
        if(iquit.eq.0) go to 160
        if(iquit.eq.1) then
! if running with x-point constraint, test if want to fix strike points
! instead.
        if(acoef(901).eq.3 .or. acoef(901).eq.4._R8) then
        write(nterm,7354)
 7354   format('to change x-point constraint to strike point constraint,   &  
     &   en                                                                &  
     &   ter 1, otherwise, enter 0')
         read(nterm,*) istrike
        if(istrike.eq.0) goto 7777
        if(istrike.eq.1) then
          write(nterm,7356)
 7356     format('enter strike point location')
          read(nterm,*) rstr,zstr
          xcon0(istart,ncnt+1)=rstr
          zcon0(istart,ncnt+1)=zstr
          ncnt=ncnt+1
          if(acoef(901).eq.4._R8) acoef(901)=2._R8
          if(acoef(901).eq.3._R8) acoef(901)=1._R8
        endif
      endif
 7777 write(nterm,7457) acoef(906),acoef(50),acoef(909)
      write(nterm,7458)
      read(nterm,*) acoef(906),acoef(50),acoef(909)
      endif
      endif
      eps1 = 0._R8
      eps2 = 0._R8
      denom = acoef(23)*acoef(50)*abs(psilim-psimin)
      if(denom.ne.0) eps1 = abs(psimin-psimin0)/denom
      denom2= acoef(24)*alz
      if(denom2.ne.0) eps2 = abs(zcurf-zcurf0)/denom2
      epsmax = max(eps1,eps2)
170   continue
      zcurf0 = zcurf
      psimin0 = psimin
      gpt = gp2
      if(ifunc.eq.2) gpt = p0*udsp
      write(nout,1001) iter,gp1,gpt,psimin,xmag,zmag,epsmax,psilim,      &  
     &    xcurf,zcurf,iplim
 1001 format(i4,1p9e12.4,i3)

        if (mod(iter,10).ne.0)   go to 15
        write (nterm, 10001) iter,gp1,gpt,psimin,psilim,xmag,zmag,       &  
     &epsmax,iplim
10001   format(i4,1p7e12.4,i5)
      if(acoef(901).eq.0._R8.or. iter.lt.acoef(907)) goto 15
      itval=acoef(910)
      if(iter/itval*itval.eq.iter) then
      write(nterm,7457) acoef(906),acoef(50),acoef(909)
7457  format('sigma-cutoff, acoef(50), acoef(909) = ',1p10e12.4)
      write(nterm,6679)
 6679 format('singular values are:')
      write(nterm,6680) (sigmasvd(k),k=1,ngscfb)
 6680 format(1p7e12.4)
      write(nterm,6681)
 6681 format('values scaled to sigmax=1. are:')
      write(nterm,6680) (sigmasvd(k)/sigmax,k=1,ngscfb)
      write(nterm,7458)
7458  format('enter new sigma-cutoff, acoef(50), acoef(909)')
      read(nterm,*) acoef(906),acoef(50),acoef(909)
      endif
      goto 15
160   continue
!
        write (nterm, 1001)   iter, gp1,gp2, psimin, psilim, xmag, zmag,  &  
!    &                                                                   &  
     &   epsmax
! after convergence with multiple limiter points, remove and
! replace with a single limiter point on inboard edge.
      if(acoef(74).ne.1._R8) goto 16
!
!..............special limiter adjustment section for acoef(74) .ne. 0 .
!
      if(itag.eq.1) goto 16
      itag=itag+1
      iold=imag
826   continue
      inew=iold-1
      pr1=psi(inew,nzp/2+1)-psilim
      pr2=psi(iold,nzp/2+1)-psilim
      if(pr1*pr2.le.0._R8) goto 827
      iold=inew
      goto 826
827   continue
      do 828 i=1,nlim
      xlima(i)=0._R8
      zlima(i)=0._R8
      ilima(i)=0._R8
      jlima(i)=0._R8
828   continue
      nlim=1
      xlima(1)=(xary(iold)*pr1-xary(inew)*pr2)/(pr1-pr2)
      zlima(1)=0._R8
      ilima(1)= int((xlima(1)-ccon)/deex   + 2.5_R8)
      jlima(1)= int((zlima(1)-zzero)/deez   + 2.5_R8)
!
!..............end of special limiter adjustment section .............
!
      goto 170
   15 continue
      if(neqmax.lt.0 .or. lrswtch.gt.0 .or. igone.ge.1) go to 16
      ineg=12
   16 continue
      if(acoef(76).eq.0) go to 17
      i1 = acoef(77)
      i2 = acoef(78)
      do 18 i=i1,i2
   18 fbfac(i) = 0._R8
   17 continue
!
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

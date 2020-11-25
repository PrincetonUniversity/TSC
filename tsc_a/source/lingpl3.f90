      subroutine lingpl3
!
!............  J.E calculation
!
      USE CLINAM
      USE SCR1
      USE SCR12A
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!
!
!.....------>   This subroutine needs a lot of work
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER n,i,j,indx
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 psiinc,ptest,fac,xmin,xmax,ymin,ymax
!============
      return
      go to 501
      entry lingin3
      psiinc = (psilim-psimin)/npts
      do 200 n=1,npts+1
      uplot2(n) = 0._R8
  200 cplot(n) = 0._R8
!
!........ Jp.Ep  ..............
!
      do 400 i=iminn,imaxx
      do 300 j=jminn,jmaxx
      if(iexv(i,j).eq.1 .or. iexs(i,j).eq.1) go to 300
      ptest = .25_R8*(psi(i-1,j)+psi(i,j)+psi(i,j-1)+psi(i-1,j-1))
      if(ptest.ge.psilim.or.ptest.lt.psimin) go to 300
      indx = (ptest-psimin)/psiinc + 1
      fac = 1.0_R8
      uplot2(indx) = uplot2(indx) + gling2(i,j)/etay(i,j)*               &  
     &   (gling1(i,j)+gling2(i,j)+gling3(i,j))/udst
      cplot(indx) = cplot(indx) + fac
  300 continue
  400 continue
      xmin = psimin*tpi
      xmax = psilim*tpi
      ymin = 0._R8
      ymax = 0._R8
      do 499 n=1,npts
      indx = npts+1-n
      if(cplot(indx).eq.0) go to 550
      uplot2(indx) = uplot2(indx)/cplot(indx)
      go to 499
  550 continue
      uplot2(indx) = uplot2(indx+1)
  499 continue
!
!.....   Jt.Et ...........................
!
      do 571 n=1,npts+1
      eplot2(n) = 0._R8
  571 cplot(n) = 0._R8
      do 590 i=iminn,imaxx
      do 580 j=jminn,jmaxx
      if(iexv(i,j).eq.1 .or. iexs(i,j).eq.1) go to 580
      ptest = psi(i,j)
      if(ptest.ge.psilim.or.ptest.lt.psimin) go to 580
      indx = (ptest-psimin)/psiinc + 1
      fac = 1.0_R8
      if(j.ne.2 .and. isym.eq.1) fac = 2.0_R8
      eplot2(indx) = eplot2(indx) + ajp2(i,j)/xary(i)*                   &  
     &   (pling5(i,j)+etay(i,j)*ajp2(i,j)-gs(i,j)*ajp4(i,j))/usdt*fac
      cplot(indx) = cplot(indx) + fac
  580 continue
  590 continue
      do 699 n=1,npts
      indx = npts+1-n
      if(cplot(indx).eq.0) go to 650
      eplot2(indx) = eplot2(indx)/cplot(indx)
      go to 699
  650 continue
      eplot2(indx) = eplot2(indx+1)
  699 continue
!
!......add smoothing
!
      do 401 n=2,npts-1
      uplot(n) = .25_R8*uplot2(n-1) + .5_R8*uplot2(n) + .25_R8*uplot2(n+  &  
     & 1)
      eplot(n) = .25_R8*eplot2(n-1) + .5_R8*eplot2(n) + .25_R8*eplot2(n+  &  
     & 1)
  401 continue
      uplot(npts) = .667_R8*uplot2(npts) + .333_R8*uplot2(npts-1)
      eplot(npts) = .667_R8*eplot2(npts) + .333_R8*eplot2(npts-1)
      uplot(1) = uplot2(1)
      eplot(1) = eplot2(1)
      return
  501 continue
!
      do 500 n=1,npts
      indx = npts+1-n
      tplot(indx) = uplot(indx)+eplot(indx)
      ymin = min(uplot(indx),ymin)
      ymin = min(eplot(indx),ymin)
      ymin = min(tplot(indx),ymin)
      ymax = max(uplot(indx),ymax)
      ymax = max(eplot(indx),ymax)
      ymax = max(tplot(indx),ymax)
      xxplot(indx) = (psimin+(indx-0.5_R8)*psiinc)*tpi
  500 continue
      call maps(xmin,xmax,ymin,ymax,.142_R8,.750_R8,.300_R8,1._R8)
      call tracec(1h1,xxplot,uplot,npts,-1,-1,0._R8,0._R8)
      call tracec(1h2,xxplot,eplot,npts,-1,-1,0._R8,0._R8)
      call tracec(1ht,xxplot,tplot,npts,-1,-1,0._R8,0._R8)
      call setld(2._R8,30._R8,0,0,2,1)
      write(s100,1001) kcycle
      call gtext(s100,80,0)
      call setold(xmax,ymax,1,0,1,0)
      write(s100,1002)
      call gtext(s100,80,0)
 1001 format(" J.E  N =",i6)
 1002 format(1x,/," 1 ...Jp.Ep",/,                                       &  
     &       " 2 ...Jt.Et     ",/,                                       &  
     &       " T ...Total" )
      call setold(xmin,ymin-.12_R8*(ymax-ymin),1,0,1,0)
      write(s100,1003)
      call gtext(s100,80,0)
 1003 format(10x," poloidal flux")
      write(nsc1,1001) kcycle
      call frscj(6)
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

      subroutine lingpl2
!
!
      USE CLINAM
      USE SCR1
      USE SCR12B
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER n,i,j,indx
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 psiinc,ptest,fac,xmin,xmax,ymin,ymax
!============
      go to 501
      entry lingin2
      psiinc = (psilim-psimin)/npts
      do 200 n=1,npts+1
      uplot2(n) = 0._R8
      eplot2(n) = 0._R8
      bplot2(n) = 0._R8
  200 cplot(n) = 0._R8
      do 400 i=iminn,imaxx
      do 300 j=jminn,jmaxx
      ptest = .25_R8*(psi(i-1,j)+psi(i,j)+psi(i,j-1)+psi(i-1,j-1))
      if(ptest.ge.psilim.or.ptest.lt.psimin) go to 300
      indx = (ptest-psimin)/psiinc + 1
      fac = 1.0_R8
      uplot2(indx) = uplot2(indx) + gling1(i,j)/udst
      eplot2(indx) = eplot2(indx) + gling2(i,j)/udst
      bplot2(indx) = bplot2(indx) + gling3(i,j)/udst
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
      eplot2(indx) = eplot2(indx)/cplot(indx)
      bplot2(indx) = bplot2(indx)/cplot(indx)
      go to 499
  550 continue
      uplot2(indx) = uplot2(indx+1)
      eplot2(indx) = eplot2(indx+1)
      bplot2(indx) = bplot2(indx+1)
  499 continue
!
!......add smoothing
      do 401 n=2,npts-1
      uplot(n) = .25_R8*uplot2(n-1) + .5_R8*uplot2(n) + .25_R8*uplot2(n+  &  
     & 1)
      eplot(n) = .25_R8*eplot2(n-1) + .5_R8*eplot2(n) + .25_R8*eplot2(n+  &  
     & 1)
      bplot(n) = .25_R8*bplot2(n-1) + .5_R8*bplot2(n) + .25_R8*bplot2(n+  &  
     & 1)
  401 continue
      uplot(npts) = .667_R8*uplot2(npts) + .333_R8*uplot2(npts-1)
      eplot(npts) = .667_R8*eplot2(npts) + .333_R8*eplot2(npts-1)
      bplot(npts) = .667_R8*bplot2(npts) + .333_R8*bplot2(npts-1)
      uplot(1) = uplot2(1)
      eplot(1) = eplot2(1)
      bplot(1) = bplot2(1)
      return
  501 continue
!
      do 500 n=1,npts
      indx = npts+1-n
      tplot(indx) = uplot(indx)+eplot(indx)+bplot(indx)
      ymin = min(uplot(indx),ymin)
      ymin = min(eplot(indx),ymin)
      ymin = min(bplot(indx),ymin)
      ymin = min(tplot(indx),ymin)
      ymax = max(uplot(indx),ymax)
      ymax = max(eplot(indx),ymax)
      ymax = max(bplot(indx),ymax)
      ymax = max(tplot(indx),ymax)
      xxplot(indx) = (psimin+(indx-0.5_R8)*psiinc)*tpi
  500 continue
      call maps(xmin,xmax,ymin,ymax,.142_R8,.750_R8,.300_R8,1._R8)
      call tracec(1hu,xxplot,uplot,npts,-1,-1,0._R8,0._R8)
      call tracec(1he,xxplot,eplot,npts,-1,-1,0._R8,0._R8)
      call tracec(1hh,xxplot,bplot,npts,-1,-1,0._R8,0._R8)
      call tracec(1ht,xxplot,tplot,npts,-1,-1,0._R8,0._R8)
      call setld(2._R8,30._R8,0,0,2,1)
      write(s100,1001) kcycle
      call gtext(s100,80,0)
      call setold(xmax,ymax,1,0,1,0)
      write(s100,1002)
      call gtext(s100,80,0)
 1001 format(" tor f ohms law  N =",i6)
 1002 format(1x,/," u ...V.Grad g",/,                                    &  
     &       " E ...eta*J     ",/,                                       &  
     &       " H ...Hyper  res" ,/,                                      &  
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

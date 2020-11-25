      subroutine ufstore
!
!
      USE CLINAM
      USE SAPROP
      USE SCR1
      USE SPECIE
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER i,imin,imax,nhp,k,ispcr,nimp,lval,lv,lvalm,nd
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 radpow,xprof,xx,xprofm,ps,pval,ppval,eval,epval,rval
      REAL*8 rpval,tival,sval,rojmin,val,v1,v2,f2
!============
      dimension radpow(21)
!============      
!
!...............................................................
!
! LAYOUT of data in UFDATA
!
! ufdata(nrecord,i,1)  te(i)
!                  2   ti(i)
!                  3   ne(i)
!                  4
!                 5-10 not now used
!                11-28 impurity charge states
!                 29   bolometry cord i
!
!...............................................................
!
!
      xprof = acoef(42)
      if(xprof.le.0) xprof = ccon
      do 40 i=3,nx
      imin = i
      xx = xary(i)
      if(xx.gt.xprof) go to 45
   40 continue
      imin = 3
   45 continue
      xprofm = acoef(43)
      if(xprofm.le.0) xprofm = alx
      do 48 i=3,nx
      imax = i
      xx = xary(i)
      if(xx.gt.xprofm) go to 47
   48 continue
      imax = nx
   47 continue
      nhp = jmag+1
      do 100 i=imin,imax
      if(isurf.eq.1) then
      ps = psi(i,jmag)
      call peval(ps,2,pval,ppval,i,jmag)
      call eeval(ps,2,eval,epval,i,jmag)
      call reval(ps,idens,1,rval,rpval,i,jmag)
      call tieval(ps,2,tival,i,jmag,pval,eval,rval)
      ufdata(nrecord,i,3) = rval*udsd
      ufdata(nrecord,i,1) = (eval*udsh)/(rval*udsd)
      ufdata(nrecord,i,2) = tival
      call seval(ps,sval,i,jmag)
      ufdata(nrecord,i,4) = sval
      else
      ufdata(nrecord,i,2) = .5_R8*(roj(i,nhp)                            &  
     &                  +  roj(i+1,nhp) )*udsd
      rojmin = fracn0*r0
      if(ufdata(nrecord,i,2) .lt. rojmin*udsd)                           &  
     &    ufdata(nrecord,i,2) = rojmin*udsd
      ufdata(nrecord,i,3) = .25_R8*(pr(i,nhp)                            &  
     &                  +  pr(i+1,nhp) )*udsh/ufdata(nrecord,i,2)
      ufdata(nrecord,i,1) = .25_R8*(pr(i,nhp)                            &  
     &                  +  pr(i+1,nhp) )*udsh/ufdata(nrecord,i,2)
      endif
!
      if(numicsuf.le.0) go to 102
      do 101 k=1,numicsuf
      ispcr = ichargst(k)
      nimp = impnum(k)
      lval = 2
      do 120 lv=3,npsit
      lval = lv
      if(xsv(lval) .gt. ps) go to 130
  120 continue
      val = 0._R8
      go to 131
  130 continue
      lvalm = lval - 1
      v1 = nq(ispcr,nimp,lval)/vp(lval)
      v2 = nq(ispcr,nimp,lvalm)/vp(lvalm)
      f2 = (xsv(lval)-ps)/(xsv(lval)-xsv(lvalm))
      if(f2.le.0) f2 = 0
      if(f2.ge.1._R8) f2 = 1._R8
      val = f2*v2 + (1._R8-f2)*v1
  131 continue
      ufdata(nrecord,i,10+k) = val
  101 continue
  102 continue
!
  100 continue
      nd = 21
      call diagbol(radpow,nd)
      do 200 k=1,nd
      ufdata(nrecord,k,29) = radpow(k)
  200 continue
!
!
      return
!
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

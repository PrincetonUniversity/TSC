      subroutine xptplt
!
!
!.....make a plot of x-point region
!
      USE CLINAM
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER nxpl,i,j,k1,k2
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 cval,xmid,zmid,dxdist,dzdist,ratiox,ratioz,xmin,xmax
      REAL*8 zmin,zmax,pmins,pmaxs,xco,zco
!============
      dimension cval(2)
!============      
!
      if(idiv.eq.0) return
      if(iminsep.ge.imaxsep) return
      if(jminsep.ge.jmaxsep) return
      if(psimin.ge.psilim) return
      if(xsep(1).le.0 .and. xsep(2).le.0) return
      if(lrswtch.gt.0) return
!
      nxpl = nx
      if(acoef(48) .ge. 1) nxpl = int(acoef(48))
!
      xmid = .5_R8*(xary(iminsep)+xary(imaxsep))
      zmid = .5_R8*(zary(jminsep)+zary(jmaxsep))
      dxdist = xary(imaxsep)-xary(iminsep)
      dzdist = zary(jmaxsep)-zary(jminsep)
      ratiox = 1._R8
      ratioz = 1._R8
      if(dxdist.gt.dzdist) ratioz = dxdist/dzdist
      if(dzdist.gt.dxdist) ratiox = dzdist/dxdist
!
      xmin = xmid - 0.5_R8*dxdist*ratiox
      xmax = xmid + 0.5_R8*dxdist*ratiox
      zmin = zmid - 0.5_R8*dzdist*ratioz
      zmax = zmid + 0.5_R8*dzdist*ratioz
      call maps(xmin,xmax,zmin,zmax,.142_R8,.858_R8,.274_R8,1.0_R8)
!
      pmins = psi(iminsep,imaxsep)
      pmaxs = psi(iminsep,imaxsep)
      do 100 i=iminsep,imaxsep
      do 100 j=jminsep,jmaxsep
      pmins = min(psi(i,j),pmins)
      pmaxs = max(psi(i,j),pmaxs)
  100 continue
      if(pmaxs .le. pmins) return
      k1 = 0
      cval(1) = psilim
      cval(2) = (pmaxs-pmins)/nxpl
!
      call rcontr(k1,cval,k2,psi,penx,xary,iminsep,imaxsep,1,            &  
     &                                zary,jminsep,jmaxsep,1 )
      xco = xmin - 0.15_R8*(xmax-xmin)
      zco = zmin
      call setold(xco,zco,1,0,1,1)
      write(s100,1001) kcycle,cval(2)
      call gtext(s100,80,0)
      write(s100,1002) xsep(1),zsep(1),iplim
      call gtext(s100,80,0)
 1002 format("(x,z)=(",1p2e10.3,"), iplim=",                             &  
     &     i4)
      write(nsc1,1000) kcycle
      call frscj(6)
 1001 format(" cycle",i7," pol flux per radian increment",1pe12.4)
 1000 format(" x-point flux,   cycle",i7)
!
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

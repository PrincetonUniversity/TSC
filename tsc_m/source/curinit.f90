      subroutine curinit
!.....5.90 curinit
!
!.....initialize current for initial psi determination
      USE CLINAM
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER i,j,ii,n
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 xzer,ax,zzer,bz,alphax,betaz,t1,t2,aj0,argx,argz,xmid
!============
      xzer = xzeric
      ax = axic
      zzer = zzeric
      bz = bzic
      alphax = 1._R8
      betaz = 0.5_R8
      psilim = 0._R8
      t1 = 4._R8*alphax*ax/(2._R8*alphax + 1._R8)
      t2 = 4._R8*betaz*bz /(2._R8*betaz  + 1._R8)
      aj0 = tcuro/(t1*t2)
      do 30 i=2,nxp
      do 30 j=2,nzp
      ajphi(i,j) = 0._R8
      vd(i,j) = 0._R8
   30 psi(i,j) = 0._R8
      do 40 i=2,nxp
      do 40 j=2,nzp
      if(xary(i) .lt. xzer-ax) go to 40
      if(xary(i) .gt. xzer+ax) go to 40
      if(zary(j) .lt. zzer-bz) go to 40
      if(zary(j) .gt. zzer+bz) go to 40
      argx = ((xary(i)-xzer)/ax)**2
      argz = ((zary(j)-zzer)/bz)**2
      ajphi(i,j) = aj0*(1._R8-argx**alphax)*(1._R8-argz**betaz)
      psi(i,j) = -ajphi(i,j)/(ax*bz)
   40 continue
      do 150 ii=1,nwire
      n = ncoil-nwire+ii
      i = iwire(ii)
      j = jwire(ii)
  150 vd(i,j) = ccoil(n)*xary(i)
      do 500 i=iminn,imaxx
      do 500 j=jminn,jmaxx
      xmid = xary(i)
  500 vd(i,j) = vd(i,j) + xmid*dxdz*ajphi(i,j)
      call icalc
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

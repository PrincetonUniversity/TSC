      function amutlrc (a1,zc1,zd1,rd1,t1,a2,zc2,zd2,rd2,t2,n)
!     mutual of tw0 c0-axial coils
!     summing flux from filaments
!     dl is the size of the subunit for discretizing
!
!===      real mutlrc
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER n,nr1,nz1,nr2,nz2,ir,iz,jr,jz
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 a1,zc1,zd1,rd1,t1,a2,zc2,zd2,rd2,t2,amutlrc,dl,dr1,dz1
      REAL*8 dr2,dz2,sflx,r2,z2,r1,z1,flx
!============
      dl=0.5_R8*sqrt((a1+a2)**2+(zc1+zc2)**2)/n
      nr1=rd1/dl+0.5_R8
      if(nr1.lt.1) nr1=1
      nz1=zd1/dl+0.5_R8
      if(nz1.lt.1) nz1=1
      dr1=rd1/nr1
      dz1=zd1/nz1
!
      nr2=rd2/dl+0.5_R8
      if(nr2.lt.1) nr2=1
      nz2=zd2/dl+0.5_R8
      if(nz2.lt.1) nz2=1
      dr2=rd2/nr2
      dz2=zd2/nz2
!     write(nterm,*) nr1,nz1,nr2,nz2
!     write(nterm,*) rd1,zd1,dl
!
      sflx=0.0_R8
!
      do 10 ir=1,nr2
      r2=(ir-0.5_R8)*dr2-0.5_R8*rd2+a2
      r2=r2*(1._R8+(dr2/r2)**2/24._R8)
      do 10 iz=1,nz2
      z2=(iz-0.5_R8)*dz2-0.5_R8*zd2+zc2
      do 10 jr=1,nr1
      r1=(jr-0.5_R8)*dr1-0.5_R8*rd1+a1
      r1=r1*(1._R8+(dr1/r1)**2/24._R8)
      do 10 jz=1,nz1
      z1=(jz-0.5_R8)*dz1-0.5_R8*zd1+zc1
      call sfilfx(r1,z1,1.0_R8,r2,z2,flx)
      sflx=sflx+flx
!     write(nterm,3) nz,sflx,flx,r1,z1,r2,z21
   10 continue
!     write(nterm,*) 'nzs = ',nzs
!     write(nterm,*) 'nzs = ',nzs
!
      amutlrc=t1*t2*sflx/(nr1*nz1*nr2*nz2)
!
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

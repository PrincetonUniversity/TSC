      function selfrc (a,zd,rd,t,n)
!     self-inductance of coil
!     summing flux from filaments
!     a/n is the subunit for discretizing
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER n,nr,nz2,nz,jp,ja,ia
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 a,zd,rd,t,selfrc,dr,dz,sr,zp,sflx,rp,ra,za,flx,flxp
!============
      REAL*8 mu0
      data mu0 / 12.56637E-7_R8/
!
      nr=rd/(a/n)+0.5_R8
      if(nr.lt.1) nr=1
      nz2=0.5_R8*zd/(a/n)+0.5_R8
      nz=2*nz2
      if(nz.lt.2) nz=2
      dr=rd/nr
      dz=zd/nz
      sr=sqrt(dr*dz/3.14159_R8)
      zp=0.5_R8*(dz-zd)
!
      sflx=0.0_R8
      do 10 jp=1,nr
      rp=(jp-0.5_R8)*dr-0.5_R8*rd+a
      do 10 ja=jp,nr
      ra=(ja-0.5_R8)*dr-0.5_R8*rd+a
 
      do 10 ia=1,nz
      za=(ia-0.5_R8)*dz-0.5_R8*zd
      if(ja.eq.jp.and.ia.eq.1) then
      flx=mu0*rp*(log(8.0_R8*rp/sr)-1.75_R8)
      else
      call sfilfx(ra,za,1.0_R8,rp,zp,flx)
      end if
      if(ia.eq.1) flx=0.5_R8*flx
      if(ja.ne.jp) flx=2.0_R8*flx
      sflx=sflx+2*(nz+1-ia)*flx
      flxp=2*(nz+1-ia)*flx
!     type 2, jp,ja,ia,flxp,ra,za,rp,zp
   10 continue
      selfrc=t**2*sflx/(nr*nz)**2
!
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

!#include "f77_dcomplx.h"
      subroutine shpcntl(ig,nl,alphac,term)
!
      USE CLINAM
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER nl,ig,nxtra,i,j,ncount,ngount,iabs
      INTEGER ii,nsg,k,kk,mn,m,n,ierr,null,ishpwrt
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 alphac,term
      REAL*8 frmshp,sumg,ans,zupp,zlow,ans1,ans2,psicnt,pinterp,gic
      REAL*8 amat1,bvec1,relerr,tau,wmax,wmin,wrat,rsq,rmsnorm
      REAL*8 wrtstp,delimax
      REAL*8 AREAL, sum
!============
!     dimension gmat(50,50),bvec(50),xvec(50),sigma(50),
!    +amat(50,50),xcon(50),zcon(50),delerr(50),
!    +nactc(50),nactg(50),nactgg(50)
!============      
      INTEGER :: istat = 0 
!     REAL*8, ALLOCATABLE, DIMENSION(:,:) :: gmat
!     REAL*8, ALLOCATABLE, DIMENSION(:) :: bvec
!     REAL*8, ALLOCATABLE, DIMENSION(:) :: xvec
!     REAL*8, ALLOCATABLE, DIMENSION(:) :: sigma
!     REAL*8, ALLOCATABLE, DIMENSION(:,:) :: amat
!     REAL*8, ALLOCATABLE, DIMENSION(:) :: xcon
!     REAL*8, ALLOCATABLE, DIMENSION(:) :: zcon
!     REAL*8, ALLOCATABLE, DIMENSION(:) :: delerr
!     INTEGER, ALLOCATABLE, DIMENSION(:) :: nactc
!     INTEGER, ALLOCATABLE, DIMENSION(:) :: nactg
!     INTEGER, ALLOCATABLE, DIMENSION(:) :: nactgg
      REAL*8, DIMENSION(50,50) :: gmat
      REAL*8, DIMENSION(50) :: bvec
      REAL*8, DIMENSION(50) :: xvec
      REAL*8, DIMENSION(50) :: sigma
      REAL*8, DIMENSION(50,50) :: amat
      REAL*8, DIMENSION(50) :: xcon
      REAL*8, DIMENSION(50) :: zcon
      REAL*8, DIMENSION(50) :: delerr
      INTEGER, DIMENSION(50) :: nactc
      INTEGER, DIMENSION(50) :: nactg
      INTEGER, DIMENSION(50) :: nactgg
!============      
!     IF(.not.ALLOCATED(gmat)) ALLOCATE( gmat(50,50), STAT=istat)
!     IF(.not.ALLOCATED(bvec)) ALLOCATE( bvec(50), STAT=istat)
!     IF(.not.ALLOCATED(xvec)) ALLOCATE( xvec(50), STAT=istat)
!     IF(.not.ALLOCATED(sigma)) ALLOCATE( sigma(50), STAT=istat)
!     IF(.not.ALLOCATED(amat)) ALLOCATE( amat(50,50), STAT=istat)
!     IF(.not.ALLOCATED(xcon)) ALLOCATE( xcon(50), STAT=istat)
!     IF(.not.ALLOCATED(zcon)) ALLOCATE( zcon(50), STAT=istat)
!     IF(.not.ALLOCATED(delerr)) ALLOCATE( delerr(50), STAT=istat)
!     IF(.not.ALLOCATED(nactc)) ALLOCATE( nactc(50), STAT=istat)
!     IF(.not.ALLOCATED(nactg)) ALLOCATE( nactg(50), STAT=istat)
!     IF(.not.ALLOCATED(nactgg)) ALLOCATE( nactgg(50), STAT=istat)
!============      
      if (istat .ne. 0) stop 'Allocation Error : shpcntl  ' 
!============      
!
      if(ncnt .le. 0) return
      term=0.0_R8
      nxtra=0
!
!...feedback system turned on or off
!
      if(tfbon(nl) .ge. 0 .and. time .lt. tfbon(nl)) return
      if(tfbons(nl) .lt. 0 .and. kcycle .lt. int(-tfbons(nl))) return
      if(tfbof(nl) .ge. 0 .and. time .gt. tfbof(nl)) return
      if(tfbofs(nl) .lt. 0 .and. kcycle .gt. int(-tfbofs(nl))) return
!
!...factor for starting up the shape control system
!
      frmshp=1.0_R8
      if((times-tfbons(nl)) .le. acoef(512)) then
      frmshp=(times-tfbons(nl))/acoef(512)
!     alphac=alphac*acoef(512)/(times-tfbons(nl))
      endif
!
!...calculate present value of the constraint points for shape
!...control
!
      do 801 i=1,ncnt
      xcon(i)=0.0_R8
      zcon(i)=0.0_R8
      do 802 j=1,ntpts
      xcon(i)=xcon(i)+fact(j)*xcon0(j,i)
      zcon(i)=zcon(i)+fact(j)*zcon0(j,i)
 802  continue
 801  continue
!
!...find coils and wires active in shape control, and their
!...home group numbers
!
      ncount=0
      ngount=0
!
      do 700 i=1,ncoil-nwire
      if(ngrvc1(i) .eq. ig .or. ngrvc2(i) .eq. ig .or. ngrvc3(i)         &  
     &.eq. ig .or. ngrvc4(i) .eq. ig .or. iabs(igroupc(i))               &  
     &.eq. ig) then
      ncount=ncount+1
      nactc(ncount)=i
      nactg(ncount)=iabs(igroupc(i))
      endif
700   continue
      do 701 i=ncoil-nwire+1,ncoil
      ii=i-(ncoil-nwire)
      if(ngrvw1(ii) .eq. ig .or. ngrvw2(ii) .eq. ig .or. ngrvw3(ii)      &  
     &.eq. ig .or. ngrvw4(ii) .eq. ig .or. iabs(igroupw(ii))             &  
     &.eq. ig) then
      ncount=ncount+1
      nactc(ncount)=i
      nactg(ncount)=iabs(igroupw(ii))
      endif
701   continue
      do 702 i=1,pngroup
      do 703 j=1,ncount
      if(nactg(j) .eq. i) then
      ngount=ngount+1
      nactgg(ngount)=nactg(j)
      go to 702
      endif
703   continue
702   continue
!
!...calculate green's functions for coils and wires to constraint
!...points, and collapse those whose home group numbers are
!...the same
!
      do 100 i=1,ncnt
      do 101 j=1,ngount
      sumg=0.0_R8
      nsg=0
      do 102 k=1,ncount
      kk=nactc(k)
      if(nactg(k) .eq. nactgg(j)) then
      call gf(ineg,0,xcon(i),zcon(i),xcoil(kk),zcoil(kk),ans)
      sumg=sumg+ans/tpi/udsi
      nsg=nsg+1
      if(isym .eq. 1 .and. zcoil(kk) .ne. 0) then
      call gf(ineg,0,xcon(i),zcon(i),xcoil(kk),-zcoil(kk),ans)
      sumg=sumg+ans/tpi/udsi
      endif
      endif
102   continue
      gmat(j,i)=sumg/AREAL(nsg)
101   continue
100   continue
!
!...calculate green's functions for coils and wires to plasma
!...magnetic axis, and collapse those whose home group numbers are
!...the same
!
      nxtra=nxtra+1
      do 900 j=1,ngount
      sumg=0.0_R8
      nsg=0
      do 901 k=1,ncount
      kk=nactc(k)
      if(nactg(k) .eq. nactgg(j)) then
      call gf(ineg,0,xmag,zmag,xcoil(kk),zcoil(kk),ans)
      sumg=sumg+ans/tpi/udsi
      nsg=nsg+1
      if(isym .eq. 1 .and. zcoil(kk) .ne. 0) then
      call gf(ineg,0,xmag,zmag,xcoil(kk),-zcoil(kk),ans)
      sumg=sumg+ans/tpi/udsi
      endif
      endif
901   continue
      gmat(j,ncnt+nxtra)=sumg/AREAL(nsg)
900   continue
!
!...calculate  poloidal flux for coils and wires above
!...and below the magnetic axis for asymmetric plasmas,
!...and collapse those whose home group numbers are the same
!
      if(isym .eq. 0) then
      nxtra=nxtra+1
      do 1100 j=1,ngount
      sumg=0.0_R8
      nsg=0
      do 1101 k=1,ncount
      kk=nactc(k)
      if(nactg(k) .eq. nactgg(j)) then
      zupp=zmag+0.5_R8
      zlow=zmag-0.5_R8
      call gf(ineg,0,xmag,zupp,xcoil(kk),zcoil(kk),ans1)
      call gf(ineg,0,xmag,zlow,xcoil(kk),zcoil(kk),ans2)
      sumg=sumg+(ans1-ans2)/tpi/udsi
      nsg=nsg+1
      endif
1101  continue
      gmat(j,ncnt+nxtra)=sumg/AREAL(nsg)
1100  continue
      endif
!
!...setup matrix terms for shape control
!
      do 120 i=1,ncnt+nxtra
      do 121 j=1,ngount
      amat(i,j)=gmat(j,i)
121   continue
!
!...setup RHS flux errors for shape control
!
!...special TPX control line
!     if(i .eq. 1) go to 1124
      if(i .eq. ncnt+1) go to 1121
      if(i .eq. ncnt+2) go to 1122
      if(i .eq. ncnt+3) go to 1123
      psicnt=pinterp(xcon(i),zcon(i),0,0)
      shpdif(i)=((psisep-psicnt)-shpdifo(i))/dts
      shpsum(i)=shpsum(i)+(psisep-psicnt)*dts
      bvec(i)=acoef(513)*(psisep-psicnt)+acoef(511)*shpsum(i)            &  
     &+acoef(514)*shpdif(i)
      shpdifo(i)=(psisep-psicnt)
      go to 120
1121  continue
      bvec(i)=0.0_R8
      go to 120
1122  continue
      bvec(i)=0.0_R8
      go to 120
1123  continue
      bvec(i)=0.0_R8
      go to 120
1124  continue
!
!...special TPX control RHS definition
      call gf(ineg,0,xcon(1),zcon(1),xcoil(8),zcoil(8),ans)
      gic=ans/tpi/udsi
      if(isym .eq. 1 .and. zcoil(8) .ne. 0.0_R8) then
      call gf(ineg,0,xcon(1),zcon(1),xcoil(8),-zcoil(8),ans)
      gic=gic+ans/tpi/udsi
      endif
      shpsum(i)=shpsum(i)+(gic*ccoil(8)*udsi)*dts
      bvec(i)=gic*ccoil(8)*udsi+acoef(511)*shpsum(i)
      go to 120
120   continue
!
!...append regularization terms to matrix
!
      do 122 i=1,ngount
      do 123 j=1,ngount
      if(i .eq. j) then
      amat(i+ncnt+nxtra,j)=alphac
      else
      amat(i+ncnt+nxtra,j)=0.0_R8
      endif
123   continue
122   continue
!
!...time delay option
!
      if(idelay(nl) .le. 0) go to 500
      indxd1(nl)=indxd1(nl)+1
      if(indxd1(nl) .gt. idelay(nl)) indxd1(nl)=indxd1(nl)-idelay(nl)
      do 150 i=1,ncnt+ngount+nxtra
      do 151 j=1,ngount
      amat1=amats(i,j,indxd1(nl))
      amats(i,j,indxd1(nl))=amat(i,j)
      amat(i,j)=amat1
151   continue
      bvec1=bvecs(i,indxd1(nl))
      bvecs(i,indxd1(nl))=bvec(i)
      bvec(i)=bvec1
150   continue
      if(kcycle .lt. idelay(nl)) return
500   continue
!
!...invert matrix to solve for the active coil and wire current
!...incremental changes
!
      mn=50
      m=ncnt+ngount+nxtra
      n=ngount
      if( m .gt. mn .or. n .gt. mn) then
      write(6,*) " error in shpcntl, m,n>mn: ",m,n,mn
      ineg=28
      return
      endif
      ierr=0
      relerr=0.0_R8
      call svd2(mn,m,n,amat,bvec,sigma,1,1,xvec,ierr,relerr,             &  
     &tau,wmax,wmin,wrat,rsq,rmsnorm,null)
      if(ierr .eq. 0) go to 600
      write(nout,601)
601   format("  error return from svd2 ")
      ineg=28
      return
600   continue
!
      do 1109 i=1,ncnt+nxtra
      sum=0.0_R8
      do 1102 j=1,ngount
      sum=sum+amat(i,j)*xvec(j)
1102  continue
      delerr(i)=bvec(i)-sum
1109  continue
!
!...periodic diagnostic output
!
      if(wrtstp .eq. 10) then
      write(nout,1105) frmshp,alphac                                       
1105  format(/,'***shape match***',5x,'frmshp = ',e10.3,                 &
     &5x,'alphac = ',e10.3,/)
      write(nout,1103) (bvec(j),j=1,ncnt+nxtra)
      write(nout,1103) (delerr(j),j=1,ncnt+nxtra)
      write(nout,1103) (xvec(j),j=1,ngount)
1103  format(7e10.3)
      wrtstp=0
      endif
      wrtstp=wrtstp+1
!
!...one time diagnostic output when control started or code restarted
!
      if(ishpwrt .ge. 1) go to 1002
      ishpwrt=1
      write(nout,1050)
1050  format(/'***shape control information***',/)
      write(nout,1049) ncnt,nxtra,ngount,m,n
1049  format(5i5,/)
      write(nout,1048) xmag,zmag
1048  format(2e10.3,/)
      write(nout,1047) (xcon(j),zcon(j),j=1,ncnt)
1047  format(2e10.3)
      write(nout,1043)
      write(nout,1046) ((amat(i,j),j=1,ngount),i=1,ncnt+nxtra)
1046  format(7e10.3)
      write(nout,1043)
      write(nout,1103) (bvec(j),j=1,ncnt+nxtra)
      write(nout,1103) (delerr(j),j=1,ncnt+nxtra)
      write(nout,1103) (xvec(j),j=1,ngount)
      write(nout,1043)
      write(nout,1044) psisep,psimin,alphac
1044  format(3e10.3)
      write(nout,1043)
      write(nout,1103) (sigma(j),j=1,ngount)
      write(nout,1051) wmax,wmin,wrat,relerr,tau
1051  format(5e10.3,/)
1043  format(/)
1002  continue
!
!...assignment of currents to coils
!
      delimax=-1.0E30_R8
      do 130 i=1,ngount
      if(abs(xvec(i)) .gt. delimax) delimax=abs(xvec(i))
130   continue
      term=frmshp*delimax
      do 131 i=1,ngount
      do 132 j=1,ncount
      if(nactg(j) .eq. nactgg(i)) then
      fadj(ig,nactc(j))=xvec(i)/delimax
      endif
132   continue
131   continue
!
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
! 15Apr2005 fgtok

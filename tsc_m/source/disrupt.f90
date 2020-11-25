      subroutine disrupt
!
      USE CLINAM
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER idstart,i,j,iw,jw,mpolc,n
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 position
      REAL*8 g1,g2,cppr1,cppr2,cppz1,cppz2,tcur,tcurmain,diamag,fac
      REAL*8 g5,g6,g8,g9,tcurhalo
!============
!     dimension bpolr(pnwire),bpolz(pnwire),cpolr(2*pnwire),
!    +btoro(2*pnwire),cpolpr(pnwire),cpolpz(pnwire),vecx(penx,penz),
!    +vecy(penx,penz),vecz(penx,penz),itap1(penx,penz),cpolz(2*pnwire),
!    +iptw(pnwire,4),jptw(pnwire,4),iwtw(pnwire,2),jwtw(pnwire,2),
!    +xpolc(2*pnwire),zpolc(2*pnwire),forcr(pnwire),forcz(pnwire),
!    +forcpr(2*pnwire),forcpz(2*pnwire)
!============      
      INTEGER :: istat = 0 
      REAL*8, ALLOCATABLE, DIMENSION(:) :: bpolr
      REAL*8, ALLOCATABLE, DIMENSION(:) :: bpolz
      REAL*8, ALLOCATABLE, DIMENSION(:) :: cpolr
      REAL*8, ALLOCATABLE, DIMENSION(:) :: btoro
      REAL*8, ALLOCATABLE, DIMENSION(:) :: cpolpr
      REAL*8, ALLOCATABLE, DIMENSION(:) :: cpolpz
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: vecx
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: vecy
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: vecz
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: itap1
      REAL*8, ALLOCATABLE, DIMENSION(:) :: cpolz
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: iptw
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: jptw
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: iwtw
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: jwtw
      REAL*8, ALLOCATABLE, DIMENSION(:) :: xpolc
      REAL*8, ALLOCATABLE, DIMENSION(:) :: zpolc
      REAL*8, ALLOCATABLE, DIMENSION(:) :: forcr
      REAL*8, ALLOCATABLE, DIMENSION(:) :: forcz
      REAL*8, ALLOCATABLE, DIMENSION(:) :: forcpr
      REAL*8, ALLOCATABLE, DIMENSION(:) :: forcpz
!============      
      IF(.not.ALLOCATED(bpolr)) ALLOCATE( bpolr(pnwire), STAT=istat)
      IF(.not.ALLOCATED(bpolz)) ALLOCATE( bpolz(pnwire), STAT=istat)
      IF(.not.ALLOCATED(cpolr)) ALLOCATE( cpolr(2*pnwire), STAT=istat)
      IF(.not.ALLOCATED(btoro)) ALLOCATE( btoro(2*pnwire), STAT=istat)
      IF(.not.ALLOCATED(cpolpr)) ALLOCATE( cpolpr(pnwire), STAT=istat)
      IF(.not.ALLOCATED(cpolpz)) ALLOCATE( cpolpz(pnwire), STAT=istat)
      IF(.not.ALLOCATED(vecx)) ALLOCATE( vecx(penx,penz), STAT=istat)
      IF(.not.ALLOCATED(vecy)) ALLOCATE( vecy(penx,penz), STAT=istat)
      IF(.not.ALLOCATED(vecz)) ALLOCATE( vecz(penx,penz), STAT=istat)
      IF(.not.ALLOCATED(itap1)) ALLOCATE( itap1(penx,penz), STAT=istat)
      IF(.not.ALLOCATED(cpolz)) ALLOCATE( cpolz(2*pnwire), STAT=istat)
      IF(.not.ALLOCATED(iptw)) ALLOCATE( iptw(pnwire,4), STAT=istat)
      IF(.not.ALLOCATED(jptw)) ALLOCATE( jptw(pnwire,4), STAT=istat)
      IF(.not.ALLOCATED(iwtw)) ALLOCATE( iwtw(pnwire,2), STAT=istat)
      IF(.not.ALLOCATED(jwtw)) ALLOCATE( jwtw(pnwire,2), STAT=istat)
      IF(.not.ALLOCATED(xpolc)) ALLOCATE( xpolc(2*pnwire), STAT=istat)
      IF(.not.ALLOCATED(zpolc)) ALLOCATE( zpolc(2*pnwire), STAT=istat)
      IF(.not.ALLOCATED(forcr)) ALLOCATE( forcr(pnwire), STAT=istat)
      IF(.not.ALLOCATED(forcz)) ALLOCATE( forcz(pnwire), STAT=istat)
      IF(.not.ALLOCATED(forcpr)) ALLOCATE( forcpr(2*pnwire), STAT=istat)
      IF(.not.ALLOCATED(forcpz)) ALLOCATE( forcpz(2*pnwire), STAT=istat)
!============      
      if (istat .ne. 0) stop 'Allocation Error : disrupt  ' 
!============      
!
      if(idstart .gt. 0) go to 1000
      if( numargs .lt. 1 ) then
         filename ='displas.dat'
      else
         filename ='displas.dat' // '.' // trim(suffix)
      end if
      open(39,file=trim(filename),status='unknown',position='append')
      if( numargs .lt. 1 ) then
         filename ='diswirt.dat'
      else
         filename ='diswirt.dat' // '.' // trim(suffix)
      end if
      open(40,file=trim(filename),status='unknown',position='append')
      if( numargs .lt. 1 ) then
         filename ='discoil.dat'
      else
         filename ='discoil.dat' // '.' // trim(suffix)
      end if
      open(41,file=trim(filename),status='unknown',position='append')
      if( numargs .lt. 1 ) then
         filename ='disjayp.dat'
      else
         filename ='disjayp.dat' // '.' // trim(suffix)
      end if 
      open(42,file=trim(filename),status='unknown',position='append')
      if( numargs .lt. 1 ) then
         filename ='diswirp.dat'
      else
         filename ='diswirp.dat' // '.' // trim(suffix)
      end if
      open(43,file=trim(filename),status='unknown',position='append')
!
      do 89 i=1,nxp
      do 89 j=1,nzp
  89  itap1(i,j)=0
!
      do 90 i=1,nwire
      iw=iwire(i)
      jw=jwire(i)
      itap1(iw,jw)=1
  90  continue
!
      do 92 i=1,nwire
      do 93 j=1,4
      iptw(i,j)=0
      jptw(i,j)=0
  93  continue
  92  continue
!
      do 91 i=1,nwire
      iw=iwire(i)
      jw=jwire(i)
      if(itap1(iw+1,jw) .eq. 1) then
      iwtw(i,1)=iw+1
      jwtw(i,1)=jw
      endif
      if(itap1(iw,jw+1) .eq. 1) then
      iwtw(i,2)=iw
      jwtw(i,2)=jw+1
      endif
      if(itap1(iw+1,jw) .ne. 1 .and. iexv(iw+1,jw) .eq. 0) then
      iptw(i,1)=iw+1
      jptw(i,1)=jw
      endif
      if(itap1(iw,jw+1) .ne. 1 .and. iexv(iw,jw+1) .eq. 0) then
      iptw(i,2)=iw
      jptw(i,2)=jw+1
      endif
      if(itap1(iw-1,jw) .ne. 1 .and. iexv(iw-1,jw) .eq. 0) then
      iptw(i,3)=iw-1
      jptw(i,3)=jw
      endif
      if(itap1(iw,jw-1) .ne. 1 .and. iexv(iw,jw-1) .eq. 0) then
      iptw(i,4)=iw
      jptw(i,4)=jw-1
      endif
  91  continue
!
      idstart=1
 1000 continue
!
      mpolc=0
      do 100 i=1,nwire
      iw=iwire(i)
      jw=jwire(i)
      n=ncoil-nwire+i
      cwires(i)=ccoil(n)*udsi
      bpolr(i)=(psi(iw,jw+1)-psi(iw,jw-1))/(2._R8*deez*xwire(i))
      bpolz(i)=-(psi(iw+1,jw)-psi(iw-1,jw))/(2._R8*deex*xwire(i))
      forcr(i)=cwires(i)*bpolz(i)*tpi*xwire(i)
      forcz(i)=-cwires(i)*bpolr(i)*tpi*xwire(i)
!
      if(iwtw(i,1) .ne. 0) then
      mpolc=mpolc+1
      xpolc(mpolc)=(xary(iw+1)+xary(iw))/2._R8
      zpolc(mpolc)=zary(jw)
      g1=g(iw+1,jw)*xsqoj(iw+1)
      g2=g(iw+1,jw+1)*xsqoj(iw+1)
      cpolr(mpolc)=-tpi*udsi*(g2-g1)
      cpolz(mpolc)=0._R8
      btoro(mpolc)=0.5_R8*(g1+g2)/xpolc(mpolc)
      forcpz(mpolc)=cpolr(mpolc)*btoro(mpolc)*deex
      endif
      if(iwtw(i,2) .ne. 0) then
      mpolc=mpolc+1
      xpolc(mpolc)=xary(iw)
      zpolc(mpolc)=(zary(jw+1)+zary(jw))/2._R8
      g1=g(iw,jw+1)*xsqoj(iw)
      g2=g(iw+1,jw+1)*xsqoj(iw+1)
      cpolz(mpolc)=tpi*udsi*(g2-g1)
      cpolr(mpolc)=0._R8
      btoro(mpolc)=0.5_R8*(g1+g2)/xpolc(mpolc)
      forcpr(mpolc)=-cpolz(mpolc)*btoro(mpolc)*deez
      endif
!
      cppr1=0._R8
      cppr2=0._R8
      cppz1=0._R8
      cppz2=0._R8
      if(iptw(i,1) .ne. 0) then
      g1=g(iw+1,jw)*xsqoj(iw+1)
      g2=g(iw+1,jw+1)*xsqoj(iw+1)
      cppr1=-tpi*udsi*(g2-g1)
      endif
      if(iptw(i,2) .ne. 0) then
      g1=g(iw,jw+1)*xsqoj(iw)
      g2=g(iw+1,jw+1)*xsqoj(iw+1)
      cppz1=tpi*udsi*(g2-g1)
      endif
      if(iptw(i,3) .ne. 0) then
      g1=g(iw,jw)*xsqoj(iw)
      g2=g(iw,jw+1)*xsqoj(iw)
      cppr2=-tpi*udsi*(g2-g1)
      endif
      if(iptw(i,4) .ne. 0) then
      g1=g(iw,jw)*xsqoj(iw)
      g2=g(iw+1,jw)*xsqoj(iw+1)
      cppz2=tpi*udsi*(g2-g1)
      endif
      cpolpr(i)=cppr1+cppr2
      cpolpz(i)=cppz1+cppz2
 100  continue
!
      tcur=0._R8
      tcurmain=0._R8
      diamag=0._R8
      fac=1._R8
      if(isym .eq. 1) fac=2._R8
      do 200 i=2,nxp
      do 201 j=2,nzp
      if(iexv(i,j) .ne. 0) go to 201
      vecx(i,j)=ajphi(i,j)*udsi*dxdz
      g5=g(i,j)*xsqoj(i)
      g6=g(i+1,j)*xsqoj(i+1)
      g8=g(i,j+1)*xsqoj(i)
      g9=g(i+1,j+1)*xsqoj(i+1)
      vecy(i,j)=-tpi*udsi*0.5_R8*(g8+g9-g6-g5)
      vecz(i,j)=tpi*udsi*0.5_R8*(g9+g6-g8-g5)
      tcur=tcur+fac*ajphi(i,j)*udsi*dxdz
      if(psi(i,j) .le. psilim) then
      tcurmain=tcurmain+fac*ajphi(i,j)*udsi*dxdz
      diamag=diamag+fac*(g(i,j)-gzero/xsqoj(i))
      endif
 201  continue
 200  continue
      tcurhalo=tcur-tcurmain
!
      write(39,599)
 599  format('         TIME           IP       IPMAIN       IPHALO',     &  
     &'         RMAG         ZMAG       DIAMAG           LI',            &  
     &'          WTH        THALO        WHALO')
      write(39,600) times,tcur,tcurmain,tcurhalo,xmag,zmag,diamag,       &  
     &2._R8*ali2,uint,thalos,whalos
 600  format(11(1pe13.5))
!
      write(40,610)
 610  format('         TIME')
      write(40,611) times
 611  format(1pe13.5)
      write(40,601)
 601  format('GRP            R            Z         ITOR',               &  
     &'        BPOLR        BPOLZ        FORCR        FORCZ',            &  
     &'       IPPOLR       IPPOLZ')
      do 801 j=1,pngroup
      do 800 i=1,nwire
      if(igroupw(i) .eq. j) then
      write(40,602) igroupw(i),xwire(i),zwire(i),cwires(i),bpolr(i),     &  
     &bpolz(i),forcr(i),forcz(i),cpolpr(i),cpolpz(i)
      endif
 800  continue
 801  continue
 602  format(i3,9(1pe13.5))
!
      write(41,610)
      write(41,611) times
      write(41,604)
 604  format('GRP            R            Z           DR           DZ',  &  
     &'            I')
      do 850 i=1,ncoil-nwire
      write(41,605) igroupc(i),xcoil(i),zcoil(i),dxcoil(i),dzcoil(i),    &  
     &ccoils(i)
 850  continue
 605  format(i3,5(1pe13.5))
!
!
      write(42,610)
      write(42,611) times
      write(42,612)
 612  format('            R            Z            I')
      do 300 i=2,nxp
      do 301 j=2,nzp
      if(iexv(i,j) .ne. 0) go to 301
      write(42,500) xary(i),zary(j),vecx(i,j)
 301  continue
 300  continue
 500  format(3(1pe13.5))
!
      write(43,610)
      write(43,611) times
      write(43,613)
 613  format('            R            Z        IPOLR       IPOLZ',      &  
     &'         BTOR        FORCR        FORCZ')
      do 410 i=1,mpolc
      write(43,510) xpolc(i),zpolc(i),cpolr(i),cpolz(i),btoro(i),        &  
     &forcpr(i),forcpz(i)
 410  continue
 510  format(5(1pe13.5))
!
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

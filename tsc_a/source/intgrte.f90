      subroutine intgrte
!
!.....10.2 avrge
!
      USE CLINAM
      USE BALCLI
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER iwr,istep,nord,nordx,mh,i,n,j,iflag,lmax
      INTEGER int
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 epsx,thpi,range,thmax,h,thet1,h0,thet2,tdpi
!============
!     dimension y(2*ppsi),yp(2*ppsi),y1past(ppsi)
      data iwr / 0 /
!============      
      INTEGER :: istat = 0 
      REAL*8, ALLOCATABLE, DIMENSION(:) :: y
      REAL*8, ALLOCATABLE, DIMENSION(:) :: yp
      REAL*8, ALLOCATABLE, DIMENSION(:) :: y1past
!============      
      IF(.not.ALLOCATED(y)) ALLOCATE( y(2*ppsi), STAT=istat)
      IF(.not.ALLOCATED(yp)) ALLOCATE( yp(2*ppsi), STAT=istat)
      IF(.not.ALLOCATED(y1past)) ALLOCATE( y1past(ppsi), STAT=istat)
!============      
      if (istat .ne. 0) stop 'Allocation Error : intgrte  ' 
!============      
!
      epsx = 1.0E-4_R8
      thpi = 10.0_R8
!
      range = thpi * pi
      istep = ( range + .0001_R8) / dthe + 1
      thmax = range
!
      h = dthe
      nord = 2*(npsit-4)
      nordx = nord
      mh = 1
!
!..... initialize the integration.
!
      do 5 i = 1, nord
      yp(i) = 0.0_R8
    5 continue
!
      thet1 = 0.0_R8
      do 6 n=1,nord,2
      j = (n+1)/2 + 2
      node1(j) = 0
      y(n) = 1._R8
      y(n+1) = 0._R8
      y1past(j) = 1._R8
    6 continue
!
!
      h0 = h
      iflag = 1
!
      do 100 int = 1, istep
!
      lmax = 0
      if ( range-thet1 .le. h0 ) lmax = 1
!
      if ( lmax .eq. 1 ) h0 = range - thet1
            thet2 = thet1 + h0
      call de(nord,y,thet1,thet2,epsx,epsx,iflag)
      if(iflag.ne.2) go to 101
      if(ifcount.gt.50*thpi*(1+isym)*nthe) go to 101
!
!..... count and place last node.
!
      do 20 n=1,nord,2
      j = (n+1)/2+2
      if(y1past(j)*y(n) .gt. 0.0_R8) go to 20
      node1(j) = node1(j) + 1
!
   20 continue
!
!
      do 21 n=1,nord,2
      j = (n+1)/2+2
      y1past(j) = y(n)
   21 continue
      if ( lmax .eq. 1 ) go to 120
!
  100 continue
  101 continue
!
      tdpi = thet1/pi
      write (nout,110 ) iflag, tdpi,lmax,int,thet1,ifcount,kcycle
  110 format ( 1x, "warning from balloon:", i4, 2x, "tdpi = ",1pe12.4 ,  &  
     &      "lmax = ",i3,"  int,thet1=",i3,1pe12.4,                      &  
     &      "  ifcount = ",i10, "  kcycle = ", i10)
      ineg=31
!
  120 continue
!
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

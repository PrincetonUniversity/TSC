      subroutine de(neqn,y,t,tout,relerr,abserr,iflag)
!
      USE PARAM
      USE DECOM

      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER iflag,neqn,maxnum,isn,iabs,nostep,kle4,l,kold,k
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 y,t,tout,relerr,abserr,psi
      REAL*8 fouru,eps,del,absdel,tend,releps,abseps,x,h,hold
!============
      logical start,crash,stiff
! shampine-gordon code
!
!     common/decom/ delsgn, told, isnold
      dimension y(neqn),psi(12)
!     dimension yy(2*ppsi),wt(2*ppsi),phi(2*ppsi,16),p(2*ppsi),
!    1          yp(2*ppsi),ypout(2*ppsi)
!
      data fouru/2.84E-14_R8/
      data maxnum/500/
!============      
      INTEGER :: istat = 0 
      REAL*8, ALLOCATABLE, DIMENSION(:) :: yy
      REAL*8, ALLOCATABLE, DIMENSION(:) :: wt
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: phi
      REAL*8, ALLOCATABLE, DIMENSION(:) :: p
      REAL*8, ALLOCATABLE, DIMENSION(:) :: yp
      REAL*8, ALLOCATABLE, DIMENSION(:) :: ypout
!============      
      IF(.not.ALLOCATED(yy)) ALLOCATE( yy(2*ppsi), STAT=istat)
      IF(.not.ALLOCATED(wt)) ALLOCATE( wt(2*ppsi), STAT=istat)
      IF(.not.ALLOCATED(phi)) ALLOCATE( phi(2*ppsi,16), STAT=istat)
      IF(.not.ALLOCATED(p)) ALLOCATE( p(2*ppsi), STAT=istat)
      IF(.not.ALLOCATED(yp)) ALLOCATE( yp(2*ppsi), STAT=istat)
      IF(.not.ALLOCATED(ypout)) ALLOCATE( ypout(2*ppsi), STAT=istat)
!============      
      if (istat .ne. 0) stop 'Allocation Error : de  ' 
!============      
      if(neqn.lt.1.or.neqn.gt.2*ppsi) goto 10
      if(t.eq.tout) goto 10
      if(relerr.lt.0._R8.or.abserr.lt.0._R8) goto 10
      eps=max(relerr,abserr)
      if(eps.lt.0._R8) goto 10
      if(iflag.eq.0) goto 10
      isn=isign(1,iflag)
      iflag=iabs(iflag)
      if(iflag.eq.1) goto 20
      if(t.ne.told) goto 10
      if(iflag.ge.2.and.iflag.le.5) goto 20
10    iflag=6
      return
20    del=tout-t
      absdel=abs(del)
      tend=t+10.0_R8*del
      if(isn.lt.0) tend=tout
      nostep=0
      kle4=0
      stiff=.false.
      releps=relerr/eps
      abseps=abserr/eps
      if(iflag.eq.1) goto 30
      if(isnold.lt.0) goto 30
      if(delsgn*del.gt.0._R8) goto 50
30    start=.true.
      x=t
      do 40 l=1,neqn
40    yy(l)=y(l)
      delsgn=sign(1.0_R8,del)
      h=sign(max(abs(tout-x),fouru*abs(x)),tout-x)
50    if(abs(x-t).lt.absdel) goto 60
      call intrp(x,yy,tout,y,ypout,neqn,kold,phi,psi)
      iflag=2
      t=tout
      told=t
      isnold=isn
      return
60    if(isn.gt.0.or.abs(tout-x).ge.fouru*abs(x)) goto 80
      h=tout-x
      call der(x,yy,yp,neqn)
      do 70 l=1,neqn
70    y(l)=yy(l)+h*yp(l)
      iflag=2
      t=tout
      told=t
      isnold=isn
      return
80    if(nostep.lt.maxnum) goto 100
      iflag=isn*4
      if(stiff) iflag=isn*5
      do 90 l=1,neqn
90    y(l)=yy(l)
      t=x
      told=t
      isnold=1
      return
100   h=sign(min(abs(h),abs(tend-x)),h)
      do 110 l=1,neqn
110   wt(l)=releps*abs(yy(l))+abseps
      call step(x,yy,neqn,h,eps,wt,start,                                &  
     &hold,k,kold,crash,phi,p,yp,psi)
      if(.not.crash) goto 130
      iflag=isn*3
      relerr=eps*releps
      abserr=eps*abseps
      do 120 l=1,neqn
120   y(l)=yy(l)
      t=x
      told=t
      isnold=1
      return
130   nostep=nostep+1
      kle4=kle4+1
      if(kold.gt.4) kle4=0
      if(kle4.ge.50) stiff=.true.
      goto 50
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

      subroutine volt2(v1,v2,v3,v4,v5,v6,v7)
!...This subroutine calculates the voltages on the PF coils from
!...the linear gain law using estimates of the state vector (currents).
!...The estimator is advanced in time driven by the gap distances
!...and their derivatives, from GA/LLNL group
!
      USE CLINAM
      USE CINTG

      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER ihumsrt,ios,i,j,neqn,ifail
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 v2,v3,v4,v5,v6,v7,v1
      REAL*8 work,tstart,tend,tol,fcnh
!============
      external fcnh
!     common/cintg/ahat(201,201),bhat(201,7),chat(7,201),fhat(201,7),
!    +gapdw(6),dgapdw(6)
      dimension work(201,7)
!============      
!
      if(ihumsrt .gt. 0) go to 1111
      if( numargs .lt. 1 ) then
         filename = 'ahatmat'
      else
         filename = 'ahatmat' // '.' // trim(suffix)
      end if
      open(67,file=trim(filename),status='unknown',iostat=ios)
      read(67,1112) ((ahat(i,j),i=1,201),j=1,201)
 1112 format(e17.10)
      close(67)
      if( numargs .lt. 1 ) then
         filename = 'bhatmat'
      else
         filename = 'bhatmat' // '.' // trim(suffix)
      end if
      open(68,file=trim(filename),status='unknown',iostat=ios)
      read(68,1113) ((bhat(i,j),i=1,201),j=1,6)
 1113 format(e17.10)
      close(68)
      if( numargs .lt. 1 ) then
         filename = 'chatmat'
      else
         filename = 'chatmat' // '.' // trim(suffix)
      end if
      open(69,file=trim(filename),status='unknown',iostat=ios)
      read(69,1114) ((chat(i,j),i=1,6),j=1,201)
 1114 format(e17.10)
      close(69)
      if( numargs .lt. 1 ) then
         filename = 'fhatmat'
      else
         filename = 'fhatmat' // '.' // trim(suffix)
      end if
      open(70,file=trim(filename),status='unknown',iostat=ios)
      read(70,1115) ((fhat(i,j),i=1,201),j=1,6)
 1115 format(e17.10)
      close(70)
      ihumsrt=1
      if( numargs .lt. 1 ) then
         filename = 'gapsout'
      else
         filename = 'gapsout' // '.' // trim(suffix)
      end if
      open(61,file=trim(filename),status='unknown',iostat=ios)
 1111 continue
!
!...construct estimated values for the state vector
!
      call gaps2(gapdw(1),gapdw(2),gapdw(3),gapdw(4),gapdw(5),           &  
     &gapdw(6),dgapdw(1),dgapdw(2),dgapdw(3),dgapdw(4),dgapdw(5),        &  
     &dgapdw(6))
      tstart=times-dts
      tend=times
      neqn=201
      tol=1.0E-6_R8
      ifail=0
      call d02bae(tstart,tend,neqn,xveh,tol,fcnh,work,ifail)
      if(ifail .eq. 0) go to 600
      write(nout,601)
601   format("  error return from d02bae ")
      ineg=28
      return
600   continue
!
!...construct coil voltages from feedback law
!
      v1=0.0_R8
      v2=0.0_R8
      v3=0.0_R8
      v4=0.0_R8
      v5=0.0_R8
      v6=0.0_R8
      v7=0.0_R8
      do 120 i=1,201
      v2=v2+chat(1,i)*xveh(i)
      v3=v3+chat(2,i)*xveh(i)
      v4=v4+chat(3,i)*xveh(i)
      v5=v5+chat(4,i)*xveh(i)
      v6=v6+chat(5,i)*xveh(i)
      v7=v7+chat(6,i)*xveh(i)
 120  continue
!
      write(61,1001) times,gapdw(1),gapdw(2),gapdw(3),gapdw(4),gapdw(5),  &  
!    &                                                                   &  
     &gapdw(6)
 1001 format(7e13.4)
      write(61,1002) times,v1,v2,v3,v4,v5,v6,v7
 1002 format(8e13.4)
!
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

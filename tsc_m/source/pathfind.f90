        subroutine  pathfind
!
!               IPATH(i,j) contains the poloidal current paths m that
!                          touch wire i.
!               IPATH(i,1) = 0  says no paths touch wire i.
!
!               IPATH(i,j+1) = 0
!               IPATH(i,j  ).ne.0   says j paths touch wire i.
!
!               IPATH(i,j) = (-,+)m  says (tail,head) of path m touches wire i
!
!               R.O. Sayer  7 Dec 89
!
      USE CLINAM
      USE SCR21
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER ibug,i,j,jhit,m
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 rcl,rtest
!============
        data  rcl /0.0001_R8/
      data ibug /0/
!
        write (nout, 10)
   10   format(' PATHFIND Poloidal Current Paths             :'       /  &  
     & '                                                                 &  
     & wire  GRP   xwire   zwire    path    xpc1    zpc1     xpc2        &  
     & zpc2                                                              &  
     & ')
!--------------------------------------------------------------------------
        do 60  i=1,nwire
        do 20  j=1,8
   20   ipath(i,j) = 0
        jhit = 0
!
        do 50  m=1,mmax
        rtest = (xpc1(m) - xwire(i))**2 + (zpc1(m) - zwire(i))**2
        if (rtest.gt.rcl)   go to 30
!                                       Tail touches
        jhit = jhit + 1
        if (jhit.lt.9)   ipath(i,jhit) = - m
!ccccccc        jgroupc(m,1) = jgroup(i)
        go to 35
   30   rtest = (xpc2(m) - xwire(i))**2 + (zpc2(m) - zwire(i))**2
        if (rtest.gt.rcl)   go to 50
!                                       Head touches
        jhit = jhit + 1
        if (jhit.lt.9)   ipath(i,jhit) = + m
!ccccccc        jgroupc(m,2) = jgroup(i)
   35   if (ibug.eq.0)   go to 50
        if (jhit.eq.1)   then
        write (nout, 40)   i, igroupw(i), xwire(i), zwire(i),            &  
     &                     m,  xpc1(m), zpc1(m), xpc2(m), zpc2(m)
   40   format (2i5, 2f8.3, i8, 2f8.3, 1x, 2f8.3)
        else
        write (nout, 42)   m,  xpc1(m), zpc1(m), xpc2(m), zpc2(m)
                         endif
   42   format (26x,i8, 2f8.3, 1x, 2f8.3)
   50   continue
        if (mod(i,200).eq.0)   write (nout,10)
   60   continue
!--------------------------------------------------------------------------
        return
        end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

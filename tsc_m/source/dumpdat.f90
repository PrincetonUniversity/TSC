      subroutine dumpdat
!             ROS    19 Dec 1996
!        Writes the standard "ographa" file if acoef(931).eq.0.
!        Write special "ograph" file for ROS postprocessing.
!
      USE CLINAM
      USE SCR14
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER ios22,iadscr,irec,l,i,imin,imax,j,izero1,izero2
      INTEGER jgrgrp,ii
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 xprof,xx,xprofm
!============
      character*8  rosograph
!-----------------------------------------------------
!..open file for writing pltsav data for postprocessing
      if (acoef(931) .eq. 0.0_R8)   then
      if( numargs .lt. 1 ) then
         filename = 'ograph' // isuffix(1:1)
      else
         filename = 'ograph' // '.' // trim(suffix)
      end if
!     ograph(1:6) = 'ograph'
!     ograph(7:7) = isuffix(1:1)
      ograph = trim(filename)
      open(n13,file=trim(ograph),status='unknown',iostat=ios22)
      write(n13,1313) nrecord,name
 1313 format(i10,8a10,4i10)
      iadscr = 0
      do 450 irec=1,nrecord
      write(n13,1314) (pltsav(iadscr+l),l=1,lenscr)
      iadscr = iadscr + lenscr
  450 continue
 1314 format(5e20.12)
      close (n13)
                  endif
!----------------------------------------------------
!....> added 04/01/97 to produce output for ufiles      SCJ
!-----------------------------------------------------
!..open file for writing ufdata data for postprocessing
      if (acoef(3001).eq. 1.0_R8)   then
      if (numargs .lt. 1) then
         filename = 'ufiles' // isuffix(1:1)
      else
         filename = 'ufiles' // '.' // trim(suffix)
      end if
      ograph = trim(filename)
!     ograph(1:6) = 'ufiles'
!     ograph(7:7) = isuffix(1:1)
      open(n13,file=trim(ograph),status='unknown',iostat=ios22)
!
      xprof = acoef(42)
      if(xprof.le.0) xprof = ccon
      do 40 i=3,nx
      imin = i
      xx = xary(i)
      if(xx.gt.xprof) go to 45
   40 continue
      imin = 3
   45 continue
      xprofm = acoef(43)
      if(xprofm.le.0) xprofm = alx
      do 48 i=3,nx
      imax = i
      xx = xary(i)
      if(xx.gt.xprofm) go to 47
   48 continue
      imax = nx
   47 continue
      write(n13,1313) nrecord,name,imin,imax,numicsuf,pngroup
      if(numicsuf.le.0) go to 147
      write(n13,1147) (ichargst(i),i=1,numicsuf)
      write(n13,1147) (impnum(i),i=1,numicsuf)
  147 continue
 1147 format(18i5)
!
      call disk(35,1,big1,big2)
      write(n13,1314) (big2(i),i=1,nrecord)
      write(n13,1314) (xary(i),i=1,pnx)
      do 460 irec=1,nrecord
      write(n13,1314) ((ufdata(irec,j,nn),j=1,pnx),nn=1,30)
  460 continue
!
!     xmag and ip
      call disk(  8,  9,big1,big2)
      write(n13,1314) (big1(i),i=1,nrecord)
      write(n13,1314) (big2(i),i=1,nrecord)
!
!     diag and vsurf
      call disk( 13, 14,big1,big2)
      write(n13,1314) (big1(i),i=1,nrecord)
      write(n13,1314) (big2(i),i=1,nrecord)
!
!     q0   and qedge
      call disk( 15, 16,big1,big2)
      write(n13,1314) (big1(i),i=1,nrecord)
      write(n13,1314) (big2(i),i=1,nrecord)
!
!     beta and dn/dt from pellet   NOTE:  changed 8/6/99
      call disk( 18,201,big1,big2)
      write(n13,1314) (big1(i),i=1,nrecord)
      write(n13,1314) (big2(i),i=1,nrecord)
!
!     betap and li/2
      call disk( 22, 23,big1,big2)
      write(n13,1314) (big1(i),i=1,nrecord)
      write(n13,1314) (big2(i),i=1,nrecord)
!
!     te(0) and ti(i)
      call disk( 28, 29,big1,big2)
      write(n13,1314) (big1(i),i=1,nrecord)
      write(n13,1314) (big2(i),i=1,nrecord)
!
!     r0   and a
      call disk( 35, 36,big1,big2)
      write(n13,1314) (big1(i),i=1,nrecord)
      write(n13,1314) (big2(i),i=1,nrecord)
!
!     ptot and n-vol
      call disk( 46, 79,big1,big2)
      write(n13,1314) (big1(i),i=1,nrecord)
      write(n13,1314) (big2(i),i=1,nrecord)
!
!     vol  and no. part
      call disk( 84, 85,big1,big2)
      write(n13,1314) (big1(i),i=1,nrecord)
      write(n13,1314) (big2(i),i=1,nrecord)
!
!     zeff and imprad
      call disk(136,122,big1,big2)
      write(n13,1314) (big1(i),i=1,nrecord)
      write(n13,1314) (big2(i),i=1,nrecord)
!
!     #imp and "part
      call disk(206,207,big1,big2)
      write(n13,1314) (big1(i),i=1,nrecord)
      write(n13,1314) (big2(i),i=1,nrecord)
!
!..PF coil currents
      izero1=pgls+3*pncoil+1
      izero2=pgls+3*pncoil+2*pngroup
      jgrgrp = 0
      do 1900 ii=izero1,izero2,2
      jgrgrp = jgrgrp + 1
      call disk(ii+1,ii+2,big1,big3)
      write(n13,1314) (big1(i),i=1,nrecord)
1900  continue
      close (n13)
                  endif
!----------------------------------------------------
        write (nout, 220)   (j, j=1,9)
  220   format(/' masking array :  iexs(i,j) for separatrix'             &  
     &  /4x, 10i10)
        write (nout, 222)
  222   format (   '  j  ', 11('        0 ') )
        write (nout, 240)
        do 300 j=nzp,2,-1
        write (nout, 240)   j, zary(j), (iexs(i,j), i=2,nxp)
  240   format (i4, f6.2, 1x, 120i1)
  300   continue
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

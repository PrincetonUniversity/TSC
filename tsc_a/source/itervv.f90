!                                        ./tsc_a/source/itervv.f90
      subroutine itervv (ivvlo,ivvhi)
 
!********************************************************************
 !
!.....special for ITER       ROS.    15 Jan 1997.         needs  acoef(1) = 4
!
!     modified for ITER-TAC6 VV, SB, FW, divertor    5 Dec 1994.
!         changed to work for expanded z grid  12 Jan 1995
!		 changed to work for expanded z grid  18 Jan 1995   (alz=7.95)
!
!			modified to work for both 1997 and 2010 structures
!			acoef(1997) = 1. says 1997
!			acoef(1997) = 0. says 2010  (default)
      USE CLINAM
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER ivvhi,ivvlo,kpexist,impbx,jbpbx,jtpbx,kiterdiv,j,i,ii
      INTEGER n,jrow,ixlow,ixhi
      integer kpexist97, kpexist10
      integer jbpbx97, jtpbx97, jbpbx10, jtpbx10
!============
      dimension ivvlo(1),ivvhi(1), kpexist(3,20), kpexist97(3,20), kpexist10(3,20)
!
      data impbx, kiterdiv / 27, 12 /
      data jbpbx97, jtpbx97/ 24, 96 /
      data jbpbx10, jtpbx10/ 24,109 /
!
      data kpexist97 /16, 34,34,  17, 33,34,  18, 32,34,  19, 31,33,      &
     &		      20, 30,33,  21, 30,32,  21, 15,17,  22, 29,31,      &
     &		      22, 15,18,  23, 16,19,  23, 28,31,  24, 18,21,      &
     &		      24, 28,31,  25, 19,22,  25, 27,31,  26, 19,31,      &
     &		      27, 20,31,  95, 28,28,    6*0/
!
      data kpexist10 /18, 27,27,  19, 26,27, 20, 26,27,  21, 26,27,       &
     &		      22, 26,27,  23, 25,27, 24, 15,17,                   &
     &		      24, 24,27,  25, 16,18,  25, 23,27,  26, 19,27,      &
!!!!!!!!     &		      27, 20,28, 110, 21,21, 111, 21,21, 112, 21,21,      &
!!!!!!!!     &		      109, 21,23,      12*0/
     &		      27, 20,28, 109, 21,23,   21*0/
!
	if (acoef(1997) .eq. 0._R8)	then
					kpexist = kpexist10
					jbpbx   = jbpbx10
					jtpbx   = jtpbx10
					write (nout, '("itervv : 2010 structure")')
					write (nterm,'("itervv : 2010 structure")')
	else if (acoef(1997) .eq. 1._R8)	then
						kpexist = kpexist97
						jbpbx   = jbpbx97
						jtpbx   = jtpbx97
						write (nout, '("itervv : 1997 structure")')
						write (nterm,'("itervv : 1997 structure")')
	else
	write (nout, '("itervv error: acoef(1997) .ne. 0 or 1")')
	write (nout, '("              acoef(1997) = ", 1pe12.3)') acoef(1997)
	write (nterm, '("itervv error: acoef(1997) .ne. 0 or 1")')
	write (nterm, '("              acoef(1997) = ", 1pe12.3)') acoef(1997)
	stop
					endif
!
!                      If  halo, use divertor top
!                      If no halo, use VV bottom
      if (acoef(97).le.0.0_R8)   jbpbx = 23
!
      write (nout, 20)  impbx, jbpbx, jtpbx, kpexist
   20 format (/'TSC 97A - ITERVV:  impbx, jbpbx, jtpbx =', 3i5/' kpexist =',    &  
     &                     /(10x, 12i5))
      if(isym.eq.1) then
                    jbpbx = 1
                    jtpbx = 51
                    go to 101
                    endif
      do 100 j=2,jbpbx
      ivvlo(j) = impbx+1
      ivvhi(j) = impbx-1
  100 continue
  101 continue
      do 200 j=jtpbx,nzp
      ivvlo(j) = impbx+1
      ivvhi(j) = impbx-1
  200 continue
!-----------------
      do 380 j=jbpbx+1,jtpbx-1
      do 340 i=impbx,nxp
      ii = i
        do 320 n=1,nwire
!ccccccc    if (igroupw(n) .eq. kiterdiv)   go to 320
        if(iwire(n).eq.ii .and. jwire(n).eq.j) go to 346
  320   continue
  340 continue
  346 ivvhi(j) = ii
      do 360 i=2,impbx
      ii = 2+impbx-i
        do 350 n=1,nwire
!ccccccc    if (igroupw(n) .eq. kiterdiv)   go to 350
        if(iwire(n).eq.ii .and. jwire(n).eq.j) go to 366
  350    continue
  360 continue
  366 ivvlo(j) = ii
  380 continue
!------------------------------
!.....define masking arrays
      do 400 j=2,nzp
      do 400 i=2,nxp
      iexv(i,j) = 0
      if(i.le.ivvlo(j)) iexv(i,j) = 1
      if(i.ge.ivvhi(j)) iexv(i,j) = 1
  400 continue
!
      do 540  j=1,20
      jrow = kpexist(1,j)
      if (jrow.le.0)   go to 550
      ixlow = kpexist(2,j)
      ixhi  = kpexist(3,j)
      do 520  i=ixlow,ixhi
      iexv(i,jrow) = 0
  520    continue
  540    continue
  550    continue
!
      do 50 n=1,nwire
   50    iexv(iwire(n),jwire(n)) = 4
      call  maskout
!
      do 650 n=1,nwire
  650    iexv(iwire(n),jwire(n)) = 1
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

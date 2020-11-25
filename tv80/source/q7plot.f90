      subroutine q7plot (bc)
 
!**************************************************************************
!
!  q7plot
!
!  calls        q8plot
!
!  description  This was designed to only be called by the contour routines.
!
!**************************************************************************
 
!
!c author: f. n. fritsch
!
!c revised: dec14.1978 to conform to 7600 version
!
!c variable declarations:
!
      USE tvgxx1
      USE tvgxx2
      USE q7quad
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER ncont,j1,nr
      INTEGER jkk,karrz
!============
! idecl:  explicitize implicit REAL declarations:
      REAL   toolrg,fdf,rh,fnp,dx,dy,dmx,dmy,pp,xn
!============
!     common /tvgxx1/ ksent, kclip, krastr, ktyp, cxfact, cxcons
!    * , cyfact, cycons, cxmi, cxma, cymi, cyma, klx, kly, kpol
!    * , knamsb, knerr, crasmx, karrsz, cxarr(256)
!    * , cyarr(256), cxlast, cylast, kspace, kcfx, kcfy
!
!     common /tvgxx2/ kxarr(511), kyarr(511), kitcnt, klchar
!    * , kntenv, kdashv
!    * , knumch, kcoff, kcpt, kcprm(7,4), klindx, klindy
!    * , kchdx, kchdy, kindyo
!    * , klabel(4,7), kltv80(9), krasmx, kact
!    * , kdstat(8), kname(8), klenfi(8), klevd(8)
!    * , kdbuf(8), kdbufs(8), knfrm(8)
!    * , kepflg(8), kmxfil(8), kminbf(8), klerad(8)
!    * , kdfact(8), knumfi(8)
!
!
      dimension bc(*)
      REAL*8 bc
!     dimension bc(2)
!
!     common /q7quad/ kp, ctest, cxc(5), cyc(5), czc(5)
!    * , int, k01, k02, c01, cdf, czmin, czmax
!
! next line revised dec14.78 to conform to 7600 version
      save toolrg
      data toolrg /100.0 /
!
!c program statements:
!
      cxc(5) = cxc(1)
      cyc(5) = cyc(1)
      czc(5) = czc(1)
!
      if (k01 .ne. 0) go to 7
!
!  plot all possible contours at spacing fdf from c01.
      fdf = cdf
!
! next 3 lines revised 14-dec-78 to conform to 7600 version
      ncont = (czmax-czmin)/fdf
      if ((czmax-czmin) .gt. (toolrg*fdf)) go to 50
!     if (toolrg .eq. 100.0) go to 200
!
      if ((czmax-czmin)/fdf-2.0 ) 200, 200, 1000
!
 1000 continue
!
!  increase fdf if 'too many' contours pass through box...  (fdfdet)
      do 201 j1 = 1, 4
      rh = abs(czc(j1+1)-czc(j1))
      if (rh .eq. 0) go to 201
      fnp = fdf/rh
!
! 1/fnp is number of contours passing through.
      if (fnp .ge. 0.5 ) go to 201
      dx = abs(cxc(j1+1)-cxc(j1))
      dy = abs(cyc(j1+1)-cyc(j1))
      dmx = dx/(cxma-cxmi)
      dmy = dy/(cyma-cymi)
      pp = 0.005 /max(dmx,dmy,0.707 *(dmx+dmy))
      if (fnp .ge. pp) go to 201
      nr = (rh*pp)/fdf
!
! find closest power-of-two approximation to nr...
      jkk = 0
  210 nr = nr/2
      if (nr .eq. 0) go to 220
      jkk = jkk+1
      go to 210
!
  220 fdf = (2**jkk)*fdf
  201 continue
!
  200 continue
      xn = (czmin-c01)/fdf
      karrz = xn
      if (xn) 21, 22, 22
!
   21 karrz = karrz-1
      int = 0
      go to 23
!
   22 int = 1
   23 xn = karrz
      ctest = c01+xn*fdf
!
! the only way to end this loop is for ctest to exceed czmax.
   24 if (ctest-czmax) 25, 25, 9999
!
   25 call q8plot
   28 ctest = ctest+fdf
      karrz = karrz+1
      if (karrz) 24, 29, 24
!
   29 int = 1
      go to 24
!
! plot from an array of ctest values.
    7 if (k02) 8, 8, 9
!
    8 int = 1
      go to 10
!
    9 int = 0
   10 do 18  j1 = 1, k01
      if (j1-k02) 12, 11, 12
!
   11 int = 1
   12 ctest = bc(j1)
      if (ctest .gt. czmax) go to 18
      if (czmin .gt. ctest) go to 18
      call q8plot
!
   18 continue
      go to 9999
!
! next 3 lines revised 14-dec-78 to conform to 7600 version
   50 print *,"q7plot error: too many contours. change bc(2)    "
      kp = 1
!
9999  return
 
      END
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

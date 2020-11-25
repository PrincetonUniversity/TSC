!     -----------------------------------------------------------------|
      SUBROUTINE pSmoUnsm(xNpar,ySpec,n)
      USE DqlBins
      USE FeBins
      USE params
      USE WkAry
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      CHARACTER*10 MyString
      INTEGER n, i, iv, ivtabl, iymjr, iymnr, nsm2
      EXTERNAL ivtabl
      REAL*8    xNpar(n), ySpec(n), ymax, vbig
      REAL*8    PofV(NVELDIM)
      REAL*8    ZERO, ONE
      DATA    ZERO, ONE/                                                 &  
     &         0.0_R8, 1.0_R8/
      call BrodCast(nv,PofV,0.00_R8)
      call BrodCast(nv,wkv ,0.00_R8)
      nsm2 = (nsmoo - 1) / 2
      ymax = 0.00_R8
      vbig = 0.50_R8
      do 10 i=1,n
        iv = ivtabl(xNpar(i))
        wkv (iv) = wkv(iv) + ySpec(i)
        if (ymax .lt. wkv(iv)) ymax = wkv(iv)
        if (vbig .lt. abs(vpar(iv))) vbig = 1.00_R8
 10   continue
      call EZrnd2  (ymax, ymax, iymjr, iymnr)
      call convolve(nv,nsm2,qlsm,PofV,wkv)
      call EZsets(600,1000,475,725, -vbig,vbig, ZERO, ymax,1)
!     call GNsets(600,1000,475,725, -vbig,vbig, ZERO, ymax,1)
      call GNsets(500,1000,375,750, -vbig,vbig, ZERO, ymax,1)
      call ezwrit(600,750,'Relative power vs v/c$',0,0)
      call GNtitl(        'Relative power vs v/c '    )
      call EZaxes(2,5,iymjr,iymnr)
      do 20 iv=1,nv
        if (wkv(iv) .gt. 0.00_R8) call EZbars(vpar(iv),wkv(iv),1,'y')
 20   continue
                               call GNbars(vpar    ,wkv   ,nv,'y')
      call EZsets(600,1000,125,375, -vbig,vbig, ZERO, ymax,1)
!     call GNsets(600,1000,125,375, -vbig,vbig, ZERO, ymax,1)
      call GNsets(500,1000,  0,375, -vbig,vbig, ZERO, ymax,1)
      call ezwrit(600,400,'Relative power smoothed$',0,0)
      call GNtitl(        'Relative power smoothed'     )
      call EZaxes(2,5,iymjr,iymnr)
      call EZcurv(vpar,PofV,nv)
      call GNcurv(vpar,PofV,nv)
      write(MyString,'(''nsm2 '',i3,''$'')') nsm2
      call ezwrit(850,350,MyString,0,0)
      write(MyString,'(''nv/2 '',i3,''$'')') (nv-1)/2
      call ezwrit(850,325,MyString,0,0)
      return
      END
!     -----------------------------------------------------------------|
!     -                                                                |
!     -                                                                |
!     SpecInit ends   -------------------------------------------------|
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"

      subroutine gaxdrw (gxp, gyp, txt)
 
 
! 929 "gaxisi.F"
 
 
 
 
 
      USE tvaxis
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!============
! idecl:  explicitize implicit REAL declarations:
!     REAL   cxb,cyb,cxe,cye,ctickx,cticky,coffx,coffy,cac,cbc,cpx
!     REAL   cpy,gxs,gys
!============
      REAL     gxp, gyp, gxs, gys
      REAL*8 dgxp,dgyp,dgxs,dgys
      character*(*) txt
!
! this subroutine will draw the text and the tick mark for gaxis
!
!     logical klr
!     common /tvaxis/ cxb, cyb, cxe, cye, ctickx, cticky
!    * , coffx, coffy, cac, cbc, kori, ksize, kn1, klr
!    * , klxsav, klysav, cpx(30), cpy(30), kdiv, kepclp
!
!c program statements:
!
      if (klxsav .ne. 0) gxp = alog10(gxp)
      if (klysav .ne. 0) gyp = alog10(gyp)
! draw the tick mark
      call line (gxp*1.0d0 , gyp*1.0d0 ,                                 &  
     &  (gxp+cac*ctickx)*1.0d0 , (gyp+cbc*cticky)*1.0d0 )
! jump if horizontal text on a horizontal axis
      if (kori .eq. 0 .and. abs(cac) .le. .01 ) go to 1
! jump if vertical text on a vertical axis
      if (kori .eq. 1 .and. abs(cbc) .le. .01 ) go to 2
! get the starting position of the text
      gxs = gxp+cac*coffx
      gys = gyp+cbc*coffy
! change the starting position if the orientation requires it
      if (klr) go to 3
      if (kori .eq. 0) gxs = gxs-kn1*coffx
      if (kori .eq. 1) gys = gys-kn1*coffy
         go to 3
!
    1    gxs = gxp-kn1*coffx/2. 
         gys = gyp+2. *sign(coffy,cbc)
         go to 3
!
    2    gxs = gxp+2. *sign(coffx,cac)
         gys = gyp-kn1*coffy/2. 
! position the beam and draw the text
    3 dgxs = gxs
      dgys = gys
      call setlch (dgxs, dgys, 0, ksize, kori, -1)
      call gtext (txt, kn1, 0)
!
 9999 return
      END
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

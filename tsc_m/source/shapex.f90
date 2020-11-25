      subroutine shapex(ellipx,delta,x1,x2,z1,z2,psitar)
!.....6.95 shape
      USE CLINAM
      USE WALLCL
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER iaxis,jaxis,n,icount,iinc,jinc,iv,jv,ia,ja
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 delta,x1,x2,z1,z2,psitar,ellipx,xe,ze,xv,zv,psiv
      REAL*8 pinterp,xa,za,psia,ellip
!============
      dimension xe(50),ze(50)
!============      
!
      ellipx = 1._R8
      delta = 0._R8
      x1 = 1._R8
      x2 = .5_R8
      z1 = 1._R8
      z2 = -.1_R8
      if(lrswtch .gt. 0) return
      if(isym.eq.0) return
      iaxis = int((xmag-ccon)/deex)+2
      jaxis = int((zmag+alz*(1-isym))/deez)+2
      do 20 n=1,12
      icount = 0
      iinc = 0
      jinc = 0
      xv = xmag
      zv = zmag
      iv = iaxis
      jv = jaxis
      psiv = psimin
      go to(1,2,3,4,5,6,7,8,9,10,11,12),n
    1 iinc = 1
      go to 15
    2 iinc = -1
      go to 15
    3 jinc = 1
      xv = 0.5_R8*(xe(1)+xe(2))
      zv = ze(2)
      iv = int((xv-ccon)/deex) + 2
      jv = int((zv+alz*(1-isym))/deez) + 2
      psiv = pinterp(xv,zv,iv,jv)
      go to 15
    4 iinc = 1
      xv = 0.5_R8*(xe(1)+xe(2))
      zv = 0.5_R8*(ze(1)+ze(3))
      iv = int((xv-ccon)/deex) + 2
      jv = int((zv+alz*(1-isym))/deez) + 2
      psiv = pinterp(xv,zv,iv,jv)
      go to 15
    5 iinc = -1
      xv = 0.5_R8*(xe(1)+xe(2))
      zv = 0.5_R8*(ze(1)+ze(3))
      iv = int((xv-ccon)/deex) + 2
      jv = int((zv+alz*(1-isym))/deez) + 2
      psiv = pinterp(xv,zv,iv,jv)
      go to 15
    6 jinc = 1
      xv = 0.5_R8*(xe(4)+xe(5))
      zv = 0.5_R8*(ze(1)+ze(3))
      iv = int((xv-ccon)/deex) + 2
      jv = int((zv+alz*(1-isym))/deez) + 2
      psiv = pinterp(xv,zv,iv,jv)
      go to 15
    7 iinc = 1
      xv = 0.5_R8*(xe(4)+xe(5))
      zv = 0.5_R8*(ze(3)+ze(6))
      iv = int((xv-ccon)/deex) + 2
      jv = int((zv+alz*(1-isym))/deez) + 2
      psiv = pinterp(xv,zv,iv,jv)
      go to 15
    8 iinc = -1
      xv = 0.5_R8*(xe(4)+xe(5))
      zv = 0.5_R8*(ze(3)+ze(6))
      iv = int((xv-ccon)/deex) + 2
      jv = int((zv+alz*(1-isym))/deez) + 2
      psiv = pinterp(xv,zv,iv,jv)
      go to 15
    9 jinc = 1
      xv = 0.5_R8*(xe(7)+xe(8))
      zv = 0.5_R8*(ze(3)+ze(6))
      iv = int((xv-ccon)/deex) + 2
      jv = int((zv+alz*(1-isym))/deez) + 2
      psiv = pinterp(xv,zv,iv,jv)
      go to 15
   10 iinc = 1
      xv = 0.5_R8*(xe(7)+xe(8))
      zv = 0.5_R8*(ze(6)+ze(9))
      iv = int((xv-ccon)/deex) + 2
      jv = int((zv+alz*(1-isym))/deez) + 2
      psiv = pinterp(xv,zv,iv,jv)
      go to 15
   11 iinc = -1
      xv = 0.5_R8*(xe(7)+xe(8))
      zv = 0.5_R8*(ze(6)+ze(9))
      iv = int((xv-ccon)/deex) + 2
      jv = int((zv+alz*(1-isym))/deez) + 2
      psiv = pinterp(xv,zv,iv,jv)
      go to 15
   12 jinc = 1
      xv = 0.5_R8*(xe(10)+xe(11))
      zv = 0.5_R8*(ze(6)+ze(9))
      iv = int((xv-ccon)/deex) + 2
      jv = int((zv+alz*(1-isym))/deez) + 2
      psiv = pinterp(xv,zv,iv,jv)
      go to 15
   15 xa = xv
      za = zv
      ia = iv
      ja = jv
      psia = psiv
      xv = xv + iinc*deex
      zv = zv + jinc*deez
      iv = iv + iinc
      jv = jv + jinc
      if(iv.lt.iminn-1 .or.iv.gt.imaxx+1) go to 20
      if(jv.lt.jminn-1.or.jv.gt.jmaxx+1) go to 20
      psiv = pinterp(xv,zv,iv,jv)
      icount = icount + 1
!     if(icount.ge.10*(nx+nz)) go to 10
      if(psiv.lt.psitar) go to 15
      xe(n) = xa + (psitar-psia)*(xv-xa)/(psiv-psia)
      ze(n) = za + (psitar-psia)*(zv-za)/(psiv-psia)
   20 continue
      if(isym.eq.1) go to 21
      if(xe(1).ne.xe(2))                                                 &  
     &ellip = (ze(3)-ze(4))/(xe(1)-xe(2))
      delta = 0._R8
      x1 = xe(1)
      x2 = xe(2)
      z1 = ze(3)
      z2 = ze(4)
      return
   21 continue
!
      if(xe(1).ne.xe(2)) then
      ellipx = 2._R8*ze(12)/(xe(1)-xe(2))
      delta = 2._R8*(0.5_R8*(xe(1)+xe(2))-xe(12))/(xe(1)-xe(2))
      endif
!
      x1 = xe(1)
      x2 = xe(2)
      z1 = ze(12)
      z2 =-ze(12)
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

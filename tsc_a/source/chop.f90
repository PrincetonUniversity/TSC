      subroutine chop(curt,NCHOP,itype,VPULL,ierr)
!
!     calculates pulling voltage of chopper given the current
!     in the coil.  if there is more than one chopper, current
!     gets divided by no of choppers.  itype=1 for hx-chopper,
!                                            2 for x-chopper
!     ierr=1 if error is found
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER nchop,itype,ierr,i
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 vpull,curt,cur,r,a,b,c,e,f,cut,vold,vnew,dum
!============
      ierr=0
      if(itype.le.0 .or. itype.gt.3) go to 100
      if(NCHOP.le.0) go to 100
!
         cur = curt/NCHOP
!
!
      go to(1,2,2),itype
    1 continue
!
!.....values for hx-choppers
      r = 0.8_R8
         a = 6._R8
         b = 415
         c = 4.3_R8
!
      go to 3
    2 continue
!
!
!.....values for x-choppers
      r = 0.9_R8
          a = 6._R8
          b = 190
          c = 4.5_R8
    3 continue
 
!
         e = 1._R8/c
         f = 1._R8/(c-1._R8)
         cut = ((1.5_R8* b/r)**c/3._R8/a)**f
!
         if(cur.le.cut) then
!
         vold = cur * r
!
!
         do 21 i=1,99
         vnew = r * (cur - a * (vold/b)**c)
         if (abs((vnew-vold)/(vnew+vold)).le. 1.E-6_R8) go to 23
!
         vnew = (0.8_R8*vnew + vold)/1.8_R8
!
!
         dum = vold * .5_R8
         vold = max(vnew,dum)
 21      continue
      go to 100
 23      continue
!
         else
!
         vold = .67_R8* r * cur
!
         do 31 i=1,99
         vnew = cur - vold/r
         vnew = vnew/a
!
         vnew = b * vnew**e
         if (abs((vnew-vold)/(vnew+vold)).le. 1.E-6_R8) go to 33
         vold = vnew
 31      continue
      go to 100
 33      continue
      endif
!
      VPULL = vnew
      return
  100 continue
      ierr=1
      VPULL=0
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

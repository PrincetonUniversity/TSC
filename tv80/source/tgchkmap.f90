        integer function tgchkmap (xl,xr,yb,yt,vl,vr,vb,vt,itype)
 
!*************************************************************************
!
!  tgchkmap - check to see if the specified mapping is legal
!
!  synopsis     call tgchkmap (xleft,xright,ybot,ytop,xmin,xmax,
!          ymin,ymax,itype)
!
!        real xleft,xright    Users x range
!        real ybot,ytop       Users y range
!        real xmin,xmax       Virtual x range
!        real ymin,ymax       Virtual y range
!        integer itype     The type of mapping
!
!  description  Checks to see if the desired mapping is legal.  Sets the
!        mapping to one if it is not.  The type of maps specified
!        by itype are 1 (linear-linear), 2 (linear-log), 3 (log-
!        linear) and 4 (log-log)
!
!*************************************************************************
 
!============
! idecl:  explicitize implicit REAL declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      REAL   denom
!============
        REAL   xl,xr,yb,yt
        REAL   vl,vr,vb,vt
        integer itype
        REAL   epsilon
 
 
        epsilon = 1.0E-6 
 
 
        tgchkmap = 1
 
! x range should not be equivalent
 
        if (xl .eq. xr) then
          print *,'left: ',xl,'  right: ',xr
          call tgerror (12)
          return
        endif
        if (vl .eq. vr) then
          print *,'normalized left: ',vl,'  normalized right: ',vr
          call tgerror (12)
          return
        endif
 
! y range should not be equivalent
 
        if (yb .eq. yt) then
          print *,'bottom: ',yb,'  top: ',yt
          call tgerror (13)
          return
        endif
        if (vb .eq. vt) then
          print *,'normalized bottom: ',vb,'  normalized top: ',vt
          call tgerror (13)
          return
        endif
 
! minimum x should be less than maximum
 
        if (xl .gt. xr .or. vl .gt. vr) then
          call tgerror (14)
          return
        endif
 
! minimum y should be less than maximun
 
        if (yb .gt. yt .or. vb .eq. vt) then
          call tgerror (15)
          return
        endif
 
! log x should be greater than zero
 
        if (itype .eq. 3 .or. itype .eq. 4) then
          if (xl .le. 0 .or. xr .le. 0) then
            call tgerror (16)
            return
          endif
        endif
 
! log y should be greater than zero
 
        if (itype .eq. 2 .or. itype .eq. 4) then
          if (yb .le. 0 .or. yt .le. 0) then
            call tgerror (17)
            return
          endif
        endif
 
! virtual x should be in the range of 0 to 1
 
        if (vl .lt. 0.0 .or. vr .gt. 1.0 ) then
          call tgerror (18)
          return
        endif
 
! virtual y should be in the range of 0 to 1
 
        if (vb .lt. 0.0 .or. vt .gt. 1.0 ) then
          call tgerror (19)
          return
        endif
 
! make sure the range has enough significant bits
 
        denom = 0.5 * (abs(xl) + abs(xr))
        if (abs(xr-xl)/denom .lt. epsilon) then
          print *,'left: ',xl,'  right: ',xr
          call tgerror (22)
          return
        endif
 
        denom = 0.5 * (abs(yb) + abs(yt))
        if (abs(yt-yb)/denom .lt. epsilon) then
          print *,'bottom: ',yb,'  top: ',yt
          call tgerror (23)
          return
        endif
 
! mapping has my blessing
 
        tgchkmap = 0
        return
        end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

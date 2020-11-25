 
      SUBROUTINE ugridEXC(vector,npts,vmin,vmax,EXCLUDED)
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER npts, ipts
      REAL*8    vector(npts), vmin, vmax, rnpts1, dv
      REAL*8                                                             &  
     &        EXCLUDED, vmaxPLUS, vminPLUS, vmaxMINU, vminMINU, v
      REAL*8    ZERO
      DATA    ZERO/                                                      &  
     &         0.0_R8/
!
!     generate a uniform grid vector(j) from vmin to vmax
!     excluding points with abs value .le. EXCLUDED
!
      if(npts .le. 1 .or. vmax .le. vmin ) return
      rnpts1=npts-1
      vmaxPLUS = max( abs(EXCLUDED) , vmax )
      vminPLUS = max( abs(EXCLUDED) , vmin )
      vmaxMINU = min(-abs(EXCLUDED) , vmax )
      vminMINU = min(-abs(EXCLUDED) , vmin )
      dv =      max( (vmaxPLUS-vminPLUS) , ZERO) +                       &  
     &          max( (vmaxMINU-vminMINU) , ZERO)
!     dv =      max( (vmaxPLUS-vminPLUS) , 0.0 ) +
!    ^          max( (vmaxMINU-vminMINU) , 0.0 )
      dv=dv/rnpts1
      if (dv .le. 0.00_R8) return
!
      v = vmin
      do 10 ipts=1,npts
      vector(ipts) = v
        v = v + dv
        if (abs(v) .lt. abs(EXCLUDED)) go to 20
 10   continue
      return
!
!     return for normal exit when the excluded zone is not an issue
!
 20   continue
      v = vmax
      do 30 ipts=npts,1,-1
      vector(ipts) = v
        v = v - dv
        if (abs(v) .lt. abs(EXCLUDED)) go to 40
 30   continue
 40   return
!
!     return for unusual exit when the positive used zone comes
!     after a negative used zone
!
 
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"

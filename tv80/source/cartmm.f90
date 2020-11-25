        subroutine cartmm (n, xmin, xmax, x, inc)
 
!**************************************************************************
!
!  tgtv80.f - Routines used in the tv80 lib that shouldn't be there at all
!
!  contents    cartmm - calculate min and max of an array
!
!  description  The routines in this file are here for backwards support only.
!               These routines should not be used by any TV80 program.
!        However, some programmers used them, so they are here.
!
!    National Energy Research Supercomputer Center
!    Lawrence Livermore National Laboratory
!    University of California
!    Livermore, California 94550, USA
!
!    (C) Copyright 1991 The Regents of the University of California.
!    All Rights Reserved.
!
!**************************************************************************
 
!
!c variable declarations:
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER inc,n,i
!============
! idecl:  explicitize implicit REAL declarations:
      REAL   xmin,xmax,tmin,tmax
!============
        REAL   x(*)
 
        tmin = x(1)
        tmax = x(1)
        do 1 i = 1, n, inc
          if (x(i) .lt. tmin) tmin = x(i)
          if (x(i) .gt. tmax) tmax = x(i)
    1   continue
        xmin = tmin
        xmax = tmax
9999    return
        END
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

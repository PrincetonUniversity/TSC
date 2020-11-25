        subroutine vvexc1(jvvlo, jvvhi)
!
!               Set the arrays JVVLO, JVVHI to contain "J" indices that
!               define a region outside the CIT vacuum vessel from which
!               plasma is to be excluded by setting ETA to ETAV in RESIST.
!
!         ****  Assume no gaps in VV. ****
!         ****  This routine will not treat multiply connected regions; e.g.
!               a "bean" shaped or "indented" VV.
!
!                               R.O. Sayer   17 Oct 1989
!
!
! ====>  Extended 7/17/91 to also treat ASDEX-U Vacuum Vessel for acoef(1)=3
!
!
!               KVVGRP(1,2,3,4) = group numbers for CIT vacuum vessel
!
      USE CLINAM
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER jvvhi,jvvlo,kvvgrp,kmxx,i,n,k,j,iw,jw
!============
        dimension  jvvlo(2), jvvhi(2), kvvgrp(4)
        data  kvvgrp /20,21,22,23/
!============      
!
!      special output to terminal
!
      if(acoef(1).eq.2) write(nterm,3011)
      if(acoef(1).eq.3) write(nterm,3012)
 3011 format(" NOTE: Special logic assumes BPX Vacuum vessel ")
 3012 format(" NOTE: Special logic assumes ASDEX-U Vacuum vessel")
!
      kmxx = 4
!
        write (nout, 20)
   20   format (/' VVEXC : Exclude plasma from CIT vacuum vessel'        &  
     &        , /'     i    xary   jvvlo   jvvhi     zlo     zhi')
        do 100  i=2,nxp
        if (isym.eq.1)   then
                         jvvlo(i) = 1
                         jvvhi(i) = 1
        else
        jvvlo(i) = nzp
        jvvhi(i) = 2
                         endif
!
        do 40  n=1,nwire
        if (iwire(n).ne.i)   go to 40
      do 38 k=1,kmxx
        if (igroupw(n).eq.kvvgrp(k))   then
        if (zwire(n).le.0.0_R8.and.jwire(n).lt.jvvlo(i))  jvvlo(i)       &  
     & =jwire(n)
        if (zwire(n).ge.0.0_R8.and.jwire(n).gt.jvvhi(i))  jvvhi(i)       &  
     & =jwire(n)
                                       endif
   38   continue
   40   continue
        write (nout,150)    i, xary(i), jvvlo(i),jvvhi(i),               &  
     &                          zary(jvvlo(i)), zary(jvvhi(i))
  150   format (i7, f8.3, 2i8, 2f8.3)
  100   continue
!
!.....define masking array
      do 200 i=2,nxp
      do 200 j=2,nzp
      iexv(i,j) = 0
      if(j.le.jvvlo(i)) iexv(i,j) = 1
      if(j.ge.jvvhi(i)) iexv(i,j) = 1
  200 continue
      do 50 n=1,nwire
      do 48 k=1,kmxx
!     if(igroupw(n).eq.kvvgrp(k))   then
                                    iw = iwire(n)
                                    jw = jwire(n)
                                    iexv(iw,jw) = 1
!                                   endif
   48   continue
   50 continue
!
        return
        end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

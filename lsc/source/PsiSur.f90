!#include "f77_dcomplx.h"
!     -----------------------------------------------------------------
      SUBROUTINE PsiSur ( PsiVal , Rary , Zary , ISIZE , Npts )
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      EXTERNAL plasma2d
!     Finds surfaces of constant Psi, such that the flux on a
!     surface is PsiVal.  The first element of the array shows
!     how to start looking.  It moves to small  r  from the given
!     value until a Psi is correct; then it moves along B
!     to find the next value.
      CHARACTER*70 ErrMsg
      INTEGER i, ISIZE, Npts, NSCALE
      REAL*8    Bpol, deltf, deltS, f, fs, f1, f2,                       &  
     &        HUGE,                                                      &  
     &        PsiVal, r0, r1, r2, Rstart, Stepr,                         &  
     &                z0, z1, z2, Zstart, Stepz
      REAL*8                                                             &  
     &                    r,z, psi,Br,Bz,RBphi,omc,Tee,pe2,pi2,aio,ael
      REAL*8 Rary(ISIZE), Zary(ISIZE)
      REAL*8 AREAL
      DATA HUGE /1.0E+30_R8/

      Rstart= Rary(1)
      Zstart= Zary(1)
!     NSCALE= ISIZE divided by something like 3 or 4; 24/7 for a while. 5 Now.
!     NSCALE= ISIZE * 7 / 24 worked for a long time.  D-III-D needed 1/5
      NSCALE= ISIZE * 1 /  5
 
      Npts  = 1
      deltS = Rstart/AREAL(NSCALE)
      r     = Rstart
      z     = Zstart
      call     plasma2d (r,z, psi,Br,Bz,RBphi,omc,Tee,pe2,pi2,aio,ael)
      fs = Psi - PsiVal
 
      do 10 i=1, NSCALE
      r     = Rstart - deltS
      call     plasma2d (r,z, psi,Br,Bz,RBphi,omc,Tee,pe2,pi2,aio,ael)
      if ( Psi .ge. HUGE ) return
      f = Psi - PsiVal
!
!     .                                 If zero crossing, kick out
      if((f.gt.0._R8.and.fs.le.0._R8).or.                                &  
     &   (f.lt.0._R8.and.fs.ge.0._R8))    go to 12
!
!     .                                 If we seem to be going away from
!     .                                 the zero crossing, then we must have
!     .                                 gone thru an extremum; kick out
      if( abs(f) .ge. abs(fs)  ) then
        write(ErrMsg,'(                                                  &  
     &''leaving 0-crossing at R,Psi,PsiVal:'',                           &  
     &  3(1x,1pe9.2))') r, psi, psival
        call LSCwarn(ErrMsg)
        go to 25
      else
        fs    = f
      endif
   10 Rstart= r
      write(ErrMsg,'(''no start point in PsiSur '',                      &  
     &               ''R,Psi,PsiVal:'', 3(1x,1pe9.2))')
      call LSCwarn(ErrMsg)
      return
 
   12 continue
      deltf = fs-f
      if (deltf .eq. 0._R8) then
        call LSCwarn( ' fs=f in PsiSur')
      else
        Rstart= Rstart - deltS * fs/(fs-f)
      endif
 25   r0    = Rstart
      z0    = Zstart
      Npts = 1
      Rary(Npts) = Rstart
      Zary(Npts) = Zstart
 
      do 30 i=2,ISIZE
      call     plasma2d (r0,z0, psi,Br,Bz,RBphi,omc,Tee,pe2,pi2,aio,ael)
      Bpol  = sqrt( Br*Br + Bz*Bz )
      Stepr = deltS * Br/Bpol
      Stepz = deltS * Bz/Bpol
!
!     .                                 Take a step along the field line
      r     = r0 + Stepr
      z     = z0 + Stepz
!     .                                 Take a tiny step perp to the last one
!     .                                 First outward from the plasma ...
!     .                                 assuming current in phi direction
      r1    = r  + Stepz/10._R8
      z1    = z  - Stepr/10._R8
!     .                                 Take a tiny step the other way..inward
      r2    = r  - Stepz/10._R8
      z2    = z  + Stepr/10._R8
!     .                                 We have 3 points in a line, all
!     .                                 pretty near the right place.
!     .                                 Find the point with the smallest
!     .                                 error
      call     plasma2d (r ,z , psi,Br,Bz,RBphi,omc,Tee,pe2,pi2,aio,ael)
        f  = abs(psi-PsiVal)
      call     plasma2d (r1,z1, psi,Br,Bz,RBphi,omc,Tee,pe2,pi2,aio,ael)
        f1 = abs(psi-PsiVal)
      call     plasma2d (r2,z2, psi,Br,Bz,RBphi,omc,Tee,pe2,pi2,aio,ael)
        f2 = abs(psi-PsiVal)
!
        fs=min(f,f1,f2)
!     .                                 The point labeled '2' is preferred be-
!     .                                 cause it tends to be in the plasma
      if (f2 .eq. fs) then
        r = r2
        z = z2
      else if (f1 .eq. fs) then
        r = r1
        z = z1
      endif
!
      Npts = i
      Rary(Npts) = r
      Zary(Npts) = z
      if( i .gt. 5 .and.                                                 &  
     & (r-Rstart)**2+(z-Zstart)**2 .lt. 4.0_R8*deltS**2) return
      r0=r
   30 z0=z
      write(ErrMsg,'(''no closed psi; i,R,Z,Psi:'',                      &  
     & i5,3(1x,1pe9.2) )') npts, r, z, psi
      call LSCwarn(ErrMsg)
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
! 24May2005 fgtok

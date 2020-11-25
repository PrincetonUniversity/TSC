!#include "f77_dcomplx.h"
!     -----------------------------------------------------------------
      SUBROUTINE RayDamp
      USE dielec
      USE Doflags
      USE FeBins
      USE params
      USE PIetc
      USE PlPr
      USE RayBins
      USE RayWrk
      USE WkAry
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      CHARACTER*70 ErrMsg
      INTEGER i, j, iry, izn, ThisRay
      INTEGER iv
      REAL*8    dum
      REAL*8    xz(NZONDIM), yX(NZONDIM), yQL(NZONDIM)
      REAL*8    MxdPdZ
      REAL*8    AREAL
      DATA    MxdPdZ / 0.00010_R8/
!
!     compute exponential decrement, deposit into yX (Maxwellian) and
!     yQL (quasilinear)
      izn = nzones
!     nzones is stored in, read from   ray.dat
!
      do 30 iry = 1, nrays
        xz(1)  = 1._R8
        yX(1)  = 1._R8
        yQL(1) = 1._R8
!       First do Maxwellian damping dlnPdsX
!
        do 10 i=1,izn-1
          xz(i) = AREAL(i)
 
            iv = ivind(i , iry )
!                                       If  ivind=0, the ray was stopped before
!                                       reaching this izone
!                                       No calculation is appropriate.
            if (iv .eq. 0) then
              yX(i+1) = yX(i)
              go to 10
            endif
 
          dum = 1._R8+ dlnPdsX(i,iry)
!                                       Limit the damping per zone to the
!                                       value MxdPdZ
          if (dum .lt. MxdPdz) then
            write(ErrMsg,'( '' Mxw damp!iry,izn,iv,'',                   &
     &                      ''dlnPdsX     :'',i4,i4,i4,                  &
     &                        1pe11.3)')                                 &
     &                        iry, i, ivind(i,iry),dlnPdsX(i,iry)
!dbg        call LSCwarn( ErrMsg )
            dum = MxdPdZ
          endif
!                                       Guard against amplification of wave!
          if (dum .gt. 1._R8+MxdPdZ) then
            write(ErrMsg,'( '' Mxw ampl!iry,izn,iv,'',                   &
     &                      ''dlnPdsX     :'',i4,i4,i4,                  &
     &                        1pe11.3)')                                 &
     &                        iry, i, ivind(i,iry),dlnPdsX(i,iry)
!dbg        call LSCwarn( ErrMsg )
            dum = 1._R8
 
          endif
!         yX(i) =  yX(i-1) * dum ! before jul 92
          yX(i+1) =  yX(i) * dum
          if(dum .eq. MxdPdZ) then
            do 9 j = i+1,izn
              yX(j) = 0.00_R8
 9          continue
            go to 11
          endif
 10     continue
 11     continue
!
!       Next do Quasilinear damping dlnPdsK with the Kernel
!
        do 20 i=1,izn-1
          rFudgDmp(i,iry) = 1.00_R8
          xz(i) = AREAL(i)
            iv = ivind(i , iry )
!                                       If  ivind=0, the ray was stopped before
!                                       reaching this izone (izn).
!                                       No calculation is appropriate.
            if (iv .eq. 0) then
              yQL(i+1) = yQL(i)
              go to 20
            endif
 
            dum = dlnPdsK(i,iry)*dfdv(ivind(i,iry),izind(i,iry),2)
 
!                                       But if we are at negative velocity
!                                       then dlnPdsK*dfdv needs a minus sign.
          if(iv .lt. IvZero) then
            dum = - dum
          endif
!
!                                       Make dum a positive number, hoped to
!                                       be small compared to 1
            dum = - dum
 
 
          if (dum .lt. 0.01_R8) then
            rFudgDmp(i,iry) = 1.0_R8- dum/2.0_R8+ dum**2/6.0_R8
          else
            rFudgDmp(i,iry) = ( 1.0_R8- exp(-dum) ) / dum
          endif
            yQL(i+1) = yQL(i) * exp (-dum)
 
          if (yQL(i+1) .lt. 1.0E-20_R8) then
            do 19 j = i+1,izn
              yQL(j) = 0.00_R8
 19         continue
            go to 21
          endif
 20     continue
 21     continue
!
        call svmult(izn, yQL, power(1, iry), power(1, iry))
!       compute power along ray
 30   continue
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"
! 24May2005 fgtok

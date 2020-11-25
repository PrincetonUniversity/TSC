!
!     -----------------------------------------------------------------
!
      REAL*8 FUNCTION ran3(idum)
!     Transcribed 1993 by D. W. Ignat
!     References:
!     W. H. Press et al, Numerical Recipes in Fortran, p199
!     D. Knuth, Seminumerical Algorithms, Vol 2 of The Art of Computer
!          Programming
!     Returns a uniform random deviate between 0.0 and 1.0.
!     Set idum to any negative value to initialize or reinitialize the
!     sequence.
!     Substitute the CCCommented lines for the ones following to
!     render the routine entirely floating point.
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER iabs
!============
      INTEGER idum, MBIG, MSEED, MZ, iff, mj,mk
      INTEGER i,ii,k, inext, inextp
!     .                                 55 dimension is special; no changes!
      INTEGER ma(55)
      REAL*8    FAC
!CC   IMPLICIT REAL*4(M)
!CC   PARAMETER (MBIG=4000000.,   MSEED=1618033., MZ=0.,FAC=1./MBIG)
      PARAMETER (MBIG=1000000000, MSEED=161803398,MZ=0, FAC=1._R8/MBIG)
      DATA iff /0/
!                                       Initialize
      if(idum .lt. 0 .or. iff .eq. 0) then
        iff=1
!     .                                 Initialize ma(55) using the seed idum
!     .                                 and the large number MSEED
        mj=MSEED-iabs(idum)
        mj=mod(mj,MBIG)
        ma(55)=mj
        mk=1
!     .                                 Now initialize the rest of the table
!     .                                 in a slightly random order
!     .                                 with numbers that are not especially
!     .                                 random
        do 11 i=1,54
          ii=mod(21*i,55)
          ma(ii)=mk
          mk=mj-mk
          if(mk .lt. MZ) mk=mk+MBIG
          mj=ma(ii)
 11     continue
!     .                                 We randomize them by
!     .                                 'warming up the generator'
        do 13 k=1,4
          do 12 i=1,55
            ma(i)=ma(i)-ma(1+mod(i+30,55))
            if(ma(i) .lt. MZ) ma(i)=ma(i)+MBIG
 12       continue
 13     continue
!     .                                 Prepare indices for our first
!     .                                 generated number.   The constant 31
!     .                                 is special...see Knuth
        inext=0
        inextp=31
        idum=1
      endif
!
!     .                                 Here is where we start usually
!     .                                 Increment inext, wrap around 56 to 1
      inext=inext+1
      if(inext .eq. 56)inext=1
      inextp=inext+1
!     .                                 Ditto for inextp
      if(inextp .eq. 56)inextp=1
!     .                                 Now generate a new random number
      mj=ma(inext)-ma(inextp)
!     .                                 Be sure it is in the range
      if(mj .lt. MZ)mj=mj+MBIG
!     .                                 Store it, and output ran3
      ma(inext)=mj
      ran3=mj*FAC
      return
      END
!
!     -----------------------------------------------------------------
!
!
!
!                                                                      |
!                                                                      |
!     PredcLSC ends                         ---------------------------|
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"

      subroutine simson(psmout,pfin,pdxin,nxin)
!
!
!  Simpson's one-third rule for integrating a function over fixed steps
!
!  S = [ f(1) + 4 * f(2) + 2 * f(3) + 4 * f(4) + 2 * f(5) + ...
!             + 2 * f(nsimp-2) + 4 * f(nsimp-1) + f(nsimp) ] * dxin / 3.
!
!  where nsimp is an odd - valued integer
!
!  usage:
!         call simson ( psmout, pfin, pdxin, nxin )
!
!  psmout = output result
!  pfin   = input integrand array
!  pdxin  = input fixed value of dx
!  nxin   = input number of values of pfin (should be odd - valued)
!
!    If nxin is given and even value, Simpson's rule is then used up to
!  the nxin - 1 value of the pfin array and then the trapezoidal rule
!  is used over the last interval.
!
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER nxin,nsimp,j
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 pfin,pdxin,psmout,zsimp
!============
      dimension pfin(1)
!============      
!
      if ( nxin .lt. 2 ) then
        psmout = 0._R8
!tss    write (nout,*) '  WARNING -- nxin .lt. 2 in sbrtn simson '
        return
      endif
!
!  Simpson's rule must use odd number of nodes (nsimp)
!
      nsimp = nxin - mod(nxin+1,2)
      zsimp = 0._R8
      do 10 j=2,nsimp-1
        if ( mod(j,2) .eq. 0 ) then
          zsimp = zsimp + 4._R8* pfin(j)
        else
          zsimp = zsimp + 2._R8* pfin(j)
        endif
  10  continue
!
      psmout = ( pfin(1) + zsimp + pfin(nsimp) ) * pdxin / 3._R8
!
!  use trapezoidal rule if an extra zone is left over
!    after Simpson's rule
!
      if ( nxin .gt. nsimp ) then
        psmout = psmout + 0.5_R8* pdxin * ( pfin(nxin-1) + pfin(nxin) )
      endif
!
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

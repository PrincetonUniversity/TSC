      subroutine boundc
!
      USE CLINAM
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!cj dir$ ivdep
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER i,j
!============
      do 401 i=2,nxp
      psi(i,nz+2) = 2._R8*psi(i,nz+1) - psi(i,nz)
      w(i,nz+2) = 0.0_R8
      g(i,nz+2) = gzero/xsqoj(i)
      r(i,nz+2) = r(i,nz+1)
      q(i,nz+2) = q(i,nz+1)
      b(i,nz+1) = 0._R8
      u(i,nz+2) = 0._R8
      psio(i,nz+2) = 2._R8*psio(i,nz+1) - psio(i,nz)
      wo(i,nz+2) = 0.0_R8
      go(i,nz+2) = gzero/xsqoj(i)
      ro(i,nz+2) = ro(i,nz+1)
      qo(i,nz+2) = qo(i,nz+1)
      bo(i,nz+1) = 0._R8
      uo(i,nz+2) = 0._R8
!
!.....inside vessel, in shaddow of separatrix
      ajphi(i,nzp+1)=ajphi(i,nz)
  401 continue
      if(isym.eq.1) go to 404
      do 403 i=2,nxp
      psi(i,1)=2.0_R8*psi(i,2)-psi(i,3)
      w(i,2) = 0.0_R8
      g(i,2) = gzero/xsqoj(i)
      r(i,2) = r(i,3)
      q(i,2) = q(i,3)
      b(i,2) = 0._R8
      u(i,2)=0.0_R8
      ajphi(i,1) = ajphi(i,3)
      psio(i,1)=2.0_R8*psio(i,2)-psio(i,3)
      wo(i,2) = 0.0_R8
      go(i,2) = gzero/xsqoj(i)
      ro(i,2) = ro(i,3)
      qo(i,2) = qo(i,3)
      bo(i,2) = 0._R8
      uo(i,2)=0.0_R8
  403 continue
      go to 405
  404 continue
!cj dir$ ivdep
      do 406 i=2,nxp
      psi(i,1) = psi(i,3)
      w(i,2) = -w(i,3)
      pr(i,2) = pr(i,3)
      g(i,2) = g(i,3)
!..boundary condition for antisymmetric toroidal field
      if(jsym.eq.-1) g(i,2) = -g(i,3)
      r(i,2) = r(i,3)
      q(i,2) = q(i,3)
      b(i,1) =-b(i,3)
      u(i,2) = u(i,3)
      ajphi(i,1) = ajphi(i,3)
      psio(i,1) = psio(i,3)
      wo(i,2) = -wo(i,3)
      go(i,2) = go(i,3)
!..boundary condition for antisymmetric toroidal field
      if(jsym.eq.-1) go(i,2) = -go(i,3)
      ro(i,2) = ro(i,3)
      qo(i,2) = qo(i,3)
!...corrected 3/12/86
      bo(i,1) = -bo(i,3)
      uo(i,2) = u(i,3)
  406 continue
  405 continue
!cj dir$ ivdep
      do 402 j=1,nzp+1
      psi(nx+2,j) = 2._R8*psi(nx+1,j) - psi(nx,j)
      w(nx+2,j) = 0.0_R8
      g(nx+2,j) = gzero/xsqoj(nx+2)
      r(nx+2,j) = r(nx+1,j)
      q(nx+2,j) = q(nx+1,j)
      b(nx+1,j) = 0._R8
      u(nx+2,j) = 0._R8
      psio(nx+2,j) = 2._R8*psio(nx+1,j) - psio(nx,j)
      wo(nx+2,j) = 0.0_R8
      go(nx+2,j) = gzero/xsqoj(nx+2)
      ro(nx+2,j) = ro(nx+1,j)
      qo(nx+2,j) = qo(nx+1,j)
      bo(nx+1,j) = 0._R8
      uo(nx+2,j) = 0._R8
!
      r(2,j) = r(3,j)*ajey(2)/ajey(3)
      q(2,j) = q(3,j)*ajey(2)/ajey(3)
      g(2,j) = gzero/xsqoj(2)
!..boundary condition for antisymmetric toroidal field
      if(jsym.eq.-1) g(2,j) = -gzero/xsqoj(2)
      w(2,j) = 0.0_R8
!
      ro(2,j) = ro(3,j)*ajey(2)/ajey(3)
      qo(2,j) = qo(3,j)*ajey(2)/ajey(3)
      go(2,j) = gzero/xsqoj(2)
!..boundary condition for antisymmetric toroidal field
      if(jsym.eq.-1) go(2,j) = -gzero/xsqoj(2)
      wo(2,j) = 0.0_R8
      psi(1,j) = 2._R8*psi(2,j) - psi(3,j)
      psio(1,j) = 2._R8*psio(2,j) - psi(3,j)
  402 continue
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

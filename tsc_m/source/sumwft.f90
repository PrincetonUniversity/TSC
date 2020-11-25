      subroutine  sumwft (n1, n2, frdum, fzdum)
!*************************************************************************
!
!           sum wire forces due to toroidal currents
!              ROS   13 June 1992
      USE CLINAM
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER n2,n1,ii,n,i,j
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 frdum,fzdum,dpsir,dpsiz
!============
      frdum = 0.0_R8
      fzdum = 0.0_R8
      if (n1.le.0 .or. n1.gt.nwire)   return
      if (n2.lt.n1.or. n2.gt.nwire)   return
        do 30  ii=n1,n2
        n = ncoil - nwire + ii
        cwires(ii) = ccoil(n) * udsi
        i = iwire(ii)
        j = jwire(ii)
        dpsir = 0.2_R8* (psi(i+1,j) - psi(i-1,j))                        &  
     &        + 0.4_R8* (psi(i+2,j) - psi(i-2,j))
!
        if (j.lt.3)   then
                      dpsiz = psi(i,j+1) - psi(i,j-1)
        else
        dpsiz = 0.2_R8* (psi(i,j+1) - psi(i,j-1))                        &  
     &        + 0.4_R8* (psi(i,j+2) - psi(i,j-2))
                      endif
        frdum = frdum -0.5_R8* cwires(ii) * dpsir / deex
        fzdum = fzdum -0.5_R8* cwires(ii) * dpsiz / deez
   30   continue
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

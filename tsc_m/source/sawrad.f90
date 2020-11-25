      subroutine sawrad(rsaw)
!
      USE CLINAM
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER i,j
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 rsaw,qval,qget,psival,qold,rout,rin
!============
      rsaw = 0._R8
      if(isurf.eq.0) return
      i = imag
      j = jmag
      qval = qget(psi(imag,jmag))
      if(qval.gt.1.0_R8) return
   10 i = i+1
      if(i.gt.nxp) go to 100
      psival = psi(i,j)
      if(psival.gt.psilim) go to 100
      qold = qval
      qval = qget(psival)
      if(qold.lt.1.0_R8.and. qval.ge.1.0_R8) go to 50
      go to 10
   50 rout = ((qval-1.0_R8)*xary(i-1) + (1.0_R8-qold)*xary(i))           &  
     &     / (qval-qold)
      i = imag
      qval = qget(psi(imag,jmag))
   20 i = i-1
      if(i.lt.2) go to 100
      psival = psi(i,j)
      if(psival.gt.psilim) go to 100
      qold = qval
      qval = qget(psival)
      if(qold.lt.1.0_R8.and. qval.ge.1.0_R8) go to 60
      go to 20
   60 rin = ((qval-1.0_R8)*xary(i+1) + (1.0_R8-qold)*xary(i))            &  
     &    / (qval-qold)
      rsaw = 0.5_R8*(rout-rin)
      if(rsaw.lt.0) rsaw = 0._R8
  100 return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

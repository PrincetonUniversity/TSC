      SUBROUTINE R8TOMSQZ(N,NA,AR,AI,BR,BI,ALFR,ALFI,BETA,ZVR,ZVI,IFAIL)    
!
!dmc      subroutine echo(lbl,arr,jmaxm)
!dmc      character*(*) lbl
!dmc      real*8 arr(0:jmaxm)
!
!dmc      write(6,1001) lbl
!dmc 1001 format(/5x,a)
!dmc      write(6,1002) (arr(j),j=0,jmaxm)
!dmc 1002 format(4(1x,1pe17.10))
!
!dmc      return
!dmc      end
!dmc -- real*8 version generated using `fgtok'.
!  names changed:  tomsqz -> r8tomsqz, cqzhes -> r8cqzhes, etc.
!  all constants converted to "D" expontent form
!  all declarations remapped to REAL*8 / COMPLEX*16
!  cpp for standardizing REAL/COMPLEX intrinsics:
!
!#include "f77_dcomplx.h"
!
!     SUBROUTINE R8TOMSQZ(N,NA,AR,AI,BR,BI,ALFR,ALFI,BETA,ZVR,ZVI,IFAIL)
!-----------------------------------------------------------------------
! TOMSQZ  written by P. Strand 27-apr-98,       elfps@elmagn.chalmers.se
!-----------------------------------------------------------------------
! CQZHES, CQZVEC, and CQZVAL: Fortran subroutines implementing the QZ
! algorithm for solving the generalized eigenvalue problem for complex
! matrices. (See B.S.C Garbow, ACM TOMS 4 (1978) pp. 404-410.).
!-----------------------------------------------------------------------
!
! ON INPUT THE GENERALIZED EIGENVALUE PROBLEM IS DEFINED THROUGH THE
! COMPLEX MATRICES
!
!       A = cmplx (AR, AI)  AND   B = cmplx (BR, BI)
!
! WHERE LEADING DIMENSION N IS AS DEFINED IN THE CALLING ROUTINE AND
! WHERE NA IS THE ROW  RANGE IN THE CURRENT PROBLEM. THE EIGENVALUE
! PROBLEM IS THEN DEFINED THROUGH
!
!       A x = w B x
!
! WHERE  THE COMPLEX EIGENVECTORS
!
!       x = cmplx (ZVR, ZVI)
!
! TOGETHER WITH THE COMPLEX EIGENVALUE
!
!        w = cmplx(alfr, alfi)/beta
!
! IS OUTPUT FROM THE ROUTINE
!
! IFAIL WILL BE NONZERO IF CONVERGENCE HAS NOT BEEN REACH WITHIN 50
! ITERATIONS
!-----------------------------------------------------------------------
! DECLARATIONS FOR INPUT VARIABLES
!-----------------------------------------------------------------------
 
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER N, NA
      REAL*8 AR(N,NA),AI(N,NA),BR(N,NA),BI(N,NA)
 
!-----------------------------------------------------------------------
! DECALRATIONS FOR OUTPUT VARIABLES
!-----------------------------------------------------------------------
 
      REAL*8 ALFR(N),ALFI(N),BETA(N)
      REAL*8 ZVR(N,NA), ZVI(N,NA)
      INTEGER IFAIL
 
!-----------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------
 
      LOGICAL WANTX
      REAL*8 EPS1
      REAL*8 ZERO,ZONE
 
!-----------------------------------------------------------------------
! START OF ACTUAL CODING
!-----------------------------------------------------------------------
 
      WANTX = .TRUE.
      EPS1  = -0.0_R8
 
      CALL R8CQZHES(N,NA,AR,AI,BR,BI,WANTX,ZVR,ZVI)
      CALL R8CQZVAL(N,NA,AR,AI,BR,BI,EPS1,ALFR,ALFI,BETA,WANTX,          &  
     &            ZVR,ZVI,IFAIL)
      CALL R8CQZVEC(N,NA,AR,AI,BR,BI,ALFR,ALFI,BETA,ZVR,ZVI)
      RETURN
      END
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

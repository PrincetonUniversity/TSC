!#include "f77_dcomplx.h"
      SUBROUTINE CQZVEC(NM,N,AR,AI,BR,BI,ALFR,ALFI,BETA,ZR,ZI)
!
!     ------------------------------------------------------------------
!
!
!============
! idecl:  explicitize implicit REAL declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      REAL*8 zero
!============
      INTEGER I,J,K,M,N,EN,II,JJ,NA,NM,NN
      REAL*8 AR(NM,N),AI(NM,N),BR(NM,N),BI(NM,N),ALFR(N),ALFI(N),        &  
     &       BETA(N),ZR(NM,N),ZI(NM,N)
      REAL*8 R,T,RI,TI,XI,ALMI,ALMR,BETM,EPSB
      REAL*8 abs
      COMPLEX*16 Z3
      COMPLEX*16 CMPLX
      REAL*8 AREAL,AIMAG
!
!
!
!
!
!     THIS SUBROUTINE IS A COMPLEX ANALOGUE OF THE FOURTH STEP OF THE
!     QZ ALGORITHM FOR SOLVING GENERALIZED MATRIX EIGENVALUE PROBLEMS,
!     SIAM J. NUMER. ANAL. 10, 241-256(1973) BY MOLER AND STEWART.
!
!     THIS SUBROUTINE ACCEPTS A PAIR OF COMPLEX MATRICES IN UPPER
!     TRIANGULAR FORM, WHERE ONE OF THEM FURTHER MUST HAVE REAL DIAGONAL
!     ELEMENTS.  IT COMPUTES THE EIGENVECTORS OF THE TRIANGULAR PROBLEM
!     AND TRANSFORMS THE RESULTS BACK TO THE ORIGINAL COORDINATE SYSTEM.
!     IT IS USUALLY PRECEDED BY  CQZHES  AND  CQZVAL.
!
!     ON INPUT-
!
!        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
!          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
!          DIMENSION STATEMENT,
!
!        N IS THE ORDER OF THE MATRICES,
!
!        A=(AR,AI) CONTAINS A COMPLEX UPPER TRIANGULAR MATRIX,
!
!        B=(BR,BI) CONTAINS A COMPLEX UPPER TRIANGULAR MATRIX WITH REAL
!          DIAGONAL ELEMENTS.  IN ADDITION, LOCATION BR(N,1) CONTAINS
!          THE TOLERANCE QUANTITY (EPSB) COMPUTED AND SAVED IN  CQZVAL,
!
!        ALFR, ALFI, AND BETA ARE VECTORS WITH COMPONENTS WHOSE
!          RATIOS ((ALFR+I*ALFI)/BETA) ARE THE GENERALIZED
!          EIGENVALUES.  THEY ARE USUALLY OBTAINED FROM  CQZVAL,
!
!        Z=(ZR,ZI) CONTAINS THE TRANSFORMATION MATRIX PRODUCED IN THE
!          REDUCTIONS BY  CQZHES  AND  CQZVAL, IF PERFORMED.
!          IF THE EIGENVECTORS OF THE TRIANGULAR PROBLEM ARE
!          DESIRED, Z MUST CONTAIN THE IDENTITY MATRIX.
!
!     ON OUTPUT-
!
!        A IS UNALTERED,
!
!        B HAS BEEN DESTROYED,
!
!        ALFR, ALFI, AND BETA ARE UNALTERED,
!
!        Z CONTAINS THE EIGENVECTORS.  EACH EIGENVECTOR IS NORMALIZED
!          SO THAT THE MODULUS OF ITS LARGEST COMPONENT IS 1.0 .
!
!     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
!     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!
!     ------------------------------------------------------------------
!
      ZERO = 0.0_R8
      IF (N .LE. 1) GO TO 1001
      EPSB = BR(N,1)
!     ********** FOR EN=N STEP -1 UNTIL 2 DO -- **********
      DO 800 NN = 2, N
         EN = N + 2 - NN
         NA = EN - 1
         ALMR = ALFR(EN)
         ALMI = ALFI(EN)
         BETM = BETA(EN)
!     ********** FOR I=EN-1 STEP -1 UNTIL 1 DO -- **********
         DO 700 II = 1, NA
            I = EN - II
            R = 0.0_R8
            RI = 0.0_R8
            M = I + 1
!
            DO 610 J = M, EN
               T = BETM * AR(I,J) - ALMR * BR(I,J) + ALMI * BI(I,J)
               TI = BETM * AI(I,J) - ALMR * BI(I,J) - ALMI * BR(I,J)
               IF (J .EQ. EN) GO TO 605
               XI = T * BI(J,EN) + TI * BR(J,EN)
               T = T * BR(J,EN) - TI * BI(J,EN)
               TI = XI
  605          R = R + T
               RI = RI + TI
  610       CONTINUE
!
            T = ALMR * BETA(I) - BETM * ALFR(I)
            TI = ALMI * BETA(I) - BETM * ALFI(I)
            IF (T .EQ. ZERO .AND. TI .EQ. ZERO) T = EPSB
            Z3 = CMPLX(R,RI,R8) / CMPLX(T,TI,R8)
            BR(I,EN) = REAL(Z3)
            BI(I,EN) = AIMAG(Z3)
  700    CONTINUE
!
  800 CONTINUE
!     ********** END BACK SUBSTITUTION.
!                TRANSFORM TO ORIGINAL COORDINATE SYSTEM.
!                FOR J=N STEP -1 UNTIL 2 DO -- **********
      DO 880 JJ = 2, N
         J = N + 2 - JJ
         M = J - 1
!
         DO 880 I = 1, N
!
            DO 860 K = 1, M
               ZR(I,J) = ZR(I,J) + ZR(I,K) * BR(K,J) - ZI(I,K) * BI(K,J)    
               ZI(I,J) = ZI(I,J) + ZR(I,K) * BI(K,J) + ZI(I,K) * BR(K,J)    
  860       CONTINUE
!
  880 CONTINUE
!     ********** NORMALIZE SO THAT MODULUS OF LARGEST
!                COMPONENT OF EACH VECTOR IS 1 **********
      DO 950 J = 1, N
         T = 0.0_R8
!
         DO 930 I = 1, N
            R = abs(CMPLX(ZR(I,J),ZI(I,J),R8))
            IF (R .GT. T) T = R
  930    CONTINUE
!
         DO 940 I = 1, N
            ZR(I,J) = ZR(I,J) / T
            ZI(I,J) = ZI(I,J) / T
  940    CONTINUE
!
  950 CONTINUE
!
 1001 RETURN
!     ********** LAST CARD OF CQZVEC **********
      END
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
! 15Apr2005 fgtok

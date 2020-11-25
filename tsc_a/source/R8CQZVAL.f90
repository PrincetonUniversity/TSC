!#include "f77_dcomplx.h"
      SUBROUTINE R8CQZVAL(NM,N,AR,AI,BR,BI,EPS1,ALFR,ALFI,BETA,          &  
     &                                       MATZ,ZR,ZI,IERR)
!
!     ------------------------------------------------------------------
!
!
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER I,J,K,L,N,EN,K1,K2,LL,L1,NA,NM,ITS,KM1,LM1,                &  
     &        ENM2,IERR,LOR1,ENORN
      REAL*8 AR(NM,N),AI(NM,N),BR(NM,N),BI(NM,N),                        &  
     &       ALFR(N),ALFI(N),                                            &  
     &       BETA(N),ZR(NM,N),ZI(NM,N)
      REAL*8 R,S,A1,A2,EP,SH,U1,U2,XI,XR,YI,YR,                          &  
     &       ANI,A1I,A33,A34,A43,A44,                                    &  
     &       BNI,B11,B33,B44,SHI,U1I,A33I,A34I,A43I,A44I,B33I,B44I,      &  
     &       EPSA,EPSB,EPS1,ANORM,BNORM,B3344,B3344I
      REAL*8 AREAL
      INTEGER max
      LOGICAL MATZ
      COMPLEX*16 Z3
!
      REAL*8 ZERO,ZONE,ZTWO
!
!
!
!
!     THIS SUBROUTINE IS A COMPLEX ANALOGUE OF STEPS 2 AND 3 OF THE
!     QZ ALGORITHM FOR SOLVING GENERALIZED MATRIX EIGENVALUE PROBLEMS,
!     SIAM J. NUMER. ANAL. 10, 241-256(1973) BY MOLER AND STEWART,
!     AS MODIFIED IN TECHNICAL NOTE NASA TN E-7305(1973) BY WARD.
!
!     THIS SUBROUTINE ACCEPTS A PAIR OF COMPLEX MATRICES, ONE OF THEM
!     IN UPPER HESSENBERG FORM AND THE OTHER IN UPPER TRIANGULAR FORM,
!     THE HESSENBERG MATRIX MUST FURTHER HAVE REAL SUBDIAGONAL ELEMENTS.
!     IT REDUCES THE HESSENBERG MATRIX TO TRIANGULAR FORM USING
!     UNITARY TRANSFORMATIONS WHILE MAINTAINING THE TRIANGULAR FORM
!     OF THE OTHER MATRIX AND FURTHER MAKING ITS DIAGONAL ELEMENTS
!     REAL AND NON-NEGATIVE.  IT THEN RETURNS QUANTITIES WHOSE RATIOS
!     GIVE THE GENERALIZED EIGENVALUES.  IT IS USUALLY PRECEDED BY
!     CQZHES  AND POSSIBLY FOLLOWED BY  CQZVEC.
!
!     ON INPUT-
!
!        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
!          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
!          DIMENSION STATEMENT,
!
!        N IS THE ORDER OF THE MATRICES,
!
!        A=(AR,AI) CONTAINS A COMPLEX UPPER HESSENBERG MATRIX
!          WITH REAL SUBDIAGONAL ELEMENTS,
!
!        B=(BR,BI) CONTAINS A COMPLEX UPPER TRIANGULAR MATRIX,
!
!        EPS1 IS A TOLERANCE USED TO DETERMINE NEGLIGIBLE ELEMENTS.
!          EPS1 = 0.0 (OR NEGATIVE) MAY BE INPUT, IN WHICH CASE AN
!          ELEMENT WILL BE NEGLECTED ONLY IF IT IS LESS THAN ROUNDOFF
!          ERROR TIMES THE NORM OF ITS MATRIX.  IF THE INPUT EPS1 IS
!          POSITIVE, THEN AN ELEMENT WILL BE CONSIDERED NEGLIGIBLE
!          IF IT IS LESS THAN EPS1 TIMES THE NORM OF ITS MATRIX.  A
!          POSITIVE VALUE OF EPS1 MAY RESULT IN FASTER EXECUTION,
!          BUT LESS ACCURATE RESULTS,
!
!        MATZ SHOULD BE SET TO .TRUE. IF THE RIGHT HAND TRANSFORMATIONS
!          ARE TO BE ACCUMULATED FOR LATER USE IN COMPUTING
!          EIGENVECTORS, AND TO .FALSE. OTHERWISE,
!
!        Z=(ZR,ZI) CONTAINS, IF MATZ HAS BEEN SET TO .TRUE., THE
!          TRANSFORMATION MATRIX PRODUCED IN THE REDUCTION
!          BY  CQZHES, IF PERFORMED, OR ELSE THE IDENTITY MATRIX.
!          IF MATZ HAS BEEN SET TO .FALSE., Z IS NOT REFERENCED.
!
!     ON OUTPUT-
!
!        A HAS BEEN REDUCED TO UPPER TRIANGULAR FORM.  THE ELEMENTS
!          BELOW THE MAIN DIAGONAL HAVE BEEN SET TO ZERO,
!
!        B IS STILL IN UPPER TRIANGULAR FORM, ALTHOUGH ITS ELEMENTS
!          HAVE BEEN ALTERED.  IN PARTICULAR, ITS DIAGONAL HAS BEEN SET
!          REAL AND NON-NEGATIVE.  THE LOCATION BR(N,1) IS USED TO
!          STORE EPS1 TIMES THE NORM OF B FOR LATER USE BY  CQZVEC,
!
!        ALFR AND ALFI CONTAIN THE REAL AND IMAGINARY PARTS OF THE
!          DIAGONAL ELEMENTS OF THE TRIANGULARIZED A MATRIX,
!
!        BETA CONTAINS THE REAL NON-NEGATIVE DIAGONAL ELEMENTS OF THE
!          CORRESPONDING B.  THE GENERALIZED EIGENVALUES ARE THEN
!          THE RATIOS ((ALFR+I*ALFI)/BETA),
!
!        Z CONTAINS THE PRODUCT OF THE RIGHT HAND TRANSFORMATIONS
!          (FOR BOTH STEPS) IF MATZ HAS BEEN SET TO .TRUE.,
!
!        IERR IS SET TO
!          ZERO       FOR NORMAL RETURN,
!          J          IF AR(J,J-1) HAS NOT BECOME
!                     ZERO AFTER 50 ITERATIONS.
!
!     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
!     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!
!     ------------------------------------------------------------------
!
      ZTWO = 2.0_R8
      ZONE = 1.0_R8
      ZERO = 0.0_R8
 
      IERR = 0
!     ********** COMPUTE EPSA,EPSB **********
      ANORM = 0.0_R8
      BNORM = 0.0_R8
!
      DO 30 I = 1, N
         ANI = 0.0_R8
         IF (I .NE. 1) ANI = ABS(AR(I,I-1))
         BNI = 0.0_R8
!
         DO 20 J = I, N
            ANI = ANI + ABS(AR(I,J)) + ABS(AI(I,J))
            BNI = BNI + ABS(BR(I,J)) + ABS(BI(I,J))
   20    CONTINUE
!
         IF (ANI .GT. ANORM) ANORM = ANI
         IF (BNI .GT. BNORM) BNORM = BNI
   30 CONTINUE
!
      IF (ANORM .EQ. ZERO) ANORM = 1.0_R8
      IF (BNORM .EQ. ZERO) BNORM = 1.0_R8
      EP = EPS1
      IF (EP .GT. ZERO) GO TO 50
!     ********** COMPUTE ROUNDOFF LEVEL IF EPS1 IS ZERO **********
      EP = 1.0_R8
   40 EP = EP / ZTWO
      IF (ZONE + EP .GT. ZONE) GO TO 40
   50 EPSA = EP * ANORM
      EPSB = EP * BNORM
!     ********** REDUCE A TO TRIANGULAR FORM, WHILE
!                KEEPING B TRIANGULAR **********
      LOR1 = 1
      ENORN = N
      EN = N
!     ********** BEGIN QZ STEP **********
   60 IF (EN .EQ. 0) GO TO 1001
      IF (.NOT. MATZ) ENORN = EN
      ITS = 0
      NA = EN - 1
      ENM2 = NA - 1
!     ********** CHECK FOR CONVERGENCE OR REDUCIBILITY.
!                FOR L=EN STEP -1 UNTIL 1 DO -- **********
   70 DO 80 LL = 1, EN
         LM1 = EN - LL
         L = LM1 + 1
         IF (L .EQ. 1) GO TO 95
         IF (ABS(AR(L,LM1)) .LE. EPSA) GO TO 90
   80 CONTINUE
!
   90 AR(L,LM1) = 0.0_R8
!     ********** SET DIAGONAL ELEMENT AT TOP OF B REAL **********
   95 B11 = abs(CMPLX(BR(L,L),BI(L,L),R8))
      IF (B11     .EQ. ZERO) GO TO 98
      U1 = BR(L,L) / B11
      U1I = BI(L,L) / B11
!
      DO 97 J = L, ENORN
         XI = U1 * AI(L,J) - U1I * AR(L,J)
         AR(L,J) = U1 * AR(L,J) + U1I * AI(L,J)
         AI(L,J) = XI
         XI = U1 * BI(L,J) - U1I * BR(L,J)
         BR(L,J) = U1 * BR(L,J) + U1I * BI(L,J)
         BI(L,J) = XI
   97 CONTINUE
!
      BI(L,L) = 0.0_R8
   98 IF (L .NE. EN) GO TO 100
!     ********** 1-BY-1 BLOCK ISOLATED **********
      ALFR(EN) = AR(EN,EN)
      ALFI(EN) = AI(EN,EN)
      BETA(EN) = B11
      EN = NA
      GO TO 60
!     ********** CHECK FOR SMALL TOP OF B **********
  100 L1 = L + 1
      IF (B11 .GT. EPSB) GO TO 120
      BR(L,L) = 0.0_R8
      S = ABS(AR(L,L)) + ABS(AI(L,L)) + ABS(AR(L1,L))
      U1 = AR(L,L) / S
      U1I = AI(L,L) / S
      U2 = AR(L1,L) / S
      R = SQRT(U1*U1+U1I*U1I+U2*U2)
      U1 = U1 / R
      U1I = U1I / R
      U2 = U2 / R
      AR(L,L) = R * S
      AI(L,L) = 0.0_R8
!
      DO 110 J = L1, ENORN
         XR = AR(L,J)
         XI = AI(L,J)
         YR = AR(L1,J)
         YI = AI(L1,J)
         AR(L,J) = U1 * XR + U1I * XI + U2 * YR
         AI(L,J) = U1 * XI - U1I * XR + U2 * YI
         AR(L1,J) = U1 * YR - U1I * YI - U2 * XR
         AI(L1,J) = U1 * YI + U1I * YR - U2 * XI
         XR = BR(L,J)
         XI = BI(L,J)
         YR = BR(L1,J)
         YI = BI(L1,J)
         BR(L1,J) = U1 * YR - U1I * YI - U2 * XR
         BR(L,J) = U1 * XR + U1I * XI + U2 * YR
         BI(L,J) = U1 * XI - U1I * XR + U2 * YI
         BI(L1,J) = U1 * YI + U1I * YR - U2 * XI
  110 CONTINUE
!
      LM1 = L
      L = L1
      GO TO 90
!     ********** ITERATION STRATEGY **********
  120 IF (ITS .EQ. 50) GO TO 1000
      IF (ITS .EQ. 10) GO TO 135
!     ********** DETERMINE SHIFT **********
      B33 = BR(NA,NA)
      B33I = BI(NA,NA)
      IF (abs(CMPLX(B33,B33I,R8)) .GE. EPSB) GO TO 122
      B33 = EPSB
      B33I = 0.0_R8
  122 B44 = BR(EN,EN)
      B44I = BI(EN,EN)
      IF (abs(CMPLX(B44,B44I,R8)) .GE. EPSB) GO TO 124
      B44 = EPSB
      B44I = 0.0_R8
  124 B3344 = B33 * B44 - B33I * B44I
      B3344I = B33 * B44I + B33I * B44
      A33 = AR(NA,NA) * B44 - AI(NA,NA) * B44I
      A33I = AR(NA,NA) * B44I + AI(NA,NA) * B44
      A34 = AR(NA,EN) * B33 - AI(NA,EN) * B33I                           &  
     &    - AR(NA,NA) * BR(NA,EN) + AI(NA,NA) * BI(NA,EN)
      A34I = AR(NA,EN) * B33I + AI(NA,EN) * B33                          &  
     &     - AR(NA,NA) * BI(NA,EN) - AI(NA,NA) * BR(NA,EN)
      A43 = AR(EN,NA) * B44
      A43I = AR(EN,NA) * B44I
      A44 = AR(EN,EN) * B33 - AI(EN,EN) * B33I - AR(EN,NA) * BR(NA,EN)
      A44I = AR(EN,EN) * B33I + AI(EN,EN) * B33 - AR(EN,NA) * BI(NA,EN)
      SH = A44
      SHI = A44I
      XR = A34 * A43 - A34I * A43I
      XI = A34 * A43I + A34I * A43
      IF (XR .EQ. ZERO .AND. XI .EQ. ZERO) GO TO 140
      YR = (A33 - SH) / 2.0_R8
      YI = (A33I - SHI) / 2.0_R8
      Z3 = sqrt(CMPLX(YR**2-YI**2+XR,2.0_R8*YR*YI+XI,R8))
      U1 = REAL(Z3)
      U1I = AIMAG(Z3)
      IF (YR * U1 + YI * U1I .GE. ZERO) GO TO 125
      U1 = -U1
      U1I = -U1I
  125 Z3 = (CMPLX(SH,SHI,R8) - CMPLX(XR,XI,R8) / CMPLX(YR+U1,YI+U1I,R8))          &  
     &   / CMPLX(B3344,B3344I,R8)
      SH = REAL(Z3)
      SHI = AIMAG(Z3)
      GO TO 140
!     ********** AD HOC SHIFT **********
  135 SH = AR(EN,NA) + AR(NA,ENM2)
      SHI = 0.0_R8
!     ********** DETERMINE ZEROTH COLUMN OF A **********
  140 A1 = AR(L,L) / B11 - SH
      A1I = AI(L,L) / B11 - SHI
      A2 = AR(L1,L) / B11
      ITS = ITS + 1
      IF (.NOT. MATZ) LOR1 = L
!     ********** MAIN LOOP **********
      DO 260 K = L, NA
         K1 = K + 1
         K2 = K + 2
         KM1 = max(K-1,L)
!     ********** ZERO A(K+1,K-1) **********
         IF (K .EQ. L) GO TO 170
         A1 = AR(K,KM1)
         A1I = AI(K,KM1)
         A2 = AR(K1,KM1)
  170    S = ABS(A1) + ABS(A1I) + ABS(A2)
         U1 = A1 / S
         U1I = A1I / S
         U2 = A2 / S
         R = SQRT(U1*U1+U1I*U1I+U2*U2)
         U1 = U1 / R
         U1I = U1I / R
         U2 = U2 / R
!
         DO 180 J = KM1, ENORN
            XR = AR(K,J)
            XI = AI(K,J)
            YR = AR(K1,J)
            YI = AI(K1,J)
            AR(K,J) = U1 * XR + U1I * XI + U2 * YR
            AI(K,J) = U1 * XI - U1I * XR + U2 * YI
            AR(K1,J) = U1 * YR - U1I * YI - U2 * XR
            AI(K1,J) = U1 * YI + U1I * YR - U2 * XI
            XR = BR(K,J)
            XI = BI(K,J)
            YR = BR(K1,J)
            YI = BI(K1,J)
            BR(K,J) = U1 * XR + U1I * XI + U2 * YR
            BI(K,J) = U1 * XI - U1I * XR + U2 * YI
            BR(K1,J) = U1 * YR - U1I * YI - U2 * XR
            BI(K1,J) = U1 * YI + U1I * YR - U2 * XI
  180    CONTINUE
!
         IF (K .EQ. L) GO TO 240
         AI(K,KM1) = 0.0_R8
         AR(K1,KM1) = 0.0_R8
         AI(K1,KM1) = 0.0_R8
!     ********** ZERO B(K+1,K) **********
  240    S = ABS(BR(K1,K1)) + ABS(BI(K1,K1)) + ABS(BR(K1,K))
         U1 = BR(K1,K1) / S
         U1I = BI(K1,K1) / S
         U2 = BR(K1,K) / S
         R = SQRT(U1*U1+U1I*U1I+U2*U2)
         U1 = U1 / R
         U1I = U1I / R
         U2 = U2 / R
         IF (K .EQ. NA) GO TO 245
         XR = AR(K2,K1)
         AR(K2,K1) = U1 * XR
         AI(K2,K1) = -U1I * XR
         AR(K2,K) = -U2 * XR
!
  245    DO 250 I = LOR1, K1
            XR = AR(I,K1)
            XI = AI(I,K1)
            YR = AR(I,K)
            YI = AI(I,K)
            AR(I,K1) = U1 * XR + U1I * XI + U2 * YR
            AI(I,K1) = U1 * XI - U1I * XR + U2 * YI
            AR(I,K) = U1 * YR - U1I * YI - U2 * XR
            AI(I,K) = U1 * YI + U1I * YR - U2 * XI
            XR = BR(I,K1)
            XI = BI(I,K1)
            YR = BR(I,K)
            YI = BI(I,K)
            BR(I,K1) = U1 * XR + U1I * XI + U2 * YR
            BI(I,K1) = U1 * XI - U1I * XR + U2 * YI
            BR(I,K) = U1 * YR - U1I * YI - U2 * XR
            BI(I,K) = U1 * YI + U1I * YR - U2 * XI
  250    CONTINUE
!
         BI(K1,K1) = 0.0_R8
         BR(K1,K) = 0.0_R8
         BI(K1,K) = 0.0_R8
         IF (.NOT. MATZ) GO TO 260
!
         DO 255 I = 1, N
            XR = ZR(I,K1)
            XI = ZI(I,K1)
            YR = ZR(I,K)
            YI = ZI(I,K)
            ZR(I,K1) = U1 * XR + U1I * XI + U2 * YR
            ZI(I,K1) = U1 * XI - U1I * XR + U2 * YI
            ZR(I,K) = U1 * YR - U1I * YI - U2 * XR
            ZI(I,K) = U1 * YI + U1I * YR - U2 * XI
  255    CONTINUE
!
  260 CONTINUE
!     ********** SET LAST A SUBDIAGONAL REAL AND END QZ STEP **********
      IF (AI(EN,NA) .EQ. ZERO) GO TO 70
      R = abs(CMPLX(AR(EN,NA),AI(EN,NA),R8))
      U1 = AR(EN,NA) / R
      U1I = AI(EN,NA) / R
      AR(EN,NA) = R
      AI(EN,NA) = 0.0_R8
!
      DO 270 J = EN, ENORN
         XI = U1 * AI(EN,J) - U1I * AR(EN,J)
         AR(EN,J) = U1 * AR(EN,J) + U1I * AI(EN,J)
         AI(EN,J) = XI
         XI = U1 * BI(EN,J) - U1I * BR(EN,J)
         BR(EN,J) = U1 * BR(EN,J) + U1I * BI(EN,J)
         BI(EN,J) = XI
  270 CONTINUE
!
      GO TO 70
!     ********** SET ERROR -- BOTTOM SUBDIAGONAL ELEMENT HAS NOT
!                BECOME NEGLIGIBLE AFTER 50 ITERATIONS **********
 1000 IERR = EN
!     ********** SAVE EPSB FOR USE BY CQZVEC **********
 1001 IF (N .GT. 1) BR(N,1) = EPSB
      RETURN
!     ********** LAST CARD OF CQZVAL **********
      END
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
! 15Apr2005 fgtok

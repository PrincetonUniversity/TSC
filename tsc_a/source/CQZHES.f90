      SUBROUTINE CQZHES(NM,N,AR,AI,BR,BI,MATZ,ZR,ZI)
 
!     ------------------------------------------------------------------
!
!
!============
! idecl:  explicitize implicit REAL declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      REAL*8 zero,zone
!============
      INTEGER I,J,K,L,N,K1,LB,L1,NM,NK1,NM1
      REAL*8 AR(NM,N),AI(NM,N),BR(NM,N),BI(NM,N),ZR(NM,N),ZI(NM,N)
      REAL*8 R,S,T,TI,U1,U2,XI,XR,YI,YR,RHO,U1I
      REAL*8, intrinsic ::  SQRT,ABS
      LOGICAL MATZ
      COMPLEX*16 CMPLX
!
!     THIS SUBROUTINE IS A COMPLEX ANALOGUE OF THE FIRST STEP OF THE
!     QZ ALGORITHM FOR SOLVING GENERALIZED MATRIX EIGENVALUE PROBLEMS,
!     SIAM J. NUMER. ANAL. 10, 241-256(1973) BY MOLER AND STEWART.
!
!     THIS SUBROUTINE ACCEPTS A PAIR OF COMPLEX GENERAL MATRICES AND
!     REDUCES ONE OF THEM TO UPPER HESSENBERG FORM WITH REAL (AND NON-
!     NEGATIVE) SUBDIAGONAL ELEMENTS AND THE OTHER TO UPPER TRIANGULAR
!     FORM USING UNITARY TRANSFORMATIONS.  IT IS USUALLY FOLLOWED BY
!     CQZVAL  AND POSSIBLY  CQZVEC.
!
!     ON INPUT-
!
!        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
!          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
!          DIMENSION STATEMENT,
!
!        N IS THE ORDER OF THE MATRICES,
!
!        A=(AR,AI) CONTAINS A COMPLEX GENERAL MATRIX,
!
!        B=(BR,BI) CONTAINS A COMPLEX GENERAL MATRIX,
!
!        MATZ SHOULD BE SET TO .TRUE. IF THE RIGHT HAND TRANSFORMATIONS
!          ARE TO BE ACCUMULATED FOR LATER USE IN COMPUTING
!          EIGENVECTORS, AND TO .FALSE. OTHERWISE.
!
!     ON OUTPUT-
!
!        A HAS BEEN REDUCED TO UPPER HESSENBERG FORM.  THE ELEMENTS
!          BELOW THE FIRST SUBDIAGONAL HAVE BEEN SET TO ZERO, AND THE
!          SUBDIAGONAL ELEMENTS HAVE BEEN MADE REAL (AND NON-NEGATIVE),
!
!        B HAS BEEN REDUCED TO UPPER TRIANGULAR FORM.  THE ELEMENTS
!          BELOW THE MAIN DIAGONAL HAVE BEEN SET TO ZERO,
!
!        Z=(ZR,ZI) CONTAINS THE PRODUCT OF THE RIGHT HAND
!          TRANSFORMATIONS IF MATZ HAS BEEN SET TO .TRUE.
!          OTHERWISE, Z IS NOT REFERENCED.
!
!     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
!     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!
!     ------------------------------------------------------------------
!
!     ********** INITIALIZE Z **********
 
      ZERO = 0.0_R8
      ZONE = 1.0_R8
      IF (.NOT. MATZ) GO TO 10
!
      DO 3 I = 1, N
!
         DO 2 J = 1, N
            ZR(I,J) = 0.0_R8
            ZI(I,J) = 0.0_R8
    2    CONTINUE
!
         ZR(I,I) = 1.0_R8
    3 CONTINUE
!     ********** REDUCE B TO UPPER TRIANGULAR FORM WITH
!                TEMPORARILY REAL DIAGONAL ELEMENTS **********
   10 IF (N .LE. 1) GO TO 170
      NM1 = N - 1
!
      DO 100 L = 1, NM1
         L1 = L + 1
         S = 0.0_R8
!
         DO 20 I = L, N
            S = S + ABS(BR(I,L)) + ABS(BI(I,L))
   20    CONTINUE
!
         IF (S .EQ. ZERO) GO TO 100
         RHO = 0.0_R8
!
         DO 25 I = L, N
            BR(I,L) = BR(I,L) / S
            BI(I,L) = BI(I,L) / S
            RHO = RHO + BR(I,L)**2 + BI(I,L)**2
   25    CONTINUE
!
         R = SQRT(RHO)
         XR = abs(CMPLX(BR(L,L),BI(L,L),R8))
         IF (XR .EQ. ZERO) GO TO 27
         RHO = RHO + XR * R
         U1 = -BR(L,L) / XR
         U1I = -BI(L,L) / XR
         YR = R / XR + 1.0_R8
         BR(L,L) = YR * BR(L,L)
         BI(L,L) = YR * BI(L,L)
         GO TO 28
!
   27    BR(L,L) = R
         U1 = -1.0_R8
         U1I = 0.0_R8
!
   28    DO 50 J = L1, N
            T = 0.0_R8
            TI = 0.0_R8
!
            DO 30 I = L, N
               T = T + BR(I,L) * BR(I,J) + BI(I,L) * BI(I,J)
               TI = TI + BR(I,L) * BI(I,J) - BI(I,L) * BR(I,J)
   30       CONTINUE
!
            T = T / RHO
            TI = TI / RHO
!
            DO 40 I = L, N
               BR(I,J) = BR(I,J) - T * BR(I,L) + TI * BI(I,L)
               BI(I,J) = BI(I,J) - T * BI(I,L) - TI * BR(I,L)
   40       CONTINUE
!
            XI = U1 * BI(L,J) - U1I * BR(L,J)
            BR(L,J) = U1 * BR(L,J) + U1I * BI(L,J)
            BI(L,J) = XI
   50    CONTINUE
!
         DO 80 J = 1, N
            T = 0.0_R8
            TI = 0.0_R8
!
            DO 60 I = L, N
               T = T + BR(I,L) * AR(I,J) + BI(I,L) * AI(I,J)
               TI = TI + BR(I,L) * AI(I,J) - BI(I,L) * AR(I,J)
   60       CONTINUE
!
            T = T / RHO
            TI = TI / RHO
!
            DO 70 I = L, N
               AR(I,J) = AR(I,J) - T * BR(I,L) + TI * BI(I,L)
               AI(I,J) = AI(I,J) - T * BI(I,L) - TI * BR(I,L)
   70       CONTINUE
!
            XI = U1 * AI(L,J) - U1I * AR(L,J)
            AR(L,J) = U1 * AR(L,J) + U1I * AI(L,J)
            AI(L,J) = XI
   80    CONTINUE
!
         BR(L,L) = R * S
         BI(L,L) = 0.0_R8
!
         DO 90 I = L1, N
            BR(I,L) = 0.0_R8
            BI(I,L) = 0.0_R8
   90    CONTINUE
!
  100 CONTINUE
!     ********** REDUCE A TO UPPER HESSENBERG FORM WITH REAL SUBDIAGONAL
!                ELEMENTS, WHILE KEEPING B TRIANGULAR **********
      DO 160 K = 1, NM1
         K1 = K + 1
!     ********** SET BOTTOM ELEMENT IN K-TH COLUMN OF A REAL **********
         IF (AI(N,K) .EQ. ZERO) GO TO 105
         R = abs(CMPLX(AR(N,K),AI(N,K),R8))
         U1 = AR(N,K) / R
         U1I = AI(N,K) / R
         AR(N,K) = R
         AI(N,K) = 0.0_R8
!
         DO 103 J = K1, N
            XI = U1 * AI(N,J) - U1I * AR(N,J)
            AR(N,J) = U1 * AR(N,J) + U1I * AI(N,J)
            AI(N,J) = XI
  103    CONTINUE
!
         XI = U1 * BI(N,N) - U1I * BR(N,N)
         BR(N,N) = U1 * BR(N,N) + U1I * BI(N,N)
         BI(N,N) = XI
  105    IF (K .EQ. NM1) GO TO 170
         NK1 = NM1 - K
!     ********** FOR L=N-1 STEP -1 UNTIL K+1 DO -- **********
         DO 150 LB = 1, NK1
            L = N - LB
            L1 = L + 1
!     ********** ZERO A(L+1,K) **********
            S = ABS(AR(L,K)) + ABS(AI(L,K)) + AR(L1,K)
            IF (S .EQ. ZERO) GO TO 150
            U1 = AR(L,K) / S
            U1I = AI(L,K) / S
            U2 = AR(L1,K) / S
            R = SQRT(U1*U1+U1I*U1I+U2*U2)
            U1 = U1 / R
            U1I = U1I / R
            U2 = U2 / R
            AR(L,K) = R * S
            AI(L,K) = 0.0_R8
            AR(L1,K) = 0.0_R8
!
            DO 110 J = K1, N
               XR = AR(L,J)
               XI = AI(L,J)
               YR = AR(L1,J)
               YI = AI(L1,J)
               AR(L,J) = U1 * XR + U1I * XI + U2 * YR
               AI(L,J) = U1 * XI - U1I * XR + U2 * YI
               AR(L1,J) = U1 * YR - U1I * YI - U2 * XR
               AI(L1,J) = U1 * YI + U1I * YR - U2 * XI
  110       CONTINUE
!
            XR = BR(L,L)
            BR(L,L) = U1 * XR
            BI(L,L) = -U1I * XR
            BR(L1,L) = -U2 * XR
!
            DO 120 J = L1, N
               XR = BR(L,J)
               XI = BI(L,J)
               YR = BR(L1,J)
               YI = BI(L1,J)
               BR(L,J) = U1 * XR + U1I * XI + U2 * YR
               BI(L,J) = U1 * XI - U1I * XR + U2 * YI
               BR(L1,J) = U1 * YR - U1I * YI - U2 * XR
               BI(L1,J) = U1 * YI + U1I * YR - U2 * XI
  120       CONTINUE
!     ********** ZERO B(L+1,L) **********
            S = ABS(BR(L1,L1)) + ABS(BI(L1,L1)) + ABS(BR(L1,L))
            IF (S .EQ. ZERO) GO TO 150
            U1 = BR(L1,L1) / S
            U1I = BI(L1,L1) / S
            U2 = BR(L1,L) / S
            R = SQRT(U1*U1+U1I*U1I+U2*U2)
            U1 = U1 / R
            U1I = U1I / R
            U2 = U2 / R
            BR(L1,L1) = R * S
            BI(L1,L1) = 0.0_R8
            BR(L1,L) = 0.0_R8
!
            DO 130 I = 1, L
               XR = BR(I,L1)
               XI = BI(I,L1)
               YR = BR(I,L)
               YI = BI(I,L)
               BR(I,L1) = U1 * XR + U1I * XI + U2 * YR
               BI(I,L1) = U1 * XI - U1I * XR + U2 * YI
               BR(I,L) = U1 * YR - U1I * YI - U2 * XR
               BI(I,L) = U1 * YI + U1I * YR - U2 * XI
  130       CONTINUE
!
            DO 140 I = 1, N
               XR = AR(I,L1)
               XI = AI(I,L1)
               YR = AR(I,L)
               YI = AI(I,L)
               AR(I,L1) = U1 * XR + U1I * XI + U2 * YR
               AI(I,L1) = U1 * XI - U1I * XR + U2 * YI
               AR(I,L) = U1 * YR - U1I * YI - U2 * XR
               AI(I,L) = U1 * YI + U1I * YR - U2 * XI
  140       CONTINUE
!
            IF (.NOT. MATZ) GO TO 150
!
            DO 145 I = 1, N
               XR = ZR(I,L1)
               XI = ZI(I,L1)
               YR = ZR(I,L)
               YI = ZI(I,L)
               ZR(I,L1) = U1 * XR + U1I * XI + U2 * YR
               ZI(I,L1) = U1 * XI - U1I * XR + U2 * YI
               ZR(I,L) = U1 * YR - U1I * YI - U2 * XR
               ZI(I,L) = U1 * YI + U1I * YR - U2 * XI
  145       CONTINUE
!
  150    CONTINUE
!
  160 CONTINUE
!
  170 RETURN
!     ********** LAST CARD OF CQZHES **********
      END
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

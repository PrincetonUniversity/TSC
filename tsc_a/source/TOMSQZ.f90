      SUBROUTINE TOMSQZ(N,NA,AR,AI,BR,BI,ALFR,ALFI,BETA,ZVR,ZVI,IFAIL)
!|
!| Note that
!| $$ \omega_{De} = 2 k_\perp T_e / e B R = 2 k_\perp \rho_s c_s / R $$
!| $$ D \propto \gamma / k^2
!|   \propto \rho_s^2 \omega_{De}
!|     \frac{ \gamma / \omega_{De} }{ k^2 \rho_s^2 }. $$
!| Note, $ \omega_{*e} = k_\perp \rho_s c_s g_{ne} $,
!| $ \rho_s = c_s / \omega_{ci} $, $ c_s = \sqrt{ 2 T_e / m_i} $,
!| and $ \omega_{ci} = Z_i e B / m_i $,
!| $ \omega_{De} R = k_y \rho_s c_s $,
!| $v_i = \sqrt{T_i / m_i}$, and $ g_n = - R \Partial{n}{r} / n $
!| in the above normalizations.
!| Note that all the diffusivities in this routine are normalized by
!| $ \omega_{De} / k_y^2 $,
!| convective velocities are normalized by $ \omega_{De} / k_y $,
!| and all the frequencies are normalized by $ \omega_{De} $.
!|
!|
!| \begin{thebibliography}{99}
!|
!| \bibitem{bate98a}
!| Glenn Bateman, Arnold~H. Kritz, Jon~E. Kinsey, Aaron~J. Redd, and Jan Weiland,
!| ``Predicting temperature and density profiles in tokamaks,''
!| {\em Physics of Plasmas,} {\bf 5} (1998) 1793--1799.
!|
!| \bibitem{weil92a} J. Weiland,
!| ``Low Frequency Modes Associated with Drift Motions in
!| Inhomogeneous Plasmas,''
!| report CTH--IEFT/PP-1992-17,
!| Institute for Electromagnetic Field Theory and Plasma Physics,
!| Chalmers University of Technology,
!| G\"{o}teborg, Sweden.
!|
!| \bibitem{froj92a} M. Fr\"{o}jdh, M. Liljestr\"{o}m, H. Nordman,
!| ``Impurity effects on $\eta_i$ mode stability and transport,''
!| Nuclear Fusion {\bf 32} (1992) 419--428.
!|
!| \bibitem{nord92a} H. Nordman and J. Weiland, ``Comments on
!| `Ion-temperature-gradient-driven
!| transport in a density modification experiment on the tokamak fusion test
!| reactor [Phys. Fluids {\bf B4} (1992) 953]' ''.
!|
!| \bibitem{weil92b} J. Weiland,
!| ``Nonlinear effects in velocity space and drift wave
!| transport in tokamaks,'' Phys. Fluids {\bf B4} (1992) 1388--1390.
!|
!| \bibitem{weil92c} J. Weiland and H. Nordman, ``Drift wave model for inward
!| energy transport in tokamak plasmas,'' Institute for Electromagnetic Field
!| Theory and Plasma Physics, Gothenburg, Sweden, (1992) CTH-IEFT/PP-1992-13 ISSN.
!|
!| \bibitem{weil92d} J.Weiland and A. Hirose, ``Electromagnetic and kinetic
!| effects on the ion temperature gradient mode,'' Nucl. Fusion {\bf 32} (1992)
!| 151--155.
!|
!| \bibitem{nord91a} H. Nordman and J. Weiland, ``The concept of marginal
!| stability and recent experimental results from the TFTR tokamak,'' Institute
!| for Electromagnetic Field Theory and Plasma Physics, Gothenburg, Sweden, (1991)
!| CTH-IEFT/PP-1991-26 ISSN.
!|
!| \bibitem{weil91a} J. Weiland and H. Nordman, ``Enhanced confinement regimes in
!| transport code simulations of toroidal drift wave transport,'' Nucl. Fusion
!| {\bf 31} (1991) 390--394.
!|
!| \bibitem{nord90a} H. Nordman, J. Weiland, and A. Jarmen, ``Simulation of
!| toroidal drift mode turbulence driven by temperature gradients and electron
!| trapping,'' Nucl. Fusion {\bf 30} (1990) 983--996.
!|
!| \bibitem{weil89a} J. Weiland, A.B. Jarm\'{e}n, and H. Nordman, ``Diffusive
!| particle and heat pinch effects in toroidal plasmas,'' Nucl. Fusion {\bf 29}
!| (1989) 1810--1814.
!|
!| \bibitem{nord89a} H. Nordman and J. Weiland, ``Transport due to toroidal
!| $\eta_i$ mode turbulence in tokamaks,'' Nucl. Fusion {\bf 29} (1989) 251--263.
!|
!| \bibitem{ande88a} P. Andersson and J. Weiland,
!| ``A fully toroidal fluid analysis
!| of the magnetohydrodynamic ballooning mode branch in tokamaks,''
!| Phys. Fluids {\bf 31} (1988) 359--365.
!|
!| \bibitem{jarm87a} A. Jarm\'{e}n, P. Andersson, and J. Weiland,
!| ``Fully toroidal ion temperature gradient driven drift modes,''
!| Nucl. Fusion {\bf 27} (1987) 941--949.
!| \end{thebibliography}
!|
!| \end{document}
!|
!|
!|
!|
!|
!|
!|
!|
!|
!|
!|
!     SUBROUTINE TOMSQZ(N,NA,AR,AI,BR,BI,ALFR,ALFI,BETA,ZVR,ZVI,IFAIL)
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
 
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER ifail
!============
      INTEGER N, NA
      REAL*8 AR(N,NA),AI(N,NA),BR(N,NA),BI(N,NA)
 
!-----------------------------------------------------------------------
! DECALRATIONS FOR OUTPUT VARIABLES
!-----------------------------------------------------------------------
 
      REAL*8 ALFR(N),ALFI(N),BETA(N)
      REAL*8 ZVR(N,NA), ZVI(N,NA)
 
!-----------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------
 
      LOGICAL WANTX
      REAL*8 EPS1
 
!-----------------------------------------------------------------------
! START OF ACTUAL CODING
!-----------------------------------------------------------------------
 
      WANTX = .TRUE.
      EPS1  = -0.0_R8
 
      CALL CQZHES(N,NA,AR,AI,BR,BI,WANTX,ZVR,ZVI)
      CALL CQZVAL(N,NA,AR,AI,BR,BI,EPS1,ALFR,ALFI,BETA,WANTX,            &  
     &            ZVR,ZVI,IFAIL)
      CALL CQZVEC(N,NA,AR,AI,BR,BI,ALFR,ALFI,BETA,ZVR,ZVI)
      RETURN
      END
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

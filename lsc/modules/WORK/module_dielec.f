C     dielec.inc:contains mostly dielectric properties of plasma ------|
C                                                                      |
C                                                                      |
c     NEQS is the number of equations integrated for the ray:
c     r; z; phi; k_r; k_z; k_phi=n_phi; and time.
c     n_phi is conserved, so this is not needed.
c     time is used in calculating E_z, thru the time spent is a zone.
c
!============
! idecl:  explicitize implicit INTEGER declarations:
!============
      INTEGER NEQS, NEQSP1
      PARAMETER (NEQS=7,NEQSP1=8)
      REAL*8    f(NEQSP1),f1(NEQS),f2(NEQS),f3(NEQS),
     ^        y(NEQSP1),y1(NEQS),y2(NEQS),y3(NEQS)
      COMMON /IGRAT1/ f,y
      COMMON /IGRAT2/ f1,f2,f3,y1,y2,y3
 
      INTEGER IDM
      PARAMETER (IDM = NZONDIM)
      INTEGER iplt
      REAL*8    d1,d2,d4
      REAL*8    denom,DKpar,DKper,wdDdw,dDdkABS,epQL, epsL, epsz
      REAL*8    yrr(IDM),yzz(IDM),ypp(IDM),
     ^        ypr(IDM),ypl(IDM),ypw(IDM),
     ^        yed(IDM),yps(IDM),
     ^        y1old   ,y2old   ,y3old
 
      REAL*8    Eper , Epar , Exy  , Aion , Aelc , OmEC2,
     ^        EparI, EparIK, cEparIK
 
      REAL*8    D11er, D33er, D12er,
     ^        D11ar, D33ar, D12ar,
     ^        D11w0, D33w0, D12w0
 
      COMMON /DNORM / d1,d2,d4
      COMMON /DAMP  /
     ^        denom,DKpar,DKper,wdDdw,dDdkABS,epQL, epsL, epsz
      COMMON/PltBlk/
     ^        yrr,yzz,ypp,
     ^        ypr,ypl,ypw,
     ^        yed,yps,
     ^        y1old   ,y2old   ,y3old   ,
     ^        iplt
 
c      yrr:r       yzz:z      ypp:\phi    yps: sqrt( psi_{norm} )
c      ypr:k_\perp ypl:n_\parallel
c      ypw:pwr(s)  power remaining
c      yed:electon damping term
C      \Im{k_\perp} / \Re{k_\perp}
 
      Common /EPS1 /
     ^        Eper , Epar , Exy  , Aion , Aelc , OmEC2,
     ^        Epari, EparIK, cEparIK
      Common /EPS2 /
     ^        D11er, D33er, D12er,
     ^        D11ar, D33ar, D12ar,
     ^        D11w0, D33w0, D12w0
 
C     Blocks EPS give dielectric tensor elements and their derivatives.
C     Using Stix notation of S (Kper) P(Kpar) D(Kxy) R L and
C     a = (wp/w)^2   b = (wce/w)^1
C     OmEC2  =  1. - b^2
C     we have
C     Lower Hybrid Case           Electron Cyclotron Case
C Eper       S                             S (1-b^2)
C Epar       P                             P
C Exy        D
C Erl                                     R L(1-b^2)
C Aion    Ion mode conversion
C Aelc    Electron term                  Electron mode conversion
C     The derivatives of the tensor elements are given in the
C     D terms.  Usage differs for ecrh and lhrh.
C
C                       K               K               K
C   ECRH                 PERP            RL              PAR
C
C  (K-PERP)**2          D11ER           DRLER           D33ER
 
C  (K-PAR )**2          D11AR           DRLAR           D33AR
 
C  (OMEGA )**2          D11w0           DRLw0           D33w0
C
C                       K               K               K
C    LHRH                PERP            XY              PAR
C
C  (K-PERP)**2          D11ER           D12ER           D33ER
 
C  (K-PAR )**2          D11AR           D12AR           D33AR
 
C  (OMEGA )**2          D11w0           D12w0           D33w0
C
C     iplt      index of plot array
C     y         function   array for integration
c             y(1   2   3   4   5   6   7   8)
c               r   z   phi kr  kz  n   wt  s
C     f         derivative of y by s, array for integration
C     f1,       derivative
C     f2,               array
C     f3,                     at old points
C     y1,       function
C     y2,                 array
C     y3                        at old points
C     d1,       terms
C     d2,             of dispersion
C     d4,                          relation
C     DKpar ==  \p D/\p k_{\parallel}^2
C     DKper ==  \p D/\p k_{\perp}^2
C     wdDdw ==  \omega \p D/\p \omega
C     dDdkABS == | \p D/{\p {\bf k}} |
C             == wdDdw * dsdwt
C                which is of interest in the damping calculation:
C     d {ln P}/ds = - 2 {\p D/\p Epar}/{dDdkABS} x Im{Epar}
C     where Im{Epar} =  PI ( \omega_p/k_{\par}v_{\par})^2
C                       (-v_{\par}^2 \p f_o/\p v) @ v=\omega/k_{\par}
C     if Maxwellian
C           Im{Epar} = PI^(1/2)(\omega_p/\omega)^2 2 \zeta^3 exp(-\zeta^2)
C           \zeta = \omega/(2^{1/2} k_{\par} v_t)^2
C     dlnPdsX == d ln(P)/ds \cdot ds if maXwellian, or at maXimum
C                                       without quasilinear burnthrough.
C     dlnPdsX = -2 (\p D/\p Epar) / dDdkABS * Im{Epar} \times ds
C                                       where Epar is K_{zz} is K_{33}
c     dlnPdsK == d ln(P)/ds \cdot ds per f_e^{\prime}, or Kernel of
c                                       the damping decrement.
c     EparI   == Im{K_{zz}} assuming Maxwellian distribution
c     EparIK  == Im{K_{zz}} for unit value of  \p f_e/\p v_{\parallel}
c                                        Note: f_e contains electron density
c                                        and v_{\parallel} is normalized to c
c
C                                                                      |
C                                                                      |
C     dielec.inc:contains mostly dielectric properties of plasma ------|
! 24May2005 fgtok -s r8_precision.sub "r8con_incl.csh conversion"

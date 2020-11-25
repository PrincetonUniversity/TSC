      subroutine rloads1d
!======================================================================
 
!  new version calculating separately the stabilizing forces of
!  passive conductors and vacuum vessel (a. krause 23.10.90)
!     -----------------
!********************************************************************
!********************************************************************
      USE CLINAM
      USE RADTAB
      USE SAPROP
      USE SPECIE
      USE WALLCL
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER j,nimp,nchrgs,n
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 zjlkev,factrx,temin,anemin,ademin,fac4,dadeo,rsum
      REAL*8 sradmax
!============
      data zjlkev / 1.6021E-16_R8/
      data factrx / 0.9_R8/
!
!  * * evaluate plasma radiant/ionization losses
!
      temin = max(0.2_R8,smallt*usdd/usdh)
      anemin = fracn0*r0*udsd
      do 10 j = 2,npsit
      if(kcycle.le.0) sradiono(j)=0._R8
!
      ademin = temin*ane(j)*vpg(j)/udsh
      fac4 = (2._R8/3._R8)*max(dtsf/usdt,0.5E-6_R8)*vpg(j)*usdp
      dadeo = (ade(j)-adeo(j))+fac4*sradion(j)
      if(ade(j).le.ademin) dadeo = 0.0_R8
!
      rsum = 0.0_R8
      do 45 nimp = 1,pimp
      nchrgs = nchrgsr(nimp)
      do 40 n = 1,nchrgs
      rsum = rsum+rad(nimp,n,j)*nq(n,nimp,j)/vp(j)
   40 continue
   45 continue
!
      sradion(j) = sradiono(j) + 0.2_R8*(rsum-sradiono(j))
      sradiono(j) = sradion(j)
      sradmax = factrx*((ade(j)-ademin)+dadeo)/fac4
      if(sradion(j).gt.0.0_R8) go to 5
      sradion(j) = 0.0_R8
      go to 10
    5 continue
      if(sradion(j).le.sradmax) go to 10
!
!
      sradion(j) = max(sradmax,0.0_R8)
!
   10 continue
!
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

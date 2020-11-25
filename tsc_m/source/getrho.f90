      subroutine getrho(itype,tcu,rho,gamma)
!.....9.31 sprop
!.....
!.....    itype - conductor type. Presently (but easily extended):
!.....            1 - OFHC Copper
!.....            2 - AL25 (Glidcop)
!.....            3 - BeCu
!.....    tcu   - temperature (K)
!.....    rho   - resistivity (Ohm-m)
!.....    gamma - density (kg/m**3)
!.....
!.....                      REFERENCE:
!. OFHC  Handbook On Materials For Superconducting Machinery, Battelle, 1977
!. AL25  GRD-47 Peter Winn, MIT
!. BeCu  Private Communication with W.Riersen
!.....
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER itype
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 tcu,rho,gamma,rho0cu,acu,bcu,rho0al25,aal25,bal25
      REAL*8 rho0becu,abecu,bbecu,rho0,arho,brho,tuse
!============
      data rho0cu  /1.724E-8_R8/,acu  /-0.2013_R8/,bcu  /0.0041_R8/
      data rho0al25/1.000_R8/,aal25/-2.5023E-9_R8/,                       &
     &             bal25/7.627907E-11_R8/
      data rho0becu/1.000_R8/,abecu/ 8.170E-9_R8/,bbecu/7.91E-11_R8/
!.....     Error Trap
      if(itype.lt.1.or.itype.gt.3) then
!........................NP ?????
      end if
      if(itype.eq.1) then
        rho0 = rho0cu
        arho = acu
        brho = bcu
        gamma= 8950._R8
      end if
      if(itype.eq.2) then
        rho0 = rho0al25
        arho = aal25
        brho = bal25
        gamma= 8830._R8
      end if
      if(itype.eq.3) then
        rho0 = rho0becu
        arho = abecu
        brho = bbecu
        gamma= 8950._R8
      end if
      tuse = tcu
      if(tuse.lt.70._R8) tuse = 70._R8
      rho = rho0*(arho+brho*tuse)
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

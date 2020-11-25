      subroutine  getemp (i, j, tempval)
!*************************************************************************
!
!     returns temperature in ev at (i,j)
!
      USE CLINAM
      USE SAPROP
      USE SCR15
      USE WALLCL
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!              Set temp to vacuum value
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER j,i
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 tempval,rval,dens,tmid,allam,thalo
!============
      tempval = tevv
      if (  i.lt.2 .or. i.gt.nxp                                         &  
     &    .or.j.lt.2 .or. j.gt.nzp)   return
      if (abs(etay(i,j)) .lt. 1.E-30_R8)   return
!
!.....            check if outside of vacuum vessel
      if(iexv(i,j).eq.1)   return
!.....            check if on wrong side of separatrix
      if(igone.eq.1)   go to 499
      if(iexs(i,j).eq.1) go to 532
      go to 498
  532 continue
      if(abs(psi(i,j)-psilim) .lt. abs(phalo-psilim)) go to 499
         return
  498 if(psi(i,j).gt.phalo)   return
      if(psi(i,j).gt.psilim) go to 499
!
!.....            point lies in plasma
      rval = .25_R8*(roj(i,j)+roj(i+1,j)+roj(i,j+1)+roj(i+1,j+1))
      dens = udsd*rval
!           Test dens  26 Apr 1993
      tmid = tevv
      if (dens.gt.1.E-30_R8)                                             &  
     & tmid=.25_R8*(.5_R8*udsh)/dens*(pr(i,j)+pr(i+1,j)+pr(i,j+1)+pr(i+  &  
     & 1,j+1))
      if(tmid.lt.tevv)      tmid = tevv
      if(tmid.gt.acoef(71)) tmid = acoef(71)
      allam = 24._R8-log(1.E-3_R8*sqrt(dens)/tmid)
      tempval = (0.5_R8* 1.03E-4_R8*allam * usdr * zeff)/etay(i,j)
      if (tempval.gt.1.E-30_R8)   tempval = tempval**(0.6666666667_R8)
      return
!              Halo point
  499    continue
      if (acoef(99).gt.0.0_R8.and. phalo.gt.psilim)   then
      thalo =acoef(98)-(acoef(98)-tevv)*(psi(i,j)-psilim)/(phalo-psilim)    
        if (thalo.gt.acoef(98))   thalo = acoef(98)
        if (thalo.lt.tevv)        thalo = tevv
                     tempval = thalo
      else
      tempval = acoef(98)
                                                    endif
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

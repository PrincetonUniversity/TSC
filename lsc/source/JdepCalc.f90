!
!     ------------------------------------------------------------------
!
      SUBROUTINE JdepCalc
      USE Jrf
      USE params
      USE PlPr
      USE ProfBody
      USE RayBins
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER iabs
!============
      INTEGER ip, ipLastRuna, iGotRuna, nGotRuna, iFillJray
      INTEGER nrLines
      CHARACTER*70 ErrMsg
      DATA nrLines / 0 /
      REAL*8                                                             &  
     &        ZERO, ONE
      DATA    ZERO, ONE/                                                 &  
     &         0.0_R8, 1.0_R8/
 
 
!
!     compute rf driven current given power dissipation vs vpar, psi
!
      call GetEdc(ZERO)
!     call GetEdc(0.0)
!                                       Fill array EdcAry
      call BrodCast(NPSIDIM * NVELDIM, Jray,    0._R8)
      call VmaxCalc
      nGotRuna = 0
      do 10 ip = 1, npsi
        call jnorm(ip)
!                                       set up normalization,
!                                       tabulate u lookup values
!
!     .
!     .                                 This is based on figure 2 of the
!     .                                 Karney Fisch paper.  It shows that
!     .                                 if you are trying to drive current in
!     .                                 cooperation with the electric field
!     .                                 at v/vrun you have had it.  If trying
!     .                                 to drive against the electric field
!     .                                 you can go up to v/vrun of 2 before
!     .                                 runaways take over.
        if (EdcAry(ip) .ge. 0._R8) then
          vnormPos(ip) = VparMaxP(ip)/vnorm /(1.0_R8)
          vnormNeg(ip) = VparMaxN(ip)/vnorm /(2.0_R8)
        else
          vnormPos(ip) = VparMaxP(ip)/vnorm /(2.0_R8)
          vnormNeg(ip) = VparMaxN(ip)/vnorm /(1.0_R8)
        endif
 
        vnormNOK(ip) = max(int(vnormPos(ip)),                            &  
     &                 iabs(int(vnormNeg(ip))))
 
!
!     .                                 compute jrf
        iFillJray = 1
        call mkj(ip,js(ip), iGotRuna, iFillJray)
!     .                                 Report if a runaway problem
!
!
        if (iGotRuna .ge. 1) then
          nGotRuna = nGotRuna + 1
          ipLastRuna = ip
        endif
 
!     .                                 compute d/dt {n_runa; J-runa}
!     .                                 as in eqn (20) and (21c) of
!     .                                 Karney Fisch paper.  Aug 93
        call ddtNrnJrn ( ip,nRunDot(ip), jRunDot(ip), vRunIdx(ip) )
 
 
 10   continue
!
        if (nGotRuna .ge. 1) then
          write(ErrMsg,                                                  &  
     &     '(1x,i4,'' Rway shells; last at ip='',i3)')                   &  
     &       nGotRuna, ipLastRuna
          call LSCwarn(ErrMsg)
        endif
!!           call JrfDiagn was once here...seemed to double the graphs
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"

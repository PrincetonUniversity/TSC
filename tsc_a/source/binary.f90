      subroutine binary(rl,zl,rr,zr,pval,rt,zt)
!
      USE CLINAM
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER il,jl,ir,jr,igt,jgt,ilt,jlt,l,itry,jtry
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 zl,rr,zr,pval,rt,zt,rl,zpl,pinterp,zpr,zpgt,rgt,zgt
      REAL*8 zplt,rlt,zlt,rtry,ztry,zptry
!============
      il = (rl-ccon)/deex + 2.0_R8
      jl = (zl-zzero)/deez + 2.0_R8
      zpl = pinterp(rl,zl,il,jl)
      ir = (rr-ccon)/deex + 2.0_R8
      jr = (zr-zzero)/deez + 2.0_R8
      zpr = pinterp(rr,zr,ir,jr)
      if(zpl.le.pval .and. zpr.le.pval) go to 400
      if(zpl.ge.pval .and. zpr.ge.pval) go to 400
      if(zpl.gt.pval) go to 100
      zpgt = zpr
      rgt = rr
      zgt = zr
      igt = ir
      jgt = jr
      zplt = zpl
      rlt = rl
      zlt = zl
      ilt = il
      jlt = jl
      go to 200
  100 continue
      zpgt = zpl
      rgt = rl
      zgt = zl
      igt = il
      jgt = jl
      zplt = zpr
      rlt = rr
      zlt = zr
      ilt = ir
      jlt = jr
  200 continue
      do 300 l=1,20
      rtry = 0.5_R8*(rgt+rlt)
      ztry = 0.5_R8*(zgt+zlt)
      itry = (rtry-ccon)/deex + 2.0_R8
      jtry = (ztry-zzero)/deez + 2.0_R8
      zptry = pinterp(rtry,ztry,itry,jtry)
      if(zptry.gt.pval) go to 310
      rlt = rtry
      zlt = ztry
      go to 300
  310 continue
      rgt = rtry
      zgt = ztry
  300 continue
      rt = 0.5_R8*(rgt+rlt)
      zt = 0.5_R8*(zgt+zlt)
      return
  400 continue
      if(abs(pval-zpl) .lt. abs(pval-zpr)) go to 450
      rt = rr
      zt = zr
      return
  450 continue
      rt = rl
      zt = zl
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

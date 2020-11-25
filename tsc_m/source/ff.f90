      function ff(itype,aa,i,j)
!
!****
      USE CLINAM
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!****
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER itype,i,j
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 aa,ff,phiarg,phix,eb,ebx,ppf
!============
      phiarg=(psilim-aa)*delpsi
      ff=0.0_R8
      if(iexv(i,j).eq.1 .or. iexs(i,j).eq.1) return
      if(phiarg.le.0.0_R8)return
!
      go to(1,2,3,4,5,6,7,8),itype
    1 continue
      ff=phiarg**(alphag+1._R8)/(delpsi*(alphag+1._R8))
      go to 10
    2 continue
      ff=phiarg**(alphag+1)*(4._R8/(alphag+1)-4._R8*phiarg/(alphag+2))/  &  
     & delpsi
      go to 10
    3 continue
      ff=-phiarg**alphag
      go to 10
    4 continue
      ff=-4.0_R8*phiarg**alphag*(1._R8-phiarg)
      go to 10
    5 continue
      phix = 1._R8-phiarg
      eb = exp(-alphag)
      ebx = exp(-alphag*phix)
      if(alphag.eq.0) go to 40
      ppf = p0*(ebx-eb)/(eb-1)*delpsi
      go to 50
   40 ppf = -p0*phiarg*delpsi
   50 continue
      ff = (1._R8/betaj-1._R8)*xplas**2*ppf
      go to 10
    6 continue
      ff = delpsi*(delg-1._R8)*alphag*phiarg**(alphag-1._R8)             &  
     &   * (1._R8+(delg-1._R8)*phiarg**alphag)
      go to 10
 7    continue
      ff = alphag*gzero*delpsi*phiarg**(alphag-1._R8)
      go to 10
 8    continue
      ff = alphag*delpsi*phiarg**(2._R8*alphag-1._R8)
   10 continue
!
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

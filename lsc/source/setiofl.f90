 
      SUBROUTINE setiofl(i)
      USE params
      USE PlPr
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER PLflst(NPLFLG), PRflst(NPRFLG), i, j, mkout
      if(i .eq. 1)then
         do 10 j = 1, nPlFlg
            PLflst(j) = PlFlg(j)
 10      continue
         do 20 j = 1, nPrFlg
            PRflst(j) = PrFlg(j)
 20      continue
      endif
      do 30 j = 1, nPlFlg
         PlFlg(j) = FALSE
 30   continue
      do 40 j = 1, nPrFlg
         PrFlg(j) = FALSE
 40   continue
      mkout = FALSE
      do 50 j = 1, NDIAG
         if(i .eq. idiag(j))then
            mkout = TRUE
            go to 100
         endif
 50   continue
!
 100  continue
!
      if(mkout .eq. TRUE)then
         do 60 j = 1, nPlFlg
            PlFlg(j) = PLflst(j)
 60      continue
         do 70 j = 1, nPrFlg
            PrFlg(j) = PRflst(j)
 70      continue
      endif
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"

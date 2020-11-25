      subroutine pelsource(p1,p2,p3,p4,pa,pb,frac)
!
!.....compute overlap of pellet trajectory and zone
!============
! idecl:  explicitize implicit REAL declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      REAL*8 p2,p3,p4,pa,pb,frac,p1,a1,a2,a3,area,a
!============
      a1 = .5_R8*(p2-p1)
      a2 =     p3-p2
      a3 = .5_R8*(p4-p3)
      area = a1+a2+a3
!.Case 1
      if(pb .le. p1) then
        frac = 0._R8
        return
      endif
!
      if(pa .le. p1) then
!.Case 2
        if(pb .le. p2) then
          a = .5_R8*(pb-p1)**2/(p2-p1)
          frac = a/area
          return
        endif
!.Case 3
        if(pb .le. p3) then
          a = a1 + pb-p2
          frac = a/area
          return
        endif
!.Case 4
        if(pb .le. p4) then
          a = area - .5_R8*(p4-pb)**2/(p4-p3)
          frac = a/area
          return
        endif
!.Case 5
        if(pb .gt. p4) then
          frac = 1._R8
          return
        endif
      endif
!
      if(pa .le. p2) then
!.Case 6
        if(pb .le. p2) then
          a = .5_R8*((pb-p1)**2-(pa-p1)**2)/(p2-p1)
          frac = a/area
          return
        endif
!.Case 7
        if(pb .le. p3) then
          a = a1  - .5_R8*(pa-p1)**2/(p2-p1)+ pb-p2
          frac = a/area
          return
        endif
!.Case 8
        if(pb .le. p4) then
          a = area - .5_R8*(p4-pb)**2/(p4-p3)                            &  
     &             - .5_R8*(pa-p1)**2/(p2-p1)
          frac = a/area
          return
        endif
!.Case 9
        if(pb .gt. p4) then
          a = area - .5_R8*(pa-p1)**2/(p2-p1)
          frac = a/area
          return
        endif
      endif
!
      if(pa .le. p3) then
!.Case 10
        if(pb .le. p3) then
          a = pb - pa
          frac = a/area
          return
        endif
!.Case 11
        if(pb .le. p4) then
          a = p3 - pa + a3 - 0.5_R8*(p4-pb)**2/(p4-p3)
          frac = a/area
          return
        endif
!.Case 12
        if(pb .gt. p4) then
          a = p3 - pa + a3
          frac = a/area
          return
        endif
      endif
      if(pa .le. p4) then
!.Case 13
        if(pb .le. p4) then
          a = 0.5_R8*((p4-pa)**2 - (p4-pb)**2)/(p4-p3)
          frac = a/area
          return
        endif
!.Case 14
        if(pb .gt. p4) then
          a = .5_R8*(p4-pa)**2/(p4-p3)
          frac = a/area
          return
        endif
      endif
!
!
!.Case 15
      frac = 0
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

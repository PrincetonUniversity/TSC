      subroutine jetsource
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      USE CLINAM
      USE SAPROP
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER j
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 smlt,dsumjet
!============
      save
!
!.....defines edge source terms for density equation for isurf=1 and idens=1
!
!
      if(times.lt.acoef(817) .or. times.gt.acoef(818)) go to 50
!
!....... jet on
      if(sjane0.le.0) then
      acoef(818) = times + (acoef(818)-acoef(817))
      acoef(817) = times
      call jeton
               endif
!
      if((times+dts).gt.acoef(818) ) then
      smlt = (acoef(818)-times)/dts
      do 20 j = 2,npsit
      sravejet(j) = smlt*sravejet(j)
   20 continue
                    endif
      dsumjet = 0._R8
      do 30 j = 2,npsit
      dsumjet = dsumjet + vp(j)*sravejet(j)*dts*udsd/udst
   30 continue
!
      dsumjet = dsumjet/(1._R8-acoef(819))
      if((sumjet+dsumjet).gt.acoef(816)*sjane0                           &  
     &     .and. dsumjet.gt.0._R8) then
      smlt = (acoef(816)*sjane0-sumjet)/dsumjet
      do 40 j = 2,npsit
      sravejet(j) = smlt*sravejet(j)
   40 continue
      dsumjet = smlt*dsumjet
                     endif
      sumjet  = sumjet + dsumjet
      return
!
!....... jet off
   50 do 60 j=2,npsi
      sravejet(j) = 0._R8
   60 continue
!
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

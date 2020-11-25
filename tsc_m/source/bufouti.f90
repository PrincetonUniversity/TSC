      subroutine bufouti(lu, ibeg, ifin)
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER ibeg,ifin,lu,i
!============
      dimension ibeg(*)
!============      
      integer loc,len,lnw
!
      lnw = loc(ibeg(2)) - loc(ibeg(1))
      len = (loc(ifin) - loc(ibeg))/lnw + 1
!ccccc      write(*,*) ' len=', len
!     write(96,1001) loc(ibeg),loc(ifin),len
!1001 format(" common block written(i)",3i12)
      write(lu) (ibeg(i),i=1,len)
!
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

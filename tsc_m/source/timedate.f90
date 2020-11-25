      subroutine timedate(itime,idate,imach,isuffix)
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      character*10 idate,itime,izone,imach,isuffix
      integer date_time(8)
      call date_and_time(idate,itime,izone,date_time)
      imach = 'a         '
      isuffix='a         '
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

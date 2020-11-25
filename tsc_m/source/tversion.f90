!./tsc_m/source/tversion.f90
!
        subroutine  tversion (lu)
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER lu
!============
      character*10 idate,itime,imach,isuffix
      character*100 svnVersionString, get_svn_version
!
        call timedate(itime,idate,imach,isuffix)
        svnVersionString = get_svn_version()
!!!!!!!!        write (lu ,10)   itime,idate,imach,isuffix
        write (lu ,10)   itime,idate
   10   format(" TSC Version UNX10.8(v152)  : 20 Dec 10", 7x, 2(a10,2x),2a10/)
        return
        end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

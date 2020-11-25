      subroutine movieerr(iplacearg,iflagarg)
      USE ezcdf
      USE CLINAM
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER iflagarg,iplacearg
!============
!     include 'netcdf.inc'
      write(nout,1001) iplacearg,iflagarg
      write(nterm,1001)iplacearg,iflagarg
      print *, nf_strerror(iflagarg)
 1001 format(" error in moviedata,  iplace=",i3,"   iflag=",i3)
      ineg=58
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

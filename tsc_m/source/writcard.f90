      subroutine writcard(nunit,itype,card,char)
!....... 5.32    writcard
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER itype,nunit,i
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 card
!============
      character*80 char
      dimension card(9)
!============      
!
!...........................................................
!
!.....writes input card data to unit = "nunit"
!........... itype .lt. 0 writes comment card "char"
!........... itype .ge. 0 writes data    card "card"
!
!..........................................................
!
      if(itype .ge. 0) go to 100
      write(nunit,1000) char(1:80)
 1000 format(a80)
      return
!
  100 continue
      write(nunit,2000) itype,(card(i),i=1,7)
 2000 format(1x,i2,1p7e11.3)
!
!
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

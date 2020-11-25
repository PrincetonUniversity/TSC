        subroutine tgtolow (string,outstr,num)
 
!*************************************************************************
!
!  tgtolow  -  convert a string to lower case
!
!  synopsis     call tgtolow (string)
!               character*(*) string    String to convert
!
!*************************************************************************
 
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER num,i
!============
        character*(*) string
        character*(*) outstr
 
        do 10 i=1,num
          if (ichar(string(i:i)) .ge. ichar('A') .and.                   &  
     &        ichar(string(i:i)) .le. ichar('Z')) then
            outstr(i:i) = char (ichar(string(i:i))+32)
          else
            outstr(i:i) = string(i:i)
          endif
10      continue
 
        return
        end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

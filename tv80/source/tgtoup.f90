        subroutine tgtoup (string,outstr,num)
 
 
!*************************************************************************
!
!  tgutil.F  -  Tv80 to gks shell utility routines
!
!  contents     tgtoup          Convert string to uppercase
!               tgtolow         Convert string to lowercase
!               clnxy           Convert linear to normalized
!               cufxy           Convert user to transformed normalized
!               culxy           Convert user to linear
!        tgerror     Print an error message
!        tgisfont    Test for validity of a font
!
!    National Energy Research Supercomputer Center
!    Lawrence Livermore National Laboratory
!    University of California
!    Livermore, California 94550, USA
!
!    (C) Copyright 1991 The Regents of the University of California.
!    All Rights Reserved.
!
!*************************************************************************
 
!*************************************************************************
!
!  tgtoup  -  convert a string to upper case
!
!  synopsis     call tgtoup (string)
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
          if (ichar(string(i:i)) .ge. ichar('a') .and.                   &  
     &        ichar(string(i:i)) .le. ichar('z')) then
            outstr(i:i) = char (ichar(string(i:i))-32)
          else
            outstr(i:i) = string(i:i)
          endif
10      continue
 
        return
        end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

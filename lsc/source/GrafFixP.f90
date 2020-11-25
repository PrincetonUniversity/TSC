 
      SUBROUTINE GrafFixP ( MyString, Iargumnt )
      USE tscunits
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!     SUBROUTINE gRafFltP ( MyString, Rargumnt )  ! given as ENTRY !
!     SUBROUTINE GrafChrS ( MyString, Cargumnt )  ! given as ENTRY !
!============
! idecl:  explicitize implicit REAL declarations:
!============
      CHARACTER*8 LSCdate
      COMMON /LSCversn/ LSCdate
      CHARACTER*8 Cargumnt
      INTEGER Iargumnt
      REAL*8    Rargumnt
      INTEGER istring
      INTEGER i,j, idelta, jdelta
      CHARACTER*10 MyString
      CHARACTER*21 PlString
      CHARACTER*40 TiString
      DATA idelta, jdelta , istring / 340, 25, 0 /
      REAL*8    ZERO, ONE
      DATA    ZERO, ONE/                                                 &  
     &         0.0_R8, 1.0_R8/
      write(PlString,'(a10,  i10,''$'')') MyString, Iargumnt
      go to 10
      ENTRY      GrafFltP ( MyString, Rargumnt )
      write(PlString,'(a10,f10.4,''$'')') MyString, Rargumnt
      go to 10
      ENTRY      GrafChrS ( MyString, Cargumnt )
      write(PlString,'(a10,2x,a8,''$'')') MyString, Cargumnt
10    continue
!
      if (MyString .eq. 'BEGIN PLOT' .or.                                &  
     &    istring  .eq. 0                 ) then
	call EZinit
        call GNinit(nTSCgraf)
        call EZsets(1, 1023, 1, 768,ZERO, ONE,ZERO, ONE, 1)
        call GNsets(1, 1023, 1, 768,ZERO, ONE,ZERO, ONE, 1)
        write(TiString,                                                  &  
     &  '(''LSC_'',a8, '' inputs and parameters:$'')') LSCdate
        call EZwrit(100, 735,TiString,0,0)
        call GNlist(100, 735,TiString,0,0)
	j = 700
	i =  20
        istring = 1
        return
      endif
!
      if (MyString .eq. 'END PLOT  ') then
        istring = 0
        call EZfini(0,0)
        call GNfinL
        return
      endif
!
      call EZwrit(i, j, PlString, 0,0)
      call GNlist(i, j, PlString, 0,0)
      istring = istring + 1
      j = j - jdelta
      if (j .le. jdelta ) then
        j = 700
        i = i + idelta
      endif
!
      if (i .ge. 3*idelta) then
        call EZfini(0,0)
        call GNfinL
        istring = 0
      endif
!
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"

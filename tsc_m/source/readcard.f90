      subroutine readcard(nin,nout,nscr,itype,card,char,inegp)
!...........................................................
!
!.....this subroutine checks for illegal data cards.
!
!...........................................................
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER nout,nscr,itype,inegp,nin,j,isum,i,i1,i2,i3,i10,i9,i8
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 card
!============
      character*8 char*80,datac,a1,c,dec,e,blk,string*80
      character*2 blk2
      dimension  card(9)
      dimension  a1(80),datac(80)
!
      c = 'c'
      dec = '.'
      e = 'e'
      blk = ' '
      blk2 = '  '
      itype = -1
!
  102 format (a80)
  108 format (i2,8x,7g10.2)
  120 format (2x,'there is an illegal decimal point on this data card.')    
  121 format(2x,"  number displaced to the left on this data card ")
  122 format(2x,"  decimal point missing on this data card ")
!
  200 continue
!
      read(nin,102) char(1:80)
!
!.....if 'c' is in column one read another card.
!
      if (char(1:1).eq.c) return
      if (char(1:1).eq.'C') return
      if (char(1:2) .eq. blk2) return
!
!.....check for illegal decimal point.
!
      do 235 j=1,7
      isum=0
      do 233 i=1,10
      i1=10*j+i
      if (char(i1:i1).eq.dec) isum=isum+1
  233 continue
      if (isum.gt.1) go to 400
  235 continue
!
!.....check and adjust the righthand alignment of e format numbers.
!
      do 239 j=1,7
      do 237 i=1,10
      i1=10*j+i
      if (char(i1:i1).eq.e) go to 238
  237 continue
      go to 239
  238 continue
      i2=10*j+10
      if (char(i2:i2).ne.blk) go to 239
      do 241 i=1,9
      i3=10*(j+1)+1-i
      char(i3:i3)=char(i3-1:i3-1)
  241 continue
      i3=10*j+1
      char(i3:i3)=blk
      go to 238
  239 continue
!
!.....check for number shifted one space into wrong field
!
      do 243 j=1,7
      i10 = 10*j+10
      i9 = i10 - 1
      i8 = i9  - 1
      if(char(i10:i10) .ne. blk .and. char(i9:i9) .eq. blk) go to 450
      if(char(i8:i8).eq.blk .and. (char(i9:i9).ne.blk .and.              &  
     &  char(i10:i10).ne.blk)                                            &  
     & .and. (char(i8:i8).ne.e .and. char(i9:i9).ne.e .and.              &  
     &  char(i10:i10).ne.e)                                              &  
     &  .and. (char(i9:i9).ne.dec .and. char(i10:i10).ne.dec)) go to     &  
     & 450
  243 continue
!
!.....check for at least 1 decimal point
      do 255 j=1,7
      do 253 i=1,10
      i1 = 10*j + i
      if(char(i1:i1).ne.blk) go to 254
  253 continue
      go to 255
  254 continue
      do 256 i=1,10
      i1 = 10*j + i
      if(char(i1:i1).eq.dec) go to 255
  256 continue
      write(nout,102) char(1:80)
      write(nout,122)
      inegp = 17
      return
  255 continue
      go to 500
!
  400 continue
      write (nout,102) char(1:80)
      write (nout,120)
      inegp = 17
      return
!
  450 continue
      write(nout,102) char(1:80)
      write(nout,121)
      inegp = 17
      return
!
  500 continue
      write(string,102) char(1:80)
      read(string,108)  itype,(card(i),i=1,7)
!
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

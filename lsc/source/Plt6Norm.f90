!     ------------------------------------------------------------------
      SUBROUTINE Plt6Norm (x, y, NumPts, chrts6, iplace, GrTyp)
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      CHARACTER*1 GrTyp
      INTEGER iplace, NumPts, MAXPts, n
      PARAMETER (MAXPts = 451)
      REAL*8     x(NumPts),  y(NumPts)
      REAL*8    xg(MAXPts), yg(MAXPts)
      CHARACTER *(*) chrts6
      INTEGER  PltsOnPg, i, iCross0, iydiv(2), iyinv(2)
      PARAMETER (PltsOnPg =6)
      INTEGER nxlo(PltsOnPg) , nxhi(PltsOnPg),                           &  
     &        nylo(PltsOnPg) , nyhi(PltsOnPg), nywri(PltsOnPg)
      INTEGER ixlo(PltsOnPg) , ixhi(PltsOnPg),                           &  
     &        iylo(PltsOnPg) , iyhi(PltsOnPg), iywri(PltsOnPg)
 
 
      REAL*8    xmin, xmax, ymin(2), ymax(2)
      REAL*8    ymnVal, ymxVal
      CHARACTER*24 MyString
      CHARACTER*24 NDString
      DATA     xmin    , xmax /                                          &  
     &           0.0_R8,  1.0_R8/
      DATA     ymin(1) , ymin(2) , ymax(1) , ymax (2) /                  &  
     &           0.0_R8,  -2.0_R8,  1.5_R8,   2.0_R8/
      DATA     iydiv(1), iydiv(2), iyinv(1), iyinv(2) /                  &  
     &            3    ,   2     ,   1     ,   2      /
      DATA    (nxlo(i) , i=1,PltsOnPg)                                   &  
     &       / 50,      400,      750,       50,      400,      750  /
      DATA    (nxhi(i), i=1,PltsOnPg)                                    &  
     &    /   300,      650,     1000,      300,      650,     1000  /
      DATA    (nylo(i) , i=1,PltsOnPg)                                   &  
     &      / 400,      400,      400,        50,      50,       50  /
      DATA   (nyhi(i) , i=1,PltsOnPg)                                    &  
     &      / 650,      650,      650,       300,     300,      300  /
      DATA   (nywri(i) , i=1,PltsOnPg)                                   &  
     &      / 625,      625,      625,       275,     275,      275  /
 
      DATA    (ixlo(i) , i=1,PltsOnPg)                                   &  
     &       /  0,      340,      680,        0,      340,      680  /
      DATA    (ixhi(i), i=1,PltsOnPg)                                    &  
     &    /   340,      680,     1020,      340,      680,     1020  /
      DATA    (iylo(i) , i=1,PltsOnPg)                                   &  
     &      / 380,      380,      380,        0,        0,        0  /
      DATA   (iyhi(i) , i=1,PltsOnPg)                                    &  
     &      / 760,      760,      760,       380,     380,      380  /
      DATA   (iywri(i) , i=1,PltsOnPg)                                   &  
     &      / 625,      625,      625,       275,     275,      275  /
 
!
!     -                                 Set x to root of normalized \psi
 
      n = NumPts
      if (n .gt. MAXPts) n = MAXPts
      do 10 i = 1, n
        xg(i) = ( x(i) - x(1) ) /                                        &  
     &     ( x(NumPts) - x(1) )
        xg(i) = sqrt ( xg(i) )
 10   continue
!
      i = iplace
      call EZnorm ( y, yg , ymnVal, ymxVal, n, iCross0 )
      write (MyString, '( a6, '' /'', 1pe8.1,''$'')') chrts6, ymxVal
      write (NDString, '( a6, '' /'', 1pe8.1,''$'')') chrts6, ymxVal
      call EZsets ( nxlo(i), nxhi(i), nylo(i), nyhi(i),                  &  
     &                xmin, xmax, ymin(iCross0),ymax(iCross0), 1 )
      call EZaxes ( 1 , 10, iydiv(iCross0), iyinv(iCross0) )
      call EZwrit ( nxlo(i), nywri(i), MyString , 0, 0 )
 
      call GNsets ( ixlo(i), ixhi(i), iylo(i), iyhi(i),                  &  
     &                xmin, xmax, ymin(iCross0),ymax(iCross0), 1 )
      call EZwrit ( nxlo(i), nywri(i), MyString , 0, 0 )
      call GNtitl (                    NDString        )
 
      if                                                                 &  
     & (GrTyp .eq. ' ' .or. GrTyp .eq. 'l' .or. GrTyp .eq. '-') then
        call EZcurv ( xg , yg , n )
        call GNcurv ( xg , yg , n )
      else if                                                            &  
     & (GrTyp .eq. 'x' .or. GrTyp .eq. 'X' .or. GrTyp .eq. '*') then
        call EZcros ( xg , yg , n )
        call GNcros ( xg , yg , n )
      else if                                                            &  
     & (GrTyp .eq. 'b' .or. GrTyp .eq. 'B' .or. GrTyp .eq. '|') then
        call EZbars ( xg , yg , n ,'y')
        call GNbars ( xg , yg , n ,'y')
      else if                                                            &  
     & (GrTyp .eq. 'p' .or. GrTyp .eq. 'P' .or. GrTyp .eq. '.') then
        call EZpnts ( xg , yg , n )
        call GNpnts ( xg , yg , n )
      else if                                                            &  
     & (GrTyp .eq. 's' .or. GrTyp .eq. 'S' .or. GrTyp .eq. ',') then
        call EZpnt2 ( xg , yg , n )
        call GNpnt2 ( xg , yg , n )
      else
        call EZcurv ( xg , yg , n )
        call GNcurv ( xg , yg , n )
      endif
!
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"

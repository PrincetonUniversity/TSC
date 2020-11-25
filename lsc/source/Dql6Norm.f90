!
!     -----------------------------------------------------------------|
!
      SUBROUTINE Dql6Norm (x, y, NumPts, chrst10, GrTyp,                 &  
     &                                               ThisPlot, TotlPlts)
      USE tscunits
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      CHARACTER*1 GrTyp
      INTEGER MAXPts, n
      INTEGER NumPts, ThisPlot, TotlPlts, PltsOnPg, PltIndx
      INTEGER FirstInFrame
      PARAMETER (PltsOnPg =6)
      PARAMETER (MAXPts = 451)
      REAL*8     x(NumPts),  y(NumPts)
      REAL*8                yg(MAXPts)
      CHARACTER*10 chrst10
      INTEGER i, iCross0, iydiv(2), iyinv(2), iyldiv, iylinv
      INTEGER xlow(PltsOnPg) , xhigh(PltsOnPg),                          &  
     &        ylow(PltsOnPg) , yhigh(PltsOnPg), ywrit(PltsOnPg)
 
      INTEGER ixlo(PltsOnPg) , ixhig(PltsOnPg),                          &  
     &        iylo(PltsOnPg) , iyhig(PltsOnPg)
 
      INTEGER ytitl
      REAL*8    xmin, xmax, ymin(2), ymax(2), ylmn, ylmx
      REAL*8    ymnVal, ymxVal
      CHARACTER*1 logtest
      CHARACTER*24 MyString
      CHARACTER*30 TSCstring
      DATA     xmin    , xmax , ymnVal, ymxVal /                         &  
     &          -1.0_R8,  1.0_R8,    0.0_R8,    1.0_R8/
      DATA     ymin(1) , ymin(2) , ymax(1) , ymax (2) , ylmn , ylmx /    &  
     &           0.0_R8,  -2.0_R8,  1.5_R8,   2.0_R8,  -5.0_R8,  1.0_R8/
      DATA     iydiv(1), iydiv(2), iyinv(1), iyinv(2), iyldiv, iylinv /  &  
     &            3    ,   2     ,   1     ,   2     ,    6  ,    1   /
      DATA    (xlow(i) , i=1,PltsOnPg)                                   &  
     &       / 50,      400,      750,       50,      400,      750  /
      DATA    (xhigh(i), i=1,PltsOnPg)                                   &  
     &    /   300,      650,     1000,      300,      650,     1000  /
      DATA    (ylow(i) , i=1,PltsOnPg)                                   &  
     &      / 400,      400,      400,        50,      50,       50  /
      DATA   (yhigh(i) , i=1,PltsOnPg)                                   &  
     &      / 650,      650,      650,       300,     300,      300  /
      DATA   (ywrit(i) , i=1,PltsOnPg)                                   &  
     &      / 660,      660,      660,       310,     310,      310  /
      DATA    ytitl   / 725 /
 
      DATA    (ixlo(i) , i=1,PltsOnPg)                                   &  
     &       /  0,      340,      680,        0,      340,      680  /
      DATA    (ixhig(i), i=1,PltsOnPg)                                   &  
     &    /   340,      680,     1020,      340,      680,     1020  /
      DATA    (iylo(i) , i=1,PltsOnPg)                                   &  
     &      / 380,      380,      380,         0,       0,        0  /
      DATA   (iyhig(i) , i=1,PltsOnPg)                                   &  
     &      / 760,      760,      760,       380,     380,      380  /
 
!
      PltIndx = mod(ThisPlot-1, PltsOnPg) + 1
      if (GrTyp .eq. ' ' .or. GrTyp .eq. '-' .or.                        &  
     &    GrTyp .eq. 'l' .or. GrTyp .eq. 'x' .or. GrTyp .eq. '|') then
          FirstInFrame = 1
        else
          FirstInFrame = 0
      endif
 
      if (ThisPlot .eq. 1 .and. FirstInFrame .eq. 1 ) then
         call EZinit
         call GNinit(nTSCgraf)
      endif
      n = NumPts
      if (n .gt. MAXPts) n = MAXPts
      logtest = chrst10(1:1)
!
 
      if (FirstInFrame .eq. 1) then
        call EZnorm ( y, yg , ymnVal, ymxVal, n, iCross0 )
        write (MyString, '( a10,''/'', 1pe8.1,''$'')') chrst10, ymxVal
        call EZwrit(xlow(PltIndx), ywrit(PltIndx), MyString , 0, 0 )
        call GNtitl(                               MyString        )
        if (logtest .ne. 'l' .and. logtest .ne. 'L' ) then
          call EZsets ( xlow(PltIndx), xhigh(PltIndx),                   &  
     &                  ylow(PltIndx), yhigh(PltIndx),                   &  
     &                  xmin,           xmax,                            &  
     &                  ymin(iCross0),   ymax(iCross0),    1 )
 
          call GNsets ( ixlo(PltIndx), ixhig(PltIndx),                   &  
     &                  iylo(PltIndx), iyhig(PltIndx),                   &  
     &                  xmin,           xmax,                            &  
     &                  ymin(iCross0),   ymax(iCross0),    1 )
 
          if(GrTyp .eq. 'x') then
            call EZcros ( x  , yg , n )
            call GNcros ( x  , yg , n )
          else if(GrTyp .eq. '|') then
            call EZbars ( x  , yg , n ,'y')
            call GNbars ( x  , yg , n ,'y')
          else
            call EZcurv ( x  , yg , n )
            call GNcurD ( x  , yg , n )
          endif
          call EZaxes ( 1 , 10, iydiv(iCross0), iyinv(iCross0) )
        else
          do 10 i = 1, n
            yg(i) = log10( abs(yg(i)) + 1.E-6_R8)
            if (yg(i) .lt. ylmn ) yg(i) = ylmn
            if (yg(i) .gt. ylmx ) yg(i) = ylmx
 10       continue
          call EZsets ( xlow(PltIndx), xhigh(PltIndx),                   &  
     &                  ylow(PltIndx), yhigh(PltIndx),                   &  
     &                  xmin,           xmax,                            &  
     &                  ylmn          , ylmx             , 1 )
 
          call GNsets ( ixlo(PltIndx), ixhig(PltIndx),                   &  
     &                  iylo(PltIndx), iyhig(PltIndx),                   &  
     &                  xmin,           xmax,                            &  
     &                  ylmn          , ylmx             , 1 )
 
          if(GrTyp .eq. 'x') then
            call EZcros ( x  , yg , n )
            call GNcros ( x  , yg , n )
          else if(GrTyp .eq. '|') then
            call EZbars ( x  , yg , n ,'y')
            call GNbars ( x  , yg , n ,'y')
          else
            call EZcurv ( x  , yg , n )
            call GNcurD ( x  , yg , n )
          endif
          call EZaxes ( 1 , 10, iyldiv, iylinv )
        endif
 
      else
        if (logtest .ne. 'l' .and. logtest .ne. 'L' ) then
          do 20 i = 1, n
            yg(i) = y(i)/ymxVal
            if (yg(i) .gt. 1.48_R8) yg(i) = 1.48_R8
            if (yg(i) .lt. 0.02_R8) yg(i) = 0.02_R8
 20       continue
          call EZpnts ( x  , yg , n )
          call GNpntD ( x  , yg , n )
        else
          do 30 i = 1, n
            yg(i) = log10( abs(y(i)/ymxVal) + 1.E-6_R8)
            if (yg(i) .lt. ylmn ) yg(i) = ylmn
            if (yg(i) .gt. ylmx ) yg(i) = ylmx
 30       continue
          call EZpnts ( x  , yg , n )
          call GNpntD ( x  , yg , n )
        endif
!
      endif
!
!
      if (( PltIndx .eq. PltsOnPg .or. ThisPlot .eq. TotlPlts ))then
      if ( FirstInFrame .eq. 0 .or. GrTyp .eq. 'x'                       &  
     &                         .or. GrTyp .eq. '|' ) then
        write(TSCstring,'( a10, '' vs v/c: page of 6  '')') chrst10
        call EZwrit(xlow(1),ytitl,                                       &  
     &                   TSCstring // '$' ,0,0)
        call MkGrfLst( TSCstring )
        call EZfini(0,0)
        call GNfini
      endif
      endif
!
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"

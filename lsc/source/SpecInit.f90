!     SpecInit begins -------------------------------------------------|
!     -                                                                |
!     -                                                                |
      SUBROUTINE SpecInit
      USE Doflags
      USE params
      USE PlPr
      USE RayBins
      USE RayWrk
      USE tscunits
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!============
! idecl:  explicitize implicit REAL declarations:
!============
      INTEGER oldNtors, oldNpols, oldnGrps, oldNslic, oldDoBram
      REAL*8                                                             &  
     &        oldPowrs(NGRPDIM),oldCents(NGRPDIM),oldWidts(NGRPDIM)
      REAL*8                                                             &  
     &        oldPhsDg(NGRPDIM)
      CHARACTER*8 oldCoupl(NGRPDIM)
      CHARACTER*80 pline
      CHARACTER*40 MyString
      INTEGER igrup, i, iry, ity, ipy, ifirst
      INTEGER ixmjr, ixmnr, iymjr, iymnr, nRyPrGr
      REAL*8                                                             &  
     &        powtot, RayGroup, GrpNorm(NGRPDIM)
      REAL*8                                                             &  
     &        eNs (NTORDIM,NGRPDIM), SofNs (NTORDIM,NGRPDIM)
      REAL*8                                                             &  
     &        ens1(NGRPDIM*NTORDIM), SofNs1(NGRPDIM*NTORDIM)
      EQUIVALENCE (  eNs (1,1),   eNs1(1) )
      EQUIVALENCE (SofNs (1,1), SofNs1(1) )
      REAL*8                                                             &  
     &        SpShape, MaxOfAry, MaxPwr, smx, xmin, xmax, EXCLUDED
      DATA    EXCLUDED, xmin, xmax, ixmjr,ixmnr,iymjr,iymnr /            &  
     &           1.20_R8,-10.0_R8, 10.0_R8,     2,    5,    2,    1 /
      DATA    ifirst / 1 /
      REAL*8    ZERO, ONE
      DATA    ZERO, ONE/                                                 &  
     &         0.0_R8, 1.0_R8/
      EXTERNAL SpShape
!     .                                 If this is not the first call, then
!     .                                 test to see if any of the parameters
!     .                                 of the spectrum have changed.
!     .                                 If stable, then return to avoid another
!     .                                 calculation, which in brambJES is
!     .                                 time consuming.
      if (ifirst .eq. 0 ) then
        i = 0
        do 5 igrup=1,nGrps
          if( oldPowrs(igrup) .ne. powers  (igrup)) i=i+1
          if( oldCents(igrup) .ne. centers (igrup)) i=i+1
          if( oldWidts(igrup) .ne. widths  (igrup)) i=i+1
          if( oldPhsDg(igrup) .ne. phaseDeg(igrup)) i=i+1
          if( oldCoupl(igrup) .ne. couplers(igrup)) i=i+1
 5      continue
          if( oldNtors .ne. ntors   .or. oldnGrps .ne. nGrps             &  
     &   .or. oldNpols .ne. npols                                        &  
     &   .or. oldNslic .ne. nslices .or. oldDoBram .ne. DoBram )i=i+1
!
        if (i .eq. 0) return
!
        iRayTrsi = 1
!
        else
        ifirst = 0
      endif
!
!
!     .                                 Make the spectrum in poloidal number
      npol(1) = (npolmax+npolmin)/2.0_R8
      if (npols .gt. 1 )  call ugrid(npol,npols,npolmin,npolmax)
 
 
 
!     .                                 IF ----------------------------|
!     .                                 If we are not to do a Brambilla calc,
!     .                                 make a model spectrum using
!a    .                                 gaussian form.
                                        if(DoBram .eq.  0 ) then
      call ugridEXC (ntor, ntors, nparmin, nparmax, EXCLUDED)
      powtot = 0._R8
 
      do 17 igrup = 1, nGrps
        GrpNorm(igrup) = 0._R8
          do 15  i = 1, ntors
            GrpNorm(igrup) = GrpNorm(igrup) +                            &  
     &        SpShape(ntor(i), centers(igrup), widths(igrup))
 15       continue
        GrpNorm(igrup) = 1._R8/ GrpNorm(igrup)
 17   continue
!
      do 30 ity = 1, ntors
         Spec(ity) = 0._R8
         do 25 igrup = 1, nGrps
            RayGroup = GrpNorm(igrup) * powers(igrup) *                  &  
     &        SpShape(ntor(ity), centers(igrup), widths(igrup))
            Spec(ity) = Spec(ity) + RayGroup
 25      continue
 30   continue
!
      call vecnorm(ntors,Spec)
 
!
!     .                                 But if DoBram =1, then
!     .                                 get a Brambillla spectrum.
!     .                                 ELSE --------------------------|
                                        else if (DoBram .eq. 1) then
!
 
!     .                                 Parcel out a ray count by groups.
!     .                                 Over distribute if a it doesnt come out
!     .                                 even
                  nRyPrGr = ntors / nGrps
                  if (nGrps * nRyPrGr .lt. ntors) nRyPrgr = nRyPrGr + 1
!     .
      if (nRyPrGr .gt. NTORDIM) nRyPrGr = NTORDIM
      call BrodCast(NTORDIM*NGRPDIM, SofNs , 0._R8)
      call BrodCast(NTORDIM*NGRPDIM,   eNs , 0._R8)
      do 50 i=1,nGrps
        call BrambJES(fghz, phaseDeg(i), nRyPrGr, eNs(1,i), SofNs(1,i),  &  
     &                      EXCLUDED, nslices, couplers(i), TurnNegs )
 50   continue
      do 55 igrup=1,nGrps
        do 55 ity=1,nRyPrGr
          SofNs(ity,igrup) = powers(igrup)*SofNs(ity,igrup)
 55   continue
      call PikSr2NR(NTORDIM*NGRPDIM, SofNs, eNs)
      do 60 i = 1, ntors
        ntor(i) = eNs1  (i + NTORDIM*NGRPDIM-ntors )
        Spec(i)     = SofNs1(i + NTORDIM*NGRPDIM-ntors )
 60   continue
      call PikSr2NR(ntors, ntor, Spec)
      call vecnorm(ntors,Spec)
                                        endif
!     .                                 ENDIF -------------------------|
!     .
!     .                                 It can happen that there are two
!     .                                 rays at the same n_//. Check and report
      do 80 i = 2, ntors
        if ( ntor(i) .ne. ntor(i-1) ) go to 80
        write(MyString,'(''identical n// at '',1pe10.3)') ntor(i)
        call LSCwarn(MyString)
 80   continue
!
!
!     Option to write out the n-tors and Powers of rays     -----------|
!     -----------------------------------------------------------------|
!     Spectrum IO
      if ( PrFlg(NPAPWRWR) .ge. TRUE ) then
         write(nLSCcom2, 600)
 600     format(/,'  Ray Index     n-tor  Power(W)')
 601     format(1x,     8x,i2,1x,1pe9.2,1x,1pe9.2  )
         do 160 ity=1,ntors
            write(nLSCcomm, 601) ity, ntor(ity), Spec(ity)
 160     continue
      endif
!     Option ends -------------------------------------------------|
!
!     Option to plot ----------------------------------------------|
       if(PlFlg(SPECPL) .ge. TRUE) then
        call EZinit
        call GNinit(nTSCgraf)
        MaxPwr = MaxOfAry(Spec,ntors)
        call EZrnd2  (MaxPwr, smx, iymjr, iymnr)
        call ezsets(100,500,125,375,xmin,xmax,ZERO,smx,1)
!       call GNsets(100,500,125,375,xmin,xmax,ZERO,smx,1)
        call GNsets(  0,500,  0,375,xmin,xmax,ZERO,smx,1)
        call ezwrit(100,400,'Relative power vs n_tor$',0,0)
        call GNtitl(        'Relative power vs n_tor'     )
        call ezbars(ntor, Spec, ntors, 'y')
        call GNbars(ntor, Spec, ntors, 'y')
        call ezaxes( ixmjr, ixmnr,  iymjr, iymnr )
        if (npols .eq. 1) then
          write(pline,'(f5.2,'' GHz w/ '',                               &  
     &     i3, '' rays in '', i2,'' groups.$'' )') fghz, ntors, nGrps
        else
          write(pline,'(f5.2,'' GHz;'',                                  &  
     &     i3, '' ntors;'', i2, '' nGrps;'',                             &  
     &     i2,'' npols$'' )') fghz, ntors, nGrps, npols
        endif
        call ezwrit(  1,750, pline, 0, 0)
        call GNlist(  1,750, pline, 0, 0)
!
        if(DoBram .eq. 1) then
          write(pline,'(''Brambilla with '', i4,                         &  
     &                 '' slices.$'' )') nslices
          call ezwrit(  1,725, pline, 0, 0)
          call GNlist(  1,725, pline, 0, 0)
          write(pline,'( a8,'' at'', f7.1,                               &  
     &                      '' deg &'', f4.1,'' pwr$'')')                &  
     &                    couplers(1),phaseDeg(1),powers(1)
          call ezwrit(  1,700, pline, 0, 0)
          call GNlist(  1,700, pline, 0, 0)
          if (nGrps .ge. 2) then
          write(pline,'( a8,'' at'', f7.1,                               &  
     &                      '' deg &'', f4.1,'' pwr$'')')                &  
     &                    couplers(2),phaseDeg(2),powers(2)
          call ezwrit(  1,675, pline, 0, 0)
          call GNlist(  1,675, pline, 0, 0)
          endif
          if (nGrps .ge. 3) then
          write(pline,'( a8,'' at'', f7.1,                               &  
     &                      '' deg &'', f4.1,'' pwr$'')')                &  
     &                    couplers(3),phaseDeg(3),powers(3)
          call ezwrit(  1,650, pline, 0, 0)
          call GNlist(  1,650, pline, 0, 0)
          endif
 
        else
          write(pline,'(                                                 &  
     &       ''centers: '', 3(f5.2))') (centers (i),i=1,nGrps)
          call DolAtEnd(pline)
          call ezwrit(  1,725, pline, 0, 0)
          call GNlist(  1,725, pline, 0, 0)
 
          write(pline,'(                                                 &  
     &       ''widths:  '', 3(f5.2))') (widths  (i),i=1,nGrps)
          call DolAtEnd(pline)
          call ezwrit(  1,700, pline, 0, 0)
          call GNlist(  1,700, pline, 0, 0)
 
          write(pline,'(                                                 &  
     &       ''powers:  '', 3(f5.2))') (powers  (i),i=1,nGrps)
          call DolAtEnd(pline)
          call ezwrit(  1,675, pline, 0, 0)
          call GNlist(  1,675, pline, 0, 0)
 
        endif
        call pSmoUnsm(ntor,Spec,ntors)
        call EZwrit(  1,30,EquTitle (1:50)//'$',0,0)
        call GNlist(  1,30,EquTitle (1:50)//'$',0,0)
        call EZwrit(  1, 5,EquTitle(51:96)//'$',0,0)
        call GNlist(  1, 5,EquTitle(51:96)//'$',0,0)
        call MkGrfLst(' Launched lwr hybrid spectrum ')
 
        call EZfini(0,0)
        call GNfini
 
      endif
! dwi save        open (nTSCunus,status='old',file = FileNa1, err=1300 )
!     Option to plot ends -----------------------------------------|
!
!     Prepare to return by saving old values
      do 90 igrup = 1,nGrps
      oldPowrs(igrup) = powers(igrup)
      oldCents(igrup) = centers(igrup)
      oldWidts(igrup) = widths(igrup)
      oldPhsDg(igrup) = phaseDeg(igrup)
      oldCoupl(igrup) = couplers(igrup)
 90   continue
!     ntors    = nInclIBW !
      nrays    = ntors*npols
      oldNtors = ntors
      oldNpols = npols
      oldnGrps = nGrps
      oldNslic = Nslices
      oldDoBram= DoBram
!
!
      return
!
!     -----------------------------------------------------------------|
      ENTRY PowrInit
      iry=1
      do 100 ity=1,ntors
      do 100 ipy=1,npols
         npar (1,iry) = ntor(ity)
         iry = iry+1
 100  continue
      call rspwr(1._R8)
!     generate launch spectrum at full power
      do 110 iry = 1, nrays
        call blockcpy(nzones - 1, power(2, iry), 1, power(1, iry), 1)
!       copy
 110   continue
!     copy results onto entire trajectory
!
!
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"

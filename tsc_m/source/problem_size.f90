      subroutine problem_size
!=======================================================================
!***********************************************************************
!                                                                      *
!.....read in all input                                                *
!                                                                      *
!***********************************************************************
!
      USE PARAM
      USE CLINAM 
      USE SCRATCH
      USE PROFCOM
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
 
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER l62z,l63z
      INTEGER i,l,l16z,l17z,l18z,l23z,l24z
      INTEGER l26z,l27z,l28z,l29z,l30z,l31z,l32z,l34z,l35z,l36z
      INTEGER l42z,l43z,l44z,l45z,l46z,l47z,l48z,l50z,l51z,l52z
      INTEGER l53z,l54z,l55z,l56z,l57z,l58z,l59z,l60z,l61z,l64z
      INTEGER l65z,l66z,l67z,l68z,l69z,l70z,l71z,l72z,l73z,l75z
      INTEGER itype,inegp,ll,iabs,nlimch,nlim0,ii,isave
      INTEGER jsave,kobs,n,ic0sv,ic01,ic0,ico,nco,icmax,indx,ic
      INTEGER igroup,inumfb,lll,nnoplot,itempi,itroub,inot,irenum
      INTEGER in,io,ictype,igr,ios5
      INTEGER nlim_max, igroupw_max,ngr
      REAL*8, EXTERNAL ::  AREAL
      INTEGER istat
!============
      REAL*8 card(9)
      character*80 char
!
!...........................................................
!
!.....initialize some variables used in inpt
!
!...........................................................
      do 1 i=1,9
    1 card(i) = 0._R8
      ncoil   = 0
      nwire   = 0
      nobs    = 0
      numfb   = 0
      nlim    = 0
      nopest = 0
      numicsuf = 0
      numsaw = 0
      ngrmax = 0
      nlim_max = 0
      igroupw_max = 0
      nsepmax = 2
!...........................................................
!
!.....read name card
!
!...........................................................
!     read(nin,2000) (name(i),i=1,8)
      read(nin,2000) 
!===========================================================
!
!.....read type 0 card
!
!===========================================================
      go to(801,802),iformat
  801 continue
      call readcard(nin,nout,nsc1,itype,card,char,inegp)
      if(inegp.gt.0) ineg=inegp
      if(itype.lt.0) go to 801
      if(itype.eq.0) go to 803
      write(nout,1803)
 1803 format(" no type 00 card")
      ineg=21
      go to 803
  802 continue
      read(nin,1003) itype,(card(i),i=1,9)
 1003 format(i2,9f8.0)
  803 continue
!...........................................................
!
!...card  0 : control
!
!...........................................................
!
      irst1  = int(card(1))
      irst2  = int(card(2))
!===========================================================
!
!.....read card types .ge. one
!
!===========================================================
      ll = 0
      ineg=0
    5 continue
      go to(806,804),iformat
  806 call readcard(nin,nout,nsc1,itype,card,char,inegp )
      if(inegp.gt.0) ineg=inegp
      go to 805
  804 read(nin,1003) itype,(card(i),i=1,9)
  805 continue
      if(itype .lt. 0) go to 5
      if(itype.eq.99) go to 500
      if(itype.ne.0 .and. itype.le.98) go to 7
      write(nout,2007) itype
 2007 format(" illegal card type",i3)
      ineg=21
      go to 999
    7 continue
      if(itype .gt. 75) go to 5
!
!.....process card according to card type
      go to (10,20,30,40,50,60,70,80,90,100,110,120,                     &  
     &   130,140,150,160,170,180,190,200,210,220,230,240,                &  
     &   250,260,270,280,290,300,310,320,330,340,350,360,                &  
     &   370,380,390,400,410,420,430,440,450,                            &  
     &   460,470,480,490,495,5101,5201,5301,5401,5501,5601,              &  
     &   5701,5801,5901,6001,6101,6201,6301,6401,6501,6601,              &  
     &   6701,6801,6901,7001,7101,7201,7301,7401,7501),itype
!...........................................................
!
!...card  1 : dimensions
!
!...........................................................
   10 continue
!
      nx      = int(card(1))
      nz      = int(card(2))
      jsym    = int(card(5))
      isym=iabs(jsym)

      if(irst1.ne.1) go to 5
      nzp = nz + 1
      if(isym.eq.0) nzp = 2*(int(AREAL(nz+1)/2._R8+.1_R8))
      nz = nzp-1
      go to 5
!...........................................................
!
!...card  2 : time step and switches
!
!...........................................................
   20 continue
      go to 5
!...........................................................
!
!...card  3 : numerical
!
!...........................................................
   30 continue
      go to 5
!...........................................................
!
!...card  4 : surface average
!
!...........................................................
   40 continue
!
      npsi    = int(card(2))
      go to 5
!...........................................................
!
!...card  5 : limiter points
!
!...........................................................
   50 continue
!
      nlim        = int(card(1))
      nlimch      = int(card(1)+.9999999_R8)
      if(nlim.eq.nlimch) go to 54
      write(nterm,1054)
      write(nout ,1054)
      ineg=21
      go to 999
 1054 format(" *** error with type 05 (limiter) card *** ",//,           &  
     &       " in version 8.00 and all future TSC versions the first",/,  &  
     &                                                                   &  
     &       " field in the type 05 input card is the index number of",/  &  
     & ,                                                                 &  
     &       " the first limiter to be defined on the card.  Up to 3 ",/  &  
     & ,                                                                 &  
     &       " limiter points can still be defined on each card. An  ",/  &  
     & ,                                                                 &  
     &       " example of a proper type 05 card is as follows",//,       &  
     &       " 05        1.0       2.13      0.12      2.13      0.24",/  &  
     & /,                                                                &  
     &       " where this means xlima(1)=2.13, zlima(1)=0.12, ",/,       &  
     &       "                  xlima(2)=2.13, zlima(2)=0.24, etc"   )
   54 continue
      nlim0       = nlim
      if(card(4).le.0) go to 55
      nlim        = nlim + 1
      if(card(6).le.0) go to 55
      nlim        = nlim+1
   55 continue
      nlim_max = max( nlim_max, nlim)
      if(ineg.ne.0) go to 999
      go to 5
!...........................................................
!
!...card  6 : divertor
!
!...........................................................
   60 continue
      nsepmax = card(7)
      if(nsepmax.le.0) nsepmax = 2
      go to 5
!...........................................................
!
!...card  7 :  impurities
!
!...........................................................
   70 continue
!
      nthe = int(card(7))
      if (nthe .eq. 0) nthe = 100
      go to 5
!...........................................................
!
!...card  8 : observation pairs
!
!...........................................................
   80 continue
!
      kobs = 2*int(card(1))-1
   81 kobs = kobs + 1
   82 if(kobs.ge.nobs) nobs = kobs
      go to 5
!...........................................................
!
!...card  9 : external coils
!
!...........................................................
   90 continue
!
      ic0 = int(card(1))
      if(ic0.le.0) ineg=90
!
      if(ic0.ge.4000) ngr=ic0-4000
      if(ic0.ge.3000) ngr=ic0-3000
      if(ic0.ge.2000) ngr=ic0-2000
      if(ic0.ge.1000) ngr=ic0-1000
    
      if(ic0 .lt. 1000 .and. ic0.ge.ncoil) ncoil = ic0
      ngrmax = max(ngrmax,ngr)
      if(ineg.gt.0) go to 999
      go to 5
!...........................................................
!
!...card 10 : internal coils
!
!...........................................................
  100 continue
!
      ic0          = abs(card(1))
      ngr = int(card(4))
!     if(ic0.gt.5000) ineg=21
      if(ic0.gt.4000 .and. card(1).gt.0) ngr = ic0 - 4000
      if(ic0.gt.3000 .and. card(1).gt.0) ngr = ic0 - 3000
      if(ic0.gt.2000 .and. card(1).gt.0) ngr = ic0 - 2000
      if(ic0.gt.1000 .and. card(1).gt.0) ngr = ic0 - 1000
      if(ic0.le.0) ineg=21
      if((ic0 .lt. 1000 .or. card(1) .lt. 0) .and. ic0.ge.nwire) nwire = ic0
      igroupw_max = max(igroupw_max,ngr)
      go to 5
101   go to 5
103   go to 5
105   go to 5
107   go to 5
!...........................................................
!
!...card 11 : acoef array
!
!...........................................................
  110 continue
      go to 5
!...........................................................
!
!...card 12 : transport
!
!...........................................................
  120 continue
      go to 5
!...........................................................
!
!...card 13 : initial conditions
!
!...........................................................
  130 continue
       go to 5
!...........................................................
!
!...card 14 : more initial conditions
!
!...........................................................
  140 continue
      go to 5
!...........................................................
!
!...card 15 : coil groups
!
!...........................................................
  150 continue
!
      igroup = int(abs(card(1)))
      if(igroup.gt.ngrmax) ngrmax = igroup
      go to 5
!...........................................................
!
!...card 16 : plasma current
!
!...........................................................
  160 continue
      go to 5
!...........................................................
!
!...card 17 : plasma pressure
!
!...........................................................
  170 continue
      go to 5
!...........................................................
!
!...card 18 : timing
!
!...........................................................
  180 continue
!
      do 181 l=1,6
      ll = l + l18z
      if(card(l+1).ne.0) ntpts = ll
  181 continue
      l18z = l18z + 6
      go to 5
!...........................................................
!
!...card 19 : radial feedback  - 1
!
!...........................................................
  190 continue
!
      if(int(card(1)).eq.1000) go to 191
      inumfb = int(card(1))
      if(inumfb.gt.numfb) numfb = inumfb
      ngr = int(card(2))
      ngrmax = max(ngrmax,ngr)
      go to 5
  191 continue
      go to 5
!...........................................................
!
!...card 20 : radial feedback - 2
!
!...........................................................
  200 continue
      inumfb = int(card(1))
      go to 5
!...........................................................
!
!...card 21 : optional contour plots
!
!...........................................................
  210 continue
      go to 5
!...........................................................
!
!...card 22 : optional vector plots
!
!...........................................................
  220 continue
      go to 5
!..............................................................
!
!...card 23 : neutral beams
!
!...............................................................
  230 continue
      l23z = l23z + 6
      if(l23z.gt.ntpts) ntpts=l23z
      go to 5
!.................................................................
!
!.....card 24 : normalized density
!
!....................................................................
  240 continue
      l24z = l24z + 6
      if(l24z.gt.ntpts) ntpts=l24z
      go to 5
!.......................................................................
!
!.....card 25 : beam deposition profile
!
!.......................................................................
  250 continue
      go to 5
!.......................................................................
!
!.....card 26 : transport multiplier
!
!.......................................................................
  260 continue
      l26z = l26z + 6
      if(l26z.gt.ntpts) ntpts=l26z
      go to 5
!.......................................................................
!
!.....card 27 : time dependent external toroidal field
!
!.......................................................................
  270 continue
      l27z = l27z + 6
      if(l27z.gt.ntpts) ntpts=l27z
      go to 5
!.......................................................................
!
!.....card 28 : pre programmed loop voltage
!
!.......................................................................
  280 continue
      l28z = l28z + 6
      if(l28z.gt.ntpts) ntpts=l28z
      go to 5
!.......................................................................
!
!.....card 29 : Pest output times
!
!.......................................................................
290   continue
      l29z=l29z+6
      go to 5
!.......................................................................
!
!.....card 30 : preprogrammed x-position of magnetic axis
!
!.......................................................................
 300  continue
      l30z=l30z+6
      if(l30z.gt.ntpts) ntpts=l30z
      go to 5
!.......................................................................
!
!.....card 31 : preprogrammed z-position of magnetic axis
!
!.......................................................................
 310  continue
      l31z=l31z+6
      if(l31z.gt.ntpts) ntpts=l31z
      go to 5
!.....................................................................
!
!.....card 32 : Divertor Plate info
!
!......................................................................
  320 continue
      ic0 = int(card(1))
      if(ic0.ge.1000) go to 321
      iplate = 1
      l32z = 0
      if(ic0.gt.nplate) nplate = ic0
      go to 5
  321 continue
!     do 322 l=1,3
!     if(card(2*l).le.0) go to 322
!     if(l+l32z .gt. pnseg+1) go to 322
! 322 continue
      l32z = l32z + 3
      go to 5
!.....................................................................
!
!.....card 33 : additional coil group informations
!
!.....................................................................
!     card(1) : coil group number
!     card(2) : gap resistance in ohms (resgs)
!     card(3) : initial energy in coil group
!
 330  continue
      go to 5
!.......................................................................
!
!.....card 34 : preprogrammed vacuum temperature tevv
!.......................................................................
  340 continue
      l34z = l34z+6
      if(l34z.gt.ntpts) ntpts=l34z
      go to 5
!.......................................................................
!
!.....card 35: preprogrammed mass enhancement ffac
!
!.......................................................................
  350 continue
      l35z = l35z+6
      if(l35z.gt.ntpts) ntpts=l35z
      go to 5
!.......................................................................
!
!.....card 36: preprogrammed z-effective zeff
!
!.......................................................................
  360 continue
      l36z = l36z+6
      if(l36z.gt.ntpts) ntpts=l36z
      go to 5
!......................................................................
!
!.....card 37: preprogrammed group voltage
!
!......................................................................
  370 continue
      igroup = int(abs(card(1)))
      go to 5
!....................................................................
!
!.....card 38: lower hybrid heating and current drive parameters
!
!....................................................................
  380 continue
      ilhcd = int(card(1))
      go to 5
!..................................................................
!
!.....card 39: additional external coil info
!
!..................................................................
  390 continue
      ic0 = int(card(1))
      go to 5
!
!..................................................................
!
!.....card 40: plot output suppression (noplot)
!
!..................................................................
!
 400  continue
      goto 5
!..................................................................
!
!.....card 41: TF ripple parameters
!
!..................................................................
  410 continue
      go to 5
!...........................................................................
!
!.....card 42 : preprogrammed center of last flux surface
!
!.......................................................................
 420  continue
      l42z=l42z+6
      if(l42z.gt.ntpts) ntpts=l42z
      go to 5
!.......................................................................
!
!.....card 43 : preprogrammed minor radius
!
!.......................................................................
 430  continue
      l43z=l43z+6
      if(l43z.gt.ntpts) ntpts=l43z
      go to 5
!.......................................................................
!
!.....card 44 : preprogrammed ellipticity
!
!.......................................................................
 440  continue
      l44z=l44z+6
      if(l44z.gt.ntpts) ntpts=l44z
      go to 5
!.......................................................................
!
!.....card 45 : preprogrammed triangularity
!
!.......................................................................
 450  continue
      l45z=l45z+6
      if(l45z.gt.ntpts) ntpts=l45z
      go to 5
!
!..................................................................
!
!.....CARD 46:  preprogramed Lower Hybrid heating powers
!
!...................................................................
  460 continue
      l46z=l46z+6
      if(l46z.gt.ntpts) ntpts=l46z
      go to 5
!..................................................................
!
!.....CARD 47:  preprogrammed density exponent --- 1
!
!...................................................................
  470 continue
      l47z=l47z+6
      if(l47z.gt.ntpts) ntpts=l47z
      go to 5
!..................................................................
!
!.....CARD 48:  preprogrammed density exponent
!
!...................................................................
  480 continue
      l48z=l48z+6
      if(l48z.gt.ntpts) ntpts=l48z
      go to 5
!...................................................................
!
!.....card 49 : Multipole "coils"
!
!...................................................................
  490 continue
      ic0 = int(card(1))
      if(ic0 .ge. nmult)nmult = ic0
      if(ic0 .le. 0) ineg=490
      if(ineg .ne. 0) go to 999
      go to 5
!..................................................................
!
!.....CARD 50:  fraction_of_beam_parallel_to_B
!
!...................................................................
  495 continue
      l50z=l50z+6
      if(l50z.gt.ntpts) ntpts=l50z
      go to 5
!..................................................................
!
!.....CARD 51:  position of the peak of the lh power profile
!
!...................................................................
 5101 continue
      l51z=l51z+6
      if(l51z.gt.ntpts) ntpts=l51z
      go to 5
!..................................................................
!
!.....CARD 52:  width of the lh power profile
!
!...................................................................
 5201 continue
      l52z=l52z+6
      if(l52z.gt.ntpts) ntpts=l52z
      go to 5
!..................................................................
!
!.....CARD 53:  exponent value near center of the lh power profile
!
!...................................................................
 5301 continue
      l53z=l53z+6
      if(l53z.gt.ntpts) ntpts=l53z
      go to 5
!..................................................................
!
!.....CARD 54:  exponent value at edge of the lh power profile
!
!...................................................................
 5401 continue
      l54z=l54z+6
      if(l54z.gt.ntpts) ntpts=l54z
      go to 5
!..................................................................
!
!.....CARD 55:  position of the peak of the current profile
!
!...................................................................
 5501 continue
      l55z=l55z+6
      if(l55z.gt.ntpts) ntpts=l55z
      go to 5
!..................................................................
!
!.....CARD 56:  width of the current profile
!
!...................................................................
 5601 continue
      l56z=l56z+6
      if(l56z.gt.ntpts) ntpts=l56z
      go to 5
!..................................................................
!
!.....CARD 57:  exponent value near center for current profile
!
!...................................................................
 5701 continue
      l57z=l57z+6
      if(l57z.gt.ntpts) ntpts=l57z
      go to 5
!..................................................................
!
!.....CARD 58:  exponent value at edge of the current profile
!
!...................................................................
 5801 continue
      l58z=l58z+6
      if(l58z.gt.ntpts) ntpts=l58z
      go to 5
!...............................................................
 5901 continue
      iicrh = 1
      l59z=l59z+6
      if(l59z.gt.ntpts) ntpts=l59z
      go to 5
!...............................................................
 6001 continue
      l60z=l60z+6
      if(l60z.gt.ntpts) ntpts=l60z
      go to 5
!...............................................................
 6101 continue
      l61z=l61z+6
      if(l61z.gt.ntpts) ntpts=l61z
      go to 5
!..............................................................
 6201 continue
      igroup = int(abs(card(1)))
      go to 5
!..............................................................
 6301 continue
      igroup = int(abs(card(1)))
      go to 5
!
!.....................................................................
!.....CARD 64:  FWCD Current
!
!...................................................................
 6401 continue
      ifwcd = 1
      l64z=l64z+6
      if(l64z.gt.ntpts) ntpts=l64z
      go to 5
!..................................................................
!
!.....CARD 65:  position of the peak of the fw power profile
!
!...................................................................
 6501 continue
      l65z=l65z+6
      if(l65z.gt.ntpts) ntpts=l65z
      go to 5
!..................................................................
!
!.....CARD 66:  width of the fw power profile
!
!...................................................................
 6601 continue
      l66z=l66z+6
      if(l66z.gt.ntpts) ntpts=l66z
      go to 5
!..................................................................
!
!.....CARD 67:  exponent value near center of the fw power profile
!
!...................................................................
 6701 continue
      l67z=l67z+6
      if(l67z.gt.ntpts) ntpts=l67z
      go to 5
!..................................................................
!
!.....CARD 68:  exponent value at edge of the fw power profile
!
!...................................................................
 6801 continue
      l68z=l68z+6
      if(l68z.gt.ntpts) ntpts=l68z
      go to 5
!..................................................................
!
!.....CARD 69:  position of the peak of the current profile
!
!...................................................................
 6901 continue
      l69z=l69z+6
      if(l69z.gt.ntpts) ntpts=l69z
      go to 5
!..................................................................
!
!.....CARD 70:  width of the current profile
!
!...................................................................
 7001 continue
      l70z=l70z+6
      if(l70z.gt.ntpts) ntpts=l70z
      go to 5
!..................................................................
!
!.....CARD 71:  exponent value near center for current profile
!
!...................................................................
 7101 continue
      l71z=l71z+6
      if(l71z.gt.ntpts) ntpts=l71z
      go to 5
!..................................................................
!
!.....CARD 72:  exponent value at edge of the current profile
!
!...................................................................
 7201 continue
      l72z=l72z+6
      if(l72z.gt.ntpts) ntpts=l72z
      go to 5
!..................................................................
!
!.....CARD 73:  helium time constant
!
!...................................................................
 7301 continue
      l73z=l73z+6
      if(l73z.gt.ntpts) ntpts=l73z
      go to 5
!
!....CARD 74 special input for writing ufiles of impurity charge states
 7401 continue
      if(card(1).le.0) go to 5
      do 7402 l=1,6
      itempi = int(card(1+l) + 0.1_R8)
      if(itempi.le.0) go to 7402
      numicsuf = numicsuf + 1
 7402 continue
      if(numicsuf .gt. 18) ineg=7401
      go to 5
!.......................................................................
!
!.....card 75 : Sawtooth times for isaw=2
!
!.......................................................................
 7501 continue
      do 7502 l=1,6
      lll=l+l75z
      if(lll.gt.500) ineg=75
7502   continue
      l75z=l75z+6
      if(ineg .ne. 0) go to 999
      go to 5
!...............................................................
!
!
!
!===========================================================
!
!....... end card read
!
!===========================================================
!
  500 continue
!
      if(isym .eq. 0) nthe = 2*(nthe/2._R8)
!
      if(nmult.le.0) go to 516
!.....add multipolar coils to regular coils
      do 515 l=1,nmult
      ncoil = ncoil + 1
  515 continue
  516 continue
!...........................................................
      if(nwire.eq.0) go to 605
      do 604 n=1,nwire
      ncoil = ncoil + 1
  604 continue
  605 continue
      pnx = nx + 4
      pnz = nz + 4
      pnwire =nwire
      pncoil = ncoil
      ppsi = npsi + 1
      pnthe = nthe + 4
      pobs = nobs
      pngroup = ngrmax
      pngroup = max( pngroup, igroupw_max)
      pnlim = nlim_max
      pnplat = nplate
      pnsep = nsepmax
      ptpts = ntpts
      pnfeed = numfb
      pncoil = pncoil + 2
      pnwire = pnwire + 2
!...........................................................
!     define size based on restart file if a restart
!     checking is not done here

      if( irst1 .eq. 1 ) then
              if( numargs .lt. 1 ) then
                   isprsi = 'sprsin' // isuffix(1:1)
              else
                   isprsi = 'sprsin' // '.' // trim(suffix)
              end if
              open(nsprsi,file=trim(isprsi),status='old',iostat=ios5,    &  
     &        form='unformatted')
              if( ios5 .ne. 0)                                           &  
     &              stop "Error opening restart file : problem_size"
              if(.not.allocated(nnp))allocate(nnp(27),STAT=istat)
              nnp=0
              call bufini(nsprsi,nnp(1),nnp(27))
              pnx = nnp(3)
              pnz = nnp(4)
              pncoil = nnp(5)
              pobs = nnp(6)
              pngroup = nnp(7)
              pnlim = nnp(8)
              ptpts = nnp(9)
              pnthe = nnp(13)
              ppsi = nnp(14)
              pnplat = nnp(22)
              pnfeed = nnp(25)
              pnsep = nnp(24)
              pnwire = nnp(27)
              nnp = 0
              close(nsprsi)
      end if
!...........................................................
      penx = pnx
      penz = pnz
      pnthew = pnthe 
      pboun=4*penx+8*penz
!     pnparts=pnx+2*pnz
      pnparts=320
      pgls=pglobs+ppsi+4*pnplat + pnplat*pnseg
      pglob=pgls+3*pncoil+10*pngroup+6
      pglobp=pglob+1
      pnxz=pnx*pnz
!...........................................................
      nbig=22*pnx*pnz
      pnxz1=  pnx*pnz+1
      pnxz2=2*pnx*pnz+1
      pnxz3=3*pnx*pnz+1
      pnxz4=4*pnx*pnz+1
      pnxz5=5*pnx*pnz+1
      pnxz6=6*pnx*pnz+1
      pnxz7=7*pnx*pnz+1
      pnxz8=8*pnx*pnz+1
      pnxz9=9*pnx*pnz+1
      pnxz10=10*pnx*pnz+1
      pnxz11=11*pnx*pnz+1
      pnxz12=12*pnx*pnz+1
      pnxz13=13*pnx*pnz+1
      pnxz14=14*pnx*pnz+1
      pnxz15=15*pnx*pnz+1
      pnxz16=16*pnx*pnz+1
      pnxz17=17*pnx*pnz+1
      pnxz18=18*pnx*pnz+1
      pnxz19=19*pnx*pnz+1
      pnxz20=20*pnx*pnz+1
      pnxz21=21*pnx*pnz+1

      pnxx1=pnx*pnx+1
!...........................................................
      pnplts=pnplt+pnplat
!...........................................................
  620 continue
!
      rewind (unit=nin)
      return
!
!
  999 continue
      ineg=21
      write(nout,1999)
 1999 format(" one or more fields on preceeding card is",                &  
     &       " outside specified range")
!
      stop
 1001 format(" type",i2,2x,1p9e11.3)
 2000 format(8a10)
 2001 format(20x,8a10)
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"
! 15Apr2005 fgtok

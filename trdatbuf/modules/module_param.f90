      MODULE PARAM
      IMPLICIT NONE
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER iformat,nordp,nptsp
!============
      integer pboun,                                                     &  
     &        pglob, pglobp, pnsave, pneq, pfour, pkw,                   &  
     &        plw, pimp, pne, pte,pnode,pklw,                            &  
     &        pdelay, pgls, p50, pspda,pnseg,pglobs,ptwir,               &  
     &        pchrg,pnframe,pnxz,pnluup,picur,psaw,                      &  
     &        pnshap,pnvvel,pmmay,pnprob,pnfob,pnbob,pnob,pngrps,        &  
     &        pchrgmx,pnparts,                                           &  
     &        pcgOx,pcgC,pcgFe,pcgBe,pcgNe,pcgKr,pcgAr,pcgW
!
      INTEGER :: penx = 110, penz = 105
      INTEGER :: pnx  = 110, pnz  = 105
      INTEGER :: pncoil=600, pnwire=999
      INTEGER :: pnfeed=30
      INTEGER :: pobs=4*26
      INTEGER :: pngroup=60 
      INTEGER :: pnlim=101
      INTEGER :: ptpts=4*6
      INTEGER :: pnsep=10
      INTEGER :: pnthew=20, pnthe=110
      INTEGER :: ppsi=801
      INTEGER :: pnplat=60
!...............description of parameters: .........................
!
!......penx-----max no of equilib. zones in x-direction + 4
!     parameter(penx=91)
!
!......penz-----max no of equilib. zones in z-direction + 4
!     parameter(penx=110,penz=105)
      
!
!......pnx------max no of eulerian zones in x-direction + 4
!     parameter(pnx=penx)
!
!......pnz------max no of eulerian zones in z-direction + 4
!     parameter(pnz=penz)
!
!......pncoil---max no of coils (type9+type10 cards) + 2
!     parameter(pncoil=600)
!
!......pnwire---max no of wires(type10 cards) + 2
!     parameter(pnwire=999)
!
!......pdelay---max no of time delay steps on type 19 card
      parameter(pdelay=200)
!
!.....pnfeed---max no of feedback systems defined on type 19 card
!     parameter(pnfeed=30)
!
!......pobs-----max no of observation points  (type 8)
!     parameter(pobs=4*26)
!
!......pngroup--max no of coil groups(type 15)
!     parameter(pngroup=60)
!
!......pngrps---max no of coil groups  with series connections
      parameter(pngrps=20)
!
!......pnlim----max no of limiter points(type 5)
!     parameter(pnlim=101)
!
!......ptpts----no of time-points (type 18 cards)
!     parameter(ptpts=4*6)
!
!......pnsave---no of points saved for time history
      parameter(pnsave=1000)
!......pnsep---max number of separatrixes searched for
!     parameter(pnsep=10)
!
!......pneq-----number of adiabatic variables
      parameter(pneq=4)
!
!......pnthew---no of theta zones for par trans(for iwall=1, type 07)
!     parameter(pnthew=20)
!
!......pnthe----number of theta zones for balloon stability
!     parameter(pnthe=110)
!
!......ppsi-----no of psi zones for surface averaged transport( type 04)
!     parameter(ppsi=801)
!
!......pfour----number of surface averaged variables advanced
      parameter(pfour=4)
!
!......pspda----number of time points read by spdata when idata=1,2
      parameter(pspda=10000)
!
!......pnplat---number of divertor plates for heat calculation(type 32)
!     parameter(pnplat=60)
!
!......pnseg---number of subdivisions for each divertor plate
      parameter(pnseg=100)
!
!......pnframe---number of movie frames when imovie .gt. 0
      parameter (pnframe=1000)
!
!......iformat..=1 regular input    ..=2 lausanne input
      parameter(iformat=1)
!
!......ptwir----max number of time points in wire file(for iwayne.gt.0 )
      parameter(ptwir=1000)
!......psaw ---- maximum number of sawtooth oscillations
      parameter(psaw=500)
!
!
!......parameters for svd routine
      parameter(nordp=7)
      parameter(nptsp=9)
!
!......pkw------no of walls
!......plw------no of segments in wall
!......pimp-----no of impurity species
!......pne------no of density locations in radiation arrays
!......pte------no of temp locations in radiation arrays
!......pchrg----charge state in radiation arrays (block data)
!......pchrgmx--maximum charge state for impurity tranpsort (nucz+1)
      parameter(pkw=1,plw=1,pnode=04,pchrg=1)
      parameter(pne=5,pte=25,pcgOx=9,pcgC=7,pcgFe=27,pcgBe=5,pcgNe=11)
      parameter(pcgKr=37,pcgAr=19,pcgW=75)
      parameter(pimp=08,pchrgmx=75)
!     parameter(pklw=pnlim)
      parameter(pklw=101)
      parameter(pnluup=50,picur=75,pnshap=25,pnvvel=50,pmmay=20)
      parameter(pnprob=50,pnfob=40,pnbob=5,pnob=50)
!...................................................................
!
!
!     parameter(pboun=4*penx+8*penz)
!     parameter(pnparts=pnx+2*pnz)
!     parameter(pglobs=250,pgls=pglobs+ppsi+4*pnplat + pnplat*pnseg)
      parameter(pglobs=250)
!     parameter(pglob=pgls+3*pncoil+10*pngroup+6,pglobp=pglob+1)
!     ------------------------------------------------------------------
!     parameter(pnxz=pnx*pnz)

!
! 15Apr2005 fgtok -s r8_precision.sub "r8con_incl.csh conversion"
      END MODULE PARAM

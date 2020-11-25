c
      integer pnx, pnz, pncoil, pboun, pnwire, pobs, pngroup, pnlim,
     +        pglob, pglobp, pnsave, pneq, pnthe, ppsi, pfour, pkw,
     +        plw, pimp, pne, pte,pnode,pklw,penx,penz,ptpts,
     +        pdelay, pgls, p50, pspda,pnplat,pnseg,pglobs,ptwir,
     +        pnthew,pchrg,pnframe,pnxz,pnfeed,pnluup,picur,psaw,
     +        pnshap,pnvvel,pmmay,pnprob,pnfob,pnbob,pnob,pngrps,    
     +        pnsep,pchrgmx,pnxz1,pnxz2,pnxz3,pnxx1,pnparts,
     +        pcgOx,pcgC,pcgFe,pcgBe,pcgNe,pcgKr,pcgAr,pcgW      
c
c...............description of parameters: .........................
c
c......penx-----max no of equilib. zones in x-direction + 4
c     parameter(penx=91)
c
c......penz-----max no of equilib. zones in z-direction + 4
      parameter(penx=150,penz=235)
c
c......pnx------max no of eulerian zones in x-direction + 4
      parameter(pnx=penx)
c
c......pnz------max no of eulerian zones in z-direction + 4
      parameter(pnz=penz)
c
c......pncoil---max no of coils (type9+type10 cards) + 2
      parameter(pncoil=999)
c
c......pnwire---max no of wires(type10 cards) + 2
      parameter(pnwire=999)
c
c......pdelay---max no of time delay steps on type 19 card
      parameter(pdelay=200)
c
c.....pnfeed---max no of feedback systems defined on type 19 card
      parameter(pnfeed=30)
c
c......pobs-----max no of observation points  (type 8)
      parameter(pobs=4*26)
c
c......pngroup--max no of coil groups(type 15)
      parameter(pngroup=60)
c
c......pngrps---max no of coil groups  with series connections
      parameter(pngrps=20)
c
c......pnlim----max no of limiter points(type 5)
      parameter(pnlim=101)
c
c......ptpts----no of time-points (type 18 cards)
      parameter(ptpts=4*6)
c
c......pnsave---no of points saved for time history
      parameter(pnsave=1000)
c......pnsep---max number of separatrixes searched for
      parameter(pnsep=10)
c
c......pneq-----number of adiabatic variables
      parameter(pneq=4)
c
c......pnthew---no of theta zones for par trans(for iwall=1, type 07)
      parameter(pnthew=06)
c
c......pnthe----number of theta zones for balloon stability
      parameter(pnthe=110)
c
c......ppsi-----no of psi zones for surface averaged transport( type 04)
      parameter(ppsi=701)
c
c......pfour----number of surface averaged variables advanced
      parameter(pfour=4)
c
c......pspda----number of time points read by spdata when idata=1,2
      parameter(pspda=10000)
c
c......pnplat---number of divertor plates for heat calculation(type 32)
      parameter(pnplat=60)
c
c......pnseg---number of subdivisions for each divertor plate
      parameter(pnseg=100)
c
c......pnframe---number of movie frames when imovie .gt. 0
      parameter (pnframe=1000)
c
c......iformat..=1 regular input    ..=2 lausanne input
      parameter(iformat=1)
c
c......ptwir----max number of time points in wire file(for iwayne.gt.0 )
      parameter(ptwir=1000)
c......psaw ---- maximum number of sawtooth oscillations
      parameter(psaw=500)
c
c
c......parameters for svd routine
      parameter(nordp=7)
      parameter(nptsp=9)
c
c......pkw------no of walls
c......plw------no of segments in wall
c......pimp-----no of impurity species
c......pne------no of density locations in radiation arrays
c......pte------no of temp locations in radiation arrays
c......pchrg----charge state in radiation arrays (block data)
c......pchrgmx--maximum charge state for impurity tranpsort (nucz+1)
      parameter(pkw=1,plw=1,pnode=04,pchrg=1)
      parameter(pne=5,pte=25,pcgOx=9,pcgC=7,pcgFe=27,pcgBe=5,pcgNe=11)
      parameter(pcgKr=37,pcgAr=19,pcgW=75)
      parameter(pimp=08,pchrgmx=75)
      parameter(pklw=pnlim)
      parameter(pnluup=50,picur=75,pnshap=25,pnvvel=50,pmmay=20)
      parameter(pnprob=50,pnfob=40,pnbob=5,pnob=50)
c...................................................................
c
c
      parameter(pboun=4*penx+8*penz)
      parameter(pnparts=pnx+2*pnz)
      parameter(pglobs=250,pgls=pglobs+ppsi+4*pnplat + pnplat*pnseg)
      parameter(pglob=pgls+3*pncoil+10*pngroup+6,pglobp=pglob+1)
c     ------------------------------------------------------------------

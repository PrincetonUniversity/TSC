      subroutine surfplt2
!
!
!
!
!
!.....plots surface source and sink terms
!
      USE CLINAM
      USE SAPROP
      USE WALLCL
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER n
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 aitot
      REAL*8 aibs,aicd,aifw,aiec,alhcd,darea,fac1,hmax,hmin,ramax,ramin
      REAL*8 vmax,vmin,ajmax,ajmin,xmin,xmax,xlabset,tplcur
!============
!     dimension thary(ppsi),psiy(ppsi),psih(ppsi),
!    1    ahary(ppsi),cycary(ppsi),alhhary(ppsi),
!    2    psit(ppsi),aimpary(ppsi),breary(ppsi),viary(ppsi),
!    3    ajtary(ppsi),absary(ppsi),acdary(ppsi),alhary(ppsi),
!    1      afwary(ppsi)
!============      
      INTEGER :: istat = 0 
      REAL*8, ALLOCATABLE, DIMENSION(:) :: thary
      REAL*8, ALLOCATABLE, DIMENSION(:) :: psiy
      REAL*8, ALLOCATABLE, DIMENSION(:) :: psih
      REAL*8, ALLOCATABLE, DIMENSION(:) :: ahary
      REAL*8, ALLOCATABLE, DIMENSION(:) :: cycary
      REAL*8, ALLOCATABLE, DIMENSION(:) :: alhhary
      REAL*8, ALLOCATABLE, DIMENSION(:) :: psit
      REAL*8, ALLOCATABLE, DIMENSION(:) :: aimpary
      REAL*8, ALLOCATABLE, DIMENSION(:) :: breary
      REAL*8, ALLOCATABLE, DIMENSION(:) :: viary
      REAL*8, ALLOCATABLE, DIMENSION(:) :: ajtary
      REAL*8, ALLOCATABLE, DIMENSION(:) :: absary
      REAL*8, ALLOCATABLE, DIMENSION(:) :: acdary
      REAL*8, ALLOCATABLE, DIMENSION(:) :: alhary
      REAL*8, ALLOCATABLE, DIMENSION(:) :: afwary
      REAL*8, ALLOCATABLE, DIMENSION(:) :: aecary
      REAL*8, ALLOCATABLE, DIMENSION(:) :: aechary
      REAL*8, ALLOCATABLE, DIMENSION(:) :: abmhary
      REAL*8, ALLOCATABLE, DIMENSION(:) :: afwhary
!============      
      IF(.not.ALLOCATED(thary)) ALLOCATE( thary(ppsi), STAT=istat)
      IF(.not.ALLOCATED(psiy)) ALLOCATE( psiy(ppsi), STAT=istat)
      IF(.not.ALLOCATED(psih)) ALLOCATE( psih(ppsi), STAT=istat)
      IF(.not.ALLOCATED(ahary)) ALLOCATE( ahary(ppsi), STAT=istat)
      IF(.not.ALLOCATED(cycary)) ALLOCATE( cycary(ppsi), STAT=istat)
      IF(.not.ALLOCATED(alhhary)) ALLOCATE( alhhary(ppsi), STAT=istat)
      IF(.not.ALLOCATED(psit)) ALLOCATE( psit(ppsi), STAT=istat)
      IF(.not.ALLOCATED(aimpary)) ALLOCATE( aimpary(ppsi), STAT=istat)
      IF(.not.ALLOCATED(breary)) ALLOCATE( breary(ppsi), STAT=istat)
      IF(.not.ALLOCATED(viary)) ALLOCATE( viary(ppsi), STAT=istat)
      IF(.not.ALLOCATED(ajtary)) ALLOCATE( ajtary(ppsi), STAT=istat)
      IF(.not.ALLOCATED(absary)) ALLOCATE( absary(ppsi), STAT=istat)
      IF(.not.ALLOCATED(acdary)) ALLOCATE( acdary(ppsi), STAT=istat)
      IF(.not.ALLOCATED(alhary)) ALLOCATE( alhary(ppsi), STAT=istat)
      IF(.not.ALLOCATED(afwary)) ALLOCATE( afwary(ppsi), STAT=istat)
      IF(.not.ALLOCATED(aecary)) ALLOCATE( aecary(ppsi), STAT=istat)
      IF(.not.ALLOCATED(aechary)) ALLOCATE( aechary(ppsi), STAT=istat)
      IF(.not.ALLOCATED(abmhary)) ALLOCATE( abmhary(ppsi), STAT=istat)
      IF(.not.ALLOCATED(afwhary)) ALLOCATE( afwhary(ppsi), STAT=istat)
!============      
      if (istat .ne. 0) stop 'Allocation Error : surfplt2  ' 
!============      
!
      if(dpsi.le.0) return
      npts = nx/2
      if(npsit.gt.1) npts = npsit-1
      aitot = 0._R8
      aibs  = 0._R8
      aicd  = 0._R8
      aifw = 0._R8
      aiec = 0._R8
      alhcd = 0._R8
!
      do 100 n=2,npts+1
!     thary(n) = (savei(n) + savee(n) + (savia(n)+savea(n))*ialpha       &  
!    &         + savefw(n) + savebm(n) + savibm(n) + savifw(n)           &  
!    &         + savelh(n) + savilh(n) )*vp(n)*udsp/udst
!     ahary(n) = (savia(n) + savea(n))*vp(n)*udsp/udst
!     afwhary(n) = (savefw(n)+savifw(n))*vp(n)*udsp/udst
!     abmhary(n) = (savei(n) + savee(n))*vp(n)*udsp/udst
      thary(n) = (savei(n) + savee(n) + (savia(n)+savea(n))*ialpha       &  
     &         + savefw(n) + savebm(n) + savibm(n) + savifw(n)           &  
     &         + savelh(n) + savilh(n) )*udsp/udst
      ahary(n) = ialpha*(savia(n) + savea(n))*udsp/udst
      afwhary(n) = (savefw(n)+savifw(n))*udsp/udst
      abmhary(n) = (savei(n) + savee(n))*udsp/udst
      aechary(n) = (savebm(n) + savibm(n))*udsp/udst
      alhhary(n) = (savelh(n) + savilh(n))*udsp/udst 
!
!     breary(n) = savebre(n)*vp(n)*udsp/udst
!     cycary(n) = savecyc(n)*vp(n)*udsp/udst
!     aimpary(n) =( sradion(n) + saveimp(n)*udsp/udst)*vp(n)
!     alhhary(n) = savelh(n)*vp(n)*udsp/udst
      breary(n) = savebre(n)*udsp/udst
      cycary(n) = savecyc(n)*udsp/udst
      aimpary(n) =( sradion(n) + saveimp(n)*udsp/udst)
!
      viary(n) = .5_R8*(as(n)+as(n-1))*udsv
!
      darea = vp(n)/(tpi*xplas)
      ajtary(n) = (gxmja2(n)-gxmja2(n-1))*rdpsi*udsi/tpi
      absary(n) = .5_R8*(ajavbs(n)+ajavbs(n-1))*udsi
      acdary(n) = .5_R8*(ajavcd(n)+ajavcd(n-1))*udsi
      afwary(n) = ajavfw(n)*udsi
      aecary(n) = ajavec(n)*udsi
      alhary(n) = ajavlh(n)*udsi
      fac1 = .5_R8*(ajbtsq(n)/ajbsq(n) + ajbtsq(n-1)/ajbsq(n-1))
      aitot = aitot + ajtary(n)*dpsi
      aibs  = aibs  + absary(n)*fac1*dpsi
      aicd  = aicd  + acdary(n)*fac1*dpsi
      aifw  = aifw  + afwary(n)*fac1*dpsi
      aiec  = aiec  + aecary(n)*fac1*dpsi
      alhcd = alhcd + alhary(n)*fac1*dpsi
! 100 psit(n) = (xsv(n)-psimin)*tpi
! 100 psit(n) = rminora(n)
  100 psit(n) = sqrt((float(n-1)*dpsi)/(float(npsit-1)*dpsi))

!     alhhary(1) = savelh(1)*udsp/udst
!     afwhary(1) = (savefw(1)+savifw(1))*udsp/udst
!     abmhary(1) = (savei(1)+savee(1))*udsp/udst
!     ahary(1) = (savia(1)+savea(1))*udsp/udst
!     thary(1) = alhhary(1)+afwhary(1)+abmhary(1)+ahary(1)
!
      hmax = thary(npts+1)
      hmin = thary(npts+1)
      ramax = breary(npts+1)
      ramin = breary(npts+1)
      vmax = viary(npts+1)
      vmin = viary(npts+1)
      ajmax = max(ajtary(2),ajtary(npts+1))
      ajmin = min(ajtary(2),ajtary(npts+1))
      xmin = 0._R8
      xmax = psit(npts+1)
      do 200 n=2,npts
      hmax = max(hmax,thary(n),ahary(n))
      hmin = min(hmin,thary(n),ahary(n))
      ramax = max(ramax,breary(n),cycary(n),aimpary(n))
      ramin = min(ramin,breary(n),cycary(n),aimpary(n))
      ajmax = max(ajmax,ajtary(n),absary(n),acdary(n),alhary(n),         &  
     &      afwary(n),aecary(n))
      ajmin = min(ajmin,ajtary(n),absary(n),acdary(n),alhary(n),         &  
     &      afwary(n),aecary(n))
      vmax = max(vmax,viary(n))
      vmin = min(vmin,viary(n))
  200 continue
      if(hmax.le.hmin) hmax = hmin+1._R8
      if(ramax.le.ramin) ramax = ramin+1._R8
      if(ajmax.le.ajmin) ajmax = ajmin+1._R8
      if(vmax .gt. 3.0_R8) vmax = 3.0_R8
      if(vmin .lt. -2._R8) vmin = -2._R8
      if(vmax.le.vmin) vmax = vmin+1._R8
!
      call mapg(xmin,xmax,hmin,hmax,.1_R8,.47_R8,.7_R8,1._R8)
      call tracec(1ht,psit(2),thary(2),npts,-1,-1,0._R8,0._R8)
      call tracec(1ha,psit(2),ahary(2),npts,-1,-1,0._R8,0._R8)
      call tracec(1hl,psit(2),alhhary(2),npts,-1,-1,0._R8,0._R8)
      call tracec(1hf,psit(2),afwhary(2),npts,-1,-1,0._R8,0._R8)
      call tracec(1hn,psit(2),abmhary(2),npts,-1,-1,0._R8,0._R8)
      call tracec(1he,psit(2),aechary(2),npts,-1,-1,0._R8,0._R8)
      xlabset = xmin-(xmax-xmin)*.20_R8
      call setold(xlabset,hmin,1,0,1,1)
      write(s100,1001)
      call gtext(s100,80,0)
 1001 format("aux heating vs psi")
!     write(nterm,9333)
!9333 format(" after aux heating")
!
      call mapg(xmin,xmax, ramin,ramax,.59_R8,.96_R8,.7_R8,1._R8)
      call tracec(1hb,psit(2),breary(2),npts,-1,-1,0._R8,0._R8)
      call tracec(1hc,psit(2),cycary(2),npts,-1,-1,0._R8,0._R8)
      call tracec(1hi,psit(2),aimpary(2),npts,-1,-1,0._R8,0._R8)
      call setold(xmax+(xmax-xmin)*.05_R8,ramin,1,0,1,1)
      write(s100,1002)
      call gtext(s100,80,0)
 1002 format("radiation vs psi")
!
      call mapg(xmin,xmax,ajmin,ajmax,.1_R8,.47_R8,.35_R8,.65_R8)
      call tracec(1ht,psit(2),ajtary(2),npts,-1,-1,0._R8,0._R8)
      call tracec(1hb,psit(2),absary(2),npts,-1,-1,0._R8,0._R8)
      call tracec(1hn,psit(2),acdary(2),npts,-1,-1,0._R8,0._R8)
      call tracec(1hf,psit(2),afwary(2),npts,-1,-1,0._R8,0._R8)
      call tracec(1he,psit(2),aecary(2),npts,-1,-1,0._R8,0._R8)
      call tracec(1hl,psit(2),alhary(2),npts,-1,-1,0._R8,0._R8)
      call setold(xlabset,ajmin,1,0,1,1)
      write(s100,1003)
      call gtext(s100,80,0)
 1003 format("Curr Dens (Amps/W-tor flx)")
!
      call setold(xlabset,ajmin-.2_R8*(ajmax-ajmin),1,0,1,0)
      write(s100,1005) kcycle,times
      call gtext(s100,80,0)
 1005 format("  kcycle=",i7,"  time=",1pe13.5)
      write(s100,1006) aitot,aibs,aicd,alhcd,aifw,aiec,ajpress
      call gtext(s100,80,0)
      write(nout,1006) aitot,aibs,aicd,alhcd,aifw,aiec,ajpress
 1006 format("SA: I",1pe9.2," BS",1pe9.2,                                &  
     &  " NB",1pe9.2," LH",1pe9.2," FW",1pe9.2," EC",1pe9.2," PS",1pe9.2)
!
      tplcur = tcurdtp*tpi*udsi
      write(s100,1007) tplcur,bsitot,cditot,alhitot,fwitot,ecitot,diatot
      call gtext(s100,86,0)
      write(nout,1007) tplcur,bsitot,cditot,alhitot,fwitot,ecitot,diatot
 1007 format("2D: I",1pe9.2," BS",1pe9.2,                                &  
     &  " NB",1pe9.2," LH",1pe9.2," FW",1pe9.2," EC",1pe9.2," PS",1pe9.2)
!
      call mapg(xmin,xmax,vmin,vmax,.59_R8,.96_R8,.35_R8,.65_R8)
      call tracec(1hi,psit(2),viary(2),npts-1,-1,-1,0._R8,0._R8)
      call setold(xmax+(xmax-xmin)*.05_R8,vmin,1,0,1,1)
      write(s100,1004)
      call gtext(s100,80,0)
 1004 format("loop voltage vs psi")
!
!     write(nterm,9332)
!9332 format("after loop voltage")
      write(nsc1,1088) kcycle
 1088 format(" sources and sinks,cycle=",i7)
      call frscj(6)
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

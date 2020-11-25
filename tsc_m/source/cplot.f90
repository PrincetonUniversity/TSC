       subroutine cplot(a,itype)
!
!***********************************************************************
!                                                                       *
!.....makes contour plots of a scalar function                          *
!                                                                       *
!           IMOVIE = 6  option added 9 Jan 88 to plot poloidal flux     *
!           contours with a specified flux increment (PSINCR) fixed     *
!           throughout the calculation.  ROS.                           *
!***********************************************************************
!	Modified 31 Mar 93  to plot jphi, psi, g on one page   ROS
!	Modified 24 Mar 94  to plot complete halo contour      ROS
!	Modified 26 Sep 94  to plot contours in colors         ROS
!
      USE CLINAM
      USE SCRATCH
      USE POLCUR
      USE SCR13, except_this => a
      USE EQRUNS
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER itype,nconlev,icallcp,nframe,nlimnstx,imnpl,imxpl
      INTEGER jmnpl,jmxpl,i,j,k1,k2,jinc,iright,icenter,ileft,ixax
      INTEGER jyax,ilcrtest,ii,itmax,nxpl,nxv,n,it,imin,imax,jmin
      INTEGER jmax,jrcontr,k2sav,k1sav,l,ipassw
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 cval,psismal,psincr,xlimnstx,zlimnstx,ytopval
      REAL*8 x858,dx858142,xco,x142,zmin,xmax,zmax,xcosv,x858p
      REAL*8 x142p,xlabval,amin,amax,fac,aa,x1,z1,x2,z2,x3,z3,x4,z4
      REAL*8 sym,cval1sav,z1neg,zlneg,asetld,plcprnt,alamda,aiprnt
      REAL*8 aloop,psiaxis,xmin,ymin,ymax
!============
      REAL*8, DIMENSION(penx,penz) :: a
        parameter (nconlev=200)
        dimension cval(nconlev)
        data  psismal, psincr, icallcp / -2.765_R8, 0.10_R8, 0/
!
!
!     common/eqruns/findex

       character*1 nchar(3),nchar2(3)
       data nchar/"p","s","e"/
       data nchar2/ "v","c","p" /
       data nframe /0/
!
!.....special data structures for NSTX
      dimension xlimnstx(6),zlimnstx(6)
      character*1 nstxchar(6)
      data nlimnstx /6/
      data xlimnstx /0.1632_R8,0.1632_R8,0.8839_R8,                      &
     &               1.5893_R8,1.6891_R8,1.8503_R8/
      data zlimnstx /0.2540_R8,1.0287_R8,1.8961_R8,                      &
     &               1.5984_R8,1.0160_R8,0.8518_R8/
      data nstxchar /"A","B","C","D","E","F" /
!============      
 
!
	if (itype.eq.9 .and. noplot(7).gt.0)	return
      imnpl = imag
      imxpl = imag
      jmnpl = jmag
      jmxpl = jmag
!
!
      do 710 i=iminn,imaxx
      do 710 j=jminn,jmaxx
      if(iexv(i,j).eq.1) go to 710
      if(iexs(i,j).eq.1) go to 710
      if(psi(i,j) .gt. psilim) go to 710
      imnpl = min(i,imnpl)
      imxpl = max(i,imxpl)
      jmnpl = min(j,jmnpl)
      jmxpl = max(j,jmxpl)
  710 continue
      imnpl = imnpl-1
      imxpl = imxpl+1
      jmxpl = jmxpl+1
      if(isym.eq.0) jmnpl = jmnpl-1
!
!
      if(itype.eq.10) then
        imnpl = 2
        imxpl = nxp
        jmnpl = 2
        jmxpl = nzp
      endif
!
        if (imovie.eq.7)   then
                           call  dplplt
                           return
                           endif
!
        if (imovie.eq.6)   then
                           if (icallcp.ne.0)   go to 420
                           k1 = nconlev
                           cval(nconlev) = -psismal
                           do 410  i=nconlev,2,-1
  410                      cval(i-1) = cval(i) - psincr
                           icallcp = 1
!                                       k2 = index in cval corresponding to
!                                            first value in vacuum
  420                      k2 = nconlev + 1
                           do 430  i=nconlev,2,-1
                           if (cval(i).lt.-psilim)   go to 440
                           k2 = k2 - 1
  430                      continue
                           endif
!
  440 continue
!
	ytopval = 0.96_R8
      x858 = .858_R8
	dx858142 = 0.716_R8
      if(alx-ccon .gt. 2._R8*alz*.716_R8/.55_R8) go to 1212
      xco = ytopval - 0.55_R8
	dx858142 =  .55_R8*(alx-ccon)/(2._R8*alz)
	if (icplgf.gt.0 .and. dx858142.gt.0.31_R8)   then
				xco = ytopval - 0.55_R8* 0.31_R8/dx858142
				dx858142 = 0.31_R8
						   endif
      	x142 = .858_R8- dx858142
      go to 1213
 1212 continue
      x142 = .142_R8
      xco = ytopval - .716_R8*(2._R8*alz/(alx-ccon))
 1213 continue
      if(itype.eq.10) then
        x142 = x142 - 0.5
        x858 = x858 - 0.5
      endif
      jinc = (nzp-2)/((2._R8-isym)*1._R8)
      zmin = -alz/1.0_R8
      xmax = alx
      zmax = alz/1.0_R8
      xcosv = xco
!
!.....special logic to put flux, current, g plot on same page if they fit
      iright=0
	icenter = 0
      ileft =0
	ixax = 0
	jyax = 0
      if(itype.eq.8.and.(imovie.eq.0 .or. imovie.ge.10).and.lrswtch.eq.0 &
     &             .and..858_R8-x142 .le.0.4_R8) icenter=1
      if(itype.eq.1.and.(imovie.eq.0 .or. imovie.ge.10).and.lrswtch.eq.0 &
     &             .and..858_R8-x142 .le.0.4_R8) iright=1
      if(itype.eq.3.and.(imovie.eq.0 .or. imovie.ge.10).and.lrswtch.eq.0 &
     &             .and..858_R8-x142 .le.0.4_R8) then
      ileft = 1
      if(icplgf.le.0) then
						icenter = 1
						ileft = 0
						jyax = 1
						endif
      endif
!
      if(iright.eq.1) then
                           x858 = 1._R8
			   if (icplgf.le.0)   x858 = 0.88_R8
                           x142 = x858 - dx858142
                      endif
      if(icenter .eq.1)	   then
                           x858 = 0.68_R8
			   if (icplgf.le.0)   x858 = 0.48_R8
                           x142 = x858 - dx858142
                           endif
      if(ileft .eq.1) then
                           x858 = 0.36_R8
                           x142 = x858 - dx858142
			   jyax = 1
                      endif
	ilcrtest = ileft + icenter + iright
!
!..logic for line plots of poloidal flux and toroidal current
!..along vertical line through magnetic axes for case of spheromak
!..and noplot(7).gt.0
      if(iright.eq.1.and.noplot(7).gt.0) then
      x858p=x858
      x142p=x142
      x858=x858p-1.1_R8*(x858p-x142p)
      x142 = x142p-1.1_R8*(x858p-x142p)
      endif
!
	call  colora("white")
	if (ilcrtest.ne.0)				then
	call map (ccon,xmax, zmin,zmax, x142,x858, xco,ytopval)
	call  setlch (ccon, zmax + 0.015_R8*(zmax-zmin), 0,1,0, -1)
	if(ileft .eq.1 .or. (icenter.eq.1 .and. icplgf.le.0))   then
        			call colora("magenta")
				write (s100, 350)   0.001_R8*apl
  350				format ('Toroidal Current   ', f6.0)
				call gtext (s100,40,0)
				xlabval = ccon - 0.08_R8*(xmax-ccon)
		call  setld  (5._R8,5._R8,1,0,1,0)
        			call colora("white")
				write (s100, 351)   1000._R8*times, kcycle
  351			format ('Time = ', f10.2,' ms       cycle ', i7)
								endif
	if(icenter .eq.1 .and. icplgf.gt.0)   	then
						call colora("green")
						write (s100, 352)
  352	format ('   Toroidal g function')
						endif
	if(iright .eq.1)	then
				call colora("red")
				write (s100, 354)   iplim
  354	format ('Poloidal Flux    IPLIM ', i2)
				endif
	call gtext (s100,40,0)
	call colora("white")
!						Border and Y-axis for contour plot
      call boxax (ccon,xmax, zmin,zmax, ixax, jyax)
!							Ticks
	call  ticdraw (ccon,xmax, zmin,zmax, 1)
	else
      call maps(ccon,xmax,zmin,zmax,x142,x858,xco,1._R8)
							 endif
      if(itype.eq.3) then
                     endif
!
!
       amin = 1.E20_R8
       amax =-1.E20_R8
       do 10 i=2,nxp
       do 10 j=2,nzp
       go to(5,5,5,5,4,4,8,7,5,5),itype
    5 vecx(i,j) = a(i,j)
       if(itype.eq.1) vecx(i,j) = -a(i,j)
       if(itype.ne.3) go to 11
       do 13 ii=nc0,nwire
   13 if(i.eq.iwire(ii) .and. j.eq.jwire(ii)) vecx(i,j) = 0._R8
       go to 11
    6 vecx(i,j) = 0._R8
      vecx(i,j) = .25_R8*(a(i,j)/ajey(i)                                 &
     &                 +a(i+1,j)/ajey(i+1)                               &
     &                 +a(i,j+1)/ajey(i)                                 &
     &                 +a(i+1,j+1)/ajey(i+1) )
      if(i.le.2) vecx(i,j)=.5_R8*(a(i+1,j)+a(i+1,j+1))/ajey(i+1)
      go to 11
    8 continue
      fac = 1.0_R8
      if(irfp.eq.1) go to 17
      vecx(i,j)   = .25_R8*(a(i,j)/(ajey(i)*xary(i))                     &
     &           +        a(i+1,j)/(ajey(i+1)*xary(i+1))                 &
     &           +        a(i,j+1)/(ajey(i)*xary(i))                     &
     &           +        a(i+1,j+1)/(ajey(i+1)*xary(i+1)) )
       go to 11
    7 vecx(i,j) = .25_R8*(a(i,j)*xsqoj(i) + a(i+1,j)*xsqoj(i+1)          &
     &            +a(i,j+1)*xsqoj(i)+a(i+1,j+1)*xsqoj(i+1)) - gzero
      if(i.le.2) vecx(i,j)=.5_R8*(a(i+1,j)+a(i+1,j+1))*xsqoj(i+1) -gzero
      go to 11
    4 vecx(i,j) = 0
       fac = (udsh)/(2._R8*udsd)
       if(itype.eq.6) fac = 1._R8
   17 continue
       if(i.eq.2) go to 3
       vecx(i,j) = .25_R8*(a(i,j)+a(i+1,j)+a(i,j+1)+a(i+1,j+1))*fac
      if(itype.eq.5) vecx(i,j) = .25_R8*fac*(                            &
     &   pr(i,j)                                                         &
     & +pr(i+1,j)                                                        &
     & +pr(i,j+1)                                                        &
     & +pr(i+1,j+1) )
      go to 11
    3 vecx(i,j) = .5_R8*fac*(a(i+1,j)+a(i+1,j+1))
      if(itype.eq.5) vecx(i,j) = .5_R8*fac*(                             &
     &   pr(i+1,j)                                                       &
     &+ pr(i+1,j+1))                                                     
      go to 11
    9 vecx(i,j) = sqrt(a(i,j))
   11 continue
      if(acoef(501) .ne. 0  .and. itype.eq.8) then
      if(psi(i,j) .lt. psilim) go to 10
      if(psi(i,j) .gt. phalo) go to 10
      if(iexv(i,j).ne.0) go to 10
      endif
      if(acoef(504).ne.0 .and. itype.eq.2) then
      if(iexv(i,j).ne.0) go to 10
      if(psi(i,j) .gt. psilim) go to 10
      if(psi(i,j) .lt. psimin) go to 10
      endif
      aa = vecx(i,j)
      if(aa.gt.amax) amax=aa
      if(aa.lt.amin) amin=aa
   10 continue
      itmax =1+isym
!
       nxpl = nx
       if(acoef(48).ge.1._R8) nxpl = int(acoef(48))
       nxv = nx
       if(acoef(49).ge.1._R8) nxv = int(acoef(49))
!
!        write (nout,887)
! 887    format (" jrcontr   it   k1   k2    cval1    cval2     amin
!    & amax  psilim   psimin")
       if (imovie.ge.6 .and. imovie.lt.10)   go to 500
!=============================================================================
      if(itype.eq.9) then
      		do 83 n=1,nobs
      		x1 = xobs(n)
      		z1 = zobs(n) + .5_R8*deez
      		x2 = xobs(n) - .5_R8*deex
      		z2 = zobs(n)
      		x3 = xobs(n)
      		z3 = zobs(n) - .5_R8*deez
      		x4 = xobs(n) + .5_R8*deex
      		z4 = zobs(n)
      		call setcrt(x1,z1)
      		call vector(x2,z2)
      		call vector(x3,z3)
      		call vector(x4,z4)
      		call vector(x1,z1)
   83 		continue
                     endif
!
	call  colora("white")
!--------------------------------------
       do 100 it=1,itmax
       sym = 1._R8
       if(it.eq.2) sym = -1._R8
!
      do 50 j=2,nzp
      vzy(j) = sym*zary(j)
   50 continue
!
      if(itype.eq.9) then
		k1 = 1
		k2 = 2
		cval(1) = 0._R8
		imin = 2
		imax = nxp
		jmin = nh-jinc*(1-isym)
		jmax = nh+jinc
      call rcontr(k1,cval,k2,vecx,penx,xary,imin,imax,1,vzy,jmin,jmax,1)
                    endif
!
	if (itype.eq.3)	call colora("magenta")
	if (itype.eq.8)	call colora("green")
      k1 = -nxpl
      k2 = 0
      if(itype.eq.1.and.igone.eq.0) k2 = -k1+1
!..logic for antisymmetric toroidal field plots
      if(jsym.eq.-1.and.itype.eq.3.and.it.eq.2) k2=-k1+1
      if(jsym.eq.-1.and.itype.eq.8.and.it.eq.2) k2=-k1+1
!
      cval(1) = amin
      cval(2) = amax
!
      if(itype.eq.1.and.igone.eq.0) cval(1) = -psilim
      if(itype.eq.1.and.igone.eq.0) cval(2) = -psimin
      imin = 2
      imax = imnpl
      jmin = nh-jinc*(1-isym)
      jmax = nh+jinc
      if(itype.eq.1.and.igone.eq.0) go to 74
      if(imovie.ge.3 .and. imovie.lt.10) call colora("red")
!
! 888   format (4x, 4i5, 6f8.4)
        jrcontr = 1
      call rcontr(k1,cval,k2,vecx,penx,xary,imin,imax,1,vzy,jmin,jmax,1)
      if(itype.ne.1.or.igone.ne.0) go to 75
   74 k1 = -nxv
      k2 = -k1+1
      cval(1) = amin
      cval(2) = -psilim
      if(imovie.ge.3 .and. imovie.lt.10) call colora("cyan")
        jrcontr = 2
!     if (it.eq.1) write (nout,888) jrcontr, it, k1,k2, cval(1),cval(2),
!    &                              amin,amax, psilim, psimin
      call rcontr(k1,cval,k2,vecx,penx,xary,imin,imax,1,vzy,jmin,jmax,1)
!------------------- 30 Apr 92  ----------------------------------
!indranil: replaced the following line with the next one
!        if (itype.ne.1 .or. acoef(97).le.0.0_R8)  go to 75
        if (itype.ne.1 .or. whalos.le.0.0_R8)  go to 75
!                                    Plot one solid contour at psi= phalo
        call colora("yellow")
        k1 = 1
        k2 = 1
        cval(1) = - phalo
      call rcontr(k1,cval,k2,vecx,penx,xary,imin,imax,1,vzy,jmin,jmax,1)
        call colora("red")
!----------------------------------------------------------------
   75 continue
	if (itype.eq.3)	call colora("magenta")
	if (itype.eq.8)	call colora("green")
      k1 = -nxpl
      k2 = 0
!
      cval(1) = amin
      cval(2) = amax
!..logic for antisymmetric toroidal field plots
      if(jsym.eq.-1.and.itype.eq.3.and.it.eq.2) k2=-k1+1
!
      if(itype.eq.1.and.igone.eq.0) cval(1) = -psilim
      if(itype.eq.1.and.igone.eq.0) cval(2) = -psimin
      imin = imnpl
      imax = imxpl
      jmin = nh-jinc*(1-isym)
      jmax = nh+jinc
      if(itype.eq.1.and.igone.eq.0) jmin = jmnpl
      if(itype.eq.1.and.igone.eq.0 .and. isym.eq.1) jmin = 2
      if(itype.eq.1.and.igone.eq.0) jmax = jmxpl
      if(imovie.ge.3 .and. imovie.lt.10) call colora("red")
!
        jrcontr = 3
!     if (it.eq.1) write (nout,888) jrcontr, it, k1,k2, cval(1),cval(2),
!    &                              amin,amax, psilim, psimin
      call rcontr(k1,cval,k2,vecx,penx,xary,imin,imax,1,vzy,jmin,jmax,1)
!------------------- 30 Apr 92  ----------------------------------
!indranil: replaced the following line with the next one
!        if (itype.eq.1 .and. acoef(97).gt.0.0)    then
        if (itype.eq.1 .and. whalos.gt.0.0)    then        
!                                    Plot one solid contour at psi= phalo
        call colora("yellow")
        k1 = 1
        k2 = 1
        cval(1) = - phalo
      call rcontr(k1,cval,k2,vecx,penx,xary,imin,imax,1,vzy,jmin,jmax,1)
        call colora("red")
                           			endif
!----------------------------------------------------------------
!....added 04/15/05
      if(itype.ne.1.or.igone.ne.0 ) go to 85
      k1 = -nxv
      k2 = -k1+1
      cval(1) = amin
      cval(2) = -psilim
      jmin = nh-jinc*(1-isym)
      jmax = nh+jinc
      if(imovie.ge.3 .and. imovie.lt.10) call colora("cyan")
         jrcontr = 4
      call rcontr(k1,cval,k2,vecx,penx,xary,imin,imax,1,vzy,jmin,jmax,1)
!------------------- 24 Mar 94  ----------------------------------
!indranil: replaced the following line with the next one
!        if (itype.eq.1 .and. acoef(97).gt.0.0_R8)    then
        if (itype.eq.1 .and. whalos.gt.0.0)    then
!                                    Plot one solid contour at psi= phalo
        call colora("yellow")
        k1 = 1
        k2 = 1
        cval(1) = - phalo
      call rcontr(k1,cval,k2,vecx,penx,xary,imin,imax,1,vzy,jmin,jmax,1)
        call colora("red")
                           			endif
!----------------------------------------------------------------
   85 continue
	if (itype.eq.3)	call colora("magenta")
	if (itype.eq.8)	call colora("green")
      k1 = -nxpl
      k2 = 0
      if(itype.eq.1.and.igone.eq.0) k2 = -k1+1
!..logic for antisymmetric toroidal field plots
      if(jsym.eq.-1.and.itype.eq.3.and.it.eq.2) k2=-k1+1
!
      cval(1) = amin
      cval(2) = amax
!
      if(itype.eq.1.and.igone.eq.0) cval(1) = -psilim
      if(itype.eq.1.and.igone.eq.0) cval(2) = -psimin
      imin = imxpl
      imax = nxp
      jmin = nh-jinc*(1-isym)
      jmax = nh+jinc
      if(itype.eq.1.and.igone.eq.0) go to 76
      if(imovie.ge.3 .and. imovie.lt.10) call colora("red")
!
      call rcontr(k1,cval,k2,vecx,penx,xary,imin,imax,1,vzy,jmin,jmax,1)
      if(itype.ne.1.or.igone.ne.0) go to 77
   76 k1 = -nxv
      k2 = -k1+1
      cval(1) = amin
      cval(2) = -psilim
      if(imovie.ge.3 .and. imovie.lt.10) call colora("cyan")
      call rcontr(k1,cval,k2,vecx,penx,xary,imin,imax,1,vzy,jmin,jmax,1)
!------------------- 30 Apr 92  ----------------------------------
!indranil: replaced the following line with the next one
!        if (itype.ne.1 .or. acoef(97).le.0.0_R8)  go to 77
        if (itype.ne.1 .or. whalos.le.0.0_R8)  go to 77
!                                    Plot one solid contour at psi= phalo
        call colora("yellow")
        k1 = 1
        k2 = 1
        cval(1) = - phalo
      call rcontr(k1,cval,k2,vecx,penx,xary,imin,imax,1,vzy,jmin,jmax,1)
        call colora("red")
!----------------------------------------------------------------
   77 continue
!
!
      if(itype.eq.1 .and. acoef(501).ne.0 .and. phalo.gt.psilim) then
      imin = iminn
      imax = imaxx
      jmin = jminn
      jmax = jmaxx
      cval(1) = -phalo
      cval(2) = -psilim
      k1 = -nxv
      k2 = -k1+1
      call rcontr(k1,cval,k2,vecx,penx,                                  &
     &     xary,imin,imax,1,vzy,jmin,jmax,1)
            endif
!
!
  100 continue
        go to 600
!=============================================================================
!
!*****************************************************************************
!                                       IMOVIE = 6
  500   continue
!
       do 590 it=1,itmax
       sym = 1._R8
       if(it.eq.2) sym = -1._R8
!
       do 530 j=2,nzp
       vzy(j) = sym*zary(j)
  530 continue
      imin = 2
      imax = imnpl
      jmin = nh-jinc*(1-isym)
      jmax = nh+jinc
      if(itype.eq.1.and.igone.eq.0) go to 574
      call colora("red")
!
      call rcontr(k1,cval,k2,vecx,penx,xary,imin,imax,1,vzy,jmin,jmax,1)
      if(itype.ne.1.or.igone.ne.0) go to 575
  574   call colora("cyan")
        jrcontr = 2
      k2sav = k2
      k1 = k2
!     if (it.eq.1) write (nout,888) jrcontr, it, k1,k2, cval(1),cval(2),
!    &                              amin,amax, psilim, psimin
      call rcontr(k1,cval,k2,vecx,penx,xary,imin,imax,1,vzy,jmin,jmax,1)
      k2 = k2sav
      k1 = k1sav
  575 continue
      imin = imnpl
      imax = imxpl
      jmin = nh-jinc*(1-isym)
      jmax = nh+jinc
      if(itype.eq.1.and.igone.eq.0) jmin = jmnpl
      if(itype.eq.1.and.igone.eq.0 .and. isym.eq.1) jmin = 2
      if(itype.eq.1.and.igone.eq.0) jmax = jmxpl
      call colora("red")
!
         jrcontr = 3
!     if (it.eq.1) write (nout,888) jrcontr, it, k1,k2, cval(1),cval(2),
!    &                              amin,amax, psilim, psimin
      k1sav = k1
      k2sav = k2
      k1 = nconlev-k2+1
      k2 = 1
      call rcontr(k1,cval(k2sav),k2,vecx,penx,xary,imin,imax,1,vzy,jmin  &
     &    ,jmax,1)
      k1 = k1sav
      k2 = k2sav
!
!                               Plot one solid contour at psi= psilim
        call colora("yellow")
        k1sav = k1
        k2sav = k2
        cval1sav = cval(1)
        k1 = 1
        k2 = 1
        cval(1) = - psilim
      call rcontr(k1,cval,k2,vecx,penx,xary,imin,imax,1,vzy,jmin,jmax,1)
        k1 = k1sav
        k2 = k2sav
        cval(1) = cval1sav
      if(itype.ne.1.or.igone.ne.0) go to 585
      jmin = nh-jinc*(1-isym)
      jmax = nh+jinc
      call colora("cyan")
         jrcontr = 4
        k2sav = k2
        k2 = nconlev + 1
      k1sav = k1
      k1 = k2sav-1
      call rcontr(k1,cval,k2,vecx,penx,xary,imin,imax,1,vzy,jmin,jmax,1)
        k2 = k2sav
      k1 = k1sav
  585 continue
      imin = imxpl
      imax = nxp
      jmin = nh-jinc*(1-isym)
      jmax = nh+jinc
      if(itype.eq.1.and.igone.eq.0) go to 576
      call colora("red")
!
      call rcontr(k1,cval,k2,vecx,penx,xary,imin,imax,1,vzy,jmin,jmax,1)
      if(itype.ne.1.or.igone.ne.0) go to 577
  576 call colora("cyan")
      k1sav = k1
      k2sav = k2
      k1  = k2sav - 1
      call rcontr(k1,cval,k2,vecx,penx,xary,imin,imax,1,vzy,jmin,jmax,1)
      k2 = k2sav
      k1 = k1sav
  577 continue
  590 continue
!*****************************************************************************
!
  600  	call colora("white")
	if(itype.ne.1) go to 99
      do 80 i=2,nxp
      do 80 j=2,nzp
   80 vecx(i,j) = -vecx(i,j)
   99 continue
      if(itype.gt.3) go to 8100
	call colora("white")
      do 81 ii=nc0,nwire
      if(acoef(1).eq.1._R8.and.igroupw(ii).eq.8) go to 81
      call coildr(iwire(ii),jwire(ii),ilwire(ii))
   81 continue
	if (itype.eq.1)		then
!							green limiters
	call colora("green")
      do 82 n=1,nlim
      x1 = xlima(n) + .25_R8*deex
      z1 = zlima(n) + .25_R8*deez
      x2 = xlima(n) - .25_R8*deex
      z2 = zlima(n) - .25_R8*deez
      x3 = xlima(n) + .25_R8*deex
      z3 = zlima(n) - .25_R8*deez
      x4 = xlima(n) - .25_R8*deex
      z4 = zlima(n) + .25_R8*deez
      call setcrt(x1,z1)
      call vector(x2,z2)
      call setcrt(x3,z3)
      call vector(x4,z4)
   82 continue
      if(idata.eq.9 .or. idata.eq.11) then
      do 182 n=1,nlimnstx
      x1 = xlimnstx(n) + .25_R8*deex
      z1 = zlimnstx(n) + .25_R8*deez
      x2 = xlimnstx(n) - .25_R8*deex
      z2 = zlimnstx(n) - .25_R8*deez
      x3 = xlimnstx(n) + .25_R8*deex
      z3 = zlimnstx(n) - .25_R8*deez
      x4 = xlimnstx(n) - .25_R8*deex
      z4 = zlimnstx(n) + .25_R8*deez
      call setcrt(x1,z1)
      call vector(x2,z2)
      call setcrt(x3,z3)
      call vector(x4,z4)
      call setlch(xlimnstx(n),zlimnstx(n),0,1,0,-1)
      call gtext(nstxchar(n),1,0)
  182 continue
                     endif
				endif
      if(iplate.eq.0 .or. nplate.eq.0) go to 8100
      do 8110 n=1,nplate
      call setcrt(xsega(n,1),zsega(n,1))
      do 8111 l=2,nseg(n)+1
 8111 call vector(xsega(n,l),zsega(n,l))
      if(isym.eq.0) go to 8110
      z1neg = -zsega(n,1)
      call setcrt(xsega(n,1),z1neg)
      do 8112 l=2,nseg(n)+1
      zlneg = -zsega(n,l)
 8112 call vector(xsega(n,l),zlneg)
 8110 continue
 8100 continue
!.....special addition to G plot for CHI simulations
      if(itype.eq.8 .and. ilower*iupper*jlower*jupper .gt.0) then
      x1 = xarh(iupper+1)
      z1 = .5_R8*(zary(jupper-1)+zary(jupper-2))
      x2 = xarh(ilower+5)
      z2 = .5_R8*(zary(jlower+2)+zary(jlower+1))
      x3 = xarh(ilower+1)
      z3 = .5_R8*(zary(jlower+1)+zary(jlower  ))
      call setcrt(x1,z1)
      call vector(x2,z2)
      call setcrt(x1,z1)
      call vector(x3,z3)
      endif
!
      if((imovie.gt.0 .and. imovie.lt.10)) go to 400
      if(acoef(901).eq.0._R8) then
      asetld = 1._R8
!     if(iright.eq.1) asetld = 24.
      if(iright.eq.1) asetld = 3.0_R8
      call setld(asetld,40._R8,1,0,1,1)
      if(iright.ne.1 .and. ileft.ne.1 .and. icenter.ne.1) then
      write(s100,1000)amin,amax
      call gtext(s100,80,0)
      endif
                                       endif
 1000 format(1x," min =",1pe12.4,"  max =",1pe12.4)
      if(irfp.eq.1) go to 200
	if (ilcrtest.ne.0 .and. itype.ne.3)   go to 110
      go to(101,102,103,104,105,106,107,108,109,110),itype
!
!
!c
  101 continue
      if(acoef(901).gt.0._R8) then
      call setld(1._R8,11._R8,1,0,1,0)
      plcprnt = tcurdtp*tpi*udsi
      write(s100,2005) gzero/xplas,plcprnt, fluxlink
 2005 format(" Bt[T]=",1pe10.3,"  Ip[A]  =",1pe10.3,                     &
     &   " Flux Linkage[W]=",1pe10.3)
      call gtext(s100,80,0)
!
      alamda = ali2+betapol
      write(s100,2167) 2._R8*ali2,betapol,beta*100._R8
 2167 format(" Li(3)=",1pe10.3,"  Betapol=",1pe10.3,                     &
     &   " Beta[%]        =",1pe10.3)
      call gtext(s100,80,0)
!
      write(s100,2168) rmajor,rminor
 2168 format(" R0[m]=   ",f7.3,"  a[m]   =   ",f7.3)
      call gtext(s100,80,0)
!
      write(s100,3168) shape5,shape3
 3168 format(" k    =   ",f7.3,"  d      =   ",f7.3)
      call gtext(s100,80,0)
!
      write(s100,21680) el95,del95
21680 format(" k95  =   ",f7.3,"  d95    =   ",f7.3)
      call gtext(s100,80,0)
!
      write(s100,2169) shape6,shape7
 2169 format(" Rx[m]=   ",f7.3,"  Zx[m]  =   ",f7.3)
      call gtext(s100,80,0)
!
      write(s100,2170) qprof2(1),global(15)
 2170 format(" q0   =   ",f7.3,"  q100   =   ",f7.3)
      call gtext(s100,80,0)
      write(s100,3170) global(61),global(62)
 3170 format(" q95  =   ",f7.3,"  qstar  =   ",f7.3)
      call gtext(s100,80,0)
!      write(s100,3171) findex
! 3171 format(" nindx=   ",f7.3)
!      call gtext(s100,80,0)
      write(nsc1,1019)
 1019 format(" poloidal flux")
!
      goto 110
      endif
      write(s100,1001) times,kcycle
      call gtext(s100,80,0)
 1001 format(1x," poloidal flux , time=",1pe10.3," cyc=",i7)
      aiprnt = tcurdtp*tpi*udsi
      write(s100,2001) aiprnt,iplim
      call gtext(s100,80,0)
 2001 format(" plasma current ",1pe12.4,"   iplim",i3)
      if(iright.ne.1) write(nsc1,1011) kcycle
 1011 format(1x," poloidal flux,  cyc=",i7)
      go to 110
  102 write(s100,1002) times,kcycle
      call gtext(s100,80,0)
 1002 format(1x,"res, time=",1pe10.3," cyc=",i7)
      write(nsc1,1012) kcycle
 1012 format(1x," res, cyc=",i7)
!     write(99,9001) kcycle
!     do j=2,nzp
!        write(99,9002) (vecx(i,j),i=2,nxp)
!     enddo
!9001 format("  resistivity,  kcycle =",i6)
!9002 format(1p8e11.3)
      go to 110
  103 if (ilcrtest.eq.0)	then
				write(s100,1003) times,kcycle
				call gtext(s100,80,0)
 1003 format(1x," tor current , time=",1pe10.3," cyc=",i7)
      write(nsc1,1013) kcycle
 1013 format(" tor current, cyc=",i7)
      go to 110
				endif
	if(ileft.eq.1)	write(nsc1,2013) kcycle
 2013			format(" current, g, flux,cyc=",i7)
	if(icenter.eq.1 .and. icplgf.eq.0)	write(nsc1,3013) kcycle
 3013	format(" current and flux,cyc=",i7)
	go to 110
  104 write(s100,1004) times,kcycle
      call gtext(s100,80,0)
 1004 format(1x," b..delstar(a) , time=",1pe10.3," cyc=",i7)
      write(nsc1,1014) kcycle
 1014 format(1x," b..delstar(a)  cyc=",i7)
      go to 110
  105 write(s100,1005) times,kcycle,beta
      call gtext(s100,80,0)
 1005 format(1x,"(ev) time=",1pe10.3," cyc=",i7,                         &
     &      " beta",1pe10.3)
      write(nsc1,1015) kcycle
 1015 format(" pressure,   cyc=",i7)
      go to 110
  106 write(s100,1006) times,kcycle
      call gtext(s100,80,0)
 1006 format(1x," u..del sq(omega) , time=",1pe10.3," cyc=",i7)
      write(nsc1,1016) kcycle
 1016 format(1x," u..del sq(omega),cyc=",i7)
      go to 110
  107 write(s100,1007) times,kcycle
      call gtext(s100,80,0)
 1007 format(1x," toroidal velocity , time=",1pe10.3," cyc=",i7)
      write(nsc1,1017) kcycle
 1017 format(1x," toroidal velocity, cyc=",i7)
      go to 110
  108 write(s100,1008) times,kcycle
      call gtext(s100,80,0)
 1008 format(1x," tor g funct, time=",1pe10.3," cyc=",i7)
      write(nsc1,1018) kcycle
 1018 format(1x," tor g funct, cyc=",i7)
      go to 110
  109 write(s100,1009) times,kcycle
      call gtext(s100,80,0)
 1009 format (1x," asymmetric flux  time=",1pe12.5," cyc=",i7)
      write(nsc1,10190) kcycle
10190 format(1x," asymmetric flux, cyc=",i7)
      go to 110
  200 continue
      go to(101,202,103,204,105,206,207,108,109,202),itype
  202 write(s100,2002) times,kcycle
      call gtext(s100,80,0)
 2002 format(1x,"eta.J, time=",1pe10.3," cyc=",i7)
      write(nsc1,2012) kcycle
 2012 format(1x," eta.J, cyc=",i7)
      go to 110
  204 write(s100,2004) times,kcycle
      call gtext(s100,80,0)
 2004 format(1x," Hyper/J , time=",1pe10.3," cyc=",i7)
      write(nsc1,2014) kcycle
 2014 format(1x," Hyper/J  cyc=",i7)
      go to 110
  206 write(s100,2006) times,kcycle
      call gtext(s100,80,0)
 2006 format(1x," (eta.J+Hyper)/eta.J , time=",1pe10.3," cyc=",i7)
      write(nsc1,2016) kcycle
 2016 format(1x," (eta.J+Hyper)/eta.J,cyc=",i7)
      go to 110
  207 write(s100,2007) times,kcycle
      call gtext(s100,80,0)
 2007 format(1x," J.B/B**2 , time=",1pe10.3," cyc=",i7)
      write(nsc1,2017) kcycle
 2017 format(1x," J.B/B**2, cyc=",i7)
      go to 110
  110 continue
      ipassw = itype
      if(acoef(901).eq.0._R8)                                            &
     &call wplot(ipassw,x142,xco,x858,iright,ileft,icenter)
      if(acoef(901).gt.0._R8) call frscj(8)
      return
!
  400 continue
      if(imovie.ge.3 .and. imovie.lt.10) call colora("white")
      call setold(xmax,zmax,1,0,1,0)
      aloop = reboun*udsv*tpi
      psiaxis = psimin*tpi
      write(s100,1010) times,kcycle,betat0, xmag, zmag, psiaxis
      call gtextm(s100,80,0,1,18)
 1010 format(/,4x,"time",/,f8.5,//,3x,"cycle",/,i8,//,4x,"beta",/,       &
     & f8.3,//,4x,"xmag",/,f8.3,//,4x,"zmag",/,f8.3,//,4x,"psim",/,      &
     &   f8.3)
      if(imovie.eq.2 .or. imovie.eq.4) go to 300
      xmin = 0
      xmax = tpros(ntpts)
      ymin = 0._R8
      ymax = 0
      do ii=1,ntpts
      ymax = max(ymax,pcur(ii)*1.05_R8/1.E6_R8)
      enddo
!
!*************       xco = 1.-.716*2./(alx-acoef(3))
      xco = xco - 0.18_R8*(xco - 0.11328_R8)
      if(imovie.ge.3 .and. imovie.lt.10) call colora("white")
      call maps(xmin,xmax,ymin,ymax,x142,x858,.11328_R8,xco)
      call setold(xmin - 0.18_R8*(xmax-xmin),ymin,1,0,1,1)
      write(s100,9020)
 9020 format(" Plasma Current (MA) ")
      call gtextm(s100,80,0,1,1)
!DIAG
      write(nterm,8888) xmin,xmax,ymin,ymax,x142,x858,xco,nframe
 8888 format( 7e12.4,i5)
      nframe = nframe+1
      xnf(nframe) = times
      ynf(nframe,1) = apl/1.E6_R8
      ynf(nframe,2) = paux/1.E7_R8
      ynf(nframe,3) = palpha/1.E7_R8
!***c***c***do 210 ii=1,3,2
!***c***c***ic = ii+1
!***c***c***if(imovie.ge.3 .and. imovie.lt.10) call color(ic)
!***210 call tracec(nchar(ii),xnf(1),ynf(1,ii),nframe)
      if(imovie.ge.3 .and. imovie.lt.10) call colora("red")
      call trace(xnf(1), ynf(1,1), nframe,-1,-1,0._R8,0._R8)
      if(istart.gt.1) call trace(tpros,pcur,istart,-1,-1,0._R8,0._R8)
      if(imovie.ge.3 .and. imovie.lt.10) call colora("yellow")
      call trace (xnf(1), ynf(1,2), nframe,-1,-1,0._R8,0._R8)
      if(imovie.ge.3 .and. imovie.lt.10) call colora("white")
      call trace (xnf(1), ynf(1,3), nframe,-1,-1,0._R8,0._R8)
      call frscj(8)
      return
  300 continue
      call frscj(8)
      return
      end
!**********************************************************************
!......6.92 wplot
! 21Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

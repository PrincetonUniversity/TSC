      Program test
c
c     implicit real*8 (a-h,o-z)
c     include 'comfile.i'
      parameter(nloopsp=100)
      parameter(ncoilsp=100)
      dimension xplot(10000),yplot(10000),zplot(10000),
     1 eplot(10000), fplot(10000),gplot(10000),hplot(1000)
      dimension pflux(nloopsp),eflux(nloopsp),
     1          error(nloopsp),icolor(nloopsp),
     2          flux_r(nloopsp),flux_z(nloopsp),
     3          iflux_active(nloopsp)
      dimension pfluxplt(10000,nloopsp),
     1          efluxplt(10000,nloopsp)
      dimension csuma(ncoilsp),ppcur(ncoilsp),fbcur(ncoilsp),div(2)
      character*4 istring4
      character*14 label,istring14
      character*2 inum(75)
      data (inum(i),i=1,75)/
     1 " 1"," 2"," 3"," 4"," 5"," 6"," 7"," 8"," 9",
     1 "10","11","12","13","14","15","16","17","18","19",
     1 "20","21","22","23","24","25","26","27","28","29",
     1 "30","31","32","33","34","35","36","37","38","39",
     1 "40","41","42","43","44","45","46","47","48","49",
     1 "50","51","52","53","54","55","56","57","58","59",
     1 "60","61","62","63","64","65","66","67","68","69",
     1 "70","71","72","73","74","75"/
c
c
c
c....1.0  initialize constants and open files
      pi = 3.1415926535
c
      no167a = 9
      open(7,file='plotout',status='unknown',iostat=ios)
      open(no167a,file='tscnstxdata_out',status='old',iostat=ios22)
      call ncarcgm(1,"out.cgm")
      call dders(-1)
c
c.....initialize NCAR GKS graphics
      iwkid = 2
      iwtype = 8
      lunit = 2
      ierrf = 6
c     call gopks(ierrf,1000)
      call gopwk(iwkid,lunit,iwtype)
      call gacwk(iwkid)
c
c.....initialize program by ploting flux loop positions
      call readmag(icolor,flux_r,flux_z,iflux_active,nloopsp,nc,nf)
c
c.....first pass through data to determine error catagories
      do i=1,nf
      error(i) = 0.
      enddo
c
      do irec=1,10000
      read(no167a,7010,end=102) times,plcurs,vespplass,vescurs,plcure,
     1                          vescure, vespplase,chicur
      read(no167a,7010) (csuma(k),k=1,nc-2)
      read(no167a,7010) (ppcur(k),k=1,nc-2)
      read(no167a,7010) (fbcur(k),k=1,nc-2)
      if(nf.gt.0) read(no167a,7010) (pflux(i),i=1,nf)
      if(nf.gt.0) read(no167a,7010) (eflux(i),i=1,nf)
      read(no167a,7011)
      write(99,7010) (pflux(i),i=1,5),(eflux(i),i=1,5)
      do i=1,nf
      error(i) = error(i) + (pflux(i)-eflux(i))**2
      enddo
c
      sumerror = 0.
      do i=1,nf
      sumerror = sumerror + abs(pflux(i)-eflux(i))
      sumerror2= sumerror2+ (pflux(i)-eflux(i))**2
      enddo
      sumerror = sumerror/nf
      sumerror2= sqrt(sumerror2/nf)
      write(7,1006) irec,times,plcurs,sumerror,sumerror2
 1006 format(i5,1p4e12.4)
      enddo
 102  continue
      do i=1,nf
      error(i) = sqrt(error(i)/irec)
      icolor(i) = 2
      if(error(i).gt. 0.090) icolor(i) = 3
      if(error(i).lt. 0.025) icolor(i) = 1
      enddo
      write(6,1005) ((i,error(i)),i=1,nf)
 1005 format(3(i5,1pe10.2))
c    
 1000 continue
      write(6,1002)
      read(5,*)nn
      call frame(0)
c
      if(nn.lt.0) go to 2000
c

c
      if(nn.gt.100 .and. nn.lt.200) then
      label = "current no"
      label(11:12) = inum(nn-100)
      endif        
c
c
      if(nn.gt.0 .and. nn.le.nf) then
      label = "flux loop "
      label(11:12) = inum(nn)
      endif
c
      if(nn.eq.0) then
      label = "vessel current"
      endif
      if(nn.eq.201) then
      label = "CHI current"
      endif
c
      write(6,1010) nn,label
 1010 format("nn = ",i3,3x,a14)
c
      rewind(no167a)
c
      ymin = 0.
      ymax = 0.
      do irec=1,10000
      read(no167a,7010,end=101) times,plcurs,vespplass,vescurs,plcure,
     1       vescure, vespplase, chicure,vescure2,chicurs
      read(no167a,7010) (csuma(k),k=1,nc-2)
      read(no167a,7010) (ppcur(k),k=1,nc-2)
      read(no167a,7010) (fbcur(k),k=1,nc-2)
      if(nf.gt.0) read(no167a,7010) (pflux(i),i=1,nf)
      if(nf.gt.0) read(no167a,7010) (eflux(i),i=1,nf)
      read(no167a,7011)
c    
c
      do i=1,nf
      pfluxplt(irec,i) = pflux(i)
      efluxplt(irec,i) = eflux(i)
      enddo
c
      xplot(irec) = times

      if(nn.eq.200) then
      yplot(irec) = csuma(1)/100000.
      zplot(irec) = csuma(2)/1000.
      eplot(irec) = csuma(3)/1000.
      fplot(irec) = csuma(5)/1000.
                    endif
c
      if(nn.gt.100 .and. nn.lt.200) then
      yplot(irec) = csuma(nn-100)
      zplot(irec) = ppcur(nn-100)
      fplot(irec) = fbcur(nn-100)
      endif
c
c
c
      if(nn.gt.0 .and. nn.le.nf) then
      yplot(irec) = pflux(nn)
      zplot(irec) = eflux(nn)
      endif
      if(nn.eq.201) then
      yplot(irec) = chicurs
      zplot(irec) = chicure
      endif
c
      if(nn.eq.0) then
c     yplot(irec) = vespplass
c     zplot(irec) = vespplase - ppcur(7)
c     zplot(irec) = vespplase
      eplot(irec) = plcure
      fplot(irec) = plcurs
      gplot(irec) = vescure
      hplot(irec) = vescurs
      endif
c
      ymin = amin1(ymin,yplot(irec),zplot(irec))
      ymax = amax1(ymax,yplot(irec),zplot(irec))
      if(nn.eq.0) then
      ymin = amin1(ymin,eplot(irec),fplot(irec))
      ymax = amax1(ymax,eplot(irec),fplot(irec))
      ymin = amin1(ymin,gplot(irec),hplot(irec))
      ymax = amax1(ymax,gplot(irec),hplot(irec))
      endif
      if(nn.eq.200) then
c
      ymin = -380.
      ymax =  180.
      endif
c
      imax = irec
      enddo
      close(no167a)
 101  continue
      write(6,3333) imax,ymin,ymax
 3333 format(" imax,ymin,ymax =",i5,1p2e12.4)
      if(nn.eq.300) go to 300
c
      bmax = .89
      bmin = .20
      if(nn.eq.200) then
      bmax = .49
      bmin = .20
      endif
      call colora("blue")
      call maps(xplot(1),xplot(imax),ymin,ymax,.20,.99,bmin,bmax)
      if(nn.eq.0) go to 198
      if(nn.eq.200) go to 200
      call tracec(1hS,xplot,yplot,imax,-1,-1,0.,0.)
      call colora("red")
      call tracec(1hE,xplot,zplot,imax,-1,-1,0.,0.)
      call colora("blue")
      call tracec(1hS,xplot,yplot,imax,-1,-1,0.,0.)
c
c      if(nn.gt.100.and. nn.lt.200) then
c      call tracec(1hF,xplot,fplot,imax,-1,-1,0.,0.)
c      call tracec(1hF,xplot,fplot,imax,-1,-1,0.,0.)
c      endif
c
  198 continue
      if(nn.eq.0) then
      call colora("red")
      call tracec(1hE,xplot,eplot,imax,-1,-1,0.,0.)
      call tracec(1hS,xplot,fplot,imax,-1,-1,0.,0.)
      call colora("green")
      call tracec(1hS,xplot,hplot,imax,-1,-1,0.,0.)
      call tracec(1hE,xplot,gplot,imax,-1,-1,0.,0.)
      call tracec(1hE,xplot,gplot,imax,-1,-1,0.,0.)
      call colora("blue")
      endif
c
c     bottom label
      x = xplot(1) + 0.3*(xplot(imax)-xplot(1))
      y = ymin     - .15*(ymax       -ymin    )
      icase = 1
      isize = 3
      iorient = 0
      istring4 = "time"
      nchar = 4
      call setlch(x,y,icase,isize,iorient,-1)
      call gtext(istring4,nchar,0)
c
c     left label
      x = xplot(1) - .15*(xplot(imax)-xplot(1))
      y = ymin     + .20*(ymax       -ymin    )
      icase = 1
      isize = 3
      iorient = 1
      istring14 = label
      nchar = 14
      call setlch(x,y,icase,isize,iorient,-1)
      call gtext(istring14,nchar,0)
      go to 201
 200  continue
c
      call trace(xplot,yplot,imax,-1,-1,0.,0.)
      call colora("red")
      call trace(xplot,zplot,imax,-1,-1,0.,0.)
      call colora("green")
      call trace(xplot,eplot,imax,-1,-1,0.,0.)
      call colora("blue")
      call trace(xplot,fplot,imax,-1,-1,0.,0.)
      call trace(xplot,fplot,imax,-1,-1,0.,0.)
      call colora("blue")
      go to 201
 300  continue
c
      rmax = 0.
      zmax = 0.
      rmin = 0.
      zmin = 0.
      do i=1,nf
      rmax = amax1(rmax,1.3*flux_r(i))
      rmin = amin1(rmin,    flux_r(i))
      zmax = amax1(zmax,    flux_z(i))
      zmin = amin1(zmin,1.4*flux_z(i))
      enddo
      rdis = 0.3*(rmax-rmin)
      zdis = 0.2*(zmax-zmin)
c

      fluxmin = 0.
      fluxmax = 0.
      do 303 j=1,nf
      if(iflux_active(j).eq.0) go to 303
      do 302 i=1,imax
      fluxmin = amin1(fluxmin,pfluxplt(i,j),efluxplt(i,j))
      fluxmax = amax1(fluxmax,pfluxplt(i,j),efluxplt(i,j))
 302  continue
 303  continue
c
      do 301 i=2,nf

      r1norm = (flux_r(i)     -rmin)/(rmax-rmin)
      r2norm = (flux_r(i)+rdis-rmin)/(rmax-rmin)
      z1norm = (flux_z(i)-zdis-zmin)/(zmax-zmin)
      z2norm = (flux_z(i)     -zmin)/(zmax-zmin)
c 
      if(iflux_active(i).eq.0) go to 301
      call map(xplot(1),xplot(imax),fluxmin,fluxmax,
     1          r1norm,r2norm,z1norm,z2norm)
      call colora("red")
      call trace(xplot,pfluxplt(1,i),imax,-1,-1,0.,0.)
      call colora("green")
      call trace(xplot,efluxplt(1,i),imax,-1,-1,0.,0.)
      call colora("blue")
      call setlch(xplot(1),fluxmax,0,1,0,-1)
      call gtext(inum(i),2,0)
      write(6,1001) i,flux_r(i),flux_z(i),error(i)
 1001 format(i3,1p4e12.4)
c
 301  continue
      call colora("blue")
      call maps(xplot(1),xplot(imax),fluxmin,fluxmax,
     1         .7,1.,.8,1.)
 201  continue
c  
c
c.....empty buffers
      call guwk(2,0)

c
      go to 1000
 2000 continue
      call gclrwk(iwkid,1)
      call gdawk(iwkid)
      call gclwk(iwkid)
      call plote

      stop
c
 7010 format(1p10e12.4)
 7011 format(1x) 
 1002 format(/,"type -1:quit",
     .       /,"      0: I(vv)",
     .       /,"      #: flux loop #",
     .       /,"  100+n: cur n",
     1       /,"    200:special plot of currents",
     2       /,"    201:chi currents",
     3       /,"    300:summary of flux loops" )
      end


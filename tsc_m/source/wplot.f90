      subroutine wplot(ipassw,x142,xco,x858,iright,ileft,icenter)
!......6.92 wplot
!
!.... plots values across the mid line
!
!
      USE CLINAM
      USE SCRATCH
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER iright,ileft,icenter,ipassw,ilcrtest,j,nz2,i,kyside,ii
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 x142,xco,x858,div,xmin,xmax,ymin,ymax,xm1,xm2,xlaboff
      REAL*8 frac
!============
      dimension div(2)
      character*4  iyaxfmt
!============      
!
      call colora("white")
      ilcrtest = ileft + icenter + iright
!
      if(noplot(7).gt.0) then
      xmin=+1.E6_R8
      xmax=-1.E6_R8
      do 602 j=1,nz
      w11(j)=vecx(imag,nzp+1-j)
      w11(2*nz+1-j)=w11(j)
      if(ipassw.eq.3.or.ipassw.eq.8) w11(2*nz+1-j)=jsym*w11(j)
      xmin=min(xmin,w11(j),w11(2*nz+1-j))
      xmax=max(xmax,w11(j),w11(2*nz+1-j))
      w12(j)=zary(nzp+1-j)
      w12(2*nz+1-j)=-w12(j)
      ymin=min(ymin,w12(j),w12(2*nz+1-j))
      ymax=max(ymax,w12(j),w12(2*nz+1-j))
602   continue
      xm1=x142+1.1_R8*(x858-x142)
      xm2=x858+1.1_R8*(x858-x142)
      xm1=xm1-1.E-6_R8
      xm2=xm2-1.E-6_R8
      call map(xmin,xmax,ymin,ymax,xm1,xm2,xco,1.0_R8)
      div(1)=xmin
      div(2)=xmax
      call gaxisf(xmin,ymin,xmax,ymin,1,1,0,'e9.2',4,div)
      ymin=ymin+1.E-6_R8*(ymax-ymin)
      xmax=xmax-1.E-6_R8*(xmax-xmin)
      call setcrt(xmin,ymin)
      call vector(xmax,ymin)
      call vector(xmax,ymax)
      call vector(xmin,ymax)
      call vector(xmin,ymin)
      nz2=2*nz
      call tracec(1hx,w11(1),w12(1),nz2,-1,-1,0._R8,0._R8)
      if(iright.ne.1) call frscj(6)
      return
      endif
!
!
!.....modified 04/27/02
      if(ipassw.eq.2) go to 30
      ymin=1.E6_R8
      ymax=-1.E6_R8
      do 601 i=1,nx
      w11(i) = vecx(i+1,nh)
      if((ipassw.eq.7.or.ipassw.eq.4.or.ipassw.eq.9).and.irfp.eq.0)      &  
     &w11(i) = (vecx(i+1,nh+1)-vecx(i+1,nh))/deez
      if(ipassw.eq.2) w11(i) = max(1.E-15_R8,vecx(i+1,nh))
      ymin = min(ymin,w11(i))
      ymax = max(ymax,w11(i))
      w12(i) = xary(i+1)
  601 continue
      if(ipassw.eq.1) then
!
!.....write midplane values of poloidal flux
!
      write(nout,1600) times, kcycle
 1600 format(" poloidal flux, time =",1pe12.4, "  cycle=",               &  
     &        i7)
      write(nout,1601)
 1601 format("  i      R          psi      ",                            &  
     &   "     i      R          psi      ",                             &  
     &   "     i      R          psi      ")
 1602 format(i3,1p2e12.4,5x,0pi3,1p2e12.4,5x,0pi3,1p2e12.4)
      do 604 i = 1,nx,3
      if(i.le.nx-2) then
      write(nout,1602) i,w12(i),w11(i),i+1,w12(i+1),w11(i+1),            &  
     &                                 i+2,w12(i+2),w11(i+2)
      go to 604
                 endif
      if(i.le.nx-1) then
      write(nout,1602) i,w12(i),w11(i),i+1,w12(i+1),w11(i+1)
      go to 604
                 endif
      write(nout,1602) i,w12(i),w11(i)
  604 continue
                      endif
!
      if(ymax.le.ymin) go to 380
      xmin = ccon
      xmax = w12(nx)
      if(ipassw.eq.2) go to 30
      call map (xmin,xmax,ymin,ymax,x142,x858, .200_R8, xco-.050_R8)
!                          yaxis
      div(1) = ymin
      div(2) = ymax
      iyaxfmt = 'e9.2'
!     if (ymin.gt.-0.0049_R8)   iyaxfmt = 'f4.2'
!     if (abs(ymax).ge.10.0_R8.or. abs(ymin).ge.10.0_R8)   iyaxfmt =     &
!    & 'e9.2'
      kyside = 0
      xlaboff = -0.1_R8
      if (icplgf.eq.0)   then
               kyside = 1
               xlaboff = -0.3_R8
               endif
      call gaxisf(xmin,ymin,xmin,ymax, 0,1,kyside, iyaxfmt, 4,div)
!
!                       Box + x-axis
      call boxax (xmin,xmax, ymin,ymax, 1,0)
!                       Ticks for xaxis
      call  ticdraw (xmin,xmax, ymin,ymax, -1)
      if(ipassw.eq.1)   call colora("red")
      if(ipassw.eq.3)   call colora("magenta")
      if(ipassw.eq.8)   call colora("green")
      call trace (w12(1),w11(1), nx, -1,-1,0._R8,0._R8)
      frac = .15_R8
      if(iright.eq.1 .or. ileft.eq.1) frac = .4375_R8
      call setold(xmin-frac*(xmax-xmin),ymin,1,0,1,1)
      if((ipassw.eq.7.or.ipassw.eq.4.or.ipassw.eq.9).and.irfp.eq.0)      &  
     &   go to 20
      if(ileft .eq.1 .or. (icenter.eq.1 .and. icplgf.le.0))   then
         call  setlch (xmin+xlaboff*(xmax-xmin), ymin, 0,1,1, -1)
               write (s100, 350)
  350             format ('z = zmag')
               call gtext (s100,40,0)
               endif
      go to 380
   20 continue
      call crtbcd('deriv   ' )
      go to 380
!
!.....special for resistivity plots
   30 continue
      ymin=1.E6_R8
      ymax=-1.E6_R8
      ii = 0
      do 603 i=2,nxp
!     if(iexv(i,nh) .ne. 0) go to 603
      ii = ii + 1
      w11(ii) = max(1.E-15_R8,vecx(i,nh))
      ymin = min(ymin,w11(ii))
      ymax = max(ymax,w11(ii))
      w12(ii) = xary(i)
  603 continue
      if(ymax.le.ymin) go to 380
      xmin = ccon
      xmax = w12(nx)
      call mapgsl(xmin,xmax,ymin,ymax,x142,x858,.200_R8,xco)
      call tracec(1hx,w12(1),w11(1),ii,-1,-1,0._R8,0._R8)
      write(nout,5556) kcycle, nh
 5556 format(" Resistivity plot at cycle, nh=",i7,i3)
  380    if (ilcrtest.eq.0)   then
               call  frscj(6)
      else
      if (icplgf.gt.0 .and. icenter.eq.1) call frscj(6)
      if (icplgf.le.0 .and. iright .ne.1) call frscj(6)
               endif
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

      subroutine frscj(iptype)
!.....6.96 frscj
!
!.....produces a graphics index at end of output file
!
      USE CLINAM
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER iptype,nframe,k,ipmes,n,nhalf,np
!============
! idecl:  explicitize implicit REAL declarations:
!     REAL*8 pltindx
!============
      character*8 ititle(2000,4)
      character*32 mes(18)
      data nframe/0/
!
        data (mes(k),k=1,17) /                                           &  
     & "  timing information            ",                               &  
     & "  input data                    ",                               &  
     & "  grid,coils, and limiters      ",                               &  
     & " filament growth rate model     ",                               &  
     & " initial coil and wire info     ",                               &  
     & "                                ",                               &  
     & " coil currents ,cycle=          ",                               &  
     & " movie frame    cycle=          ",                               &  
     & "  summary plot                  ",                               &  
     & " flux contours in plasma        ",                               &  
     & " switch and time step info      ",                               &  
     & " wires, CIT plasma-VV boundary  ",                               &  
     & " density diffusion and sources  ",                               &  
     & "  alpha and total pressure      ",                               &  
     & "                                ",                               &  
     & " plasma impurity plots          ",                               &  
     & " coil current info              "/
!
      nframe = nframe + 1
      if(nframe.ge.2000) ineg=40
!
        ipmes = iptype
        if (ipmes.gt.16)   ipmes = 16
        if (iptype.eq.6)   go to 100
        if (iptype.ne.7 .and. iptype.ne.8)   then
        write (nsc1,1001)  mes(ipmes)
 1001   format (a32)
        else
        write (nsc1,1002)  mes(ipmes), kcycle
 1002   format (a24, i7)
                                             endif
  100 rewind nsc1
      read(nsc1,1100) (ititle(nframe,n),n=1,4)
!
!.....debug
!     write(nterm,1100)  (ititle(nframe,n),n=1,4)
!.end debug
 1100 format(4a8)
      rewind nsc1
!
!.....write identifier
      call setld(77._R8,6._R8,1,0,1,0)
      write(s100,8485) nframe
      call gtext(s100,80,0)
 8485 format(1x,i8)
      call setld(77._R8,5._R8,1,0,1,0)
      write(s100,8484) name(1)
      call gtext(s100,80,0)
      call setld(77._R8,4._R8,1,0,1,0)
      write(s100,8484) idate
      call gtext(s100,80,0)
!     call setld(77.,4.,1,0,1,0)
!     write(s100,8484) itime
!     call gtext(s100,80,0)
 8484 format(1x,a10 )
!
      call frame(0)
      return
!
      entry pltindx
      if(nframe.le.0) return
      write(nout,2000)
 2000 format(1h1,"   *** plot index ***",//,                             &  
     &" frame         description",12x,                                  &  
     &" frame         description",/)
      nhalf = (nframe+1)/2
      do 200 nn=1,nhalf
      np = nn + nhalf
      if(np.gt.nframe) go to 210
      write(nout,3000) nn,(ititle(nn,n),n=1,4),                          &  
     &              np,(ititle(np,n),n=1,4)
      go to 200
  210 write(nout,3000) nn,(ititle(nn,n),n=1,4)
 3000 format(i4,4a8,2x,i4,4a8)
  200 continue
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

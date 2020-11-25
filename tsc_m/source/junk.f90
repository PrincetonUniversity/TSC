      subroutine rstrt1
!......5.40 rstrt1
!
!.....special initialization for restart run
!
      USE CLINAM
      USE FVV1
      USE NONCOR
      USE RADTAB
      USE RUNAWAY
      USE SAPROP
      USE SPECIE
      USE TCVCOM
      USE WALLCL
      USE WALLP
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER ifirstr,nnplocal,i,ios5,iflag,ios15
!============
! idecl:  explicitize implicit REAL declarations:
!     REAL*8 rstrt2
!============
      data ifirstr/1/
      character*31 mess1
      character*17 mess2
!
      dimension nnplocal(27)
!============      
!
!
      do 10 i=1,27
   10 nnplocal(i) = nnp(i)
!
!
!................................................................
!
!.....open disk file sprsin
!
!................................................................
!     isprsi(1:6) = 'sprsin'
!     isprsi(7:7) = isuffix(1:1)
      if( numargs .lt. 1 ) then
         isprsi = 'sprsin' // isuffix(1:1)
      else
         isprsi = 'sprsin' // '.' // trim(suffix)
      end if
      open(nsprsi,file=trim(isprsi),status='old',iostat=ios5,            &  
     &    form='unformatted')
!................................................................
!
!.....read in all common blocks
!
!................................................................
      call bufini(nsprsi,nnp(1),nnp(27))
!
      iflag = 0
      do 20 i=1,27
      if(nnp(i) .ne. nnplocal(i)) iflag = 1
   20 continue
      if(iflag.eq.0) go to 30
      ineg=1
      write(nout,1005)
 1005 format(" mismatch:  parameters on restart file printed first")
      write(nout,1001)
 1001 format(10x," iformat       pnx  pnz coil pobs pngr nlim tpts",     &  
     &" glob nsav pneq nthe ppsi  pkw four  plw pimp  pne  pte node")
      write(nout,1000) (nnp(i),i=1,27)
      write(nout,1000) (nnplocal(i),i=1,27)
 1000 format(" parameters:",27i5)
      return
!
   30 continue
      call buf_ini(nsprsi,"ibgcm1","ifncm1aa")
      call buf_ini(nsprsi,"ibgcm1aa","ifncm1")
      call buf_ini(nsprsi,"ibgc1a","ifnc1a")
      call buf_in(nsprsi,"begcm2","fincm2aa")
      call buf_in(nsprsi,"begc2aa","fincm2")
      call buf_in(nsprsi,"begc2a","finc2a")
      call buf_in(nsprsi,"begc2b","finc2b")
      call buf_in(nsprsi,"begc2c","finc2c")
!
      call buf_in(nsprsi,"begcm3","fincm3")
      call buf_in(nsprsi,"begc3a","finc3a")
      call buf_in(nsprsi,"begc3b","finc3b")
      call buf_in(nsprsi,"begc3c","finc3c")
      call buf_in(nsprsi,"begc3d","finc3d")
      call buf_in(nsprsi,"begcm4","fincm4")
      call buf_in(nsprsi,"begc4a","finc4a")
      call buf_in(nsprsi,"begc4b","finc4ba")
      call buf_in(nsprsi,"begc4ba","finc4b")
      write(6,*) "plhamp :", plhamp
      write(6,*) "picrh  :", picrh
      write(6,*) "frcparv:", frcparv
      write(6,*) "betarv:", betarv
      call buf_in(nsprsi,"begc4c","finc4c")
      call buf_in(nsprsi,"begcm5","fincm5")
      call buf_in(nsprsi,"begc5a","finc5a")
      call buf_in(nsprsi,"begc5b","finc5b")
      call buf_in(nsprsi,"begc5c","finc5c")
      call buf_in(nsprsi,"begcap","fincap")
      call buf_in(nsprsi,"begc6a","finc6a")
      call buf_in(nsprsi,"begcot","fincot")
!
      call buf_in(nsprsi,"begcsa","fincsa")
      call buf_in(nsprsi,"begcim","fincim")
      call buf_in(nsprsi,"begcsp","fincsp")
      call buf_in(nsprsi,"begcrt","fincrt")
      call buf_in(nsprsi,"plascur","vvfin")
      call buf_in(nsprsi,"begtcv","fintcv")
!
      call buf_in(nsprsi,"begufc","finufc")
!
      call buf_in(nsprsi,"begrun","finrun")
      call buf_in(nsprsi,"begcwa","fincwa")
!
!...reset some ifrst values :
!
      ifrst(3) = 1
      ifrst(4) = 1
      ifrst(5) = 0
      ifrst(6) = 0
      ifrst(9) = 1
!
      write(nout,1112) kcycle,times
      write(nterm,1112) kcycle,times
 1112 format(" * * * restart file read n=",i7," t=",1pe14.6," * * *")
!
      if(ilhcd.ge.0 .and. ifk.eq.2 .and. bpowsmw .gt. 0)                 &  
     &  call getlscrestart
      close(nsprsi)
      return
!................................................................
!
!.....write out common blocks for restart run
!
!................................................................
      entry rstrt2
!
!
!.....create disk file sprsou and change name of old file to osprsou
!
!................................................................
!     isprsp(1:7) = 'osprsou'
!     isprsp(8:8) = isuffix(1:1)
!     mess2(1:9) = 'rm `pwd`/'
!     mess2(10:17)=isprsp(1:8)
!     istat = ishell(mess2)
!     if(istat.lt.0) write(nterm,6331)istat, mess2
!6331 format(" istat=",i3," for ishell with ",a17)
!     isprso(1:6) = 'sprsou'
!     isprso(7:7) = isuffix(1:1)
!     mess1(1:9) = 'mv `pwd`/'
!     mess1(10:16)= isprso(1:7)
!     mess1(17:23) = ' `pwd`/'
!     mess1(24:31) = isprsp(1:8)
!     if(ifirstr.eq.0) istat=ishell(mess1)
!     if(istat.lt.0 .and. ifirstr.eq.0) write(nterm,6332) istat,mess1
!6332 format(" istat=",i3," for ishell with ",a31)
      if( numargs .lt. 1 ) then
         isprsp = 'osprsou' // isuffix(1:1)
         isprso = 'sprsou'  // isuffix(1:1)
      else
         isprsp = 'osprsou' // '.' // trim(suffix)
         isprso = 'sprsou'  // '.' // trim(suffix)
      endif
         
!     if(ifirstr.eq.0) call rename(isprso(1:7) , isprsp(1:8))
      if(ifirstr.eq.0) call rename(trim(isprso) , trim(isprsp))
      open(nsprso,file=trim(isprso),status='unknown',form='unformatted',  &  
     &     iostat=ios15)
      ifirstr=0
!................................................................
!
!.....output all common blocks
!
!................................................................
      call bufouti(nsprso,nnp(1),nnp(27))
!
      call buf_outi(nsprso,"ibgcm1","ifncm1aa")
      call buf_outi(nsprso,"ibgcm1aa","ifncm1")
      call buf_outi(nsprso,"ibgc1a","ifnc1a")
      call buf_out(nsprso,"begcm2","fincm2aa")
      call buf_out(nsprso,"begc2aa","fincm2")
      call buf_out(nsprso,"begc2a","finc2a")
      call buf_out(nsprso,"begc2b","finc2b")
      call buf_out(nsprso,"begc2c","finc2c")
!
      call buf_out(nsprso,"begcm3","fincm3")
      call buf_out(nsprso,"begc3a","finc3a")
      call buf_out(nsprso,"begc3b","finc3b")
      call buf_out(nsprso,"begc3c","finc3c")
      call buf_out(nsprso,"begc3d","finc3d")
      call buf_out(nsprso,"begcm4","fincm4")
      call buf_out(nsprso,"begc4a","finc4a")
      call buf_out(nsprso,"begc4b","finc4ba")
      call buf_out(nsprso,"begc4ba"",finc4b")
      write(6,*) "plhamp :", plhamp
      write(6,*) "picrh  :", picrh
      write(6,*) "frcparv:", frcparv
      write(6,*) "betarv:", betarv
      call buf_out(nsprso,"begc4c","finc4c")
      call buf_out(nsprso,"begcm5","fincm5")
      call buf_out(nsprso,"begc5a","finc5a")
      call buf_out(nsprso,"begc5b","finc5b")
      call buf_out(nsprso,"begc5c","finc5c")
      call buf_out(nsprso,"begcap","fincap")
      call buf_out(nsprso,"begc6a","finc6a")
      call buf_out(nsprso,"begcot","fincot")
!
      call buf_out(nsprso,"begcsa","fincsa")
      call buf_out(nsprso,"begcim","fincim")
      call buf_out(nsprso,"begcsp","fincsp")
      call buf_out(nsprso,"begcrt","fincrt")
      call buf_out(nsprso,"plascur","vvfin")
      call buf_out(nsprso,"begtcv","fintcv")
!
      call buf_out(nsprso,"begufc","finufc")
!
      call buf_out(nsprso,"begrun","finrun")
      call buf_out(nsprso,"begcwa","fincwa")
!
      close(nsprso)
      write(nout,1111) kcycle,times
      write(nterm,1111) kcycle,times
 1111 format(" * * restart file written n=",i7," t=",1pe14.6," * *")
!
      if(ilhcd.ge.0 .and. ifk.eq.2 .and. bpowsmw .gt. 0)                 &  
     &  call putlscrestart
!
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

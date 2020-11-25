 
 
 
 
 
 
 
 
 
 
 
 
 
 
!
!     level  -- explanation
!     0      plot not called yet
!     1      a plot file is opened
!     2      limits for plotting set
!     3      a plot is ready, or finished
 
!
!----------------------------------------------------------------------
!
      SUBROUTINE GNinit(iarg)
      USE gnuI
      USE tscunits
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
 
!     status: NEW        OLD    UNKNOWN     SCRATCH
!     access: SEQUENTIAL APPEND TRANSPARENT DIRECT
      INTEGER iarg
      INTEGER ifirst
      DATA    ifirst /                                                   &  
     &        1      /
 
      CHARACTER*19 myfile
      CHARACTER*18 myout
 
      CHARACTER*8  zzdate
      CHARACTER*10 zztime
      CHARACTER*5  zzzone
!
!     zzdate  ccyymmdd    cc=century yy=year mm=month dd=day
!     zztime  hhmmss.ttt  hh=hour mm=minute ss=second ttt=millisecond
!     zzzone  Shhmm       S=sign (-=west, +=east) hh=hour mm=minute
!
      INTEGER ivals(8)
!                             EXAMPLE
!     ivals ( 1 ) year        "2000"
!     ivals ( 2 ) month       "10"
!     ivals ( 3 ) day         "12"
!     ivals ( 4 ) t - UTC     "-0400"
!     ivals ( 5 ) hour        "8"
!     ivals ( 6 ) minute      "21"
!     ivals ( 7 ) second      "23"
!     ivals ( 8 ) millisecond "621"
!
 
      if (ifirst .eq. 1) then
          level = 0
          ifirst = 0
      endif
      iUnit = iarg
 
 
      if (level .eq. 0)      then
         call date_and_time(zzdate,zztime,zzzone,ivals)
         write(myfile,'(a)') zzdate(3:4)//'-'//zzdate(5:6)//'-'//        &  
     &                       zzdate(7:8)//'-'//zztime(1:6)//             &  
     &                       '.gnp'
 
         write(myout ,'(a)') zzdate(3:4)//'-'//zzdate(5:6)//'-'//        &  
     &                       zzdate(7:8)//'-'//zztime(1:6)//             &  
     &                       '.ps'
!
!        myfile
!        yy-mm-dd-hhmmss.gnp
!        123456789 123456789
!
!        myout
!        yy-mm-dd-hhmmss.ps
!        123456789 12345678
!
         open (iUnit, status='unknown',                                  &  
     &                file = myfile,                                     &  
     &                err=1300 )
         write (iUnit,'(a)') 'set terminal postscript'
         write (iUnit,'(a)') 'set output "'//myout//'"'
         write (iUnit,'(a)') 'set nokey'
         write (iUnit,'(a)') '# terminal opened by GNinit'
         level = 1
      else if (level .eq. 1) then
         write (iUnit,'(a)') 'set multiplot # by GNinit'
         write (iUnit,'(a)') 'set nolabel   # by GNinit'
         level = 2
      else if (level .eq. 2) then
         level = 3
      else
         level = 3
      endif
 
      return
 1300 write(nTSCscrn,'('' cannot open GNUplot file: '', a )') myfile
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"

      subroutine abexit
!......1.20 abexit
!                               R.O.S.  7 Dec 1987
!
!......writes an error message if ineg .ne. 0
!......then returns to main program for normal exit
!
      USE CLINAM
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!----------------------------------------------------------------------------
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER j,infeed,i,nr1
!============
      character*56 mes(70)
!
        data  (mes(j),j=1,10) /                                          &  
     & " error-parameters do not match at restart time          ",       &  
     & " error-coil or wire input fields not defined            ",       &  
     & " error-two wires with same indices                      ",       &  
     & " error--time step too small                             ",       &  
     & " error-- zero resistivity in wire                       ",       &  
     & " error-- parameters too small                           ",       &  
     & " error-observation point lies outside grid              ",       &  
     & " error-wire lies near or outside grid boundary          ",       &  
     & " error-limiter lies outside grid                        ",       &  
     & " error -no plasma current, check xlim and zlim          "/
        data  (mes(j),j=11,20) /                                         &  
     & " error-plasma current centroid lies outside grid        ",       &  
     & " error-initial equilibrium not converged                ",       &  
     & " error-out of range in g evaluation                     ",       &  
     & " trouble in flxvol                                      ",       &  
     & " ier ne 0 after call to nag library                     ",       &  
     & " ran out of toroidal flux zones, tfmult too small       ",       &  
     & "  error in readcard                                     ",       &  
     & " mag axis out of range or equal a wire                  ",       &  
     & "  trouble in facr                                       ",       &  
     & " error in feedback system definition                    "/
        data  (mes(j),j=21,30) /                                         &  
     & " input variable out of range                            ",       &  
     & " ******* negative jacobian in metriw ********           ",       &  
     & " *******  hyper did not converge *******                ",       &  
     & " error in externally supplied transport sub             ",       &  
     & "   error...parameter spda too small in spdata           ",       &  
     & " coil(type 09, lies inside grid boundary                ",       &  
     & " error in divertor plate specification (type32)         ",       &  
     & " error in svd analysis                                  ",       &  
     & " udsd not defined (type 12 card)                        ",       &  
     & " coil param not defined(type 39 card) for icirc=1       "/
        data  (mes(j),j=31,40) /                                         &  
     & " trouble in balloon eqn                                 ",       &  
     & " error in ripple calculation                            ",       &  
     & " error--too little poloidal flux in plasma              ",       &  
     & " psilim.lt.psimin.. check lims(or x-points search)      ",       &  
     & " psilim.lt.psimin..lower acoef(50),(51)                 ",       &  
     & " bc trouble in bounda                                   ",       &  
     & " enin tape in spdtftr wrong format                      ",       &  
     & " alphar and/or betar zero on type 04 card               ",       &  
     & " error in multipolar coil definition  or evaluation     ",       &  
     & " too many graphics frames                               "/
      data (mes(j),j=41,50) /                                            &  
     & " iteration  in subroutine iterate did not converge      ",       &  
     & " error in values of acoef(801)-acoef(805) for ffac<0    ",       &  
     & " error .... negative temperatures in trcdef             ",       &  
     & " error, v-vessel (type 10) incomplete with acoef(1).ne.0",       &  
     & " error in subroutine lsc:  ierr .ne. 0 on exit          ",       &  
     & " error in subroutine missionc                           ",       &  
     & " error in q-iteration converging in subroutine flxvol   ",       &  
     & " error in bootstrap calculation                         ",       &  
     & " g**2 negative in defg routine                          ",       &  
     & " one or more subroutines missing for options specified  "/
        data  (mes(j),j=51,60) /                                         &  
     & " error in tieval                                        ",       &  
     & " error in peval                                         ",       &  
     & " error in pevalo                                        ",       &  
     & " error in eeval                                         ",       &  
     & " error in reval                                         ",       &  
     & "  error in mmm95, nerr > 0 upon exit                    ",       &  
     & "  Max number of sawteeth exceeded...increase psaw       ",       &  
     & "  error in subroutine moviedata                         ",       &  
     & "  same acoef defined twice                              ",       &  
     & "  negative hydrogen density                             "/
      data (mes(j),j=61,70) /                                            &  
     & "  error in reading or writing to the plasma state       ",       &
     & "  error in trxpl profile fetching                       ",       &  
     & "  array size exceeded                                   ",       &
     & "  error 64                                              ",       &  
     & "  error 65                                              ",       &  
     & "  error 66                                              ",       &  
     & "  error 67                                              ",       &  
     & "  error 68                                              ",       &  
     & "  error 69                                              ",       &  
     & "  error 70                                              "/
!----------------------------------------------------------------------------
!
        if(ineg.eq.0)   return
        if(ineg.lt.0 .or. ineg.gt.70)  then
                                       write (nout,10)  ineg
   10   format (" abnormal exit...ineg= ",i4,"  out of range !!")
                                       return
                                       endif
        write (nout, 20)  ineg, mes(ineg)
        write (nterm,20)  ineg, mes(ineg)
   20   format (" abnormal exit...ineg= ",i4, 1x, a56)
        if (ineg.ne.6)   return
!
      write(nout,30)
   30 format(1h ,//," param  required   actual error ")
      infeed = 0
      do 40 i=1,numfb
   40 if(idelay(i).gt.infeed) infeed = idelay(i)
!
        call  parwrt (nout, " pdelay ", infeed, nnp(2))
        call  parwrt (nout, " pnx    ", nx+2,   nnp(3))
        call  parwrt (nout, " pnz    ", nz+2,   nnp(4))
      nr1 = ncoil
      if(nwire.gt.ncoil) nr1 = nwire
        call  parwrt (nout, " pncoil ", nr1,    nnp(5))
        call  parwrt (nout, " pnwire ", nwire, nnp(27))
        call  parwrt (nout, " ppsi   ", npsi+1, nnp(14))
        call  parwrt (nout, " pimp   ", 2,   nnp(18))
        call  parwrt (nout, " pobs   ", nobs,   nnp(6))
        call  parwrt (nout, " pngroup", ngrmax, nnp(7))
        call  parwrt (nout, " pnlim  ", nlim,   nnp(8))
        call  parwrt (nout, " ptpts  ", ntpts,  nnp(9))
        call  parwrt (nout, " pnplat ", nplate, nnp(22))
        call  parwrt (nout, " pnsep  ", nsepmax,nnp(24))
        call  parwrt (nout, " pnfeed ", numfb  ,nnp(25))
        call  parwrt (nout, " pngrps ", ngroups,nnp(26))
        if(ibalsw.le.0 .and. iwall.le.0 .and. irippl.le.0)   return
        call  parwrt (nout, " pnthe  ", nthe+4, nnp(13))
        return
        end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

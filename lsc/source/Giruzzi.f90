 
!     pbxio.f(or)   starts      ---------------------------------------|
!     -                                                                |
!     -                                                                |
!
      SUBROUTINE Giruzzi
      USE DqlBins
      USE FeBins
      USE params
      USE ProfBody
      USE RayBins
      USE RayWrk
      USE WkAry
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER iv, ip, iNotSmoo, iYesSmoo, ipstep, ThisPlot, NumPlts
      INTEGER iFullFe, iOrigMaxwl
      DATA iNotSmoo, iYesSmoo / 1 , 2 /
 
      CHARACTER*1 curve, dots, altdots
      CHARACTER*10 chip
      DATA curve, dots,altdots, NumPlts  / ' ', '.', ',', 12 /
      ipstep = (npsi - 2)/NumPlts
      iv = 1
        iFullFe    = iITR
        iOrigMaxwl = mod(iITR,2) + 1
      ThisPlot = 1
      do 10 ip = ipstep, ipstep*NumPlts, ipstep
!     .                                 if chip starts with a 'l' or 'L'
!     .                                 then the logarithm is plotted
        write(chip,'(''Dql('',i2,'')'')') ip
!
!     .                                 call the curve line first... this
!     .                                 sets the scale of the graph
        call Dql6Norm(Vpar(iv),Dql(iv,ip,iYesSmoo),nv,                   &  
     &                chip,curve, ThisPlot, NumPlts)
!     .                                 call the dots second....previous scale
!     .                                 will be used, and curve clipped if
!     .                                 necessary
        call Dql6Norm(Vpar(iv),Dql(iv,ip,iNotSmoo),nv,                   &  
     &                chip,dots , ThisPlot, NumPlts)
!
        ThisPlot = ThisPlot + 1
!
 10   continue
 
 
      iv = 1
      ThisPlot = 1
      do 20 ip = ipstep, ipstep*NumPlts, ipstep
!     .                                 if chip starts with a 'l' or 'L'
!     .                                 then the logarithm is plotted
        write(chip,'(''logFe ('',i2,'')'')') ip
!
!     .                                 call the curve line first... this
!     .                                 sets the scale of the graph
        call Dql6Norm(Vpar(iv),fe(iv,ip,iFullFe),nv,                     &  
     &                chip,curve, ThisPlot, NumPlts)
!     .                                 call the dots second....previous scale
!     .                                 will be used, and curve clipped if
!     .                                 necessary
        call Dql6Norm(Vpar(iv),fe(iv,ip,iOrigMaxwl),nv,                  &  
     &                chip,dots , ThisPlot, NumPlts)
!
        ThisPlot = ThisPlot + 1
!
 20   continue
 
 
      return
      END
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"

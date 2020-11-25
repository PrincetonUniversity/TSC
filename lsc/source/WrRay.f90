!
!     Rayini.F  ends                           ------------------------|
!     Rayio.F  begins                       ---------------------------|
!     -                                                                |
!     -                                                                |
!
!     WrRay     begins                      ---------------------------|
!     -                                                                |
!     -                                                                |
      SUBROUTINE WrRay
      USE Doflags
      USE params
      USE RayBins
      USE TSCgrap
      USE tscunits
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      CHARACTER*16 FileNa
      INTEGER i,j, nzThisRay
!
      if (DoTRAN .eq. 1) then
        FileNa = Equhd(1) // Equhd(2)            ! Transp Run ID // 'RAY.DAT'
        If (FileNa .Eq. ' ') FileNa = 'ray.dat'  ! NonTransp name
      else
        FileNa = 'ray.dat'
      endif
!
      open (nTSCunus,status='unknown',file = FileNa, err=1300)
      write(nTSCunus,'('' nzones, nrays, npsi, iLastRy are: '')' )
      write(nTSCunus,'( 1x, 4i5 )' )  nzones, nrays, npsi, iLastRy
 
      do 10 j=1,nrays
        nzThisRay=nzones
        do 5 i=1,nzones
          if(ivind(i,j) .ne. 0) goto 5
          nzThisRay=i
          goto 6
 5      continue
 6      continue
      write(nTSCunus,                                                    &  
     &     '(1x, i4, '' izind for ray number '', i4 )') nzThisRay,j
      write(nTSCunus,600)                (izind(i,j)  ,i=1,nzThisRay)
      write(nTSCunus,                                                    &  
     &     '(1x, i4, '' ivind for ray number '', i4 )') nzThisRay, j
      write(nTSCunus,600)                (ivind(i,j)  ,i=1,nzThisRay)
      write(nTSCunus,'(                                                  &  
     & '' ezsq; npar; dlnPdsX; dlnPdsK for ray number '', i4 )') j
      write(nTSCunus,601)                 (ezsq(i,j)  ,i=1,nzThisRay)
      write(nTSCunus,601)                 (npar(i,j)  ,i=1,nzThisRay)
      write(nTSCunus,601)              (dlnPdsX(i,j)  ,i=1,nzThisRay)
      write(nTSCunus,601)              (dlnPdsK(i,j)  ,i=1,nzThisRay)
 10   continue
 600  format( (20i4)    )
 601  format( (1x, 5e15.8) )
 
      close(nTSCunus)
      return
 1300 continue
      call LSCstop (' error opening file in WrRay')
      return
      END
!     -                                                                |
!     -                                                                |
!     WrRay     ends                        ---------------------------|
! 24May2005 fgtok -s r8_precision.sub "r8con.csh conversion"

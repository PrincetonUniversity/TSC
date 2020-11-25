      subroutine map12(sary1d,sary2d)
!
!
!.....interpolates from 1D aray SARY1D(PPSI) (cell centered)
!                       (corresponding to flux values XSV(PPSI)
!.....defines values on 2D spatial array SARY2D(PNX,PNZ) (cell centered)
!
      USE CLINAM
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER i,j,l
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 sary2d,sary1d,psival,denom
!============
!     dimension sary1d(ppsi),sary2d(pnx,pnz),xsvm(ppsi),xsv2m(ppsi)
      dimension sary1d(ppsi),sary2d(pnx,pnz)
!============      
      INTEGER :: istat = 0 
      REAL*8, ALLOCATABLE, DIMENSION(:) :: xsvm
      REAL*8, ALLOCATABLE, DIMENSION(:) :: xsv2m
!============      
      IF(.not.ALLOCATED(xsvm)) ALLOCATE( xsvm(ppsi), STAT=istat)
      IF(.not.ALLOCATED(xsv2m)) ALLOCATE( xsv2m(ppsi), STAT=istat)
!============      
      if (istat .ne. 0) stop 'Allocation Error : map12  ' 
!============      
!
!.....initialize to zero
      do 100 i=3,nxp
      do 100 j=3,nzp
  100 sary2d(i,j) = 0
!
      do 105 j = 1,npsi
      xsvm (j) = xsv (j) - (xsv2(1)-psimin)
      xsv2m(j) = xsv2(j) - (xsv2(1)-psimin)
  105 continue
!
!.....use psi value at cell center to interpolate from sary1d aray
!
      do 201 i=iminn+1,imaxx
      do 200 j=jminn+1,jmaxx
!
!.....exclude cells not in plasma
      if(iexvc(i,j).gt.0 .or. iexs(i,j).eq.1) go to 200
      psival = .25_R8*(psi(i,j)+psi(i-1,j)+psi(i,j-1)+psi(i-1,j-1))
      if(psival.gt.xsv2m(npsit)) go to 200
      if(psival.gt.xsvm (npsit)) go to 190
      if(psival.lt.xsvm (2)    ) go to 180
!
      l = 2
  150 l = l+1
      if(l.gt.npsit) go to 200
      if(psival.lt.xsvm(l)) go to 170
      go to 150
  170 continue
!
!.....  xsvm(l-1)  < psival < xsvm(l)
!
      denom = xsvm(l) - xsvm(l-1)
      if(denom.eq.0) go to 200
      sary2d(i,j) = ((psival-xsvm(l-1))*sary1d(l)                        &  
     &            +  (xsvm(l)  -psival)*sary1d(l-1))/denom
      go to 200
  180 continue
!
!.....  xsv2m(1)   < psival < xsvm(2)
      sary2d(i,j) = sary1d(2)
      go to 200
  190 continue
!
!.....xsv2m(npsit) < psival < xsvm(npsit)
      sary2d(i,j) = sary1d(npsit)
  200 continue
      if(isym.ne.1) go to 201
      sary2d(i,2) = sary2d(i,3)
  201 continue
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

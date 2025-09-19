      subroutine iterate(itype,aa1)
!.....* * * number 4.60 * * *
!
      USE CLINAM
      USE SCRATCH
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER itype,imin2,i,j,imaxa,jmaxa,ibna,ibnb,ibnc
      INTEGER ibnd,iskip,ni,ict,nisave
!============
! idecl:  explicitize implicit REAL declarations:
!     REAL*8 aa1,vec,aa2,aa3,grsb,grsd,facconv,gtsa,gtsc,denom,facm
      REAL*8 aa1
      REAL*8 grsb,grsd,facconv,gtsa,gtsc,denom,facm
      REAL*8 usum,umaxa,omsum,umax,omcor,anum1,term2,diffmx,error
      REAL*8 sum
!============
!     equivalence (vec(1,1,1),bigcom(1)),
!    3   ( aa2(1,1),bigcom(pnxz7)),
!    4   ( aa3(1,1),bigcom(pnxz8))
!     dimension aa2(pnx,pnz),aa3(pnx,pnz)
!     dimension vec(7,pnz,pnx)
!!!!  vec, aa2, aa3 may be locally allocated and deallocated

      dimension aa1(pnx,pnz)
!               ,jmin2y(pnx)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: vec
      REAL*8, ALLOCATABLE, DIMENSION(:,:)   :: aa2, aa3
!============      
      INTEGER :: istat = 0 
      INTEGER, ALLOCATABLE, DIMENSION(:) :: jmin2y
!============      
!
      if(itype.le.0 .or. itype.gt.3) return

      IF(.not.ALLOCATED(jmin2y)) ALLOCATE( jmin2y(pnx), STAT=istat)
      if(.NOT.ALLOCATED(vec)) ALLOCATE (vec(7,pnz,pnx), STAT=istat)
      if(.NOT.ALLOCATED(aa2)) ALLOCATE (aa2(pnx,pnz), STAT=istat)
      if(.NOT.ALLOCATED(aa3)) ALLOCATE (aa3(pnx,pnz), STAT=istat)
      if (istat .ne. 0) stop "Allocation Error : iterate"
!============      

      go to(701,702,703),itype
  701 continue
!
!.....space reserved for poloidal flux operator
!
!     DEALLOCATE(vec,aa2,aa3)
      return
  702 continue
      imin2 = iminn
!
!
!.....itype=2 for velocity stream function
!
      grsb = deex/deez
      grsd = deex/deez
      phi2 = acoef(831)
      sf = acoef(832)
      facconv = acoef(833)
      nimax = acoef(834)
!
!
      do 100 i=2,nxp
      jmin2y(i) = jminny(i)
      gtsa = (xary(i)/xarh(i+1))*deez/deex
      gtsc = (xary(i)/xarh(i  ))*deez/deex
!
      do 150 j=2,nzp
      denom = 1._R8+ phi2*sf/nx
      facm = sf**2/(grsb+grsd+gtsa+gtsc)/denom
!
      if(iexv(i,j).eq.1) go to 140
!
      aa2(i,j) = aa1(i,j)
      aa3(i,j) = aa1(i,j)
      vec(1,j,i) = -(1._R8-phi2*sf/nx)/denom
      vec(4,j,i) =-facm*(grsb+grsd+gtsa+gtsc) + 2._R8/denom
      vec(2,j,i) = facm*(grsd)
      vec(3,j,i) = facm*(gtsc)
      vec(5,j,i) = facm*(gtsa)
      vec(6,j,i) = facm*(grsb)
      vec(7,j,i) =-facm*xary(i)*b(i,j)
      go to 150
  140 continue
      aa2(i,j) = 0._R8
      aa3(i,j) = 0._R8
      vec(1,j,i)=0._R8
      vec(4,j,i)=0._R8
      vec(2,j,i)=0._R8
      vec(3,j,i)=0._R8
      vec(5,j,i)=0._R8
      vec(6,j,i)=0._R8
      vec(7,j,i)=0._R8
  150 continue
  100 continue
      go to 800
  703 continue
      imin2 = iminn+1
!
!.....itype=3 for velocity potential
!
!
      grsb = deex/deez
      grsd = deex/deez
      phi2 = acoef(841)
      sf = acoef(842)
      facconv = acoef(843)
      nimax = acoef(844)
      usum=0._R8
      umaxa = 0._R8
      sum = 0._R8
      omsum = 0
      do 750 i=3,nxp
      jmin2y(i) = jminny(i)+1
      do 750 j=3,nzp
      if(iexvc(i,j).gt.0) go to 750
      if(abs(u(i,j)) .lt. umaxa) go to 749
      umax = u(i,j)
      umaxa = abs(u(i,j))
      imaxa = i
      jmaxa = j
  749 continue
      usum = usum + u(i,j)
      sum = sum + 1._R8
      omsum = omsum + omeg(i,j)
  750 continue
!     write(nterm,5959)
!5959 format(".phi2     .sf       .facconv  .nimax    .")
!     read(59,5960) phi2,sf,facconv,nimax
!5960 format(3e10.2,i10)
      ucor = usum/sum
      omcor = omsum/sum
!
      do 101 i=3,nxp
      gtsa = (xary(i  )/xarh(i))*deez/deex
      gtsc = (xary(i-1)/xarh(i))*deez/deex
!
      do 151 j=3,nzp
!
      ibna = 1
      ibnb = 1
      ibnc = 1
      ibnd = 1
!
      if(iexvc(i,j).gt.0) go to 141
      u(i,j) = u(i,j)-ucor
      uo(i,j) = uo(i,j) - ucor
      aa1(i,j)=aa1(i,j)-omcor
!
      if(i.ne.nxp) go to 91
      ibna = 0
   91 if(i.ne. 3 ) go to 90
      ibnc = 0
   90 continue
      if(iexv(i,j).eq.0.or.iexv(i,j-1).eq.0) go to 94
      ibna=0
   94 continue
      if(iexv(i-1,j).eq.0.or.iexv(i,j).eq.0) go to 95
      ibnb=0
   95 continue
      if(iexv(i-1,j-1).eq.0.or.iexv(i-1,j).eq.0) go to 96
      ibnc=0
   96 continue
      if(iexv(i,j-1).eq.0.or.iexv(i-1,j-1).eq.0) go to 97
      ibnd=0
   97 continue
      if(j.eq.3 .and. isym.eq.0) ibnd=0
      if(j.eq.nzp)ibnb=0
!
      anum1 = sf**2/(grsb+grsd+gtsa+gtsc)
      term2 = .5_R8*anum1*(gtsa*(1._R8-ibna)+gtsc*(1._R8-ibnc)           &  
     &                 +grsb*(1._R8-ibnb)+grsd*(1._R8-ibnd) )
      denom = 1._R8+ phi2*sf/nx - term2
      facm = anum1/denom
!
      aa2(i,j) = aa1(i,j)
      aa3(i,j) = aa1(i,j)
      vec(1,j,i) = -(1._R8-phi2*sf/nx-term2)/denom
      vec(4,j,i)=-facm*(grsb+grsd+gtsa+gtsc)+2._R8/denom
      vec(2,j,i) = facm*(grsd*ibnd)
      vec(3,j,i) = facm*(gtsc*ibnc)
      vec(5,j,i) = facm*(gtsa*ibna)
      vec(6,j,i) = facm*(grsb*ibnb)
      vec(7,j,i) =-facm*u(i,j)/xarh(i)
      go to 151
  141 continue
      aa2(i,j) = 0._R8
      aa3(i,j) = 0._R8
      vec(1,j,i)=0._R8
      vec(4,j,i)=0._R8
      vec(2,j,i)=0._R8
      vec(3,j,i)=0._R8
      vec(5,j,i)=0._R8
      vec(6,j,i)=0._R8
      vec(7,j,i)=0._R8
  151 continue
      if(isym.eq.0) go to 101
      aa1(i,2) = aa1(i,3)
      aa2(i,2) = aa2(i,3)
      aa3(i,2) = aa3(i,3)
      u(i,2) = u(i,3)
      uo(i,2) = uo(i,3)
  101 continue
!     write(nterm,5101) usum
 5101 format(" usum = ",1pe12.4)
!
  800 continue
      iskip=0
      ni = 0
  300 ni = ni + 1
!
      ict = 0
      diffmx = 0
!
!
!
      iskip = iskip + 1
!
      do 20 i=imin2,imaxx
      if(jmaxxy(i).lt.jmin2y(i)) go to 20
      do 10 j=jmin2y(i),jmaxxy(i)
!
      aa3(i,j) = vec(1,j,i)*aa1(i,j)                                     &  
     &        + vec(2,j,i)*aa2(i  ,j-1)                                  &  
     &        + vec(3,j,i)*aa2(i-1,j  )                                  &  
     &        + vec(4,j,i)*aa2(i  ,j  )                                  &  
     &        + vec(5,j,i)*aa2(i+1,j  )                                  &  
     &        + vec(6,j,i)*aa2(i  ,j+1)                                  &  
     &        + vec(7,j,i)
!
   10 continue
      if(isym.eq.0) go to 20
      if(itype.eq.2) then
                     aa3(i,1)=-aa3(i,3)
                     aa3(i,2)=0._R8
                     endif
      if(itype.eq.3) then
                     aa3(i,2) = aa3(i,3)
                     endif
   20 continue
!
!
!
      do 40 i=imin2,imaxx
      if(jmaxxy(i).lt.jmin2y(i)) go to 40
      do 30 j=jmin2y(i),jmaxxy(i)
!
      aa1(i,j) = vec(1,j,i)*aa2(i,j)                                     &  
     &        + vec(2,j,i)*aa3(i  ,j-1)                                  &  
     &        + vec(3,j,i)*aa3(i-1,j  )                                  &  
     &        + vec(4,j,i)*aa3(i  ,j  )                                  &  
     &        + vec(5,j,i)*aa3(i+1,j  )                                  &  
     &        + vec(6,j,i)*aa3(i  ,j+1)                                  &  
     &        + vec(7,j,i)
!
   30 continue
      if(isym.eq.0) go to 40
      if(itype.eq.2) then
                     aa1(i,1)=-aa1(i,3)
                     aa1(i,2)=0._R8
                     endif
      if(itype.eq.3) then
                     aa1(i,2) = aa1(i,3)
                     endif
   40 continue
!
      do 60 i=imin2,imaxx
      if(jmaxxy(i).lt.jmin2y(i)) go to 60
      do 50 j=jmin2y(i),jmaxxy(i)
!
      aa2(i,j) = vec(1,j,i)*aa3(i,j)                                     &  
     &        + vec(2,j,i)*aa1(i  ,j-1)                                  &  
     &        + vec(3,j,i)*aa1(i-1,j  )                                  &  
     &        + vec(4,j,i)*aa1(i  ,j  )                                  &  
     &        + vec(5,j,i)*aa1(i+1,j  )                                  &  
     &        + vec(6,j,i)*aa1(i  ,j+1)                                  &  
     &        + vec(7,j,i)
!
   50 continue
      if(isym.eq.0) go to 60
      if(itype.eq.2) then
                     aa2(i,1)=-aa2(i,3)
                     aa2(i,2)=0._R8
                     endif
      if(itype.eq.3) then
                     aa2(i,2) = aa2(i,3)
                     endif
   60 continue
!
      if(iskip.lt.10) go to 300
      iskip = 0
!
      sumb = 0
      sumdb = 0
      do 70 i=imin2,imaxx
      if(jmaxxy(i).lt.jmin2y(i)) go to 70
      do 80 j=jmin2y(i),jmaxxy(i)
!
      sumdb= sumdb+ abs(aa1(i,j)-aa2(i,j))
!
      sumb = sumb + abs(aa1(i,j))
!
!
   80 continue
   70 continue
!
      if(sumb.gt.1.E10_R8) go to 250
      if(sumdb .lt. facconv*sumb) go to 801
      if(sumdb .lt. 1.E-13_R8*nx*nz) go to 801
!
  200 continue
!
      if(ni.gt.nimax) go to 700
      go to 300
  250 continue
      write(nterm,9250) sf
      write(nout,9250) sf
 9250 format(" iteration diverging in iterate",                          &  
     &      /," lower sf...present value",f8.3 )
!
!     abnormal end of loop
700   continue
      ineg=41
      error = sumdb/sumb
      write(nterm,9260) error,facconv
      write(nout, 9260) error,facconv
      write(nterm,*) "nimax",nimax
      write(nout,*) "nimax",nimax
 9260 format(" iteration not converged in iterate",/,                    &  
     &       " error,faccon= ",1p2e12.4)
!
!     normal end of loop
  801 continue
      nisave = ni
      if(itype.le.2) then
!       DEALLOCATE(vec,aa2,aa3)
        return
      end if
!
!.....ITYPE=3 only
!
      do 660 i=imin2,imaxx
      if(jmaxxy(i).lt.jmin2y(i)) go to 660
      do 650 j=jmin2y(i),jmaxxy(i)
 
      ibna = 1
      ibnb = 1
      ibnc = 1
      ibnd = 1
!
      if(iexvc(i,j).gt.0) go to 650
!
      if(i.ne.nxp) go to 691
      ibna = 0
  691 if(i.ne. 3 ) go to 690
      ibnc = 0
  690 continue
      if(iexv(i,j).eq.0.or.iexv(i,j-1).eq.0) go to 694
      ibna=0
  694 continue
      if(iexv(i-1,j).eq.0.or.iexv(i,j).eq.0) go to 695
      ibnb=0
  695 continue
      if(iexv(i-1,j-1).eq.0.or.iexv(i-1,j).eq.0) go to 696
      ibnc=0
  696 continue
      if(iexv(i,j-1).eq.0.or.iexv(i-1,j-1).eq.0) go to 697
      ibnd=0
  697 continue
      if(j.eq.3 .and. isym.eq.0) ibnd=0
      if(j.eq.nzp) ibnb=0
      if(ibna.eq.0) aa1(i+1,j) = aa1(i,j)
      if(ibnb.eq.0) aa1(i,j+1) = aa1(i,j)
      if(ibnc.eq.0) aa1(i-1,j) = aa1(i,j)
      if(ibnd.eq.0) aa1(i,j-1) = aa1(i,j)
  650 continue
  660 continue
!
 
!     write(nterm,3059) itype,ni,usum,imaxa,jmaxa,umax
 3059 format(" itype,ni =",2i4," usum=",1pe10.2,2i3,1pe10.2)
      nisave=ni
!     DEALLOCATE (vec,aa2,aa3)
      return
!
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

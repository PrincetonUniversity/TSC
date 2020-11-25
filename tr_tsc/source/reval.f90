      subroutine reval(psv,itype,jsurf,rval,rpval,ig,jg)
!.....number 8.60
!
      USE trtsc
      USE CLINAM
      USE SAPROP
      USE SCR3
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!.....for itype=1,2  evaluate from analytic function
!.......  itype=0,3  evaluate from mu profile
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER itype,ig,jg,n1,n2,jv,j,lval,lv,lm,itypeselect
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 rval,rpval,psv,phiarg,phitil,phiargc,fac1,frac,psve,psl
      REAL*8 phiarg2,a3013,a3014
!============
      REAL*8, ALLOCATABLE, DIMENSION(:) :: xi_bdy, rmajmp, bmidp, ner
      INTEGER :: i, jsurf 
      REAL*8 :: xi, qtmp
!============
      logical :: first_read_ne=.true.
      logical :: ex_ne=.false. 
      real*8 :: tint, rint, rval1, rval2, rvalp, rvalm
      real*8, dimension(:),allocatable :: time_ne, rad_ne
      real*8, dimension(:,:),allocatable :: rdn
      integer :: j1, j2, kk, ii
      integer :: ntime_ne, nrad_ne
      real*8 :: tfluxd, tflux1, tflux2, denompsi,sq_tflx
      real*8 :: tfluxn, fluxrat, flfactor, rhox


      if(iexvc(ig,jg).gt.0 .or. iexs(ig,jg).eq.1)go to 40
      if(psv.ge.psilim) go to 40
!     modified following 3 lines  11/20/2015   (SCJ)
      if(itype.gt.3) go to 100   
      if((itype.eq.0 .or.itype.eq.3) .and. jsurf.eq.1) go to 10
!      if(jsurf.eq.1) go to 10
      if(ifunc.eq.6) go to 110

      itypeselect = itype
      if(itype.eq.0) itypeselect = 1
      select case (itypeselect)
    
      case (2)
!
!     using experimental profile from trdatbuf
!
      if(.not.allocated(xi_bdy))allocate(xi_bdy(npsi))
      if(.not.allocated(rmajmp))allocate(rmajmp(2*npsi-1))
      if(.not.allocated(bmidp))allocate(bmidp(2*npsi-1))
      if(.not.allocated(ner))allocate(ner(npsi))

      if( npsit .gt. 1 ) then
      xi_bdy(1) = 0.0_R8
      do i = 1, npsit-1
      xi_bdy(i+1) = xi_bdy(i) + 1.0_R8/float(npsit-1)
      enddo
      xi = 0._R8
      do i = 1, npsit-1
      if(psv .gt. xsv2(i) .and. psv .le. xsv2(i+1)) then
!     if(psv .gt. xsv(i) .and. psv .le. xsv(i+1)) then
!     xi = xi_bdy(i)+ (psv-xsv(i))/(xsv(i+1)-xsv(i))               &
!                     *(xi_bdy(i+1)-xi_bdy(i)) 
      xi = xi_bdy(i)*dpsi*float(npsit-1)                               &
                      +tpi*qprof2(i)*(psv-xsv2(i))                     &
                      +tpi*0.5_R8*(qprof2(i+1)-qprof2(i))/             &
                       (xsv2(i+1)-xsv2(i))*(psv-xsv2(i))**2
      xi = xi/(dpsi*float(npsit-1))
      exit
      endif
      enddo

      xi = sqrt(xi)
      do i= 1, npsit
      xi_bdy(i) = sqrt( xi_bdy(i) )
      enddo

      do i = 1, npsit-1
!     rmajmp(npsit-i) = rmajor + (1._R8-xi_bdy(i+1)**2)*(xmag-rmajor)  &
!                              - xi_bdy(i+1)*rminor
!     rmajmp(npsit+i) = rmajor + (1._R8-xi_bdy(i+1)**2)*(xmag-rmajor)  &
!                              + xi_bdy(i+1)*rminor

      rmajmp(npsit-i) = rmajora(i+1) - rminora(i+1)
      rmajmp(npsit+i) = rmajora(i+1) + rminora(i+1)

      bmidp(npsit-i) = gzero/rmajmp(npsit-i)
      bmidp(npsit+i) = gzero/rmajmp(npsit+i)
      enddo
      rmajmp(npsit) = xmag
      bmidp(npsit) = gzero/xmag
      rmajmp(:2*npsit-1) = rmajmp(:2*npsit-1)*100.0_R8

      call ne_from_trdatbuf(times,xi_bdy,rmajmp,bmidp,ner,npsit-1,ierr)
!     if(ierr .ne. 0) stop " Error in ne_from_trdatbuf. Program stopped"
      do i = 1, npsit-1
      if(xi .gt. xi_bdy(i) .and. xi .le. xi_bdy(i+1)) then
      rval = ner(i) + (xi-xi_bdy(i))/(xi_bdy(i+1)-xi_bdy(i))             &
                    * (ner(i+1)-ner(i))
      rpval = (ner(i+1)-ner(i))/(xsv(i+1)-xsv(i))
      exit
      endif
      enddo
      rval = rval*1.0e+06*usdd
      rpval = rpval*1.0e+06*usdd
      if( psv .le. psimin ) then
      rval = ner(1)*1.0e+06*usdd
      rpval = 0._R8
      endif
      else
      rval = 5.0e+17*usdd
      rpval = 0.0_R8
      endif
! --------------------------------------------------------------

      case (1)

      if(use_user_ne .and. first_read_ne) then
      inquire(file="user_ne_data", exist=ex_ne)
      if (ex_ne) then
      open(55,file="user_ne_data",form="formatted",status="old")
      read(55,*) ntime_ne, nrad_ne
      allocate(time_ne(ntime_ne), rad_ne(nrad_ne))
      allocate(rdn(nrad_ne,ntime_ne))
      read(55,*) time_ne(1:ntime_ne)
      read(55,*) rad_ne(1:nrad_ne)
      read(55,*) rdn(1:nrad_ne,1:ntime_ne)
      close (55)
      endif
      first_read_ne=.false.
      endif

      if(use_user_ne .and. ex_ne) then
      do  kk = 1,ntime_ne-1
      if(times .gt. time_ne(kk) .and. times .le. time_ne(kk+1)) then
      j1 = kk
      j2 = kk + 1
      tint = (times-time_ne(j1))/(time_ne(j2)-time_ne(j1))
      exit
      endif
      enddo
      if(times .le. time_ne(1)) then
      j1 = 1
      j2 = 2
      tint = 0.
      endif
      if(times .gt. time_ne(ntime_ne)) then
      j1 = ntime_ne-1
      j2 = ntime_ne
      tint = 1.
      endif

      if(psv .le. psimin) then
      rval = rdn(1,j1)+(rdn(1,j2)-rdn(1,j1))*tint
      rval = rval*1.e20*usdd/1.0
      rpval = 0.
      else
      do i=1,npsit-1
      if(psv .gt. xsv(i) .and. psv .le. xsv(i+1)) then
      tflux1 = float(i-1)*dpsi
      tflux2 = float(i)*dpsi
      tfluxd = tflux1 + ((tflux2-tflux1)/(xsv(i+1)-xsv(i)))             &
     & *(psv-xsv(i))
      tfluxd = tfluxd/(float(npsit-1)*dpsi)
      denompsi = (xsv(i+2)-xsv(i))
      exit
      endif
      enddo

      do i=1,nrad_ne 
         sq_tflx = sqrt(tfluxd)
!      if(tfluxd .gt. rad_ne(i) .and. tfluxd .le. rad_ne(i+1)) then
      if(sq_tflx .gt. rad_ne(i) .and. sq_tflx .le. rad_ne(i+1)) then
!      rint = (tfluxd-rad_ne(i))/(rad_ne(i+1)-rad_ne(i))
      rint = (sq_tflx-rad_ne(i))/(rad_ne(i+1)-rad_ne(i))
      rval1 = rdn(i,j1) + (rdn(i+1,j1)-rdn(i,j1))*rint
      rval2 = rdn(i,j2) + (rdn(i+1,j2)-rdn(i,j2))*rint
      rval = rval1 + (rval2-rval1)*tint
      rvalp = rdn(i+2,j1)+(rdn(i+2,j2)-rdn(i+2,j1))*tint
      rvalm = rdn(i  ,j1)+(rdn(i  ,j2)-rdn(i  ,j1))*tint
      rval = rval*1.e20*usdd/1.0
      rpval = (rvalp-rvalm)/denompsi
      rpval = rpval*1.e20*usdd/1.0
      exit
      endif
      enddo
      endif

      return
      endif
!============================================================================
!	modified by F.M.Poli on 08/12/2011
!	removed coding that uses constant n1, n2 and acoef(883) 
!	rpval is re-calculated from rval

      phiarg = 1.0_R8
      if(psv.gt.psimin)                                                 &
     &phiarg = 1._R8- ((psv-psimin)*delpsi)**betar
      if(nflag .gt. 0 .and. expn2 .ne. 0.)  &
     &   phiarg2 = 1._R8 - ((psv-psimin)*delpsi)**expn2
      phitil =  (psv-psimin)*delpsi
      rval = r0*(1._R8-fracn0-newden)*phiarg**alphar             &
     &  + fracn0*r0
      if(nflag .gt. 0 .and. expn2 .ne. 0.) rval = rval                  &
     & + newden*r0*phiarg2**expn1
!
       phiargc = 0._R8
      if(psv.gt.psimin) phiargc = (psv-psimin)*delpsi
       rpval = (-r0*(1.-fracn0-newden)*alphar*betar*              &
     &  (phiargc**(betar-1.))*                                          &
     &  ((1.-phiargc**betar)**(alphar-1.)))*delpsi
      if(nflag .gt. 0 .and. expn2 .ne. 0.) rpval = rpval                &
     & + ((-r0*newden)*expn1*expn2*(phiargc**(expn2-1._R8))*            &
     & (1.-phiargc**expn2)**(expn1-1._R8))*delpsi
!========================================================================
!     special coding for ITER
      if(int(acoef(4977)+0.1) .eq. 1) then
      tfluxn = 0._R8
      do 9000 i=1,npsit-1
      if(psv .gt. xsv(i) .and. psv .le. xsv(i+1)) then
      tflux1 = float(i-1)*dpsi
      tflux2 = float(i)*dpsi
      tfluxd = tflux1 + ((tflux2-tflux1)/(xsv(i+1)-xsv(i)))              &
     & *(psv-xsv(i))
      tfluxn = tfluxd/(float(npsit-1)*dpsi)
      fluxrat = (xsv(i+1)-xsv(i))/(tflux2-tflux1)
      flfactor = fluxrat*2.*sqrt(tfluxd*float(npsit-1)*dpsi)
      go to 9001
      endif
 9000 continue
 9001 continue
      rhox = sqrt(tfluxn)
      rval = r0*(1.-fracn0-newden)*((1.-rhox**betar)**alphar)    &
     &  + r0*fracn0
      if(nflag .gt. 0 .and. expn2 .ne. 0.) rval = rval                  &
     & + r0*newden*(1.-rhox**expn2)**expn1
      if(rhox .gt. 1.0) rval = r0*fracn0
      rpval = (-r0*(1.-fracn0-newden)*alphar*betar*              &
     &  (rhox**(betar-1.))*                                             &
     &  ((1.-rhox**betar)**(alphar-1.)))/flfactor
      if(nflag .gt. 0 .and. expn2 .ne. 0.) rpval = rpval                &
     & + ((-r0*newden)*expn1*expn2*(rhox**(expn2-1._R8))*               &
     & (1.-rhox**expn2)**(expn1-1._R8))/flfactor
!     if(rhox .le. 0.925) then
!     rval = r0
!     rpval = 0.
!     else
!     rval = r0 + ((r0*fracn0-r0)/(1.0-0.925))*(rhox-0.925)
!     if(rhox .ge. 1.0) rval = r0*fracn0
!     rpval = ((r0*fracn0-r0)/(1.0-0.925))/flfactor
!     endif
      endif

  case default
      stop "error in reval : invalid iden type"

  end select

      if(iimp.le.1) return
!
      jv = 2
      frac = 1.0_R8
      if(psv.gt.xsv(2)) then
      do 15 j=3,npsit+1
      frac = (psv-xsv(j-1))/(xsv(j)-xsv(j-1))
      jv = j
      if(frac.ge.0.0_R8.and. frac.le.1.0_R8) go to 21
   15 continue
      frac = 1.0_R8
                        endif
   21 rval = frac*aneimp(jv)+(1._R8-frac)*aneimp(jv-1) + rval
!
      return
   10 continue
!
!.....get r from interpolating polynomial for surface averaged transport
      psve = max(psv,xsv(2))
      lval = 2
      do 20 lv=3,npsit
      lval = lv
      if(xsv(lval).gt.psve) go to 30
   20 continue
      go to 40
   30 continue
      lm = lval-1
      rval=asv(lm,4)+psve*(bsv(lm,4)+psve*(csv(lm,4)+psve*dsv(lm,4)))
      if(rval.lt. fracn0*r0) rval = fracn0*r0
      rpval = bsv(lm,4)+psve*(2._R8*csv(lm,4)+3._R8*dsv(lm,4)*psve)
      return
  110 continue
      psl =  (psv-psimin)*delpsi
!
!
!.....get r from interpolating polynomial for ifunc=6
      lval = 2
      do 120 lv=3,npsit6+1
      lval = lv
      if(xs6(lval).gt.psl) go to 130
  120 continue
      go to 40
  130 continue
      lm = lval-1
      rval=as6(lm,4)+psl*(bs6(lm,4)+psl*(cs6(lm,4)+psl*ds6(lm,4)))
      if(rval.lt. fracn0*r0) rval = fracn0*r0
      rpval = bs6(lm,4)+psl*(2._R8*cs6(lm,4)+3._R8*ds6(lm,4)*psl)
      return
   40 continue
!
!
!.....vacuum region
      rval = fracn0*r0
      rpval = 0._R8
      return
  100 continue
!.....error exit
      ineg=55
      rval = 0._R8
      rpval = 0._R8
      write(nout,1000) psv,psimin,psilim
 1000 format(" error in reval,psv,psimin,psilim=",1p3e12.4)
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

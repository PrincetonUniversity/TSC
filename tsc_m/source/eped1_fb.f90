      subroutine eped1_fb(t_eped)
!
      USE CLINAM
      USE SAPROP
!     USE COMWOY
!     USE SVDCOM
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER l,ll,llsav,li,ii,ngr,iabs,ngrvwl,ngrvcl,ig,iigr,igr,i
      INTEGER icnt8,icnt9,iw,iswtch,n,nl,ie,np,iem10,iem20,iii
      INTEGER index, lsav
!============
!============      
!
      REAL*8, DIMENSION(20) :: eped_cur,eped_den,eped_tem
      REAL*8 ne_ped,te_ped,te_ped1,te_ped2,t_eped,cur_t0,tped_corr,AREAL
      REAL*8 eped_den_min,eped_den_max,eped_cur_min,eped_cur_max,tint,zip
      REAL*8, DIMENSION(5) :: te_dum
      INTEGER j,j1,j2,npsihm,jd1,jd2,jc1,jc2


!     hardwire values from EPED1 calculations      
      eped_cur(1:5) = 7.0E+06
      eped_cur(6:10) = 8.0E+06
      eped_cur(11:15) = 9.0E+06
      eped_cur(16:20) = 1.0E+07
      eped_den(1) = 4.00E+19
      eped_den(6) = 4.00E+19
      eped_den(11) = 4.00E+19
      eped_den(16) = 4.00E+19
      eped_den(2) = 5.50E+19
      eped_den(7) = 5.50E+19
      eped_den(12) = 5.50E+19
      eped_den(17) = 5.50E+19
      eped_den(3) = 7.00E+19
      eped_den(8) = 7.00E+19
      eped_den(13) = 7.00E+19
      eped_den(18) = 7.00E+19
      eped_den(4) = 8.50E+19
      eped_den(9) = 8.50E+19
      eped_den(14) = 8.50E+19
      eped_den(19) = 8.50E+19
      eped_den(5) = 1.00E+20
      eped_den(10) = 1.00E+20
      eped_den(15) = 1.00E+20
      eped_den(20) = 1.00E+20
      eped_tem(1) = 4.02E+03
      eped_tem(2) = 3.42E+03
      eped_tem(3) = 3.08E+03
      eped_tem(4) = 3.01E+03
      eped_tem(5) = 2.84E+03
      eped_tem(6) = 4.35E+03
      eped_tem(7) = 3.74E+03
      eped_tem(8) = 3.27E+03
      eped_tem(9) = 3.14E+03
      eped_tem(10) = 2.98E+03
      eped_tem(11) = 4.43E+03
      eped_tem(12) = 4.02E+03
      eped_tem(13) = 3.54E+03
      eped_tem(14) = 3.32E+03
      eped_tem(15) = 3.17E+03
      eped_tem(16) = 5.39E+03
      eped_tem(17) = 4.09E+03
      eped_tem(18) = 3.91E+03
      eped_tem(19) = 3.34E+03
      eped_tem(20) = 3.3E+03


      eped_den_min = 4.0E+19
      eped_den_max = 1.0E+20
      eped_cur_min = 7.0E+06
      eped_cur_max = 1.0E+07
!-------------------------------------------
!	interpolate 
! needed: total current at time t
!         pedestal density
!         a3003
!	  pedestal location (for pedestal density)

!  values are tabulated for current 7-10MA, nped=4-10e19
!  for lower currents we need to extrapolate pedestal values otherwise
!  the pedestal height will be too large - for example in the ramp-up phase
!  this will be changed as soon as a more complete table of values is available             
      
      npsihm = int(pwidthc * AREAL(npsit))
      ne_ped = ane(npsihm)
      cur_t0 = apld0*udsi
      zip = tcurdtp*tpi*udsi*1.E-6_R8
      cur_t0 = zip*1.e6
      te_ped = 0.5*(te(npsihm)+ti(npsihm))
!      if (ne_ped .lt. eped_den_min) ne_ped = eped_den_min
!      if (ne_ped .gt. eped_den_max) ne_ped = eped_den_max
!      if (cur_t0 .lt. eped_cur_min) cur_t0 = eped_cur_min
!      if (cur_t0 .gt. eped_cur_max) cur_t0 = eped_cur_max
!      if (ne_ped .ge. eped_den_min .and. ne_ped .le. eped_den_max) then 
      
      if(cur_t0 .le. eped_cur_min) j1=1
      if(cur_t0 .gt. eped_cur_max) j1=16
      if(cur_t0 .le. eped_cur_min) then 
         te_dum(1)=eped_tem(1)-1.e-6*(eped_tem(6)-eped_tem(1))*(eped_cur_min-cur_t0)
         te_dum(2)=eped_tem(2)-1.e-6*(eped_tem(7)-eped_tem(2))*(eped_cur_min-cur_t0)
         te_dum(3)=eped_tem(3)-1.e-6*(eped_tem(8)-eped_tem(3))*(eped_cur_min-cur_t0)
         te_dum(4)=eped_tem(4)-1.e-6*(eped_tem(9)-eped_tem(4))*(eped_cur_min-cur_t0)
         te_dum(5)=eped_tem(5)-1.e-6*(eped_tem(10)-eped_tem(5))*(eped_cur_min-cur_t0)
         do i=j1,j1+3
           if(ne_ped .gt. eped_den(i) .and. ne_ped .le. eped_den(i+1)) then
             jd1 = i
             jd2 = i+1
             tint = (ne_ped-eped_den(jd1))/(eped_den(jd2)-eped_den(jd1))
             exit
           end if
         end do
         if(ne_ped .le. eped_den(j1)) then
            jd1 = j1
            jd2 = j1+1
            tint = 0.0
         end if
         if(ne_ped .gt. eped_den(j1+4)) then
            jd1 = j1+3
            jd2 = j1+4
            tint = 1.0
         end if
         t_eped = te_dum(jd1) + (te_dum(jd2)-te_dum(jd1))*tint
      endif 
      if(cur_t0 .gt. eped_cur_max) then 
         te_dum(1)=eped_tem(16)+1.e-6*(eped_tem(16)-eped_tem(11))*(cur_t0-eped_cur_max)
         te_dum(2)=eped_tem(17)-1.e-6*(eped_tem(17)-eped_tem(12))*(cur_t0-eped_cur_max)
         te_dum(3)=eped_tem(18)-1.e-6*(eped_tem(18)-eped_tem(13))*(cur_t0-eped_cur_max)
         te_dum(4)=eped_tem(19)-1.e-6*(eped_tem(19)-eped_tem(14))*(cur_t0-eped_cur_max)
         te_dum(5)=eped_tem(20)-1.e-6*(eped_tem(20)-eped_tem(15))*(cur_t0-eped_cur_max)
         do i=j1,j1+3
           if(ne_ped .gt. eped_den(i) .and. ne_ped .le. eped_den(i+1)) then
             jd1 = i-15
             jd2 = i-14
             tint = (ne_ped-eped_den(i))/(eped_den(i+1)-eped_den(i))
             exit
           end if
         end do
         if(ne_ped .le. eped_den(j1)) then
            jd1 = 1
            jd2 = 2
            tint = 0.0
         end if
         if(ne_ped .gt. eped_den(j1+4)) then
            jd1 = 4
            jd2 = 5
            tint = 1.0
         end if
         t_eped = te_dum(jd1) + (te_dum(jd2)-te_dum(jd1))*tint
      end if
      if(cur_t0 .gt. eped_cur_min .and. cur_t0 .le. eped_cur_max) then 
         do j=1,19
            if(cur_t0 .gt. eped_cur(j) .and. cur_t0 .le. eped_cur(j+1)) then
               j1=j-4
               j2=j+1
            end if
         end do
         do i=j1,j1+3
           if(ne_ped .gt. eped_den(i) .and. ne_ped .le. eped_den(i+1)) then
             jd1 = i
             jd2 = i+1
             tint = (ne_ped-eped_den(jd1))/(eped_den(jd2)-eped_den(jd1))
             exit
           end if
         end do
         if(ne_ped .le. eped_den(j1)) then
           jd1 = j1
           jd2 = j1+1
           tint = 0.0
         end if
         if(ne_ped .gt. eped_den(j1+4)) then
           jd1 = j1+3
           jd2 = j1+4
           tint = 1.0
         end if
         te_ped1 = eped_tem(jd1) + (eped_tem(jd2)-eped_tem(jd1))*tint
         do i=j2,j2+3
           if(ne_ped .gt. eped_den(i) .and. ne_ped .le. eped_den(i+1)) then
             jd1 = i
             jd2 = i+1
             tint = (ne_ped-eped_den(jd1))/(eped_den(jd2)-eped_den(jd1))
             exit
           end if
         end do
         if(ne_ped .le. eped_den(j1)) then
           jd1 = j2
           jd2 = j2+1
           tint = 0.0
         end if
         if(ne_ped .gt. eped_den(j1+4)) then
           jd1 = j2+3
           jd2 = j2+4
           tint = 1.0
         end if
         te_ped2 = eped_tem(jd1) + (eped_tem(jd2)-eped_tem(jd1))*tint
         tint = (cur_t0-eped_cur(j1))/(eped_cur(j2)-eped_cur(j1))
         t_eped = te_ped1 + (te_ped2-te_ped1)*tint
      end if 
!     write(*,*) 'Tped = ',t_eped
!     tped_corr = te_ped/t_eped

      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

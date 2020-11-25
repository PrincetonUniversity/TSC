      subroutine advanp
!......2.30 advanp
!
!
!......calculates time derivative term for point centered variables
!          b = (j/x**2)*delstar(a) . . .  curl of velocity
!          psi = poloidal flux
!
      USE CLINAM
      USE SCR1
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER i,j
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 xa,x5,rhxsq5,rhxsqa,ajmid,denom,th3mmax
      REAL*8 th3m,th3,befo,t28,t4,t6,t5,bt,fac3,fac4,gis,gjs,rjf,oi
      REAL*8 oj
!============
!     dimension aboune(pnz),abounf(pnz)
!============      
      INTEGER :: istat = 0 
      REAL*8, ALLOCATABLE, DIMENSION(:) :: aboune
      REAL*8, ALLOCATABLE, DIMENSION(:) :: abounf
!============      
      IF(.not.ALLOCATED(aboune)) ALLOCATE( aboune(pnz), STAT=istat)
      IF(.not.ALLOCATED(abounf)) ALLOCATE( abounf(pnz), STAT=istat)
!============      
      if (istat .ne. 0) stop 'Allocation Error : advanp  ' 
!============      
!
      i = icol
      do 7 j=2,nzp
      abounf(j) = 1
      if(iexv(i  ,j) .eq.1) abounf(j)=0.0_R8
      aboune(j) = 1.0_R8
      if(iexv(i+1,j) .eq. 1) aboune(j)=0.0_R8
      if(iexv(i-1,j) .eq. 1) aboune(j)=0.0_R8
      if(iexv(i,j+1) .eq. 1) aboune(j)=0.0_R8
      if(iexv(i,j-1) .eq. 1) aboune(j)=0.0_R8
    7 if(iexv(i ,j ) .eq. 1) aboune(j)=0.0_R8
!
!
!
      if(i.eq.nxp) go to 60
      if(i.eq.2 .and. lrswtch .gt. 0) go to 60
      if(lrswtch.gt.0) go to 32
!
!
!
!.....define left fluxes
      do 8 j=3-isym,nz
      bf2pa(j) = bf1pa(j)
    8 continue
!
      xa = xarh(i+1)
      x5 = xary(i)
      rhxsq5 = .5_R8/x5**2
      rhxsqa = .5_R8/xa**2
!
      do 9 j=2,nz
!
!.....calculate the top fluxes
!
      ajmid = (ajphi(i,j)+ajphi(i,j+1))*x5
      bf1ta(j) = (-ajmid*.25_R8*(pip(j+1)+pip(j)-pim(j+1)-pim(j))        &  
     & - (gip(j+1)*xsqoj(i+1))**2+(giz(j+1)*xsqoj(i))**2)*rhxsq5         &  
     & *.5_R8*(rdens(i,j)+rdens(i,j+1))
!
!.....calculate the right fluxes
!
      ajmid = (ajphi(i,j)+ajphi(i+1,j))*xa
      bf1pa(j) = (ajmid*.25_R8*(piz(j+1)+pip(j+1)-piz(j-1)-pip(j-1))     &  
     & + (gip(j+1)*xsqoj(i+1))**2-(gip(j)*xsqoj(i+1))**2)*rhxsqa         &  
     & *.5_R8*(rdens(i,j)+rdens(i+1,j))
!
    9 continue
      bf1ta(1) = bf1ta(2)
!...................................................................
!
!.....advance b   (leap frog method)
!
!.......................................................................
!
      if(i.eq.2) go to 60
      denom = (8._R8*(dt+dtold)*acoef(10)*amux)
      th3mmax = 0.5_R8
      if(denom.gt.0.0_R8)                                                &  
     &th3mmax = (deex**2+deez**2)/denom
      th3m = min(0.5_R8,th3mmax)
!
!....temp fix 6/28/98
      th3m = 0._R8
      th3 = 1._R8- th3m
!
      befo = acoef(10)*amux
      t28= befo*dzisq
      t4 = befo*dxisq*xary(i-1)/xarh(i)
      t6 = befo*dxisq*xary(i+1)/xarh(i+1)
      t5 = th3*befo*(2._R8*dzisq + dxisq*xary(i)*(1._R8/xarh(i)+1._R8/   &  
     & xarh(i+1)))
!
      do 20 j=3-isym,nz
      bt = (bf1pa(j) - bf2pa(j) + bf1ta(j) - bf1ta(j-1))
!
      denom = 1._R8+ .5_R8*(dt+dtold)*(t5 + acoef(92) )
      fac3 = (2._R8-denom)*abounf(j)/denom
      fac4 = (dt+dtold)*abounf(j)/denom
!
      b(i,j) = fac3*bo(i,j) + fac4*(bt                                   &  
     & +t28*(th3*biz(j+1)+th3m*(bo2(i,j+1)-bo2(i,j))                     &  
     &      +th3*biz(j-1)+th3m*(bo2(i,j-1)-bo2(i,j)))                    &  
     & + t4*(th3*bim(j)  +th3m*(bo2(i-1,j)-xary(i)/xary(i-1)*bo2(i,j)))  &  
     & + t6*(th3*bip(j)  +th3m*(bo2(i+1,j)-xary(i)/xary(i+1)*bo2(i,j))))     
!
      bo(i,j) = biz(j)
   20 continue
!
      gis = .25_R8*dxisq
      gjs = .25_R8*dzisq
      rjf = .25_R8/(xary(i)*dxdz)
      do 30 j=3-isym,nz
 
      oi = (omeg(i+1,j+1)+omeg(i+1,j)-omeg(i,j+1)-omeg(i,j))
      oj = (omeg(i+1,j+1)+omeg(i,j+1)-omeg(i,j)-omeg(i+1,j))
      ptsave(i,j) = -((pip(j)-pim(j))*((aiz(j+1)-aiz(j-1))*rjf+oi*gis)   &  
     &               -(piz(j+1)-piz(j-1))*((abip(j)-aim(j))*rjf-oj*gjs))     
!
!     vx = .5*vecx(i,j) + .125*(vecx(i-1,j) + vecx(i+1,j)
!    1                        + vecx(i,j+1) + vecx(i,j-1))
!     vz = .5*vecz(i,j) + .125*(vecz(i-1,j) + vecz(i+1,j)
!    1                        + vecz(i,j+1) + vecz(i,j-1))
!     xpt = xary(i) - (dt+dtold)*vx
!     zpt = zary(j) - (dt+dtold)*vz
!     call grap3(zpt,xpt,psistar)
!     if(dt+dtold .gt. 0)
!    1ptsave(i,j) =(psistar - pso2(i,j))/(dt+dtold)
   30 continue
   32 continue
!
      b(i,nzp) = 0._R8
      ptsave(i,nzp) = reboun
      if(isym.ne.0) return
      b(i,2) = 0._R8
      ptsave(i,2) = reboun
      return

   60 continue
      do 70 j=2,nzp
      b(i,j) = 0._R8
      ptsave(i,j) = reboun
   70 continue
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

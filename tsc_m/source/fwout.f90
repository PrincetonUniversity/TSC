       subroutine  fwout
!                                Write forces on wires to l.u. nout
!                                Revised 18 Feb 88 to condense output
      USE CLINAM
      USE FVV1
      USE SAPROP
      USE SCR3
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER ipr,ii,n,i,j,k
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 fcpr,t1,dpsir,dpsiz
!============
        dimension  fcpr(6,2)
!============      
!============      
!       dimension  fcpr(6,2), fcpr1(6), fcpr2(6)
!       equivalence (fcpr(1,1), fcpr1(1)), (fcpr(1,2), fcpr2(1))
!
      write(nout,10)
   10 format(/'                                                          &  
     & ________________________________________________________          &  
     &________________________________________________________________'/  &  
     &                                                                   &  
     & '   nw  cwire(a)     psi     dpsir    dpsiz  fr(kNt)  fz(kNt)'4x,  &  
     &                                                                   &  
     & '   nw  cwire(a)     psi     dpsir    dpsiz  fr(kNt)  fz(kNt)')
      write(nout,15)
   15 format( '                                                          &  
     & ________________________________________________________          &  
     &____     _______________________________________________________')     
!
        ipr = 0
      do 85 ii=nc0,nwire
      n = ncoil-nwire+ii
        i = iwire(ii)
        j = jwire(ii)
        t1 = ccoil(n)*udsi
!
        dpsir = 0.2_R8* (psi(i+1,j) - psi(i-1,j))                        &  
     &        + 0.4_R8* (psi(i+2,j) - psi(i-2,j))
!
        if (j.lt.3)   then
                      dpsiz = psi(i,j+1) - psi(i,j-1)
        else
        dpsiz = 0.2_R8* (psi(i,j+1) - psi(i,j-1))                        &  
     &        + 0.4_R8* (psi(i,j+2) - psi(i,j-2))
                      endif
!****        bpolr(ii) = + 0.5* dpsiz / (deez * xwire(ii))
!****        bpolz(ii) = - 0.5* dpsir / (deex * xwire(ii))
        fwirer(ii) = - 0.5_R8* t1 * dpsir / deex
        fwirez(ii) = - 0.5_R8* t1 * dpsiz / deez
!****      bpol = sqrt(bpolr(ii)**2 + bpolz(ii)**2)
!****      fwire= sqrt(fwirer(ii)**2 + fwirez(ii)**2)
!****      write(nout,5520) ii,i,j,xwire(ii),zwire(ii),t1,psi(i,j),dpsir
!****     &         , dpsiz, bpolr(ii), bpolz(ii), bpol, fwirer(ii),
!****     &        fwirez(ii), fwire
!
        ipr = ipr + 1
        fcpr(1,ipr) = t1
        fcpr(2,ipr) = psi(i,j)
        fcpr(3,ipr) = dpsir
        fcpr(4,ipr) = dpsiz
        fcpr(5,ipr) = fwirer(ii) * tpi * 0.001_R8
        fcpr(6,ipr) = fwirez(ii) * tpi * 0.001_R8
        if (ipr.eq.2)   then
                        ipr = 0
                        write (nout,20)   ii-1, (fcpr(k,1),k=1,6),       &  
     &                                    ii,   (fcpr(k,2),k=1,6)
   20 format(1x, 2(i4, f9.1, f9.4, 1x, 2f9.4, 2x, 2f8.1, 5x))
                        endif
   85 continue
      write(nout,15)
        write (nout,110)   vvdum
  110   format(' betapol    li/2     q95   xcurf   zcurf                 &  
     & ceg1                                                              &  
     &    ceg2    ceg3    ceg4    ceg5    ceg6    ceg7    ceg8    ceg9   &  
     & ceg10' / 5f8.3, 11f8.1)
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

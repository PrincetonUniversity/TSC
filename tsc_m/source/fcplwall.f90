      subroutine  fcplwall (iwopen)
!*************************************************************************
!
!           Calculate poloidal current between plasma & wall
!
      USE CLINAM
      USE FVV1
      USE SCR15
      USE SCR16
      USE SCR21
      USE SCR22
      USE SCR23
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER iwopen,m,n,iw,jw,msegmax
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 etawirev,etaplasv
!============
        if (iwopen.gt.0)   go to 150
      write (nout, 105)
      m = 0
      do 100 n=1,nwire
      iw = iwire(n)
      jw = jwire(n)
!              Exclude outer shell, port, & belly band of SSAT vessel
!              Exclude connecting groups 55 and 56
!ccccccc    if (igroupw(n).eq.19 .or. igroupw(n).eq.18 .or.igroupw(n).eq.50
!ccccccc     &    .or.igroupw(n).eq.55 .or. igroupw(n).eq.56)   go to 100
!                       Exclude top part of inner shell
!ccccccc    if (igroupw(n).eq.20 .and. abs(xwire(n)).gt.1.60)   go to 100
      etawirev = etay(iw,jw)
      if (itap1(iw-1,jw).ne.1) then
               m = m + 1
               kgrouplw(m) = igroupw(n)
               jplexv(m) = iexv(iw-1,jw)
               xplw1(m) = xary(iw-1)
               zplw1(m) = zary(jw)
               xplw2(m) = xary(iw)
               zplw2(m) = zary(jw)
               etaplasv = etay(iw-1,jw)
         write (nout, 110)   m, kgrouplw(m), xplw1(m),                   &  
     &       xplw2(m), zplw1(m), zplw2(m), jplexv(m), etawirev,etaplasv
               endif
      if (itap1(iw+1,jw).ne.1) then
               m = m + 1
               jplexv(m) = iexv(iw+1,jw)
               kgrouplw(m) = igroupw(n)
               xplw1(m) = xary(iw)
               zplw1(m) = zary(jw)
               xplw2(m) = xary(iw+1)
               zplw2(m) = zary(jw)
               etaplasv = etay(iw+1,jw)
         write (nout, 110)   m, kgrouplw(m), xplw1(m),                   &  
     &       xplw2(m), zplw1(m), zplw2(m), jplexv(m), etawirev,etaplasv
               endif
!
      if (itap1(iw,jw-1).ne.1) then
               m = m+1
               jplexv(m) = iexv(iw,jw-1)
               kgrouplw(m) = igroupw(n)
               xplw1(m) = xary(iw)
               zplw1(m) = zary(jw-1)
               xplw2(m) = xary(iw)
               zplw2(m) = zary(jw)
               etaplasv = etay(iw,jw-1)
         write (nout, 110)   m, kgrouplw(m), xplw1(m),                   &  
     &       xplw2(m), zplw1(m), zplw2(m), jplexv(m), etawirev,etaplasv
               endif
      if (itap1(iw,jw+1).ne.1) then
               m = m+1
               jplexv(m) = iexv(iw,jw+1)
               kgrouplw(m) = igroupw(n)
               xplw1(m) = xary(iw)
               zplw1(m) = zary(jw)
               xplw2(m) = xary(iw)
               zplw2(m) = zary(jw+1)
               etaplasv = etay(iw,jw+1)
         write (nout, 110)   m, kgrouplw(m), xplw1(m),                   &  
     &       xplw2(m), zplw1(m), zplw2(m), jplexv(m), etawirev,etaplasv
               endif
  100    continue
      mplmax = m
      write (nout, 105)
  105    format (' FCPLWALL:',                                           &  
     &        /'    m  grp          x1        x2          z1        z2   &  
     &jplexv   etawirev   etaplasv')
  110    format (2i5, 2x, 2f10.3, 2x, 2f10.3, i8, 2x, 1p2e11.3)
      return
!----------------------------------------------------------------------
  150    continue
!ccccccc    write (nout, 196)  times, udsr
!ccccccc  196  format (' FCPLWALL:     time =', 1pe11.3, '  udsr =', 1pe11.3,
!ccccccc     &        /'    m  grp  polcur(A)  polvolt(V)   etawire   etaplas
!ccccccc     &       x1        x2        z1        z2  jplexv')
      m = 0
      do 200 n=1,nwire
      iw = iwire(n)
      jw = jwire(n)
!              Exclude outer shell, port, & belly band of SSAT vessel
!              Exclude connecting groups 55 and 56
!ccccccc    if (igroupw(n).eq.19 .or. igroupw(n).eq.18 .or.igroupw(n).eq.50
!ccccccc     &    .or.igroupw(n).eq.55 .or. igroupw(n).eq.56)   go to 200
!                       Exclude top & right part of inner shell
!ccccccc    if (igroupw(n).eq.20 .and. abs(xwire(n)).gt.1.60)   go to 200
      if (itap1(iw-1,jw).ne.1) then
               m = m + 1
         polsegpl(m) = -tpi*udsi * ( g(iw,jw+1) * xsqoj(iw)              &  
     &                       -g(iw,jw  ) * xsqoj(iw))
      vpolsegpl(m) = -polsegpl(m) * (etay(iw,jw) + etay(iw-1,jw)) *      &  
     & 0.5_R8                                                            &  
     &                * deex * udsr / (tpi*xary(iw)*deez)
!ccccccc    write (nout,198) m, kgrouplw(m), polsegpl(m), vpolsegpl(m),
!ccccccc     &                etay(iw,jw), etay(iw-1,jw),
!ccccccc     &      xplw1(m), xplw2(m), zplw1(m), zplw2(m), jplexv(m)
  198    format (1x,2i4, 1x, 1p4e11.3, 2x, 0p2f9.3, 2x, 2f9.3, i4)
               endif
      if (itap1(iw+1,jw).ne.1) then
               m = m + 1
         polsegpl(m) = -tpi*udsi * ( g(iw+1,jw+1) * xsqoj(iw+1)          &  
     &                     -g(iw+1,jw  ) * xsqoj(iw+1))
      vpolsegpl(m) = -polsegpl(m) * (etay(iw,jw) + etay(iw+1,jw)) *      &  
     & 0.5_R8                                                            &  
     &                * deex * udsr / (tpi*xary(iw)*deez)
!ccccccc    write (nout,198) m, kgrouplw(m), polsegpl(m), vpolsegpl(m),
!ccccccc     &                etay(iw,jw), etay(iw+1,jw),
!ccccccc     &      xplw1(m), xplw2(m), zplw1(m), zplw2(m), jplexv(m)
               endif
!
      if (itap1(iw,jw-1).ne.1) then
               m = m+1
         polsegpl(m) =  tpi*udsi * ( g(iw+1,jw) * xsqoj(iw+1)            &  
     &                       -g(iw  ,jw) * xsqoj(iw)  )
      vpolsegpl(m) = -polsegpl(m) * (etay(iw,jw) + etay(iw,jw-1)) *      &  
     & 0.5_R8                                                            &  
     &                * deez * udsr / (tpi*xary(iw)*deex)
!ccccccc    write (nout,198) m, kgrouplw(m), polsegpl(m), vpolsegpl(m),
!ccccccc     &                etay(iw,jw), etay(iw,jw-1),
!ccccccc     &      xplw1(m), xplw2(m), zplw1(m), zplw2(m), jplexv(m)
               endif
      if (itap1(iw,jw+1).ne.1) then
               m = m+1
         polsegpl(m) =  tpi*udsi * ( g(iw+1,jw+1) * xsqoj(iw+1)          &  
     &                       -g(iw  ,jw+1) * xsqoj(iw)  )
      vpolsegpl(m) = -polsegpl(m) * (etay(iw,jw) + etay(iw,jw+1)) *      &  
     & 0.5_R8                                                            &  
     &                * deez * udsr / (tpi*xary(iw)*deex)
!ccccccc    write (nout,198) m, kgrouplw(m), polsegpl(m), vpolsegpl(m),
!ccccccc     &                etay(iw,jw), etay(iw,jw+1),
!ccccccc     &      xplw1(m), xplw2(m), zplw1(m), zplw2(m), jplexv(m)
               endif
  200    continue
      mplmax = m
!
      plsegmax = 0.0_R8
      kgrplmax = 0
      msegmax = 0
      vpolsegmax = 0.0_R8
      do 300  m=1,mplmax
      if (abs(polsegpl(m)) .gt. abs(plsegmax))  then
                     plsegmax = polsegpl(m)
                     vpolsegmax = vpolsegpl(m)
                     kgrplmax = kgrouplw(m)
                     msegmax = m
                     endif
  300    continue
      write (nterm,214)  msegmax, kgrplmax, plsegmax, vpolsegmax
  214 format('FCPLWALL: m, grp, plsegmax, volts =', 5x,2i4, 1p2e11.3)
!ccccccc    if (times.lt.0.00019 .or. times.gt.0.00021)   return
!ccccccc    write (nout, 303)  udsr
!ccccccc  303  format (' udsr =', 1pe11.3)
!ccccccc    write (nout, 305)  times
!ccccccc  305  format (' FCPLWALL:     time =', 1pe11.3
!ccccccc     &        /'    m  grp  polcur(kA)  polvolt(kV)         x1        x2
!ccccccc     &        z1        z2  jplexv')
!ccccccc    write (nout, 310)   (m, kgrouplw(m), 1.e-3*polsegpl(m),
!ccccccc     &                     1.e-3*vpolsegpl(m),
!ccccccc     &      xplw1(m), xplw2(m), zplw1(m), zplw2(m), jplexv(m), m=1,mplmax)
!ccccccc  310  format (2i5, 2x, 1p2e11.3, 2x, 0p2f10.3, 2x, 2f10.3, i8)
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

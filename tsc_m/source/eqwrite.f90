      subroutine eqwrite
!===================================================================
!
!.....writes equilibrium functions at end of run
!
      USE CLINAM
      USE SAPROP
      USE SCR3
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER iparam,ios16,l,i,j,ios6,numfbo,jmin
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 aparam,alxo,alzo,ccono,r0o,gzeroo
!============
      dimension aparam(15)
!============      
!     dimension iparam(15),aparam(15)
!     equivalence (iparam(1),aparam(1))
!
!.....create disk file eqflou
!
      if( numargs .lt. 1 ) then
         filename = 'eqflou' // isuffix(1:1)
      else
         filename = 'eqflou' // '.' // trim(suffix)
      end if
!     ieqflou(1:6) = 'eqflou'
!     ieqflou(7:7) = isuffix(1:1)
      ieqflou = trim(filename)
      open(neqou,file=trim(ieqflou),status='unknown',form='unformatted',  &  
     &     iostat=ios16)
!
      aparam(1) = nx
      aparam(2) = nz
      aparam(3) = isym
      aparam(4) = isurf
      aparam(5) = tfmax
      aparam(6) = psib
      aparam(7) = psiob
      aparam(8) = npsit
      aparam(9) = rebouno
      aparam(10)= numfb
      aparam(11) = alx
      aparam(12) = alz
      aparam(13) = ccon
      aparam(14) = r0
      aparam(15) = gzero
!
      write(neqou) (aparam(l),l=1,15)
      write(neqou) ((psi(i,j),i=2,nxp),j=1,nzp)
      write(neqou) ((ajphi(i,j),i=2,nxp),j=1,nzp)
      write(neqou) ((q(i,j),i=2,nxp),j=2,nzp)
      write(neqou) ((g(i,j),i=2,nxp),j=2,nzp)
      write(neqou) (savfeed(l),l=1,numfb)
      write(neqou) (sumfeed(l),l=1,numfb)
      write(neqou) (psibbo(i),i=2,nxp)
      write(neqou) (psibto(i),i=2,nxp)
      write(neqou) (psiblo(j),j=2,nzp)
      write(neqou) (psibro(j),j=2,nzp)
!
      if(isurf.ne.1) go to 209
      write(neqou) (adn(l),l=1,npsi+1)
      write(neqou) (ade(l),l=1,npsi+1)
      write(neqou) (adp(l),l=1,npsi+1)
      write(neqou) (adi(l),l=1,npsi+1)
  209 continue
      close(neqou)
      return
!     ----->  read equilibrium from disk file eqflin for irst1=2  <-----
!
      entry eqread
      if( numargs .lt. 1 ) then
         filename = 'eqflin' // isuffix(1:1)
      else 
         filename = 'eqflin' // '.' // trim(suffix)
      end if
      ieqflin = trim(filename)
!     ieqflin(1:6) = 'eqflin'
!     ieqflin(7:7) = isuffix(1:1)
      open(neqin,file=trim(ieqflin),status='old',iostat=ios6,            &  
     &    form='unformatted')
!
      read(neqin) (aparam(l),l=1,15)
      tfmax = aparam(5)
      if(npsi.gt.1) dpsi = tfmax/(npsi-1)
      if(dpsi.gt.0) rdpsi = 1._R8/dpsi
      psib = aparam(6)
      psiob = aparam(7)
      npsit = int(aparam(8))
      rebouno = aparam(9)
      numfbo = int(aparam(10))
      alxo = aparam(11)
      alzo = aparam(12)
      ccono = aparam(13)
      r0o = aparam(14)
      gzeroo = aparam(15)
!
      write(nout,1209) (int(aparam(l)),l=1,4),nx,nz,isym,isurf
 1209 format(" nx,nz,isym,isurf .. old file,new file  ",8i5)
      write(nout,1210) tfmax
 1210 format(" tfmax from restart file =",1pe12.4)
      write(nout,1211) psib,psiob,npsit,rebouno
      write(nout,1212) numfbo,alxo,alzo,ccono,r0o,gzeroo
      write(nout,1213) numfb,alx,alz,ccon,r0,gzero
      if(alz .ne. alzo) ineg=1
      if(alx .ne. alxo) ineg=1
      if(ccon.ne.ccono) ineg=1
      if(r0.ne.r0o .or.gzero.ne.gzeroo) write(nout,1214)
 1214 format(" *** warning..gzero and/or r0 do not match ***")
 1211 format(" psib,psiob,npsit,rebouno  from restart file",             &  
     &       5x,1p2e12.4,i3,1pe12.4  )
 1212 format(" numfb,alx,alz,ccon,r0,gzero from restart file",           &  
     &       5x,i3,1p5e12.4 )
 1213 format(" numfb,alx,alz,ccon,r0,gzero in new start file",           &  
     &       5x,i3,1p5e12.4)
      if(int(aparam(1)) .ne. nx) ineg=1
      if(isym.eq.int(aparam(3)) .and.                                    &  
     &     nz.ne.int(aparam(2))) ineg=1
      if(int(aparam(3)) .eq. 1 .and. isym .eq. 0 .and.                   &  
     &     nz.ne.2*int(aparam(2))-1  ) ineg=1
      if(int(aparam(4)) .eq.0 .and. isurf .eq. 1) ineg=1
      if(ineg.ne.0) return
!
      jmin = nh
      if(int(aparam(3)).eq.0) jmin = 2
      read(neqin) ((psi(i,j),i=2,nxp),j=jmin-1,nzp)
      read(neqin) ((ajphi(i,j),i=2,nxp),j=jmin-1,nzp)
      read(neqin) ((q(i,j),i=2,nxp),j=jmin,nzp)
      read(neqin) ((g(i,j),i=2,nxp),j=jmin,nzp)
      read(neqin) (savfeed(l),l=1,numfbo)
      read(neqin) (sumfeed(l),l=1,numfbo)
      read(neqin) (psibbo(i),i=2,nxp)
      read(neqin) (psibto(i),i=2,nxp)
      read(neqin) (psiblo(j),j=jmin,nzp)
      read(neqin) (psibro(j),j=jmin,nzp)
      if(int(aparam(4)) .ne. 1) go to 309
      read(neqin) (adn(l),l=1,npsi+1)
      read(neqin) (ade(l),l=1,npsi+1)
      read(neqin) (adp(l),l=1,npsi+1)
      read(neqin) (adi(l),l=1,npsi+1)
  309 continue
      if(isym.ne.0 .or. int(aparam(3)) .ne. 1) return
      do 409 j=1,nh-1
      psiblo(nh-j) =-psiblo(nh+j)
      psibro(nh-j) =-psibro(nh+j)
      do 409 i=2,nxp
      psi(i,nh-j) = psi(i,nh+j)
      q(i,nh-j+1) = q(i,nh+j)
      g(i,nh-j+1) = g(i,nh+j)
  409 continue
      do 410 i=2,nxp
  410 psibbo(i) = -psibto(i)
!
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

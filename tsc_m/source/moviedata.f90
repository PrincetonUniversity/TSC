      subroutine moviedata
!
      USE ezcdf
      USE CLINAM
      USE SAPROP
      USE wallcl
      USE SCR10
      USE POLCUR
      USE SPDMOD
      USE SCR11
      IMPLICIT NONE
!
!   C. Kessel added large number of scalars and profiles to netcdf output 11/2009
!  (C. Kessel & DMC 17 Oct 2009): added EFIT-units ggprime, P, and P'
!       Jardin converted to F90 and added flux output 2/2010
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER ifirstmd,i,n,imin
      INTEGER imax,numr
      INTEGER j,ii
      INTEGER netcdfth
      LOGICAL :: ex
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 eps,ps,gval,qminim
      REAL*8 gpval,gppval,vloop,tflux,pval,ppval,t1,t3,t5,psilimw
      REAL*8 psisol,psival,xedge,eval,epval,rval,rpval,tival,qget
      REAL*8 psimaxx
      REAL*8 sumnb, sumlh, sumic, sumec, dvol, pnb, plh, pic, pec, rlin, rdis, zden, zr, za, zip, zk, zd
      REAL*8 zp, zb, zq,zmu
      REAL*8 taue98, taue89, ah98, ah89
!============
      REAL*8, DIMENSION(ppsi) :: te2, ane2,ti2,avez2,zeff2,gary, pres, presp, ggprime,etparas
      REAL*8, DIMENSION(2*pnz) :: zgridcdf
      REAL*8,  DIMENSION(pnx) :: xgridcdf, rcoord, ter, tir, aner, ajphir, qr
      REAL*8,  DIMENSION(2*pnx*pnz) ::packcdf
      REAL*8, DIMENSION(pnx,pnz) :: psiaux
      INTEGER, DIMENSION(pnx,pnz) :: ipackcdf
      REAL*8, DIMENSION(ppsi) :: chie, chii, anz1, anz2, anz3,           &
              anz4,anz5,anz6,anz7,anz8,                                  &
              ajbs,ajnb,ajlh,ajic,ajec,sab,                              &
              henb,hinb,heic,hiic,helh,                                  &
              heec,hiec,heal,hial,hbre,hcyc,                             &
              vl,tflx,ajto,hilh,hlin
      REAL*8 gcur01,gcur02,gcur03,gcur04,gcur05,gcur06,gcur07,gcur08,gcur09,gcur10
      REAL*8 gcur11,gcur12,gcur13,gcur14,gcur15,gcur16,gcur17,gcur18,gcur19,gcur20,gcur21,gcur22
!
!     Boundary pts (Kumar)
      REAL*8, DIMENSION(2*pnparts) :: rbnd, zbnd

      integer :: ncid, status, ier, netcdfx, netcdfz, netcdft, netcdfp, nzdimcdf, netcdfr
      integer, dimension(2) ::           start2, count2, oneddims
      integer, dimension(3) :: dimlens3, start3, count3, psidims 
!
!     Boundary point related (Kumar)
!
      integer, dimension(2) :: count2n
!
!     scalar variables (control and limits)
      integer :: timeid, psiminid, psilimid, psisolid, phaloid, npsitid
!
!     scalar energy variables (physics)
!     global     47        48      46        119    120    121
      integer :: pohmicid, pauxid, palphaid, breid, cycid, impid,  &
                 pnbid, plhid, picid, pecid
!
!     scalar current variables (physics)
!     global     8      141    142    146     198    213    242
      integer :: aplid, bsiid, cdiid, alhiid, fwiid, tniid, eciid
!
!     scalar shape variables (geometry)
!     global      34      35       36      37
      integer ::  rmajorid,rminorid,deltaid,kappaid,xmagid,zmagid, &
                  rq1id
!
!
!.....other scalar physics variables 
!                13                                  135
      integer :: vloopid, xplasid, gzeroid, tfluxid, zeffid, qzeroid,  &
        qminid, q95id, qcylid, ali1id, ali3id, betatid, betapid,       &
        wthid, taue1id, taue2id, ah98id, ah89id, rzeroid, rlineid,     &
        rvolid, rorgrid, te0id, ti0id, flxstid, vstaxid,               &
        vsiaxid, vsraxid, vsipoyid, vsrpoyid, thaloid, whaloid,        &
        tpedid, psipedid,fz1id, fz2id, fz3id, fz4id, fz5id, fz6id,     &
        fz7id, fz8id, injcurrid, injcur2id, helicityid,                &
        gcur01id, gcur02id, gcur03id, gcur04id, gcur05id, gcur06id,    &
        gcur07id, gcur08id, gcur09id, gcur10id,                        &
        gcur11id, gcur12id, gcur13id, gcur14id, gcur15id, gcur16id,    &
        gcur17id, gcur18id, gcur19id, gcur20id, gcur21id, gcur22id
!
!  (DMC 17 Oct 2009): added EFIT-units ggprime, P, and P'
!
!     1D variables (vs poloidal flux)
      integer :: xsv2id, te2id, ti2id, ane2id, ajpest2id, prpest2id,   &
                 pppest2id, avez2id, zeff2id,etparaid, qprof2id, garyid,        &
                 presid, prespid, ggprimeid,                           &
                 chieid, chiiid,anhyid,anheid,anz1id,anz2id,anz3id,    &
                 anz4id,anz5id,ajtoid,ajbsid,ajnbid,ajlhid,ajicid,     &
                 ajecid,sabid,henbid,hinbid,heicid,hiicid,helhid,      &
                 hilhid,heecid,healid,hialid,hbreid,hcycid,hlinid,     &
                 vlid,tflxid,anz6id,anz7id,anz8id
!
!     1D variables (vs major radius)
      integer :: rcoordid, terid, tirid, anerid, ajphirid, qrid
!
!     1D variables (vs theta index)
!          (Kumar)
!
      integer :: rbndid,zbndid
!
!     2D variables
      integer :: psiid,divflxid,ajphiid, prid, gsid,imaskid
!
!.....variables added for idata=9 or idata=11 option
      integer :: atotid, csumid, vesppeid, vescureid, totcureid, pfluxid, efluxid, &
                 ppcurid, fbcurid, csumaid, netcdfc, netcdffl
              
!
      save
      data ifirstmd/0/
      data eps/1.e-6/
!
!
      call globe
      do i=2,nxp
        xgridcdf(i-1) = xary(i)
      enddo
      do n=1,npsi
        ps = xsv2(n)
        call geval(ps,2,gval,gpval,gppval,imag,jmag)
        call peval(ps,2,pval,ppval,imag,jmag)
        pres(n) = udsp*pval
        presp(n) = udsp*ppval
        gary(n) = gval
        ggprime(n) = gppval
      enddo
!....removed interpolations for testing:   scj 07/14/2010
!     ggprime(2) = 2.*ggprime(3)-ggprime(4)
!     presp(2) = 2.*presp(3)-presp(4)
!     ggprime(1) = 2.*ggprime(2)-ggprime(3)
!     presp(1) = 2.*presp(2)-presp(3)
      vloop = global(13)
      tflux = dpsi*(npsit-1)
!
      qprof2(1) = 2.*qprof2(2)-qprof2(3)
      qminim = min(qprof2(1),qprof2(npsit))
!
      sumnb = 0.
      sumlh = 0.
      sumic = 0.
      sumec = 0.
      do j=2,npsit
        dvol = vary(j)-vary(j-1)
        sumnb = sumnb+(savei(j)+savee(j))*dvol
        sumlh = sumlh+(savelh(j)+savilh(j))*dvol
        sumic = sumic+(savefw(j)+savifw(j))*dvol
        sumec = sumec+(savebm(j)+savibm(j))*dvol
      enddo
      pnb = sumnb*udsp/udst
      plh = sumlh*udsp/udst
      pic = sumic*udsp/udst
      pec = sumec*udsp/udst
!
      rlin = 0.
      rdis = 0.
      do i=iminn,imaxx
        if(iexv(i,nh).eq.1 .or. iexs(i,nh).eq.1) cycle
        psival = psi(i,nh)
        if(psival.gt.psilim) cycle
        call reval(psival,idens,isurf,rval,rpval,i,nh)
        rlin = rlin + rval*deex*udsd
        rdis = rdis + deex
      enddo
      zden = (rlin/rdis)*1.e-19
      zr = rmajor
      za = rminor
      zip = tcurdtp*tpi*udsi*1.e-6
      zk = shape5
      zd = shape3
      zp = (pohmic+paux+prad+palpha)*1.e-6
      zb = gzero/zr
      zq = 2.5*za**2*zb*(1.+zk**2)/(zip*zr)
      zmu = 2.5
      taue98 = 0.0562*zip**0.93*zb**0.15*zden**0.41*zr**1.97*  &
               (za/zr)**0.58*zk**0.78*zmu**0.19/zp**0.69
      taue89 = 0.0480*zip**0.85*zb**0.2*(zden/10.)**0.1*zr**1.5* &
               (za/zr)**0.3*zk**0.5*zmu**0.5/zp**0.5
      ah98 = (global(26)/1000.)/taue98
      ah89 = (global(26)/1000.)/taue89
!
      do i=2,npsit
        chie(i) = chiesec(i)
        chii(i) = chiisec(i)
        anz1(i) = fraci(1)*ane(i)
        anz2(i) = fraci(2)*ane(i)
        anz3(i) = fraci(3)*ane(i)
        anz4(i) = fraci(4)*ane(i)
        anz5(i) = fraci(5)*ane(i)
        anz6(i) = fraci(6)*ane(i)
        anz7(i) = fraci(7)*ane(i)
        anz8(i) = fraci(8)*ane(i)
        ajto(i) = (gxmja2(i)-gxmja2(i-1))*rdpsi*udsi/tpi
        ajbs(i) = 0.5*(ajavbs(i)+ajavbs(i-1))*udsi
        ajnb(i) = 0.5*(ajavcd(i)+ajavcd(i-1))*udsi
        ajic(i) = 0.5*(ajavfw(i)+ajavfw(i-1))*udsi
        ajlh(i) = 0.5*(ajavlh(i)+ajavlh(i-1))*udsi
        ajec(i) = 0.5*(ajavec(i)+ajavec(i-1))*udsi
        sab(i) = gzero/rmajora(i)
        henb(i) = savee(i)*udsp/udst
        hinb(i) = savei(i)*udsp/udst
        heic(i) = savefw(i)*udsp/udst
        hiic(i) = savifw(i)*udsp/udst
        helh(i) = savelh(i)*udsp/udst
        hilh(i) = savilh(i)*udsp/udst
        heec(i) = savebm(i)*udsp/udst
	hiec(i) = savibm(i)*udsp/udst
        heal(i) = savea(i)*udsp/udst
        hial(i) = savia(i)*udsp/udst
        hbre(i) = savebre(i)*udsp/udst
        hcyc(i) = savecyc(i)*udsp/udst
        hlin(i) = (sradion(i)+saveimp(i)*udsp/udst)
        vl(i) = 0.5*(as(i)+as(i-1))*udsv
        tflx(i) = (float(i-1)*dpsi)/(float(npsit-1)*dpsi)
      enddo
      chie(1) = 2.*chie(2)-chie(3)
      chii(1) = 2.*chii(2)-chii(3)
      anz1(1) = 2.*anz1(2)-anz1(3)
      anz2(1) = 2.*anz2(2)-anz2(3)
      anz3(1) = 2.*anz3(2)-anz3(3)
      anz4(1) = 2.*anz4(2)-anz4(3)
      anz5(1) = 2.*anz5(2)-anz5(3)
      anz6(1) = 2.*anz6(2)-anz6(3)
      anz7(1) = 2.*anz7(2)-anz7(3)
      anz8(1) = 2.*anz8(2)-anz8(3)
      ajto(1) = 2.*ajto(2)-ajto(3)
      ajbs(1) = 2.*ajbs(2)-ajbs(3)
      ajnb(1) = 2.*ajnb(2)-ajnb(3)
      ajic(1) = 2.*ajic(2)-ajic(3)
      ajlh(1) = 2.*ajlh(2)-ajlh(3)
      ajec(1) = 2.*ajec(2)-ajec(3)
      if(ajec(1) .lt. 0.) ajec(1) = 0.
      sab(1) = gzero/xmag
      henb(1) = 2.*henb(2)-henb(3)
      if(henb(1) .lt. 0.) henb(1) = 0.
      hinb(1) = 2.*hinb(2)-hinb(3)
      if(hinb(1) .lt. 0.) hinb(1) = 0.
      heic(1) = 2.*heic(2)-heic(3)
      if(heic(1) .lt. 0.) heic(1) = 0.
      hiic(1) = 2.*hiic(2)-hiic(3)
      if(hiic(1) .lt. 0.) hiic(1) = 0.
      helh(1) = 2.*helh(2)-helh(3)
      if(helh(1) .lt. 0.) helh(1) = 0.
      hilh(1) = 2.*hilh(2)-hilh(3)
      if(hilh(1) .lt. 0.) hilh(1) = 0.
      heec(1) = 2.*heec(2)-heec(3)
      if(heec(1) .lt. 0.) heec(1) = 0.
      heal(1) = 2.*heal(2)-heal(3)
      hial(1) = 2.*hial(2)-hial(3)
      hbre(1) = 2.*hbre(2)-hbre(3)
      hcyc(1) = 2.*hcyc(2)-hcyc(3)
      hlin(1) = 2.*hlin(2)-hlin(3)
      vl(1) = 2.*vl(2)-vl(3)
      tflx(1) = 0.0

      nzdimcdf = nz + isym*(nz-1)
      if(isym.eq.0) then
        do i=1,nzdimcdf
          zgridcdf(i) = zary(i+1)
        enddo
      else
        do i=3,nzp
          zgridcdf(i-2)    =   -zary(nzp+3-i)
          zgridcdf(nz-2+i) =    zary(i)
        enddo
        zgridcdf(nz) = 0
      endif
!      first 10 group currents
      call groupcur(grsum,grsum0,gvsum,gvsum0,gcur0ka,gcurfka)
      gcur01 = 0
      gcur02 = 0
      gcur03 = 0
      gcur04 = 0
      gcur05 = 0
      gcur06 = 0
      gcur07 = 0
      gcur08 = 0
      gcur09 = 0
      gcur10 = 0
      gcur11 = 0
      gcur12 = 0
      gcur13 = 0
      gcur14 = 0
      gcur15 = 0
      gcur16 = 0
      gcur17 = 0
      gcur18 = 0
      gcur19 = 0
      gcur20 = 0
      gcur21 = 0
      gcur22 = 0
      if(ngroupt.ge.1) gcur01 = grsum(1)
      if(ngroupt.ge.2) gcur02 = grsum(2)
      if(ngroupt.ge.3) gcur03 = grsum(3)
      if(ngroupt.ge.4) gcur04 = grsum(4)
      if(ngroupt.ge.5) gcur05 = grsum(5)
      if(ngroupt.ge.6) gcur06 = grsum(6)
      if(ngroupt.ge.7) gcur07 = grsum(7)
      if(ngroupt.ge.8) gcur08 = grsum(8)
      if(ngroupt.ge.9) gcur09 = grsum(9)
      if(ngroupt.ge.10) gcur10 = grsum(10)
      if(ngroupt.ge.11) gcur11 = grsum(11)
      if(ngroupt.ge.12) gcur12 = grsum(12)
      if(ngroupt.ge.13) gcur13 = grsum(13)
      if(ngroupt.ge.14) gcur14 = grsum(14)
      if(ngroupt.ge.15) gcur15 = grsum(15)
      if(ngroupt.ge.16) gcur16 = grsum(16)
      if(ngroupt.ge.17) gcur17 = grsum(17)
      if(ngroupt.ge.18) gcur18 = grsum(18)
      if(ngroupt.ge.19) gcur19 = grsum(19)
      if(ngroupt.ge.20) gcur20 = grsum(20)
      if(ngroupt.ge.21) gcur21 = grsum(21)
      if(ngroupt.ge.22) gcur22 = grsum(22)
!
      if(ifirstmd .eq. 0 ) then
!     create netcdf file...first check if start or restart run
      inquire(file='movie.cdf', exist=ex)
       if(irst1.eq.0 .or. .not. ex) then
!
       call cdf_open( ncid, 'movie.cdf', 'w' , ier)
!
          if(ier.ne.0) then
            call movieerr(1,ier)
            write(*,*) "error after cdf_open", ncid,ex,ier
          endif
!
       dimlens3 = (/nx,1,1/)
       call cdf_define(ncid,'xary',dimlens3,'R8' ,ier)
           if(ier.ne.0) call movieerr(2,ier)
       call cdf_setatt(ncid,'xary','x-coordinate of grid','m',ier)
           if(ier.ne.0) call movieerr(3,ier)

!
       dimlens3 = (/nzdimcdf,1,1/)
       call cdf_define(ncid,'zary',dimlens3,'R8' ,ier)
           if(ier.ne.0) call movieerr(5,ier)
       call cdf_setatt(ncid,'zary','z-coordinate of grid','m',ier)
           if(ier.ne.0) call movieerr(6,ier)
!
!     define time array variables using netCDF
!    0. Dimension variables
       status = nf_def_dim(ncid,'x'  ,nx          ,netcdfx)
           if(status.ne.0) call movieerr(9,status)
       status = nf_def_dim(ncid,'z'  ,nzdimcdf    ,netcdfz)
           if(status.ne.0) call movieerr(10,status)
       status = nf_def_dim(ncid,'t'  ,NF_UNLIMITED,netcdft)
           if(status.ne.0) call movieerr(11,status)
       status = nf_def_dim(ncid,'psi',ppsi        ,netcdfp)
           if(status.ne.0) call movieerr(85,status)
!
!     Following related to boundary pts. (Kumar)
!
       status =nf_def_dim(ncid,'ntheta',2*pnparts ,netcdfth)
           if(status.ne.0) call movieerr(90,status)
!
       call xlimits(imin,imax)
       numr = imax-imin+1
       status = nf_def_dim(ncid,'R  ',numr        ,netcdfr)
           if(status.ne.0) call movieerr(86,status)
!
!    1.1 Scalar variables (control and limits)     
       status = nf_def_var(ncid,'times' ,NF_DOUBLE,1,netcdft,timeid)
           if(status.ne.0) call movieerr(12,status)
       status = nf_def_var(ncid,'psimin',NF_DOUBLE,1,netcdft,psiminid)
           if(status.ne.0) call movieerr(14,status)
       status = nf_def_var(ncid,'psilim',NF_DOUBLE,1,netcdft,psilimid)
           if(status.ne.0) call movieerr(34,status)
       status = nf_def_var(ncid,'psisol',NF_DOUBLE,1,netcdft,psisolid)
           if(status.ne.0) call movieerr(34,status)
       status = nf_def_var(ncid,'phalo',NF_DOUBLE,1,netcdft,phaloid)
           if(status.ne.0) call movieerr(36,status)
       status = nf_def_var(ncid,'npsit',NF_INT,1,netcdft,npsitid)
           if(status.ne.0) call movieerr(35,status)
!
       status=nf_put_att_text(ncid,timeid,'name',16,'problem time (s)')
           if(status.ne.0) call movieerr(12,status)
       status = nf_put_att_text(ncid,psiminid,'name',21,       &
                'minimum psi in plasma')
           if(status.ne.0) call movieerr(14,status)
       status = nf_put_att_text(ncid,psilimid,'name',21,       &
                'maximum psi in plasma')
           if(status.ne.0) call movieerr(34,status)
       status = nf_put_att_text(ncid,psisolid,'name',21,       &
                'psi at edge of SOL   ')
           if(status.ne.0) call movieerr(35,status)
       status = nf_put_att_text(ncid,phaloid,'name',25,        &
                'psi value at SOL boundary')
           if(status.ne.0) call movieerr(37,status)
       status = nf_put_att_text(ncid,npsitid,'name',20,        &
                'number of flux zones')
           if(status.ne.0) call movieerr(35,status)
!   1.2 Scalar physics energy variables
       status = nf_def_var(ncid,'pohmic' ,NF_DOUBLE,1,netcdft,pohmicid)
           if(status.ne.0) call movieerr(12,status)
       status = nf_def_var(ncid,'paux'   ,NF_DOUBLE,1,netcdft,pauxid)
           if(status.ne.0) call movieerr(12,status)
       status = nf_def_var(ncid,'palpha' ,NF_DOUBLE,1,netcdft,palphaid)
           if(status.ne.0) call movieerr(12,status)
       status = nf_def_var(ncid,'bre  '  ,NF_DOUBLE,1,netcdft,breid)
           if(status.ne.0) call movieerr(12,status)
       status = nf_def_var(ncid,'cyc  '  ,NF_DOUBLE,1,netcdft,cycid)
           if(status.ne.0) call movieerr(12,status)
       status = nf_def_var(ncid,'imp  '  ,NF_DOUBLE,1,netcdft,impid)
           if(status.ne.0) call movieerr(12,status)
       status = nf_def_var(ncid,'pnb  '  ,NF_DOUBLE,1,netcdft,pnbid)
           if(status.ne.0) call movieerr(12,status)
       status = nf_def_var(ncid,'plh  '  ,NF_DOUBLE,1,netcdft,plhid)
           if(status.ne.0) call movieerr(12,status)
       status = nf_def_var(ncid,'pic  '  ,NF_DOUBLE,1,netcdft,picid)
           if(status.ne.0) call movieerr(12,status)
       status = nf_def_var(ncid,'pec  '  ,NF_DOUBLE,1,netcdft,pecid)
           if(status.ne.0) call movieerr(12,status)
!
       status = nf_put_att_text(ncid,pohmicid,'name',24,        &
                'Ohmic Heating Power (MW)')
           if(status.ne.0) call movieerr(12,status)
       status = nf_put_att_text(ncid,pauxid,'name',30,          &
                'Auxialliary Heating Power (MW)')
           if(status.ne.0) call movieerr(12,status)
       status = nf_put_att_text(ncid,palphaid,'name',16,        &
                'Alpha Power (MW)')
           if(status.ne.0) call movieerr(12,status)
       status = nf_put_att_text(ncid,breid,'name',29,           &
                'Bremsstrahlung Radiation (MW)')
           if(status.ne.0) call movieerr(12,status)
       status = nf_put_att_text(ncid,cycid,'name',24,           &
                'Cyclotron Radiation (MW)')
           if(status.ne.0) call movieerr(12,status)
       status = nf_put_att_text(ncid,impid,'name',28,           &
                'Impurity Line Radiation (MW)')
           if(status.ne.0) call movieerr(12,status)
       status = nf_put_att_text(ncid,pnbid,'name',23,           &
                'Neutral beam power (MW)')
           if(status.ne.0) call movieerr(12,status)
       status = nf_put_att_text(ncid,plhid,'name',23,           &
                'Lower hybrid power (MW)')
           if(status.ne.0) call movieerr(12,status)
       status = nf_put_att_text(ncid,picid,'name',24,           &
                'Ion cyclotron power (MW)')
           if(status.ne.0) call movieerr(12,status)
       status = nf_put_att_text(ncid,pecid,'name',24,           &
                'Electron cyclotron power (MW)')
           if(status.ne.0) call movieerr(12,status)
!
!   1.3 Scalar physics current variables
       status = nf_def_var(ncid,'apl' ,NF_DOUBLE,1,netcdft,aplid)
           if(status.ne.0) call movieerr(12,status)
       status = nf_def_var(ncid,'bsi' ,NF_DOUBLE,1,netcdft,bsiid)
           if(status.ne.0) call movieerr(12,status)
       status = nf_def_var(ncid,'cdi' ,NF_DOUBLE,1,netcdft,cdiid)
           if(status.ne.0) call movieerr(12,status)
       status = nf_def_var(ncid,'alhi',NF_DOUBLE,1,netcdft,alhiid)
           if(status.ne.0) call movieerr(12,status)
       status = nf_def_var(ncid,'fwi' ,NF_DOUBLE,1,netcdft,fwiid)
           if(status.ne.0) call movieerr(12,status)
       status = nf_def_var(ncid,'eci' ,NF_DOUBLE,1,netcdft,eciid)
           if(status.ne.0) call movieerr(12,status)
       status = nf_def_var(ncid,'tni' ,NF_DOUBLE,1,netcdft,tniid)
       status = nf_def_var(ncid,'vloop' ,NF_DOUBLE,1,netcdft,vloopid)
       status = nf_def_var(ncid,'xplas' ,NF_DOUBLE,1,netcdft,xplasid)
       status = nf_def_var(ncid,'gzero' ,NF_DOUBLE,1,netcdft,gzeroid)
       status = nf_def_var(ncid,'tflux' ,NF_DOUBLE,1,netcdft,tfluxid)
       status = nf_def_var(ncid,'rmajor' ,NF_DOUBLE,1,netcdft,rmajorid)
       status = nf_def_var(ncid,'rminor' ,NF_DOUBLE,1,netcdft,rminorid)
       status = nf_def_var(ncid,'delta' ,NF_DOUBLE,1,netcdft,deltaid)
       status = nf_def_var(ncid,'kappa' ,NF_DOUBLE,1,netcdft,kappaid)
       status = nf_def_var(ncid,'xmag' ,NF_DOUBLE,1,netcdft,xmagid)
       status = nf_def_var(ncid,'zmag' ,NF_DOUBLE,1,netcdft,zmagid)
       status = nf_def_var(ncid,'rq1' ,NF_DOUBLE,1,netcdft,rq1id)
       status = nf_def_var(ncid,'zeff' ,NF_DOUBLE,1,netcdft,zeffid)
       status = nf_def_var(ncid,'qzero' ,NF_DOUBLE,1,netcdft,qzeroid)
       status = nf_def_var(ncid,'qmin' ,NF_DOUBLE,1,netcdft,qminid)
       status = nf_def_var(ncid,'q95' ,NF_DOUBLE,1,netcdft,q95id)
       status = nf_def_var(ncid,'qcyl' ,NF_DOUBLE,1,netcdft,qcylid)
       status = nf_def_var(ncid,'ali1' ,NF_DOUBLE,1,netcdft,ali1id)
       status = nf_def_var(ncid,'ali3' ,NF_DOUBLE,1,netcdft,ali3id)
       status = nf_def_var(ncid,'betat' ,NF_DOUBLE,1,netcdft,betatid)
       status = nf_def_var(ncid,'betap' ,NF_DOUBLE,1,netcdft,betapid)
       status = nf_def_var(ncid,'wth' ,NF_DOUBLE,1,netcdft,wthid)
       status = nf_def_var(ncid,'taue1' ,NF_DOUBLE,1,netcdft,taue1id)
       status = nf_def_var(ncid,'taue2' ,NF_DOUBLE,1,netcdft,taue2id)
       status = nf_def_var(ncid,'ah98' ,NF_DOUBLE,1,netcdft,ah98id)
       status = nf_def_var(ncid,'ah89' ,NF_DOUBLE,1,netcdft,ah89id)
       status = nf_def_var(ncid,'rzero' ,NF_DOUBLE,1,netcdft,rzeroid)
       status = nf_def_var(ncid,'rline' ,NF_DOUBLE,1,netcdft,rlineid)
       status = nf_def_var(ncid,'rvol' ,NF_DOUBLE,1,netcdft,rvolid)
       status = nf_def_var(ncid,'rorgr' ,NF_DOUBLE,1,netcdft,rorgrid)
       status = nf_def_var(ncid,'te0' ,NF_DOUBLE,1,netcdft,te0id)
       status = nf_def_var(ncid,'ti0' ,NF_DOUBLE,1,netcdft,ti0id)
       status = nf_def_var(ncid,'flxst' ,NF_DOUBLE,1,netcdft,flxstid)
       status = nf_def_var(ncid,'vstax' ,NF_DOUBLE,1,netcdft,vstaxid)
       status = nf_def_var(ncid,'vsiax' ,NF_DOUBLE,1,netcdft,vsiaxid)
       status = nf_def_var(ncid,'vsrax' ,NF_DOUBLE,1,netcdft,vsraxid)
       status = nf_def_var(ncid,'vsipoy' ,NF_DOUBLE,1,netcdft,vsipoyid)
       status = nf_def_var(ncid,'vsrpoy' ,NF_DOUBLE,1,netcdft,vsrpoyid)
       status = nf_def_var(ncid,'thalo' ,NF_DOUBLE,1,netcdft,thaloid)
       status = nf_def_var(ncid,'whalo' ,NF_DOUBLE,1,netcdft,whaloid)
       status = nf_def_var(ncid,'tped' ,NF_DOUBLE,1,netcdft,tpedid)
       status = nf_def_var(ncid,'psiped' ,NF_DOUBLE,1,netcdft,psipedid)
       status = nf_def_var(ncid,'fz1' ,NF_DOUBLE,1,netcdft,fz1id)
       status = nf_def_var(ncid,'fz2' ,NF_DOUBLE,1,netcdft,fz2id)
       status = nf_def_var(ncid,'fz3' ,NF_DOUBLE,1,netcdft,fz3id)
       status = nf_def_var(ncid,'fz4' ,NF_DOUBLE,1,netcdft,fz4id)
       status = nf_def_var(ncid,'fz5' ,NF_DOUBLE,1,netcdft,fz5id)
       status = nf_def_var(ncid,'fz6' ,NF_DOUBLE,1,netcdft,fz6id)
       status = nf_def_var(ncid,'fz7' ,NF_DOUBLE,1,netcdft,fz7id)
       status = nf_def_var(ncid,'fz8' ,NF_DOUBLE,1,netcdft,fz8id)
       status = nf_def_var(ncid,'injcurr' ,NF_DOUBLE,1,netcdft,injcurrid)
       status = nf_def_var(ncid,'injcur2' ,NF_DOUBLE,1,netcdft,injcur2id)
       status = nf_def_var(ncid,'helicity' ,NF_DOUBLE,1,netcdft,helicityid)
       status = nf_def_var(ncid,'gcur01' ,NF_DOUBLE,1,netcdft,gcur01id)
       status = nf_def_var(ncid,'gcur02' ,NF_DOUBLE,1,netcdft,gcur02id)
       status = nf_def_var(ncid,'gcur03' ,NF_DOUBLE,1,netcdft,gcur03id)
       status = nf_def_var(ncid,'gcur04' ,NF_DOUBLE,1,netcdft,gcur04id)
       status = nf_def_var(ncid,'gcur05' ,NF_DOUBLE,1,netcdft,gcur05id)
       status = nf_def_var(ncid,'gcur06' ,NF_DOUBLE,1,netcdft,gcur06id)
       status = nf_def_var(ncid,'gcur07' ,NF_DOUBLE,1,netcdft,gcur07id)
       status = nf_def_var(ncid,'gcur08' ,NF_DOUBLE,1,netcdft,gcur08id)
       status = nf_def_var(ncid,'gcur09' ,NF_DOUBLE,1,netcdft,gcur09id)
       status = nf_def_var(ncid,'gcur10' ,NF_DOUBLE,1,netcdft,gcur10id)
       status = nf_def_var(ncid,'gcur11' ,NF_DOUBLE,1,netcdft,gcur11id)
       status = nf_def_var(ncid,'gcur12' ,NF_DOUBLE,1,netcdft,gcur12id)
       status = nf_def_var(ncid,'gcur13' ,NF_DOUBLE,1,netcdft,gcur13id)
       status = nf_def_var(ncid,'gcur14' ,NF_DOUBLE,1,netcdft,gcur14id)
       status = nf_def_var(ncid,'gcur15' ,NF_DOUBLE,1,netcdft,gcur15id)
       status = nf_def_var(ncid,'gcur16' ,NF_DOUBLE,1,netcdft,gcur16id)
       status = nf_def_var(ncid,'gcur17' ,NF_DOUBLE,1,netcdft,gcur17id)
       status = nf_def_var(ncid,'gcur18' ,NF_DOUBLE,1,netcdft,gcur18id)
       status = nf_def_var(ncid,'gcur19' ,NF_DOUBLE,1,netcdft,gcur19id)
       status = nf_def_var(ncid,'gcur20' ,NF_DOUBLE,1,netcdft,gcur20id)
       status = nf_def_var(ncid,'gcur21' ,NF_DOUBLE,1,netcdft,gcur21id)
       status = nf_def_var(ncid,'gcur22' ,NF_DOUBLE,1,netcdft,gcur22id)
           if(status.ne.0) call movieerr(12,status)
!
       status = nf_put_att_text(ncid,aplid,'name',18,    &
               'Plasma Current (A)')
           if(status.ne.0) call movieerr(12,status)      
       status = nf_put_att_text(ncid,bsiid,'name',21,    &
                'Bootstrap Current (A)')
           if(status.ne.0) call movieerr(12,status)
       status = nf_put_att_text(ncid,cdiid,'name',23,    &
                'Beam Driven Current (A)')
           if(status.ne.0) call movieerr(12,status)
       status = nf_put_att_text(ncid,alhiid,'name',24,   &
                'Lower Hybrid Current (A)')
           if(status.ne.0) call movieerr(12,status) 
       status = nf_put_att_text(ncid,fwiid,'name',21,    &
                'Fast Wave Current (A)')
           if(status.ne.0) call movieerr(12,status) 
       status = nf_put_att_text(ncid,eciid,'name',30,    &
                'Electron Cyclotron Current (A)')
           if(status.ne.0) call movieerr(12,status) 
       status = nf_put_att_text(ncid,tniid,'name',31,    &
                'Total Non-Inductive Current (A)')
       status = nf_put_att_text(ncid,vloopid,'name',16,  &
               'Loop Voltage (V)')
       status = nf_put_att_text(ncid,xplasid,'name',32,  &
                'Geometrical center of plasma (m)')
       status = nf_put_att_text(ncid,gzeroid,'name',35,  &
                'R times vacuum toroidal field (T-m)') 
       status = nf_put_att_text(ncid,tfluxid,'name',34,  &
                'Total Toroidal Flux in Plasma ( W)')
       status = nf_put_att_text(ncid,rmajorid,'name',26, &
                'Major Radius of Plasma (m)')
       status = nf_put_att_text(ncid,rminorid,'name',26, &
                'Minor Radius of PLasma (m)')
       status = nf_put_att_text(ncid,deltaid,'name',23,  &
                'Triangularity of Plasma')
       status = nf_put_att_text(ncid,kappaid,'name',21,  &
                'Ellipticity of Plasma')
       status = nf_put_att_text(ncid,xmagid,'name',23,   &
                'Radius of magnetic axis')
       status = nf_put_att_text(ncid,zmagid,'name',27,   &
                'Z position of magnetic axis')
       status = nf_put_att_text(ncid,rq1id,'name',21,    &
                'Radius of q=1 surface')
       status = nf_put_att_text(ncid,zeffid,'name',26,   &
                'Effective charge of Plasma')
       status = nf_put_att_text(ncid,qzeroid,'name',21,  &
                'On-axis safety factor')
       status = nf_put_att_text(ncid,qminid,'name',21,   &
                'Minimum safety factor')
       status = nf_put_att_text(ncid,q95id,'name',26,    &
               '95% pol flux safety factor')
       status = nf_put_att_text(ncid,qcylid,'name',25,   &
                'Cylindrical safety factor')
       status = nf_put_att_text(ncid,ali1id,'name',25,   &
                'Li(1) int self inductance')
       status = nf_put_att_text(ncid,ali3id,'name',25,   &
                'Li(3) int self inductance')
       status = nf_put_att_text(ncid,betatid,'name',13,  &
                'Toroidal beta')
       status = nf_put_att_text(ncid,betapid,'name',13,  &
                'Beta poloidal')
       status = nf_put_att_text(ncid,wthid,'name',26,    &
                'Thermal stored energy (MJ)')
       status = nf_put_att_text(ncid,taue1id,'name',33,  &
                'Energy confinement time w rad (s)')
       status = nf_put_att_text(ncid,taue2id,'name',34,  &
                'Energy confinement time wo rad (s)')
       status = nf_put_att_text(ncid,ah98id,'name',25,   &
                'H98(y,2) confinement mult')
       status = nf_put_att_text(ncid,ah89id,'name',20,   &
                'H89 confinement mult')
       status = nf_put_att_text(ncid,rzeroid,'name',21,  &
                'Central density (/m3)')
       status = nf_put_att_text(ncid,rlineid,'name',22,  &
                'Line ave density (/m3)')
       status = nf_put_att_text(ncid,rvolid,'name',24,   &
                'Volume ave density (/m3)')
       status = nf_put_att_text(ncid,rorgrid,'name',26,  &
                'Ratio density to Greenwald') 
       status = nf_put_att_text(ncid,te0id,'name',28,    &
                'Central electron temperature (eV)')
       status = nf_put_att_text(ncid,ti0id,'name',23,    &
                'Central ion temperature (eV)')
       status = nf_put_att_text(ncid,flxstid,'name',15,  &
                'Flux state (Wb)')
       status = nf_put_att_text(ncid,vstaxid,'name',27,  &
                'Axial total flux cons (V-s)')
       status = nf_put_att_text(ncid,vsiaxid,'name',30,  &
                'Axial internal flux cons (V-s)') 
       status = nf_put_att_text(ncid,vsraxid,'name',31,  &
                'Axial resistive flux cons (V-s)')
       status = nf_put_att_text(ncid,vsipoyid,'name',34, &
                'Poynting inductive flux cons (V-s)')
       status = nf_put_att_text(ncid,vsrpoyid,'name',34, &
                'Poynting resistive flux cons (V-s)')
       status = nf_put_att_text(ncid,thaloid,'name',21,  &
                'Halo temperature (eV)')
       status = nf_put_att_text(ncid,whaloid,'name',28,  &
                'Halo width pol flux fraction')
       status = nf_put_att_text(ncid,tpedid,'name',25,   &
                'Pedestal temperature (eV)')
       status = nf_put_att_text(ncid,psipedid,'name',35, &
                'Pedestal location pol flux fraction')
       status = nf_put_att_text(ncid,fz1id,'name',24,    &
                'Oxygen impurity fraction')
       status = nf_put_att_text(ncid,fz2id,'name',24,    &
                'Carbon impurity fraction')
       status = nf_put_att_text(ncid,fz3id,'name',22,    &
                'Iron impurity fraction')
       status = nf_put_att_text(ncid,fz4id,'name',27,    &
                'Beryllium impurity fraction')
       status = nf_put_att_text(ncid,fz5id,'name',22,    &
                'Neon impurity fraction')
       status = nf_put_att_text(ncid,fz6id,'name',25,    &
                'Krypton impurity fraction')
       status = nf_put_att_text(ncid,fz7id,'name',23,    &
                'Argon impurity fraction')
       status = nf_put_att_text(ncid,fz8id,'name',26,    &
                'Tungsten impurity fraction')
       status = nf_put_att_text(ncid,injcurrid,'name',26,    &
                'Injector Current..........')
       status = nf_put_att_text(ncid,injcur2id,'name',26,    &
                'Injector Current (2)......')
       status = nf_put_att_text(ncid,helicityid,'name',26,    &
                'Magnetic Helicity.........')
       status = nf_put_att_text(ncid,gcur01id,'name',26,    &
                'Group 01 Current (kA).....')
       status = nf_put_att_text(ncid,gcur02id,'name',26,    &
                'Group 02 Current (kA).....')
       status = nf_put_att_text(ncid,gcur03id,'name',26,    &
                'Group 03 Current (kA).....')
       status = nf_put_att_text(ncid,gcur04id,'name',26,    &
                'Group 04 Current (kA).....')
       status = nf_put_att_text(ncid,gcur05id,'name',26,    &
                'Group 05 Current (kA).....')
       status = nf_put_att_text(ncid,gcur06id,'name',26,    &
                'Group 06 Current (kA).....')
       status = nf_put_att_text(ncid,gcur07id,'name',26,    &
                'Group 07 Current (kA).....')
       status = nf_put_att_text(ncid,gcur08id,'name',26,    &
                'Group 08 Current (kA).....')
       status = nf_put_att_text(ncid,gcur09id,'name',26,    &
                'Group 09 Current (kA).....')
       status = nf_put_att_text(ncid,gcur10id,'name',26,    &
                'Group 10 Current (kA).....')
       status = nf_put_att_text(ncid,gcur11id,'name',26,    &
                'Group 11 Current (kA).....')
       status = nf_put_att_text(ncid,gcur12id,'name',26,    &
                'Group 12 Current (kA).....')
       status = nf_put_att_text(ncid,gcur13id,'name',26,    &
                'Group 13 Current (kA).....')
       status = nf_put_att_text(ncid,gcur14id,'name',26,    &
                'Group 14 Current (kA).....')
       status = nf_put_att_text(ncid,gcur15id,'name',26,    &
                'Group 15 Current (kA).....')
       status = nf_put_att_text(ncid,gcur16id,'name',26,    &
                'Group 16 Current (kA).....')
       status = nf_put_att_text(ncid,gcur17id,'name',26,    &
                'Group 17 Current (kA).....')
       status = nf_put_att_text(ncid,gcur18id,'name',26,    &
                'Group 18 Current (kA).....')
       status = nf_put_att_text(ncid,gcur19id,'name',26,    &
                'Group 19 Current (kA).....')
       status = nf_put_att_text(ncid,gcur20id,'name',26,    &
                'Group 20 Current (kA).....')
       status = nf_put_att_text(ncid,gcur21id,'name',26,    &
                'Group 21 Current (kA).....')
       status = nf_put_att_text(ncid,gcur22id,'name',26,    &
                'Group 22 Current (kA).....')
           if(status.ne.0) call movieerr(12,status)
!
!. 2.1 1D variables (poloidal flux)
!
       oneddims(1) = netcdfp
       oneddims(2) = netcdft
       status = nf_def_var(ncid,'xsv2',NF_DOUBLE,2,oneddims,xsv2id)
           if(status.ne.0) call movieerr(35,status)
       status = nf_def_var(ncid,'te2',NF_DOUBLE,2,oneddims,te2id)
           if(status.ne.0) call movieerr(36,status)
       status = nf_def_var(ncid,'ti2',NF_DOUBLE,2,oneddims,ti2id)
           if(status.ne.0) call movieerr(37,status)
       status = nf_def_var(ncid,'ane2',NF_DOUBLE,2,oneddims,ane2id)
           if(status.ne.0) call movieerr(38,status)
       status=nf_def_var(ncid,'ajpest2',NF_DOUBLE,2,oneddims,ajpest2id)
           if(status.ne.0) call movieerr(39,status)
       status=nf_def_var(ncid,'prpest2',NF_DOUBLE,2,oneddims,prpest2id)
           if(status.ne.0) call movieerr(40,status)
       status=nf_def_var(ncid,'pppest2',NF_DOUBLE,2,oneddims,pppest2id)
           if(status.ne.0) call movieerr(41,status)
       status = nf_def_var(ncid,'avez2',NF_DOUBLE,2,oneddims,avez2id)
           if(status.ne.0) call movieerr(42,status)
       status = nf_def_var(ncid,'zeff2',NF_DOUBLE,2,oneddims,zeff2id)
       status = nf_def_var(ncid,'etpara',NF_DOUBLE,2,oneddims,etparaid)
           if(status.ne.0) call movieerr(43,status)
       status = nf_def_var(ncid,'qprof2',NF_DOUBLE,2,oneddims,qprof2id)
           if(status.ne.0) call movieerr(44,status)
!
       status = nf_def_var(ncid,'gary',NF_DOUBLE,2,oneddims,garyid)
           if(status.ne.0) call movieerr(44,status)
       status = nf_def_var(ncid,'ggprime',NF_DOUBLE,2,oneddims,ggprimeid)
           if(status.ne.0) call movieerr(44,status)
       status = nf_def_var(ncid,'pres',NF_DOUBLE,2,oneddims,presid)
           if(status.ne.0) call movieerr(44,status)
       status = nf_def_var(ncid,'presp',NF_DOUBLE,2,oneddims,prespid)
       status = nf_def_var(ncid,'chie',NF_DOUBLE,2,oneddims,chieid)
       status = nf_def_var(ncid,'chii',NF_DOUBLE,2,oneddims,chiiid)
       status = nf_def_var(ncid,'anhy',NF_DOUBLE,2,oneddims,anhyid)
       status = nf_def_var(ncid,'anhe',NF_DOUBLE,2,oneddims,anheid)
       status = nf_def_var(ncid,'anz1',NF_DOUBLE,2,oneddims,anz1id)
       status = nf_def_var(ncid,'anz2',NF_DOUBLE,2,oneddims,anz2id)
       status = nf_def_var(ncid,'anz3',NF_DOUBLE,2,oneddims,anz3id)
       status = nf_def_var(ncid,'anz4',NF_DOUBLE,2,oneddims,anz4id)
       status = nf_def_var(ncid,'anz5',NF_DOUBLE,2,oneddims,anz5id)
       status = nf_def_var(ncid,'anz6',NF_DOUBLE,2,oneddims,anz6id)
       status = nf_def_var(ncid,'anz7',NF_DOUBLE,2,oneddims,anz7id)
       status = nf_def_var(ncid,'anz8',NF_DOUBLE,2,oneddims,anz8id)
       status = nf_def_var(ncid,'ajto',NF_DOUBLE,2,oneddims,ajtoid)
       status = nf_def_var(ncid,'ajbs',NF_DOUBLE,2,oneddims,ajbsid)
       status = nf_def_var(ncid,'ajnb',NF_DOUBLE,2,oneddims,ajnbid)
       status = nf_def_var(ncid,'ajlh',NF_DOUBLE,2,oneddims,ajlhid)
       status = nf_def_var(ncid,'ajic',NF_DOUBLE,2,oneddims,ajicid)
       status = nf_def_var(ncid,'ajec',NF_DOUBLE,2,oneddims,ajecid)
       status = nf_def_var(ncid,'sab',NF_DOUBLE,2,oneddims,sabid)
       status = nf_def_var(ncid,'henb',NF_DOUBLE,2,oneddims,henbid)
       status = nf_def_var(ncid,'hinb',NF_DOUBLE,2,oneddims,hinbid)
       status = nf_def_var(ncid,'heic',NF_DOUBLE,2,oneddims,heicid)
       status = nf_def_var(ncid,'hiic',NF_DOUBLE,2,oneddims,hiicid)
       status = nf_def_var(ncid,'helh',NF_DOUBLE,2,oneddims,helhid)
       status = nf_def_var(ncid,'hilh',NF_DOUBLE,2,oneddims,hilhid)
       status = nf_def_var(ncid,'heec',NF_DOUBLE,2,oneddims,heecid)
       status = nf_def_var(ncid,'heal',NF_DOUBLE,2,oneddims,healid)
       status = nf_def_var(ncid,'hial',NF_DOUBLE,2,oneddims,hialid)
       status = nf_def_var(ncid,'hbre',NF_DOUBLE,2,oneddims,hbreid)
       status = nf_def_var(ncid,'hcyc',NF_DOUBLE,2,oneddims,hcycid)
       status = nf_def_var(ncid,'hlin',NF_DOUBLE,2,oneddims,hlinid)
       status = nf_def_var(ncid,'vl',NF_DOUBLE,2,oneddims,vlid)
       status = nf_def_var(ncid,'tflx',NF_DOUBLE,2,oneddims,tflxid)
           if(status.ne.0) call movieerr(44,status)
!
       status = nf_put_att_text(ncid,xsv2id,'name',28,     &
                'poloidal flux per radian (w)')
           if(status.ne.0) call movieerr(35,status)        
       status = nf_put_att_text(ncid,te2id,'name',25,      &
                'electron temperature (ev)')
           if(status.ne.0) call movieerr(36,status)
       status = nf_put_att_text(ncid,ti2id,'name',20,      &
                'ion temperature (ev)')
           if(status.ne.0) call movieerr(37,status)
       status = nf_put_att_text(ncid,ane2id,'name',22,     &
                'electron density (m-3)')
           if(status.ne.0) call movieerr(38,status)
       status = nf_put_att_text(ncid,ajpest2id,'name',32,  &
                'surface averaged current density')
           if(status.ne.0) call movieerr(39,status)
       status = nf_put_att_text(ncid,prpest2id,'name',8,   &
                'pressure')
           if(status.ne.0) call movieerr(40,status)
       status = nf_put_att_text(ncid,pppest2id,'name',40,  &
                'derivative of pressure wrt poloidal flux')
           if(status.ne.0) call movieerr(41,status)
       status = nf_put_att_text(ncid,avez2id,'name',9,     &
                'average z')
           if(status.ne.0) call movieerr(42,status)
       status = nf_put_att_text(ncid,zeff2id,'name',11,    &
                'z-effective')
           if(status.ne.0) call movieerr(43,status) 
       status = nf_put_att_text(ncid,etparaid,'name',20,    &
                'parallel resistivity')
           if(status.ne.0) call movieerr(43,status) 
       status = nf_put_att_text(ncid,qprof2id,'name',9,    &
                'q-profile') 
           if(status.ne.0) call movieerr(44,status)
!
       status = nf_put_att_text(ncid,garyid,'name',18,     &
                'RB_T-profile (T*m)')
           if(status.ne.0) call movieerr(44,status)
       status = nf_put_att_text(ncid,ggprimeid,'name',35,  &
                'RB_T*d(RB_T)/dPsi (T*m)**2/(Wb/Rad)')
           if(status.ne.0) call movieerr(44,status)
       status = nf_put_att_text(ncid,presid,'name',18,     &
                'EFIT Pressure (Pa)') 
           if(status.ne.0) call movieerr(44,status) 
       status = nf_put_att_text(ncid,prespid,'name',26,    &
                'EFIT dP/dPsi (Pa)/(Wb/Rad)') 
       status = nf_put_att_text(ncid,chieid,'name',32,     &
                'Elect thermal diffusivity (m2/s)')
       status = nf_put_att_text(ncid,chiiid,'name',30,     &
                'Ion thermal diffusivity (m2/s)') 
       status = nf_put_att_text(ncid,anhyid,'name',24,     &
                'Hydrogenic density (/m3)') 
       status = nf_put_att_text(ncid,anheid,'name',20,     &
                'Helium density (/m3)')
       status = nf_put_att_text(ncid,anz1id,'name',20,     &
                'Oxygen density (/m3)')
       status = nf_put_att_text(ncid,anz2id,'name',20,     &
                'Carbon density (/m3)')
       status = nf_put_att_text(ncid,anz3id,'name',18,     &
                'Iron density (/m3)')
       status = nf_put_att_text(ncid,anz4id,'name',23,     &
                'Beryllium density (/m3)')
       status = nf_put_att_text(ncid,anz5id,'name',21,     &
                'Neon density (/m3)')
       status = nf_put_att_text(ncid,anz6id,'name',19,     &
                'Krypton density (/m3)') 
       status = nf_put_att_text(ncid,anz7id,'name',24,     &
                'Argon density (/m3)')   
       status = nf_put_att_text(ncid,anz8id,'name',24,     &
                'Tungsten density (/m3)')
       status = nf_put_att_text(ncid,ajtoid,'name',34,     &
                'Total par current density (A/m2-T)') 
       status = nf_put_att_text(ncid,ajbsid,'name',38,     &
                'Bootstrap par current density (A/m2-T)')
       status = nf_put_att_text(ncid,ajnbid,'name',41,     &
                'Neutral beam par current density (A/m2-T)')
       status = nf_put_att_text(ncid,ajlhid,'name',41,     &
                'Lower hybrid par current density (A/m2-T)')
       status = nf_put_att_text(ncid,ajicid,'name',36,     &
                'Ion cyc par current density (A/m2-T)')
       status = nf_put_att_text(ncid,ajecid,'name',37,     &
                'Elec cyc par current density (A/m2-T)') 
       status = nf_put_att_text(ncid,sabid,'name',20,      &
                'Surface ave Btor (T)') 
       status = nf_put_att_text(ncid,henbid,'name',31,     &
                'Elec heating from neut beam (W)')
       status = nf_put_att_text(ncid,hinbid,'name',30,     &
                'Ion heating from neut beam (W)')
       status = nf_put_att_text(ncid,heicid,'name',29,     &
                'Elec heating from ion cyc (W)') 
       status = nf_put_att_text(ncid,hiicid,'name',28,     &
                'Ion heating from ion cyc (W)') 
       status = nf_put_att_text(ncid,helhid,'name',27,     &
                'Elec heating from lower hyb (W)')
       status = nf_put_att_text(ncid,hilhid,'name',30,     &
                'Ion heating from lower hyb (W)') 
       status = nf_put_att_text(ncid,heecid,'name',30,     &
                'Elec heating from elec cyc (W)')
       status = nf_put_att_text(ncid,healid,'name',27,     &
                'Elec heating from alphas (W)') 
       status = nf_put_att_text(ncid,hialid,'name',27,     &
                'Ion heating from alphas (W)') 
       status = nf_put_att_text(ncid,hbreid,'name',28,     &
                'Bremsstrahlung radiation (W)') 
       status = nf_put_att_text(ncid,hcycid,'name',23,     &
                'Cyclotron radiation (W)') 
       status = nf_put_att_text(ncid,hlinid,'name',18,     &
                'Line radiation (W)') 
       status = nf_put_att_text(ncid,vlid,'name',16,       &
                'Loop voltage (V)') 
       status = nf_put_att_text(ncid,tflxid,'name',19,     &
                'Normalized tor flux')
           if(status.ne.0) call movieerr(44,status)
!
!. 2.2 1D variables (major radius)
!
       oneddims(1) = netcdfr
       oneddims(2) = netcdft
       status = nf_def_var(ncid,'rcoord',NF_DOUBLE,2,oneddims,rcoordid)
           if(status.ne.0) call movieerr(35,status)
       status = nf_def_var(ncid,'ter',NF_DOUBLE,2,oneddims,terid)
           if(status.ne.0) call movieerr(36,status)
       status = nf_def_var(ncid,'tir',NF_DOUBLE,2,oneddims,tirid)
           if(status.ne.0) call movieerr(37,status)
       status = nf_def_var(ncid,'aner',NF_DOUBLE,2,oneddims,anerid)
           if(status.ne.0) call movieerr(38,status)
       status=nf_def_var(ncid,'ajphir',NF_DOUBLE,2,oneddims,ajphirid)
           if(status.ne.0) call movieerr(39,status)
       status=nf_def_var(ncid,'qr',NF_DOUBLE,2,oneddims,qrid)
           if(status.ne.0) call movieerr(40,status)
!
       status = nf_put_att_text(ncid,rcoordid,'name',16,   &
                'major radius (m)')
           if(status.ne.0) call movieerr(35,status)       
       status = nf_put_att_text(ncid,terid,'name',25,      &
                'electron temperature (ev)')
           if(status.ne.0) call movieerr(36,status)
       status = nf_put_att_text(ncid,tirid,'name',20,      &
                'ion temperature (ev)') 
           if(status.ne.0) call movieerr(37,status)
       status = nf_put_att_text(ncid,anerid,'name',22,     &
                'electron density (m-3)')
           if(status.ne.0) call movieerr(38,status)
       status = nf_put_att_text(ncid,ajphirid,'name',24,   &
                'toroidal current density')
           if(status.ne.0) call movieerr(39,status)
       status = nf_put_att_text(ncid,qrid,'name',9,        &
                'q-profile')
           if(status.ne.0) call movieerr(40,status)

!
!. 2.3 1D variables (theta index)
!.      (added by Kumar)
!
         oneddims(1) = netcdfth
         oneddims(2) = netcdft
         status=nf_def_var(ncid,'rbnd',NF_DOUBLE,2,oneddims,rbndid)
            if(status.ne.0) call movieerr(35,status)
         status=nf_def_var(ncid,'zbnd',NF_DOUBLE,2,oneddims,zbndid)
            if(status.ne.0) call movieerr(36,status)

         status = nf_put_att_text(ncid,rbndid,'name',16,     &
                'Boundary coord R')
         if(status.ne.0) call movieerr(35,status)
            status = nf_put_att_text(ncid,zbndid,'name',16,'Boundary coord Z')
         if(status.ne.0) call movieerr(36,status)

!
!   3.   2D variables
         psidims(1) = netcdfx
         psidims(2) = netcdfz
         psidims(3) = netcdft
         status = nf_def_var(ncid,'psi',NF_DOUBLE,3,psidims,psiid)
           if(status.ne.0) call movieerr(13,status)
         status = nf_put_att_text(ncid,psiid,'name',28,           &
                'poloidal flux per radian (w)') 
           if(status.ne.0) call movieerr(13,status) 
!
         status = nf_def_var(ncid,'imask',NF_INT,3,psidims,imaskid) 
           if(status.ne.0) call movieerr(15,status) 
         status = nf_put_att_text(ncid,imaskid,'name',36,         &
                'masking array: 0 plot, 1 do not plot')
           if(status.ne.0) call movieerr(15,status)
! 
         status = nf_def_var(ncid,'divflx',NF_DOUBLE,3,psidims,divflxid)
           if(status.ne.0) call movieerr(14,status)
         status = nf_put_att_text(ncid,divflxid,'name',28,        &
                'poloidal flux in divertor(w)')
           if(status.ne.0) call movieerr(14,status)
!
         status = nf_def_var(ncid,'ajphi',NF_DOUBLE,3,psidims,ajphiid)
           if(status.ne.0) call movieerr(70,status)
         status = nf_put_att_text(ncid,ajphiid,'name',24,         &
               'toroidal current density')
           if(status.ne.0) call movieerr(70,status)
!
         status = nf_def_var(ncid,'pr',NF_DOUBLE,3,psidims,prid)
           if(status.ne.0) call movieerr(71,status)
         status = nf_put_att_text(ncid,prid,'name',15,            &
                'plasma pressure')
           if(status.ne.0) call movieerr(71,status)
!
         status = nf_def_var(ncid,'gs',NF_DOUBLE,3,psidims,gsid)
           if(status.ne.0) call movieerr(72,status)
         status = nf_put_att_text(ncid,gsid,'name',29,            &
                'toroidal field function (T-m)') 
           if(status.ne.0) call movieerr(72,status)
!
!
!........added by SCJ for idata=9 or idata=11 on 2/10
         if(idata.eq.9 .or. idata.eq.11) then
           status = nf_def_var(ncid,'atot' ,NF_DOUBLE,1,netcdft,atotid)
           status = nf_put_att_text(ncid,atotid,'name',27,'Vessel + Plasma current (S)')
!
           status = nf_def_var(ncid,'csum' ,NF_DOUBLE,1,netcdft,csumid)
           status = nf_put_att_text(ncid,csumid,'name',18,   'Vessel current (S)')
!
           status = nf_def_var(ncid,'vesppe' ,NF_DOUBLE,1,netcdft,vesppeid)
           status = nf_put_att_text(ncid,vesppeid,'name',18, 'Plasma current (e)')
!
           status = nf_def_var(ncid,'vescure' ,NF_DOUBLE,1,netcdft,vescureid)
           status = nf_put_att_text(ncid,vescureid,'name',27,'Vessel current (e)')
!
           status = nf_def_var(ncid,'totcure' ,NF_DOUBLE,1,netcdft,totcureid)
           status = nf_put_att_text(ncid,totcureid,'name',27,'Vessel + Plasma current (e)')
!
           status = nf_def_dim(ncid,'current groups',ncurrentg, netcdfc)
           status = nf_def_dim(ncid,'flux loops',numfluxloopsg, netcdffl)
           oneddims(1) = netcdffl
           oneddims(2) = netcdft
!
           status = nf_def_var(ncid,'pflux' ,NF_DOUBLE,2,oneddims,pfluxid)
           status = nf_put_att_text(ncid,pfluxid,'name',17,'poloidal flux (S)')
!
           status = nf_def_var(ncid,'eflux' ,NF_DOUBLE,2,oneddims,efluxid)
           status = nf_put_att_text(ncid,efluxid,'name',17,'poloidal flux (e)')
           oneddims(1) = netcdfc
!
           status = nf_def_var(ncid,'ppcur' ,NF_DOUBLE,2,oneddims,ppcurid)
           status = nf_put_att_text(ncid,ppcurid,'name',27, 'experimental group currents')
!
           status = nf_def_var(ncid,'fbcur' ,NF_DOUBLE,2,oneddims,fbcurid)
           status = nf_put_att_text(ncid,fbcurid,'name',23, 'feedback group currents')
!
           status = nf_def_var(ncid,'csuma' ,NF_DOUBLE,2,oneddims,csumaid)
           status = nf_put_att_text(ncid,csumaid,'name',21,'actual group currents')
         endif
!
!
         status = nf_enddef(ncid)
          if(status.ne.0) call movieerr(21,status)   
!
!
!........write coordinate variables using ezcdf..[note:  2 to nx+1, 
!                                                      1 to nzdimcdf]
         call cdf_write (ncid,'xary', xgridcdf,ier)
          if(ier.ne.0) call movieerr(4,ier)
         call cdf_write (ncid,'zary', zgridcdf,ier)
          if(ier.ne.0) call movieerr(7,ier)
         nreccdf = 0
!
      else     !.....on restart
!
!.....initialize restart run
      call cdf_open( ncid, 'movie.cdf', 'a' , ier)
          if(ier.ne.0) call movieerr(15,ier)
!
      status = nf_inq_varid(ncid,'times',timeid)
          if(status.ne.0) call movieerr(17,status)
      status = nf_inq_varid(ncid,'psimin',psiminid)
          if(status.ne.0) call movieerr(18,status)
      status = nf_inq_varid(ncid,'psilim',psilimid)
          if(status.ne.0) call movieerr(15,status)
      status = nf_inq_varid(ncid,'psisol',psisolid)
          if(status.ne.0) call movieerr(16,status)
      status = nf_inq_varid(ncid,'phalo',phaloid)
          if(status.ne.0) call movieerr(15,status)
      status = nf_inq_varid(ncid,'npsit',npsitid)
          if(status.ne.0) call movieerr(60,status)
!
      status = nf_inq_varid(ncid,'pohmic',pohmicid)
          if(status.ne.0) call movieerr(17,status)
      status = nf_inq_varid(ncid,'paux',pauxid)
          if(status.ne.0) call movieerr(17,status)
      status = nf_inq_varid(ncid,'palpha',palphaid)
          if(status.ne.0) call movieerr(17,status)
      status = nf_inq_varid(ncid,'bre',breid)
          if(status.ne.0) call movieerr(17,status)
      status = nf_inq_varid(ncid,'cyc',cycid)
          if(status.ne.0) call movieerr(17,status)
      status = nf_inq_varid(ncid,'imp',impid)
          if(status.ne.0) call movieerr(17,status)
      status = nf_inq_varid(ncid,'pnb',pnbid)
          if(status.ne.0) call movieerr(17,status)
      status = nf_inq_varid(ncid,'plh',plhid)
          if(status.ne.0) call movieerr(17,status)
      status = nf_inq_varid(ncid,'pic',picid)
          if(status.ne.0) call movieerr(17,status)
      status = nf_inq_varid(ncid,'pec',pecid)
          if(status.ne.0) call movieerr(17,status)
!
      status = nf_inq_varid(ncid,'apl',aplid)
          if(status.ne.0) call movieerr(17,status)
      status = nf_inq_varid(ncid,'bsi',bsiid)
          if(status.ne.0) call movieerr(17,status)
      status = nf_inq_varid(ncid,'cdi',cdiid)
          if(status.ne.0) call movieerr(17,status)
      status = nf_inq_varid(ncid,'alhi',alhiid)
          if(status.ne.0) call movieerr(17,status)
      status = nf_inq_varid(ncid,'fwi',fwiid)
          if(status.ne.0) call movieerr(17,status)
      status = nf_inq_varid(ncid,'eci',eciid)
          if(status.ne.0) call movieerr(17,status)
      status = nf_inq_varid(ncid,'tni',tniid)
      status = nf_inq_varid(ncid,'vloop',vloopid)
      status = nf_inq_varid(ncid,'xplas',xplasid)
      status = nf_inq_varid(ncid,'gzero',gzeroid)
      status = nf_inq_varid(ncid,'tflux',tfluxid)
      status = nf_inq_varid(ncid,'rmajor',rmajorid)
      status = nf_inq_varid(ncid,'rminor',rminorid)
      status = nf_inq_varid(ncid,'delta',deltaid)
      status = nf_inq_varid(ncid,'kappa',kappaid)
      status = nf_inq_varid(ncid,'xmag',xmagid)
      status = nf_inq_varid(ncid,'zmag',zmagid)
      status = nf_inq_varid(ncid,'rq1',rq1id)
      status = nf_inq_varid(ncid,'zeff',zeffid)
      status = nf_inq_varid(ncid,'qzero',qzeroid)
      status = nf_inq_varid(ncid,'qmin',qminid)
      status = nf_inq_varid(ncid,'q95',q95id)
      status = nf_inq_varid(ncid,'qcyl',qcylid)
      status = nf_inq_varid(ncid,'ali1',ali1id)
      status = nf_inq_varid(ncid,'ali3',ali3id)
      status = nf_inq_varid(ncid,'betat',betatid)
      status = nf_inq_varid(ncid,'betap',betapid)
      status = nf_inq_varid(ncid,'wth',wthid)
      status = nf_inq_varid(ncid,'taue1',taue1id)
      status = nf_inq_varid(ncid,'taue2',taue2id)
      status = nf_inq_varid(ncid,'ah98',ah98id)
      status = nf_inq_varid(ncid,'ah89',ah89id)
      status = nf_inq_varid(ncid,'rzero',rzeroid)
      status = nf_inq_varid(ncid,'rline',rlineid)
      status = nf_inq_varid(ncid,'rvol',rvolid)
      status = nf_inq_varid(ncid,'rorgr',rorgrid)
      status = nf_inq_varid(ncid,'te0',te0id)
      status = nf_inq_varid(ncid,'ti0',ti0id)
      status = nf_inq_varid(ncid,'flxst',flxstid)
      status = nf_inq_varid(ncid,'vstax',vstaxid)
      status = nf_inq_varid(ncid,'vsiax',vsiaxid)
      status = nf_inq_varid(ncid,'vsrax',vsraxid)
      status = nf_inq_varid(ncid,'vsipoy',vsipoyid)
      status = nf_inq_varid(ncid,'vsrpoy',vsrpoyid)
      status = nf_inq_varid(ncid,'thalo',thaloid)
      status = nf_inq_varid(ncid,'whalo',whaloid)
      status = nf_inq_varid(ncid,'tped',tpedid)
      status = nf_inq_varid(ncid,'psiped',psipedid)
      status = nf_inq_varid(ncid,'fz1',fz1id)
      status = nf_inq_varid(ncid,'fz2',fz2id)
      status = nf_inq_varid(ncid,'fz3',fz3id)
      status = nf_inq_varid(ncid,'fz4',fz4id)
      status = nf_inq_varid(ncid,'fz5',fz5id)
      status = nf_inq_varid(ncid,'fz6',fz6id)
      status = nf_inq_varid(ncid,'fz7',fz7id)
      status = nf_inq_varid(ncid,'fz8',fz8id)
      status = nf_inq_varid(ncid,'injcurr',injcurrid)
      status = nf_inq_varid(ncid,'injcur2',injcur2id)
      status = nf_inq_varid(ncid,'helicity',helicityid)
      status = nf_inq_varid(ncid,'gcur01',gcur01id)
      status = nf_inq_varid(ncid,'gcur02',gcur02id)
      status = nf_inq_varid(ncid,'gcur03',gcur03id)
      status = nf_inq_varid(ncid,'gcur04',gcur04id)
      status = nf_inq_varid(ncid,'gcur05',gcur05id)
      status = nf_inq_varid(ncid,'gcur06',gcur06id)
      status = nf_inq_varid(ncid,'gcur07',gcur07id)
      status = nf_inq_varid(ncid,'gcur08',gcur08id)
      status = nf_inq_varid(ncid,'gcur09',gcur09id)
      status = nf_inq_varid(ncid,'gcur10',gcur10id)
      status = nf_inq_varid(ncid,'gcur11',gcur11id)
      status = nf_inq_varid(ncid,'gcur12',gcur12id)
      status = nf_inq_varid(ncid,'gcur13',gcur13id)
      status = nf_inq_varid(ncid,'gcur14',gcur14id)
      status = nf_inq_varid(ncid,'gcur15',gcur15id)
      status = nf_inq_varid(ncid,'gcur16',gcur16id)
      status = nf_inq_varid(ncid,'gcur17',gcur17id)
      status = nf_inq_varid(ncid,'gcur18',gcur18id)
      status = nf_inq_varid(ncid,'gcur19',gcur19id)
      status = nf_inq_varid(ncid,'gcur20',gcur20id)
      status = nf_inq_varid(ncid,'gcur21',gcur21id)
      status = nf_inq_varid(ncid,'gcur22',gcur22id)
          if(status.ne.0) call movieerr(17,status)
!    
      status = nf_inq_varid(ncid,'xsv2',xsv2id)
          if(status.ne.0) call movieerr(61,status)
      status = nf_inq_varid(ncid,'te2',te2id)
          if(status.ne.0) call movieerr(62,status)
      status = nf_inq_varid(ncid,'ti2',ti2id)
          if(status.ne.0) call movieerr(63,status)
      status = nf_inq_varid(ncid,'ane2',ane2id)
          if(status.ne.0) call movieerr(64,status)
      status = nf_inq_varid(ncid,'ajpest2',ajpest2id)
          if(status.ne.0) call movieerr(65,status)
      status = nf_inq_varid(ncid,'prpest2',prpest2id)
          if(status.ne.0) call movieerr(66,status)
      status = nf_inq_varid(ncid,'pppest2',pppest2id)
          if(status.ne.0) call movieerr(67,status)
      status = nf_inq_varid(ncid,'avez2',avez2id)
          if(status.ne.0) call movieerr(68,status)
      status = nf_inq_varid(ncid,'zeff2',zeff2id)
      status = nf_inq_varid(ncid,'etpara',etparaid)
          if(status.ne.0) call movieerr(69,status)
      status = nf_inq_varid(ncid,'qprof2',qprof2id)
          if(status.ne.0) call movieerr(70,status)
!
      status = nf_inq_varid(ncid,'gary',garyid)
          if(status.ne.0) call movieerr(70,status)
      status = nf_inq_varid(ncid,'ggprime',ggprimeid)
          if(status.ne.0) call movieerr(70,status)
      status = nf_inq_varid(ncid,'pres',presid)
          if(status.ne.0) call movieerr(70,status)
      status = nf_inq_varid(ncid,'presp',prespid)
      status = nf_inq_varid(ncid,'chie',chieid)
      status = nf_inq_varid(ncid,'chii',chiiid)
      status = nf_inq_varid(ncid,'anhy',anhyid)
      status = nf_inq_varid(ncid,'anhe',anheid)
      status = nf_inq_varid(ncid,'anz1',anz1id)
      status = nf_inq_varid(ncid,'anz2',anz2id)
      status = nf_inq_varid(ncid,'anz3',anz3id)
      status = nf_inq_varid(ncid,'anz4',anz4id)
      status = nf_inq_varid(ncid,'anz5',anz5id)
      status = nf_inq_varid(ncid,'anz6',anz6id)
      status = nf_inq_varid(ncid,'anz7',anz7id)
      status = nf_inq_varid(ncid,'anz8',anz8id)
      status = nf_inq_varid(ncid,'ajto',ajtoid)
      status = nf_inq_varid(ncid,'ajbs',ajbsid)
      status = nf_inq_varid(ncid,'ajnb',ajnbid)
      status = nf_inq_varid(ncid,'ajlh',ajlhid)
      status = nf_inq_varid(ncid,'ajic',ajicid)
      status = nf_inq_varid(ncid,'ajec',ajecid)
      status = nf_inq_varid(ncid,'sab',sabid)
      status = nf_inq_varid(ncid,'henb',henbid)
      status = nf_inq_varid(ncid,'hinb',hinbid)
      status = nf_inq_varid(ncid,'heic',heicid)
      status = nf_inq_varid(ncid,'hiic',hiicid)
      status = nf_inq_varid(ncid,'helh',helhid)
      status = nf_inq_varid(ncid,'hilh',hilhid)
      status = nf_inq_varid(ncid,'heec',heecid)
      status = nf_inq_varid(ncid,'heal',healid)
      status = nf_inq_varid(ncid,'hial',hialid)
      status = nf_inq_varid(ncid,'hbre',hbreid)
      status = nf_inq_varid(ncid,'hcyc',hcycid)
      status = nf_inq_varid(ncid,'hlin',hlinid)
      status = nf_inq_varid(ncid,'vl',vlid)
      status = nf_inq_varid(ncid,'tflx',tflxid)
          if(status.ne.0) call movieerr(70,status)
!    
      status = nf_inq_varid(ncid,'rcoord',rcoordid)
          if(status.ne.0) call movieerr(61,status)
      status = nf_inq_varid(ncid,'ter',terid)
          if(status.ne.0) call movieerr(62,status)
      status = nf_inq_varid(ncid,'tir',tirid)
          if(status.ne.0) call movieerr(63,status)
      status = nf_inq_varid(ncid,'aner',anerid)
          if(status.ne.0) call movieerr(64,status)
      status = nf_inq_varid(ncid,'ajphir',ajphirid)
          if(status.ne.0) call movieerr(65,status)
      status = nf_inq_varid(ncid,'qr',qrid)
          if(status.ne.0) call movieerr(66,status)
!
!      Boundary pts (Kumar)
!
      status = nf_inq_varid(ncid,'rbnd',rbndid)
      if(status.ne.0) call movieerr(61,status)
      status = nf_inq_varid(ncid,'zbnd',zbndid)
      if(status.ne.0) call movieerr(62,status)
!
      status = nf_inq_varid(ncid,'psi',psiid)
          if(status.ne.0) call movieerr(19,status)
      status = nf_inq_varid(ncid,'imask',imaskid)
          if(status.ne.0) call movieerr(21,status)
      status = nf_inq_varid(ncid,'divflx',divflxid)
          if(status.ne.0) call movieerr(20,status)
      status = nf_inq_varid(ncid,'ajphi',ajphiid)
          if(status.ne.0) call movieerr(19,status)
      status = nf_inq_varid(ncid,'pr',prid)
          if(status.ne.0) call movieerr(19,status)
      status = nf_inq_varid(ncid,'gs',gsid)
          if(status.ne.0) call movieerr(19,status)
!
!
!.....idata=9 or idata=11 data SCJ (02/10)
      if(idata.eq.9 .or. idata.eq.11) then
        status = nf_inq_varid(ncid,'atot',atotid)
        status = nf_inq_varid(ncid,'csum',csumid)
        status = nf_inq_varid(ncid,'vesppe',vesppe)
        status = nf_inq_varid(ncid,'vescure',vescureid)
        status = nf_inq_varid(ncid,'totcure',totcureid)
        status = nf_inq_varid(ncid,'pflux',pfluxid)
        status = nf_inq_varid(ncid,'eflux',efluxid)
        status = nf_inq_varid(ncid,'ppcur',ppcurid)
        status = nf_inq_varid(ncid,'fbcur',fbcurid)
        status = nf_inq_varid(ncid,'csuma',csumaid)
      endif

!
      status = nf_enddef(ncid)
          if(status.ne.0) call movieerr(21,status)   
      endif
!
!   
      ifirstmd = 1
      endif
!
!-----------------------------------------------------------------------
!  START OF SECTION THAT WRITES DATA FOR ONE TIME POINT
!
!.....first define some auxialliary arrays
!
      do 100 j=1,npsit
      ps = xsv2(j)
      call peval(ps,1+isurf,pval,ppval,imag,jmag)
      prpest2(j) = pval
      pppest2(j) = ppval
      te2(j)    = .5*(te(j)   +te(j+1))
      ane2(j)   = .5*(ane(j)  +ane(j+1))
      ti2(j)    = .5*(ti(j)   +ti(j+1))
      avez2(j)  = .5*(avez(j) +avez(j+1))
      zeff2(j) = .5*(zeffa(j)+zeffa(j+1))
      etparas(j) = etpara(j)*udsr
      if(j.eq.1 .or. j.eq.npsit) go to 100
      t1 = rdpsi*tpi*(qprof2(j))**2/(xmja2(j))**2
      t3 = gxmja(j+1)*xmja(j+1)/(.5*(qprof2(j)+qprof2(j+1)))
      t5 = gxmja(j)*xmja(j)/(.5*(qprof2(j)+qprof2(j-1)))
      ajpest2(j) = t1*(t3-t5)
  100 continue
      write(*,9100) kcycle,times,etparas(1),etparas(2),etparas(3)
 9100 format("debug:  kcycle,times,etparas(1),etparas(2),etparas(3)",i6,1p4e12.4)
      ajpest2(1) = 2.*ajpest2(2) - ajpest2(3)
      ajpest2(npsit) = 0.
      te2(1)    = 1.5*te(2)    - 0.5*te(3)
      ane2(1)   = 1.5*ane(2)   - 0.5*ane(3)
      ti2(1)    = 1.5*ti(2)    - 0.5*ti(3)
      avez2(1)  = 1.5*avez(2)  - 0.5*avez(3)
      zeff2(1) = 1.5*zeffa(2) - 0.5*zeffa(3)
!
!......Boundary data
!    (Following by Kumar)
!

!     Initialize rbnd and zbnd arrays

       do j=1,2*pnparts
         rbnd(j)=1.e+6
         zbnd(j)=1.e+6
       enddo
!
       do j=2,kmax+1
         rbnd(j-1)=xplot(1,j)
         zbnd(j-1)=zplot(1,j)
       enddo
!
!.....warning, this need to be checked out
       if (isym.eq.1) then
        do j=3,kmax
         rbnd(kmax+j-2)=xplot(1,kmax-j+3)
         zbnd(kmax+j-2)=-zplot(1,kmax-j+3)
        enddo
       endif
!
!
!.....increment the size of the time dimension
      nreccdf = nreccdf + 1
!
      start2(1) = 1
      start2(2) = nreccdf    
!
      count2(1) = npsit
      count2(2) = 1

! Boundary point related (Kumar)

      count2n(1) = 2*pnparts
      count2n(2) = 1
!
      start3(1) = 1
      start3(2) = 1
      start3(3) = nreccdf
!
      count3(1) = nx
      count3(2) = nzdimcdf
      count3(3) = 1
!
!.....write variables for this time point:
! 1. Scalar variables
      status = nf_put_var1_double(ncid,timeid,nreccdf,times)
          if(status.ne.0) call movieerr(22,status)
      status = nf_put_var1_double(ncid,psiminid,nreccdf,psimin)
          if(status.ne.0) call movieerr(32,status)
      psilimw = psilim
      if(lrswtch.gt.0) psilimw = psilim - 0.002
      status = nf_put_var1_double(ncid,psilimid,nreccdf,psilimw)
          if(status.ne.0) call movieerr(42,status)
!
!     define psisol here
      psisol = psilim+eps
      if(iplim .lt. 0) then
        do i=imag,nx
          ii = i
          psival = psi(i,nh)
          if(psival .ge. psilim) go to 407
        enddo
        ineg=59
        return
 407    continue
        xedge = xary(ii-1) + deex*(psilim    -psi(ii-1,nh))     &
                         /     (psi(ii,nh)-psi(ii-1,nh))
        call grap3(zmag,xedge+.05,psisol)
      endif
      if(lrswtch.gt.0) psisol = psilim+.002
!
      status = nf_put_var1_double(ncid,psisolid,nreccdf,psisol)
          if(status.ne.0) call movieerr(43,status)
!
!.....changed 10/1/2011   scj
!     phalo should be defined only in subroutine delpdef
!      phalo = psilim + whalos*(psilim - psimin)
      status = nf_put_var1_double(ncid,phaloid,nreccdf,phalo)
          if(status.ne.0) call movieerr(84,status)
      status = nf_put_var1_int(ncid,npsitid,nreccdf,npsit)
          if(status.ne.0) call movieerr(83,status)
!
      status = nf_put_var1_double(ncid,pohmicid,nreccdf,global(47))
          if(status.ne.0) call movieerr(22,status)
      status = nf_put_var1_double(ncid,pauxid,nreccdf,global(48))
          if(status.ne.0) call movieerr(22,status)
      status = nf_put_var1_double(ncid,palphaid,nreccdf,global(46))
          if(status.ne.0) call movieerr(22,status)
      status = nf_put_var1_double(ncid,breid,nreccdf,global(119))
          if(status.ne.0) call movieerr(22,status)
      status = nf_put_var1_double(ncid,cycid,nreccdf,global(120))
          if(status.ne.0) call movieerr(22,status)
      status = nf_put_var1_double(ncid,impid,nreccdf,global(121))
          if(status.ne.0) call movieerr(22,status)
!
      status = nf_put_var1_double(ncid,aplid,nreccdf,global(8))
          if(status.ne.0) call movieerr(22,status)
      status = nf_put_var1_double(ncid,bsiid,nreccdf,global(141))
          if(status.ne.0) call movieerr(22,status)
      status = nf_put_var1_double(ncid,cdiid,nreccdf,global(142))
          if(status.ne.0) call movieerr(22,status)
      status = nf_put_var1_double(ncid,alhiid,nreccdf,global(146))
          if(status.ne.0) call movieerr(22,status)
      status = nf_put_var1_double(ncid,fwiid,nreccdf,global(198))
          if(status.ne.0) call movieerr(22,status)
      status = nf_put_var1_double(ncid,eciid,nreccdf,global(242))
          if(status.ne.0) call movieerr(22,status)
      status = nf_put_var1_double(ncid,tniid,nreccdf,global(213))
      status = nf_put_var1_double(ncid,vloopid,nreccdf,global(13))
      status = nf_put_var1_double(ncid,xplasid,nreccdf,xplas)
      status = nf_put_var1_double(ncid,gzeroid,nreccdf,gzero)
      tflux = dpsi*(npsit-1)
      status = nf_put_var1_double(ncid,tfluxid,nreccdf,tflux)
      status = nf_put_var1_double(ncid,rmajorid,nreccdf,global(34))
      status = nf_put_var1_double(ncid,rminorid,nreccdf,global(35))
      status = nf_put_var1_double(ncid,deltaid,nreccdf,global(36))
      status = nf_put_var1_double(ncid,kappaid,nreccdf,global(37))
      status = nf_put_var1_double(ncid,pnbid,nreccdf,pnb)
      status = nf_put_var1_double(ncid,plhid,nreccdf,plh)
      status = nf_put_var1_double(ncid,picid,nreccdf,pic)
      status = nf_put_var1_double(ncid,pecid,nreccdf,pec)
      status = nf_put_var1_double(ncid,xmagid,nreccdf,global(7))
      status = nf_put_var1_double(ncid,zmagid,nreccdf,global(6))
      status = nf_put_var1_double(ncid,rq1id,nreccdf,global(151))
      status = nf_put_var1_double(ncid,zeffid,nreccdf,global(135))
      status = nf_put_var1_double(ncid,qzeroid,nreccdf,global(14))
      status = nf_put_var1_double(ncid,qminid,nreccdf,qminim)
      status = nf_put_var1_double(ncid,q95id,nreccdf,global(61))
      status = nf_put_var1_double(ncid,qcylid,nreccdf,global(63))
      status = nf_put_var1_double(ncid,ali1id,nreccdf,global(152))
      status = nf_put_var1_double(ncid,ali3id,nreccdf,2.*global(22))
      status = nf_put_var1_double(ncid,betatid,nreccdf,global(75))
      status = nf_put_var1_double(ncid,betapid,nreccdf,global(21))
      status = nf_put_var1_double(ncid,wthid,nreccdf,global(80))
      status = nf_put_var1_double(ncid,taue1id,nreccdf,global(26))
      status = nf_put_var1_double(ncid,taue2id,nreccdf,global(25))
      status = nf_put_var1_double(ncid,ah98id,nreccdf,ah98)
      status = nf_put_var1_double(ncid,ah89id,nreccdf,ah89)
      status = nf_put_var1_double(ncid,rzeroid,nreccdf,global(79))
      status = nf_put_var1_double(ncid,rlineid,nreccdf,global(77))
      status = nf_put_var1_double(ncid,rvolid,nreccdf,global(78))
      status = nf_put_var1_double(ncid,rorgrid,nreccdf,global(149))
      status = nf_put_var1_double(ncid,te0id,nreccdf,global(27))
      status = nf_put_var1_double(ncid,ti0id,nreccdf,global(28))
      status = nf_put_var1_double(ncid,flxstid,nreccdf,global(42))
      status = nf_put_var1_double(ncid,vstaxid,nreccdf,global(42))
      status = nf_put_var1_double(ncid,vsiaxid,nreccdf,global(40))
      status = nf_put_var1_double(ncid,vsraxid,nreccdf,global(41))
      status = nf_put_var1_double(ncid,vsipoyid,nreccdf,global(168))
      status = nf_put_var1_double(ncid,vsrpoyid,nreccdf,global(169))
      status = nf_put_var1_double(ncid,thaloid,nreccdf,global(164))
      status = nf_put_var1_double(ncid,whaloid,nreccdf,global(165))
      status = nf_put_var1_double(ncid,tpedid,nreccdf,global(218))
      status = nf_put_var1_double(ncid,psipedid,nreccdf,pwidthc)
      status = nf_put_var1_double(ncid,fz1id,nreccdf,fraci(1))
      status = nf_put_var1_double(ncid,fz2id,nreccdf,fraci(2))
      status = nf_put_var1_double(ncid,fz3id,nreccdf,fraci(3))
      status = nf_put_var1_double(ncid,fz4id,nreccdf,fraci(4))
      status = nf_put_var1_double(ncid,fz5id,nreccdf,fraci(5))
      status = nf_put_var1_double(ncid,fz6id,nreccdf,fraci(6))
      status = nf_put_var1_double(ncid,fz7id,nreccdf,fraci(7))
      status = nf_put_var1_double(ncid,fz8id,nreccdf,fraci(8))
      status = nf_put_var1_double(ncid,injcurrid,nreccdf,global(220))
      status = nf_put_var1_double(ncid,injcur2id,nreccdf,global(221))
      status = nf_put_var1_double(ncid,helicityid,nreccdf,global(12))
      status = nf_put_var1_double(ncid,gcur01id,nreccdf,gcur01)
      status = nf_put_var1_double(ncid,gcur02id,nreccdf,gcur02)
      status = nf_put_var1_double(ncid,gcur03id,nreccdf,gcur03)
      status = nf_put_var1_double(ncid,gcur04id,nreccdf,gcur04)
      status = nf_put_var1_double(ncid,gcur05id,nreccdf,gcur05)
      status = nf_put_var1_double(ncid,gcur06id,nreccdf,gcur06)
      status = nf_put_var1_double(ncid,gcur07id,nreccdf,gcur07)
      status = nf_put_var1_double(ncid,gcur08id,nreccdf,gcur08)
      status = nf_put_var1_double(ncid,gcur09id,nreccdf,gcur09)
      status = nf_put_var1_double(ncid,gcur10id,nreccdf,gcur10)
      status = nf_put_var1_double(ncid,gcur11id,nreccdf,gcur11)
      status = nf_put_var1_double(ncid,gcur12id,nreccdf,gcur12)
      status = nf_put_var1_double(ncid,gcur13id,nreccdf,gcur13)
      status = nf_put_var1_double(ncid,gcur14id,nreccdf,gcur14)
      status = nf_put_var1_double(ncid,gcur15id,nreccdf,gcur15)
      status = nf_put_var1_double(ncid,gcur16id,nreccdf,gcur16)
      status = nf_put_var1_double(ncid,gcur17id,nreccdf,gcur17)
      status = nf_put_var1_double(ncid,gcur18id,nreccdf,gcur18)
      status = nf_put_var1_double(ncid,gcur19id,nreccdf,gcur19)
      status = nf_put_var1_double(ncid,gcur20id,nreccdf,gcur20)
      status = nf_put_var1_double(ncid,gcur21id,nreccdf,gcur21)
      status = nf_put_var1_double(ncid,gcur22id,nreccdf,gcur22)
          if(status.ne.0) call movieerr(22,status)
!
! 2. 1D Variables
      status = nf_put_vara_double(ncid,xsv2id,start2,count2,xsv2)
          if(status.ne.0) call movieerr(43,status)
      status = nf_put_vara_double(ncid,te2id,start2,count2,te2)
          if(status.ne.0) call movieerr(44,status)
      status = nf_put_vara_double(ncid,ti2id,start2,count2,ti2)
          if(status.ne.0) call movieerr(45,status)
      status = nf_put_vara_double(ncid,ane2id,start2,count2,ane2)
          if(status.ne.0) call movieerr(46,status)
      status = nf_put_vara_double(ncid,ajpest2id,start2,count2,ajpest2)
          if(status.ne.0) call movieerr(47,status)
      status = nf_put_vara_double(ncid,prpest2id,start2,count2,prpest2)
          if(status.ne.0) call movieerr(48,status)
      status = nf_put_vara_double(ncid,pppest2id,start2,count2,pppest2)
          if(status.ne.0) call movieerr(49,status)
      status = nf_put_vara_double(ncid,avez2id,start2,count2,avez2)
          if(status.ne.0) call movieerr(50,status)
      status = nf_put_vara_double(ncid,zeff2id,start2,count2,zeff2)
          if(status.ne.0) call movieerr(51,status)
      status = nf_put_vara_double(ncid,etparaid,start2,count2,etparas)
          if(status.ne.0) call movieerr(51,status)
      status = nf_put_vara_double(ncid,qprof2id,start2,count2,qprof2)
          if(status.ne.0) call movieerr(52,status)
!
      status = nf_put_vara_double(ncid,garyid,start2,count2,gary)
          if(status.ne.0) call movieerr(52,status)
      status = nf_put_vara_double(ncid,ggprimeid,start2,count2,ggprime)
          if(status.ne.0) call movieerr(52,status)
      status = nf_put_vara_double(ncid,presid,start2,count2,pres)
          if(status.ne.0) call movieerr(52,status)
      status = nf_put_vara_double(ncid,prespid,start2,count2,presp)
      status = nf_put_vara_double(ncid,chieid,start2,count2,chie)
      status = nf_put_vara_double(ncid,chiiid,start2,count2,chii)
      status = nf_put_vara_double(ncid,anhyid,start2,count2,anhy)
      status = nf_put_vara_double(ncid,anheid,start2,count2,anhe)
      status = nf_put_vara_double(ncid,anz1id,start2,count2,anz1)
      status = nf_put_vara_double(ncid,anz2id,start2,count2,anz2)
      status = nf_put_vara_double(ncid,anz3id,start2,count2,anz3)
      status = nf_put_vara_double(ncid,anz4id,start2,count2,anz4)
      status = nf_put_vara_double(ncid,anz5id,start2,count2,anz5)
      status = nf_put_vara_double(ncid,anz6id,start2,count2,anz6)
      status = nf_put_vara_double(ncid,anz7id,start2,count2,anz7)
      status = nf_put_vara_double(ncid,anz8id,start2,count2,anz8)
      status = nf_put_vara_double(ncid,ajtoid,start2,count2,ajto)
      status = nf_put_vara_double(ncid,ajbsid,start2,count2,ajbs)
      status = nf_put_vara_double(ncid,ajnbid,start2,count2,ajnb)
      status = nf_put_vara_double(ncid,ajlhid,start2,count2,ajlh)
      status = nf_put_vara_double(ncid,ajicid,start2,count2,ajic)
      status = nf_put_vara_double(ncid,ajecid,start2,count2,ajec)
      status = nf_put_vara_double(ncid,sabid,start2,count2,sab)
      status = nf_put_vara_double(ncid,henbid,start2,count2,henb)
      status = nf_put_vara_double(ncid,hinbid,start2,count2,hinb)
      status = nf_put_vara_double(ncid,heicid,start2,count2,heic)
      status = nf_put_vara_double(ncid,hiicid,start2,count2,hiic)
      status = nf_put_vara_double(ncid,helhid,start2,count2,helh)
      status = nf_put_vara_double(ncid,hilhid,start2,count2,hilh)
      status = nf_put_vara_double(ncid,heecid,start2,count2,heec)
      status = nf_put_vara_double(ncid,healid,start2,count2,heal)
      status = nf_put_vara_double(ncid,hialid,start2,count2,hial)
      status = nf_put_vara_double(ncid,hbreid,start2,count2,hbre)
      status = nf_put_vara_double(ncid,hcycid,start2,count2,hcyc)
      status = nf_put_vara_double(ncid,hlinid,start2,count2,hlin)
      status = nf_put_vara_double(ncid,vlid,start2,count2,vl)
      status = nf_put_vara_double(ncid,tflxid,start2,count2,tflx)
          if(status.ne.0) call movieerr(52,status)
!
      call xlimits(imin,imax)
      count2(1) = imax-imin+1
      do i=imin,imax
        ii = i + 1 - imin
        rcoord(ii) = xary(i)
        ps = psi(i,jmag)
        call peval(ps,2,pval,ppval,i,jmag)
        call eeval(ps,2,eval,epval,i,jmag)
        call reval(ps,idens,isurf,rval,rpval,i,jmag)
        call tieval(ps,2,tival,i,jmag,pval,eval,rval)
        ter(ii) = (eval*udsh)/(rval*udsd)
        tir(ii) = tival
        aner(ii) = rval*udsd
        ajphir(ii) = ajphi(i,jmag)*udsi
        qr(ii) =  qget(ps)
      enddo   
!
      status = nf_put_vara_double(ncid,rcoordid,start2,count2,rcoord)
          if(status.ne.0) call movieerr(43,status)
      status = nf_put_vara_double(ncid,terid,start2,count2,ter)
          if(status.ne.0) call movieerr(44,status)
      status = nf_put_vara_double(ncid,tirid,start2,count2,tir)
          if(status.ne.0) call movieerr(45,status)
      status = nf_put_vara_double(ncid,anerid,start2,count2,aner)
          if(status.ne.0) call movieerr(46,status)
      status = nf_put_vara_double(ncid,ajphirid,start2,count2,ajphir)
          if(status.ne.0) call movieerr(47,status)
      status = nf_put_vara_double(ncid,qrid,start2,count2,qr)
          if(status.ne.0) call movieerr(48,status)

!  Boundary points (by Kumar)

      status = nf_put_vara_double(ncid,rbndid,start2,count2n,rbnd)
      if(status.ne.0) call movieerr(43,status)
      status = nf_put_vara_double(ncid,zbndid,start2,count2n,zbnd)
      if(status.ne.0) call movieerr(44,status)
!
! 3. 2D Variables
!
!.....define psiaux array for plotting psi in plasma
      do i=2,nxp
        do j=2,nzp
          psiaux(i,j) = psi(i,j)
          psimaxx = psilim+eps
          if(phalo.ne.0) psimaxx = phalo+eps
!         if(iexv(i,j) .eq. 1 .or. iexs(i,j) .eq. 1) 
!    1               psiaux(i,j) = psimaxx
        enddo
      enddo
!
!.......................................
      call pack(psiaux,packcdf,isym,nx,nz,pnx,pnz)
      status = nf_put_vara_double(ncid,psiid,start3,count3,packcdf)
          if(status.ne.0) call movieerr(23,status)
      call packi(iexv,ipackcdf,isym,nx,nz,pnx,pnz)
      status = nf_put_vara_int(ncid,imaskid,start3,count3,ipackcdf)
          if(status.ne.0) call movieerr(25,status)
!
!
!.....define psiaux array for plotting divertor flux
      do i=2,nxp
        do j=2,nzp
          psiaux(i,j) = psi(i,j)
        enddo
      enddo
!
!....special diagnostics...............
      if(lrswtch.gt.0) then
        call oneplot(psiaux,psilimw,psisol,iminn,imaxx,jminn,jmaxx )
        write(nterm,6003) iminn,imaxx
 6003   format("iminn,imaxx =",2i5)
        write(nterm,6001) psilim,psimaxx,phalo,whalos,psimin
        write(nterm,6001) (psiaux(i,nh+2),i=2,nx)
        do j=2,nzp
          write(nterm,6002) (iexv(i,j),i=2,nx)
        enddo
 6002   format(0p55i2)
 6001   format(1p5e12.4)
      endif
!
      call pack(psiaux,packcdf,isym,nx,nz,pnx,pnz)
      status = nf_put_vara_double(ncid,divflxid,start3,count3,packcdf)
          if(status.ne.0) call movieerr(24,status)
!
      call pack(ajphi,packcdf,isym,nx,nz,pnx,pnz)
      status = nf_put_vara_double(ncid,ajphiid,start3,count3,packcdf)
          if(status.ne.0) call movieerr(80,status)
!
      call pack(pr,packcdf,isym,nx,nz,pnx,pnz)
      status = nf_put_vara_double(ncid,prid,start3,count3,packcdf)
          if(status.ne.0) call movieerr(81,status)
!
      call pack(gs,packcdf,isym,nx,nz,pnx,pnz)
      status = nf_put_vara_double(ncid,gsid,start3,count3,gs)
          if(status.ne.0) call movieerr(82,status)
!
      if(idata.eq.9 .or. idata.eq.11) then
        status = nf_put_var1_double(ncid,atotid,nreccdf,atot)
        status = nf_put_var1_double(ncid,csumid,nreccdf,csum)
        status = nf_put_var1_double(ncid,vesppeid,nreccdf,vesppe)
        status = nf_put_var1_double(ncid,vescureid,nreccdf,vescure)
        status = nf_put_var1_double(ncid,totcureid,nreccdf,totcure)
        count2(1) = numfluxloopsg
        status = nf_put_vara_double(ncid,pfluxid,start2,count2,pflux)
        status = nf_put_vara_double(ncid,efluxid,start2,count2,eflux)
        count2(1) = ncurrentg
        status = nf_put_vara_double(ncid,ppcurid,start2,count2,ppcur)
        status = nf_put_vara_double(ncid,fbcurid,start2,count2,fbcur)
        status = nf_put_vara_double(ncid,csumaid,start2,count2,csuma)
      endif
!
!.....end of writing variables
!
!
!
!
      return
!
      entry movieclose
      if(ncid .ne. 0) call cdf_close(ncid,status)
          if(status.ne.0) then
                    call movieerr(30,status)
                    write(*,*) "ncid", ncid
          endif
      return
      end
      subroutine pack(psi,packcdf,isym,nx,nz,pnx,pnz)
      IMPLICIT NONE
      integer :: pnx,pnz,isym,nx,nz, i,j,nzdimcdf
      real*8 :: psi, packcdf
      dimension psi(pnx,pnz), packcdf(*)
      nzdimcdf = nz + isym*(nz-1)      
!
!     pack variables 
      do i=2,nx+1
         if(isym.eq.0) then
             do j=1,nzdimcdf
              packcdf( (i-1)+ (j-1)*nx) = psi(i,j+1)
             enddo
         else
             do j=2,nz
              packcdf( (i-1)+ (j+nz-2)*nx) = psi(i,j+1)
              packcdf( (i-1)+ (nz-j  )*nx) = psi(i,j+1)
             enddo
              packcdf( (i-1)+ (nz-1)*nx)  = psi(i,2)
         endif
      enddo
!
      return
      end
      subroutine packi(imask,ipackcdf,isym,nx,nz,pnx,pnz)
      IMPLICIT NONE
      integer :: pnx,pnz, imask, ipackcdf,isym,nx,nz
      integer :: nzdimcdf, i,j
      dimension imask(pnx,pnz), ipackcdf(*)
      nzdimcdf = nz + isym*(nz-1)      
!
!     pack variables 
      do i=2,nx+1
         if(isym.eq.0) then
             do j=1,nzdimcdf
              ipackcdf( (i-1)+ (j-1)*nx) = imask(i,j+1)
             enddo
         else
             do j=2,nz
              ipackcdf( (i-1)+ (j+nz-2)*nx) = imask(i,j+1)
              ipackcdf( (i-1)+ (nz-j  )*nx) = imask(i,j+1)
             enddo
              ipackcdf( (i-1)+ (nz-1)*nx)  = imask(i,2)
         endif
      enddo
!
      return
      end

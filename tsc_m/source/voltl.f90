      function voltl(ipass,jpass)
!.....6.935 voltl
      USE CLINAM
      USE SAPROP
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER ipass,jpass,k,ks
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 voltl,ps,avm,avk
!============
      if(ipass.le.iminn .or. ipass.ge.imaxx) go to 200
      if(jpass.lt.jminn .or. jpass.gt.jmaxx) go to 200
      if(iexvc(ipass,jpass).gt.0 .or. iexs(ipass,jpass).eq.1) go to 200
      ps = psi(ipass,jpass)
      if(ps .ge. psilim .or. lrswtch.gt.0)                               &  
     &      go to 200
      if(isurf.le.0) go to 200
      do 240 k=3,npsit-1
      ks = k
      if(ps.lt.xsv2(k)) go  to 230
  240 continue
      go to 200
  230 continue
      avm = as(ks-1)*udsv
      avk = as(ks  )*udsv
      if(ps.gt.xsv2(2)) go to 250
      avm = 2._R8*avm - avk
      avk = avm
      ks = 2
  250 continue
      voltl = avm + (ps-xsv2(ks-1))*(avk-avm)/(xsv2(ks)-xsv2(ks-1))
      return
  200 continue
      voltl = 0
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

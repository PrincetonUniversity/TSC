      subroutine lookup2(jq,mte,mne,altei,alnei,alte,alne,alinzr,alradr  &  
     &,                  alrecr,nd1,nd2,nd3,ainzrjq,radrjq,recrjq)
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER mte,mne,nd1,nd2,nd3,jq,i,itep,item,inep,inem
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 altei,alnei,alte,alne,alinzr,alradr,alrecr,ainzrjq
      REAL*8 radrjq,recrjq,ff,gf
!============
      dimension altei(*), alnei(*),  alinzr(nd1,nd2)                     &  
     &,         alradr(nd1,nd2,nd3), alrecr(nd1,nd2,nd3)
!
      if(alte.lt.altei(1))   alte = altei(1)
      if(alte.gt.altei(mte)) alte = altei(mte)
 
      if(alne.lt.alnei(1))   alne = alnei(1)
      if(alne.gt.alnei(mne)) alne = alnei(mne)
 
      do 10 i = 1,mte
      itep = i
      item = i-1
      if(altei(i).gt.alte) go to 20
   10 continue
!
   20 ff = (alte-altei(item))                                            &  
     &              /(altei(itep)-altei(item))
!
      do 30 i = 1,mne
      inep = i
      inem = i-1
      if(alnei(i).gt.alne) go to 40
   30 continue
!
   40 gf = (alne-alnei(inem))                                            &  
     &              /(alnei(inep)-alnei(inem))
!
      ainzrjq = alinzr(jq,item)                                          &  
     &         +ff*(alinzr(jq,itep)-alinzr(jq,item))
 
      radrjq = (1._R8-ff)*(1._R8-gf)*alradr(jq,item,inem)                &  
     &       + (1._R8-ff)*    gf *alradr(jq,item,inep)                   &  
     &       +     ff *(1._R8-gf)*alradr(jq,itep,inem)                   &  
     &       +     ff *    gf *alradr(jq,itep,inep)
      recrjq = (1._R8-ff)*(1._R8-gf)*alrecr(jq,item,inem)                &  
     &       + (1._R8-ff)*    gf *alrecr(jq,item,inep)                   &  
     &       +     ff *(1._R8-gf)*alrecr(jq,itep,inem)                   &  
     &       +     ff *    gf *alrecr(jq,itep,inep)
 
!     write (nout,1000) recrjq,radrjq,ainzrjq
!1000 format (10x,"recrjq=",1pe11.4," radrjq=",1pe11.4," ainzejq=",
!    -1pe11.4," ainzrjq=",1pe11.4)
      recrjq = 10**recrjq
      radrjq = 10**radrjq
      ainzrjq = 10**ainzrjq
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

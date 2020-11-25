      subroutine lookup (nimp,jq,tekev,ane,ainzrjq,radrjq,recrjq)
!
      USE CLINAM
      USE RADTAB
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER jq,nimp,mte,mne
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 tekev,ane,ainzrjq,radrjq,recrjq,alte,alogxx,alne
!============
      alte = alogxx(tekev)
      alne = alogxx(ane)
      mte  = nte(nimp)
      mne  = nne(nimp)
!
      go to(10,20,30,40,50,60,70,80),nimp
   10 call lookup2(jq,mte,mne,altei(1,nimp),alnei(1,nimp),alte,alne      &  
     &,    alinzOx,alradOx,alrecOx,pcgOx,pte,pne,ainzrjq,radrjq,recrjq)
      return
!
   20 call lookup2(jq,mte,mne,altei(1,nimp),alnei(1,nimp),alte,alne      &  
     &,    alinzC,alradC,alrecC,pcgC,pte,pne,ainzrjq,radrjq,recrjq)
      return
!
   30 call lookup2(jq,mte,mne,altei(1,nimp),alnei(1,nimp),alte,alne      &  
     &,    alinzFe,alradFe,alrecFe,pcgFe,pte,pne,ainzrjq,radrjq,recrjq)
      return
!
   40 call lookup2(jq,mte,mne,altei(1,nimp),alnei(1,nimp),alte,alne      &  
     &,    alinzBe,alradBe,alrecBe,pcgBe,pte,pne,ainzrjq,radrjq,recrjq)
      return
!
   50 call lookup2(jq,mte,mne,altei(1,nimp),alnei(1,nimp),alte,alne      &  
     &,    alinzNe,alradNe,alrecNe,pcgNe,pte,pne,ainzrjq,radrjq,recrjq)
      return
!
   60 call lookup2(jq,mte,mne,altei(1,nimp),alnei(1,nimp),alte,alne      &  
     &,    alinzKr,alradKr,alrecKr,pcgKr,pte,pne,ainzrjq,radrjq,recrjq)
      return
!
   70 call lookup2(jq,mte,mne,altei(1,nimp),alnei(1,nimp),alte,alne      &  
     &,    alinzAr,alradAr,alrecAr,pcgAr,pte,pne,ainzrjq,radrjq,recrjq)
      return
!
   80 call lookup2(jq,mte,mne,altei(1,nimp),alnei(1,nimp),alte,alne      &  
     &,    alinzW,alradW,alrecW,pcgW,pte,pne,ainzrjq,radrjq,recrjq)
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

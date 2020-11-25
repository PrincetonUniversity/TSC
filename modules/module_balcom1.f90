      MODULE BALCOM1
      USE PARAM
      IMPLICIT NONE

      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::  bsqi,xbal,zbal,delt,deltp,  &  
     &                                 deltm,                            &  
     &                                 delp1,delm1,delp,del
      REAL*8, ALLOCATABLE, DIMENSION(:) ::  qpmh,fpe,fmh,fph,fm3,ppmh,   &  
     &                                 gmh,gpmh,                         &  
     &                           sum1,sum2,sum3,sum4,sum5,sum6


      END MODULE BALCOM1

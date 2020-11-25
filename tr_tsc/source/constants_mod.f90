        MODULE constants_mod
        USE spec_kind_mod
        IMPLICIT NONE


        
        REAL(kind=r8), PARAMETER      :: pi = 3.14159265358979_r8
        REAL(kind=r8), PARAMETER      :: mu0= pi *4E-7_r8
        REAL(kind=r8), PARAMETER      :: eCharge=1.602176487E-19_r8
        REAL(kind=r8), PARAMETER      :: vacuumPerm = 8.8452E-12_r8 
        REAL(kind=r8), PARAMETER      :: eV2Joule=eCharge
        REAL(kind=r8), PARAMETER      :: notgiven=1.234567890123456789E300_r8
        REAL(kind=r8), PARAMETER      :: zero=0.0_r8


        END MODULE constants_mod

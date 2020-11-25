      MODULE tgcblk

        LOGICAL ::  dotrans
        REAL ::  transx,transy
        REAL ::  scalex,scaley
        REAL ::  rotat
        REAL ::  centrx,centry
        REAL, DIMENSION(3,3) ::   trnmat
        LOGICAL :: matmade
        CHARACTER*8 :: tgname = "NONE    "
        LOGICAL :: doinit = .TRUE.


      END MODULE tgcblk

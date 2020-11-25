      MODULE tgchar
      
      REAL :: chx,chy,chrot
      REAL, PARAMETER :: chrht1 = .009 , chexp1 = 0.87
      REAL, PARAMETER :: chrht2 = .013 , chexp2 = 0.90 
      REAL, PARAMETER :: chrht3 = .016 , chexp3 = 0.98
      REAL, PARAMETER :: chrht4 = .018 , chexp4 = 1.28
      REAL, DIMENSION(4,4) :: chparm 
      DATA chparm                      /chrht1,.0156250 ,                &
     &                                  chexp1,.00781250 ,               &
     &                                  chrht2,.0235294 ,                &
     &                                  chexp2,.0117647 ,                &
     &                                  chrht3,.0312500 ,                &
     &                                  chexp3,.0156250 ,                &
     &                                  chrht4,.0476190 ,                &
     &                                  chexp4,.0238095 / 

      INTEGER :: ichcase,ichindx,ichfont,ichangle
      REAL :: chupx,chupy,chaddx,chaddy
      INTEGER, PARAMETER :: minfnt = 1, maxfnt = 20
      INTEGER :: ichclip
      INTEGER :: minfont = minfnt , maxfont = maxfnt
      LOGICAL :: autofeed

      END MODULE tgchar

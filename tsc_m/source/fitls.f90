      subroutine fitls(x,y,nrecord,aave,afin,a1,b1,a2,b2,lpts)
!     ----------------
!.....compute least square fits
!     linear regression is performed:  y = a*x + b
!     and the regression coefficients are calculated
!
!     a1,b1: final regression coeff. calculated over all points
!     a2,b2:                                    the last (lpts) points
!     aave(i): contains the regression coefficients (a) for: 1 <= ii <= i
!     afin(i):                                          i-lpts <= ii <= i
!
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER nrecord,lpts,i,ii
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 y,aave,afin,a1,b1,a2,b2,x,avex,avey,avexy,avexx,count
      REAL*8 yval,xval,denom,bvex,bvey,bvexy,bvexx
!============
      dimension x(*),y(*)
      dimension aave(*),afin(*)
!============      
!
!.....fit over all points (1 <= ii <=i):
      do 9006 i=1,nrecord
       aave(i) = 0._R8
       avex = 0._R8
       avey = 0._R8
       avexy = 0._R8
       avexx = 0._R8
       count = 0._R8
       if(i.le.10) go to 8006
       do 8010 ii=1,i
        yval = y(ii)
        xval = x(ii)
        avex = avex + xval
        avey = avey + yval
        avexy = avexy + xval*yval
        avexx = avexx + xval*xval
 8010  continue
       count = i
      if(count.le.0) go to 8006
       avex = avex/count
       avey = avey/count
       avexy = avexy/count
       avexx = avexx/count
       denom =(avexx-avex*avex)
      if(denom.eq.0) go to 8006
       aave(i) = (avexy-avex*avey)/denom
 8006  continue
!.....fit over last lpts points ( i-lpts <= ii <= i):
       afin(i) = 0._R8
       bvex = 0._R8
       bvey = 0._R8
       bvexy = 0._R8
       bvexx = 0._R8
       count = 0._R8
       if(i.le.lpts) go to 9006
       do 9010 ii=i-lpts,i
        count = count+1
        yval = y(ii)
        xval = x(ii)
        bvex = bvex + xval
        bvey = bvey + yval
        bvexy = bvexy + xval*yval
        bvexx = bvexx + xval*xval
 9010  continue
      if(count.le.0) go to 9006
       bvex = bvex/count
       bvey = bvey/count
       bvexy = bvexy/count
       bvexx = bvexx/count
      denom = bvexx-bvex*bvex
      if(denom.ne.0) afin(i) = (bvexy-bvex*bvey)/denom
!      afin(i) = (bvexy-bvex*bvey)/(bvexx-bvex*bvex)
 9006 continue
!
!.....construct least squares solution
!.....all points:
      a1 = aave(nrecord)
      b1 = avey - a1*avex
!.....last ten points:
      a2 = afin(nrecord)
      b2 = bvey - a2*bvex
!
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

      subroutine grap3(zz,xx,psval)
!.....number 8.60
!***************************************************************
!                                                              *
!...local bivariate interpolation                              *
!                                                              *
!***************************************************************
!
      USE CLINAM
      USE SCRATCH
      USE SAPROP
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER jval,ival,isave,jsave,ii,i,j,m,l,k
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 xx,psval,zz,vmat,fmat,poly,xpa,zpa,rhs1,rhs2,rhs3,rhs4
      REAL*8 xp,zp,sum1,f1
!============
      dimension vmat(4,4),fmat(0:3,0:3),poly(0:3,0:3)
      dimension xpa(0:3),zpa(0:3)
      data xpa(0),zpa(0)/ 1._R8, 1._R8/
!============      
!
      jval  = (zz-zzero)/deez+2
      ival  = (xx-ccon )/deex+2
!
      if(jval .lt. 2) jval = 2
      if(jval .gt. nz) jval = nz
      if(ival .lt. 2) ival = 2
      if(ival .gt. nx) ival = nx
      if(ival.eq.isave .and. jval.eq.jsave) go to 200
      isave = ival
      jsave = jval
!..........................................................
!
! calculate function values and function derivatives
! at four corners of the reference cell
!
!..........................................................
      do 10 ii=1,4
      go to(101,102,103,104),ii
101   i      = ival
      j      = jval
      go to 105
102   i      = ival
      j      = jval+1
      go to 105
103   i      = ival+1
      j      = jval
      go to 105
104   i      = ival+1
      j      = jval+1
105   continue
!
!...psi :
!
      vmat(1,ii) =      pso2(i,j)     - pso2(ival,jval)
      vmat(2,ii) = .50_R8*(pso2(i+1,j  ) - pso2(i-1,j  )   )
      vmat(3,ii) = .50_R8*(pso2(i  ,j+1) - pso2(i  ,j-1)   )
      vmat(4,ii) = .25_R8*(pso2(i+1,j+1) - pso2(i-1,j+1)                 &  
     &                  - pso2(i+1,j-1) + pso2(i-1,j-1)   )
!
10    continue
!..........................................................
!
!...determine coefficients by evaluating polynomial
!...at four corner points
!
!..........................................................
!
!...point (i,j) :
!
      fmat(0,0) = 0._R8
      fmat(1,0) = vmat(2,1)
      fmat(0,1) = vmat(3,1)
      fmat(1,1) = vmat(4,1)
!
!...point (i,j+1) :
!
      fmat(0,2) = 3._R8*vmat(1,2) - 2._R8*fmat(0,1)                      &  
     &            -    vmat(3,2)
      fmat(0,3) =    vmat(3,2) +    fmat(0,1)                            &  
     &            - 2._R8*vmat(1,2)
      fmat(1,2) = 3._R8*vmat(2,2) - 3._R8*fmat(1,0)                      &  
     &            -    vmat(4,2) - 2._R8*fmat(1,1)
      fmat(1,3) =    vmat(4,2) + 2._R8*fmat(1,0)                         &  
     &            - 2._R8*vmat(2,2) +    fmat(1,1)
!
!...point (i+1,j) :
!
      fmat(2,0) = 3._R8*vmat(1,3) - 2._R8*fmat(1,0)                      &  
     &            -    vmat(2,3)
      fmat(3,0) =    vmat(2,3) +    fmat(1,0)                            &  
     &            - 2._R8*vmat(1,3)
      fmat(2,1) = 3._R8*vmat(3,3) - 3._R8*fmat(0,1)                      &  
     &            -    vmat(4,3) - 2._R8*fmat(1,1)
      fmat(3,1) =    vmat(4,3) + 2._R8*fmat(0,1)                         &  
     &            - 2._R8*vmat(3,3) +    fmat(1,1)
!
!...point (i+1,j+1) :
!
      rhs1=   vmat(1,4)-vmat(1,3)-vmat(3,3)                              &  
     &    -  (fmat(0,2)+fmat(1,2)                                        &  
     &    +   fmat(0,3)+fmat(1,3) )
      rhs2=   vmat(2,4)-vmat(2,3)-vmat(4,3)                              &  
     &    -  (fmat(1,2)+fmat(1,3) )
      rhs3=   vmat(3,4)-vmat(3,3)                                        &  
     &    -2*(fmat(0,2)+fmat(1,2) )                                      &  
     &    -3*(fmat(0,3)+fmat(1,3) )
      rhs4=   vmat(4,4)-vmat(4,3)                                        &  
     &    -2*(fmat(1,2) )                                                &  
     &    -3*(fmat(1,3) )
!
      fmat(2,2) =  9*rhs1-3*rhs2-3*rhs3+rhs4
      fmat(3,2) = -6*rhs1+3*rhs2+2*rhs3-rhs4
      fmat(2,3) = -6*rhs1+2*rhs2+3*rhs3-rhs4
      fmat(3,3) =  4*rhs1-2*rhs2-2*rhs3+rhs4
  200 continue
!...........................................................
!
!.....evaluate function and derivatives at (xx,zz)
!
!...........................................................
      xp     = (xx-xary(ival))/deex
      zp     = (zz-zary(jval))/deez
!
      do 32 m=1,3
      xpa(m) = xp*xpa(m-1)
   32 zpa(m) = zp*zpa(m-1)
      do 33 l=0,3
      do 33 k=0,3
   33 poly(k,l) = xpa(k)*zpa(l)
!
!
      sum1   = 0._R8
      do 40 k = 0,3
      do 45 l = 0,3
      f1     = fmat(k,l)
      sum1   = sum1 + f1*poly(k,l)
   45 continue
   40 continue
      psval  = sum1 + pso2(ival,jval)
!
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

      subroutine svbksb(u,w,v,m,n,mp,np,b,x)
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER m,n,mp,np,nmax,j,i,jj
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 w,v,b,x,u,tmp,s
!============
      parameter (nmax=100)
      dimension u(mp,np),w(np),v(np,np),b(mp),x(np),tmp(nmax)
!============      
      do 12 j=1,n
        s=0._R8
        if(w(j).ne.0._R8)then
          do 11 i=1,m
            s=s+u(i,j)*b(i)
11        continue
          s=s/w(j)
        endif
        tmp(j)=s
12    continue
      do 14 j=1,n
        s=0._R8
        do 13 jj=1,n
          s=s+v(j,jj)*tmp(jj)
13      continue
        x(j)=s
14    continue
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

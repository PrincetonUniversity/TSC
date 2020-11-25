      REAL*8 function pythag(a,b)
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      REAL*8 a,b
 
!     finds sqrt(a**2+b**2) without overflow or destructive underflow
 
      REAL*8 p,r,s,t,u
!rew changed dmax1 to max
      p = max(abs(a),abs(b))
      if (p .eq. 0.000_R8) go to 20
!rew changed dmin1 to min
      r = (min(abs(a),abs(b))/p)**2
   10 continue
         t = 4.000_R8+ r
!        write(*,*) 't = ',t
         if (abs(t-4.000_R8) .lt. 1.E-5_R8) go to 20
         s = r / t
         u = 1.000_R8+ 2.000_R8*s
         p = u*p
         r = (s/u)**2 * r
      go to 10
   20 pythag = p
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

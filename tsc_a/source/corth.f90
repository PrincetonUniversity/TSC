      subroutine corth(nm,n,low,igh,ar,ai,ortr,orti)
 
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      integer i,j,m,n,ii,jj,la,mp,nm,igh,kp1,low
      REAL*8 ar(nm,n),ai(nm,n),ortr(igh),orti(igh)
      REAL*8 f,g,h,fi,fr,scale,pythag,dlapy3gf
 
!     this subroutine is a translation of a complex analogue of
!     the algol procedure orthes, num. math. 12, 349-368(1968)
!     by martin and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 339-358(1971).
 
!     given a complex general matrix, this subroutine
!     reduces a submatrix situated in rows and columns
!     low through igh to upper hessenberg form by
!     unitary similarity transformations.
 
!     on input
 
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
 
!        n is the order of the matrix.
 
!        low and igh are integers determined by the balancing
!          subroutine  cbal.  if  cbal  has not been used,
!          set low=1, igh=n.
 
!        ar and ai contain the real and imaginary parts,
!          respectively, of the complex input matrix.
 
!     on output
 
!        ar and ai contain the real and imaginary parts,
!          respectively, of the hessenberg matrix.  information
!          about the unitary transformations used in the reduction
!          is stored in the remaining triangles under the
!          hessenberg matrix.
 
!        ortr and orti contain further information about the
!          transformations.  only elements low through igh are used.
 
!     calls pythag for  sqrt(a*a + b*b) .
 
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
 
!     this version dated august 1983.
 
!     ------------------------------------------------------------------
 
      la = igh - 1
      kp1 = low + 1
      if (la .lt. kp1) go to 200
 
      do 180 m = kp1, la
         h = 0.000_R8
         ortr(m) = 0.000_R8
         orti(m) = 0.000_R8
         scale = 0.000_R8
!     .......... scale column (algol tol then not needed) ..........
         do 90 i = m, igh
   90    scale = scale + abs(ar(i,m-1)) + abs(ai(i,m-1))
 
         if (scale .eq. 0.000_R8) go to 180
         mp = m + igh
!     .......... for i=igh step -1 until m do -- ..........
         do 100 ii = m, igh
            i = mp - ii
            ortr(i) = ar(i,m-1) / scale
            orti(i) = ai(i,m-1) / scale
            h = h + ortr(i) * ortr(i) + orti(i) * orti(i)
  100    continue
 
         g = sqrt(h)
         f = dlapy3gf(ortr(m),orti(m))
         if (f .eq. 0.000_R8) go to 103
         h = h + f * g
         g = g / f
         ortr(m) = (1.000_R8+ g) * ortr(m)
         orti(m) = (1.000_R8+ g) * orti(m)
         go to 105
 
  103    ortr(m) = g
         ar(m,m-1) = scale
!     .......... form (i-(u*ut)/h) * a ..........
  105    do 130 j = m, n
            fr = 0.000_R8
            fi = 0.000_R8
!     .......... for i=igh step -1 until m do -- ..........
            do 110 ii = m, igh
               i = mp - ii
               fr = fr + ortr(i) * ar(i,j) + orti(i) * ai(i,j)
               fi = fi + ortr(i) * ai(i,j) - orti(i) * ar(i,j)
  110       continue
 
            fr = fr / h
            fi = fi / h
 
            do 120 i = m, igh
               ar(i,j) = ar(i,j) - fr * ortr(i) + fi * orti(i)
               ai(i,j) = ai(i,j) - fr * orti(i) - fi * ortr(i)
  120       continue
 
  130    continue
!     .......... form (i-(u*ut)/h)*a*(i-(u*ut)/h) ..........
         do 160 i = 1, igh
            fr = 0.000_R8
            fi = 0.000_R8
!     .......... for j=igh step -1 until m do -- ..........
            do 140 jj = m, igh
               j = mp - jj
               fr = fr + ortr(j) * ar(i,j) - orti(j) * ai(i,j)
               fi = fi + ortr(j) * ai(i,j) + orti(j) * ar(i,j)
  140       continue
 
            fr = fr / h
            fi = fi / h
 
            do 150 j = m, igh
               ar(i,j) = ar(i,j) - fr * ortr(j) - fi * orti(j)
               ai(i,j) = ai(i,j) + fr * orti(j) - fi * ortr(j)
  150       continue
 
  160    continue
 
         ortr(m) = scale * ortr(m)
         orti(m) = scale * orti(m)
         ar(m,m-1) = -g * ar(m,m-1)
         ai(m,m-1) = -g * ai(m,m-1)
  180 continue
 
  200 return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

      subroutine comqr(nm,n,low,igh,hr,hi,wr,wi,ierr)
 
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      integer i,j,l,n,en,ll,nm,igh,itn,its,low,lp1,enm1,ierr
      REAL*8 hr(nm,n),hi(nm,n),wr(n),wi(n)
      REAL*8 si,sr,ti,tr,xi,xr,yi,yr,zzi,zzr,norm,tst1,tst2,             &  
     &       pythag,dlapy3gf
 
!     this subroutine is a translation of a unitary analogue of the
!     algol procedure  comlr, num. math. 12, 369-376(1968) by martin
!     and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 396-403(1971).
!     the unitary analogue substitutes the qr algorithm of francis
!     (comp. jour. 4, 332-345(1962)) for the lr algorithm.
 
!     this subroutine finds the eigenvalues of a complex
!     upper hessenberg matrix by the qr method.
 
!     on input
 
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
 
!        n is the order of the matrix.
 
!        low and igh are integers determined by the balancing
!          subroutine  cbal.  if  cbal  has not been used,
!          set low=1, igh=n.
 
!        hr and hi contain the real and imaginary parts,
!          respectively, of the complex upper hessenberg matrix.
!          their lower triangles below the subdiagonal contain
!          information about the unitary transformations used in
!          the reduction by  corth, if performed.
 
!     on output
 
!        the upper hessenberg portions of hr and hi have been
!          destroyed.  therefore, they must be saved before
!          calling  comqr  if subsequent calculation of
!          eigenvectors is to be performed.
 
!        wr and wi contain the real and imaginary parts,
!          respectively, of the eigenvalues.  if an error
!          exit is made, the eigenvalues should be correct
!          for indices ierr+1,...,n.
 
!        ierr is set to
!          zero       for normal return,
!          j          if the limit of 30*n iterations is exhausted
!                     while the j-th eigenvalue is being sought.
 
!     calls cdiv for complex division.
!     calls csroot for complex square root.
!     calls pythag for  sqrt(a*a + b*b) .
 
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
 
!     this version dated august 1983.
 
!     ------------------------------------------------------------------
 
      ierr = 0
      if (low .eq. igh) go to 180
!     .......... create real subdiagonal elements ..........
      l = low + 1
 
      do 170 i = l, igh
         ll = min(i+1,igh)
         if (hi(i,i-1) .eq. 0.000_R8) go to 170
         norm = dlapy3gf(hr(i,i-1),hi(i,i-1))
!rew inserted norm+1.d-100
         yr = hr(i,i-1) / (norm+1.E-100_R8)
         yi = hi(i,i-1) / (norm+1.E-100_R8)
         hr(i,i-1) = norm
         hi(i,i-1) = 0.000_R8
 
         do 155 j = i, igh
            si = yr * hi(i,j) - yi * hr(i,j)
            hr(i,j) = yr * hr(i,j) + yi * hi(i,j)
            hi(i,j) = si
  155    continue
 
         do 160 j = low, ll
            si = yr * hi(j,i) + yi * hr(j,i)
            hr(j,i) = yr * hr(j,i) - yi * hi(j,i)
            hi(j,i) = si
  160    continue
 
  170 continue
!     .......... store roots isolated by cbal ..........
  180 do 200 i = 1, n
         if (i .ge. low .and. i .le. igh) go to 200
         wr(i) = hr(i,i)
         wi(i) = hi(i,i)
  200 continue
 
      en = igh
      tr = 0.000_R8
      ti = 0.000_R8
      itn = 30*n
!     .......... search for next eigenvalue ..........
  220 if (en .lt. low) go to 1001
      its = 0
      enm1 = en - 1
!     .......... look for single small sub-diagonal element
!                for l=en step -1 until low d0 -- ..........
  240 do 260 ll = low, en
         l = en + low - ll
         if (l .eq. low) go to 300
         tst1 = abs(hr(l-1,l-1)) + abs(hi(l-1,l-1))                      &  
     &            + abs(hr(l,l)) + abs(hi(l,l))
         tst2 = tst1 + abs(hr(l,l-1))
         if (tst2 .eq. tst1) go to 300
  260 continue
!     .......... form shift ..........
  300 if (l .eq. en) go to 660
      if (itn .eq. 0) go to 1000
      if (its .eq. 10 .or. its .eq. 20) go to 320
      sr = hr(en,en)
      si = hi(en,en)
      xr = hr(enm1,en) * hr(en,enm1)
      xi = hi(enm1,en) * hr(en,enm1)
      if (xr .eq. 0.000_R8.and. xi .eq. 0.000_R8) go to 340
      yr = (hr(enm1,enm1) - sr) / 2.000_R8
      yi = (hi(enm1,enm1) - si) / 2.000_R8
      call csroot(yr**2-yi**2+xr,2.000_R8*yr*yi+xi,zzr,zzi)
      if (yr * zzr + yi * zzi .ge. 0.000_R8) go to 310
      zzr = -zzr
      zzi = -zzi
  310 call cdiv(xr,xi,yr+zzr,yi+zzi,xr,xi)
      sr = sr - xr
      si = si - xi
      go to 340
!     .......... form exceptional shift ..........
  320 sr = abs(hr(en,enm1)) + abs(hr(enm1,en-2))
      si = 0.000_R8
 
  340 do 360 i = low, en
         hr(i,i) = hr(i,i) - sr
         hi(i,i) = hi(i,i) - si
  360 continue
 
      tr = tr + sr
      ti = ti + si
      its = its + 1
      itn = itn - 1
!     .......... reduce to triangle (rows) ..........
      lp1 = l + 1
 
      do 500 i = lp1, en
         sr = hr(i,i-1)
         hr(i,i-1) = 0.000_R8
         norm = dlapy3gf(dlapy3gf(hr(i-1,i-1),hi(i-1,i-1)),sr)
!rew inserted norm+1.d-100
         xr = hr(i-1,i-1) / (norm+1.E-100_R8)
         wr(i-1) = xr
         xi = hi(i-1,i-1) / (norm+1.E-100_R8)
         wi(i-1) = xi
         hr(i-1,i-1) = norm
         hi(i-1,i-1) = 0.000_R8
         hi(i,i-1) = sr / (norm+1.E-100_R8)
 
         do 490 j = i, en
            yr = hr(i-1,j)
            yi = hi(i-1,j)
            zzr = hr(i,j)
            zzi = hi(i,j)
            hr(i-1,j) = xr * yr + xi * yi + hi(i,i-1) * zzr
            hi(i-1,j) = xr * yi - xi * yr + hi(i,i-1) * zzi
            hr(i,j) = xr * zzr - xi * zzi - hi(i,i-1) * yr
            hi(i,j) = xr * zzi + xi * zzr - hi(i,i-1) * yi
  490    continue
 
  500 continue
 
      si = hi(en,en)
      if (si .eq. 0.000_R8) go to 540
      norm = dlapy3gf(hr(en,en),si)
!rew inserted norm+1.d-100
      sr = hr(en,en) / (norm+1.E-100_R8)
      si = si / (norm+1.E-100_R8)
      hr(en,en) = norm
      hi(en,en) = 0.000_R8
!     .......... inverse operation (columns) ..........
  540 do 600 j = lp1, en
         xr = wr(j-1)
         xi = wi(j-1)
 
         do 580 i = l, j
            yr = hr(i,j-1)
            yi = 0.000_R8
            zzr = hr(i,j)
            zzi = hi(i,j)
            if (i .eq. j) go to 560
            yi = hi(i,j-1)
            hi(i,j-1) = xr * yi + xi * yr + hi(j,j-1) * zzi
  560       hr(i,j-1) = xr * yr - xi * yi + hi(j,j-1) * zzr
            hr(i,j) = xr * zzr + xi * zzi - hi(j,j-1) * yr
            hi(i,j) = xr * zzi - xi * zzr - hi(j,j-1) * yi
  580    continue
 
  600 continue
 
      if (si .eq. 0.000_R8) go to 240
 
      do 630 i = l, en
         yr = hr(i,en)
         yi = hi(i,en)
         hr(i,en) = sr * yr - si * yi
         hi(i,en) = sr * yi + si * yr
  630 continue
 
      go to 240
!     .......... a root found ..........
  660 wr(en) = hr(en,en) + tr
      wi(en) = hi(en,en) + ti
      en = enm1
      go to 220
!     .......... set error -- all eigenvalues have not
!                converged after 30*n iterations ..........
 1000 ierr = en
 1001 return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

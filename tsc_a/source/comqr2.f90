      subroutine comqr2(nm,n,low,igh,ortr,orti,hr,hi,wr,wi,zr,zi,ierr)
!  MESHED overflow control WITH vectors of isolated roots (10/19/89 BSG)
!  MESHED overflow control WITH triangular multiply (10/30/89 BSG)
 
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      integer i,j,k,l,m,n,en,ii,jj,ll,nm,nn,igh,ip1,                     &  
     &        itn,its,low,lp1,enm1,iend,ierr
      REAL*8 hr(nm,n),hi(nm,n),wr(n),wi(n),zr(nm,n),zi(nm,n),            &  
     &       ortr(igh),orti(igh)
      REAL*8 si,sr,ti,tr,xi,xr,yi,yr,zzi,zzr,norm,tst1,tst2,             &  
     &       pythag, dlapy3gf
 
!     this subroutine is a translation of a unitary analogue of the
!     algol procedure  comlr2, num. math. 16, 181-204(1970) by peters
!     and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 372-395(1971).
!     the unitary analogue substitutes the qr algorithm of francis
!     (comp. jour. 4, 332-345(1962)) for the lr algorithm.
 
!     this subroutine finds the eigenvalues and eigenvectors
!     of a complex upper hessenberg matrix by the qr
!     method.  the eigenvectors of a complex general matrix
!     can also be found if  corth  has been used to reduce
!     this general matrix to hessenberg form.
 
!     on input
 
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
 
!        n is the order of the matrix.
 
!        low and igh are integers determined by the balancing
!          subroutine  cbal.  if  cbal  has not been used,
!          set low=1, igh=n.
 
!        ortr and orti contain information about the unitary trans-
!          formations used in the reduction by  corth, if performed.
!          only elements low through igh are used.  if the eigenvectors
!          of the hessenberg matrix are desired, set ortr(j) and
!          orti(j) to 0.000 for these elements.
 
!        hr and hi contain the real and imaginary parts,
!          respectively, of the complex upper hessenberg matrix.
!          their lower triangles below the subdiagonal contain further
!          information about the transformations which were used in the
!          reduction by  corth, if performed.  if the eigenvectors of
!          the hessenberg matrix are desired, these elements may be
!          arbitrary.
 
!     on output
 
!        ortr, orti, and the upper hessenberg portions of hr and hi
!          have been destroyed.
 
!        wr and wi contain the real and imaginary parts,
!          respectively, of the eigenvalues.  if an error
!          exit is made, the eigenvalues should be correct
!          for indices ierr+1,...,n.
 
!        zr and zi contain the real and imaginary parts,
!          respectively, of the eigenvectors.  the eigenvectors
!          are unnormalized.  if an error exit is made, none of
!          the eigenvectors has been found.
 
!        ierr is set to
!          zero       for normal return,
!          j          if the limit of 30*n iterations is exhausted
!                     while the j-th eigenvalue is being sought.
 
!     calls cdiv for complex division.
!     calls csroot for complex square root.
!     calls pythag for  sqrt(a*a + b*b) .
 
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
 
!     this version dated october 1989.
 
!     ------------------------------------------------------------------
 
      ierr = 0
!     .......... initialize eigenvector matrix ..........
      do 101 j = 1, n
 
         do 100 i = 1, n
            zr(i,j) = 0.000_R8
            zi(i,j) = 0.000_R8
  100    continue
         zr(j,j) = 1.000_R8
  101 continue
!     .......... form the matrix of accumulated transformations
!                from the information left by corth ..........
      iend = igh - low - 1
      if (iend) 180, 150, 105
!     .......... for i=igh-1 step -1 until low+1 do -- ..........
  105 do 140 ii = 1, iend
         i = igh - ii
         if (ortr(i) .eq. 0.000_R8.and. orti(i) .eq. 0.000_R8) go to     &  
     & 140
         if (hr(i,i-1) .eq. 0.000_R8.and. hi(i,i-1) .eq. 0.000_R8) go    &  
     & to 140
!     .......... norm below is negative of h formed in corth ..........
         norm = hr(i,i-1) * ortr(i) + hi(i,i-1) * orti(i)
         ip1 = i + 1
 
         do 110 k = ip1, igh
            ortr(k) = hr(k,i-1)
            orti(k) = hi(k,i-1)
  110    continue
 
         do 130 j = i, igh
            sr = 0.000_R8
            si = 0.000_R8
            do 115 k = i, igh
               sr = sr + ortr(k) * zr(k,j) + orti(k) * zi(k,j)
               si = si + ortr(k) * zi(k,j) - orti(k) * zr(k,j)
  115       continue
!
!rew inserted norm+1.d-100
            sr = sr / (norm+1.E-100_R8)
            si = si / (norm+1.E-100_R8)
 
            do 120 k = i, igh
               zr(k,j) = zr(k,j) + sr * ortr(k) - si * orti(k)
               zi(k,j) = zi(k,j) + sr * orti(k) + si * ortr(k)
  120       continue
 
  130    continue
 
  140 continue
!     .......... create real subdiagonal elements ..........
  150 l = low + 1
 
      do 170 i = l, igh
         ll = min(i+1,igh)
         if (hi(i,i-1) .eq. 0.000_R8) go to 170
         norm = dlapy3gf(hr(i,i-1),hi(i,i-1))
!rew     inserted norm+1.d-100
         yr = hr(i,i-1) / (norm+1.E-100_R8)
         yi = hi(i,i-1) / (norm+1.E-100_R8)
         hr(i,i-1) = norm
         hi(i,i-1) = 0.000_R8
 
         do 155 j = i, n
            si = yr * hi(i,j) - yi * hr(i,j)
            hr(i,j) = yr * hr(i,j) + yi * hi(i,j)
            hi(i,j) = si
  155    continue
 
         do 160 j = 1, ll
            si = yr * hi(j,i) + yi * hr(j,i)
            hr(j,i) = yr * hr(j,i) - yi * hi(j,i)
            hi(j,i) = si
  160    continue
 
         do 165 j = low, igh
            si = yr * zi(j,i) + yi * zr(j,i)
            zr(j,i) = yr * zr(j,i) - yi * zi(j,i)
            zi(j,i) = si
  165    continue
 
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
  220 if (en .lt. low) go to 680
      its = 0
      enm1 = en - 1
!     .......... look for single small sub-diagonal element
!                for l=en step -1 until low do -- ..........
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
 
         do 490 j = i, n
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
      if (en .eq. n) go to 540
      ip1 = en + 1
 
      do 520 j = ip1, n
         yr = hr(en,j)
         yi = hi(en,j)
         hr(en,j) = sr * yr + si * yi
         hi(en,j) = sr * yi - si * yr
  520 continue
!     .......... inverse operation (columns) ..........
  540 do 600 j = lp1, en
         xr = wr(j-1)
         xi = wi(j-1)
 
         do 580 i = 1, j
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
 
         do 590 i = low, igh
            yr = zr(i,j-1)
            yi = zi(i,j-1)
            zzr = zr(i,j)
            zzi = zi(i,j)
            zr(i,j-1) = xr * yr - xi * yi + hi(j,j-1) * zzr
            zi(i,j-1) = xr * yi + xi * yr + hi(j,j-1) * zzi
            zr(i,j) = xr * zzr + xi * zzi - hi(j,j-1) * yr
            zi(i,j) = xr * zzi - xi * zzr - hi(j,j-1) * yi
  590    continue
 
  600 continue
 
      if (si .eq. 0.000_R8) go to 240
 
      do 630 i = 1, en
         yr = hr(i,en)
         yi = hi(i,en)
         hr(i,en) = sr * yr - si * yi
         hi(i,en) = sr * yi + si * yr
  630 continue
 
      do 640 i = low, igh
         yr = zr(i,en)
         yi = zi(i,en)
         zr(i,en) = sr * yr - si * yi
         zi(i,en) = sr * yi + si * yr
  640 continue
 
      go to 240
!     .......... a root found ..........
  660 hr(en,en) = hr(en,en) + tr
      wr(en) = hr(en,en)
      hi(en,en) = hi(en,en) + ti
      wi(en) = hi(en,en)
      en = enm1
      go to 220
!     .......... all roots found.  backsubstitute to find
!                vectors of upper triangular form ..........
  680 norm = 0.000_R8
 
      do 720 i = 1, n
 
         do 720 j = i, n
            tr = abs(hr(i,j)) + abs(hi(i,j))
            if (tr .gt. norm) norm = tr
  720 continue
 
      if (n .eq. 1 .or. norm .eq. 0.000_R8) go to 1001
!     .......... for en=n step -1 until 2 do -- ..........
      do 800 nn = 2, n
         en = n + 2 - nn
         xr = wr(en)
         xi = wi(en)
         hr(en,en) = 1.000_R8
         hi(en,en) = 0.000_R8
         enm1 = en - 1
!     .......... for i=en-1 step -1 until 1 do -- ..........
         do 780 ii = 1, enm1
            i = en - ii
            zzr = 0.000_R8
            zzi = 0.000_R8
            ip1 = i + 1
 
            do 740 j = ip1, en
               zzr = zzr + hr(i,j) * hr(j,en) - hi(i,j) * hi(j,en)
               zzi = zzi + hr(i,j) * hi(j,en) + hi(i,j) * hr(j,en)
  740       continue
 
            yr = xr - wr(i)
            yi = xi - wi(i)
            if (yr .ne. 0.000_R8.or. yi .ne. 0.000_R8) go to 765
               tst1 = norm
               yr = tst1
  760          yr = 0.01_R8* yr
               tst2 = norm + yr
               if (tst2 .gt. tst1) go to 760
  765       continue
            call cdiv(zzr,zzi,yr,yi,hr(i,en),hi(i,en))
!     .......... overflow control ..........
            tr = abs(hr(i,en)) + abs(hi(i,en))
            if (tr .eq. 0.000_R8) go to 780
            tst1 = tr
            tst2 = tst1 + 1.000_R8/tst1
            if (tst2 .gt. tst1) go to 780
            do 770 j = i, en
               hr(j,en) = hr(j,en)/tr
               hi(j,en) = hi(j,en)/tr
  770       continue
 
  780    continue
 
  800 continue
!     .......... end backsubstitution ..........
!     .......... vectors of isolated roots ..........
      do  840 i = 1, N
         if (i .ge. low .and. i .le. igh) go to 840
 
         do 820 j = I, n
            zr(i,j) = hr(i,j)
            zi(i,j) = hi(i,j)
  820    continue
 
  840 continue
!     .......... multiply by transformation matrix to give
!                vectors of original full matrix.
!                for j=n step -1 until low do -- ..........
      do 880 jj = low, N
         j = n + low - jj
         m = min(j,igh)
 
         do 880 i = low, igh
            zzr = 0.000_R8
            zzi = 0.000_R8
 
            do 860 k = low, m
               zzr = zzr + zr(i,k) * hr(k,j) - zi(i,k) * hi(k,j)
               zzi = zzi + zr(i,k) * hi(k,j) + zi(i,k) * hr(k,j)
  860       continue
 
            zr(i,j) = zzr
            zi(i,j) = zzi
  880 continue
 
      go to 1001
!     .......... set error -- all eigenvalues have not
!                converged after 30*n iterations ..........
 1000 ierr = en
 1001 return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

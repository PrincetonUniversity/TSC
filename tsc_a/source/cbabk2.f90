      subroutine cbabk2(nm,n,low,igh,scale,m,zr,zi)
 
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      integer i,j,k,m,n,ii,nm,igh,low
      REAL*8 scale(n),zr(nm,m),zi(nm,m)
      REAL*8 s
 
!     this subroutine is a translation of the algol procedure
!     cbabk2, which is a complex version of balbak,
!     num. math. 13, 293-304(1969) by parlett and reinsch.
!     handbook for auto. comp., vol.ii-linear algebra, 315-326(1971).
 
!     this subroutine forms the eigenvectors of a complex general
!     matrix by back transforming those of the corresponding
!     balanced matrix determined by  cbal.
 
!     on input
 
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
 
!        n is the order of the matrix.
 
!        low and igh are integers determined by  cbal.
 
!        scale contains information determining the permutations
!          and scaling factors used by  cbal.
 
!        m is the number of eigenvectors to be back transformed.
 
!        zr and zi contain the real and imaginary parts,
!          respectively, of the eigenvectors to be
!          back transformed in their first m columns.
 
!     on output
 
!        zr and zi contain the real and imaginary parts,
!          respectively, of the transformed eigenvectors
!          in their first m columns.
 
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
 
!     this version dated august 1983.
 
!     ------------------------------------------------------------------
 
      if (m .eq. 0) go to 200
      if (igh .eq. low) go to 120
 
      do 110 i = low, igh
         s = scale(i)
!     .......... left hand eigenvectors are back transformed
!                if the foregoing statement is replaced by
!                s=1.000/scale(i). ..........
         do 100 j = 1, m
            zr(i,j) = zr(i,j) * s
            zi(i,j) = zi(i,j) * s
  100    continue
 
  110 continue
!     .......... for i=low-1 step -1 until 1,
!                igh+1 step 1 until n do -- ..........
  120 do 140 ii = 1, n
         i = ii
         if (i .ge. low .and. i .le. igh) go to 140
         if (i .lt. low) i = low - ii
         k = scale(i)
         if (k .eq. i) go to 140
 
         do 130 j = 1, m
            s = zr(i,j)
            zr(i,j) = zr(k,j)
            zr(k,j) = s
            s = zi(i,j)
            zi(i,j) = zi(k,j)
            zi(k,j) = s
  130    continue
 
  140 continue
 
  200 return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

      subroutine cgg_glf(nm,n,ar,ai,wr,wi,matz,zr,zi,fv1,fv2,fv3,ierr)
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cgg
!---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
!
 
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      integer n,nm,is1,is2,ierr,matz
      REAL*8 ar(nm,n),ai(nm,n),wr(n),wi(n),zr(nm,n),zi(nm,n),            &  
     &       fv1(n),fv2(n),fv3(n)
 
!     this subroutine calls the recommended sequence of
!     subroutines from the eigensystem subroutine package (eispack)
!     to find the eigenvalues and eigenvectors (if desired)
!     of a complex general matrix.
 
!     on input
 
!        nm  must be set to the row dimension of the two-dimensional
!        array parameters as declared in the calling program
!        dimension statement.
 
!        n  is the order of the matrix  a=(ar,ai).
 
!        ar  and  ai  contain the real and imaginary parts,
!        respectively, of the complex general matrix.
 
!        matz  is an integer variable set equal to zero if
!        only eigenvalues are desired.  otherwise it is set to
!        any non-zero integer for both eigenvalues and eigenvectors.
 
!     on output
 
!        wr  and  wi  contain the real and imaginary parts,
!        respectively, of the eigenvalues.
 
!        zr  and  zi  contain the real and imaginary parts,
!        respectively, of the eigenvectors if matz is not zero.
 
!        ierr  is an integer output variable set equal to an error
!           completion code described in the documentation for comqr
!           and comqr2.  the normal completion code is zero.
 
!        fv1, fv2, and  fv3  are temporary storage arrays.
 
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
 
!     this version dated august 1983.
 
!     ------------------------------------------------------------------
 
      if (n .le. nm) go to 10
      ierr = 10 * n
      go to 50
 
   10 call  cbal(nm,n,ar,ai,is1,is2,fv1)
      call  corth(nm,n,is1,is2,ar,ai,fv2,fv3)
      if (matz .ne. 0) go to 20
!     .......... find eigenvalues only ..........
      call  comqr(nm,n,is1,is2,ar,ai,wr,wi,ierr)
      go to 50
!     .......... find both eigenvalues and eigenvectors ..........
   20 call  comqr2(nm,n,is1,is2,fv2,fv3,ar,ai,wr,wi,zr,zi,ierr)
      if (ierr .ne. 0) go to 50
      call  cbabk2(nm,n,is1,is2,fv1,n,zr,zi)
   50 return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

      subroutine pack(psi,packcdf,isym,nx,nz,pnx,pnz)
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER isym,nx,nz,nzdimcdf,i,j
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 packcdf,psi
!============
      integer pnx,pnz
      dimension psi(pnx,pnz), packcdf(*)
!============      
      nzdimcdf = nz + isym*(nz-1)
!
!     pack variables
      do i=2,nx+1
         if(isym.eq.0) then
             do j=1,nzdimcdf
              packcdf( (i-1)+ (j-1)*nx) = psi(i,j+1)
             enddo
         else
             do j=2,nz
              packcdf( (i-1)+ (j+nz-2)*nx) = psi(i,j+1)
              packcdf( (i-1)+ (nz-j  )*nx) = psi(i,j+1)
             enddo
              packcdf( (i-1)+ (nz-1)*nx)  = psi(i,2)
         endif
      enddo
!
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

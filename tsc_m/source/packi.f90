      subroutine packi(imask,ipackcdf,isym,nx,nz,pnx,pnz)
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER ipackcdf,isym,nx,nz,imask,nzdimcdf,i,j
!============
      integer pnx,pnz
      dimension imask(pnx,pnz), ipackcdf(*)
!============      
      nzdimcdf = nz + isym*(nz-1)
!
!     pack variables
      do i=2,nx+1
         if(isym.eq.0) then
             do j=1,nzdimcdf
              ipackcdf( (i-1)+ (j-1)*nx) = imask(i,j+1)
             enddo
         else
             do j=2,nz
              ipackcdf( (i-1)+ (j+nz-2)*nx) = imask(i,j+1)
              ipackcdf( (i-1)+ (nz-j  )*nx) = imask(i,j+1)
             enddo
              ipackcdf( (i-1)+ (nz-1)*nx)  = imask(i,2)
         endif
      enddo
!
      return
      end
! 15Apr2005 fgtok -s r8_precision.sub "r8con.csh conversion"

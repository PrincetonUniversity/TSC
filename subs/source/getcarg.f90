      subroutine getcarg(index, arg, numargs)
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer :: index, numargs
      character*(*), intent(out) :: arg
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: numchars, ier
!-----------------------------------------------
      integer iargc, getarg
      numargs = iargc()
      numchars = getarg(index, arg)

      end subroutine getcarg


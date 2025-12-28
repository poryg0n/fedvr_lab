      program reading
      implicit none
      integer :: i, j, s, ut0, nmax
      real*8 :: aux
      real*8, dimension(:), allocatable :: wfc

      include 'param'

      nmax = lnbr*(nnbr-1)-1
      ut0 = 15
      allocate(wfc(5))
!     open(ut0,file='data/basis_vect.dat', status='unknown') 
      open(ut0,file='dummyfile.dat', status='unknown') 

      read(ut0,*) s, i, wfc

      write(*,*) wfc

      end 

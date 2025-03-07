program test
  use sphpack
  use omp_lib
  implicit none
  
  integer, parameter :: jcut = 49
  
  integer, parameter :: jmcut  = jcut*(jcut+1)/2+jcut+1
  
  real(kind=dbl)                 :: start, end
  real(kind=dbl),    pointer     :: grid(:,:,:)
  type(c_ptr)                    :: c_grid
  complex(kind=dbl), allocatable :: c1(:), c2(:), cout(:)
  
  !*****************************************************************************!
  !** Init sphpack  ************************************************************!
  !*****************************************************************************!
  call init_sphpack( jcut )
  
  !*****************************************************************************!
  !** harmsy/harman test  ******************************************************!
  !*****************************************************************************!
  allocate( c1(jmcut)   ) ; call fill_scalar_sub( c1 )
  allocate( cout(jmcut) ) ; cout = cmplx(0._dbl, 0._dbl, kind=dbl)
  
  c_grid = malloc( 32, int(2 * nth * nph * c_sizeof(0._dbl), kind=4) )
  call c_f_pointer( c_grid,  grid,  [nth,2,nph] )
  
  start = omp_get_wtime()
    call harmsy( c1, grid )
    call harman( grid, cout )
  end = omp_get_wtime()
  
  call free( c_grid)
  
  write(*,*) 'harmsy/harman:'
  write(*,*) 'exec time: ',end-start
  write(*,*) 'supreme error: ', maxval( abs( c1 - cout ) )
  
  deallocate( c1, cout )
  
  !*****************************************************************************!
  !** Cleaning  ****************************************************************!
  !*****************************************************************************!
  call clean_sphpack()
  
  contains
  
  subroutine fill_scalar_sub(cs)
    complex(kind=dbl), intent(out) :: cs(:)
    integer                        :: j, m, jm
    real(kind=dbl)                 :: val
    
    do j = 0, jcut
      m = 0
        jm = j*(j+1)/2+m+1
        
        call random_number( val ); cs(jm)%re = val+2; cs(jm)%im = 0._dbl
      
      do m = 1, j
        jm = j*(j+1)/2+m+1
        
        call random_number( val ); cs(jm)%re = val+2
        call random_number( val ); cs(jm)%im = val+2
      end do
    end do
    
  end subroutine fill_scalar_sub
    
end program test

program test
  use sphpack
  use omp_lib
  implicit none
  
  integer, parameter :: jcut  = 497
  integer, parameter :: jmcut = jcut*(jcut+1)/2+jcut+1
  
  real(kind=dbl)                 :: start, end
  complex(kind=dbl), allocatable :: c1(:), c2(:), cout(:)
  
  type(T_grid) :: grid
  
  !*****************************************************************************!
  !** Init sphpack  ************************************************************!
  !*****************************************************************************!
  call init_sphpack( jcut )
  
  !*****************************************************************************!
  !** harmsy/harman test  ******************************************************!
  !*****************************************************************************!
  allocate( c1(jmcut)   ) ; call fill_scalar_sub( c1 )
  allocate( cout(jmcut) ) ; cout = cmplx(0._dbl, 0._dbl, kind=dbl)
  
  call alloc_grid( grid )
  
  start = omp_get_wtime()
    call harmsy( c1,      grid%tp )
    call harman( grid%tp, cout    )
  end = omp_get_wtime()
  
  call free_grid( grid )
  
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

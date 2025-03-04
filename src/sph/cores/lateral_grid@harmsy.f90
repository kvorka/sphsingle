submodule (lateral_grid) harmsy
  implicit none; contains
  
  module procedure harmsy_sub
    type(c_ptr)                         :: c_rcc
    real(kind=dbl), pointer, contiguous :: rcc(:)
    
    !Prepare output array
    call zero_rarray_sub( 2*this%lgp%n*this%fourtrans%n, grid )
    
    !Transform to suitable real input
    call alloc_aligned_sub( 32, 4*this%lgp%nrma, c_rcc, rcc )
    call this%lgp%index_bwd_sub( cin, rcc )
    
    !Transform
    call this%lgp%bwd_legesum_sub( rcc, grid )
    call this%fourtrans%fft_c2r_sub( 2*this%lgp%n, grid )
    
    !Cleaning
    call free( c_rcc )
    
  end procedure harmsy_sub
  
end submodule harmsy

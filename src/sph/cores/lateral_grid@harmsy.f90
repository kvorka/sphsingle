submodule (lateral_grid) harmsy
  implicit none; contains
  
  module procedure harmsy_sub
    type(c_ptr)                         :: c_rcc
    real(kind=dbl), pointer, contiguous :: rcc(:)
    
    !Transform to suitable real input
    call this%lgp%allocate_lgp_arr_sub( c_rcc, rcc )
    call this%lgp%index_bwd_sub( cin, rcc )
    
    call zero_rarray_sub( 2*this%lgp%n*this%fourtrans%n, grid )
    
    call this%lgp%bwd_legesum_sub( rcc, grid(1,1,1), grid(1,1,2) )
    
    call this%fourtrans%fft_c2r_sub( this%lgp%n, grid(1,1,1) )
    call this%fourtrans%fft_c2r_sub( this%lgp%n, grid(1,1,2) )
    
    !Cleaning
    call free( c_rcc )
    
  end procedure harmsy_sub
  
end submodule harmsy

submodule (lateral_grid) harman
  implicit none; contains
  
  module procedure harman_sub
    type(c_ptr)                         :: c_rcr
    real(kind=dbl), pointer, contiguous :: rcr(:)
    
    !Allocate input array
    call alloc_aligned_sub( 32, 4*this%lgp%nrma, c_rcr, rcr )
    
    call this%fourtrans%fft_r2c_sub( this%lgp%n, grid(1,1,1) )
    call this%fourtrans%fft_r2c_sub( this%lgp%n, grid(1,1,2) )
    
    call this%lgp%fwd_legesum_sub( grid(1,1,1), grid(1,1,2), rcr )
    
    !Reindex output array
    call this%lgp%index_fwd_sub( rcr, cout )
    
    !Cleaning
    call free( c_rcr )
    
  end procedure harman_sub
  
end submodule harman

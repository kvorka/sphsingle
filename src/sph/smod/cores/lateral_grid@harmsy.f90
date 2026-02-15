submodule (lateral_grid) harmsy
  implicit none; contains
  
  module procedure harmsy_sub
    type(c_ptr)                         :: c_rcc
    real(kind=dbl), pointer, contiguous :: rcc(:)
    
    !! Reindex the array into suitable form and store the input into an aligned array.
    call alloc_aligned_sub( 4*this%lgp%nrma, c_rcc, rcc )
    call this%lgp%index_bwd_sub( cin, rcc )
    
    !! Sum of associated Legendre polynomials. High grid frequencies are used to dealias FFT 
    !! and not needed anywhere else, therefore used as a work array for summation of associated
    !! Legendre polynomials.
    call this%lgp%bwd_legesum_sub( rcc, grid(1), grid(this%lgp%nFreq+1) )
    
    !! Zero the high frequencies, which were used as work arrays during the previous operations.
    call zero_rarray_sub( this%nGrid-this%lgp%nFreq, grid(this%lgp%nFreq+1) )
    
    !! Fourier transform into the physical space. FFT package vectorizes with a simd length of ndbl. 
    !! Therefore, the total number of independent FFTs is 2*this%lgp%n_dbl. FFT package leverages
    !! the parity of total number of simd transforms as well as additional factor of 2 comming from
    !! construnction of nLege.
    call this%fft%fft_c2r_sub( 2*this%lgp%n_dbl, grid )
    
    !! Cleaning
    call free( c_rcc )
    
  end procedure harmsy_sub
  
end submodule harmsy

submodule (lateral_grid) harmsy
  implicit none; contains
  
  module procedure harmsy_sub
    
    !! Sum of associated Legendre polynomials. High grid frequencies are used to dealias FFT 
    !! and not needed anywhere else, therefore used as a work array for summation of associated
    !! Legendre polynomials.
    call this%lgp%bwd_legesum_sub( cjm, grid )
    
    !! Zero the high frequencies, which were used as work arrays during the previous operations.
    call zero_rarray_sub( this%nGrid-this%lgp%nFreq, grid(this%lgp%nFreq+1) )
    
    !! Fourier transform into the physical space. FFT package vectorizes with a simd length of ndbl. 
    !! Therefore, the total number of independent FFTs is 2*this%lgp%n/ndbl. FFT package leverages
    !! the parity of total number of simd transforms as well as additional factor of 2 comming from
    !! construnction of nLege.
    call this%fft%fft_bwd_sub( this%lgp%n_dbl_2, grid )
    
  end procedure harmsy_sub
  
end submodule harmsy

submodule (lateral_grid) harman
  implicit none; contains
  
  module procedure harman_sub
    type(c_ptr)                         :: c_rcr
    real(kind=dbl), pointer, contiguous :: rcr(:)
    
    !! Fourier transform from the physical space. FFT package vectorizes with a simd length of ndbl. 
    !! Therefore, the total number of independent FFTs is 2*this%lgp%n_dbl. FFT package leverages
    !! the parity of total number of simd transforms.
    call this%fft%fft_r2c_sub( 2*this%lgp%n_dbl, grid )
    
    !! Allocate aligned temporary input array for sph coeffs and fill it with zeros.
    call alloc_aligned_sub( 4*this%lgp%nrma, c_rcr, rcr )
    
    !! Gauss-Legendre quadrature into the associated Legendre polynomials. High grid frequencies are 
    !! used to dealias fft transform and not needed afterwards, therefore used as a work array for
    !! integration.
    call this%lgp%fwd_legesum_sub( grid(1), rcr, grid(this%lgp%nFreq+1) )
    
    !! Reindex into the output array.
    call this%lgp%index_fwd_sub( rcr, cout )
    
    !! Clean.
    call free( c_rcr )
    
  end procedure harman_sub
  
end submodule harman

submodule (lateral_grid) harman
  implicit none; contains
  
  module procedure harman_sub
    real(kind=dbl), allocatable :: rcr(:)
    
    !! Fourier transform from the physical space. FFT package vectorizes with a simd length of ndbl. 
    !! Therefore, the total number of independent FFTs is 2*this%lgp%n_dbl. FFT package leverages
    !! the parity of total number of simd transforms as well as additional factor of 2 comming from
    !! construnction of nLege.
    call this%fft%fft_r2c_sub( this%lgp%n_dbl_2, grid )
    
    !! Allocate aligned temporary input array for sph coeffs.
    allocate( rcr(4*this%lgp%nrma) )
    
    !! The sum works as an accumulator, so set it to zero.
    call zero_rarray_sub( 4*this%lgp%nrma, rcr(1) )
    
    !! Gauss-Legendre quadrature into the associated Legendre polynomials. High grid frequencies are 
    !! used to dealias fft transform and not needed afterwards, therefore used as a work array for
    !! integration.
    call this%lgp%fwd_legesum_sub( grid(1), rcr )
    
    !! Reindex into the output array.
    call this%lgp%index_fwd_sub( rcr, cout )
    
    !! Clean.
    deallocate( rcr )
    
  end procedure harman_sub
  
end submodule harman

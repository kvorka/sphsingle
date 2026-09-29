submodule (lateral_grid) harman
  implicit none; contains
  
  module procedure harman_sub
    real(kind=dbl), allocatable :: cr(:)
    
    !! Allocate output array for sph coeffs.
    allocate( cr(4*this%lgp%nrma) )
    
    !! Fourier transform from the physical space. FFT package vectorizes with a simd length of ndbl. 
    !! Therefore, the total number of independent FFTs is 2*this%lgp%n/ndbl. FFT package leverages
    !! the parity of total number of simd transforms as well as additional factor of 2 comming from
    !! construnction of nLege.
    call this%fft%fft_r2c_sub( this%lgp%n_dbl_2, grid )
    
    !! Gauss-Legendre quadrature into the associated Legendre polynomials. High grid frequencies are 
    !! used to dealias fft transform and not needed afterwards, therefore used as a work array for
    !! integration.
    call this%lgp%fwd_legesum_sub( grid, cr )
    
    !! Reindex into the output array.
    call this%lgp%index_fwd_sub( cr, cout )
    
    !! Clean.
    deallocate( cr )
    
  end procedure harman_sub
  
end submodule harman

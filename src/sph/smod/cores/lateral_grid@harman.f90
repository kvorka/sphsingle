submodule (lateral_grid) harman
  implicit none; contains
  
  module procedure harman_sub
    integer                             :: i1
    type(c_ptr)                         :: c_rcr
    real(kind=dbl), pointer, contiguous :: rcr(:)
    
    !! Allocate aligned temporary input array and fill it with zeros
    call alloc_aligned_sub( 4*this%lgp%nrma, c_rcr, rcr )
    
    !$omp simd aligned (rcr:alig)
    do i1 = 1, 4*this%lgp%nrma
      rcr(i1) = 0._dbl
    end do
    
    !! Fourier transform from the physical space. FFT package leverages
    !! the fact that lgp%n is multiple of 16 for vectorization.
    call this%fft%fft_r2c_sub( this%lgp%n/16, grid )
    
    !! Gauss-Legendre quadrature into the associated Legendre polynomials.
    call this%lgp%fwd_legesum_sub( grid, rcr )
    
    !! Reindex into the output array.
    call this%lgp%index_fwd_sub( rcr, cout )
    
    !! Clean.
    call free( c_rcr )
    
  end procedure harman_sub
  
end submodule harman

submodule (lateral_grid) harmsy
  implicit none; contains
  
  module procedure harmsy_sub
    integer                             :: nGrid32, nLege32, i1, i2
    type(c_ptr)                         :: c_rcc
    real(kind=dbl), pointer, contiguous :: rcc(:)
    
    !! Size of the whole grid and size of the space used by the  backword Legendre transform. High frequencies are used
    !! to dealias fft transform and will be set to zero. Definition incorporates factor of 32 to leverage that lgp%n is
    !! divisible by 16 by construction.
    nGrid32 = 2 * this%lgp%n * this%fourtrans%n    / 32
    nLege32 = 4 * this%lgp%n * ( this%lgp%jmax+1 ) / 32
    
    !! Reindex the array into suitable form and store the input into an aligned array.
    call alloc_aligned_sub( 4*this%lgp%nrma, c_rcc, rcc )
    call this%lgp%index_bwd_sub( cin, rcc )
    
    !! Backword associated Legendre polynomials sum.
    call this%lgp%bwd_legesum_sub( rcc, grid )
    
    !! Zero the high frequencies, which were used as work arrays during the previous operations. Leverage the fact that
    !! number of grid points is divisible by 32 by construction.
    do i2 = nLege32+1, nGrid32
      !$omp simd
      do i1 = 1, 32
        grid(i1,i2) = 0._dbl
      end do
    end do
    
    !! Fourier transform into the physical space.
    call this%fourtrans%fft_c2r_sub( this%lgp%n/16, grid )
    
    !! Cleaning
    call free( c_rcc )
    
  end procedure harmsy_sub
  
end submodule harmsy

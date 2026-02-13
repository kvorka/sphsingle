submodule (lateral_grid) harmsy
  implicit none; contains
  
  module procedure harmsy_sub
    integer                             :: nGrid, nLege
    type(c_ptr)                         :: c_rcc
    real(kind=dbl), pointer, contiguous :: rcc(:)
    
    !! Size of the whole grid and size of the space used by the  backword Legendre transform. High frequencies are used
    !! to dealias fft transform and will be set to zero.
    nGrid = 2 * this%lgp%n * this%fft%n
    nLege = 4 * this%lgp%n * ( this%lgp%jmax+1 )
    
    !! Reindex the array into suitable form and store the input into an aligned array.
    call alloc_aligned_sub( 4*this%lgp%nrma, c_rcc, rcc )
    call this%lgp%index_bwd_sub( cin, rcc )
    
    !! Sum of associated Legendre polynomials.
    call this%lgp%bwd_legesum_sub( rcc, grid )
    
    !! Zero the high frequencies, which were used as work arrays during the previous operations.
    call zero_rarray_sub( nGrid-nLege, grid(nLege+1) )
    
    !! Fourier transform into the physical space. FFT package leverages the fact that lgp%n is multiple of 2*n_dbl. 
    !! Therefore, the total number of independent FFTs is 2*n_dbl*this%lgp%n_dbl. The default simd length is set to 
    !! step = 2*n_dbl within the package.
    call this%fft%fft_c2r_sub( this%lgp%n_dbl, grid )
    
    !! Cleaning
    call free( c_rcc )
    
  end procedure harmsy_sub
  
end submodule harmsy

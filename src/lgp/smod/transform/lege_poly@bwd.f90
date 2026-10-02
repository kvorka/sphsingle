submodule (lege_poly) bwd
  implicit none; contains
  
  module procedure bwd_legesum_sub
    integer                             :: im, ima1, ima2, nma
    real(kind=dbl), allocatable         :: cc(:)
    real(kind=dbl), pointer, contiguous :: swork(:), pmm(:), pmj(:), pmj1(:)
    
    !! Reindexing of real/imaginary parts of the original sequence from order-fast jm 
    !! to degree-fast mj indexing. Afterwards, reindexing and rescaling for transform.
    !! Quadruplets of the rescaled coefficietns order real/imaginary, odd/even degrees
    !! components, respectively, for cache friendly behaviour.
    allocate( cc(0:4*this%nrma-1) ); call bwd_rxd_sub( this%jmax, cjm, this%emj, this%amj, cc )
    
    !! Memory preparation, the grid frequencies above nFreq required only for 
    !! the FFT are now used as temporal storage for polynomials and partial sums
    pmm   => grid( this%nFreq + 0 * this%n : this%nFreq + 1 * this%n - 1 )
    pmj   => grid( this%nFreq + 1 * this%n : this%nFreq + 2 * this%n - 1 )
    pmj1  => grid( this%nFreq + 2 * this%n : this%nFreq + 3 * this%n - 1 )
    swork => grid( this%nFreq + 3 * this%n : this%nFreq + 7 * this%n - 1 )
    
    !! Initialization of this weird iterator
    ima2 = 0
    
    !! Cycle over the harmonic orders
    do im = 0, this%jmax-1
      ima1 = ima2
      ima2 = this%mamj(im+1)
      nma  = ima2-ima1-1
      
      call bwd_sum_m_sub( this%n_dbl_2,     &
                        & im,               &
                        & nma,              &
                        & this%fmj(2*ima1), &
                        & this%cosx,        &
                        & this%sinx,        &
                        & this%cosx2,       &
                        & pmm,              &
                        & pmj1,             &
                        & pmj,              &
                        & cc(4*(ima1)),     &
                        & swork,            &
                        & grid(4*this%n*im) )
    end do
    
    !! im == this%jmax
    call bwd_sum_jmax_sub( this%n_dbl_2,            &
                         & this%fmj(1+2*ima2),      &
                         & this%cosx,               &
                         & this%sinx,               &
                         & pmm,                     &
                         & pmj1,                    &
                         & pmj,                     &
                         & cc(4*ima2),              &
                         & swork,                   &
                         & grid(4*this%n*this%jmax) )
    
    !! Cleaning
    deallocate( cc )
    
  end procedure bwd_legesum_sub
  
end submodule bwd
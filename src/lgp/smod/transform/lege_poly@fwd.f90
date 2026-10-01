submodule (lege_poly) fwd
  implicit none; contains

  module procedure fwd_legesum_sub
    integer                             :: im, ima1, ima2
    real(kind=dbl), pointer, contiguous :: swork(:), pmm(:), pmj(:), pmj1(:)
    
    !! Memory preparation, the grid frequencies above nFreq required only for 
    !! the FFT are now used as temporal storage for polynomials and partial sums
    pmm   => grid( this%nFreq + 0 * this%n : this%nFreq + 1 * this%n - 1 )
    pmj   => grid( this%nFreq + 1 * this%n : this%nFreq + 2 * this%n - 1 )
    pmj1  => grid( this%nFreq + 2 * this%n : this%nFreq + 3 * this%n - 1 )
    swork => grid( this%nFreq + 3 * this%n : this%nFreq + 7 * this%n - 1 )
    
    !! Initialization of this weird iterator
    ima2 = 1
    
    !! Cycle over the harmonic orders
    do im = 0, this%jmax-1
      ima1 = ima2
      ima2 = this%mamj(im+1)
      
      call fwd_sum_m_sub( this%n_dbl_2,     &
                        & ima1,             &
                        & ima2-1,           &
                        & this%fmj(1,ima1), &
                        & this%cosx,        &
                        & this%sinx,        &
                        & this%cosx2,       &
                        & this%wght,        &
                        & pmm,              &
                        & pmj1,             &
                        & pmj,              &
                        & swork,            &
                        & cr(1,ima1),       &
                        & grid(4*this%n*im) )
      
    end do
    
    !! im == this%jmax
    call fwd_sum_jmax_sub( this%n_dbl_2,            &
                         & this%fmj(2,ima2),        &
                         & this%cosx,               &
                         & this%sinx,               &
                         & this%wght,               &
                         & pmm,                     &
                         & pmj1,                    &
                         & pmj,                     &
                         & swork,                   &
                         & cr(1,ima2),              &
                         & grid(4*this%n*this%jmax) )
    
  end procedure fwd_legesum_sub
  
end submodule fwd
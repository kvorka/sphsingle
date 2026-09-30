submodule (lege_poly) bwd
  implicit none; contains
  
  module procedure bwd_legesum_sub
    integer                             :: im, ima1, ima2
    real(kind=dbl), pointer, contiguous :: swork(:), pmm(:), pmj(:), pmj1(:)
    
    !! Memory preparation, the grid frequencies above nFreq required only for 
    !! the FFT are now used as temporal storage for polynomials and partial sums
    pmm   => grid( this%nFreq + 0 * this%n : this%nFreq + 1 * this%n - 1 )
    pmj   => grid( this%nFreq + 1 * this%n : this%nFreq + 2 * this%n - 1 )
    pmj1  => grid( this%nFreq + 2 * this%n : this%nFreq + 3 * this%n - 1 )
    swork => grid( this%nFreq + 3 * this%n : this%nFreq + 7 * this%n - 1 )
    
    !! Cycle over the harmonic orders
    do im = 0, this%jmax
      ima1 = this%mamj(im)
      ima2 = this%mamj(im+1)-1
      
      call bwd_sum_m_sub( this%n_dbl_2,     &
                        & ima1,             &
                        & ima2,             &
                        & this%fmj(1,ima1), &
                        & this%cosx,        &
                        & this%sinx,        &
                        & this%cosx2,       &
                        & pmm,              &
                        & pmj1,             &
                        & pmj,              &
                        & cc(1,ima1),       &
                        & swork,            &
                        & grid(4*this%n*im) )
    end do
    
  end procedure bwd_legesum_sub
  
end submodule bwd
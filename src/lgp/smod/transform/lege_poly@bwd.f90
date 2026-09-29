submodule (lege_poly) bwd
  implicit none; contains
  
  module procedure bwd_legesum_sub
    integer                             :: im, ima1, ima2
    real(kind=dbl), pointer, contiguous :: swork(:), pmj(:), pmj1(:)
    
    !! Memory preparation, the grid frequencies required only for the FFT are now used
    !! as temporal storage for polynomials and partial sums
    pmj   => grid( ( this%jmax + 1 ) * this%n * 4 : ( this%jmax + 1 ) * this%n *  5 - 1 )
    pmj1  => grid( ( this%jmax + 1 ) * this%n * 5 : ( this%jmax + 1 ) * this%n *  6 - 1 )
    swork => grid( ( this%jmax + 1 ) * this%n * 6 : ( this%jmax + 1 ) * this%n * 10 - 1 )
    
    !! Cycle over the harmonic orders
    do im = 0, this%jmax
      ima1 = this%mamj(im)
      ima2 = this%mamj(im+1)-1
      
      call bwd_sum_m_sub( this%n_dbl_2,     &
                        & ima1,             &
                        & ima2,             &
                        & this%fmj(1,ima1), &
                        & this%cosx,        &
                        & this%cosx2,       &
                        & this%pmm(:,im+1), &
                        & pmj1,             &
                        & pmj,              &
                        & cc(1,ima1),       &
                        & swork,            &
                        & grid(4*this%n*im) )
    end do
    
  end procedure bwd_legesum_sub
  
end submodule bwd
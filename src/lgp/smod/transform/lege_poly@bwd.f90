submodule (lege_poly) bwd
  implicit none; contains
  
  module procedure bwd_legesum_sub
    integer                             :: im, ima1, ima2
    type(c_ptr)                         :: c_work
    real(kind=dbl), pointer, contiguous :: work(:), swork(:), pmj(:), pmj1(:)
    
    call alloc_aligned_sub( 6*this%n, c_work, work )
    
    pmj   => work(          1 :   this%n )
    pmj1  => work(   this%n+1 : 2*this%n )
    swork => work( 2*this%n+1 : 6*this%n )
    
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
                        & grid(1,im)        )
    end do
    
    call free(c_work)
    
  end procedure bwd_legesum_sub
  
end submodule bwd
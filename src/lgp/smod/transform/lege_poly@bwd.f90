submodule (lege_poly) bwd
  implicit none; contains
  
  module procedure bwd_legesum_sub
    integer                             :: im, ima1
    type(c_ptr)                         :: c_work
    real(kind=dbl), pointer, contiguous :: work(:), acc(:), pmj(:,:), pmj1(:,:)
    
    !$omp parallel private (ima1,c_work,work,acc,pmj1,pmj)
    call alloc_aligned_sub( 6*this%n, c_work, work )
    
    pmj(1:ndbl,1:this%n_dbl)  => work(          1 :   this%n )
    pmj1(1:ndbl,1:this%n_dbl) => work(   this%n+1 : 2*this%n )
    acc                       => work( 2*this%n+1 : 6*this%n )
    
    !$omp do schedule (dynamic)
    do im = 0, this%jmax
      ima1 = this%mamj(im)
      
      call bwd_sum_m_sub( this%n_dbl,         &
                        & ima1,               &
                        & this%mamj(im+1)-1,  &
                        & this%fmj(1,ima1),   &
                        & this%cosx,          &
                        & this%cosx2,         &
                        & this%pmm(:,:,im+1), &
                        & pmj1,               &
                        & pmj,                &
                        & cc(1,ima1),         &
                        & acc,                &
                        & grid(1,im)          )
    end do
    !$omp end do
    
    call free(c_work)
    !$omp end parallel
    
  end procedure bwd_legesum_sub
  
end submodule bwd
submodule (lege_poly) fwd
  implicit none; contains

  module procedure fwd_legesum_sub
    integer                             :: im, ima1
    type(c_ptr)                         :: c_work
    real(kind=dbl), pointer, contiguous :: work(:), pmj1(:,:), pmj(:,:), acc(:), acc1(:), acc2(:)
    
    !$omp parallel private (ima1,c_work,work,acc,acc1,acc2,pmj1,pmj)
    call alloc_aligned_sub( 6*this%n+8*ndbl, c_work, work )
    
    pmj(1:ndbl,1:this%n_dbl)  => work(          1 :   this%n )
    pmj1(1:ndbl,1:this%n_dbl) => work(   this%n+1 : 2*this%n )
    acc                       => work( 2*this%n+1 : 6*this%n )
    
    acc1 => work( 6*this%n+       1 : 6*this%n + 4*ndbl )
    acc2 => work( 6*this%n+4*ndbl+1 : 6*this%n + 8*ndbl )
    
    !$omp do schedule (dynamic)
    do im = 0, this%jmax
      ima1 = this%mamj(im)
      
      call fwd_sum_m_sub( this%n_dbl,         &
                        & ima1,               &
                        & this%mamj(im+1)-1,  &
                        & this%fmj(1,ima1),   &
                        & this%cosx,          &
                        & this%cosx2,         &
                        & this%wght,          &
                        & this%pmm(:,:,im+1), &
                        & pmj1,               &
                        & pmj,                &
                        & acc,                &
                        & cr(1,ima1),         &
                        & acc1,               &
                        & acc2,               &               
                        & grid(1,im)          )

    end do
    !$omp end do
    
    call free(c_work)
    !$omp end parallel
    
  end procedure fwd_legesum_sub
  
end submodule fwd
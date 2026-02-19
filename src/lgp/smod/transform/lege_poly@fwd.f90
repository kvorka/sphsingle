submodule (lege_poly) fwd
  implicit none; contains

  module procedure fwd_legesum_sub
    integer                             :: im, ima, ima1
    real(kind=dbl), pointer, contiguous :: pmm(:), pmj1(:), pmj(:), pmj2(:), acc(:), acc2(:)
    
    pmm   => work(                     1 :   this%n         )
    pmj   => work(   this%n +          1 : 2*this%n         )
    pmj1  => work( 2*this%n +          1 : 3*this%n         )
    acc   => work( 3*this%n +          1 : 3*this%n + 4*ndbl)
    acc2  => work( 3*this%n + 4*ndbl + 1 : 3*this%n + 8*ndbl)
    
    do im = 0, this%jmax
      ima1 = this%mamj(im)
      
      call fwd_shuffle_sub( this%n_dbl,   &
                          & this%cosx,    &
                          & this%wght,    &
                          & grid(1,1,im), &
                          & grid(1,2,im), &
                          & grid(1,3,im), &
                          & grid(1,4,im)  )
      
      !! ima = ima1
        call fwd_sum1_sub( this%n_dbl,       &
                         & ima1,             &
                         & this%fmj(2,ima1), &
                         & this%sinx,        &
                         & this%cosx,        &
                         & pmm,              &
                         & pmj1,             &
                         & pmj,              &
                         & grid(1,1,im),     &
                         & grid(1,2,im),     &
                         & grid(1,3,im),     &
                         & grid(1,4,im),     &
                         & cr(1,ima1),       &
                         & acc,              &
                         & acc2              )
      
      do ima = ima1+1, this%mamj(im+1)-1
        pmj2 => pmj1
        pmj1 => pmj
        pmj  => pmj2
        
        call fwd_sum2_sub( this%n_dbl,      &
                         & this%fmj(1,ima), &
                         & this%cosx2,      &
                         & pmj1,            &
                         & pmj,             &
                         & grid(1,1,im),    &
                         & grid(1,2,im),    &
                         & grid(1,3,im),    &
                         & grid(1,4,im),    &
                         & cr(1,ima),       &
                         & acc,             &
                         & acc2             )
      end do
    end do
    
  end procedure fwd_legesum_sub
  
end submodule fwd
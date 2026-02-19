submodule (lege_poly) bwd
  implicit none; contains
  
  module procedure bwd_legesum_sub
    integer                             :: im, ima, ima1
    real(kind=dbl), pointer, contiguous :: pmm(:), pmj(:), pmj1(:), pmj2(:)
    
    pmm  => work(          1 :   this%n )
    pmj  => work(   this%n+1 : 2*this%n )
    pmj1 => work( 2*this%n+1 : 3*this%n )
    
    do im = 0, this%jmax
      ima1 = this%mamj(im)
      
      !! ima = ima1
        call bwd_sum1_sub( this%n_dbl,       &
                         & ima1,             &
                         & this%fmj(2,ima1), &
                         & this%sinx,        &
                         & this%cosx,        &
                         & pmm,              & 
                         & pmj1,             &
                         & pmj,              &
                         & cc(1,ima1),       &
                         & grid(1,1,im),     &
                         & grid(1,2,im),     &
                         & grid(1,3,im),     &
                         & grid(1,4,im)      )
      
      do ima = ima1+1, this%mamj(im+1)-1
        pmj2 => pmj1
        pmj1 => pmj
        pmj  => pmj2
        
        call bwd_sum2_sub( this%n_dbl,      &
                         & this%fmj(1,ima), &
                         & this%cosx2,      &
                         & pmj1,            &
                         & pmj,             &
                         & cc(1,ima),       &
                         & grid(1,1,im),    &
                         & grid(1,2,im),    &
                         & grid(1,3,im),    &
                         & grid(1,4,im)     )
      end do
      
      call bwd_shuffle_sub( this%n_dbl,   &
                          & this%cosx,    &
                          & grid(1,1,im), &
                          & grid(1,2,im), &
                          & grid(1,3,im), &
                          & grid(1,4,im)  )
    end do
    
  end procedure bwd_legesum_sub
  
end submodule bwd
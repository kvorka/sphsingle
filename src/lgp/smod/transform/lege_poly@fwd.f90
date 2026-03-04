submodule (lege_poly) fwd
  implicit none; contains

  module procedure fwd_legesum_sub
    integer                             :: im, ima1
    real(kind=dbl), pointer, contiguous :: pmj1(:,:), pmj(:,:)
    
    pmj(1:ndbl,1:this%n_dbl)  => work(        1 :   this%n )
    pmj1(1:ndbl,1:this%n_dbl) => work( this%n+1 : 2*this%n )
    
    do im = 0, this%jmax
      ima1 = this%mamj(im)
      
      call fwd_sum_m_sub( this%n_dbl,        &
                        & ima1,              &
                        & this%mamj(im+1)-1, &
                        & this%fmj(1,ima1),  &
                        & this%sinx,         &
                        & this%cosx,         &
                        & this%cosx2,        &
                        & this%wght,         &
                        & work(2*this%n+1),  &
                        & pmj1,              &
                        & pmj,               &
                        & work(3*this%n+1),  &
                        & cr(1,ima1),        &
                        & work(7*this%n+1),  &
                        & work(8*this%n+1),  &
                        & grid(1,im)         )
    end do
    
  end procedure fwd_legesum_sub
  
end submodule fwd
submodule (lege_poly) fwd
  implicit none; contains

  module procedure fwd_legesum_sub
    integer                             :: im, ima, ima1
    real(kind=dbl), pointer, contiguous :: swork(:), pmj1(:), pmj(:), pmj2(:), acc(:), acc2(:)
    
    pmj1  => work(                   1 :   this%n         )
    pmj   => work(   this%n         +1 : 2*this%n         )
    swork => work( 2*this%n         +1 : 6*this%n         )
    acc   => work( 6*this%n         +1 : 6*this%n + 4*ndbl)
    acc2  => work( 6*this%n + 4*ndbl+1 : 6*this%n + 8*ndbl)
    
    do im = 0, this%jmax
      ima1 = this%mamj(im)
      
      call this%fwd_shuffle_sub( grid(1,im), swork )
      
      !ima = ima1
        call this%mm_set_sub( im, pmj1, pmj )
        call this%fwd_sum_sub( pmj, swork, cr(1,ima1), acc, acc2 )
      
      do ima = ima1+1, this%mamj(im+1)-1
        pmj2 => pmj1
        pmj1 => pmj
        pmj  => pmj2
        
        call this%mj_rec_sub( ima, pmj1, pmj )
        call this%fwd_sum_sub( pmj, swork, cr(1,ima), acc, acc2 )
      end do
    end do
    
  end procedure fwd_legesum_sub
  
end submodule fwd
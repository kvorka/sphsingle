submodule (lege_poly) fwd
  implicit none; contains

  module procedure fwd_legesum_sub
    integer                             :: im, ima, ima1
    real(kind=dbl), pointer, contiguous :: swork(:), pmm(:), pmj1(:), pmj(:), pmj2(:), acc(:), acc2(:)
    
    pmm   => work(                   1 :   this%n         )
    pmj   => work(   this%n         +1 : 2*this%n         )
    pmj1  => work( 2*this%n         +1 : 3*this%n         )
    swork => work( 3*this%n         +1 : 7*this%n         )
    acc   => work( 7*this%n         +1 : 7*this%n + 4*ndbl)
    acc2  => work( 7*this%n + 4*ndbl+1 : 7*this%n + 8*ndbl)
    
    do im = 0, this%jmax
      ima1 = this%mamj(im)
      
      call this%fwd_shuffle_sub( grid(1,im), swork )
      
      !ima = ima1
        call this%mm_rec_sub( ima1, pmj1, pmj, pmm )
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
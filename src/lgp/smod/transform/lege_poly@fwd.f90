submodule (lege_poly) fwd
  implicit none; contains

  module procedure fwd_legesum_sub
    integer                             :: im, ima, ima1, ima2
    real(kind=dbl), pointer, contiguous :: swork(:), pmj1(:), pmj(:), pmj2(:), acc(:)
    
    pmj1  => grid( 1:  this%n, this%jmax+2 )
    pmj   => grid( 1:  this%n, this%jmax+3 )
    swork => grid( 1:4*this%n, this%jmax+4 )
    acc   => grid( 1:16,       this%jmax+5 )
    
    ima1 = 0
    ima2 = this%mamj(0)-1
    
    do im = 0, this%jmax
      ima1 = ima2+1
      ima2 = this%mamj(im+1)-1
      
      call this%fwd_shuffle_sub( grid(1,im), swork )
      
      !ima = ima1
        call this%mm_set_sub( im, pmj1, pmj )
        call this%fwd_sum_sub( pmj, swork, cr(1,ima1), acc )
      
      do ima = ima1+1, ima2
        pmj2 => pmj1
        pmj1 => pmj
        pmj  => pmj2
        
        call this%mj_rec_sub( ima, pmj1, pmj )
        call this%fwd_sum_sub( pmj, swork, cr(1,ima), acc )
      end do
    end do
    
  end procedure fwd_legesum_sub
  
end submodule fwd
submodule (lege_poly) bwd
  implicit none; contains
  
  module procedure bwd_legesum_sub
    integer                             :: im, ima, ima1, ima2
    real(kind=dbl), pointer, contiguous :: swork(:), pmj1(:), pmj(:), pmj2(:)
    
    pmj1  => work(          1 :   this%n )
    pmj   => work(   this%n+1 : 2*this%n )
    swork => work( 2*this%n+1 : 6*this%n )
    
    do im = 0, this%jmax
      ima1 = this%mamj(im)
      ima2 = this%mamj(im+1)-1
      
      call zero_rarray_sub( 4*this%n, swork )
      
      !ima = ima1
        call this%mm_set_sub( im, pmj1, pmj )
        call this%bwd_sum_sub( pmj, cc(1,ima1), swork )
      
      do ima = ima1+1, ima2
        pmj2 => pmj1
        pmj1 => pmj
        pmj  => pmj2
        
        call this%mj_rec_sub( ima, pmj1, pmj )
        call this%bwd_sum_sub( pmj, cc(1,ima), swork )
      end do
      
      call this%bwd_shuffle_sub( swork, grid(1,im) )
    end do
    
  end procedure bwd_legesum_sub
  
end submodule bwd
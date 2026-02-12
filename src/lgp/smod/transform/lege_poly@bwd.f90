submodule (lege_poly) bwd
  implicit none; contains
  
  module procedure bwd_legesum_sub
    integer                             :: i1, i2, im, ima, ima1, ima2
    real(kind=dbl), pointer, contiguous :: swork(:), pmj1(:), pmj(:), pmj2(:)
    
    pmj1  => grid( 1:  this%n, this%jmax+2 )
    pmj   => grid( 1:  this%n, this%jmax+3 )
    swork => grid( 1:4*this%n, this%jmax+4 )
    
    ima1 = 0
    ima2 = this%mamj(0)-1
    
    do im = 0, this%jmax
      ima1 = ima2+1
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
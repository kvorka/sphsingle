submodule (lege_poly) bwd
  implicit none; contains
  
  module procedure bwd_legesum_sub
    integer                             :: im, ima, ima1
    real(kind=dbl), pointer, contiguous :: swork(:), pmm(:), pmj(:), pmj1(:), pmj2(:)
    
    pmm   => work(          1 :   this%n )
    pmj   => work(   this%n+1 : 2*this%n )
    pmj1  => work( 2*this%n+1 : 3*this%n )
    swork => work( 3*this%n+1 : 7*this%n )
    
    do im = 0, this%jmax
      ima1 = this%mamj(im)
      
      !! ima = ima1
        call this%bwd_sum1_sub( ima1, pmm, pmj1, pmj, cc(1,ima1), swork )
      
      do ima = ima1+1, this%mamj(im+1)-1
        pmj2 => pmj1
        pmj1 => pmj
        pmj  => pmj2
        
        call this%bwd_sum2_sub( ima, pmj1, pmj, cc(1,ima), swork )
      end do
      
      call this%bwd_shuffle_sub( swork, grid(1,im) )
    end do
    
  end procedure bwd_legesum_sub
  
end submodule bwd
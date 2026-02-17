submodule (lege_poly) fwd
  implicit none; contains

  module procedure fwd_legesum_sub
    integer                             :: im, ima, ima1
    real(kind=dbl), pointer, contiguous :: swork(:), pmm(:), pmj1(:), pmj(:), pmj2(:), acc(:)
    
    pmm   => work(            1 :   this%n         )
    pmj   => work(   this%n + 1 : 2*this%n         )
    pmj1  => work( 2*this%n + 1 : 3*this%n         )
    swork => work( 3*this%n + 1 : 7*this%n         )
    acc   => work( 7*this%n + 1 : 7*this%n + 4*ndbl)
    
    do im = 0, this%jmax
      ima1 = this%mamj(im)
      
      call fwd_shuffle_sub( this%n_dbl, this%cosx, this%wght, grid(1,im), swork )
      
      !! ima = ima1
        call fwd_sum1_sub( this%n_dbl, ima1, this%fmj(2,ima1), this%sinx, this%cosx, pmm, pmj1, pmj, swork, cr(1,ima1), acc )
      
      do ima = ima1+1, this%mamj(im+1)-1
        pmj2 => pmj1
        pmj1 => pmj
        pmj  => pmj2
        
        call fwd_sum2_sub( this%n_dbl, this%fmj(1,ima), this%cosx2, pmj1, pmj, swork, cr(1,ima), acc )
      end do
    end do
    
  end procedure fwd_legesum_sub
  
end submodule fwd
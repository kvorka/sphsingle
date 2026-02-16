submodule (lege_poly) fwd_shuffle
  implicit none; contains
  
  module procedure fwd_shuffle_sub
    integer                             :: i1, i2
    real(kind=dbl)                      :: w1, w2, cx1, cx2, g11, g21, g31, g41, g12, g22, g32, g42
    real(kind=dbl), pointer, contiguous :: wght(:,:,:), cosx(:,:,:)
    
    cosx(1:ndbl,1:2,1:this%n_dbl_2) => this%cosx
    wght(1:ndbl,1:2,1:this%n_dbl_2) => this%wght
    
    do i2 = 1, this%n_dbl_2
      !$omp simd aligned (wght,cosx:alig)
      do i1 = 1, ndbl
        w1  = wght(i1,1,i2)
        w2  = wght(i1,2,i2)
        cx1 = cosx(i1,1,i2)
        cx2 = cosx(i1,2,i2)
        
        g11 = grid(i1,1,i2,1,1)
        g21 = grid(i1,1,i2,2,1)
        g31 = grid(i1,1,i2,1,2)
        g41 = grid(i1,1,i2,2,2)
        g12 = grid(i1,2,i2,1,1)
        g22 = grid(i1,2,i2,2,1)
        g32 = grid(i1,2,i2,1,2)
        g42 = grid(i1,2,i2,2,2)
        
        swork(i1,1,i2,1,1) = ( g11 - g21 ) * w1
        swork(i1,1,i2,1,2) = ( g11 + g21 ) * w1 * cx1
        swork(i1,1,i2,2,1) = ( g31 - g41 ) * w1
        swork(i1,1,i2,2,2) = ( g31 + g41 ) * w1 * cx1
        swork(i1,2,i2,1,1) = ( g12 - g22 ) * w2
        swork(i1,2,i2,1,2) = ( g12 + g22 ) * w2 * cx2
        swork(i1,2,i2,2,1) = ( g32 - g42 ) * w2
        swork(i1,2,i2,2,2) = ( g32 + g42 ) * w2 * cx2
      end do
    end do
    
  end procedure fwd_shuffle_sub
  
end submodule fwd_shuffle
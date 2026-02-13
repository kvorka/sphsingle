submodule (lege_poly) fwd_shuffle
  implicit none; contains
  
  module procedure fwd_shuffle_sub
    integer                             :: i1, i2
    real(kind=dbl), pointer, contiguous :: wght(:,:), cosx(:,:)
    
    cosx(1:ndbl,1:this%n_dbl) => this%cosx
    wght(1:ndbl,1:this%n_dbl) => this%wght
    
    do i2 = 1, this%n_dbl
      !$omp simd aligned (wght,cosx:alig)
      do i1 = 1, ndbl
        swork(i1,i2,1,1) = ( grid(i1,i2,1,1) - grid(i1,i2,2,1) ) * wght(i1,i2)
        swork(i1,i2,1,2) = ( grid(i1,i2,1,1) + grid(i1,i2,2,1) ) * wght(i1,i2) * cosx(i1,i2)
        swork(i1,i2,2,1) = ( grid(i1,i2,1,2) - grid(i1,i2,2,2) ) * wght(i1,i2)
        swork(i1,i2,2,2) = ( grid(i1,i2,1,2) + grid(i1,i2,2,2) ) * wght(i1,i2) * cosx(i1,i2)
      end do
    end do
    
  end procedure fwd_shuffle_sub
  
end submodule fwd_shuffle
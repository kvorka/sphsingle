submodule (lege_poly) fwd_shuffle
  implicit none; contains
  
  module procedure fwd_shuffle_sub
    integer                             :: i1, i2
    real(kind=dbl), pointer, contiguous :: w(:,:,:), cosx(:,:,:)
    
    cosx(1:4,1:2,1:this%n_8) => this%cosx
    w(1:4,1:2,1:this%n_8)    => this%wght
    
    do i2 = 1, this%n_8
      !$omp simd aligned (w,cosx:alig)
      do i1 = 1, 4
        swork(i1,1,i2,1,1) = ( grid(i1,1,i2,1,1) - grid(i1,1,i2,2,1) ) * w(i1,1,i2)
        swork(i1,1,i2,1,2) = ( grid(i1,1,i2,1,1) + grid(i1,1,i2,2,1) ) * w(i1,1,i2) * cosx(i1,1,i2)
        swork(i1,2,i2,1,1) = ( grid(i1,2,i2,1,1) - grid(i1,2,i2,2,1) ) * w(i1,2,i2)
        swork(i1,2,i2,1,2) = ( grid(i1,2,i2,1,1) + grid(i1,2,i2,2,1) ) * w(i1,2,i2) * cosx(i1,2,i2)
        swork(i1,1,i2,2,1) = ( grid(i1,1,i2,1,2) - grid(i1,1,i2,2,2) ) * w(i1,1,i2)
        swork(i1,1,i2,2,2) = ( grid(i1,1,i2,1,2) + grid(i1,1,i2,2,2) ) * w(i1,1,i2) * cosx(i1,1,i2)
        swork(i1,2,i2,2,1) = ( grid(i1,2,i2,1,2) - grid(i1,2,i2,2,2) ) * w(i1,2,i2)
        swork(i1,2,i2,2,2) = ( grid(i1,2,i2,1,2) + grid(i1,2,i2,2,2) ) * w(i1,2,i2) * cosx(i1,2,i2)
      end do
    end do
    
  end procedure fwd_shuffle_sub
  
end submodule fwd_shuffle
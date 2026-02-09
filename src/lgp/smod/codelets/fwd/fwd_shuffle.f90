submodule (lege_poly) fwd_shuffle
  implicit none; contains
  
  module procedure fwd_shuffle_sub
    integer                             :: i1, i2
    real(kind=dbl), pointer, contiguous :: w(:,:), cosx(:,:)
    
    w(1:16,1:this%n16)    => this%wght
    cosx(1:16,1:this%n16) => this%cosx
    
    do i2 = 1, this%n16
      !$omp simd
      do i1 = 1, 16
        swork(i1,i2,1,1) = ( grid(i1,i2,1,1) - grid(i1,i2,2,1) ) * w(i1,i2)
        swork(i1,i2,2,1) = ( grid(i1,i2,1,2) - grid(i1,i2,2,2) ) * w(i1,i2)
        swork(i1,i2,1,2) = ( grid(i1,i2,1,1) + grid(i1,i2,2,1) ) * w(i1,i2) * cosx(i1,i2)
        swork(i1,i2,2,2) = ( grid(i1,i2,1,2) + grid(i1,i2,2,2) ) * w(i1,i2) * cosx(i1,i2)
      end do
    end do
    
  end procedure fwd_shuffle_sub
  
end submodule fwd_shuffle
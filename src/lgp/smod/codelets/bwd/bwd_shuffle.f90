submodule (lege_poly) bwd_shuffle
  implicit none; contains
  
  module procedure bwd_shuffle_sub
    integer                             :: i1, i2
    real(kind=dbl), pointer, contiguous :: cosx(:,:,:)
    
    cosx(1:4,1:2,1:this%n_8) => this%cosx
    
    do i2 = 1, this%n_8
      !$omp simd aligned (cosx:alig)
      do i1 = 1, 4
        grid(i1,1,i2,1,1) = swork(i1,1,i2,1,2) * cosx(i1,1,i2) + swork(i1,1,i2,1,1)
        grid(i1,1,i2,2,1) = swork(i1,1,i2,1,2) * cosx(i1,1,i2) - swork(i1,1,i2,1,1)
        grid(i1,2,i2,1,1) = swork(i1,2,i2,1,2) * cosx(i1,2,i2) + swork(i1,2,i2,1,1)
        grid(i1,2,i2,2,1) = swork(i1,2,i2,1,2) * cosx(i1,2,i2) - swork(i1,2,i2,1,1)
        grid(i1,1,i2,1,2) = swork(i1,1,i2,2,2) * cosx(i1,1,i2) + swork(i1,1,i2,2,1)
        grid(i1,1,i2,2,2) = swork(i1,1,i2,2,2) * cosx(i1,1,i2) - swork(i1,1,i2,2,1)
        grid(i1,2,i2,1,2) = swork(i1,2,i2,2,2) * cosx(i1,2,i2) + swork(i1,2,i2,2,1)
        grid(i1,2,i2,2,2) = swork(i1,2,i2,2,2) * cosx(i1,2,i2) - swork(i1,2,i2,2,1)
      end do
    end do
    
  end procedure bwd_shuffle_sub

end submodule bwd_shuffle
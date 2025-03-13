submodule (lege_poly) bwd_shuffle
  implicit none; contains
  
  module procedure bwd_shuffle_c
    integer :: i1, i2
    
    do i2 = 1, (n/32)*32, 32
      !$omp simd
      do i1 = 0, 31
        grid(i1+i2,1,1) = swork(i1+i2,1,2) * cosx(i1+i2) + swork(i1+i2,1,1)
        grid(i1+i2,2,1) = swork(i1+i2,1,2) * cosx(i1+i2) - swork(i1+i2,1,1)
        grid(i1+i2,1,2) = swork(i1+i2,2,2) * cosx(i1+i2) + swork(i1+i2,2,1)
        grid(i1+i2,2,2) = swork(i1+i2,2,2) * cosx(i1+i2) - swork(i1+i2,2,1)
      end do
    end do
    
    do i2 = (n/32)*32+1, n, 8
      !$omp simd
      do i1 = 0, 7
        grid(i1+i2,1,1) = swork(i1+i2,1,2) * cosx(i1+i2) + swork(i1+i2,1,1)
        grid(i1+i2,2,1) = swork(i1+i2,1,2) * cosx(i1+i2) - swork(i1+i2,1,1)
        grid(i1+i2,1,2) = swork(i1+i2,2,2) * cosx(i1+i2) + swork(i1+i2,2,1)
        grid(i1+i2,2,2) = swork(i1+i2,2,2) * cosx(i1+i2) - swork(i1+i2,2,1)
      end do
    end do
    
  end procedure bwd_shuffle_c

end submodule bwd_shuffle
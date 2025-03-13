submodule (lege_poly) fwd_shuffle
  implicit none; contains
  
  module procedure fwd_shuffle_c
    integer :: i1, i2
    
    do i2 = 1, (n/32)*32, 32
      !$omp simd
      do i1 = 0, 31
        swork(i1+i2,1,1) = ( grid(i1+i2,1,1) - grid(i1+i2,2,1) ) * w(i1+i2)
        swork(i1+i2,2,1) = ( grid(i1+i2,1,2) - grid(i1+i2,2,2) ) * w(i1+i2)
        swork(i1+i2,1,2) = ( grid(i1+i2,1,1) + grid(i1+i2,2,1) ) * w(i1+i2) * cosx(i1+i2)
        swork(i1+i2,2,2) = ( grid(i1+i2,1,2) + grid(i1+i2,2,2) ) * w(i1+i2) * cosx(i1+i2)
      end do
    end do
    
    do i2 = (n/32)*32+1, n, 8
      !$omp simd
      do i1 = 0, 7
        swork(i1+i2,1,1) = ( grid(i1+i2,1,1) - grid(i1+i2,2,1) ) * w(i1+i2)
        swork(i1+i2,2,1) = ( grid(i1+i2,1,2) - grid(i1+i2,2,2) ) * w(i1+i2)
        swork(i1+i2,1,2) = ( grid(i1+i2,1,1) + grid(i1+i2,2,1) ) * w(i1+i2) * cosx(i1+i2)
        swork(i1+i2,2,2) = ( grid(i1+i2,1,2) + grid(i1+i2,2,2) ) * w(i1+i2) * cosx(i1+i2)
      end do
    end do
    
  end procedure fwd_shuffle_c
  
end submodule fwd_shuffle
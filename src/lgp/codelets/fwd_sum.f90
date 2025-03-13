submodule (lege_poly) fwd_sum
  implicit none; contains
  
  module procedure fwd_sum_c
    integer :: i1, i2
    
    do i1 = 1, (n/32)*32, 32
      do i2 = 0, 31
        cr(1) = cr(1) + pmj(i2+i1) * swork(i2+i1,1)
        cr(2) = cr(2) + pmj(i2+i1) * swork(i2+i1,2)
        cr(3) = cr(3) + pmj(i2+i1) * swork(i2+i1,3)
        cr(4) = cr(4) + pmj(i2+i1) * swork(i2+i1,4)
      end do
    end do
    
    do i1 = (n/32)*32+1, n, 8
      do i2 = 0, 7
        cr(1) = cr(1) + pmj(i2+i1) * swork(i2+i1,1)
        cr(2) = cr(2) + pmj(i2+i1) * swork(i2+i1,2)
        cr(3) = cr(3) + pmj(i2+i1) * swork(i2+i1,3)
        cr(4) = cr(4) + pmj(i2+i1) * swork(i2+i1,4)
      end do
    end do
    
  end procedure fwd_sum_c
  
end submodule fwd_sum

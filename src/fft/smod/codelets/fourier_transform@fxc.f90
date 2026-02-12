submodule (fourier_transform) fxc
  implicit none; contains
  
  module procedure fxc2r
    integer        :: i1, i2
    real(kind=dbl) :: t1, t2, addre, subre, addim, subim
    
    t1 = t(1)
    t2 = t(2)
    
    do i2 = 1, m
      !$omp simd
      do i1 = 1, 16
        addre = x11(i1,i2) + x21(i1,i2)
        subre = x11(i1,i2) - x21(i1,i2)
        addim = x12(i1,i2) + x22(i1,i2)
        subim = x12(i1,i2) - x22(i1,i2)
        
        x11(i1,i2) = addre - subre * t2 - addim * t1
        x12(i1,i2) = subim - addim * t2 + subre * t1
        
        x21(i1,i2) = -x11(i1,i2) + 2 * addre
        x22(i1,i2) = +x12(i1,i2) - 2 * subim
      end do
    end do
    
  end procedure fxc2r
  
end submodule fxc
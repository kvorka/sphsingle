submodule (fourier_transform) fxr
  implicit none; contains
  
  module procedure fxr2c
    integer        :: i1, i2
    real(kind=dbl) :: t1, t2, addre, subre, addim, subim
    
    t1 = t(1)
    t2 = t(2)
    
    !GCC$ unroll 4
    !DIR$ unroll (4)
    do i2 = 1, m
      !$omp simd
      do i1 = 1, ndbl
        addre = x11(i1,i2) + x21(i1,i2)
        subre = x11(i1,i2) - x21(i1,i2)
        addim = x22(i1,i2) + x12(i1,i2)
        subim = x22(i1,i2) - x12(i1,i2)
        
        x11(i1,i2) = ( addre + subre * t2 + addim * t1 ) / 2
        x12(i1,i2) = ( subim - addim * t2 + subre * t1 ) / 2
        
        x21(i1,i2) = -x11(i1,i2) + addre
        x22(i1,i2) = +x12(i1,i2) - subim
      end do
    end do
    
  end procedure fxr2c
  
end submodule fxr
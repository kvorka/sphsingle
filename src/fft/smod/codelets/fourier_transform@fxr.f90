submodule (fourier_transform) fxr
  implicit none; contains
  
  module procedure fxr2c
    integer        :: i1, i2
    real(kind=dbl) :: t1, t2, addre1, subre1, addim1, subim1, addre2, subre2, addim2, subim2
    
    t1 = t(1)
    t2 = t(2)
    
    do i2 = 1, m, 2
      !$omp simd
      do i1 = 1, ndbl
        addre1 = x11(i1,i2  ) + x21(i1,i2  )
        subre1 = x11(i1,i2  ) - x21(i1,i2  )
        addim1 = x22(i1,i2  ) + x12(i1,i2  )
        subim1 = x22(i1,i2  ) - x12(i1,i2  )
        addre2 = x11(i1,i2+1) + x21(i1,i2+1)
        subre2 = x11(i1,i2+1) - x21(i1,i2+1)
        addim2 = x22(i1,i2+1) + x12(i1,i2+1)
        subim2 = x22(i1,i2+1) - x12(i1,i2+1)
        
        x11(i1,i2  ) = ( addre1 + subre1 * t2 + addim1 * t1 ) / 2
        x12(i1,i2  ) = ( subim1 - addim1 * t2 + subre1 * t1 ) / 2
        x11(i1,i2+1) = ( addre2 + subre2 * t2 + addim2 * t1 ) / 2
        x12(i1,i2+1) = ( subim2 - addim2 * t2 + subre2 * t1 ) / 2
        
        x21(i1,i2  ) = -x11(i1,i2  ) + addre1
        x22(i1,i2  ) = +x12(i1,i2  ) - subim1
        x21(i1,i2+1) = -x11(i1,i2+1) + addre2
        x22(i1,i2+1) = +x12(i1,i2+1) - subim2
      end do
    end do
    
  end procedure fxr2c
  
end submodule fxr
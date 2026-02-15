submodule (fourier_transform) fxad
  implicit none; contains
  
  module procedure fxaddsub
    integer        :: i1, i2
    real(kind=dbl) :: add1, add2, add3, add4
    
    do i2 = 1, m, 4
      !$omp simd
      do i1 = 1, ndbl
        add1 = arr1(i1,i2  )
        add2 = arr1(i1,i2+1)
        add3 = arr1(i1,i2+2)
        add4 = arr1(i1,i2+3)
        
        arr1(i1,i2  ) = arr1(i1,i2  ) + arr2(i1,i2  )
        arr1(i1,i2+1) = arr1(i1,i2+1) + arr2(i1,i2+1)
        arr1(i1,i2+2) = arr1(i1,i2+2) + arr2(i1,i2+2)
        arr1(i1,i2+3) = arr1(i1,i2+3) + arr2(i1,i2+3)

        arr2(i1,i2  ) = add1 - arr2(i1,i2  )
        arr2(i1,i2+1) = add2 - arr2(i1,i2+1)
        arr2(i1,i2+2) = add3 - arr2(i1,i2+2)
        arr2(i1,i2+3) = add4 - arr2(i1,i2+3)
      end do
    end do
    
  end procedure fxaddsub
  
end submodule fxad
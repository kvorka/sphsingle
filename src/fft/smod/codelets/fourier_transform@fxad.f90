submodule (fourier_transform) fxad
  implicit none; contains
  
  module procedure fxaddsub
    integer        :: i1, i2
    real(kind=dbl) :: add1, add2
    
    do i2 = 1, m, 2
      !$omp simd
      do i1 = 1, ndbl
        add1 = arr1(i1,i2  )
        add2 = arr1(i1,i2+1)
        
        arr1(i1,i2  ) = arr1(i1,i2  ) + arr2(i1,i2  )
        arr1(i1,i2+1) = arr1(i1,i2+1) + arr2(i1,i2+1)

        arr2(i1,i2  ) = add1 - arr2(i1,i2  )
        arr2(i1,i2+1) = add2 - arr2(i1,i2+1)
      end do
    end do
    
  end procedure fxaddsub
  
end submodule fxad
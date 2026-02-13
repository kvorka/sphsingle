submodule (fourier_transform) fxad
  implicit none; contains
  
  module procedure fxaddsub
    integer        :: i1, i2
    real(kind=dbl) :: add
    
    do i2 = 1, m
      !$omp simd
      do i1 = 1, step
        add         =               arr1(i1,i2)
        arr1(i1,i2) = arr1(i1,i2) + arr2(i1,i2)
        arr2(i1,i2) = add         - arr2(i1,i2)
      end do
    end do
    
  end procedure fxaddsub
  
end submodule fxad
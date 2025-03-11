submodule (math) arrcopy
  implicit none; contains
  
  module procedure copy_rarray_c
    integer :: i
    
    do i = 1, n
      arrto(i) = arrfrom(i)
    end do
    
  end procedure copy_rarray_c
  
end submodule arrcopy
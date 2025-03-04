submodule (fourier_transform) fxshf
  implicit none; contains
  
  module procedure fxzshf
    integer                             :: iv, j, isj, isj2
    type(c_ptr)                         :: c_y
    real(kind=dbl), pointer, contiguous :: y(:)
    
    call alloc_aligned_sub( default_alig, 2*howmany, c_y, y )
    
    j = 1  
      do while (j <= this%n/2-2)
        isj = this%it(j)
        
        if (isj < 0) then
          j = j + 1
          
        else
          !$omp simd
          do iv = 1, 2*howmany
            y(iv) = x(iv,isj)
          end do
          
          do
            j    = j + 1
            isj2 = this%it(j)
            
            if ( isj2 < 0 ) then
              !$omp simd
              do iv = 1, 2*howmany
                x(iv,isj)      = x(iv,isj2-imm)
                x(iv,isj2-imm) = y(iv)
              end do
              
              j = j + 1
              exit
            
            else
              !$omp simd
              do iv = 1, 2*howmany
                x(iv,isj) = x(iv,isj2)
              end do
              
              isj = isj2
            end if
          end do
        end if
      end do
    
    call free( c_y )
    
  end procedure fxzshf
  
end submodule fxshf

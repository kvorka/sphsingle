submodule (fourier_transform) fxshf
  implicit none; contains
  
  module procedure fxzshf
    integer                             :: iv, iv1, j, isj, isj2
    type(c_ptr)                         :: c_y
    real(kind=dbl), pointer, contiguous :: y(:)
    
    call alloc_aligned_sub( 2*howmany, c_y, y )
    
    j = 1  
      do while (j <= this%n/2-2)
        isj = this%it(j)
        
        if (isj < 0) then
          j = j + 1
          
        else
          do iv = 1, 2*howmany, 32
            !$omp simd
            do iv1 = 0, 31
              y(iv1+iv) = x(iv1+iv,isj)
            end do
          end do
          
          do
            j    = j + 1
            isj2 = this%it(j)
            
            if ( isj2 < 0 ) then
              do iv = 1, 2*howmany, 32
                !$omp simd
                do iv1 = 0, 31
                  x(iv1+iv,isj)      = x(iv1+iv,isj2-imm)
                  x(iv1+iv,isj2-imm) = y(iv1+iv)
                end do
              end do
              
              j = j + 1
              exit
            
            else
              do iv = 1, 2*howmany, 32
                !$omp simd
                do iv1 = 0, 31
                  x(iv1+iv,isj) = x(iv1+iv,isj2)
                end do
              end do
              
              isj = isj2
            end if
          end do
        end if
      end do
    
    call free( c_y )
    
  end procedure fxzshf
  
end submodule fxshf
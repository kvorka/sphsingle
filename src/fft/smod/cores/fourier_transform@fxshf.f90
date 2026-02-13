submodule (fourier_transform) fxshf
  implicit none; contains
  
  module procedure fxzshf
    integer                             :: j, i30, i31
    real(kind=dbl), pointer, contiguous :: y(:)
    type(c_ptr)                         :: c_y
    
    call alloc_aligned_sub( m*step, c_y, y )
    
    j = 1
    
    do while ( j <= n/2-2 )
      i30 = it(j)
      
      if ( i30 < 0 ) then
        j = j + 1
        
      else
        call fxcpy( m, x(1,i30), y )
        
        do
          j   = j + 1
          i31 = it(j)
          
          if ( i31 < 0 ) then
            call fxcpy( m, x(1,i31-imm), x(1,i30    ) )
            call fxcpy( m, y,            x(1,i31-imm) )
            
            j = j + 1
            exit
          
          else
            call fxcpy( m, x(1,i31), x(1,i30) )
            
            i30 = i31
          end if
        end do
      end if
    end do
    
    call free( c_y )
    
  end procedure fxzshf
  
end submodule fxshf
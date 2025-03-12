module sphpack
  use lateral_grid
  implicit none; public
  
  type(T_lateralGrid), private   :: sph
  integer,             protected :: nth, nph
  
  contains
  
  subroutine init_sphpack(jcut)
    integer, intent(in) :: jcut
    
    call sph%init_sub( jcut )
    
    nth = 2 * sph%lgp%n
    nph = sph%fourtrans%n
    
  end subroutine init_sphpack
  
  subroutine alloc_grid(c_grid, f_grid)
    type(c_ptr),                         intent(out) :: c_grid
    real(kind=dbl), pointer, contiguous, intent(out) :: f_grid(:,:)
    
    call sph%alloc_grid_sub( c_grid, f_grid )
    
  end subroutine alloc_grid
  
  subroutine harmsy(cajm, grid)
    complex(kind=dbl), intent(in)  :: cajm(*)
    real(kind=dbl),    intent(out) :: grid(*)
    
    call sph%harmsy_sub( cajm, grid )
    
  end subroutine harmsy
  
  subroutine harman(grid, cajm)
    real(kind=dbl),    intent(inout) :: grid(*)
    complex(kind=dbl), intent(out)   :: cajm(*)
    
    call sph%harman_sub( grid, cajm )
    
  end subroutine harman
  
  subroutine clean_sphpack()
    
    call sph%deallocate_sub()
    
  end subroutine clean_sphpack
  
end module sphpack
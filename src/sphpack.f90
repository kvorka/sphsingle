module sphpack
  use lateral_grid
  implicit none; public
  
  type(T_lateralGrid), private   :: sph
  integer,             protected :: nth, nph
  
  contains
  
  subroutine init_sphpack(jcut)
    integer, intent(in) :: jcut
    
    call sph%init_sub( jcut )
    
    nth = sph%lgp%n
    nph = sph%fourtrans%n
    
  end subroutine init_sphpack
  
  subroutine harmsy(cajm, grid)
    complex(kind=dbl), intent(in)  :: cajm(*)
    real(kind=dbl),    intent(out) :: grid(*)
    
    call sph%harmsy_sub( cajm, grid )
    
  end subroutine harmsy
  
  subroutine harman(grid, cajm)
    real(kind=dbl),    intent(in)  :: grid(*)
    complex(kind=dbl), intent(out) :: cajm(*)
    
    call sph%harman_sub( grid, cajm )
    
  end subroutine harman
  
  subroutine clean_sphpack()
    
    call sph%deallocate_sub()
    
  end subroutine clean_sphpack
  
end module sphpack
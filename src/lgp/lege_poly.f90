module lege_poly
  use math
  implicit none
  
  type, public :: T_legep
    integer                             :: jmax, jms, n, nrma
    integer,        allocatable         :: mma(:)
    real(kind=dbl), allocatable         :: emj(:), fmj(:,:), amj(:)
    real(kind=dbl), pointer, contiguous :: cosx(:), cosx2(:), wght(:), pmm(:,:)
    type(c_ptr)                         :: c_cosx, c_cosx2, c_wght, c_pmm
    
    contains
    
    procedure, public,  pass :: init_sub       => init_lege_sub
    procedure, private, pass :: roots_sub      => find_roots_sub
    procedure, private, pass :: coeffs_sub     => compute_coeffs_sub
    procedure, private, pass :: pmm_sub        => compute_pmm_sub
    procedure, public,  pass :: deallocate_sub => deallocate_lege_sub
    
    procedure, public, pass :: index_bwd_sub, bwd_legesum_sub
    procedure, public, pass :: index_fwd_sub, fwd_legesum_sub
    
  end type T_legep
  
  real(kind=dbl), parameter :: deps = 1.0d-15
  real(kind=qbl), parameter :: qeps = 1.0d-28
  real(kind=qbl), parameter :: qpi  = acos(-1._qbl)
  
  interface
    module subroutine init_lege_sub(this, jmax, n, wfac)
      class(T_legep), intent(inout) :: this
      integer,        intent(in)    :: jmax, n
      real(kind=dbl), intent(in)    :: wfac
    end subroutine init_lege_sub
    
    module subroutine deallocate_lege_sub(this)
      class(T_legep), intent(inout) :: this
    end subroutine deallocate_lege_sub
    
    module subroutine find_roots_sub(this)
      class(T_legep), intent(inout) :: this
    end subroutine find_roots_sub
    
    module subroutine compute_coeffs_sub(this)
      class(T_legep), intent(inout) :: this
    end subroutine compute_coeffs_sub
    
    module subroutine compute_pmm_sub(this)
      class(T_legep), intent(inout) :: this
    end subroutine compute_pmm_sub
    
    module subroutine index_bwd_sub(this, cjm, rcab)
      class(T_legep),    intent(in)  :: this
      complex(kind=dbl), intent(in)  :: cjm(this%jms)
      real(kind=dbl),    intent(out) :: rcab(4*this%nrma)
    end subroutine index_bwd_sub
    
    module subroutine index_fwd_sub(this, rcab, cjm)
      class(T_legep),    intent(in)    :: this
      real(kind=dbl),    intent(inout) :: rcab(4*this%nrma)
      complex(kind=dbl), intent(out)   :: cjm(this%jms)
    end subroutine index_fwd_sub
    
    module subroutine bwd_legesum_sub(this, cc, grid)
      class(T_legep), intent(in)  :: this
      real(kind=dbl), intent(out) :: grid(4*this%n,0:this%jmax)
      real(kind=dbl), intent(in)  :: cc(4,this%nrma)
    end subroutine bwd_legesum_sub
    
    module subroutine fwd_legesum_sub(this, grid, cr)
      class(T_legep), intent(in)    :: this
      real(kind=dbl), intent(in)    :: grid(4*this%n,0:this%jmax)
      real(kind=dbl), intent(inout) :: cr(4,this%nrma)
    end subroutine fwd_legesum_sub
  end interface
  
  interface
    module pure subroutine bwd_indx_c(jmax, emj, icab, ocab)
      integer,        intent(in)  :: jmax
      real(kind=dbl), intent(in)  :: emj(*), icab(2,*)
      real(kind=dbl), intent(out) :: ocab(2,2,*)
    end subroutine bwd_indx_c
    
    module pure subroutine fwd_indx_c(jmax, emj, ocab, icab)
      integer,        intent(in)  :: jmax
      real(kind=dbl), intent(in)  :: emj(*), ocab(2,2,*)
      real(kind=dbl), intent(out) :: icab(2,*)
    end subroutine fwd_indx_c
    
    module pure subroutine is_rescale_c(n, cff, rcab)
      integer,        intent(in)    :: n
      real(kind=dbl), intent(in)    :: cff(n)
      real(kind=dbl), intent(inout) :: rcab(4,n)
    end subroutine is_rescale_c
    
    module pure subroutine mm_set_c(n, pmm, pmj1, pmj)
      integer,        intent(in)  :: n
      real(kind=dbl), intent(in)  :: pmm(n)
      real(kind=dbl), intent(out) :: pmj1(n), pmj(n)
    end subroutine mm_set_c
    
    module pure subroutine mj_rec_c(n, cff, cosx2, pmj1, pmj)
      integer,        intent(in)    :: n
      real(kind=dbl), intent(in)    :: cff(2), cosx2(n)
      real(kind=dbl), intent(inout) :: pmj1(n), pmj(n)
    end subroutine mj_rec_c
    
    module pure subroutine bwd_sum_c(n, pmj, cc, swork)
      integer,        intent(in)  :: n
      real(kind=dbl), intent(in)  :: pmj(n)
      real(kind=dbl), intent(in)  :: cc(4)
      real(kind=dbl), intent(out) :: swork(n,4)
    end subroutine bwd_sum_c
    
    module pure subroutine bwd_shuffle_c(n, cosx, swork, grid)
      integer,        intent(in)    :: n
      real(kind=dbl), intent(in)    :: cosx(n)
      real(kind=dbl), intent(inout) :: swork(n,2,2)
      real(kind=dbl), intent(out)   :: grid(n,2,2)
    end subroutine bwd_shuffle_c
    
    module pure subroutine fwd_shuffle_c(n, cosx, w, grid, swork)
      integer,        intent(in)  :: n
      real(kind=dbl), intent(in)  :: w(n), cosx(n)
      real(kind=dbl), intent(in)  :: grid(n,2,2)
      real(kind=dbl), intent(out) :: swork(n,2,2)
    end subroutine fwd_shuffle_c
    
    module pure subroutine fwd_sum_c(n, pmj, swork, cr)
      integer,        intent(in)    :: n
      real(kind=dbl), intent(in)    :: pmj(n)
      real(kind=dbl), intent(in)    :: swork(n,4)
      real(kind=dbl), intent(inout) :: cr(4)
    end subroutine fwd_sum_c
  end interface
  
end module lege_poly

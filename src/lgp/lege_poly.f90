module lege_poly
  use math
  implicit none
  
  type, public :: T_legep
    integer                     :: jmax, jms, n, nrma
    real(kind=dbl), allocatable :: emj(:), fmj(:,:), amj(:), rw(:,:)
    
    contains
    
    procedure, public,  pass :: init_sub       => init_lege_sub
    procedure, private, pass :: roots_sub      => find_roots_sub
    procedure, private, pass :: coeffs_sub     => compute_coeffs_sub
    procedure, public,  pass :: deallocate_sub => deallocate_lege_sub
    
    procedure, public, pass :: allocate_lgp_arr_sub
    
    procedure, private, pass :: is_rescale_sub
    procedure, public,  pass :: index_bwd_sub
    procedure, public,  pass :: index_fwd_sub
    
    procedure, public,  pass :: bwd_legesum_sub
    procedure, public,  pass :: fwd_legesum_sub
    
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
    
    module subroutine allocate_lgp_arr_sub(this, arr)
      class(T_legep),              intent(in)  :: this
      real(kind=dbl), allocatable, intent(out) :: arr(:)
    end subroutine allocate_lgp_arr_sub
    
    module subroutine is_rescale_sub(this, rcab)
      class(T_legep), intent(in)    :: this
      real(kind=dbl), intent(inout) :: rcab(4,this%nrma)
    end subroutine is_rescale_sub
    
    module subroutine index_bwd_sub(this, cjm, rcab)
      class(T_legep),    intent(in)  :: this
      complex(kind=dbl), intent(in)  :: cjm(*)
      real(kind=dbl),    intent(out) :: rcab(2,2,this%nrma)
    end subroutine index_bwd_sub
    
    module  subroutine index_fwd_sub(this, rcab, cjm)
      class(T_legep),    intent(in)    :: this
      real(kind=dbl),    intent(inout) :: rcab(2,2,this%nrma)
      complex(kind=dbl), intent(out)   :: cjm(*)
    end subroutine index_fwd_sub
    
    module  subroutine bwd_legesum_sub(this, cc, sumN, sumS, cosx, sinx, cosx2, pmm, pmj, pmj1, swork)
      class(T_legep), intent(in)  :: this
      real(kind=dbl), intent(in)  :: cosx(step), sinx(step), cosx2(step)
      real(kind=dbl), intent(out) :: swork(4*step), pmm(step), pmj(step), pmj1(step)
      real(kind=dbl), intent(out) :: sumN(2*step,0:this%jmax), sumS(2*step,0:this%jmax)
      real(kind=dbl), intent(in)  :: cc(4,this%nrma)
    end subroutine bwd_legesum_sub
    
    module  subroutine fwd_legesum_sub(this, sumN, sumS, cr, cosx, sinx, cosx2, weight, pmm, pmj, pmj1, swork)
      class(T_legep), intent(in)    :: this
      real(kind=dbl), intent(in)    :: cosx(step), sinx(step), cosx2(step), weight(step)
      real(kind=dbl), intent(out)   :: swork(4*step), pmm(step), pmj(step), pmj1(step)
      real(kind=dbl), intent(in)    :: sumN(2*step,0:this%jmax), sumS(2*step,0:this%jmax)
      real(kind=dbl), intent(inout) :: cr(4,this%nrma)
    end subroutine fwd_legesum_sub
  end interface
  
  interface
    module pure subroutine mmset_sub(ma, cff, cosx, sinx, pmm, pmj1, pmj)
      integer,        intent(in)    :: ma
      real(kind=dbl), intent(in)    :: cff, cosx(step), sinx(step)
      real(kind=dbl), intent(inout) :: pmm(step)
      real(kind=dbl), intent(out)   :: pmj1(step), pmj(step)
    end subroutine mmset_sub
    
    module pure subroutine mjrec_sub(cff, cosx2, pmj1, pmj)
      real(kind=dbl), intent(in)    :: cff(2), cosx2(step)
      real(kind=dbl), intent(inout) :: pmj1(step), pmj(step)
    end subroutine mjrec_sub
    
    module pure subroutine bwd_sum_sub(pmj, cc, swork)
      real(kind=dbl), intent(in)  :: pmj(step)
      real(kind=dbl), intent(in)  :: cc(4)
      real(kind=dbl), intent(out) :: swork(step,4)
    end subroutine bwd_sum_sub
    
    module pure subroutine bwd_shuffle_sub(cosx, swork, sumN, sumS)
      real(kind=dbl), intent(in)    :: cosx(step)
      real(kind=dbl), intent(inout) :: swork(step,2,2)
      real(kind=dbl), intent(out)   :: sumN(step,2), sumS(step,2)
    end subroutine bwd_shuffle_sub
    
    module pure subroutine fwd_sum_sub(pmj, swork, cr)
      real(kind=dbl), intent(in)    :: pmj(step)
      real(kind=dbl), intent(in)    :: swork(step,4)
      real(kind=dbl), intent(inout) :: cr(4)
    end subroutine fwd_sum_sub
    
    module pure subroutine fwd_shuffle_sub(w, cosx, sumN, sumS, swork)
      real(kind=dbl), intent(in)  :: w(step), cosx(step)
      real(kind=dbl), intent(in)  :: sumN(step,2), sumS(step,2)
      real(kind=dbl), intent(out) :: swork(step,2,2)
    end subroutine fwd_shuffle_sub
  end interface
  
end module lege_poly

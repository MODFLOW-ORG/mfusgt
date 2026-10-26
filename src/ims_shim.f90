!> @file ims_shim.f90
!! Phase-1 isolation shims for the IMS linear solver extracted from
!! MODFLOW 6 tag 6.4.4. These modules replace the MF6 core modules that the
!! four IMS files (ims8base/linear/misc/reordering) `use`, so the solver can
!! be compiled and tested *in isolation* under the production gfortran flags
!! (-fdec -fallow-argument-mismatch -std=legacy).
!!
!! Each shim is trimmed to exactly the symbols the IMS files reference (see the
!! Phase-1 report). NONE of these shims pull additional MF6 modules -- they are
!! flat and self-contained. In the USG-T integration phase the SimModule shim is
!! where IMS failure calls route to USG-T's USTOP.

! ---------------------------------------------------------------------------
!  KindModule -- kept (tiny), matches MF6 KindModule for DP / I4B / LGP
! ---------------------------------------------------------------------------
module KindModule
  implicit none
  public
  integer, parameter :: DP = KIND(1.0D0)            !< double precision real kind
  integer, parameter :: I4B = SELECTED_INT_KIND(8)  !< 4-byte integer kind
  integer, parameter :: LGP = KIND(.true.)          !< default logical kind
end module KindModule

! ---------------------------------------------------------------------------
!  ConstantsModule -- trimmed: only the constants the 4 IMS files import.
!  Values are verbatim from MF6 6.4.4 src/Utilities/Constants.f90.
! ---------------------------------------------------------------------------
module ConstantsModule
  use KindModule, only: I4B, DP
  implicit none
  public
  ! -- string lengths
  integer(I4B), parameter :: LINELENGTH = 300
  integer(I4B), parameter :: LENSOLUTIONNAME = 16
  integer(I4B), parameter :: LENMEMPATH = 256
  ! -- integer constants
  integer(I4B), parameter :: IZERO = 0
  ! -- real constants
  real(DP), parameter :: DZERO = 0.0_DP
  real(DP), parameter :: DHALF = 0.5_DP
  real(DP), parameter :: DONE = 1.0_DP
  real(DP), parameter :: DTWO = 2.0_DP
  real(DP), parameter :: DEM1 = 1.0e-1_DP
  real(DP), parameter :: DEM2 = 1.0e-2_DP
  real(DP), parameter :: DEM3 = 1.0e-3_DP
  real(DP), parameter :: DEM4 = 1.0e-4_DP
  real(DP), parameter :: DEM5 = 1.0e-5_DP
  real(DP), parameter :: DEM6 = 1.0e-6_DP
  real(DP), parameter :: DEM8 = 1.0e-8_DP
  real(DP), parameter :: DPREC = EPSILON(1.0_DP)
  real(DP), parameter :: DSAME = 1.0e2_DP * DPREC
  ! -- enumerated verbosity level used by IMS (VDEBUG = 2)
  integer(I4B), parameter :: VDEBUG = 2
end module ConstantsModule

! ---------------------------------------------------------------------------
!  MathUtilModule -- is_close ported verbatim from MF6 6.4.4 src/Utilities/MathUtil.f90
! ---------------------------------------------------------------------------
module MathUtilModule
  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: DZERO, DSAME
  implicit none
  private
  public :: is_close
contains
  pure logical function is_close(a, b, rtol, atol, symmetric)
    real(DP), intent(in) :: a
    real(DP), intent(in) :: b
    real(DP), intent(in), optional :: rtol
    real(DP), intent(in), optional :: atol
    logical(LGP), intent(in), optional :: symmetric
    real(DP) :: lrtol, latol
    logical(LGP) :: lsymmetric
    ! check for exact equality
    if (a == b) then
      is_close = .true.
      return
    end if
    ! process optional arguments
    if (.not. present(rtol)) then
      lrtol = DSAME
    else
      lrtol = rtol
    end if
    if (.not. present(atol)) then
      latol = DZERO
    else
      latol = atol
    end if
    if (.not. present(symmetric)) then
      lsymmetric = .true.
    else
      lsymmetric = symmetric
    end if
    if (lsymmetric) then
      is_close = abs(a - b) <= max(lrtol * max(abs(a), abs(b)), latol)
    else
      is_close = (abs(a - b) <= (latol + lrtol * abs(b)))
    end if
  end function is_close
end module MathUtilModule

! ---------------------------------------------------------------------------
!  SimModule -- guardrail shim. store_error accumulates + prints; ustop prints
!  and STOPs. For Phase-1 isolation a print+STOP is correct; the integration
!  phase routes ustop() to USG-T's USTOP (utl7u1.f:2660).
! ---------------------------------------------------------------------------
!  PHASE 2 (USG-T wiring): store_error writes to the USG-T listing unit IOUT and
!  accumulates; ustop routes to USG-T's external USTOP (utl7u1.f:2660) via the
!  fixed-form bridge IMS_USTOP (a module procedure named `ustop` cannot call the
!  case-insensitively-identical external USTOP directly). The listing unit is
!  registered once by IMS7U1AR via set_sim_iout; reset_errors clears the counter
!  at the start of each solve so count_errors() reflects only the current call.
module SimModule
  use KindModule, only: I4B
  implicit none
  private
  public :: store_error, count_errors, ustop, deprecation_warning
  public :: set_sim_iout, reset_errors
  integer(I4B), save :: nerrors = 0
  integer(I4B), save :: sim_iout = 0   !< USG-T listing unit (0 => stdout only)
contains
  subroutine set_sim_iout(iu)
    integer(I4B), intent(in) :: iu
    sim_iout = iu
  end subroutine set_sim_iout

  subroutine reset_errors()
    nerrors = 0
  end subroutine reset_errors

  subroutine store_error(msg, terminate)
    character(len=*), intent(in) :: msg
    logical, intent(in), optional :: terminate
    nerrors = nerrors + 1
    if (sim_iout > 0) write (sim_iout, '(2a)') ' IMS ERROR: ', trim(msg)
    write (*, '(2a)') ' IMS ERROR: ', trim(msg)
    if (present(terminate)) then
      if (terminate) call ustop('IMS terminating on stored error')
    end if
  end subroutine store_error

  integer(I4B) function count_errors()
    count_errors = nerrors
  end function count_errors

  subroutine ustop(msg)
    character(len=*), intent(in), optional :: msg
    external :: IMS_USTOP   !< fixed-form bridge -> USG-T USTOP (utl7u1.f)
    if (sim_iout > 0 .and. present(msg)) &
      write (sim_iout, '(2a)') ' IMS STOP: ', trim(msg)
    if (present(msg)) then
      call IMS_USTOP(trim(msg))
    else
      call IMS_USTOP(' ')
    end if
  end subroutine ustop

  ! Imported by ims8linear (use SimModule, only: ..., deprecation_warning).
  ! Not exercised by the isolation driver; provided so the `only:` import resolves.
  subroutine deprecation_warning(cblock, cvar, cver, cmsg, iunit)
    character(len=*), intent(in) :: cblock
    character(len=*), intent(in) :: cvar
    character(len=*), intent(in) :: cver
    character(len=*), intent(in), optional :: cmsg
    integer(I4B), intent(in), optional :: iunit
    write (*, '(5a)') ' IMS DEPRECATION: ', trim(cblock), ' / ', trim(cvar), &
      ' '//trim(cver)
  end subroutine deprecation_warning
end module SimModule

! ---------------------------------------------------------------------------
!  MemoryHelperModule -- only create_mem_path is referenced (AR builds a path
!  string). With the memory manager removed the path is cosmetic; return a
!  simple "name/subcomponent" string.
! ---------------------------------------------------------------------------
module MemoryHelperModule
  implicit none
  private
  public :: create_mem_path
contains
  function create_mem_path(name, subcomponent) result(mem_path)
    character(len=*), intent(in) :: name
    character(len=*), intent(in), optional :: subcomponent
    character(len=:), allocatable :: mem_path
    if (present(subcomponent)) then
      mem_path = trim(name)//'/'//trim(subcomponent)
    else
      mem_path = trim(name)
    end if
  end function create_mem_path
end module MemoryHelperModule

! ---------------------------------------------------------------------------
!  BlockParserModule -- IMS only declares a private `parser` of this type and
!  calls parser%StoreErrorUnit() inside two error branches. The option-reading
!  path is bypassed in Phase 1, so a no-op stub type is sufficient.
! ---------------------------------------------------------------------------
module BlockParserModule
  use KindModule, only: I4B
  implicit none
  private
  public :: BlockParserType
  type :: BlockParserType
    integer(I4B) :: iuext = 0
  contains
    procedure :: StoreErrorUnit => bp_store_error_unit
  end type BlockParserType
contains
  subroutine bp_store_error_unit(this, terminate)
    class(BlockParserType), intent(inout) :: this
    logical, intent(in), optional :: terminate
    ! no-op: in isolation there is no input unit to report
  end subroutine bp_store_error_unit
end module BlockParserModule

! ---------------------------------------------------------------------------
!  ImsLinearSettingsModule -- data holder for the IMS linear options. In MF6
!  these are read from the IMS block via BlockParser; in Phase 1 the driver
!  sets them directly. Fields are plain scalars (no memory manager).
! ---------------------------------------------------------------------------
module ImsLinearSettingsModule
  use KindModule, only: DP, I4B
  implicit none
  private
  integer(I4B), public, parameter :: CG_METHOD = 1
  integer(I4B), public, parameter :: BCGS_METHOD = 2
  type, public :: ImsLinearSettingsType
    real(DP) :: dvclose = 0.0_DP    !< dependent variable closure criterion
    real(DP) :: rclose = 0.0_DP     !< residual closure criterion
    integer(I4B) :: icnvgopt = 0    !< convergence option
    integer(I4B) :: iter1 = 0       !< max. inner iterations
    integer(I4B) :: ilinmeth = 0    !< linear solver method (1=CG, 2=BCGS)
    integer(I4B) :: iscl = 0        !< scaling method
    integer(I4B) :: iord = 0        !< reordering method
    integer(I4B) :: north = 0       !< number of orthogonalizations
    real(DP) :: relax = 0.0_DP      !< relaxation factor
    integer(I4B) :: level = 0       !< nr. of preconditioner levels
    real(DP) :: droptol = 0.0_DP    !< drop tolerance for preconditioner
    integer(I4B) :: ifdparam = 0    !< complexity option
  end type ImsLinearSettingsType
end module ImsLinearSettingsModule

! ---------------------------------------------------------------------------
!  ConvergenceSummaryModule -- data holder the accelerators write per-model
!  convergence info into. Fields/shape match MF6 6.4.4 ConvergenceSummary.f90,
!  but arrays are ALLOCATABLE (no memory manager). The driver calls init() to
!  size them for (convnmod, nitermax).
! ---------------------------------------------------------------------------
module ConvergenceSummaryModule
  use KindModule, only: DP, I4B
  implicit none
  private
  public :: ConvergenceSummaryType
  type :: ConvergenceSummaryType
    integer(I4B) :: iter_cnt = 0                                  !< iteration counter within the timestep
    integer(I4B) :: convnmod = 1                                  !< number of models in the solution
    integer(I4B) :: nitermax = 0                                  !< max. nr. of inner iterations in a timestep
    integer(I4B), allocatable :: itinner(:)                       !< inner iteration number within each picard iteration
    integer(I4B), allocatable :: locdv(:)                         !< location of max dependent-variable change
    real(DP), allocatable :: dvmax(:)                             !< max dependent-variable change
    integer(I4B), allocatable :: locdr(:)                         !< location of max flow change
    real(DP), allocatable :: drmax(:)                             !< max flow change
    integer(I4B), allocatable :: convlocdv(:, :)                  !< per-model max dv-change location
    real(DP), allocatable :: convdvmax(:, :)                      !< per-model max dv-change
    integer(I4B), allocatable :: convlocdr(:, :)                  !< per-model max flow-change location
    real(DP), allocatable :: convdrmax(:, :)                      !< per-model max flow-change
  contains
    procedure :: init => cs_init
    procedure :: destroy => cs_destroy
  end type ConvergenceSummaryType
contains
  subroutine cs_init(this, convnmod, nitermax)
    class(ConvergenceSummaryType), intent(inout) :: this
    integer(I4B), intent(in) :: convnmod
    integer(I4B), intent(in) :: nitermax
    this%convnmod = convnmod
    this%nitermax = nitermax
    this%iter_cnt = 0
    allocate (this%itinner(nitermax))
    allocate (this%locdv(convnmod))
    allocate (this%dvmax(convnmod))
    allocate (this%locdr(convnmod))
    allocate (this%drmax(convnmod))
    allocate (this%convlocdv(convnmod, nitermax))
    allocate (this%convdvmax(convnmod, nitermax))
    allocate (this%convlocdr(convnmod, nitermax))
    allocate (this%convdrmax(convnmod, nitermax))
    this%itinner = 0
    this%locdv = 0
    this%dvmax = 0.0_DP
    this%locdr = 0
    this%drmax = 0.0_DP
    this%convlocdv = 0
    this%convdvmax = 0.0_DP
    this%convlocdr = 0
    this%convdrmax = 0.0_DP
  end subroutine cs_init

  subroutine cs_destroy(this)
    class(ConvergenceSummaryType), intent(inout) :: this
    if (allocated(this%itinner)) deallocate (this%itinner)
    if (allocated(this%locdv)) deallocate (this%locdv)
    if (allocated(this%dvmax)) deallocate (this%dvmax)
    if (allocated(this%locdr)) deallocate (this%locdr)
    if (allocated(this%drmax)) deallocate (this%drmax)
    if (allocated(this%convlocdv)) deallocate (this%convlocdv)
    if (allocated(this%convdvmax)) deallocate (this%convdvmax)
    if (allocated(this%convlocdr)) deallocate (this%convlocdr)
    if (allocated(this%convdrmax)) deallocate (this%convdrmax)
  end subroutine cs_destroy
end module ConvergenceSummaryModule

! ---------------------------------------------------------------------------
!  MatrixBaseModule -- IMS only uses this to receive the CSR system in
!  IMSLINEAR_AR via `class(MatrixBaseType), pointer :: matrix` and one call
!  `matrix%get_aij(ia, ja, amat)`. A concrete type holding pointers to the
!  caller's CSR arrays, with a get_aij that associates the pointers, reproduces
!  exactly that behaviour without SparseModule / VectorBaseModule.
! ---------------------------------------------------------------------------
module MatrixBaseModule
  use KindModule, only: DP, I4B
  implicit none
  private
  public :: MatrixBaseType
  type :: MatrixBaseType
    integer(I4B), dimension(:), pointer, contiguous :: ia => null()
    integer(I4B), dimension(:), pointer, contiguous :: ja => null()
    real(DP), dimension(:), pointer, contiguous :: amat => null()
  contains
    procedure :: get_aij => mat_get_aij
  end type MatrixBaseType
contains
  subroutine mat_get_aij(this, ia, ja, amat)
    class(MatrixBaseType), intent(in) :: this
    integer(I4B), dimension(:), pointer, contiguous :: ia
    integer(I4B), dimension(:), pointer, contiguous :: ja
    real(DP), dimension(:), pointer, contiguous :: amat
    ia => this%ia
    ja => this%ja
    amat => this%amat
  end subroutine mat_get_aij
end module MatrixBaseModule

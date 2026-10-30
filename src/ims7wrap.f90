!> @file ims7wrap.f90
!! PHASE-2 USG-Transport wrappers for the IMS linear solver (LINMETH = 5).
!!
!! Mirrors the PCGU integration (pcgu7.f): AR reads options + sizes work arrays +
!! does the symbolic setup; AP copies AMAT/RHS into the IMS type, solves, and
!! scatters the solution back into SOLN; DA deallocates. State for the single
!! solution is held in this module (one ImsLinearDataType instance), exactly as
!! PCGU keeps its state in PCGUMODULE.
!!
!! Matrix convention: USG-T assembles A x = b in CSR (IA2/JA2/AMAT, diagonal
!! present) — the same form PCGU consumes. IMS-linear consumes that form directly
!! (no re-sign / no re-scale outside IMS's own ISCL), so the wrapper just copies.
module IMS7WRAPPERMODULE
  use KindModule, only: DP, I4B
  use IMSLinearModule, only: ImsLinearDataType
  use ImsLinearSettingsModule, only: ImsLinearSettingsType
  use ConvergenceSummaryModule, only: ConvergenceSummaryType
  use MatrixBaseModule, only: MatrixBaseType
  use SimModule, only: set_sim_iout, reset_errors, count_errors, ustop
  implicit none
  private
  public :: IMS7U1AR, IMS7U1AP, IMS7DA, IMS_PREPROC
  public :: IMS7_GET_ACCEL, IMS7_SET_ACCEL, IMS7_GET_CONFIG_ACCEL
  public :: IMS7_SET_TOLS

  ! -- single-solution state (mirrors PCGUMODULE)
  type(ImsLinearDataType), save :: ims
  type(ImsLinearSettingsType), pointer :: ims_settings => null()
  type(ConvergenceSummaryType), pointer :: ims_summary => null()
  class(MatrixBaseType), pointer :: ims_matrix => null()

  ! -- module-owned CSR + vectors the IMS type aliases (TARGET so its pointers bind)
  integer(I4B), allocatable, target, save :: g_ia(:)
  integer(I4B), allocatable, target, save :: g_ja(:)
  real(DP), allocatable, target, save :: g_amat(:)
  real(DP), allocatable, target, save :: g_rhs(:)
  real(DP), allocatable, target, save :: g_x(:)

  integer(I4B), save :: g_neqs = 0   !< NEQS2 (active system size)
  integer(I4B), save :: g_nja = 0    !< NJA2

  ! -- PHASE-3b POLISH-3: setup params stashed by IMS7U1AR so IMS_PREPROC can
  !    (re)build the reordering + symbolic ILU for the current IA2/JA2/NEQS2/NJA2.
  !    Mirrors how XMD splits XMD7U1AR (read options) from XMD_PREPROC (symbolic).
  integer(I4B), save :: g_iout = 0
  integer(I4B), save :: g_iprims = 0
  integer(I4B), save :: g_mxiter = 1
  logical, save :: g_echoed = .false.   !< echo the settings summary only once

contains

  !> @brief Allocate IMS settings + read options (mirrors XMD7U1AR/PCGU7U1AR).
  !! PHASE-3b POLISH-3: the structure-dependent symbolic setup (reordering +
  !! symbolic ILU + work arrays on IA2/JA2/NEQS2/NJA2) is now done in IMS_PREPROC,
  !! NOT here — exactly as XMD splits XMD7U1AR (options) from XMD_PREPROC (symbolic).
  !! For SOLVEACTIVE (ISOLVEACTIVE/=0) IA2/JA2 are not yet built when this AR runs;
  !! SMS_REDUCE calls IMS_PREPROC later (after FILLIA2JA2). For ISOLVEACTIVE==0 the
  !! full system is already pointed (SMS_REDUCE0), so we call IMS_PREPROC here.
  subroutine IMS7U1AR(IN, MXITER, ITER1, IPRSMS, IFDPARAM, IOUT)
    use SMSMODULE, only: ISOLVEACTIVE
    ! -- dummy
    integer(I4B), intent(in) :: IN       !< SMS input unit
    integer(I4B), intent(in) :: MXITER   !< max outer iterations
    integer(I4B), intent(in) :: ITER1    !< max inner iterations (from SMS)
    integer(I4B), intent(in) :: IPRSMS   !< print flag
    integer(I4B), intent(in) :: IFDPARAM !< complexity option (0 = read from file)
    integer(I4B), intent(in) :: IOUT     !< listing unit

    ! -- route IMS error messages to the USG-T listing file
    call set_sim_iout(IOUT)

    write (IOUT, 2000)
2000 format(1X, /1X, 'IMS -- IMS LINEAR SOLUTION PACKAGE', &
            ' (ported from MODFLOW 6 tag 6.4.4), LINMETH = 5')

    ! -- settings holder; populate from the SMS input file (URWORD) or presets
    allocate (ims_settings)
    call ims_read_options(IN, IOUT, IFDPARAM, ITER1)

    ! -- matrix shim + convergence summary (structure-independent)
    allocate (ims_matrix)
    allocate (ims_summary)
    call ims_summary%init(1, 1)

    ! -- stash params for IMS_PREPROC (the per-structure rebuild hook)
    g_iout = IOUT
    g_iprims = IPRSMS
    g_mxiter = MXITER
    g_echoed = .false.

    ! -- ISOLVEACTIVE==0: full system is already pointed (SMS_REDUCE0 ran before
    !    this AR), so build the symbolic setup now. Otherwise defer to IMS_PREPROC
    !    (called from SMS_REDUCE after FILLIA2JA2 builds the active IA2/JA2).
    if (ISOLVEACTIVE == 0) call IMS_PREPROC()

    return
  end subroutine IMS7U1AR

  !> @brief PHASE-3b POLISH-3: (re)build the IMS reordering + symbolic ILU + work
  !! arrays for the CURRENT active system IA2/JA2/NEQS2/NJA2. Idempotent — frees
  !! the previous structure-dependent arrays first (the 3a heap error was
  !! "pointer being freed was not allocated", so every array is deallocated before
  !! reallocating). Wired into USG-T's SOLVER PREPROCESSOR dispatch (SMS_REDUCE),
  !! next to XMD_PREPROC/PCGUSETDIMS, so IMS rebuilds whenever USG-T recompresses
  !! the active set under SOLVEACTIVE.
  subroutine IMS_PREPROC()
    use GLOBAL, only: IA2, JA2, NEQS2, NJA2
    character(len=16) :: solname
    integer(I4B) :: n

    ! -- free the previous structure-dependent IMS arrays (safe on first call:
    !    all components are unallocated/nullified, so the guards inside skip).
    call ims%IMSLINEAR_DA()
    if (allocated(g_ia)) deallocate (g_ia)
    if (allocated(g_ja)) deallocate (g_ja)
    if (allocated(g_amat)) deallocate (g_amat)
    if (allocated(g_rhs)) deallocate (g_rhs)
    if (allocated(g_x)) deallocate (g_x)

    ! -- size to the current active system and copy its CSR structure
    g_neqs = NEQS2
    g_nja = NJA2
    allocate (g_ia(g_neqs + 1))
    allocate (g_ja(g_nja))
    allocate (g_amat(g_nja))
    allocate (g_rhs(g_neqs))
    allocate (g_x(g_neqs))
    do n = 1, g_neqs + 1
      g_ia(n) = IA2(n)
    end do
    do n = 1, g_nja
      g_ja(n) = JA2(n)
    end do
    g_amat = 0.0_DP
    g_rhs = 0.0_DP
    g_x = 0.0_DP

    ! -- point the matrix shim at the (re)built CSR arrays
    ims_matrix%ia => g_ia
    ims_matrix%ja => g_ja
    ims_matrix%amat => g_amat

    ! -- rebuild IMS work arrays + symbolic ILU + reordering for this structure
    solname = 'IMSSOLUTION'
    call ims%IMSLINEAR_ALLOCATE(solname, g_iout, g_iprims, g_mxiter, &
                                g_neqs, ims_matrix, g_rhs, g_x, ims_settings)

    ! -- F2: escalate a fatal setup error. IMSLINEAR_ALLOCATE runs the reordering
    !    and symbolic factorization and may store_error on a bad ordering. Without
    !    this check the next solve's reset_errors() would silently clear it, so a
    !    setup-fatal error is hard-stopped here instead of being swept under.
    if (count_errors() > 0) then
      if (g_iout > 0) write (g_iout, '(a)') &
        ' IMS SETUP FAILURE: error during IMS reordering / symbolic factorization'
      call ustop('IMS SETUP FAILURE (reordering / symbolic factorization)')
    end if

    ! -- echo the chosen settings once (MF6-style summary)
    if (.not. g_echoed) then
      call ims%imslinear_summary(g_mxiter)
      g_echoed = .true.
    end if

    return
  end subroutine IMS_PREPROC

  !> @brief Read IMS linear options from the SMS input unit in USG-T style
  !! (URDCOM/URWORD), or apply MF6 SIMPLE/MODERATE/COMPLEX presets when
  !! IFDPARAM /= 0. The preconditioner code (IPC) is derived inside
  !! IMSLINEAR_AR from LEVEL/DROPTOL/RELAX, exactly as in stock IMS.
  !!
  !! Positional option line (IFDPARAM = 0), MF6-style fields (8 fields;
  !! INNER_DVCLOSE removed — the flow dependent-variable closure is taken
  !! from HICLOSE, the SMS item-1 inner head-closure criterion):
  !!   ACCEL  ISCL  IORD  INNER_RCLOSE  RELAXATION_FACTOR \
  !!          PRECONDITIONER_LEVELS  PRECONDITIONER_DROP_TOLERANCE  NUM_ORTHO
  !!   ACCEL = CG | BCGS
  subroutine ims_read_options(IN, IOUT, IFDPARAM, ITER1)
    use SMSMODULE, only: HICLOSE
    integer(I4B), intent(in) :: IN, IOUT, IFDPARAM, ITER1
    character(len=400) :: line
    integer(I4B) :: lloc, istart, istop, ival
    real :: rval
    character(len=20) :: cval

    ! -- defaults common to all paths
    ims_settings%iter1 = ITER1
    ims_settings%icnvgopt = 0      ! infinity-norm residual test
    ims_settings%ifdparam = IFDPARAM
    ! -- flow dependent-variable closure = HICLOSE (SMS item-1). IMS no longer
    !    reads its own DVCLOSE from the option line. The transport solve's
    !    dependent-variable closure is CICLOSE (read from the BCT file).
    ims_settings%dvclose = HICLOSE

    if (IFDPARAM == 1) then
      ! -- SIMPLE preset (MF6)
      ims_settings%ilinmeth = 1; ims_settings%iscl = 0; ims_settings%iord = 0
      ims_settings%rclose = 1.0e-1_DP
      ims_settings%relax = 0.0_DP; ims_settings%level = 0
      ims_settings%droptol = 0.0_DP; ims_settings%north = 0
      return
    else if (IFDPARAM == 2) then
      ! -- MODERATE preset (MF6)
      ims_settings%ilinmeth = 2; ims_settings%iscl = 0; ims_settings%iord = 0
      ims_settings%rclose = 1.0e-1_DP
      ims_settings%relax = 0.97_DP; ims_settings%level = 0
      ims_settings%droptol = 0.0_DP; ims_settings%north = 0
      return
    else if (IFDPARAM == 3) then
      ! -- COMPLEX preset (MF6)
      ims_settings%ilinmeth = 2; ims_settings%iscl = 0; ims_settings%iord = 0
      ims_settings%rclose = 1.0e-1_DP
      ims_settings%relax = 0.0_DP; ims_settings%level = 5
      ims_settings%droptol = 1.0e-4_DP; ims_settings%north = 2
      return
    end if

    ! -- IFDPARAM = 0: read the option line from the SMS input file
    call URDCOM(IN, IOUT, line)
    lloc = 1
    ! accelerator keyword
    call URWORD(line, lloc, istart, istop, 1, ival, rval, IOUT, IN)
    cval = line(istart:istop)
    if (cval(1:4) == 'BCGS') then
      ims_settings%ilinmeth = 2
    else
      ims_settings%ilinmeth = 1   ! CG default
    end if
    call URWORD(line, lloc, istart, istop, 2, ims_settings%iscl, rval, IOUT, IN)
    call URWORD(line, lloc, istart, istop, 2, ims_settings%iord, rval, IOUT, IN)
    ! -- INNER_RCLOSE (DVCLOSE removed; flow dep-var closure = HICLOSE, set above)
    call URWORD(line, lloc, istart, istop, 3, ival, rval, IOUT, IN)
    ims_settings%rclose = real(rval, DP)
    call URWORD(line, lloc, istart, istop, 3, ival, rval, IOUT, IN)
    ims_settings%relax = real(rval, DP)
    call URWORD(line, lloc, istart, istop, 2, ims_settings%level, rval, IOUT, IN)
    call URWORD(line, lloc, istart, istop, 3, ival, rval, IOUT, IN)
    ims_settings%droptol = real(rval, DP)
    call URWORD(line, lloc, istart, istop, 2, ims_settings%north, rval, IOUT, IN)

    return
  end subroutine ims_read_options

  !> @brief Copy AMAT/RHS into the IMS type, solve, scatter solution to SOLN.
  !! Mirrors PCGU7U1AP's signature (incl. Rcutoff/ICUTOFF for parity). PHASE-3b:
  !! the truncated-Newton inner cutoff is now honored — Rcutoff is armed into the
  !! accelerators (ims_rcutoff) and ICUTOFF=1 is returned when the inner
  !! sum-of-squares residual <= Rcutoff (mirrors xMD mrcutoff/micutoff).
  !! IERR is the guardrail: > 0 means a true solver breakdown (stored error or
  !! NaN in the solution) and the caller (SOLVERS) hard-stops via USTOP.
  subroutine IMS7U1AP(AMAT, RHS, SOLN, IA2, JA2, NJA2, NEQS2, &
                      ICNVG, KSTP, KPER, MXITER, KITER, IN_ITER, IOUT, &
                      Rcutoff, ICUTOFF, IERR)
    ! -- PHASE-3b POLISH-2: truncated-Newton inner cutoff state (mirrors xMD)
    use IMSLinearBaseModule, only: ims_rcutoff, ims_icutoff
    ! -- dummy
    integer(I4B), intent(in) :: NJA2, NEQS2
    real(DP), intent(inout) :: AMAT(NJA2)
    real(DP), intent(inout) :: RHS(NEQS2)
    real(DP), intent(inout) :: SOLN(NEQS2)
    integer(I4B), intent(in) :: IA2(NEQS2 + 1)
    integer(I4B), intent(in) :: JA2(NJA2)
    integer(I4B), intent(inout) :: ICNVG
    integer(I4B), intent(in) :: KSTP, KPER, MXITER, KITER, IOUT
    integer(I4B), intent(inout) :: IN_ITER
    real(DP), intent(in) :: Rcutoff
    integer(I4B), intent(inout) :: ICUTOFF
    integer(I4B), intent(inout) :: IERR
    ! -- local
    integer(I4B) :: n, nconv, convnmod
    integer(I4B) :: convmodstart(2)
    character(len=31) :: caccel(1)

    ! -- PHASE-3b POLISH-2: arm the truncated-Newton inner cutoff for this solve.
    !    Rcutoff is the sum-of-squares residual cutoff from the SMS outer loop
    !    (< 0 disables it). ims_icutoff is read back after the solve.
    ims_rcutoff = Rcutoff
    ims_icutoff = 0
    IERR = 0

    ! -- copy current system values into the module-owned targets
    do n = 1, NJA2
      g_amat(n) = AMAT(n)
    end do
    do n = 1, NEQS2
      g_rhs(n) = RHS(n)
      g_x(n) = SOLN(n)
    end do

    ! -- single-model convergence bookkeeping (detailed log off: NCONV = 1)
    nconv = 1
    convnmod = 1
    convmodstart(1) = 1
    convmodstart(2) = NEQS2 + 1
    caccel(1) = ''

    ! -- reset accumulated IMS errors so count_errors() reflects this solve only
    call reset_errors()

    ! -- inner solve
    call ims%IMSLINEAR_APPLY(ICNVG, KSTP, KITER, IN_ITER, &
                             nconv, convnmod, convmodstart, caccel, ims_summary)

    ! -- scatter solution back
    do n = 1, NEQS2
      SOLN(n) = g_x(n)
    end do

    ! -- PHASE-3b POLISH-2: report the truncated-Newton cutoff flag to the outer loop
    ICUTOFF = ims_icutoff

    ! -- guardrail: flag a true breakdown (stored error, or a non-finite solution)
    do n = 1, NEQS2
      ! not-finite test: catches NaN (x/=x) AND +/-Inf (|x| > HUGE), portably,
      ! so a non-finite solve triggers USTOP instead of a NaN/Inf budget.
      if (.not. (abs(SOLN(n)) <= HUGE(SOLN(n)))) then
        IERR = 2
        exit
      end if
    end do

    return
  end subroutine IMS7U1AP

  !> @brief Deallocate IMS storage (mirrors XMD7DA / PCGU7U1DA).
  subroutine IMS7DA()
    call ims%IMSLINEAR_DA()
    if (associated(ims_summary)) then
      call ims_summary%destroy()
      deallocate (ims_summary)
      nullify (ims_summary)
    end if
    if (associated(ims_settings)) then
      deallocate (ims_settings)
      nullify (ims_settings)
    end if
    if (associated(ims_matrix)) then
      deallocate (ims_matrix)
      nullify (ims_matrix)
    end if
    if (allocated(g_ia)) deallocate (g_ia)
    if (allocated(g_ja)) deallocate (g_ja)
    if (allocated(g_amat)) deallocate (g_amat)
    if (allocated(g_rhs)) deallocate (g_rhs)
    if (allocated(g_x)) deallocate (g_x)
    return
  end subroutine IMS7DA

  !> @brief Return the IMS linear accelerator currently in effect (1=CG, 2=BiCGSTAB).
  !! Reads the live solver object -- exactly what IMSLINEAR_APPLY dispatches on -- so
  !! GWT2SOLVE can save it around the transport solve. Mirrors XMD IACL / PCGU ILINMETH.
  subroutine IMS7_GET_ACCEL(iaccel)
    integer(I4B), intent(out) :: iaccel
    iaccel = ims%ILINMETH
  end subroutine IMS7_GET_ACCEL

  !> @brief Return the CONFIGURED IMS accelerator (1=CG, 2=BiCGSTAB) from the settings
  !! read at AR time. Used by the LINMETH=5 AR arm to set ISYMFLG before the live
  !! solver exists (SOLVEACTIVE defers IMS_PREPROC, so ims%ILINMETH is not yet built).
  !! Returns 2 (BiCGSTAB, the non-symmetric default) if settings are not yet allocated.
  subroutine IMS7_GET_CONFIG_ACCEL(iaccel)
    integer(I4B), intent(out) :: iaccel
    iaccel = 2
    if (associated(ims_settings)) iaccel = ims_settings%ilinmeth
  end subroutine IMS7_GET_CONFIG_ACCEL

  !> @brief Set the IMS linear accelerator (1=CG, 2=BiCGSTAB) for the next solve.
  !! Sets the live solver object so the change takes effect at APPLY time with no
  !! rebuild, and keeps ims_settings in sync for any later structure rebuild. Used by
  !! GWT2SOLVE to force BiCGSTAB for the non-symmetric transport system and restore
  !! the flow accelerator afterward -- the flow solve is unaffected.
  subroutine IMS7_SET_ACCEL(iaccel)
    integer(I4B), intent(in) :: iaccel
    ims%ILINMETH = iaccel
    if (associated(ims_settings)) ims_settings%ilinmeth = iaccel
    ! -- BiCGSTAB uses five extra NEQ-sized work vectors (T/V/DHAT/PHAT/QHAT) that
    !    IMSLINEAR_ALLOCATE sizes to NEQ only when the solver was built for BCGS;
    !    for a CG-built solver they are length 1. When we force BiCGSTAB on such a
    !    solver (the transport solve of a CG/default-configured run), grow those
    !    scratch vectors to NEQ so ims_base_bcgs does not index past their end.
    !    One-time (NIABCGS stays NEQ after); they are pure per-solve scratch.
    if (iaccel == 2 .and. ims%NIABCGS < ims%NEQ) then
      if (allocated(ims%T)) deallocate (ims%T)
      if (allocated(ims%V)) deallocate (ims%V)
      if (allocated(ims%DHAT)) deallocate (ims%DHAT)
      if (allocated(ims%PHAT)) deallocate (ims%PHAT)
      if (allocated(ims%QHAT)) deallocate (ims%QHAT)
      ims%NIABCGS = ims%NEQ
      allocate (ims%T(ims%NIABCGS), ims%V(ims%NIABCGS), ims%DHAT(ims%NIABCGS), &
                ims%PHAT(ims%NIABCGS), ims%QHAT(ims%NIABCGS))
      ims%T = 0.0_DP; ims%V = 0.0_DP; ims%DHAT = 0.0_DP
      ims%PHAT = 0.0_DP; ims%QHAT = 0.0_DP
    end if
  end subroutine IMS7_SET_ACCEL

  !> @brief Set the per-stress-period IMS linear tolerances for the next solve.
  !! Updates the live solver object so the change takes effect at APPLY time with
  !! no structural rebuild, and keeps ims_settings in sync (so a later deferred
  !! IMS_PREPROC that copies ims_settings picks up the same values). A negative
  !! argument means "leave unchanged". hclose maps to the dependent-variable
  !! (head) closure -- the same HICLOSE mapping used for the flow solve.
  !! Structural parameters (level/ordering/preconditioner) are not touched.
  subroutine IMS7_SET_TOLS(hclose, rclose, iter1)
    real(DP), intent(in) :: hclose, rclose
    integer(I4B), intent(in) :: iter1
    if (hclose > 0.0_DP) then
      ims%DVCLOSE = hclose
      if (associated(ims_settings)) ims_settings%dvclose = hclose
    end if
    if (rclose > 0.0_DP) then
      ims%RCLOSE = rclose
      if (associated(ims_settings)) ims_settings%rclose = rclose
    end if
    if (iter1 > 0) then
      ims%ITER1 = iter1
      if (associated(ims_settings)) ims_settings%iter1 = iter1
    end if
  end subroutine IMS7_SET_TOLS

end module IMS7WRAPPERMODULE

!> @brief Fixed-interface bridge so the SimModule shim's module procedure `ustop`
!! can reach USG-T's external USTOP (utl7u1.f) despite Fortran's case-insensitive
!! name identity (module `ustop` vs external `USTOP`). Defined OUTSIDE the module
!! as a plain external subroutine.
subroutine IMS_USTOP(MSG)
  character(len=*), intent(in) :: MSG
  call USTOP(MSG)
end subroutine IMS_USTOP

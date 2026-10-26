MODULE IMSLinearModule

  use KindModule, only: DP, I4B
  use ConstantsModule, only: LINELENGTH, LENSOLUTIONNAME, LENMEMPATH, &
                             IZERO, DZERO, DPREC, DSAME, &
                             DEM8, DEM6, DEM5, DEM4, DEM3, DEM2, DEM1, &
                             DHALF, DONE, DTWO, &
                             VDEBUG
  use IMSLinearBaseModule, only: ims_base_cg, ims_base_bcgs, &
                                 ims_base_pccrs, ims_base_calc_order, &
                                 ims_base_scale, ims_base_pcu, &
                                 ims_base_residual
  use BlockParserModule, only: BlockParserType
  use MatrixBaseModule
  use ConvergenceSummaryModule
  use ImsLinearSettingsModule

  IMPLICIT NONE
  private

  TYPE, PUBLIC :: ImsLinearDataType
    character(len=LENMEMPATH) :: memoryPath !< the path for storing variables in the memory manager
    ! -- PHASE-1 ISOLATION: memory-manager POINTER scalars converted to plain scalars
    integer(I4B) :: iout = 0 !< simulation listing file unit
    integer(I4B) :: IPRIMS = 0 !< print flag
    ! input variables (copied from settings structure)
    real(DP) :: DVCLOSE = 0.0_DP !< dependent variable closure criterion
    real(DP) :: RCLOSE = 0.0_DP !< residual closure criterion
    integer(I4B) :: ICNVGOPT = 0 !< convergence option
    integer(I4B) :: ITER1 = 0 !< max. iterations
    integer(I4B) :: ILINMETH = 0 !< linear solver method
    integer(I4B) :: iSCL = 0 !< scaling method
    integer(I4B) :: IORD = 0 !< reordering method
    integer(I4B) :: NORTH = 0 !< number of orthogonalizations
    real(DP) :: RELAX = 0.0_DP !< relaxation factor
    integer(I4B) :: LEVEL = 0 !< nr. of preconditioner levels
    real(DP) :: DROPTOL = 0.0_DP !< drop tolerance for preconditioner
    !
    integer(I4B) :: IPC = 0 !< preconditioner flag
    integer(I4B) :: IACPC = 0 !< preconditioner CRS row pointers
    integer(I4B) :: NITERC = 0 !<
    integer(I4B) :: NIABCGS = 0 !< size of working vectors for BCGS linear accelerator
    integer(I4B) :: NIAPC = 0 !< preconditioner number of rows
    integer(I4B) :: NJAPC = 0 !< preconditioner number of non-zero entries
    real(DP) :: EPFACT = 0.0_DP !< factor for decreasing convergence criteria in seubsequent Picard iterations
    real(DP) :: L2NORM0 = 0.0_DP !< initial L2 norm
    ! -- ilut variables
    integer(I4B) :: NJLU = 0 !< length of jlu work vector
    integer(I4B) :: NJW = 0 !< length of jw work vector
    integer(I4B) :: NWLU = 0 !< length of wlu work vector
    ! -- pointers to solution variables (aliases to caller storage -> stay POINTER)
    integer(I4B) :: NEQ = 0 !< number of equations (rows in matrix)
    integer(I4B) :: NJA = 0 !< number of non-zero values in amat
    integer(I4B), dimension(:), pointer, contiguous :: IA => NULL() !< position of start of each row
    integer(I4B), dimension(:), pointer, contiguous :: JA => NULL() !< column pointer
    real(DP), dimension(:), pointer, contiguous :: AMAT => NULL() !< coefficient matrix
    real(DP), dimension(:), pointer, contiguous :: RHS => NULL() !< right-hand side of equation
    real(DP), dimension(:), pointer, contiguous :: X => NULL() !< dependent variable
    ! VECTORS -- owned work arrays: POINTER -> ALLOCATABLE
    real(DP), ALLOCATABLE, DIMENSION(:) :: DSCALE !< scaling factor
    real(DP), ALLOCATABLE, DIMENSION(:) :: DSCALE2 !< unscaling factor
    integer(I4B), ALLOCATABLE, DIMENSION(:) :: IAPC !< position of start of each row in preconditioner matrix
    integer(I4B), ALLOCATABLE, DIMENSION(:) :: JAPC !< preconditioner matrix column pointer
    real(DP), ALLOCATABLE, DIMENSION(:) :: APC !< preconditioner coefficient matrix
    integer(I4B), ALLOCATABLE, DIMENSION(:) :: LORDER !< reordering mapping
    integer(I4B), ALLOCATABLE, DIMENSION(:) :: IORDER !< mapping to restore reordered matrix
    ! reordered system: aliased by IA0/JA0/A0 -> stay POINTER (allocated with ALLOCATE)
    integer(I4B), POINTER, DIMENSION(:), CONTIGUOUS :: IARO => NULL() !< position of start of each row in reordered matrix
    integer(I4B), POINTER, DIMENSION(:), CONTIGUOUS :: JARO => NULL() !< reordered matrix column pointer
    real(DP), POINTER, DIMENSION(:), CONTIGUOUS :: ARO => NULL() !< reordered coefficient matrix
    ! WORKING ARRAYS -- owned: POINTER -> ALLOCATABLE
    integer(I4B), ALLOCATABLE, DIMENSION(:) :: IW !< integer working array
    real(DP), ALLOCATABLE, DIMENSION(:) :: W !< real working array
    integer(I4B), ALLOCATABLE, DIMENSION(:) :: ID !< integer working array
    real(DP), ALLOCATABLE, DIMENSION(:) :: D !< real working array
    real(DP), ALLOCATABLE, DIMENSION(:) :: P !< real working array
    real(DP), ALLOCATABLE, DIMENSION(:) :: Q !< real working array
    real(DP), ALLOCATABLE, DIMENSION(:) :: Z !< real working array
    ! BICGSTAB WORKING ARRAYS -- owned: POINTER -> ALLOCATABLE
    real(DP), ALLOCATABLE, DIMENSION(:) :: T !< BICGSTAB real working array
    real(DP), ALLOCATABLE, DIMENSION(:) :: V !< BICGSTAB real working array
    real(DP), ALLOCATABLE, DIMENSION(:) :: DHAT !< BICGSTAB real working array
    real(DP), ALLOCATABLE, DIMENSION(:) :: PHAT !< BICGSTAB real working array
    real(DP), ALLOCATABLE, DIMENSION(:) :: QHAT !< rBICGSTAB eal working array
    ! POINTERS FOR USE WITH BOTH ORIGINAL AND RCM ORDERINGS (aliases -> stay POINTER)
    integer(I4B), POINTER, DIMENSION(:), CONTIGUOUS :: IA0 => NULL() !< pointer to current CRS row pointers
    integer(I4B), POINTER, DIMENSION(:), CONTIGUOUS :: JA0 => NULL() !< pointer to current CRS column pointers
    real(DP), POINTER, DIMENSION(:), CONTIGUOUS :: A0 => NULL() !< pointer to current coefficient matrix
    ! ILUT WORKING ARRAYS -- owned: POINTER -> ALLOCATABLE
    integer(I4B), ALLOCATABLE, DIMENSION(:) :: JLU !< ilut integer working array
    integer(I4B), ALLOCATABLE, DIMENSION(:) :: JW !< ilut integer working array
    real(DP), ALLOCATABLE, DIMENSION(:) :: WLU !< ilut real working array

    ! PROCEDURES (METHODS)
  CONTAINS
    PROCEDURE :: IMSLINEAR_ALLOCATE => imslinear_ar
    procedure :: imslinear_summary
    PROCEDURE :: IMSLINEAR_APPLY => imslinear_ap
    procedure :: IMSLINEAR_DA => imslinear_da
    procedure, private :: allocate_scalars
    ! -- PRIVATE PROCEDURES
    PROCEDURE, PRIVATE :: SET_IMSLINEAR_INPUT => imslinear_set_input
  END TYPE ImsLinearDataType

CONTAINS

  !> @ brief Allocate storage and read data
    !!
    !!  Allocate storage for linear accelerators and read data
    !!
  !<
  SUBROUTINE imslinear_ar(this, NAME, IOUT, IPRIMS, MXITER, &
                          NEQ, matrix, RHS, X, linear_settings)
    ! -- modules
    use MemoryHelperModule, only: create_mem_path
    use SimModule, only: store_error, count_errors, &
                         deprecation_warning
    ! -- dummy variables
    CLASS(ImsLinearDataType), INTENT(INOUT) :: this !< ImsLinearDataType instance
    CHARACTER(LEN=LENSOLUTIONNAME), INTENT(IN) :: NAME !< solution name
    integer(I4B), INTENT(IN) :: IOUT !< simulation listing file unit
    integer(I4B), TARGET, INTENT(IN) :: IPRIMS !< print option
    integer(I4B), INTENT(IN) :: MXITER !< maximum outer iterations
    integer(I4B), TARGET, INTENT(IN) :: NEQ !< number of equations
    class(MatrixBaseType), pointer :: matrix
    real(DP), DIMENSION(NEQ), TARGET, INTENT(INOUT) :: RHS !< right-hand side
    real(DP), DIMENSION(NEQ), TARGET, INTENT(INOUT) :: X !< dependent variables
    type(ImsLinearSettingsType), pointer :: linear_settings !< the settings form the IMS file
    ! -- local variables
    character(len=LINELENGTH) :: errmsg
    integer(I4B) :: i, n
    integer(I4B) :: i0
    integer(I4B) :: iscllen, iolen
    integer(I4B) :: ijlu
    integer(I4B) :: ijw
    integer(I4B) :: iwlu
    integer(I4B) :: iwk
    !
    ! -- DEFINE NAME
    this%memoryPath = create_mem_path(name, 'IMSLINEAR')
    !
    ! -- COPY IMS settings (plain scalars; PHASE-1 isolation)
    this%DVCLOSE = linear_settings%dvclose
    this%RCLOSE = linear_settings%rclose
    this%ICNVGOPT = linear_settings%icnvgopt
    this%ITER1 = linear_settings%iter1
    this%ILINMETH = linear_settings%ilinmeth
    this%iSCL = linear_settings%iscl
    this%IORD = linear_settings%iord
    this%NORTH = linear_settings%north
    this%RELAX = linear_settings%relax
    this%LEVEL = linear_settings%level
    this%DROPTOL = linear_settings%droptol
    !
    ! -- SET POINTERS TO SOLUTION STORAGE
    this%IPRIMS = IPRIMS
    this%NEQ = NEQ
    call matrix%get_aij(this%IA, this%JA, this%AMAT)
    this%NJA = size(this%AMAT)
    this%RHS => RHS
    this%X => X
    !
    ! -- ALLOCATE SCALAR VARIABLES
    call this%allocate_scalars()
    !
    ! -- initialize iout
    this%iout = iout
    !
    ! -- DEFAULT VALUES
    this%IPC = 0
    !
    this%IACPC = 0
    !
    ! -- PRINT A MESSAGE IDENTIFYING IMSLINEAR SOLVER PACKAGE
    write (iout, 2000)
02000 FORMAT(1X, /1X, 'IMSLINEAR -- UNSTRUCTURED LINEAR SOLUTION', &
           ' PACKAGE, VERSION 8, 04/28/2017')
    !
    ! -- DETERMINE PRECONDITIONER
    IF (this%LEVEL > 0 .OR. this%DROPTOL > DZERO) THEN
      this%IPC = 3
    ELSE
      this%IPC = 1
    END IF
    IF (this%RELAX > DZERO) THEN
      this%IPC = this%IPC + 1
    END IF
    !
    ! -- ERROR CHECKING FOR OPTIONS
    IF (this%ISCL < 0) this%ISCL = 0
    IF (this%ISCL > 2) THEN
      WRITE (errmsg, '(A)') 'IMSLINEAR7AR ISCL MUST BE <= 2'
      call store_error(errmsg)
    END IF
    IF (this%IORD < 0) this%IORD = 0
    IF (this%IORD > 2) THEN
      WRITE (errmsg, '(A)') 'IMSLINEAR7AR IORD MUST BE <= 2'
      call store_error(errmsg)
    END IF
    IF (this%NORTH < 0) THEN
      WRITE (errmsg, '(A)') 'IMSLINEAR7AR NORTH MUST >= 0'
      call store_error(errmsg)
    END IF
    IF (this%RCLOSE == DZERO) THEN
      IF (this%ICNVGOPT /= 3) THEN
        WRITE (errmsg, '(A)') 'IMSLINEAR7AR RCLOSE MUST > 0.0'
        call store_error(errmsg)
      END IF
    END IF
    IF (this%RELAX < DZERO) THEN
      WRITE (errmsg, '(A)') 'IMSLINEAR7AR RELAX MUST BE >= 0.0'
      call store_error(errmsg)
    END IF
    IF (this%RELAX > DONE) THEN
      WRITE (errmsg, '(A)') 'IMSLINEAR7AR RELAX MUST BE <= 1.0'
      call store_error(errmsg)
    END IF
    !
    ! -- INITIALIZE IMSLINEAR VARIABLES
    this%NITERC = 0
    !
    ! -- ALLOCATE AND INITIALIZE MEMORY FOR IMSLINEAR
    iscllen = 1
    IF (this%ISCL .NE. 0) iscllen = NEQ
    ALLOCATE(this%DSCALE(iscllen))
    ALLOCATE(this%DSCALE2(iscllen))
    !
    ! -- ALLOCATE MEMORY FOR PRECONDITIONING MATRIX
    ijlu = 1
    ijw = 1
    iwlu = 1
    !
    ! -- ILU0 AND MILU0
    this%NIAPC = this%NEQ
    this%NJAPC = this%NJA
    !
    ! -- ILUT AND MILUT
    IF (this%IPC == 3 .OR. this%IPC == 4) THEN
      this%NIAPC = this%NEQ
      IF (this%LEVEL > 0) THEN
        iwk = this%NEQ * (this%LEVEL * 2 + 1)
      ELSE
        iwk = 0
        DO n = 1, NEQ
          i = this%IA(n + 1) - this%IA(n)
          IF (i > iwk) THEN
            iwk = i
          END IF
        END DO
        iwk = this%NEQ * iwk
      END IF
      this%NJAPC = iwk
      ijlu = iwk
      ijw = 2 * this%NEQ
      iwlu = this%NEQ + 1
    END IF
    this%NJLU = ijlu
    this%NJW = ijw
    this%NWLU = iwlu
    !
    ! -- ALLOCATE BASE PRECONDITIONER VECTORS
    ALLOCATE(this%IAPC(this%NIAPC + 1))
    ALLOCATE(this%JAPC(this%NJAPC))
    ALLOCATE(this%APC(this%NJAPC))
    !
    ! -- ALLOCATE MEMORY FOR ILU0 AND MILU0 NON-ZERO ROW ENTRY VECTOR
    ALLOCATE(this%IW(this%NIAPC))
    ALLOCATE(this%W(this%NIAPC))
    !
    ! -- ALLOCATE MEMORY FOR ILUT VECTORS
    ALLOCATE(this%JLU(ijlu))
    ALLOCATE(this%JW(ijw))
    ALLOCATE(this%WLU(iwlu))
    !
    ! -- GENERATE IAPC AND JAPC FOR ILU0 AND MILU0
    IF (this%IPC == 1 .OR. this%IPC == 2) THEN
      CALL ims_base_pccrs(this%NEQ, this%NJA, this%IA, this%JA, &
                          this%IAPC, this%JAPC)
    END IF
    !
    ! -- ALLOCATE SPACE FOR PERMUTATION VECTOR
    i0 = 1
    iolen = 1
    IF (this%IORD .NE. 0) THEN
      i0 = this%NEQ
      iolen = this%NJA
    END IF
    ALLOCATE(this%LORDER(i0))
    ALLOCATE(this%IORDER(i0))
    ALLOCATE(this%IARO(i0 + 1))
    ALLOCATE(this%JARO(iolen))
    ALLOCATE(this%ARO(iolen))
    !
    ! -- ALLOCATE WORKING VECTORS FOR IMSLINEAR SOLVER
    ALLOCATE(this%ID(this%NEQ))
    ALLOCATE(this%D(this%NEQ))
    ALLOCATE(this%P(this%NEQ))
    ALLOCATE(this%Q(this%NEQ))
    ALLOCATE(this%Z(this%NEQ))
    !
    ! -- ALLOCATE MEMORY FOR BCGS WORKING ARRAYS
    this%NIABCGS = 1
    IF (this%ILINMETH == 2) THEN
      this%NIABCGS = this%NEQ
    END IF
    ALLOCATE(this%T(this%NIABCGS))
    ALLOCATE(this%V(this%NIABCGS))
    ALLOCATE(this%DHAT(this%NIABCGS))
    ALLOCATE(this%PHAT(this%NIABCGS))
    ALLOCATE(this%QHAT(this%NIABCGS))
    !
    ! -- INITIALIZE IMSLINEAR VECTORS
    DO n = 1, iscllen
      this%DSCALE(n) = DONE
      this%DSCALE2(n) = DONE
    END DO
    DO n = 1, this%NJAPC
      this%APC(n) = DZERO
    END DO
    !
    ! -- WORKING VECTORS
    DO n = 1, this%NEQ
      this%ID(n) = IZERO
      this%D(n) = DZERO
      this%P(n) = DZERO
      this%Q(n) = DZERO
      this%Z(n) = DZERO
    END DO
    DO n = 1, this%NIAPC
      this%IW(n) = IZERO
      this%W(n) = DZERO
    END DO
    !
    ! -- BCGS WORKING VECTORS
    DO n = 1, this%NIABCGS
      this%T(n) = DZERO
      this%V(n) = DZERO
      this%DHAT(n) = DZERO
      this%PHAT(n) = DZERO
      this%QHAT(n) = DZERO
    END DO
    !
    ! -- ILUT AND MILUT WORKING VECTORS
    DO n = 1, ijlu
      this%JLU(n) = DZERO
    END DO
    DO n = 1, ijw
      this%JW(n) = DZERO
    END DO
    DO n = 1, iwlu
      this%WLU(n) = DZERO
    END DO
    !
    ! -- REORDERING VECTORS
    DO n = 1, i0 + 1
      this%IARO(n) = IZERO
    END DO
    DO n = 1, iolen
      this%JARO(n) = IZERO
      this%ARO(n) = DZERO
    END DO
    !
    ! -- REVERSE CUTHILL MCKEE AND MINIMUM DEGREE ORDERING
    IF (this%IORD .NE. 0) THEN
      CALL ims_base_calc_order(this%IORD, this%NEQ, this%NJA, this%IA, &
                               this%JA, this%LORDER, this%IORDER)
    END IF
    !
    ! -- ALLOCATE MEMORY FOR STORING ITERATION CONVERGENCE DATA
    !
    ! -- RETURN
    RETURN
  END SUBROUTINE imslinear_ar

  !> @ brief Write summary of settings
    !!
    !!  Write summary of linear accelerator settings.
    !!
  !<
  subroutine imslinear_summary(this, mxiter)
    ! -- dummy variables
    class(ImsLinearDataType), intent(inout) :: this !< ImsLinearDataType instance
    integer(I4B), intent(in) :: mxiter !< maximum number of outer iterations
    ! -- local variables
    CHARACTER(LEN=10) :: clin(0:2)
    CHARACTER(LEN=31) :: clintit(0:2)
    CHARACTER(LEN=20) :: cipc(0:4)
    CHARACTER(LEN=20) :: cscale(0:2)
    CHARACTER(LEN=25) :: corder(0:2)
    CHARACTER(LEN=16), DIMENSION(0:4) :: ccnvgopt
    CHARACTER(LEN=15) :: clevel
    CHARACTER(LEN=15) :: cdroptol
    integer(I4B) :: i
    integer(I4B) :: j
    ! -- data
    DATA clin/'UNKNOWN   ', &
             &'CG        ', &
             &'BCGS      '/
    DATA clintit/'             UNKNOWN           ', &
                &'       CONJUGATE-GRADIENT      ', &
                &'BICONJUGATE-GRADIENT STABILIZED'/
    DATA cipc/'UNKNOWN             ', &
             &'INCOMPLETE LU       ', &
             &'MOD. INCOMPLETE LU  ', &
             &'INCOMPLETE LUT      ', &
             &'MOD. INCOMPLETE LUT '/
    DATA cscale/'NO SCALING          ', &
               &'SYMMETRIC SCALING   ', &
               &'L2 NORM SCALING     '/
    DATA corder/'ORIGINAL ORDERING        ', &
               &'RCM ORDERING             ', &
               &'MINIMUM DEGREE ORDERING  '/
    DATA ccnvgopt/'INFINITY NORM   ', &
                 &'INFINITY NORM S ', &
                 &'L2 NORM         ', &
                 &'RELATIVE L2NORM ', &
                 &'L2 NORM W. REL. '/
    ! -- formats
02010 FORMAT(1X, /, 7X, 'SOLUTION BY THE', 1X, A31, 1X, 'METHOD', &
           /, 1X, 66('-'), /, &
           ' MAXIMUM OF ', I0, ' CALLS OF SOLUTION ROUTINE', /, &
           ' MAXIMUM OF ', I0, &
           ' INTERNAL ITERATIONS PER CALL TO SOLUTION ROUTINE', /, &
           ' LINEAR ACCELERATION METHOD            =', 1X, A, /, &
           ' MATRIX PRECONDITIONING TYPE           =', 1X, A, /, &
           ' MATRIX SCALING APPROACH               =', 1X, A, /, &
           ' MATRIX REORDERING APPROACH            =', 1X, A, /, &
           ' NUMBER OF ORTHOGONALIZATIONS          =', 1X, I0, /, &
           ' HEAD CHANGE CRITERION FOR CLOSURE     =', E15.5, /, &
           ' RESIDUAL CHANGE CRITERION FOR CLOSURE =', E15.5, /, &
           ' RESIDUAL CONVERGENCE OPTION           =', 1X, I0, /, &
           ' RESIDUAL CONVERGENCE NORM             =', 1X, A, /, &
           ' RELAXATION FACTOR                     =', E15.5)
02015 FORMAT(' NUMBER OF LEVELS                      =', A15, /, &
           ' DROP TOLERANCE                        =', A15, //)
2030 FORMAT(1X, A20, 1X, 6(I6, 1X))
2040 FORMAT(1X, 20('-'), 1X, 6(6('-'), 1X))
2050 FORMAT(1X, 62('-'),/) !
! -- -----------------------------------------------------------
    !
    ! -- initialize clevel and cdroptol
    clevel = ''
    cdroptol = ''
    !
    ! -- write common variables to all linear accelerators
    write (this%iout, 2010) &
      clintit(this%ILINMETH), MXITER, this%ITER1, &
      clin(this%ILINMETH), cipc(this%IPC), &
      cscale(this%ISCL), corder(this%IORD), &
      this%NORTH, this%DVCLOSE, this%RCLOSE, &
      this%ICNVGOPT, ccnvgopt(this%ICNVGOPT), &
      this%RELAX
    if (this%level > 0) then
      write (clevel, '(i15)') this%level
    end if
    if (this%droptol > DZERO) then
      write (cdroptol, '(e15.5)') this%droptol
    end if
    IF (this%level > 0 .or. this%droptol > DZERO) THEN
      write (this%iout, 2015) trim(adjustl(clevel)), &
        trim(adjustl(cdroptol))
    ELSE
      write (this%iout, '(//)')
    END IF

    if (this%iord /= 0) then
      !
      ! -- WRITE SUMMARY OF REORDERING INFORMATION TO LIST FILE
      if (this%iprims == 2) then
        DO i = 1, this%neq, 6
          write (this%iout, 2030) 'ORIGINAL NODE      :', &
            (j, j=i, MIN(i + 5, this%neq))
          write (this%iout, 2040)
          write (this%iout, 2030) 'REORDERED INDEX    :', &
            (this%lorder(j), j=i, MIN(i + 5, this%neq))
          write (this%iout, 2030) 'REORDERED NODE     :', &
            (this%iorder(j), j=i, MIN(i + 5, this%neq))
          write (this%iout, 2050)
        END DO
      END IF
    end if
    !
    ! -- return
    return
  end subroutine imslinear_summary

  !> @ brief Allocate and initialize scalars
    !!
    !!  Allocate and inititialize linear accelerator scalars
    !!
  !<
  subroutine allocate_scalars(this)
    ! -- dummy variables
    class(ImsLinearDataType), intent(inout) :: this !< ImsLinearDataType instance
    !
    ! -- PHASE-1 ISOLATION: scalars are plain type components (no memory manager);
    !    only initialization remains.
    ! -- initialize scalars
    this%iout = 0
    this%ipc = 0
    this%iacpc = 0
    this%niterc = 0
    this%niabcgs = 0
    this%niapc = 0
    this%njapc = 0
    this%epfact = DZERO
    this%l2norm0 = 0
    this%njlu = 0
    this%njw = 0
    this%nwlu = 0
    !
    ! -- return
    return
  end subroutine allocate_scalars

  !> @ brief Deallocate memory
    !!
    !!  Deallocate linear accelerator memory.
    !!
  !<
  subroutine imslinear_da(this)
    ! -- dummy variables
    class(ImsLinearDataType), intent(inout) :: this !< linear datatype instance
    !
    ! -- PHASE-1 ISOLATION: owned ALLOCATABLE arrays
    if (allocated(this%dscale)) deallocate (this%dscale)
    if (allocated(this%dscale2)) deallocate (this%dscale2)
    if (allocated(this%iapc)) deallocate (this%iapc)
    if (allocated(this%japc)) deallocate (this%japc)
    if (allocated(this%apc)) deallocate (this%apc)
    if (allocated(this%iw)) deallocate (this%iw)
    if (allocated(this%w)) deallocate (this%w)
    if (allocated(this%jlu)) deallocate (this%jlu)
    if (allocated(this%jw)) deallocate (this%jw)
    if (allocated(this%wlu)) deallocate (this%wlu)
    if (allocated(this%lorder)) deallocate (this%lorder)
    if (allocated(this%iorder)) deallocate (this%iorder)
    if (allocated(this%id)) deallocate (this%id)
    if (allocated(this%d)) deallocate (this%d)
    if (allocated(this%p)) deallocate (this%p)
    if (allocated(this%q)) deallocate (this%q)
    if (allocated(this%z)) deallocate (this%z)
    if (allocated(this%t)) deallocate (this%t)
    if (allocated(this%v)) deallocate (this%v)
    if (allocated(this%dhat)) deallocate (this%dhat)
    if (allocated(this%phat)) deallocate (this%phat)
    if (allocated(this%qhat)) deallocate (this%qhat)
    !
    ! -- owned POINTER arrays (reordered system; aliased by IA0/JA0/A0)
    if (associated(this%iaro)) deallocate (this%iaro)
    if (associated(this%jaro)) deallocate (this%jaro)
    if (associated(this%aro)) deallocate (this%aro)
    !
    ! -- nullify alias pointers (caller-owned storage; not deallocated here)
    nullify (this%ia)
    nullify (this%ja)
    nullify (this%amat)
    nullify (this%rhs)
    nullify (this%x)
    nullify (this%ia0)
    nullify (this%ja0)
    nullify (this%a0)
    !
    ! -- return
    return
  end subroutine imslinear_da

  !> @ brief Set default settings
    !!
    !!  Set default linear accelerator settings.
    !!
  !<
  SUBROUTINE imslinear_set_input(this, IFDPARAM)
    ! -- dummy variables
    CLASS(ImsLinearDataType), INTENT(INOUT) :: this !< ImsLinearDataType instance
    integer(I4B), INTENT(IN) :: IFDPARAM !< complexity option
    ! -- code
    SELECT CASE (IFDPARAM)
      !
      ! -- Simple option
    CASE (1)
      this%ITER1 = 50
      this%ILINMETH = 1
      this%IPC = 1
      this%ISCL = 0
      this%IORD = 0
      this%DVCLOSE = DEM3
      this%RCLOSE = DEM1
      this%RELAX = DZERO
      this%LEVEL = 0
      this%DROPTOL = DZERO
      this%NORTH = 0
      !
      ! -- Moderate
    CASE (2)
      this%ITER1 = 100
      this%ILINMETH = 2
      this%IPC = 2
      this%ISCL = 0
      this%IORD = 0
      this%DVCLOSE = DEM2
      this%RCLOSE = DEM1
      this%RELAX = 0.97D0
      this%LEVEL = 0
      this%DROPTOL = DZERO
      this%NORTH = 0
      !
      ! -- Complex
    CASE (3)
      this%ITER1 = 500
      this%ILINMETH = 2
      this%IPC = 3
      this%ISCL = 0
      this%IORD = 0
      this%DVCLOSE = DEM1
      this%RCLOSE = DEM1
      this%RELAX = DZERO
      this%LEVEL = 5
      this%DROPTOL = DEM4
      this%NORTH = 2
    END SELECT
    !
    ! -- return
    RETURN
  END SUBROUTINE imslinear_set_input

  !> @ brief Base linear accelerator subroutine
    !!
    !!  Base linear accelerator subroutine that scales and reorders
    !!  the system of equations, if necessary, updates the preconditioner,
    !!  and calls the appropriate linear accelerator.
    !!
  !<
  SUBROUTINE imslinear_ap(this, ICNVG, KSTP, KITER, IN_ITER, &
                          NCONV, CONVNMOD, CONVMODSTART, &
                          CACCEL, summary)
    ! -- modules
    USE SimModule
    ! -- dummy variables
    CLASS(ImsLinearDataType), INTENT(INOUT) :: this !< ImsLinearDataType instance
    integer(I4B), INTENT(INOUT) :: ICNVG !< convergence flag (1) non-convergence (0)
    integer(I4B), INTENT(IN) :: KSTP !< time step number
    integer(I4B), INTENT(IN) :: KITER !< outer iteration number
    integer(I4B), INTENT(INOUT) :: IN_ITER !< inner iteration number
    ! -- convergence information dummy variables
    integer(I4B), INTENT(IN) :: NCONV !<
    integer(I4B), INTENT(IN) :: CONVNMOD !<
    integer(I4B), DIMENSION(CONVNMOD + 1), INTENT(INOUT) :: CONVMODSTART !<
    character(len=31), DIMENSION(NCONV), INTENT(INOUT) :: CACCEL !<
    type(ConvergenceSummaryType), pointer, intent(in) :: summary !< Convergence summary report
    ! -- local variables
    integer(I4B) :: n
    integer(I4B) :: innerit
    integer(I4B) :: irc
    integer(I4B) :: itmax
    real(DP) :: dnrm2
    !
    ! -- set epfact based on timestep
    IF (this%ICNVGOPT == 2) THEN
      IF (KSTP == 1) THEN
        this%EPFACT = 0.01
      ELSE
        this%EPFACT = 0.10
      END IF
    ELSE IF (this%ICNVGOPT == 4) THEN
      this%EPFACT = DEM4
    ELSE
      this%EPFACT = DONE
    END IF
    !
    ! -- SCALE PROBLEM
    IF (this%ISCL .NE. 0) THEN
      CALL ims_base_scale(0, this%ISCL, &
                          this%NEQ, this%NJA, this%IA, this%JA, &
                          this%AMAT, this%X, this%RHS, &
                          this%DSCALE, this%DSCALE2)
    END IF
    !
    ! -- PERMUTE ROWS, COLUMNS, AND RHS
    IF (this%IORD /= 0) THEN
      CALL ims_dperm(this%NEQ, this%AMAT, this%JA, this%IA, &
                 this%ARO, this%JARO, this%IARO, &
                 this%LORDER, this%ID, 1)
      CALL dvperm(this%NEQ, this%X, this%LORDER)
      CALL dvperm(this%NEQ, this%RHS, this%LORDER)
      this%IA0 => this%IARO
      this%JA0 => this%JARO
      this%A0 => this%ARO
    ELSE
      this%IA0 => this%IA
      this%JA0 => this%JA
      this%A0 => this%AMAT
    END IF
    !
    ! -- UPDATE PRECONDITIONER
    CALL ims_base_pcu(this%iout, this%NJA, this%NEQ, this%NIAPC, this%NJAPC, &
                      this%IPC, this%RELAX, this%A0, this%IA0, this%JA0, &
                      this%APC, this%IAPC, this%JAPC, this%IW, this%W, &
                      this%LEVEL, this%DROPTOL, this%NJLU, this%NJW, &
                      this%NWLU, this%JLU, this%JW, this%WLU)
    !
    ! -- F-4: stop at the store site. A hard ILUT zero-pivot in ims_base_pcu calls
    !    store_error; without this guard the accelerator would run on a partially
    !    built JLU/APC factorization before the wrapper's post-solve
    !    count_errors()>0 guard fires. Return immediately so IMS7U1AP raises IERR=1
    !    (fail-clean) right away. Inert on a healthy solve (count_errors()==0).
    if (count_errors() > 0) then
      ICNVG = 0
      IN_ITER = 0
      return
    end if
    !
    ! -- INITIALIZE SOLUTION VARIABLE AND ARRAYS
    IF (KITER == 1) then
      this%NITERC = 0
      summary%iter_cnt = 0
    end if
    irc = 1
    ICNVG = 0
    DO n = 1, this%NEQ
      this%D(n) = DZERO
      this%P(n) = DZERO
      this%Q(n) = DZERO
      this%Z(n) = DZERO
    END DO
    !
    ! -- CALCULATE INITIAL RESIDUAL
    call ims_base_residual(this%NEQ, this%NJA, this%X, this%RHS, this%D, &
                           this%A0, this%IA0, this%JA0)
    this%L2NORM0 = dnrm2(this%NEQ, this%D, 1)
    !
    ! -- CHECK FOR EXACT SOLUTION
    itmax = this%ITER1
    IF (this%L2NORM0 == DZERO) THEN
      itmax = 0
      ICNVG = 1
    END IF
    !
    ! -- SOLUTION BY THE CONJUGATE GRADIENT METHOD
    IF (this%ILINMETH == 1) THEN
      CALL ims_base_cg(ICNVG, itmax, innerit, &
                       this%NEQ, this%NJA, this%NIAPC, this%NJAPC, &
                       this%IPC, this%ICNVGOPT, this%NORTH, &
                       this%DVCLOSE, this%RCLOSE, this%L2NORM0, &
                       this%EPFACT, this%IA0, this%JA0, this%A0, &
                       this%IAPC, this%JAPC, this%APC, &
                       this%X, this%RHS, this%D, this%P, this%Q, this%Z, &
                       this%NJLU, this%IW, this%JLU, &
                       NCONV, CONVNMOD, CONVMODSTART, &
                       CACCEL, summary)
      !
      ! -- SOLUTION BY THE BICONJUGATE GRADIENT STABILIZED METHOD
    ELSE IF (this%ILINMETH == 2) THEN
      CALL ims_base_bcgs(ICNVG, itmax, innerit, &
                         this%NEQ, this%NJA, this%NIAPC, this%NJAPC, &
                         this%IPC, this%ICNVGOPT, this%NORTH, &
                         this%ISCL, this%DSCALE, &
                         this%DVCLOSE, this%RCLOSE, this%L2NORM0, &
                         this%EPFACT, this%IA0, this%JA0, this%A0, &
                         this%IAPC, this%JAPC, this%APC, &
                         this%X, this%RHS, this%D, this%P, this%Q, &
                         this%T, this%V, this%DHAT, this%PHAT, this%QHAT, &
                         this%NJLU, this%IW, this%JLU, &
                         NCONV, CONVNMOD, CONVMODSTART, &
                         CACCEL, summary)
    END IF
    !
    ! -- BACK PERMUTE AMAT, SOLUTION, AND RHS
    IF (this%IORD /= 0) THEN
      CALL ims_dperm(this%NEQ, this%A0, this%JA0, this%IA0, &
                 this%AMAT, this%JA, this%IA, &
                 this%IORDER, this%ID, 1)
      CALL dvperm(this%NEQ, this%X, this%IORDER)
      CALL dvperm(this%NEQ, this%RHS, this%IORDER)
    END IF
    !
    ! -- UNSCALE PROBLEM
    IF (this%ISCL .NE. 0) THEN
      CALL ims_base_scale(1, this%ISCL, &
                          this%NEQ, this%NJA, this%IA, this%JA, &
                          this%AMAT, this%X, this%RHS, &
                          this%DSCALE, this%DSCALE2)
    END IF
    !
    ! -- SET IMS INNER ITERATION NUMBER (IN_ITER) TO NUMBER OF
    !       IMSLINEAR INNER ITERATIONS (innerit)
    IN_ITER = innerit
    !
    ! -- RETURN
    RETURN
  END SUBROUTINE imslinear_ap

END MODULE IMSLinearModule

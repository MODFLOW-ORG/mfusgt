
c             for XMD package (version 1.3)
c
c                     M. Ibaraki
c
c             Thu Oct 29 10:03:25 EDT 2009
c
c     variable definitions:
c
c      iacl               choice of acceleration method
c                         = 0; conjugate gradient
c                         = 1; ORTHOMIN
c                         = 2; CGSTAB
c      n                  number of unknowns
c      norder             = 0; original ordering
c                         = 1; RCM ordering
c                         = 2; Minimum Degree ordering
c      nja                size of ja, a, arrays
c      njaf               size of af, jaf arrays
c      level              level of ILU
c      itmax              number of maximum allowable iterations
c      north              number of orthogonalization for the ORTHOMIN
c      liwrk              size of integer work array
c      lrwrk              size of real work array
c
c      ia(n+1),ja(nja)    usual ia, ja arrays for coefficient matrix
c      lorder(n)          ordering vector: lorder( new_order ) = old_order
c      iwork(liwrk)      temporary work array
c
c      dptol              flag for the drop tolerance
c                         =.true. perform the drop tolerance
c                         =.false. do NOT perform the drop tolerance
c
c      epsrn              drop tolerance
c      ctol               absolute convergence criteria
c      rrctol             residual reduction convergence criteria
c
c      rwork(lrwrk)       temporary work array
c      a(nja)             matrix stored as linear array
c      af(njaf)           factored matrix (each row of af contains a row L\U)
c                         where A = LU
c      b(n)               right hand side vector
c      x(n)               solution vector
c
c
c      nx,ny,nz           graph of matrix is regular rectangular grid
c                         of size nx * ny * nz
c
cmi
cmi   MODULE XMDMODULE
cmi   IMPLICIT NONE
cmi   LOGICAL, SAVE, POINTER ::  REDSYS,LDCOMB
cmi   DOUBLE PRECISION, SAVE, POINTER ::  EPSRN,RRCTOL
cmi   DOUBLE PRECISION, SAVE, DIMENSION(:),ALLOCATABLE::RWORK,AF,DGSCAL
cmi   INTEGER, SAVE, DIMENSION(:), ALLOCATABLE ::  LORDER,IWORK,MSINDX
cmi   INTEGER, SAVE, POINTER :: IACL,NORDER,NJAF,LEVEL,NORTH,LIWRK,
cmi  *  LRWRK,IDROPTOL,ICOMB,NBLACK,IERR,IDSCALE,MBLACK
cmi   END MODULE XMDMODULE
cmi

      MODULE XMDMODULE
      IMPLICIT NONE
      LOGICAL, SAVE, POINTER ::  REDSYS
      LOGICAL, SAVE, POINTER :: ILUREUSE
      DOUBLE PRECISION, SAVE, POINTER ::  EPSRN,RRCTOL
      DOUBLE PRECISION, SAVE, DIMENSION(:),ALLOCATABLE::DGSCAL
      INTEGER, SAVE, DIMENSION(:), ALLOCATABLE ::  LORDER
      INTEGER, SAVE, POINTER :: IACL,NORDER,LEVEL,NORTH,
     *  IDROPTOL,IERR,IDSCALE
      END MODULE XMDMODULE




C------------------------------------------------------------------
      SUBROUTINE XMD7U1AR(IN,IFDPARAM)

      USE GLOBAL, ONLY: NEQS,IOUT,STRT,IBOUND,AMAT,RHS,HNEW,NJA,IA,JA
      USE XMDMODULE
      USE SMSMODULE, ONLY: IPRSMS,ISOLVEACTIVE
cmi
      use xmdcmn
cmi
      IMPLICIT NONE
!     ------------------------------------------------------------------
!     SPECIFICATIONS:
!     ------------------------------------------------------------------
      INTRINSIC INT
      EXTERNAL URDCOM, URWORD
cmi
cmi   include 'xmdcmn.com'
cmi
!     ------------------------------------------------------------------
!     ARGUMENTS
!     ------------------------------------------------------------------
      INTEGER IN, IFDPARAM
!     ------------------------------------------------------------------
!     LOCAL VARIABLES
!     ------------------------------------------------------------------
      INTEGER lloc, istart, istop, i, n,IREDSYS
      CHARACTER(LEN=200) line
      REAL R,RRCTOLS,EPSRNS
!     LOCAL VARIABLES FOR XMD SOLVER
!     ------------------------------------------------------------------
!
!1------IDENTIFY PACKAGE AND INITIALIZE.
      WRITE(IOUT,1)IN
    1 FORMAT(1X,'XMD -- LINEAR SOLUTION BY XMD PACKAGE VERSION',
     &       1X,'1.30',
     & /5X,'BY MOTOMU IBARAKI, THE OHIO STATE UNIVERSITY, COLUMBUS, OH',
     & /1X,'                INPUT READ FROM UNIT',I3)
C
cmi
c     ALLOCATE (IACL,NORDER,NJAF,LEVEL,NORTH,LIWRK,LRWRK,IDROPTOL,
c    *  ICOMB,NBLACK,IERR,IDSCALE,MBLACK)
c     ALLOCATE (EPSRN,RRCTOL)
c     ALLOCATE (REDSYS,LDCOMB)
cmi
      ALLOCATE (IACL,NORDER,LEVEL,NORTH,IDROPTOL,IERR,IDSCALE)
      ALLOCATE (EPSRN,RRCTOL)
      ALLOCATE (REDSYS)
      ALLOCATE (ILUREUSE)
cmi

!-----XMD INPUT
      IF ( IFDPARAM.EQ.0 ) THEN
      CALL URDCOM(In, Iout, line)
      lloc = 1
      i = 1
      CALL URWORD(line, lloc, istart, istop, 2, IACL, r, Iout, In)
      CALL URWORD(line, lloc, istart, istop, 2, NORDER, r, Iout, In)
      CALL URWORD(line, lloc, istart, istop, 2, LEVEL, r, Iout, In)
      CALL URWORD(line, lloc, istart, istop, 2, NORTH, r, Iout, In)
csp      CALL URWORD(line, lloc, istart, istop, 2, IDSCALE, r, Iout, In)
      CALL URWORD(line, lloc, istart, istop, 2, IREDSYS, r, Iout, In)
      CALL URWORD(line, lloc, istart, istop, 3, I, RRCTOLS, Iout, In)
      CALL URWORD(line, lloc, istart, istop, 2, IDROPTOL, r, Iout, In)
      CALL URWORD(line, lloc, istart, istop, 3, I, EPSRNS, Iout, In)
C
      RRCTOL = RRCTOLS
      EPSRN = EPSRNS
      MIUNIT = IOUT
      MIOUT = IPRSMS - 1
      if(north.eq.0)north = 7
      ELSE
        CALL SET_XMDINPUT(IFDPARAM,IREDSYS)
      END IF
C-------WRITE WARNING IF REDUCED SYSTEM IS USED BY USER
      IF(IREDSYS.EQ.1)THEN
        WRITE(IOUT,15)
15      FORMAT(/5X,'*** WARNING: REDUCED SYSTEM (IREDSYS=1) SHOULD NOT',
     1    1X,'BE USED WITH MODFLOW-USG. ***'/5X,'IREDSYS=1  REQUIRES'
     2    ,1X,'ASCENDING ORDER OF JA COLUMNS AND GNC CONNECTIONS')
      ENDIF
C
      REDSYS = .FALSE.
      idscale = 0
      IF(IREDSYS.EQ.1) THEN
        REDSYS = .TRUE.
        IDSCALE = 1   ! NEED DIAGONAL SCALING FOR REDUCED SYSTEMS
      ENDIF
C
      WRITE(IOUT,23) IACL,NORDER,LEVEL,NORTH,IREDSYS,RRCTOL,
     *  IDROPTOL,EPSRN
   23 FORMAT(1X,'ACCELERATION METHOD                    (IACL) = ',I9/
     *      1X,'EQUATION ORDERING FLAG               (NORDER) = ',I9/
     *      1X,'LEVEL OF FILL                         (LEVEL) = ',I9/
     *      1X,'MAXIMUM NUMBER OF ORTHOGONALIZATIONS  (NORTH) = ',I9/
csp     *      1X,'INDEX FOR DIAGONAL SCALING          (IDSCALE) = ',I9/
     *      1X,'INDEX FOR USING REDUCED SYSTEM      (IREDSYS) = ',I9/
     *      1X,'RESID. REDUCTION CONVERGE CRITERION  (RRCTOL) = ',E13.6/
     *      1X,'INDEX FOR USING DROP TOLERANCE     (IDROPTOL) = ',I9/
     *      1X,'DROP TOLERANCE VALUE                 ((EPSRN) = ',E13.6)
!
!4-----ALLOCATE SPACE USED BY SOLVER
cmi
c      MBLACK = NEQS
c      IF(IREDSYS.EQ.1)MBLACK = NEQS * 0.5 + 1
c      NJAF = ISTORXMD * NJA
c      IF(LDCOMB)THEN
c        liwrk = 3*NEQS + 4*mblack + 2*njaf + 1  !        = .true.
c      ELSE
c        liwrk = 3*NEQS + 3*mblack + njaf + 1    ! ldcomb = .false.
c      ENDIF
c      lrwrk = ISTORXMD/2*mblack + 2*(north+1)*mblack + north
cmi



C
cmi   ALLOCATE(AF(NJAF),RWORK(LRWRK),DGSCAL(NEQS))
cmi   ALLOCATE(LORDER(NEQS),IWORK(LIWRK),MSINDX(30))
cmi   AF = 0.0
cmi   RWORK = 0.0
cmi   DGSCAL = 0.0
cmi   LORDER = 0
cmi   IWORK = 0
cmi   MSINDX = 0
C
      IF(ISOLVEACTIVE.EQ.0) CALL XMD_PREPROC()
C
      RETURN
      END
C-----------------------------------------------------------------------------------
      SUBROUTINE XMD_PREPROC()
C     ******************************************************************
C     CALL SUBROUTINES FOR SETTING XMD ARRAYS
C     ******************************************************************
      USE GLOBAL,   ONLY:IA,JA,NEQS,NJA,IA2,JA2,NEQS2,NJA2,IA1IA2,
     1  JA1JA2,IBOUND
      USE SMSMODULE,ONLY:LINMETH
      USE XMDMODULE,ONLY:REDSYS,IDROPTOL,level,DGSCAL,NORDER,ILUREUSE
C
C-----CALL SOLVER SUBROUTINES TO RESET INTERNAL ARRAYS
      IF(ALLOCATED(DGSCAL)) DEALLOCATE(DGSCAL)
      ALLOCATE(DGSCAL(NEQS2))
      DGSCAL = 0.0
      call xmdprpc(ia2, ja2, nja2, neqs2, norder, ierr, redsys)
      call xmdcheck(ia2, ja2, neqs2, nja2, ierr)
      ILUREUSE = .FALSE.
c  ---------------------------------------------------
c     PERFORM SYMBOLIC FACTORIZATION FOR LEVEL BASED PRECONDITIONING
c  ---------------------------------------------------
      IF(IDROPTOL.EQ.0)THEN
c  --------------------------------
c     level based preconditioning
c  --------------------------------
        call xmdprecl(ia2, ja2, level, nja2, neqs2, ierr)
      ENDIF
C
C-----RETURN
      RETURN
      END
C -----------------------------------------------------------------------
C
      SUBROUTINE SET_XMDINPUT(IFDPARAM,IREDSYS)
      USE XMDMODULE, ONLY:  IACL,NORDER,LEVEL,NORTH,IDROPTOL,
     +                      RRCTOL,EPSRN
      INTEGER IFDPARAM, IREDSYS
C Simple option
      SELECT CASE ( IFDPARAM )
      CASE(1)
        IACL = 1
        NORDER = 0
        LEVEL = 3
        NORTH = 5
        IREDSYS = 1
        IDROPTOL = 0
        RRCTOL = 0.0
        EPSRN = 1.0E-3
C Moderate
      CASE(2)
        IACL = 2
        NORDER = 0
        LEVEL = 3
        NORTH = 5
        IREDSYS = 1
        IDROPTOL = 1
        RRCTOL = 0.0
        EPSRN = 1.0E-3
C Complex
      CASE(3)
        IACL = 2
        NORDER = 1
        LEVEL = 5
        NORTH = 7
        IREDSYS = 1
        IDROPTOL = 1
        RRCTOL = 0.0
        EPSRN = 1.0E-5
      END SELECT
      RETURN
      END
      SUBROUTINE XMD7DA
      USE GLOBAL, ONLY: NEQS,IOUT,STRT,IBOUND,AMAT,RHS,HNEW,NJA,IA,JA
      USE XMDMODULE
      USE SMSMODULE, ONLY: IPRSMS
      use xmdcmn
      use xmdmatrix
      INTEGER ALLOC_ERR
      DEALLOCATE (IACL,NORDER,LEVEL,NORTH,IDROPTOL,IERR,IDSCALE
     & , STAT = ALLOC_ERR)
      DEALLOCATE (EPSRN,RRCTOL, STAT = ALLOC_ERR)
      DEALLOCATE (REDSYS, STAT = ALLOC_ERR)
      DEALLOCATE (ILUREUSE, STAT = ALLOC_ERR)
      DEALLOCATE(DGSCAL, STAT = ALLOC_ERR)
      RETURN
      END


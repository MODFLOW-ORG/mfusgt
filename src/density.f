C -------------------------------------------------------------------------------------
      MODULE DDFMODULE
C
        INTEGER, SAVE, POINTER :: IDDF, ITHICKAV, IMPHDD, ISHARP,ISWIUN 
        DOUBLE PRECISION, SAVE, POINTER :: RHOFRESH,RHOSTD,CSTD
        DOUBLE PRECISION, SAVE, POINTER :: DRHODC,RELRHODIFF,HSEA
        DOUBLE PRECISION, SAVE,DIMENSION(:),ALLOCATABLE::
     1         RHONORM, ZEECELL, THICKCELL, ZETASWI
C        DOUPLE PRECISION, SAVE, DIMENSION (:), ALLOCATABLE :: AMATDD
c -------amatdd to save matrix before implicit dd term; but no need since we recompute (good to debug)
C 
      END MODULE DDFMODULE
C
C -------------------------------------------------------------------------------------
C
      SUBROUTINE DDF1AR(IUDDF)
C     ******************************************************************
C     ALLOCATE SPACE AND READ INFORMATION FOR DENSITY DEPENDENT FLOW
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE DDFMODULE
      USE GLOBAL,ONLY:IUNIT,IOUT,NEQS,NODES,IFREFM,IUNSTR,INDDF,HNEW,NJA
      USE GWTBCTMODULE, ONLY: CONC
      DOUBLE PRECISION ZEE,RNORM,THICK
      REAL TEMPVAR
      CHARACTER*400 LINE
C
C     ------------------------------------------------------------------
C
C1------IDENTIFY PACKAGE.
        INDDF = IUDDF
        WRITE(IOUT,1)INDDF
    1   FORMAT(1X,/1X,'DDF -- DENSITY DEPENDENT FLOW MODULE ',
     1    'VERSION 1, 4/4/2016 INPUT READ FROM UNIT ',I4)
C
C2------ALLOCATE SCALAR VARIABLES VECTORS, AND INITIALIZE.
      ALLOCATE(IDDF, ITHICKAV, IMPHDD, ISHARP)
        IDDF = 1  ! INDEX THAT DENSITY DEPENDENT FLOW IS ACTIVE ON ONE SPECIES (HARDWIRE)
      ALLOCATE(RHOFRESH,RHOSTD,CSTD)
      ALLOCATE (DRHODC)
C
C3------READ FRESHWATER DENSITY, STANDARD SOLUTION DENSITY,  STANDARD CONCENTRATION, AND OPTIONS
      CALL URDCOM(INDDF,IOUT,LINE)
      IF(IFREFM.EQ.0) THEN
        READ(LINE,'(5F10.4)')RHOFRESH,RHOSTD,CSTD,ITHICKAV,IMPHDD,ISHARP
        LLOC=61
      ELSE
        LLOC=1
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,I,TEMPVAR,IOUT,INDDF)
        RHOFRESH = TEMPVAR
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,I,TEMPVAR,IOUT,INDDF)
        RHOSTD = TEMPVAR
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,I,TEMPVAR,IOUT,INDDF)
        CSTD = TEMPVAR 
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,ITHICKAV,R,IOUT,INDDF)
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IMPHDD,R,IOUT,INDDF)
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,ISHARP,R,IOUT,INDDF)
      END IF      
C---------------------------------------------------------------------------
C3A-----REFLECT INPUT IN OUTPUT LISTING FILE
      WRITE(IOUT,3) RHOFRESH,RHOSTD,CSTD,ITHICKAV,IMPHDD,ISHARP
    3 FORMAT(1X,'DENSITY OF FRESHWATER (RHOFRESH) =',F15.3,
     1  /1X,'DENSITY OF STANDARD SOLUTION (RHOSTD) =', F15.3,
     1  /1X,'CONCENTRATION OF STANDARD SOLUTION (CSTD) =',F15.3
     1  /1X,'FLAG FOR USING THICKNESS WEIGHTED AVERAGE (ITHICKAV) =',I2
     1  /1X,'FLAG FOR IMPLICIT TREATMENT OF HEAD TERM (IMPHDD) =',I2 
     1  /1X,'FLAG FOR SHARP INTERFACE MODEL (ISHARP) =',I2 )
C---------------------------------------------------------------------------
C4------FOR SHARP INTERFACE OPTION, ALLOCATE AND READ  VARIABLES
      IF(ISHARP.NE.0) THEN 
        CALL SHARPINT1AR
      ELSE
C---------------------------------------------------------------------------
C4------FOR DENSITY DEPENDENT SIMULATION, ALLOCATE AND READ REMAINING VARIABLES            
        ALLOCATE(RHONORM(NEQS),ZEECELL(NEQS))            
        IF(ITHICKAV.NE.0) THEN
          ALLOCATE(THICKCELL(NEQS))  
        ENDIF    
c        IF(IMPHDD.NE.0) THEN 
c          ALLOCATE(AMATDD(NJA))  
c        ENDIF             
C ---------------------------------------------------------------------------
C4------FILL ARRAY  ZEECELL, INITIALIZE RHO/RHOFRESH, AND 
C4------INITIALIZE HNEW TO BE THE NORMALIZED POTENTIAL (RHO/RHOFRESH*HEAD)
        DRHODC = (RHOSTD - RHOFRESH)/CSTD 
        DO N=1,NEQS
          CALL RHONORMCALC(N,RNORM)  
          RHONORM(N) = RNORM  
C        
          CALL ZCELLCALC(N,ZEE,THICK)
          ZEECELL(N) = ZEE
          IF(ITHICKAV.NE.0) THEN 
            THICKCELL(N) = THICK
          ENDIF  
C        
        ENDDO
      ENDIF
C---------------------------------------------------------------------------
C5----RETURN
      RETURN
      END
C--------------------------------------------------------------------------- 
      SUBROUTINE DDF1AD(IUDDF)
C     ******************************************************************
C     ADVANCE DENSITY OF ALL CELLS FROM CONCENTRATION AFTER TRANSPORT 
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE DDFMODULE, ONLY: DRHODC,RHONORM,ISHARP
      USE GLOBAL,ONLY: IUNIT,IOUT,NEQS,NODES,IFREFM,IUNSTR,INDDF,HNEW,SN
      USE GWTBCTMODULE, ONLY: CONC
      DOUBLE PRECISION ZEE,RNORM
      CHARACTER*400 LINE
C
C     ------------------------------------------------------------------
      IF (ISHARP. NE. 0) RETURN 
C1----LOOP OVER ALL CELLS
      DO N=1,NEQS
C2------COMPUTE RHO/RHOFRESH
c        IF(CONC(N,1).GT.CSTD) CONC(N,1) = CSTD
c        IF(CONC(N,1).LT.0.0) CONC(N,1) = 0.0
        CALL RHONORMCALC(N,RNORM)  
        RHONORM(N) = RNORM  
      ENDDO
C---------------------------------------------------------------------------
C3----RETURN
      RETURN
      END
C---------------------------------------------------------------------------
      SUBROUTINE RHONORMCALC(N,RNORM)
C     ******************************************************************
C     CALCULATE THE NORMALIZED RHO FOR A GIVEN CONC AT CELL N
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE DDFMODULE
      USE GLOBAL,ONLY: IUNIT,IOUT,NEQS,NODES,IFREFM,IUNSTR,INDDF,HNEW,Sn
      USE GWTBCTMODULE, ONLY: CONC
      DOUBLE PRECISION RNORM,CNC
C     ------------------------------------------------------------------
C
C1----COMPUTE RNORM
      CNC = CONC(N,1)
      IF(CNC.LT.0.0) CNC = 0.0
      IF(CNC.GT.CSTD) CNC = CSTD
      RNORM = (RHOFRESH + Sn(N) * DRHODC * CNC)/RHOFRESH
C---------------------------------------------------------------------------
C2----RETURN
      RETURN
      END      
C---------------------------------------------------------------------------
      SUBROUTINE ZCELLCALC(N,ZEE,THICK)
C     ******************************************************************
C     CALCULATE CELL ELEVATION FOR DENSITY DEPENDENT FLOW
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE DDFMODULE
      USE GLOBAL, ONLY: IUNIT,IOUT,NEQS,NODES,IFREFM,IUNSTR,INDDF,HNEW,
     1            TOP,BOT
      USE CLN1MODULE, ONLY: ACLNNDS
      USE GWTBCTMODULE, ONLY: CONC
      DOUBLE PRECISION ZEE,FRAD,FLENG,FANGLE,FDPTH,THICK
C     ------------------------------------------------------------------
C
C1----CHECK IF CELL IS BCF OR CLN
      IF(N.LE.NODES)THEN
C2------GROUNDWATER CELL ELEVATION IS AVERAGE OF TOP AND BOT
        ZEE = 0.5 * (TOP(N) + BOT(N))
        THICK = TOP(N) - BOT(N)
      ELSE 
C3------CLN CELL ELEVATION SET AT BOTTOM OF THE CLN CELL (LOCAL CELL INDEX)
        I = N - NODES
        ZEE = ACLNNDS(I,5) 
C4------ADJUST CLN CELL ELEVATION TO CENTER DEPENDING ON ORITENTATION
        IFDIR = ACLNNDS(I,3) 
        IF(IFDIR.EQ.2)THEN 
C5--------ANGLED PIPE            
          FLENG = ACLNNDS(I,4) 
          FANGLE = ACLNNDS(I,6)   
          FDPTH = FLENG * SIN(FANGLE)
        ELSEIF(IFDIR.EQ.1)THEN
C6--------HORIZONTAL PIPE
          IFTYP =  ACLNNDS(I,2)
          CALL CLNR(IFTYP,FRAD)
C          FDPTH = 2.0 * FRAD
          FDPTH = FRAD
        ELSEIF(IFDIR.EQ.0)THEN  
C7--------VERTICAL PIPE
         FLENG = ACLNNDS(I,4) 
         FDPTH = FLENG
        ENDIF  
        ZEE = ZEE + 0.5 * FDPTH 
        THICK = FDPTH
      ENDIF
C---------------------------------------------------------------------------
C8----RETURN
      RETURN
      END      
C---------------------------------------------------------------------------
      SUBROUTINE DDF1FM(KSTP,KPER)
C     ******************************************************************
C     COMPUTE FLOW BETWEEN ADJACENT CELLS DUE TO DENSITY TERM
C     FOR TRANSIENT SIMULATION ALSO ADD DENSITY STORAGE TERM      
C     ******************************************************************
C
C     SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL, ONLY:IBOUND,NEQS,IOUT,NODES,NJA,IA,JA,JAS,IUNSTR,ISYM,
     1                AMAT,RHS,Sn,ISSFLG,AREA,TOP,BOT,HNEW,HOLD
      USE GWFBASMODULE,ONLY:DELT
      USE SMSMODULE, ONLY: AMATFL
      USE GWTBCTMODULE, ONLY: PRSITY,CONC,CONCO
      USE CLN1MODULE, ONLY: ACLNNDS
      USE DDFMODULE, ONLY: ZEECELL,RHONORM,DRHODC,RHOFRESH,CSTD,
     1  THICKCELL, ITHICKAV, IMPHDD, ISHARP                                     !  , AMATDD
C
      DOUBLE PRECISION QCON,ZEENJJ,RHO,VOL,HPHI,RHOTERM,AMAT_TERM,
     1  CNC,CNCO,GRHO,OMEGA
C     ------------------------------------------------------------------
C     
      IF (ISHARP. NE. 0) RETURN 
C0 -----SAVE AMAT INTO AMATDD FOR IMPLICIT UPDATE  (NO NEED, CAN BACK OUT THE TERMS)    
c      IF(IMPHDD.EQ.1) THEN 
c        DO IJA = 1,NJA 
c          AMATDD(IJA) = AMAT(IJA)  
c        ENDDO
c      ENDIF  
C1------FILL DENSITY GRADIENT TERM FOR EVERY CELL
      DO N=1,NEQS
C
C2------IF CELL IS NOT ACTIVE GO ON TO NEXT CELL.
        IF (IBOUND(N).LE.0) CYCLE
C
C3------CALCULATE FOR ALL CONNECTING FACES.
        DO II = IA(N)+1,IA(N+1)-1
          JJ = JA(II)
          IIS = JAS(II)
          IF(IBOUND(JJ).LE.0) CYCLE
C4--------COMPUTE AVERAGE Z AND H OF THE TWO CELLS 
          IF(ITHICKAV.EQ.0) THEN 
            ZEENJJ = (ZEECELL(N) + ZEECELL(JJ)) * 0.5
            HPHI = (HNEW(N) + HNEW(JJ)) * 0.5  
            OMEGA = 0.5
          ELSE
            OMEGA = THICKCELL(N) / (THICKCELL(N) + THICKCELL(JJ))  
            ZEENJJ = (1.0 - OMEGA) * ZEECELL(N) + OMEGA * ZEECELL(JJ) 
            HPHI = (1.0 - OMEGA) * HNEW(N) + OMEGA * HNEW(JJ) 
          ENDIF  
C5--------CALCULATE FLOW DUE TO DENSITY TERM         
          IF(IMPHDD.EQ.0) THEN 
            QCON = AMAT(II) *(HPHI - ZEENJJ) * (RHONORM(N)-RHONORM(JJ))
          ELSE 
             QCON = AMAT(II) *(- ZEENJJ) * (RHONORM(N)-RHONORM(JJ))  
          ENDIF  
C6--------ADD THIS FLOW TERM ON RHS OF N        
          RHS(N) = RHS(N) + QCON 
        ENDDO
C
      ENDDO
C ------------------------------------------------------------------------------------
C7------UPDATE FLOW MATRIX WITH DENSITY TERM 
      DO N=1,NEQS
C
C8------IF CELL IS NOT ACTIVE GO ON TO NEXT CELL.
        IF (IBOUND(N).LE.0) CYCLE
C
C9------CALCULATE FOR ALL CONNECTING FACES.
        DO II = IA(N)+1,IA(N+1)-1
          JJ = JA(II)
          IIS = JAS(II)
          IF(JJ.LT.N. OR. IBOUND(JJ).LE.0) CYCLE
C10--------COMPUTE AVERAGE DENSITY OF THE TWO CELLS 
          IF(ITHICKAV.EQ.0) THEN 
            RHOTERM = (RHONORM(N) + RHONORM(JJ)) * 0.5
          ELSE
            OMEGA = THICKCELL(N) / (THICKCELL(N) + THICKCELL(JJ))   
            RHOTERM = OMEGA * RHONORM(N) + (1.0 - OMEGA) * RHONORM(JJ)
          ENDIF
C11 --------ADJUST MATRIX FOR DENSITY TERM DEPENDING ON IMPLICIT OR EXPLICIT         
          IF(IMPHDD. EQ. 0) THEN 
C11A -------FOR EXPLICIT UPDATE ON RHS ONLY ADJUST BY * RHOTERM TO ORIGINAL K-VALUE              
            AMAT_TERM = AMAT(II) 
            AMAT(II) = AMAT_TERM * RHOTERM
            AMAT(IA(N)) = AMAT(IA(N)) + AMAT_TERM * (1.0 - RHOTERM)
C            
            AMAT_TERM = AMAT(ISYM(II))            
            AMAT(ISYM(II)) = AMAT_TERM * RHOTERM
            AMAT(IA(JJ)) = AMAT(IA(JJ)) + AMAT_TERM * (1.0 - RHOTERM)
          ELSE  
C11B -------FOR IMPLICIT UPDATE NEED TO DESTROY MATRIX SYMMETRY
            GRHO = RHONORM(N)-RHONORM(JJ)
            AMAT_TERM = AMAT(II) 
            AMAT(II) = AMAT_TERM * (RHOTERM - OMEGA * GRHO)
            AMAT(IA(N)) = AMAT(IA(N)) + 
     1         AMAT_TERM * (1.0 - RHOTERM - GRHO*(1.0-OMEGA))
C            
            GRHO = -GRHO
            AMAT_TERM = AMAT(ISYM(II))            
            AMAT(ISYM(II)) = AMAT_TERM * (RHOTERM - GRHO*(1.0-OMEGA))
            AMAT(IA(JJ)) = AMAT(IA(JJ)) + 
     1         AMAT_TERM * (1.0 - RHOTERM - GRHO*OMEGA) 
          ENDIF     
        ENDDO
C
      ENDDO
C----------------------------------------------------------------------
C12------IF THE SIMULATION IS TRANSIENT ADD STORAGE TO RHS
      ISS=ISSFLG(KPER)
c      IF(ISS.NE.0) RETURN  ! need this term even if flow is steady-state
c      RETURN    !skipping small transient density term
C13------FILL TRANSIENT DENSITY TERM FOR EACH CELL      
      DO N=1,NEQS      
C
C14------IF CELL IS NOT ACTIVE GO ON TO NEXT CELL.
        IF (IBOUND(N).EQ.0) CYCLE
C
C15------COMPUTE TERM AND PUT ON RHS
        IF(N.LE.NODES)THEN
          VOL = AREA(N) * (TOP(N) - BOT(N))*PRSITY(N)
        ELSE
          I=N-NODES  
          VOL = AREA(N) * ACLNNDS(I,4)   
        ENDIF
csp        RHO = VOL*PRSITY(N)*Sn(N)/RHONORM(N)*DRHODC
        RHO = VOL*Sn(N)/RHONORM(N)/RHOFRESH*DRHODC
        CNC = CONC(N,1)
        IF(CNC.LT.0.0) CNC = 0.0
        IF(CNC.GT.CSTD) CNC = CSTD
        CNCO = CONCO(N,1)
        IF(CNCO.LT.0.0) CNCO = 0.0
        IF(CNCO.GT.CSTD) CNCO = CSTD        
        RHO = RHO * (CNC - CNCO) / DELT
        RHS(N) = RHS(N) + RHO
      ENDDO
c16----return
      RETURN
      END
C---------------------------------------------------------------------------
      SUBROUTINE DDF1ADJ(KSTP,KPER)
C     ******************************************************************
C     COMPUTE FLOW BETWEEN ADJACENT CELLS DUE TO DENSITY TERM
C     FOR TRANSIENT SIMULATION ALSO ADD DENSITY STORAGE TERM      
C     ******************************************************************
C
C     SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL, ONLY:IBOUND,NEQS,IOUT,NODES,NJA,IA,JA,JAS,IUNSTR,ISYM,
     1                AMAT,FLOWJA,Sn,ISSFLG,AREA,TOP,BOT,BUFF,HNEW
      USE GWFBASMODULE,ONLY:MSUM,ICBCFL,VBVL,VBNM,DELT,PERTIM,TOTIM
      USE SMSMODULE, ONLY: AMATFL
      USE GWTBCTMODULE, ONLY: PRSITY,CONC,CONCO
      USE DDFMODULE, ONLY: ZEECELL,RHONORM,DRHODC,RHOFRESH,CSTD, 
     1  THICKCELL, ITHICKAV, IMPHDD, ISHARP                                     !  , AMATDD 
      USE GWFBCFMODULE,ONLY:IBCFCB
      USE CLN1MODULE, ONLY: ACLNNDS
C
       double precision,    DIMENSION(:),    ALLOCATABLE ::temp
      DOUBLE PRECISION QCON,ZEENJJ,RHO,VOL,STOIN,STOUT,SSTRG,STIN,SOUT,
     1  ZERO,HPHI,RHOTERM,CNC,CNCO,STRG,GRHO,AMAT_TERM
      CHARACTER*16 TEXT
      DATA TEXT /' DENSITY STORAGE'/
C     ------------------------------------------------------------------
      IF (ISHARP. NE. 0) RETURN 
C1------INITIALIZE BUDGET ACCUMULATORS AND 1/DELT.
      ZERO=0.
      STOIN=ZERO
      STOUT=ZERO      
C2------FILL DENSITY GRADIENT TERM FOR EVERY CELL
      DO N=1,NEQS
C
C3------IF CELL IS NOT ACTIVE GO ON TO NEXT CELL.
        IF (IBOUND(N).LE.0) CYCLE
C
C4------CALCULATE FOR ALL CONNECTING FACES.
        DO II = IA(N)+1,IA(N+1)-1
          JJ = JA(II)
          IIS = JAS(II)
          IF(JJ.LT.N. OR. IBOUND(JJ).LE.0) CYCLE 
C5--------COMPUTE AVERAGE DENSITY, Z AND H OF THE TWO CELLS
          IF(ITHICKAV.EQ.0) THEN 
            ZEENJJ = (ZEECELL(N) + ZEECELL(JJ)) * 0.5
            HPHI = (HNEW(N) + HNEW(JJ)) * 0.5  
            RHOTERM = (RHONORM(N) + RHONORM(JJ)) * 0.5
          ELSE
            OMEGA = THICKCELL(N) / (THICKCELL(N) + THICKCELL(JJ))  
            ZEENJJ = (1.0 - OMEGA) * ZEECELL(N) + OMEGA * ZEECELL(JJ) 
            HPHI = (1.0 - OMEGA) * HNEW(N) + OMEGA * HNEW(JJ) 
            RHOTERM = OMEGA * RHONORM(N) + (1.0 - OMEGA) * RHONORM(JJ) 
          ENDIF  
C6 -------BACK CALCULATE ORIGINAL K VALUE FOR CELLS N AND JJ (amat has K * rho/rho_o)  
            AMAT_TERM = AMATFL(II) / RHOTERM
C 
C7--------CALCULATE FLOW DUE TO DENSITY TERM          
          QCON = AMAT_TERM * (HPHI - ZEENJJ) 
     *         * (RHONORM(N)-RHONORM(JJ))
C8--------SUBTRACT THIS FLOW TERM ON FLOWJA (MINUS ON LHS, PLUS ON RHS)         
c          temp(ii) = qcon
c          temp(isym(ii)) = -qcon
          FLOWJA(II) = FLOWJA(II) + QCON
          FLOWJA(ISYM(II)) = FLOWJA(ISYM(II)) - QCON 
        ENDDO
C
      ENDDO
C9------IF THE SIMULATION IS TRANSIENT ADD STORAGE TO RHS
      ISS=ISSFLG(KPER)
c      IF(ISS.NE.0) GO TO 400   ! need this term even if flow is steady-state
c      GO TO 400  !skipping small transient density term
C
C10------IF CELL-BY-CELL FLOWS WILL BE SAVED, SET FLAG IBD.
      IBD=0
      IF(IBCFCB.GT.0) IBD=ICBCFL 
C
C12------FILL TRANSIENT DENSITY TERM FOR EACH CELL      
      DO N=1,NEQS      
C
C13------IF CELL IS NOT ACTIVE GO ON TO NEXT CELL.
        IF (IBOUND(N).EQ.0) CYCLE
C
C14------COMPUTE TERM AND PUT ON RHS
        IF(N.LE.NODES)THEN
          VOL = AREA(N) * (TOP(N) - BOT(N))*PRSITY(N)
        ELSE
          I=N-NODES  
          VOL = AREA(N) * ACLNNDS(I,4)   
        ENDIF
csp        RHO = VOL*PRSITY(N)*Sn(N)/RHONORM(N)*DRHODC
        RHO = VOL*Sn(N)/RHONORM(N)/RHOFRESH*DRHODC
        CNC = CONC(N,1)
        IF(CNC.LT.0.0) CNC = 0.0
        IF(CNC.GT.CSTD) CNC = CSTD
        CNCO = CONCO(N,1)
        IF(CNCO.LT.0.0) CNCO = 0.0
        IF(CNCO.GT.CSTD) CNCO = CSTD        
        RHO = RHO * (CNC - CNCO) / DELT
        STRG = -RHO
C
C15-----STORE CELL-BY-CELL FLOW IN BUFFER AND ADD TO ACCUMULATORS.
        FLOWJA(IA(N)) = FLOWJA(IA(N)) - STRG
      ENDDO  
c18----return
      RETURN
      END
C---------------------------------------------------------------------------
      SUBROUTINE DDF1BD(KSTP,KPER)
C     ******************************************************************
C     COMPUTE FLOW BETWEEN ADJACENT CELLS DUE TO DENSITY TERM
C     FOR TRANSIENT SIMULATION ALSO ADD DENSITY STORAGE TERM      
C     ******************************************************************
C
C     SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL, ONLY:IBOUND,NEQS,IOUT,NODES,NJA,IA,JA,JAS,IUNSTR,ISYM,
     1                AMAT,FLOWJA,Sn,ISSFLG,AREA,TOP,BOT,BUFF,HNEW
      USE GWFBASMODULE,ONLY:MSUM,ICBCFL,VBVL,VBNM,DELT,PERTIM,TOTIM
      USE SMSMODULE, ONLY: AMATFL
      USE GWTBCTMODULE, ONLY: PRSITY,CONC,CONCO
      USE DDFMODULE, ONLY: ZEECELL,RHONORM,DRHODC,RHOFRESH,CSTD, 
     1  THICKCELL, ITHICKAV, IMPHDD, ISHARP                                     !  , AMATDD 
      USE GWFBCFMODULE,ONLY:IBCFCB
      USE CLN1MODULE, ONLY: ACLNNDS
C
       double precision,    DIMENSION(:),    ALLOCATABLE ::temp
      DOUBLE PRECISION QCON,ZEENJJ,RHO,VOL,STOIN,STOUT,SSTRG,STIN,SOUT,
     1  ZERO,HPHI,RHOTERM,CNC,CNCO,STRG,GRHO,AMAT_TERM
      CHARACTER*16 TEXT
      DATA TEXT /' DENSITY STORAGE'/
C     ------------------------------------------------------------------
      IF (ISHARP. NE. 0) RETURN 
C1------INITIALIZE BUDGET ACCUMULATORS AND 1/DELT.
      ZERO=0.
      STOIN=ZERO
      STOUT=ZERO
C9------IF THE SIMULATION IS TRANSIENT ADD STORAGE TO RHS
      ISS=ISSFLG(KPER)
c      IF(ISS.NE.0) GO TO 400   ! need this term even if flow is steady-state
c      GO TO 400  !skipping small transient density term
C
C10------IF CELL-BY-CELL FLOWS WILL BE SAVED, SET FLAG IBD.
      IBD=0
      IF(IBCFCB.GT.0) IBD=ICBCFL 
C
C12------FILL TRANSIENT DENSITY TERM FOR EACH CELL      
      DO N=1,NEQS      
C
C13------IF CELL IS NOT ACTIVE GO ON TO NEXT CELL.
        IF (IBOUND(N).EQ.0) CYCLE
C
C14------COMPUTE TERM AND PUT ON RHS
        IF(N.LE.NODES)THEN
          VOL = AREA(N) * (TOP(N) - BOT(N))*PRSITY(N)
        ELSE
          I=N-NODES  
          VOL = AREA(N) * ACLNNDS(I,4)   
        ENDIF
csp        RHO = VOL*PRSITY(N)*Sn(N)/RHONORM(N)*DRHODC
        RHO = VOL*Sn(N)/RHONORM(N)/RHOFRESH*DRHODC
        CNC = CONC(N,1)
        IF(CNC.LT.0.0) CNC = 0.0
        IF(CNC.GT.CSTD) CNC = CSTD
        CNCO = CONCO(N,1)
        IF(CNCO.LT.0.0) CNCO = 0.0
        IF(CNCO.GT.CSTD) CNCO = CSTD        
        RHO = RHO * (CNC - CNCO) / DELT
        STRG = -RHO
C
C15-----ADD CELL-BY-CELL FLOW TO BUFF AND ACCUMULATORS.
        BUFF(N)= BUFF(N) + STRG
CCSP        FLOWJA(IA(N)) = FLOWJA(IA(N)) - STRG ! DONE IN DDF1ADJ
        SSTRG=STRG
        IF(STRG.LT.ZERO) THEN
          STOUT=STOUT-SSTRG
        ELSE
          STOIN=STOIN+SSTRG
        END IF
      ENDDO  
C16------record contents of buffer for structured and unstructured grids
CSP      IF(IUNSTR.EQ.0)THEN
CSP
C16A-----IF IBD FLAG IS SET RECORD THE CONTENTS OF THE BUFFER.
CSP        IF(IBD.EQ.1) CALL UBUDSV(KSTP,KPER,TEXT,
CSP     1                       IBCFCB,BUFF,NCOL,NROW,NLAY,IOUT)
CSP        IF(IBD.EQ.2) CALL UBDSV1(KSTP,KPER,TEXT,IBCFCB,
CSP     1            BUFF,NCOL,NROW,NLAY,IOUT,DELT,PERTIM,TOTIM,IBOUND)
CSP      ELSE
C
C16B-----IF IBD FLAG IS SET RECORD THE CONTENTS OF THE BUFFER.
CSP        IF(IBD.EQ.1) CALL UBUDSVU(KSTP,KPER,TEXT,IBCFCB,BUFF(1),NODES,
CSP     1         IOUT,PERTIM,TOTIM)
CSP        IF(IBD.EQ.2) CALL UBDSV1U(KSTP,KPER,TEXT,IBCFCB,BUFF(1),NODES,
CSP     1     IOUT,DELT,PERTIM,TOTIM,IBOUND,NODES)
CSP      ENDIF
C
C17-----ADD TOTAL RATES AND VOLUMES TO VBVL & PUT TITLE IN VBNM.
  400 CONTINUE
      STIN=STOIN
      SOUT=STOUT
      VBVL(1,MSUM)=VBVL(1,MSUM)+STIN*DELT
      VBVL(2,MSUM)=VBVL(2,MSUM)+SOUT*DELT
      VBVL(3,MSUM)=STIN
      VBVL(4,MSUM)=SOUT
      VBNM(MSUM)=TEXT
      MSUM=MSUM+1
c18----return
      RETURN
      END      
C-----------------------------------------------------------------
      SUBROUTINE DDF1DA
C  Deallocate DDF data 
      USE DDFMODULE
C
        DEALLOCATE(IDDF)
        DEALLOCATE(RHOFRESH,RHOSTD,CSTD,DRHODC)
        IF(ISHARP.NE.0) THEN 
          DEALLOCATE(ZETASWI,RELRHODIFF,ISWIUN)
        ELSE
          DEALLOCATE(RHONORM,ZEECELL)   
          IF(ITHICKAV.NE.0)DEALLOCATE(THICKCELL)
        ENDIF  
        DEALLOCATE (ISHARP,ITHICKAV)  
C
      RETURN
      END      
C---------------------------------------------------------------------------
      SUBROUTINE SDDF1BDADJ
C     ******************************************************************
C     RESET MATRIX TO REMOVE IMPLICIT DENSITY TERM EFFECT
C     FLOW TERM OF MATRIX IS KEPT AS K * (RHO/RHO_o) FOR BCF BUDGET ROUTINES      
C     ******************************************************************
C
C     SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL, ONLY:IBOUND,NEQS,IOUT,NODES,NJA,IA,JA,JAS,IUNSTR,ISYM,
     1                AMAT,Sn,AREA,TOP,BOT,BUFF,HNEW
      USE SMSMODULE, ONLY: AMATFL
      USE DDFMODULE, ONLY: ZEECELL,RHONORM,DRHODC,RHOFRESH,CSTD, 
     1  THICKCELL, ITHICKAV, IMPHDD, ISHARP                                     !  , AMATDD 
      USE GWFBCFMODULE,ONLY:IBCFCB
      USE CLN1MODULE, ONLY: ACLNNDS
C
      DOUBLE PRECISION QCON,ZEENJJ,OMEGA,HPHI,RHOTERM,GRHO,AMAT_TERM
C     ------------------------------------------------------------------
      IF (ISHARP. NE. 0) RETURN       
C1------FILL DENSITY GRADIENT TERM FOR EVERY CELL
      DO N=1,NEQS
C
C2------IF CELL IS NOT ACTIVE GO ON TO NEXT CELL.
        IF (IBOUND(N).LE.0) CYCLE
C
C3------CALCULATE FOR ALL CONNECTING FACES.
        DO II = IA(N)+1,IA(N+1)-1
          JJ = JA(II)
          IIS = JAS(II)
          IF(JJ.LT.N. OR. IBOUND(JJ).LE.0) CYCLE 
C4--------COMPUTE AVERAGE DENSITY, Z AND H OF THE TWO CELLS
          IF(ITHICKAV.EQ.0) THEN 
            ZEENJJ = (ZEECELL(N) + ZEECELL(JJ)) * 0.5
            HPHI = (HNEW(N) + HNEW(JJ)) * 0.5  
C            HPHI = (HOLD(N) + HOLD(JJ)) * 0.5  !TIME LAG THIS TERM FOR STABILITY? USE HOLD? 
            RHOTERM = (RHONORM(N) + RHONORM(JJ)) * 0.5
            OMEGA = 0.5
          ELSE
            OMEGA = THICKCELL(N) / (THICKCELL(N) + THICKCELL(JJ))  
            ZEENJJ = (1.0 - OMEGA) * ZEECELL(N) + OMEGA * ZEECELL(JJ) 
            HPHI = (1.0 - OMEGA) * HNEW(N) + OMEGA * HNEW(JJ) 
C            HPHI = (1.0 - OMEGA) * HOLD(N) + OMEGA * HOLD(JJ)  !TIME LAG THIS TERM FOR STABILITY? USE HOLD?  
            RHOTERM = OMEGA * RHONORM(N) + (1.0 - OMEGA) * RHONORM(JJ) 
          ENDIF  
C5 -------BACK CALCULATE ORIGINAL K VALUE FOR CELLS N AND JJ INTO A CONSTANT VARIABLE 
          IF(IMPHDD.EQ.0) THEN 
            AMAT_TERM = AMATFL(II) / RHOTERM
          ELSE
            GRHO = RHONORM(N)-RHONORM(JJ)  
            AMAT_TERM = AMATFL(II) / (RHOTERM - OMEGA * GRHO) 
          ENDIF                    
C          
C6 --------ADJUST MATRIX FOR DENSITY TERM DEPENDING ON IMPLICIT OR EXPLICIT         
          IF(IMPHDD. EQ. 0) THEN 
C7----------NOTHING TO UPDATE FOR EXPLICIT SINCE FLOW COEFFICIENT IS 
C----------K * (RHO/RHO_o) AND LEAVE THAT FOR BCF BUDGET PACKAGE TO SOLVE
          ELSE 
C8----------BACK CALCULATE K * (RHO/RHO_o) FOR MATRIX TERMS FOR IMPLICIT DENSITY TERM
            GRHO = RHONORM(N)-RHONORM(JJ)
            AMATFL(II) = AMAT_TERM * RHOTERM 
            AMATFL(IA(N)) = AMATFL(IA(N)) + 
     1         AMAT_TERM * GRHO*(1.0-OMEGA) 
C            
            AMATFL(ISYM(II)) = AMAT_TERM * RHOTERM 
            AMATFL(IA(JJ)) = AMATFL(IA(JJ)) + 
     1         AMAT_TERM * GRHO * OMEGA 
          ENDIF               
        ENDDO
C
      ENDDO
C---------------------------------------------------------------------------
C9----RETURN
      RETURN
      END
C -------------------------------------------------------------------------------
      SUBROUTINE SHARPINT1AR
C     ******************************************************************
C     ALLOCATE SPACE AND READ INFORMATION FOR SHARP INTERFACE FLOW MODEL
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------c
      USE DDFMODULE
      USE GLOBAL,ONLY:IUNIT,IOUT,NEQS,NODES,IFREFM,IUNSTR,INDDF,HNEW,
     1  NJA,AREA,BOT,TOP,NODLAY,NROW,NCOL,IBOUND,INCLN,SN,SO,NLAY
      USE GWTBCTMODULE, ONLY: CONC
      USE GWFBCFMODULE, ONLY: SC2,LAYCON
      USE CLN1MODULE, ONLY: ACLNNDS,NCLNNDS,IFLINCLN      
      DOUBLE PRECISION ZEE,RNORM,THICK
      DOUBLE PRECISION,    DIMENSION(:,:),    ALLOCATABLE ::TEMP 
      DOUBLE PRECISION ZETA,FRAD,BBOT,TTOP,HD,THCK,TOTTHICK
      CHARACTER*400 LINE
      CHARACTER*24 ANAME(2)
      DATA ANAME(1) /'            SWI POROSITY'/
      DATA ANAME(2) /'      SWI INTERFACE ELEV'/
C
C     ------------------------------------------------------------------
C1A---------ALLOCATE ARRAYS AND COMPUTE BASIC VARIABLES       
      ALLOCATE (RELRHODIFF,HSEA,ISWIUN)  
      ALLOCATE(ZETASWI(NEQS))
      RELRHODIFF = (RHOSTD-RHOFRESH)/RHOFRESH     
C---------------------------------------------------------------------------         
C2------------READ FLAGS AND OPTIONS
        CALL URDCOM(INDDF,IOUT,LINE)          
        IF(IFREFM.EQ.0) THEN
          READ(LINE,'(F10.4,I10)')HSEA,ISWIUN
          LLOC=11
        ELSE
          LLOC=1
          CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,I,TEMPVAR,IOUT,INDDF)    
          CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,ISWIUN,R,IOUT,INDDF)
          HSEA = TEMPVAR
        END IF   
      IF(ISWIUN.GT.0) WRITE(IOUT,11) ISWIUN
11    FORMAT(1X,'IMMOBILE DOMAIN DRAWDOWN WILL BE SAVED ON UNIT',1X,
     1 '(IDPNDD)         =',I4)        
C---------------------------------------------------------------------------
C3-----------READ INTERFACE ELEVATION FOR CLN CELLS
        IF(INCLN.NE.0) THEN 
          CALL U1DREL8(ZETASWI(NODES+1),ANAME(2),NCLNNDS,0,INDDF,IOUT)
        ENDIF 
C4---------------COMPUTE INTERFACE ELEVATION AT EQUILIBRIUM WITH FRESHWATER HEAD       
        CALL ZETACALC              
C---------------------------------------------------------------------------        
C5--------INITIALIZE SN AND SO FOR FRESHWATER
        DO K=1,NLAY
          IF(LAYCON(K).GE.4) THEN
C---------LOOP THROUGH EACH CELL IN LAYER
            NNDLAY = NODLAY(K)
            NSTRT = NODLAY(K-1)+1
            DO N=NSTRT,NNDLAY
              IF(IBOUND(N).NE.0) THEN
C-------------CALCULATE SATURATED THICKNESS.
                HD=HNEW(N)
                BBOT=BOT(N)
                TTOP=TOP(N)
                TOTTHICK = TTOP - BBOT
                CALL SAT_THIK(N,HD,TOTTHICK,BBOT,THCK,K,TTOP)
                Sn(N)=THCK
                So(N) = Sn(N)
              ENDIF
            ENDDO
          ENDIF
        ENDDO
C ------------------------------------------------------------------------          
        IF(INCLN.GT.0) THEN   
C6--------CALCULATE SATURATED THICKNESS FOR CLN-NODES 
          DO  IFN=1,NCLNNDS
            N = ACLNNDS(IFN,1)
            IFLIN = IFLINCLN(IFN)
            IF(IBOUND(N).NE.0.AND.IFLIN.LE.0) THEN
              HD=HNEW(N)
              BBOT = ACLNNDS(IFN,5)
              CALL CLN_THIK(IFN,HD,BBOT,THCK)
              Sn(N)=THCK
              So(N) = Sn(N)
            ENDIF
          ENDDO
        ENDIF  
C---------------------------------------------------------------------------
C7-----RETURN
      RETURN
      END
C -------------------------------------------------------------------------------        
      SUBROUTINE ZETACALC 
C     ******************************************************************
C     CALCULATE INTERFACE ELEVATION FOR SHARP INTERFACE FLOW MODEL AT ALL NODES
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------c
      USE DDFMODULE
      USE GLOBAL,ONLY:NEQS,NODES,INDDF,HNEW,BOT,TOP,IBOUND,INCLN
      USE CLN1MODULE, ONLY: ACLNNDS,NCLNNDS,IFLINCLN      
      DOUBLE PRECISION ZETA,FRAD,BBOT,TTOP,HD 
C
C     ------------------------------------------------------------------  
      IF(ISHARP. EQ.0) RETURN
C1A------------FOR GROUNDWATER CELLS            
      DO N=1,NODES
        IF(IBOUND(N).EQ.0) CYCLE   
        HD = HNEW(N)  
        CALL ZETANODE (N,ZETA,HD)
        ZETASWI(N) = ZETA         
      ENDDO
C1B-----------FOR CLN CELLS          
      IF(INCLN.NE.0) THEN 
        DO IFN=1,NCLNNDS
          N = ACLNNDS(IFN,1)
          IFLIN = IFLINCLN(IFN)
          IF(IBOUND(N).EQ.0.AND.IFLIN.GT.0) CYCLE
          HD = HNEW(N)  
          CALL ZETANODE (N,ZETA,HD) 
          ZETASWI(N) = ZETA           
        ENDDO    
      ENDIF      
C---------------------------------------------------------------------------
C2-----RETURN
      RETURN
      END          
C---------------------------------------------------------------------------          
      SUBROUTINE ZETANODE (N,ZETA,HD)
C     ******************************************************************
C     CALCULATE INTERFACE ELEVATION FOR SHARP INTERFACE FLOW MODEL FOR A NODE
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------c
      USE DDFMODULE
      USE GLOBAL,ONLY:NEQS,NODES,INDDF,HNEW,BOT,TOP,IBOUND,INCLN
      USE CLN1MODULE, ONLY: ACLNNDS,NCLNNDS,IFLINCLN      
      DOUBLE PRECISION ZETA,FRAD,BBOT,TTOP,HD 
C
C     ------------------------------------------------------------------      
C1A------------FOR GROUNDWATER CELLS            
      IF(N. LE. NODES) THEN
        IF(IBOUND(N).NE.0) THEN  
          ZETA = (RHOSTD/RHOFRESH*HSEA - HD) / RELRHODIFF    
          IF(ZETA.LT.BOT(N)) ZETA = BOT(N)                         ! COULD SMOOTHEN AT BOTTOM AND TOP
          IF(ZETA.GT.TOP(N)) ZETA = TOP(N)                         ! COULD SMOOTHEN AT BOTTOM AND TOP
        ENDIF
      ELSE
C1B-----------FOR CLN CELLS          
        IFN = N - NODES
        IFLIN = IFLINCLN(IFN)
        IF(IBOUND(N).NE.0) THEN
C2---------GET TOP AND BOTTOM OF CLN CELL DEPENDING ON ORIENTATION
          BBOT = ACLNNDS(IFN,5)            
          IFDIR = ACLNNDS(IFN,3)
          IF(IFDIR.EQ.0)THEN
C2A-------VERTICAL LINE SEGMENT
            TTOP = ACLNNDS(ICLN,4) + BBOT 
          ELSEIF(IFDIR.EQ.1) THEN 
C2B-------HORIZONTAL LINE SEGMENT 
            IC = ACLNNDS(ICLN,2)
            CALL CLNR(IC,FRAD)
            TOP = 2.0 * FRAD + BBOT
          ELSEIF(IFDIR.EQ.2)THEN  
C2C-------ANGLED LINE SEGMENT                 
            FANGLE = ACLNNDS(ICLN,6)
            TTOP = ACLNNDS(ICLN,4) * SIN(FANGLE) +BBOT
          ENDIF  
C3--------CALCULATE INITIAL INTERFACE ELEVATION IN CLN CELL        
          ZETA = (RHOSTD/RHOFRESH*HSEA - HD) / RELRHODIFF        
          IF(ZETA.LT.BBOT) ZETA = BBOT                        ! COULD SMOOTHEN AT BOTTOM AND TOP
          IF(ZETA.GT.TTOP) ZETA = TTOP                        ! COULD SMOOTHEN AT BOTTOM AND TOP         
        ENDIF            
      ENDIF    
C---------------------------------------------------------------------------
C4-----RETURN
      RETURN
      END          
C---------------------------------------------------------------------------          
      SUBROUTINE SWI1OT (KSTP,KPER,ICNVG,ISA)
C     ******************************************************************
C     OUTPUT INTERFACE ELEVATION FOR SHARP INTERFACE FLOW MODEL
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------c
      USE DDFMODULE,   ONLY: ZETASWI
      USE GLOBAL,      ONLY:ITMUNI,IOUT,IUNSTR
      USE GWFBASMODULE,ONLY:DELT,PERTIM,TOTIM,IHDDFL,IBUDFL,
     1                      MSUM,VBVL,VBNM
C     ------------------------------------------------------------------
C
C
C1------CLEAR PRINTOUT FLAG (IPFLG)
      IPFLG=0
C
C2------IF ITERATIVE PROCEDURE FAILED TO CONVERGE PRINT MESSAGE
      IF(ICNVG.EQ.0) THEN
         WRITE(IOUT,17) KSTP,KPER
   17    FORMAT(1X,/11X,'****FAILED TO CONVERGE IN TIME STEP',I3,
     1      ' OF STRESS PERIOD ',I4,'****')
         IPFLG=1
      END IF
C
C3------IF HEAD FLAG (IHDDFL) IS SET WRITE SWI ELEVATION
C3------IN ACCORDANCE WITH FLAGS IN IOFLG.
      IF(IHDDFL.EQ.0) GO TO 100
C3A-----FOR POROUS MATRIX NODES
      IF(IUNSTR.EQ.0)THEN ! WRITE M2K5 STYLE FOR STRUCTURED GRID
        CALL SGWF2SWI1OT(KSTP,KPER,IPFLG,ISA)
      ELSE
        CALL SGWF2SWI1OTU(KSTP,KPER,IPFLG,ISA)
      ENDIF
  100 CONTINUE
      CALL SGWF2BAS7T(KSTP,KPER,DELT,PERTIM,TOTIM,ITMUNI,IOUT)
      WRITE(IOUT,101)
  101 FORMAT('1')
C
C6------RETURN
      RETURN
      END
C----------------------------------------------------------------------    

      SUBROUTINE SGWF2SWI1OT(KSTP,KPER,IPFLG,ISA)
C     ******************************************************************
C     PRINT AND RECORD SWI INTERFACE ELEVATION FOR STRUCTURED GWF GRID
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,      ONLY:NCOL,NROW,NLAY,IXSEC,HNEW,NODLAY,
     1                      IBOUND,IOUT,IDPOUT
      USE GWFBASMODULE,ONLY:PERTIM,TOTIM,IHEDFM,IHEDUN,LBHDSV,
     2                      CHEDFM,IOFLG
      USE DDFMODULE, ONLY: ISWIUN,ZETASWI
C
      REAL,          SAVE,    DIMENSION(:,:,:),    ALLOCATABLE ::BUFF 
      REAL*8,        SAVE,    DIMENSION(:,:,:),    ALLOCATABLE ::BUFF8 
      CHARACTER*16 TEXT
      DATA TEXT /' INTERFACE ELEV.'/
C     ------------------------------------------------------------------
      ALLOCATE(BUFF(NCOL,NROW,NLAY))
      ALLOCATE(BUFF8(NCOL,NROW,NLAY))
C
C1------FOR EACH LAYER MOVE HNEWIM TO BUFF IF PRINT OR SAVE IS REQUESTED.
      DO 59 K=1,NLAY
C
C2------IS HEAD NEEDED FOR THIS LAYER?
      KL=K
      IF(IXSEC.NE.0) KL=1
      IF(IOFLG(KL,1).EQ.0 .AND. IOFLG(KL,3).EQ.0) GO TO 59
C
C3------MOVE HNEWIM TO BUFF FOR THE LAYER.
      DO 58 I=1,NROW
      DO 58 J=1,NCOL
      N = (K-1)*NROW*NCOL + (I-1)*NCOL + J
      BUFF8(J,I,K)=ZETASWI(N)
      BUFF(J,I,K) = BUFF8(J,I,K)
   58 CONTINUE
   59 CONTINUE
C
C4------FOR EACH LAYER: DETERMINE IF HEAD SHOULD BE PRINTED.
C4------IF SO THEN CALL ULAPRS OR ULAPRW TO PRINT HEAD.
      IF(ISA.NE.0) THEN
         IF(IXSEC.EQ.0) THEN
           DO 69 K=1,NLAY
           KK=K
           IF(IOFLG(K,1).EQ.0) GO TO 69
           IF(IHEDFM.LT.0) CALL ULAPRS(BUFF(1,1,K),TEXT,KSTP,KPER,
     1               NCOL,NROW,KK,-IHEDFM,IOUT)
           IF(IHEDFM.GE.0) CALL ULAPRW(BUFF(1,1,K),TEXT,KSTP,KPER,
     1               NCOL,NROW,KK,IHEDFM,IOUT)
           IPFLG=1
   69      CONTINUE
C
C4A-----PRINT HEAD FOR CROSS SECTION.
         ELSE
           IF(IOFLG(1,1).NE.0) THEN
             IF(IHEDFM.LT.0) CALL ULAPRS(BUFF,TEXT,KSTP,KPER,
     1                 NCOL,NLAY,-1,-IHEDFM,IOUT)
             IF(IHEDFM.GE.0) CALL ULAPRW(BUFF,TEXT,KSTP,KPER,
     1                 NCOL,NLAY,-1,IHEDFM,IOUT)
             IPFLG=1
           END IF
         END IF
      END IF
C
C5------FOR EACH LAYER: DETERMINE IF HEAD SHOULD BE SAVED ON DISK.
C5------IF SO THEN CALL ULASAV OR ULASV2 TO SAVE HEAD.
      IFIRST=1
      IF(ISWIUN.LE.0) GO TO 80
      IF(IXSEC.EQ.0) THEN
        DO 79 K=1,NLAY
        NSTRT = NODLAY(K-1)+1
        KK=K
        IF(IOFLG(K,3).EQ.0) GO TO 79
        IF(IFIRST.EQ.1) WRITE(IOUT,74) ISWIUN,KSTP,KPER
   74   FORMAT(1X,/1X,'INTERFACE ELEV WILL BE SAVED ON UNIT ',I4,
     1      ' AT END OF TIME STEP ',I3,', STRESS PERIOD ',I4)
        IFIRST=0
        IF(IDPOUT.EQ.1) THEN
           WRITE(IOUT,*)
     1      '  INTERFACE ELEV AND TIME VARIABLES ARE SAVED AS REAL*8' 
           CALL ULASAV8(BUFF8(1,1,K),TEXT,KSTP,KPER,PERTIM,TOTIM,NCOL,
     1                NROW,KK,ISWIUN)        
        ELSEIF(CHEDFM.EQ.' ') THEN
           CALL ULASAV(BUFF(1,1,K),TEXT,KSTP,KPER,PERTIM,TOTIM,NCOL,
     1                NROW,KK,ISWIUN)
        ELSE
           CALL ULASV2(BUFF(1,1,K),TEXT,KSTP,KPER,PERTIM,TOTIM,NCOL,
     1                NROW,KK,ISWIUN,CHEDFM,LBHDSV,IBOUND(NSTRT))
        END IF
   79   CONTINUE
C
C5A-----SAVE HEAD FOR CROSS SECTION.
      ELSE
        IF(IOFLG(1,3).NE.0) THEN
          WRITE(IOUT,74) ISWIUN,KSTP,KPER
          IF(IDPOUT.EQ.1)THEN
           WRITE(IOUT,*)
     1       '  INTERRACE ELEV AND TIME VARIABLES ARE SAVED AS REAL*8' 
             CALL ULASAV8(BUFF8,TEXT,KSTP,KPER,PERTIM,TOTIM,NCOL,
     1                NLAY,-1,ISWIUN)          
          ELSEIF(CHEDFM.EQ.' ') THEN
             CALL ULASAV(BUFF,TEXT,KSTP,KPER,PERTIM,TOTIM,NCOL,
     1                NLAY,-1,ISWIUN)
          ELSE
             CALL ULASV2(BUFF,TEXT,KSTP,KPER,PERTIM,TOTIM,NCOL,
     1                  NLAY,-1,ISWIUN,CHEDFM,LBHDSV,IBOUND)
          END IF
        END IF
      END IF
80    CONTINUE
      DEALLOCATE(BUFF)
      DEALLOCATE(BUFF8)
C
C6------RETURN.
      RETURN
      END
C----------------------------------------------------------------------

      SUBROUTINE SGWF2SWI1OTU(KSTP,KPER,IPFLG,ISA)
C     ******************************************************************
C     PRINT AND RECORD SWI INTERFACE ELEVATION FOR UNSTRUCTURED GWF GRID
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,      ONLY:NCOL,NROW,NLAY,IXSEC,HNEW,NODLAY,
     1                      IBOUND,IOUT,NODES,BUFF,NEQS,IDPOUT
      USE GWFBASMODULE,ONLY:PERTIM,TOTIM,IHEDFM,IHEDUN,LBHDSV,
     2                      CHEDFM,IOFLG
      USE DDFMODULE, ONLY: ISWIUN,ZETASWI,ISHARP
C
      REAL*8,        SAVE,    DIMENSION(:),    ALLOCATABLE ::BUFF8
      CHARACTER*16 TEXT
      DATA TEXT /' INTERFACE ELEV.'/
C     ------------------------------------------------------------------
C
      IF(ISHARP.EQ.0) RETURN
C      
      ALLOCATE(BUFF8(NEQS))  
C1------FOR EACH LAYER MOVE HNEW TO BUFF IF PRINT OR SAVE IS REQUESTED.
      DO 59 K=1,NLAY
C
C2------IS HEAD NEEDED FOR THIS LAYER?
      KL=K
      IF(IXSEC.NE.0) KL=1
      IF(IOFLG(KL,1).EQ.0 .AND. IOFLG(KL,3).EQ.0) GO TO 59
C
C3------MOVE HNEW TO BUFF FOR THE LAYER.
      NNDLAY = NODLAY(K)
      NSTRT = NODLAY(K-1)+1
      DO 58 N=NSTRT,NNDLAY
      BUFF8(N)=ZETASWI(N)
      BUFF(N)=BUFF8(N)
   58 CONTINUE
   59 CONTINUE
C
C4------FOR EACH LAYER: DETERMINE IF HEAD SHOULD BE PRINTED.
C4------IF SO THEN CALL ULAPRU TO PRINT HEAD.
      IF(ISA.NE.0) THEN
         IF(IXSEC.EQ.0) THEN
           DO 69 K=1,NLAY
           KK=K
           IF(IOFLG(K,1).EQ.0) GO TO 69
           NNDLAY = NODLAY(K)
           NSTRT = NODLAY(K-1)+1
           CALL ULAPRU(BUFF,TEXT,KSTP,KPER,
     1           NSTRT,NNDLAY,KK,IABS(IHEDFM),IOUT,PERTIM,TOTIM,NODES)
           IPFLG=1
   69      CONTINUE
C
C4A-----PRINT HEAD FOR CROSS SECTION.
         ELSE
           IF(IOFLG(1,1).NE.0) THEN
           CALL ULAPRU(BUFF,TEXT,KSTP,KPER,
     1           NSTRT,NNDLAY,-1,IABS(IHEDFM),IOUT,PERTIM,TOTIM,NODES)
             IPFLG=1
C
           END IF
         END IF
      END IF
C
C5------FOR EACH LAYER: DETERMINE IF HEAD SHOULD BE SAVED ON DISK.
C5------IF SO THEN CALL ULASAV OR ULASV2 TO SAVE HEAD.
      IFIRST=1
      IF(ISWIUN.LE.0) GO TO 80
      IF(IXSEC.EQ.0) THEN
        DO 79 K=1,NLAY
        KK=K
        IF(IOFLG(K,3).EQ.0) GO TO 79
        NNDLAY = NODLAY(K)
        NSTRT = NODLAY(K-1)+1
        IF(IFIRST.EQ.1) WRITE(IOUT,74) ISWIUN,KSTP,KPER
   74   FORMAT(1X,/1X,'INTERFACE ELEV WILL BE SAVED ON UNIT ',I8,
     1      ' AT END OF TIME STEP ',I8,', STRESS PERIOD ',I8)
        IFIRST=0
        IF(IDPOUT.EQ.1)THEN 
          WRITE(IOUT,*)
     1     '  INTERFACE ELEV AND TIME VARIABLES ARE SAVED AS REAL*8'
          CALL ULASAVU8(BUFF8,TEXT,KSTP,KPER,PERTIM,TOTIM,NSTRT,
     1                NNDLAY,KK,ISWIUN,NODES)
        ELSEIF(CHEDFM.EQ.' ') THEN
           CALL ULASAVU(BUFF,TEXT,KSTP,KPER,PERTIM,TOTIM,NSTRT,
     1                NNDLAY,KK,ISWIUN,NODES)
        ELSE
           CALL ULASV2U(BUFF,TEXT,KSTP,KPER,PERTIM,TOTIM,NSTRT,
     1             NNDLAY,KK,ISWIUN,CHEDFM,LBHDSV,IBOUND(NSTRT),NODES)
        END IF
        IPFLG=1
   79   CONTINUE
C
C5A-----SAVE HEAD FOR CROSS SECTION.
      ELSE
        IF(IOFLG(1,3).NE.0) THEN
          WRITE(IOUT,74) ISWIUN,KSTP,KPER
          IF(IDPOUT.EQ.1) THEN 
            WRITE(IOUT,*)
     1       '  INTERFACE ELEV AND TIME VARIABLES ARE SAVED AS REAL*8'
            CALL ULASAVU8(BUFF8,TEXT,KSTP,KPER,PERTIM,TOTIM,NSTRT,
     1                NNDLAY,-1,ISWIUN,NODES)
          ELSEIF(CHEDFM.EQ.' ') THEN
             CALL ULASAVU(BUFF,TEXT,KSTP,KPER,PERTIM,TOTIM,NSTRT,
     1                NNDLAY,-1,ISWIUN,NODES)
          ELSE
             CALL ULASV2U(BUFF,TEXT,KSTP,KPER,PERTIM,TOTIM,NSTRT,
     1                  NNDLAY,-1,ISWIUN,CHEDFM,LBHDSV,IBOUND,NODES)
          END IF
          IPFLG=1
        END IF
      END IF
C
C6------RETURN.
   80 CONTINUE
      DEALLOCATE(BUFF8)
      RETURN
C
      END
C-----------------------------------------------------------------------




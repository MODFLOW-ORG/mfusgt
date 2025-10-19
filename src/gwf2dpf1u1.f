      MODULE GWFDPFMODULE
        INTEGER,SAVE,POINTER::IDPFCB,IDPFHD,IDPFDD
        INTEGER, SAVE, ALLOCATABLE, DIMENSION (:) :: IBOUNDIM
        REAL, SAVE, ALLOCATABLE, DIMENSION(:) ::SC1IM,SC2IM,PHIF,DDFTR
        REAL, SAVE,ALLOCATABLE,DIMENSION(:)::alphaIM,betaIM,srIM,
     1        brookIM,BPIM
         INTEGER, SAVE, ALLOCATABLE, DIMENSION(:) :: IUZONTABIM
C
        DOUBLE PRECISION, SAVE,ALLOCATABLE,DIMENSION(:) :: HNEWIM,HOLDIM
        DOUBLE PRECISION, SAVE,DIMENSION(:), POINTER :: SNIM,SOIM
        DOUBLE PRECISION, SAVE,ALLOCATABLE,DIMENSION(:) :: DIADDF,RDDF
        DOUBLE PRECISION, SAVE,ALLOCATABLE,DIMENSION(:) :: OFFDDFM
        DOUBLE PRECISION, SAVE,ALLOCATABLE,DIMENSION(:) :: OFFDDFIM
        DOUBLE PRECISION, SAVE,ALLOCATABLE,DIMENSION(:) :: AKRCIM
        DOUBLE PRECISION, SAVE,ALLOCATABLE,DIMENSION(:) :: CBCFIM
      END MODULE GWFDPFMODULE
C
C-----------------------------------------------------------------------
      SUBROUTINE GWF2DPFU1AR(IN)
C     ******************************************************************
C     INITIALIZE VARIABLES AND READ DATA FOR DUAL POROSITY FLOW PROCESS
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL, ONLY:IFREFM,IOUT,NODES,NLAY,IUNSAT,NODLAY,IDPF,IBOUND,
     1            TOP,BOT,PGF,IA,JA,JAS,ITRNSP,IDPIN
      USE GWFBASMODULE, ONLY: HNOFLO
      USE GWFDPFMODULE
      USE GWFBCFMODULE, ONLY: ISFAC,SC1,SC2,LAYCON,IBPN,HK,ITABRICH
      CHARACTER*400 LINE
      INTEGER IFRAHK
C
      REAL, DIMENSION(:),ALLOCATABLE  ::HTMP1
      REAL*8, DIMENSION(:),ALLOCATABLE  ::HTMP18
      DOUBLE PRECISION HD,TOTTHICK,TTOP,BBOT,THCK
      CHARACTER*24 ANAME(12)
      DATA ANAME(1) /'      DPF BOUNDARY ARRAY'/
      DATA ANAME(2) /'        DPF INITIAL HEAD'/
      DATA ANAME(3) /'       FRACTURE POROSITY'/
      DATA ANAME(4) /'DUALDOMAIN TRANSFER RATE'/
      DATA ANAME(5) /' DD PRIMARY STORAGE COEF'/
      DATA ANAME(6) /'DD SECONDARY STORGE COEF'/
      DATA ANAME(7) /'               DPF alpha'/
      DATA ANAME(8) /'                DPF beta'/
      DATA ANAME(9) /'                  DPF sr'/
      DATA ANAME(10) /'               DPF brook'/
      DATA ANAME(11) /' DPF BUBBLING POINT HEAD'/
      DATA ANAME(12) /'DPF UNSAT PARAM ZONE MAP'/
C     ------------------------------------------------------------------
C
C1------IDENTIFY PACKAGE
      IDPF = 1 ! FLAG FOR DUAL POROSITY FLOW IS ON
      WRITE(IOUT,1) IN
    1 FORMAT(1X,/1X,'DPF -- DUAL POROSITY FLOW PACKAGE, VERSION 1',
     1',12/26/2012',/,9X,'INPUT READ FROM UNIT',I3)
C
C2------ALLOCATE VARIABLES AND ARRAYS AND INITIALIZE
      ALLOCATE(IDPFCB,IDPFHD,IDPFDD,IBOUNDIM(NODES))
      ALLOCATE(SC1IM(NODES),SC2IM(NODES),PHIF(NODES),DDFTR(NODES))
      ALLOCATE(ALPHAIM(NODES),BETAIM(NODES),SRIM(NODES),BROOKIM(NODES))
      ALLOCATE(HNEWIM(NODES),HOLDIM(NODES),SNIM(NODES),SOIM(NODES))
      ALLOCATE(DIADDF(NODES),RDDF(NODES),OFFDDFM(NODES),OFFDDFIM(NODES))
      IF(IBPN.GT.0)THEN
        ALLOCATE(BPIM(NODES))
      ENDIF
csp      IF(ITRNSP.GT.0)THEN
        ALLOCATE(CBCFIM(NODES))
        DO N=1,NODES
          CBCFIM(N) = 0.0
        ENDDO
csp      ENDIF
C
C3------READ DPF INFORMATION
      IF(IFREFM.EQ.0) THEN
        READ(IN,2)IDPFCB,IDPFHD,IDPFDD
        LLOC=31
      ELSE
        CALL URDCOM(IN,IOUT,LINE)
        LLOC=1
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IDPFCB,R,IOUT,IN)
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IDPFHD,R,IOUT,IN)
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IDPFDD,R,IOUT,IN)
      ENDIF
C3A--------GET OPTIONS
      IFRAHK=0
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,I,R,IOUT,IN)
      IF(LINE(ISTART:ISTOP).EQ.'FRAHK') THEN
C5D1-------SET FLAG FOR HK TO BE FRACTURE-DOMAIN VALUE.
        IFRAHK = 1
      ENDIF
C-------------------------------------------------------------------
2     FORMAT(3I10)
      IF(IDPFCB.GT.0) WRITE(IOUT,9) IDPFCB
    9 FORMAT(1X,'DUAL POROSITY CELL-BY-CELL FLOWS WILL BE SAVED ON',
     1 1X,'UNIT (IDPNCB) =',I4)
      IF(IDPFHD.GT.0) WRITE(IOUT,10) IDPFHD
10    FORMAT(1X,'IMMOBILE DOMAIN HEAD WILL BE SAVED ON UNIT (IDPNHD)'
     1 ,12X,' =',I4)
      IF(IDPFDD.GT.0) WRITE(IOUT,11) IDPFDD
11    FORMAT(1X,'IMMOBILE DOMAIN DRAWDOWN WILL BE SAVED ON UNIT',1X,
     1 '(IDPNDD)         =',I4)
      IF(IFRAHK.EQ.1) WRITE(IOUT,12) IFRAHK
12    FORMAT(1X,'CONDUCTANCE ARE FOR FRACTURE VOLUME; (IFRAHK)',12X,
     1 '=',I4)
      IF(IFRAHK.EQ.0) WRITE(IOUT,13) IFRAHK
13    FORMAT(1X,'CONDUCTANCE ARE FOR TOTAL (FRACTURE + MATRIX)',
     1            1X,'DOMAIN; (IFRAHK)  =',I4)
C
C-------READ STUFF TABULAR INPUT OF RICHARDS EQUATION
      IF(ITABRICH.NE.0) THEN
C-------DIMENSION AND READ ZONE MAP
        ALLOCATE(IUZONTABIM(NODES))
        CALL U1DINT(IUZONTABIM,ANAME(12),NODES,0,IN,IOUT)
      ENDIF
C
C4------READ BOUNDARY ARRAY(IBOUNDIM).
      DO K = 1,NLAY
        KK = K
        NNDLAY = NODLAY(K)
        NSTRT = NODLAY(K-1)+1
        NDSLAY = NNDLAY - NODLAY(K-1)
        CALL U1DINT(IBOUNDIM(NSTRT),ANAME(1),NDSLAY,K,IN,IOUT)
      ENDDO
C-----------------------------------------------------------------------
      DO N=1,NODES
        IF(IBOUND(N).EQ.0) IBOUNDIM(N) = 0
      ENDDO
C-----------------------------------------------------------------------
C5------READ INITIAL HEADS IN IMMOBILE DOMAIN.
      IF(IDPIN.EQ.0) THEN !-----------------------------SINGLE PRECISION READ
      ALLOCATE(HTMP1(Nodes))
      DO K=1,NLAY
        NNDLAY = NODLAY(K)
        NSTRT = NODLAY(K-1)+1
        CALL U1DREL(Htmp1(NSTRT),ANAME(2),NNDLAY-NSTRT+1,K,IN,IOUT)
      ENDDO
      DO N=1,NODES
        HNEWIM(N) = HTMP1(N)
        IF(IBOUNDIM(N).EQ.0) HNEWIM(N)=HNOFLO
        HOLDIM(N) = HNEWIM(N)
      ENDDO
      DEALLOCATE(HTMP1)
      ELSE  !-----------------------------DOUBLE PRECISION READ
      ALLOCATE(HTMP18(Nodes))
      DO K=1,NLAY
        NNDLAY = NODLAY(K)
        NSTRT = NODLAY(K-1)+1
        CALL U1DREL8(Htmp18(NSTRT),ANAME(2),NNDLAY-NSTRT+1,K,IN,IOUT)
      ENDDO
      DO N=1,NODES
        HNEWIM(N) = HTMP18(N)
        IF(IBOUNDIM(N).EQ.0) HNEWIM(N)=HNOFLO
        HOLDIM(N) = HNEWIM(N)
      ENDDO
      DEALLOCATE(HTMP18)
      ENDIF
C
C-----------------------------------------------------------------------
C6------READ MATERIAL PROPERTIES FOR IMMOBILE DOMAIN.
      NCNVRT=0
      DO 200 K = 1,NLAY
        KK = K
        NNDLAY = NODLAY(K)
        NSTRT = NODLAY(K-1)+1
        NDSLAY = NNDLAY - NODLAY(K-1)
C6A-----READ FRACTURE POROSITY (VOLUME OF FRACTURES/TOTAL VOLUME)
        CALL U1DREL(PHIF(NSTRT),ANAME(3),NDSLAY,K,IN,IOUT)
C6B-----READ DUAL DOMAIN MASS TRANSFER RATE (IS K*A/L).
        CALL U1DREL(DDFTR(NSTRT),ANAME(4),NDSLAY,K,IN,IOUT)
C6C-----READ PRIMARY STORAGE COEFFICIENT INTO ARRAY SC1.
        CALL U1DREL(SC1IM(NSTRT),ANAME(5),NDSLAY,K,IN,IOUT)
C6D-----READ SECONDARY STORAGE COEFFICIENT INTO ARRAY SC2.
        IF(LAYCON(K).NE.0)THEN
          CALL U1DREL(SC2IM(NSTRT),ANAME(6),NDSLAY,K,IN,IOUT)
          NCNVRT=1
        ENDIF
C
C6E-----READ alpha, beta, sr, brook
        IF(LAYCON(K).EQ.5)THEN
        IF(ITABRICH.EQ.0) THEN
          CALL U1DREL(alphaIM(NSTRT),ANAME(7),NDSLAY,K,IN,IOUT)
          CALL U1DREL(betaIM(NSTRT),ANAME(8),NDSLAY,K,IN,IOUT)
          CALL U1DREL(srIM(NSTRT),ANAME(9),NDSLAY,K,IN,IOUT)
          CALL U1DREL(brookIM(NSTRT),ANAME(10),NDSLAY,K,IN,IOUT)
          IF(IBPN.GT.0)THEN
            CALL U1DREL(bPIM(NSTRT),ANAME(11),NDSLAY,K,IN,IOUT)
          ENDIF
        ENDIF
        ENDIF
200   CONTINUE
C
C-----------------------------------------------------------------------
C7------CONVERT STORAGE COEFFICIENTS TO STORAGE CAPACITIES
      IF(ISFAC.EQ.0) THEN
        CALL SGWF2LPFU1SC(SC1IM(1),1)
      ELSE
        CALL SGWF2LPFU1SC(SC1IM(1),0)
      END IF
      IF(NCNVRT.GT.0) CALL SGWF2LPFU1SC(SC2IM(1),0)
C7A-----MULTIPLY DUAL DOMAIN RATE TERM BY VOLUME OF GRID BLOCK TO GIVE TOTAL TRANSFER FOR CELL
      CALL SGWF2LPFU1SC(DDFTR(1),1)
C
C-----------------------------------------------------------------------
C8------MULTIPLY SC1 AND SC2 BY FRACTURE POROSITY TO GIVE APPROPRIATE CAPACITY
C8------MULTIPLY SC1IM AND SC2IM BY MATRIX POROSITY (1-PHIF) TO GIVE APPROPRIATE CAPACITY
      DO N=1,NODES
        SC1(N) = SC1(N) * PHIF(N)
        SC1IM(N) = SC1IM(N) * (1.0 - PHIF(N))
      ENDDO
      IF(NCNVRT.GT.0)THEN
        DO N=1,NODES
          SC2(N) = SC2(N) * PHIF(N)
          SC2IM(N) = SC2IM(N) * (1.0 - PHIF(N))
        ENDDO
      ENDIF
C--------------------------------------------------------------------------------
C9------SET INITIAL SATURATIONS IN IMMOBILE DOMAIN.
      DO K=1,NLAY
C9A-------LOOP THROUGH EACH CELL IN LAYER
        NNDLAY = NODLAY(K)
        NSTRT = NODLAY(K-1)+1
        DO N=NSTRT,NNDLAY
          IF(IBOUNDIM(N).NE.0) THEN
            IF(LAYCON(K).NE.0.OR.LAYCON(K).NE.2) THEN
C9B-------------CALCULATE SATURATION OR SATURATED THICKNESS
              HD=HNEWIM(N)
              BBOT=BOT(N)
              TTOP=TOP(N)
              TOTTHICK = TTOP - BBOT
              CALL SAT_THIKIM(N,HD,TOTTHICK,BBOT,THCK,K,TTOP)
              SnIM(N)=THCK
              SoIM(N) = SnIM(N)
            ELSE
              SnIM(N) = 1.0
              SoIM(N) = 1.0
            ENDIF
          ENDIF
        ENDDO
      ENDDO
C10-----------MULTIPLY PROPERTY-GEOMETRY FACTOR BY FRACTURE POROSITY IF IFRAHK FLAG IS SET
      IF(IFRAHK.EQ.1)THEN
        DO K=1,NLAY
          NNDLAY = NODLAY(K)
          NSTRT = NODLAY(K-1)+1
C1B-------loop over all nodes within each layer
          DO N=NSTRT,NNDLAY
C2----------loop over all connections of node N and fill upper triangle with PGF term
            DO II = IA(N)+1,IA(N+1)-1
              JJ = JA(II)
C3------------only for upper triangle of porous medium nodes
              IF(JJ.LE.N.OR.JJ.GT.NODES) CYCLE
              IIS = JAS(II)
              PGF(IIS) = PGF(IIS) * 0.5 * (PHIF(N) + PHIF(JJ))
            ENDDO
          ENDDO
          ENDDO
          DO N = 1,NODES
            HK(N) = HK(N) * PHIF(N)
          ENDDO
      ENDIF
C
C
C10------RETURN
      RETURN
      END
C
C-----------------------------------------------------------------------
      SUBROUTINE GWF2DPFU1FM(KPER)
C     ******************************************************************
C     FORMULATE DUAL POROSITY TERMS FOR STORAGE AND FLOW
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,  ONLY:NODES,NLAY,IUNSAT,NODLAY,AMAT,IA,HNEW,AKR,
     1  TOP,BOT,IBOUND,ISSFLG
      USE GWFDPFMODULE
      USE GWFBASMODULE, ONLY: DELT
      USE GWFBCFMODULE, ONLY: ISFAC,SC1,SC2,LAYCON,IALTSTO
      USE SMSMODULE, ONLY: NONMETH,DKDH
      DOUBLE PRECISION DDFLOW,RHO,TLED,TTOP,BBOT, TOTTHICK,EPS,HD,
     *  DS,RHO2,RHO1,THCK,SATN,SATO,FTERM,DFTERM,SW,EKR,AKRCIMN
C
C     ------------------------------------------------------------------
      IF(ISSFLG(KPER).EQ.1) RETURN
C1--------INITALIZE ARRAYS
      IF(NONMETH.GT.0)THEN
        ALLOCATE(AKRCIM(NODES))
      ENDIF
      DO N=1,NODES
        DIADDF(N) = 0.0
        RDDF(N) = 0.0
        OFFDDFM(N) = 0.0
        OFFDDFIM(N) = 0.0
      ENDDO
C-----------------------------------------------------------------------------
C2------FILL STORAGE TERM ON DIAGONAL AND RHS OF IMMOBILE DOMAIN EQUATION
C-----------------------------------------------------------------------------
      TLED=1.0D0/DELT
      DO 200 K=1,NLAY
        NNDLAY = NODLAY(K)
        NSTRT = NODLAY(K-1)+1
        IF(LAYCON(K).EQ.0)THEN
C-----------------------------------------------------------------------------
C3-------NON-CONVERTIBLE LAYER, SO USE PRIMARY STORAGE
          DO 140 N=NSTRT,NNDLAY
            IF(IBOUNDIM(N).LE.0) GO TO 140
            RHO=SC1IM(N)*TLED
            DIADDF(N)=DIADDF(N)-RHO
            RDDF(N)=RDDF(N)-RHO*HOLDIM(N)
  140     CONTINUE
        ELSE
C-----------------------------------------------------------------------------
C4-------COMPUTE STORAGE TERM AS PER NEWTON RAPHSON
          DO 190 N=NSTRT,NNDLAY
            IF(IBOUNDIM(N).LE.0) GO TO 190
C4A --------CALCULATE SATURATION OR SATURATED THICKNESS
            HD=HNEWIM(N)
            BBOT=BOT(N)
            TTOP=TOP(N)
            TOTTHICK = TTOP - BBOT
            CALL SAT_THIKIM(N,HD,TOTTHICK,BBOT,THCK,K,TTOP)
            SnIM(N)=THCK
C4B-----------COMPUTE PORE STORAGE TERM
            SATO = SoIM(N)
            BBOT=BOT(N)
            TTOP=TOP(N)
            TOTTHICK = TTOP - BBOT
            SATN = SnIM(N)
            EPS = 1.0E-4
            HD=HNEWIM(N)+ EPS
            CALL SAT_THIKIM(N,HD,TOTTHICK,BBOT,THCK,K,TTOP)
            DS = (THCK - SATN)/EPS
            IF(DS.LT.1.0E-20) DS = 1.0E-20
            RHO2  = SC2IM(N) * DS * TOTTHICK * TLED
            DIADDF(N) = DIADDF(N) - RHO2
            RDDF(N) = RDDF(N) - RHO2*HNEWIM(N) +
     1                SC2IM(N)*TOTTHICK*TLED*(SATN-SATO)
C4B-----COMPUTE COMPRESSIBLE STORAGE TERM
c            RHO1 = SATN * SC1IM(N) * TLED
c            DIADDF(N) = DIADDF(N) - RHO1
c            RDDF(N) = RDDF(N) - RHO1 * HOLDIM(N)
C7C --------CHANGE SATN FOR ALTERNATE COMPRESSIBLE STORAGE FORMULATION
            IF(IALTSTO.EQ.1) THEN  ! COMPRESSIBLE STORAGE KICKS ON ONLY NEAR TOP OF CELL
              EPS = TOTTHICK * 0.001     !SMOOTH AT 0.1 PERCENT OF CELL THICKNESS
              EX = HNEWIM(N) - (TTOP - EPS)
              IF(EX .LT. 0) THEN ! BELOW (TOP - SMOOTH) SO NO COMPRESSIBLE STORAGE
                SATN = 0.0
                DS = 0.0
              ELSEIF(EX. LT. EPS) THEN ! WITHIN SMOOTHENING FUNCTION
                SATN = -2.0/EPS**3 * EX**3 + 3.0/EPS**2 * EX**2
                EX = EX + 0.0001
                DS = -2.0/EPS**3 * EX**3 + 3.0/EPS**2 * EX**2
                DS = (DS - SATN) / 0.0001
              ELSEIF(EX. GE. EPS) THEN  ! APPLY FULL COMPRESSIBLE STORAGE TERM
                SATN = 1.0
                DS = 0.0
              ENDIF
            ENDIF
C4B----COMPUTE COMPRESSIBLE STORAGE TERM VIA NEWTON RAPHSON
            RHO1 = SATN * SC1IM(N) * TLED
            FTERM = RHO1 * (HNEWIM(N) - HOLDIM(N))
            DFTERM = RHO1 + SC1IM(N) * TLED * (HNEWIM(N)-HOLDIM(N)) * DS
            DIADDF(N) = DIADDF(N) - DFTERM
            RDDF(N) = RDDF(N) - DFTERM * HNEWIM(N) + FTERM
  190     CONTINUE
        ENDIF
200   CONTINUE
C
C-----------------------------------------------------------------------------
C5------FILL FLOW TERM ON OFF-DIAGONALS AND DIAGONALS
C-----------------------------------------------------------------------------
      DO 300 K=1,NLAY
        NNDLAY = NODLAY(K)
        NSTRT = NODLAY(K-1)+1
        IF(LAYCON(K).EQ.0)THEN
C-----------------------------------------------------------------------------
C6-------NON-CONVERTIBLE LAYER, SO COPY CONSTANT TERM TO OFF-DIAGONALS AND ADD TO DIAGONALS
          DO 240 N=NSTRT,NNDLAY
            IF(IBOUNDIM(N).NE.0.AND.IBOUND(N).NE.0) THEN
              DDFLOW = DDFTR(N)
              OFFDDFM(N) = DDFLOW
              OFFDDFIM(N) =DDFLOW
              AMAT(IA(N)) =  AMAT(IA(N)) - DDFLOW
              DIADDF(N) = DIADDF(N) - DDFLOW
            ENDIF
240       CONTINUE
        ELSE
C-----------------------------------------------------------------------------
C7-------CONVERTIBLE LAYER IS NONLINEAR, COMPUTE NONLINEAR TERMS
          DO 250 N=NSTRT,NNDLAY
            IF(IBOUNDIM(N).NE.0.AND.IBOUND(N).NE.0) THEN
C
C8------------COMPUTE UPSTREAM KR TERM BETWEEN M AND IM DOMAINS
              IF(HNEW(N).GT.HNEWIM(N))THEN
C8A-------------MOBILE DOMAIN IS UPSTREAM
                AKRCIMN = AKR(N)
              ELSE
C8B-------------IMMOBILE DOMAIN IS UPSTREAM
                BBOT=BOT(N)
                TTOP=TOP(N)
                TOTTHICK = TTOP - BBOT
                HD=HNEWIM(N)
                CALL SAT_THIKIM(N,HD,TOTTHICK,BBOT,SW,K,TTOP)
                CALL KR_CALIM(N,SW,EKR,K,HD)
                AKRCIMN = EKR
              ENDIF
C
C9--------------FILL NONLINEAR TERMS IN OFF-DIAGONALS AND ADD TO DIAGONALS
              DDFLOW = DDFTR(N) * AKRCIMN
              OFFDDFM(N) = DDFLOW
              OFFDDFIM(N) = DDFLOW
              AMAT(IA(N)) =  AMAT(IA(N)) - DDFLOW
              DIADDF(N) = DIADDF(N) - DDFLOW
C10-------------SAVE AKRCM FOR USE IN DERIVATIVE COMPUTATION FOR NEWTON
              IF(NONMETH.GT.0)THEN
                AKRCIM(N) = AKRCIMN
              ENDIF
            ENDIF
250       CONTINUE
        ENDIF
300   CONTINUE
C
C11------RETURN
      RETURN
      END
C
C-----------------------------------------------------------------------
      SUBROUTINE SGLO2SMS1NDPF(KPER)
C     ******************************************************************
C     FORMULATE NEWTON RAPHSON TERMS FOR DUAL POROSITY FLOW
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,  ONLY:IFREFM,IOUT,NODES,NLAY,IUNSAT,NODLAY,RHS,AMAT,
     1    IA,HNEW,IBOUND,BOT,TOP,ISSFLG
      USE GWFDPFMODULE
      USE SMSMODULE, ONLY: DKDH,EPSILON
      USE GWFBCFMODULE, ONLY: ISFAC,SC1,SC2,LAYCON
      DOUBLE PRECISION DDFLOW,RHO,TLED,TTOP,BBOT, TOTTHICK,EPS,HD,
     *  DS,RHO2,RHO1,THCK,SATN,SATO,FTERM,DFTERM,SW,EKR,DKRCIM
C
C     ------------------------------------------------------------------
      IF(ISSFLG(KPER).EQ.1) RETURN
C-----------------------------------------------------------------------------
C1------LOOP OVER NODES TO FILL NEWTON TERMS
C-----------------------------------------------------------------------------
      DO 300 K=1,NLAY
        NNDLAY = NODLAY(K)
        NSTRT = NODLAY(K-1)+1
        IF(LAYCON(K).NE.0)THEN
C-----------------------------------------------------------------------------
C2-------CONVERTIBLE LAYER IS NONLINEAR, COMPUTE NONLINEAR TERMS
          DO 250 N=NSTRT,NNDLAY
            IF(IBOUNDIM(N).GT.0.AND.IBOUND(N).GT.0) THEN
C
C3------------COMPUTE UPSTREAM DKR TERM BETWEEN M AND IM DOMAINS
              IF(HNEW(N).GT.HNEWIM(N))THEN
C3A-------------MOBILE DOMAIN IS UPSTREAM
                DKRCIM = DKDH(N)
              ELSE
C3B-------------IMMOBILE DOMAIN IS UPSTREAM
                BBOT=BOT(N)
                TTOP=TOP(N)
                TOTTHICK = TTOP - BBOT
                HD=HNEWIM(N)+ EPSILON
                CALL SAT_THIKIM(N,HD,TOTTHICK,BBOT,SW,K,TTOP)
                CALL KR_CALIM(N,SW,EKR,K,HD)
                DKRCIM = (EKR - AKRCIM(N))/(EPSILON)
              ENDIF
C
C4--------------FILL NONLINEAR TERMS IN OFF-DIAGONALS AND ADD TO DIAGONALS
              CONSTERM = DDFTR(N) * (HNEWIM(N) - HNEW(N))
              IF(HNEW(N).GT.HNEWIM(N))THEN
C4A-------------MOBILE DOMAIN IS UPSTREAM
                TERM = CONSTERM*DKRCIM
                AMAT(IA(N)) = AMAT(IA(N)) + TERM
                OFFDDFIM(N) = OFFDDFIM(N) - TERM
                RHS(N) = RHS(N) + TERM* HNEW(N)
                RDDF(N) = RDDF(N) - TERM* HNEW(N)
              ELSE
C4B-------------IMMOBILE DOMAIN IS UPSTREAM
                TERM = CONSTERM*DKRCIM
                OFFDDFM(N) = OFFDDFM(N) + TERM
                DIADDF(N) = DIADDF(N) - TERM
                RHS(N) = RHS(N) + TERM*HNEWIM(N)
                RDDF(N) = RDDF(N) - TERM*HNEWIM(N)
              ENDIF
            ENDIF
250       CONTINUE
        ENDIF
300   CONTINUE
C
C5----------DEALLOCATE UNWANTED ARRAYS
      DEALLOCATE(AKRCIM)
C
C6------RETURN
      RETURN
      END
C
C-----------------------------------------------------------------------
      SUBROUTINE GWF2DPFU1BDS(KSTP,KPER)
C     ******************************************************************
C     COMPUTE STORAGE BUDGET TERM FOR IMMOBILE DOMAIN
C     ******************************************************************
C
C     SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,  ONLY:NCOL,NROW,NLAY,ISSFLG,IBOUND,ITRNSP,FLOWJA,IA,
     1      BUFF,TOP,IOUT,NODES,NODLAY,IUNSTR,Sn,So,TOP,BOT,iunsat
      USE GWFBASMODULE,ONLY:MSUM,ICBCFL,VBVL,VBNM,DELT,PERTIM,TOTIM
      USE GWFBCFMODULE,ONLY:LAYCON,IALTSTO
      USE GWFDPFMODULE, ONLY:IDPFCB,SC1IM,SC2IM,HNEWIM,HOLDIM,
     1    IBOUNDIM,SOIM,CBCFIM
C
      CHARACTER*16 TEXT
      DOUBLE PRECISION STOIN,STOUT,SSTRG,SIN,SOUT,TLED,HSING,STRG,
     *  RHO,RHO1,RHO2,SNEW,SOLD,ONE,BBOT,TTOP,TOTTHICK,TP
C
      DATA TEXT /'IMMOBILE STORAGE'/
C     ------------------------------------------------------------------
      ISS=ISSFLG(KPER)
C
C1------INITIALIZE BUDGET ACCUMULATORS AND 1/DELT.
      ZERO=0.
      STOIN=ZERO
      STOUT=ZERO
C2------IF STEADY STATE, STORAGE TERM IS ZERO
      IF(ISS.NE.0) GOTO 400
      ONE=1.0
      TLED=ONE/DELT
C
C3------IF CELL-BY-CELL FLOWS WILL BE SAVED, SET FLAG IBD.
      IBD=0
      IF(IDPFCB.GT.0) IBD=ICBCFL
C
C4------CLEAR BUFFER.
      DO 210 N=1,NODES
      BUFF(N)=ZERO
210   CONTINUE
C
C5------LOOP THROUGH EVERY CELL IN THE GRID.
      KT=0
      DO 300 K=1,NLAY
        LC=LAYCON(K)
        NNDLAY = NODLAY(K)
        NSTRT = NODLAY(K-1)+1
        DO 301 N=NSTRT,NNDLAY
C
C6--------SKIP NO-FLOW AND CONSTANT-HEAD CELLS.
          IF(IBOUNDIM(N).LE.0) GO TO 301
          HSING=HNEWIM(N)
C
C7---------COMPUTE STORAGE
          IF(LC.EQ.0) THEN
C7A---------ONE STORAGE CAPACITY.
            RHO=SC1IM(N)*TLED
            STRG=RHO*HOLDIM(N) - RHO*HSING
          ELSE
C7B---------UNCONFINED AND CONFINED STORAGE
            HSING=HNEWIM(N)
            SOLD = SoIM(N)
            BBOT=BOT(N)
            TTOP=TOP(N)
            TOTTHICK = TTOP - BBOT
            CALL SAT_THIKIM(N,HSING,TOTTHICK,BBOT,SNEW,K,TTOP)
            RHO2 = SC2IM(N) * TLED * (SOLD - SNEW) * TOTTHICK
C7C --------CHANGE SATN FOR ALTERNATE COMPRESSIBLE STORAGE FORMULATION
            IF(IALTSTO.EQ.1) THEN  ! COMPRESSIBLE STORAGE KICKS ON ONLY NEAR TOP OF CELL
              EPS = TOTTHICK * 0.001     !SMOOTH AT 0.1 PERCENT OF CELL THICKNESS
              EX = HNEWIM(N) - (TTOP - EPS)
              IF(EX .LT. 0) THEN ! BELOW (TOP - SMOOTH) SO NO COMPRESSIBLE STORAGE
                SATN = 0.0
              ELSEIF(EX. LT. EPS) THEN ! WITHIN SMOOTHENING FUNCTION
                SATN = -2.0/EPS**3 * EX**3 + 3.0/EPS**2 * EX**2
              ELSEIF(EX. GE. EPS) THEN  ! APPLY FULL COMPRESSIBLE STORAGE TERM
                SATN = 1.0
              ENDIF
            ENDIF
C
            RHO1 = SC1IM(N) * TLED * SNEW * (HOLDIM(N) - HSING)
            STRG = RHO1 + RHO2
          ENDIF
C
C8--------STORE CELL-BY-CELL FLOW IN BUFFER AND ADD TO ACCUMULATORS.
          BUFF(N)=STRG
          IF(ITRNSP.GT.0)THEN
            CBCFIM(N) = STRG
            FLOWJA(IA(N)) = FLOWJA(IA(N)) - STRG
          ENDIF
          SSTRG=STRG
          IF(STRG.LT.ZERO) THEN
            STOUT=STOUT-SSTRG
          ELSE
            STOIN=STOIN+SSTRG
          END IF
C
  301 CONTINUE
  300 CONTINUE
C
C9------record contents of buffer for structured and unstructured grids
      IF(IUNSTR.EQ.0)THEN
C9A-----IF IBD FLAG IS SET RECORD THE CONTENTS OF THE BUFFER.
        IF(IBD.EQ.1) CALL UBUDSV(KSTP,KPER,TEXT,
     1                       IDPFCB,BUFF,NCOL,NROW,NLAY,IOUT)
        IF(IBD.EQ.2) CALL UBDSV1(KSTP,KPER,TEXT,IDPFCB,
     1            BUFF,NCOL,NROW,NLAY,IOUT,DELT,PERTIM,TOTIM,IBOUND)
      ELSE
C9B-----IF IBD FLAG IS SET RECORD THE CONTENTS OF THE BUFFER.
        IF(IBD.EQ.1) CALL UBUDSVU(KSTP,KPER,TEXT,IDPFCB,BUFF(1),NODES,
     1         IOUT,PERTIM,TOTIM)
        IF(IBD.EQ.2) CALL UBDSV1U(KSTP,KPER,TEXT,IDPFCB,BUFF(1),NODES,
     1     IOUT,DELT,PERTIM,TOTIM,IBOUND,NODES)
      ENDIF
C
C10-----ADD TOTAL RATES AND VOLUMES TO VBVL & PUT TITLE IN VBNM.
  400 CONTINUE
      SIN=STOIN
      SOUT=STOUT
      VBVL(1,MSUM)=VBVL(1,MSUM)+SIN*DELT
      VBVL(2,MSUM)=VBVL(2,MSUM)+SOUT*DELT
      VBVL(3,MSUM)=SIN
      VBVL(4,MSUM)=SOUT
      VBNM(MSUM)=TEXT
      MSUM=MSUM+1
C
C11----RETURN.
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE SSMS2DPF1RED(KPER)
C     ******************************************************************
C     REDUCE THE IMMOBILE DOMAIN EQUATION INTO THE MOBILE DOMAIN EQUATION
C     AT EACH DUAL POROSITY CELL
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,       ONLY: NODES,AMAT,RHS,IA,ISSFLG
      USE GWFDPFMODULE, ONLY: DIADDF, RDDF, OFFDDFM,OFFDDFIM,IBOUNDIM
C
C     ------------------------------------------------------------------
       IF(ISSFLG(KPER).EQ.1) RETURN
C-----------------------------------------------------------------------------
C1------LOOP OVER NODES TO FILL NEWTON TERMS
C-----------------------------------------------------------------------------
      DO N=1,NODES
        IF(IBOUNDIM(N).EQ.0) CYCLE
        IPIV = IA(N)
        AMAT(IPIV) = AMAT(IPIV) - OFFDDFIM(N)/DIADDF(N)*OFFDDFM(N)
        RHS(N) = RHS(N) - RDDF(N)/DIADDF(N)*OFFDDFM(N)
      ENDDO
C
C11----RETURN.
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE SSMS2DPF1BKS(KPER)
C     ******************************************************************
C     BACK-SUBSTITUTE INTO IMMOBILE DOMAIN EQUATION TO GET ITS HEAD
C     AT EACH DUAL POROSITY CELL
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL, ONLY:NODES,AMAT,RHS,IA,HNEW,ISSFLG
      USE GWFDPFMODULE,ONLY:DIADDF,RDDF,OFFDDFM,OFFDDFIM,HNEWIM,IBOUNDIM
C
C     ------------------------------------------------------------------
C-----------------------------------------------------------------------------
C1------LOOP OVER NODES TO FILL HNWEIM
C-----------------------------------------------------------------------------
      IF(ISSFLG(KPER).EQ.0)THEN
      DO N=1,NODES
        IF(IBOUNDIM(N).EQ.0) CYCLE
C1B-----COMPUTE HNEWIM FOR TRANSIENT STRESS PERIOD
        HNEWIM(N) = RDDF(N)/DIADDF(N) - OFFDDFIM(N)*HNEW(N)/DIADDF(N)
      ENDDO
      ELSE
C1B-----COPY HNEW INTO HNEWIM FOR STEADY-STATE STRESS PERIOD
        DO N=1,NODES
          HNEWIM(N) = HNEW(N)
        ENDDO
      ENDIF
C
C2----RETURN.
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE RES_FUNCDPF(RES)
C     ******************************************************************
C     COMPUTE RESIDUAL - USING MEAN SQUARE RESIDUAL
C     AND INCLUDING IMMOBILE DOMAIN
C     ******************************************************************
C
C      SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL, ONLY:IBOUND,HNEW,ISYM,IOUT,NEQS,AMAT,RHS,IA,JA,NODES
      USE SMSMODULE, ONLY: IBFLAG,BTOL,BREDUC,NUMTRACK,HCLOSE,HTEMP
      USE GWFDPFMODULE, ONLY: DIADDF,RDDF,OFFDDFM,OFFDDFIM,HNEWIM
      DOUBLE PRECISION ROWSUM,RESIDUAL,RES
C     ------------------------------------------------------------------
C
C1------COMPUTE Q FOR ALL NODES
      RESIDUAL = 0.0
      DO N=1,NEQS
        IF(IBOUND(N).GT.0)THEN
          ROWSUM = 0.0
          DO J = IA(N),IA(N+1)-1
            JJ = JA(J)
            ROWSUM = ROWSUM + AMAT(J) * HNEW(JJ)
          ENDDO
C1A--------ADD RESIDUAL FROM IMMOBILE DOMAIN FLOW TERM
          IF(N.LE.NODES)THEN
            ROWSUM = ROWSUM + OFFDDFM(N)*HNEWIM(N)
          ENDIF
C2----------COMPUTE MEAN SQUARE RESIDUAL FROM Q OF EACH NODE
          RESIDUAL = RESIDUAL +  (ROWSUM - RHS(N))**2
        ENDIF
      ENDDO
C3------COMPUTE FOR IMMOBILE DOMAIN OF DUAL POROSITY
      DO N=1,NODES
        ROWSUM = DIADDF(N)*HNEWIM(N) + OFFDDFIM(N)*HNEW(N) - RDDF(N)
        RESIDUAL = RESIDUAL + ROWSUM**2
      ENDDO
C----------------------------------------------------------------------
      RES = RESIDUAL
C3------RETURN
      RETURN
      END SUBROUTINE RES_FUNCDPF
C------------------------------------------------------------------
      SUBROUTINE GWF2DPF1AD
C     ******************************************************************
C     COPY NEW INTO OLD VARIABLES FOR DUAL POROSITY FLOW
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,      ONLY:NODES
      USE GWFDPFMODULE,ONLY:SOIM,SNIM,HOLDIM,HNEWIM
C     ------------------------------------------------------------------
C
C1----COPY NEW TO OLD VARIABLES
      DO 10 N=1,NODES
      SoIM(N) = SnIM(N)
   10 HOLDIM(N)=HNEWIM(N)
C
C4------RETURN
      RETURN
      END
C
      SUBROUTINE ATS1DPFCT
C     ******************************************************************
C     RESET DUAL POROSITY FLOW PARAMETER VECTORS
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL, ONLY:NODES
      USE GWFDPFMODULE, ONLY: HNEWIM,HOLDIM,SNIM,SOIM
C------------------------------------------------------------------------
C1------RESET IMMOBILE DOMAIN HEAD AND SATURATION VECTORS FROM OLD VALUES
      DO N=1,NODES
        HNEWIM(N) = HOLDIM(N)
        SNIM(N) = SOIM(N)
      ENDDO
C
C2------RETURN.
      RETURN
      END
C -----------------------------------------------------------------------
      SUBROUTINE GWF2DPF1OT(KSTP,KPER,ICNVG,ISA)
C     ******************************************************************
C     OUTPUT HEAD, AND DRAWDOWN OF IMMOBILE DOMAIN
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
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
C3------IF HEAD AND DRAWDOWN FLAG (IHDDFL) IS SET WRITE HEAD,
C3------DRAWDOWN, AND IBOUND IN ACCORDANCE WITH FLAGS IN IOFLG.
      IF(IHDDFL.EQ.0) GO TO 100
C3A-----FOR POROUS MATRIX NODES
      IF(IUNSTR.EQ.0)THEN ! WRITE M2K5 STYLE FOR STRUCTURED GRID
        CALL SGWF2DPF1H(KSTP,KPER,IPFLG,ISA)
        CALL SGWF2DPF1D(KSTP,KPER,IPFLG,ISA)
      ELSE
        CALL SGWF2DPF1HU(KSTP,KPER,IPFLG,ISA)
        CALL SGWF2DPF1DU(KSTP,KPER,IPFLG,ISA)
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
      SUBROUTINE SGWF2DPF1D(KSTP,KPER,IPFLG,ISA)
C     ******************************************************************
C     CALCULATE, PRINT, AND SAVE DRAWDOWNS FOR STRUCTURED GWF GRID
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,      ONLY:NCOL,NROW,NLAY,IXSEC,HNEW,STRT,NODLAY,
     1                      IBOUND,IOUT,SN,iunsat,IDPOUT
      USE GWFBASMODULE,ONLY:PERTIM,TOTIM,IDDNFM,IDDNUN,LBDDSV,
     2                      CDDNFM,IOFLG
      USE GWFDPFMODULE, ONLY: HNEWIM,IDPFDD,IBOUNDIM
      USE GWFBCFMODULE, ONLY: LAYCON
C
      CHARACTER*16 TEXT
      DOUBLE PRECISION SSTRT
      REAL,          SAVE,    DIMENSION(:,:,:),    ALLOCATABLE ::BUFF
      REAL*8,        SAVE,    DIMENSION(:,:,:),    ALLOCATABLE ::BUFF8
C
      DATA TEXT /'IMMOBILEDRAWDOWN'/
C     ------------------------------------------------------------------
C
C1-------ALLOCATE TEMPORARY SPACE
      ALLOCATE(BUFF(NCOL,NROW,NLAY))
      ALLOCATE(BUFF8(NCOL,NROW,NLAY))
C
C2------FOR EACH LAYER CALCULATE DRAWDOWN IF PRINT OR SAVE IS REQUESTED.
      DO 59 K=1,NLAY
C
C3------IS DRAWDOWN NEEDED FOR THIS LAYER?
      KL=K
      IF(IXSEC.NE.0) KL=1
      IF(IOFLG(KL,2).EQ.0 .AND. IOFLG(KL,4).EQ.0) GO TO 59
C
C4------CALCULATE DRAWDOWN FOR THE LAYER.
      IF(LAYCON(K).NE.5)THEN
        DO 58 I=1,NROW
        DO 58 J=1,NCOL
        N = (K-1)*NROW*NCOL + (I-1)*NCOL + J
        BUFF8(J,I,K)=HNEWIM(N)
        SSTRT=STRT(N)
        IF(IBOUNDIM(N).NE.0) BUFF8(J,I,K)=SSTRT-HNEWIM(N)
        BUFF(J,I,K) = BUFF8(J,I,K)
   58   CONTINUE
      ELSE
        DO 57 I=1,NROW
        DO 57 J=1,NCOL
        N = (K-1)*NROW*NCOL + (I-1)*NCOL + J
        BUFF8(J,I,K) = SN(N)
        BUFF(J,I,K) = BUFF8(J,I,K)
   57   CONTINUE
      ENDIF
   59 CONTINUE
C
C5------FOR EACH LAYER: DETERMINE IF DRAWDOWN SHOULD BE PRINTED.
C5------IF SO THEN CALL ULAPRS OR ULAPRW TO PRINT DRAWDOWN.
      IF(ISA.NE.0) THEN
         IF(IXSEC.EQ.0) THEN
           DO 69 K=1,NLAY
           KK=K
           IF(IOFLG(K,2).EQ.0) GO TO 69
           IF(IDDNFM.LT.0) CALL ULAPRS(BUFF(1,1,K),TEXT,KSTP,KPER,
     1                  NCOL,NROW,KK,-IDDNFM,IOUT)
           IF(IDDNFM.GE.0) CALL ULAPRW(BUFF(1,1,K),TEXT,KSTP,KPER,
     1                  NCOL,NROW,KK,IDDNFM,IOUT)
           IPFLG=1
   69      CONTINUE
C
C5A-----PRINT DRAWDOWN FOR CROSS SECTION.
         ELSE
           IF(IOFLG(1,2).NE.0) THEN
             IF(IDDNFM.LT.0) CALL ULAPRS(BUFF,TEXT,KSTP,KPER,
     1                 NCOL,NLAY,-1,-IDDNFM,IOUT)
             IF(IDDNFM.GE.0) CALL ULAPRW(BUFF,TEXT,KSTP,KPER,
     1                 NCOL,NLAY,-1,IDDNFM,IOUT)
             IPFLG=1
           END IF
         END IF
      END IF
C
C6------FOR EACH LAYER: DETERMINE IF DRAWDOWN SHOULD BE SAVED.
C6------IF SO THEN CALL A ULASAV OR ULASV2 TO RECORD DRAWDOWN.
      IFIRST=1
      IF(IDPFDD.LE.0) GO TO 80
      IF(IXSEC.EQ.0) THEN
        DO 79 K=1,NLAY
        NSTRT = NODLAY(K-1)+1
        KK=K
        IF(IOFLG(K,4).EQ.0) GO TO 79
        IF(IFIRST.EQ.1) WRITE(IOUT,74) IDPFDD,KSTP,KPER
   74   FORMAT(1X,/1X,'IMMOBILE DOMAIN DRAWDOWN WILL BE SAVED ON UNIT ',
     1      I4,' AT END OF TIME STEP ',I3,', STRESS PERIOD ',I4)
        IFIRST=0
        IF(IDPOUT. EQ. 1) THEN
           WRITE(IOUT,*)'  DDN AND TIME VARIABLES ARE SAVED AS REAL*8'
           CALL ULASAV8(BUFF8(1,1,K),TEXT,KSTP,KPER,PERTIM,TOTIM,NCOL,
     1                NROW,KK,IDPFDD)
        ELSEIF(CDDNFM.EQ.' ') THEN
           CALL ULASAV(BUFF(1,1,K),TEXT,KSTP,KPER,PERTIM,TOTIM,NCOL,
     1                NROW,KK,IDPFDD)
        ELSE
           CALL ULASV2(BUFF(1,1,K),TEXT,KSTP,KPER,PERTIM,TOTIM,NCOL,
     1                NROW,KK,IDPFDD,CDDNFM,LBDDSV,IBOUND(NSTRT))
        END IF
   79   CONTINUE
C
C6A-----SAVE DRAWDOWN FOR CROSS SECTION.
      ELSE
        IF(IOFLG(1,4).NE.0) THEN
          WRITE(IOUT,74) IDPFDD,KSTP,KPER
          IF(IDPOUT. EQ. 1) THEN
            WRITE(IOUT,*)'  DDN AND TIME VARIABLES ARE SAVED AS REAL*8'
            CALL ULASAV8(BUFF8,TEXT,KSTP,KPER,PERTIM,TOTIM,NCOL,
     1                NLAY,-1,IDPFDD)
          ELSEIF(CDDNFM.EQ.' ') THEN
             CALL ULASAV(BUFF,TEXT,KSTP,KPER,PERTIM,TOTIM,NCOL,
     1                NLAY,-1,IDPFDD)
          ELSE
             CALL ULASV2(BUFF,TEXT,KSTP,KPER,PERTIM,TOTIM,NCOL,
     1                  NLAY,-1,IDPFDD,CDDNFM,LBDDSV,IBOUND)
          END IF
        END IF
      END IF
80    CONTINUE
C
C7------DEALLOCATE TEMPORARY SPACE
      DEALLOCATE(BUFF)
      DEALLOCATE(BUFF8)

C
C8------RETURN.
      RETURN
      END
      SUBROUTINE SGWF2DPF1H(KSTP,KPER,IPFLG,ISA)
C     ******************************************************************
C     PRINT AND RECORD HEADS FOR STRUCTURED GWF GRID
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,      ONLY:NCOL,NROW,NLAY,IXSEC,HNEW,NODLAY,
     1                      IBOUND,IOUT,IDPOUT
      USE GWFBASMODULE,ONLY:PERTIM,TOTIM,IHEDFM,IHEDUN,LBHDSV,
     2                      CHEDFM,IOFLG
      USE GWFDPFMODULE, ONLY: HNEWIM,IDPFHD,IBOUNDIM
C
      REAL,          SAVE,    DIMENSION(:,:,:),    ALLOCATABLE ::BUFF
      REAL*8,        SAVE,    DIMENSION(:,:,:),    ALLOCATABLE ::BUFF8
      CHARACTER*16 TEXT
      DATA TEXT /'   IMMOBILE HEAD'/
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
      BUFF8(J,I,K)=HNEWIM(N)
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
      IF(IDPFHD.LE.0) GO TO 80
      IF(IXSEC.EQ.0) THEN
        DO 79 K=1,NLAY
        NSTRT = NODLAY(K-1)+1
        KK=K
        IF(IOFLG(K,3).EQ.0) GO TO 79
        IF(IFIRST.EQ.1) WRITE(IOUT,74) IDPFHD,KSTP,KPER
   74   FORMAT(1X,/1X,'IMMOBILE DOMAIN HEAD WILL BE SAVED ON UNIT ',I4,
     1      ' AT END OF TIME STEP ',I3,', STRESS PERIOD ',I4)
        IFIRST=0
        IF(IDPOUT.EQ.1) THEN
           WRITE(IOUT,*)'  HEAD AND TIME VARIABLES ARE SAVED AS REAL*8'
           CALL ULASAV8(BUFF8(1,1,K),TEXT,KSTP,KPER,PERTIM,TOTIM,NCOL,
     1                NROW,KK,IDPFHD)
        ELSEIF(CHEDFM.EQ.' ') THEN
           CALL ULASAV(BUFF(1,1,K),TEXT,KSTP,KPER,PERTIM,TOTIM,NCOL,
     1                NROW,KK,IDPFHD)
        ELSE
           CALL ULASV2(BUFF(1,1,K),TEXT,KSTP,KPER,PERTIM,TOTIM,NCOL,
     1                NROW,KK,IDPFHD,CHEDFM,LBHDSV,IBOUND(NSTRT))
        END IF
   79   CONTINUE
C
C5A-----SAVE HEAD FOR CROSS SECTION.
      ELSE
        IF(IOFLG(1,3).NE.0) THEN
          WRITE(IOUT,74) IDPFHD,KSTP,KPER
          IF(IDPOUT.EQ.1)THEN
           WRITE(IOUT,*)'  HEAD AND TIME VARIABLES ARE SAVED AS REAL*8'
             CALL ULASAV8(BUFF8,TEXT,KSTP,KPER,PERTIM,TOTIM,NCOL,
     1                NLAY,-1,IDPFHD)
          ELSEIF(CHEDFM.EQ.' ') THEN
             CALL ULASAV(BUFF,TEXT,KSTP,KPER,PERTIM,TOTIM,NCOL,
     1                NLAY,-1,IDPFHD)
          ELSE
             CALL ULASV2(BUFF,TEXT,KSTP,KPER,PERTIM,TOTIM,NCOL,
     1                  NLAY,-1,IDPFHD,CHEDFM,LBHDSV,IBOUND)
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
      SUBROUTINE SGWF2DPF1DU(KSTP,KPER,IPFLG,ISA)
C     ******************************************************************
C     CALCULATE, PRINT, AND SAVE DRAWDOWNS FOR UNSTRUCTURED GWF GRID
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,      ONLY:NCOL,NROW,NLAY,IXSEC,HNEW,STRT,NODLAY,
     1                 IBOUND,IOUT,BUFF,NODES,IUNSAT,SN,NEQS,IDPOUT
      USE GWFBASMODULE,ONLY:PERTIM,TOTIM,IDDNFM,IDDNUN,LBDDSV,
     2                      CDDNFM,IOFLG
       USE GWFDPFMODULE, ONLY: HNEWIM,IDPFDD,IBOUNDIM
       USE GWFBCFMODULE, ONLY: LAYCON
C
      CHARACTER*16 TEXT
      DOUBLE PRECISION SSTRT
      REAL*8,        SAVE,    DIMENSION(:),    ALLOCATABLE ::BUFF8
C
      DATA TEXT /'IMMOBILDRAWDOWNU'/
C     ------------------------------------------------------------------
C
      ALLOCATE(BUFF8(NEQS))
C1------FOR EACH LAYER CALCULATE DRAWDOWN IF PRINT OR SAVE IS REQUESTED.
      DO 59 K=1,NLAY
C
C2------IS DRAWDOWN NEEDED FOR THIS LAYER?
      KL=K
      IF(IXSEC.NE.0) KL=1
      IF(IOFLG(KL,2).EQ.0 .AND. IOFLG(KL,4).EQ.0) GO TO 59
C
C3------CALCULATE DRAWDOWN FOR THE LAYER.
      NNDLAY = NODLAY(K)
      NSTRT = NODLAY(K-1)+1
      IF(LAYCON(K).NE.5)THEN
        DO 58 N=NSTRT,NNDLAY
          BUFF8(N)=HNEWIM(N)
          SSTRT=STRT(N)
          IF(IBOUND(N).NE.0) BUFF8(N)=SSTRT-HNEWIM(N)
          BUFF(N) = BUFF8(N)
   58   CONTINUE
      ELSE
        DO 57 N=NSTRT,NNDLAY
          BUFF8(N)=SN(N)
          BUFF(N) = BUFF8(N)
   57   CONTINUE
      ENDIF
   59 CONTINUE
C
C4------FOR EACH LAYER: DETERMINE IF DRAWDOWN SHOULD BE PRINTED.
C4------IF SO THEN CALL ULAPRU TO PRINT DRAWDOWN.
      IF(ISA.NE.0) THEN
         IF(IXSEC.EQ.0) THEN
           DO 69 K=1,NLAY
           KK=K
           IF(IOFLG(K,2).EQ.0) GO TO 69
           NNDLAY = NODLAY(K)
           NSTRT = NODLAY(K-1)+1
           CALL ULAPRU(BUFF,TEXT,KSTP,KPER,
     1           NSTRT,NNDLAY,KK,IABS(IDDNFM),IOUT,PERTIM,TOTIM,NODES)
           IPFLG=1
   69      CONTINUE
C
C4A-----PRINT DRAWDOWN FOR CROSS SECTION.
         ELSE
           IF(IOFLG(1,2).NE.0) THEN
           CALL ULAPRU(BUFF,TEXT,KSTP,KPER,
     1           NSTRT,NNDLAY,-1,IABS(IDDNFM),IOUT,PERTIM,TOTIM,NODES)
             IPFLG=1
           END IF
         END IF
      END IF
C
C5------FOR EACH LAYER: DETERMINE IF DRAWDOWN SHOULD BE SAVED.
C5------IF SO THEN CALL A ULASAV OR ULASV2 TO RECORD DRAWDOWN.
      IFIRST=1
      IF(IDPFDD.LE.0) GO TO 80
      IF(IXSEC.EQ.0) THEN
        DO 79 K=1,NLAY
        KK=K
        IF(IOFLG(K,4).EQ.0) GO TO 79
        NNDLAY = NODLAY(K)
        NSTRT = NODLAY(K-1)+1
        IF(IFIRST.EQ.1) WRITE(IOUT,74) IDPFDD,KSTP,KPER
   74   FORMAT(1X,/1X,'IMMOBILE DRAWDOWN WILL BE SAVED ON UNIT ',I4,
     1      ' AT END OF TIME STEP ',I3,', STRESS PERIOD ',I4)
        IFIRST=0
        IF(IDPOUT.EQ.1) THEN
          WRITE(IOUT,*)'  DDN AND TIME VARIABLES ARE SAVED AS REAL*8'
           CALL ULASAVU8(BUFF8,TEXT,KSTP,KPER,PERTIM,TOTIM,NSTRT,
     1                NNDLAY,KK,IDPFDD,NODES)
        ELSEIF(CDDNFM.EQ.' ') THEN
           CALL ULASAVU(BUFF,TEXT,KSTP,KPER,PERTIM,TOTIM,NSTRT,
     1                NNDLAY,KK,IDPFDD,NODES)
        ELSE
           CALL ULASV2U(BUFF,TEXT,KSTP,KPER,PERTIM,TOTIM,NSTRT,
     1             NNDLAY,KK,IDPFDD,CDDNFM,LBDDSV,IBOUND(NSTRT),NODES)
        END IF
   79   CONTINUE
C
C5A-----SAVE DRAWDOWN FOR CROSS SECTION.
      ELSE
        IF(IOFLG(1,4).NE.0) THEN
          WRITE(IOUT,74) IDPFDD,KSTP,KPER
          IF(IDPOUT.EQ.1) THEN
           WRITE(IOUT,*)'  DDN AND TIME VARIABLES ARE SAVED AS REAL*8'
           CALL ULASAVU8(BUFF8,TEXT,KSTP,KPER,PERTIM,TOTIM,NSTRT,
     1                NNDLAY,-1,IDPFDD,NODES)
          ELSEIF(CDDNFM.EQ.' ') THEN
           CALL ULASAVU(BUFF,TEXT,KSTP,KPER,PERTIM,TOTIM,NSTRT,
     1                NNDLAY,-1,IDPFDD,NODES)
        ELSE
           CALL ULASV2U(BUFF,TEXT,KSTP,KPER,PERTIM,TOTIM,NSTRT,
     1             NNDLAY,-1,IDPFDD,CDDNFM,LBDDSV,IBOUND(NSTRT),NODES)
          END IF
        END IF
      END IF
80    CONTINUE
      DEALLOCATE(BUFF8)
C
C6------RETURN.
      RETURN
      END
      SUBROUTINE SGWF2DPF1HU(KSTP,KPER,IPFLG,ISA)
C     ******************************************************************
C     PRINT AND RECORD HEADS FOR UNSTRUCTURED GWF GRID
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,      ONLY:NCOL,NROW,NLAY,IXSEC,HNEW,NODLAY,
     1                      IBOUND,IOUT,NODES,BUFF,NEQS,IDPOUT
      USE GWFBASMODULE,ONLY:PERTIM,TOTIM,IHEDFM,IHEDUN,LBHDSV,
     2                      CHEDFM,IOFLG
      USE GWFDPFMODULE, ONLY: HNEWIM,IDPFHD,IBOUNDIM
C
      REAL*8,        SAVE,    DIMENSION(:),    ALLOCATABLE ::BUFF8
      CHARACTER*16 TEXT
      DATA TEXT /'  IMMOBILE HEADU'/
C     ------------------------------------------------------------------
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
      BUFF8(N)=HNEWIM(N)
      BUFF(N)=HNEWIM(N)
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
      IF(IDPFHD.LE.0) GO TO 80
      IF(IXSEC.EQ.0) THEN
        DO 79 K=1,NLAY
        KK=K
        IF(IOFLG(K,3).EQ.0) GO TO 79
        NNDLAY = NODLAY(K)
        NSTRT = NODLAY(K-1)+1
        IF(IFIRST.EQ.1) WRITE(IOUT,74) IDPFHD,KSTP,KPER
   74   FORMAT(1X,/1X,'IMMOBILE HEAD WILL BE SAVED ON UNIT ',I8,
     1      ' AT END OF TIME STEP ',I8,', STRESS PERIOD ',I8)
        IFIRST=0
        IF(IDPOUT.EQ.1)THEN
          WRITE(IOUT,*)'  HEAD AND TIME VARIABLES ARE SAVED AS REAL*8'
          CALL ULASAVU8(BUFF8,TEXT,KSTP,KPER,PERTIM,TOTIM,NSTRT,
     1                NNDLAY,KK,IDPFHD,NODES)
        ELSEIF(CHEDFM.EQ.' ') THEN
           CALL ULASAVU(BUFF,TEXT,KSTP,KPER,PERTIM,TOTIM,NSTRT,
     1                NNDLAY,KK,IDPFHD,NODES)
        ELSE
           CALL ULASV2U(BUFF,TEXT,KSTP,KPER,PERTIM,TOTIM,NSTRT,
     1             NNDLAY,KK,IDPFHD,CHEDFM,LBHDSV,IBOUND(NSTRT),NODES)
        END IF
        IPFLG=1
   79   CONTINUE
C
C5A-----SAVE HEAD FOR CROSS SECTION.
      ELSE
        IF(IOFLG(1,3).NE.0) THEN
          WRITE(IOUT,74) IDPFHD,KSTP,KPER
          IF(IDPOUT.EQ.1) THEN
            WRITE(IOUT,*)'  HEAD AND TIME VARIABLES ARE SAVED AS REAL*8'
            CALL ULASAVU8(BUFF8,TEXT,KSTP,KPER,PERTIM,TOTIM,NSTRT,
     1                NNDLAY,-1,IDPFHD,NODES)
          ELSEIF(CHEDFM.EQ.' ') THEN
             CALL ULASAVU(BUFF,TEXT,KSTP,KPER,PERTIM,TOTIM,NSTRT,
     1                NNDLAY,-1,IDPFHD,NODES)
          ELSE
             CALL ULASV2U(BUFF,TEXT,KSTP,KPER,PERTIM,TOTIM,NSTRT,
     1                  NNDLAY,-1,IDPFHD,CHEDFM,LBHDSV,IBOUND,NODES)
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
      SUBROUTINE SAT_THIKIM(N,HD,TOTTHICK,BBOT,THCK,K,TTOP)
C     ******************************************************************
C     COMPUTE Kr FOR LAYCON OF 4 AND 5
C     ******************************************************************
C
C      SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,     ONLY:ISSFLG,NODES,iunsat
      USE GWFDPFMODULE,ONLY:alphaIM,betaIM,brookIM,srIM,BPIM,IUZONTABIM
      USE GWFBCFMODULE, ONLY: IBPN,LAYCON,IDRY,RETCRVS,NUTABROWS,NUZONES
     1  ,ITABRICH
      DOUBLE PRECISION THCK,HD,BBOT,TTOP,X,S,V,COF1,COF2,FACTOR1,FACTOR2
     *  ,eps,acof,y,TOTTHICK,seff,gamma, tabS,pc, topS, botS,cellS,
     *   bott,mul,div
C     ------------------------------------------------------------------
C
C ----DEPENDING ON THE METHOD, FIND THE SATURATED THICKNESS
C
C--------------------------------------------------------
      ZERO=0.
      METHOD = 3
      IF(N.LE.NODES)THEN
        IF(LAYCON(K).EQ.5) THEN
          METHOD = 4
          IF(ITABRICH.NE.0)THEN !METHOD 5 USES TABULAR INPUT
            METHOD = 5
            IF(INTRICH.NE.0)THEN  !METHOD 6 USES INTEGRATED FUNCTIONS
              METHOD = 6
            ENDIF
          ENDIF
        ENDIF
      ENDIF
C      if(IUNSAT.EQ.1.AND.N.LE.NODES) method = 4
      IF(METHOD.EQ.1)THEN
C-------STRAIGHT LINE, NO SMOOTHING
        TTOP = BBOT + TOTTHICK
        IF(HD.GT.TTOP) HD=TTOP
        THCK = (HD - BBOT) / TOTTHICK
        IF(THCK.LT.ZERO) THCK=0.0
      ELSEIF(METHOD.EQ.2)THEN
C-------STRAIGHT LINE WITH CUBIC SMOOTHING
        Thickfact = 0.01
        x = (HD-bbot)
        s = Thickfact*TOTTHICK
        v = TOTTHICK
        cof1 = (1.0D0/s**2.0D0)
        cof2 = (2.0D0/s)
        factor1 = -cof1*x**3.0D0+cof2*x**2.0D0
        factor2 = 1.0D0 + cof1*(v-x)**3.0D0-
     +      cof2*(v-x)**2.0D0
        THCK = 0.0D0
        IF ( x.LT.0.0D0 ) THEN
          THCK = 0.0D0
        ELSEIF ( x.LT.s ) THEN
          THCK = factor1
        ELSEIF ( x.LT.v-s ) THEN
          THCK = x/v
        ELSEIF ( x.LT.v ) THEN
          THCK = factor2
        ELSEIF ( x.GE.v ) THEN
          THCK = 1.0D0
        END IF
      ELSEIF(METHOD.EQ.3)THEN
C-------STRAIGHT LINE WITH PARABOLIC SMOOTHING
        EPS = 0.01
        eps = 1.0e-6
        ACOF = 1.0 / (1.0 - EPS)
        x = (HD-bbot)/TOTTHICK
        IF(X.LT.0)THEN
          Y = 0.0
        ELSEIF(X.LT.EPS)THEN
          Y = ACOF *0.5/EPS * X**2
        ELSEIF(X.LT.1.0-EPS)THEN
          Y = ACOF * X + (1.0-ACOF)*0.5
        ELSEIF(X.LT.1.0)THEN
          X = 1.0 - X
          Y = ACOF *0.5/EPS * X**2
          Y = 1.0-Y
        ELSE
          Y = 1.0
        ENDIF
        THCK = Y
      ELSEIF(METHOD.EQ.4)THEN
C---------vanG FUNCTION WITH MODIFICATIONS FOR BUBBLEPT AND FULLYDRY
        bpn = 0.0
        if(ibpn.eq.1)then
          bpn = bpIM(N)
        endif
        TTOP = BBOT + TOTTHICK
        pc = 0.5*(ttop+bbot) - (hd-bpn)
        if(pc.le.0)then
          thck = 1.0
        else
          gamma = 1.-1./betaIM(n)
          Seff = (1. + (alphaIM(n)*pc)**betaIM(n))**gamma
          Seff = 1.0 / Seff
          if(idry.eq.0) then
            thck = seff * (1-srIM(n)) + srIM(n)
          else
            thck = seff
          endif
        endif
C------------------------------------
      ELSEIF(METHOD.EQ.5) THEN
C-------FOR TABULAR INPUT FIND SATURATION FROM TABLE
        TTOP = BBOT + TOTTHICK
        pc = 0.5*(ttop+bbot) - hd
        if(pc.le.0)then
          thck = 1.0
        else
          izon = iuzontabIM(n)
          iv = 2  !variable is saturation
          thck = tabS(pc,retcrvs,izon,iv,nutabrows,nuzones,ttop)
        endif
      ELSEIF(METHOD.EQ.6) THEN
C---------FOR USE WITH TABULAR INPUT OF INTERPOLATED TABLE
        iv = 2  !variable is saturation
        izon = iuzontabIM(n)
        bott = bbot   ! bott changes to water table if below water table in tabS
        topS = tabS((ttop-hd),retcrvs,izon,iv,nutabrows,nuzones,ttop)
        botS = tabS((bbot-hd),retcrvs,izon,iv,nutabrows,nuzones,bott)
        if(hd.lt. bbot) then
c----------water level is below bottom so integrate curve from top to bottom
          if((ttop-bott).gt.1.0e-5)then
            cellS=-(topS - botS)/(ttop - bott)
            if(cellS.lt.1.0e-20) cellS = 1.0e-20
          else ! past bottom of table so use last two entries to compute value
           mul = retcrvs(iv,nutabrows,izon)-retcrvs(iv,nutabrows-1,izon)
           div = retcrvs(1,nutabrows,izon)-retcrvs(1,nutabrows-1,izon)
           cellS = mul/div
          endif
        elseif(hd.gt. ttop) then
c----------water level is above top so integrated saturation of cell is 1
          cellS = 1.0
        else
c----------water level is within cell so integrate
          cellS = (botS - topS + 1.0 * (hd -bbot))/(ttop - bbot)
        endif
        thck = cellS
C------------------------------------
      ENDIF
C
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE KR_CALIM(N,SW,EKR,K,HD)
C     ******************************************************************
C     COMPUTE KR FOR RICHARDS EQUATION FORMULATION
C     ******************************************************************
C
C      SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,     ONLY:ISSFLG,iunsat,TOP,BOT,NODES
      USE GWFBCFMODULE,ONLY: LAYCON,retcrvs,nutabrows,nuzones,
     1   ITABRICH
      USE GWFDPFMODULE,ONLY:alphaIM,betaIM,srIM,brookIM, iuzontabIM
      DOUBLE PRECISION EKR,SW,seff,one, tabS,pc, topS, botS,cellS,
     *    BBOT,TTOP,HD,bott,mul,div
C     ------------------------------------------------------------------
      METHOD = 3  !  METHOD 3 IS LAYCON(K) = 4
      IF(N.LE.NODES)THEN
        IF(LAYCON(K).EQ.5) THEN ! METHOD 4 USES VAN GENUCHTEN AND BROOKS-COREY FUNCTIONS
          METHOD = 4
          IF(ITABRICH.NE.0)THEN !METHOD 5 USES TABULAR INPUT
            METHOD = 5
            IF(INTRICH.NE.0)THEN  !METHOD 6 USES INTEGRATED FUNCTIONS
              METHOD = 6
            ENDIF
          ENDIF
        ENDIF
      ENDIF
C      IF(LAYCON(K).EQ.4)THEN
      IF(METHOD.EQ.3)THEN
        EKR = SW
cc        IF(EKR.LT.1.0E-7) EKR = 1.0E-7
      ELSEIF(METHOD.EQ.4)THEN
C---------BROOKS COREY FUNCTION
        IF(SW.LT.SRIM(N)) THEN
          EKR = 0.0
        ELSE
          one = 1.0
          SEFF = (SW - SrIM(N))/(one - SrIM(N))
cc ......use van genuchten kr when brookimn is less than zero
          IF(BROOKIM(N). LT. 0.0) THEN
           IF(Seff. GT. one) Seff = one
cc ........van genuchten is below
            EKR = (Seff**0.5)*(one-(one-Seff**
     1        (one/(one-one/betaim(n))))**(one-one/betaim(n)))**2
          ELSE
Cc ........brooks corey is below
            EKR = SEFF ** BROOKIM(N)
          ENDIF
        ENDIF
C------------------------------------
      ELSEIF(METHOD.EQ.5) THEN
C-------FOR TABULAR INPUT FIND Kr FROM TABLE
        BBOT = BOT(N)
        TTOP = TOP(N)
        pc = 0.5*(ttop+bbot) - hd
        if(pc.le.0)then
          EKR = 1.0
        else
          izon = iuzontabIM(n)
          iv = 3  !variable is Kr
          EKR = tabS(pc,retcrvs,izon,iv,nutabrows,nuzones,ttop)
        endif
      ELSEIF(METHOD.EQ.6) THEN
C---------FOR USE WITH TABULAR INPUT OF INTERPOLATED TABLE
        iv = 3  !variable is Kr
        izon = iuzontabIM(n)
        bott = bbot   ! bott changes to water table if below water table in tabS
        topS = tabS((ttop-hd),retcrvs,izon,iv,nutabrows,nuzones,ttop)
        botS = tabS((bbot-hd),retcrvs,izon,iv,nutabrows,nuzones,bott)
        if(hd.lt. bbot) then
c----------water level is below bottom so integrate curve from top to bottom
          if((ttop-bott).gt.1.0e-5)then
            cellS=-(topS - botS)/(ttop - bott)
            if(cellS.lt.1.0e-20) cellS = 1.0e-20
          else ! past bottom of table so use last two entries to compute value
           mul = retcrvs(iv,nutabrows,izon)-retcrvs(iv,nutabrows-1,izon)
           div = retcrvs(1,nutabrows,izon)-retcrvs(1,nutabrows-1,izon)
           cellS = mul/div
          endif
        elseif(hd.gt. ttop) then
c----------water level is above top so integrated saturation of cell is 1
          cellS = 1.0
        else
c----------water level is within cell so integrate
          cellS = (botS - topS + 1.0 * (hd -bbot))/(ttop - bbot)
        endif
        EKR = cellS
      ENDIF
C----------------------------------------------------------------
      RETURN
      END

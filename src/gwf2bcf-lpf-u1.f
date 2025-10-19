      SUBROUTINE GWF2BCFU1AR(INBCF,INLAK,INLPF,INTVM)
C     ******************************************************************
C     INITIALIZE VARIABLES AND READ DATA IN BCF OR LPF FORMATS
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,  ONLY:NLAY,IOUT,IWADI
      USE GWFBCFMODULE,ONLY:IBCFCB,IWDFLG,IWETIT,IHDWET,IKVFLAG,IKCFLAG,
     *  WETFCT,HDRY,LAYCON,LAYAVG,LAYWET,IHANISO,ISFAC,IBPN,LAYAVGV,
     *  IDRY,ITABRICH,INTRICH,IALTSTO
C     ------------------------------------------------------------------
C1------ALLOCATE SCALAR VARIABLES IN FORTRAN MODULE.
      ALLOCATE(IBCFCB,IWDFLG,IWETIT,IHDWET,IKVFLAG,IKCFLAG,IBPN,IDRY,
     *  ITABRICH,INTRICH)
      ALLOCATE(WETFCT,HDRY,IHANISO,ISFAC,IALTSTO)
      ALLOCATE(LAYCON(NLAY),LAYAVG(NLAY),laywet(nlay),LAYAVGV(NLAY))
      IALTSTO = 0 
      IHANISO = 0
      IKCFLAG = 0
      IBPN = 0
      IDRY = 0
      INTRICH = 0
      ITABRICH = 0
      IWADI = 1
      ISFAC = 1 ! STORAGE COEFFICIENT READ WHEN BCF IS USED
      DO K=1,NLAY
        LAYAVGV(K) = 0  ! FINITE DIFFERENCE HARMONIC AVERAGING IN VERTICAL DIRECTION DEFAULT
      ENDDO
C
C2------READ BCF INFORMATION USING BCF OR LPF FORMATS
      IF(INBCF.NE.0) THEN
        CALL READBCF(INBCF,INLAK,INTVM)
      ELSEIF(INLPF.NE.0) THEN
        CALL READLPF(INLPF,INLAK)
        INBCF = INLPF
        IKVFLAG = 1
      ELSE
        WRITE(IOUT,*)'BCF OR LPF SHOULD BE ON, STOPPING.'
        STOP
      ENDIF
C3------RETURN
      RETURN
      END
C -------------------------------------------------------------------------------------
      SUBROUTINE READBCF(IN,INLAK,INTVM)
C     ******************************************************************
C     ALLOCATE ARRAYS AND READ DATA FOR BLOCK-CENTERED FLOW PACKAGE
C     FILL PGF FOR FLOW TERMS AND DEALLOCATE UNNECESSARY ARRAYS
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,  ONLY:IOUT,NCOL,NROW,NLAY,ITRSS,LAYHDT,LAYHDS,NODLAY,
     1             IFREFM,IUNSTR,NODES,AREA,CL1,CL2,HNEW,So,Sn,NJA,NJAS,
     2             IDEALLOC_HY,INCLN,ISSFLG,ARAD,itrnsp,
     3             IBOUND,TOP,BOT,iconcv,iunsat,INSGB,IOUTNORMAL
      USE CLN1MODULE, ONLY: ACLNNDS,NCLNNDS
      USE GWFBCFMODULE,ONLY:IBCFCB,IWDFLG,IWETIT,IHDWET,WETFCT,HDRY,CV,
     1   HK,LAYCON,LAYAVG,SC1,SC2,WETDRY,CHANI,IHANISO,
     2   IKCFLAG,IKVFLAG,LAYWET,alpha,beta,sr,brook,BP,
     3   IBPN,LAYAVGV,IDRY,
     4   ITABRICH,INTRICH,IUZONTAB,RETCRVS,NUTABROWS,NUZONES 
C
      DOUBLE PRECISION HD,THCK,TOTTHICK,BBOT,TTOP
C
      CHARACTER*24 ANAME(2)
      DATA ANAME(1) /'COLUMN TO ROW ANISOTROPY'/
      DATA ANAME(2) /'UNSAT PARAMETER ZONE MAP'/
      CHARACTER*12 AVGNAM(4)
      CHARACTER*400 LINE
      DATA AVGNAM/'HARMONIC    ','ARITHMETIC  ',
     1            'LOGARITHMIC ','*UNCONFINED*'/
C
C     ------------------------------------------------------------------
C
C1------IDENTIFY PACKAGE
      WRITE(IOUT,1) IN
    1 FORMAT(1X,/1X,'BCF -- BLOCK-CENTERED FLOW PACKAGE, VERSION 7',
     1', 5/2/2005',/,9X,'INPUT READ FROM UNIT',I3)
C
C2------DETERMINE ISS FROM ITRSS
      IF(ITRSS.EQ.0) THEN
         ISS=1
      ELSE
         ISS=0
      END IF
C3------READ AND PRINT IBCFCB (FLAG FOR PRINTING
C3------OR UNIT# FOR RECORDING CELL-BY-CELL FLOW TERMS), HDRY
C3------(HEAD AT CELLS THAT CONVERT TO DRY), AND WETTING PARAMETERS.
      IF(IFREFM.EQ.0) THEN
        IF(IUNSTR.NE.0)THEN
          READ(IN,2)IBCFCB,HDRY,IWDFLG,WETFCT,IWETIT,IHDWET,IKVFLAG,
     *      IKCFLAG
            LLOC = 81
        ELSE
          READ(IN,2)IBCFCB,HDRY,IWDFLG,WETFCT,IWETIT,IHDWET,IKVFLAG
          LLOC = 71
        ENDIF
      ELSE
        CALL URDCOM(IN,IOUT,LINE)
        LLOC=1
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IBCFCB,R,IOUT,INDIS)
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,I,HDRY,IOUT,INDIS)
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IWDFLG,R,IOUT,INDIS)
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,I,WETFCT,IOUT,INDIS)
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IWETIT,R,IOUT,INDIS)
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IHDWET,R,IOUT,INDIS)
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IKVFLAG,R,IOUT,INDIS)
        IF(IUNSTR.NE.0)
     1   CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IKCFLAG,R,IOUT,INDIS)
      ENDIF
2     FORMAT(I10,F10.3,I10,F10.3,4I10)
C-----SEARCH FOR OPTIONS 
   20 CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,I,R,IOUT,IN)
      IF(LINE(ISTART:ISTOP).EQ.'BUBBLEPT') THEN
         IBPN=1
         WRITE(IOUT,30)
   30    FORMAT(1X,'BUBBLEPT  OPTION:',/,1X,
     1    'For Richards equation the bubble point head is also ',
     2        'input at all cells')
         allocate(bp(nodes))
         DO N=1,NODES
           BP(NODES) = 0.0
         ENDDO
      ELSE IF(LINE(ISTART:ISTOP).EQ.'FULLYDRY') THEN
         IDRY=1
         WRITE(IOUT,32)
   32    FORMAT(1X,'FULLYDRY  OPTION:',/,1X,
     1    'For Richards equation, residual saturation is only applied ',
     2    'to relative permeability and soil saturation can vary from ',
     3    '0 to 1')    
      ELSE IF(LINE(ISTART:ISTOP).EQ.'TABRICH') THEN
         ITABRICH=1
         ALLOCATE(NUZONES,NUTABROWS)
         WRITE(IOUT,33)
   33    FORMAT(1X,'TABULAR INPUT OPTION:',/,1X,
     1    'For Richards equation the retention ',
     2    'and relative permeability curves are provided ',
     3    'as tabular input for different soil types') 
         CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NUZONES,R,IOUT,IN)
         WRITE(IOUT,136) NUZONES
136      FORMAT(1X,'NUMBER OF SOIL TYPE ZONES = ', I10)
         CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NUTABROWS,R,IOUT,IN)
         WRITE(IOUT,137) NUTABROWS
137      FORMAT(1X,'NUMBER OF ROWS OF TABULAR INPUT = ', I10)         
      ELSE IF(LINE(ISTART:ISTOP).EQ.'INTRICH') THEN
         INTRICH=1
         WRITE(IOUT,34)
   34    FORMAT(1X,'INTEGRATED CURVE OPTION:',/,1X,
     1    'For Richards equation with tabular input, the retention ',
     2    'and relative permeability curves are integrated over the ',
     3    'grid block thickness') 
      ELSE IF(LINE(ISTART:ISTOP).EQ.'ALTSTO') THEN
         IALTSTO=1
         WRITE(IOUT,235)
235      FORMAT(1X,'ALTERNATIVE STORAGE FORMULATION (IALTSTO = 1)') 
      END IF
C            
      IF(LLOC.LT.200) GO TO 20
C
C-------READ STUFF TABULAR INPUT OF RICHARDS EQUATION   
      IF(ITABRICH.NE.0) THEN
C-------DIMENSION AND READ ZONE MAP          
        ALLOCATE(IUZONTAB(NODES))
        CALL U1DINT(IUZONTAB,ANAME(2),NODES,0,IN,IOUT)  
C-------DIMENSION AND READ TABLES        
        ALLOCATE(RETCRVS(3,NUTABROWS,NUZONES))
C      
        DO IZON = 1,NUZONES
          WRITE(IOUT,65) IZON  
65        FORMAT(/10X,'RETENTION TABLE NUMBER',I5/10X,27('-')/
     *  10X,'CAPILLARY HEAD',7X,'SATURATION',3X,'RELATIVE PERMEABILITY')
          DO ITROWS = 1,NUTABROWS  
            READ (IN,*) (RETCRVS(I,ITROWS,IZON),I=1,3)  
            WRITE(IOUT,66)(RETCRVS(I,ITROWS,IZON),I=1,3) 
66          FORMAT(10X,E14.6,3X,E14.6,7X,E14.6) 
          ENDDO  
        ENDDO    
      ENDIF   
C
C-------SET LAYWET FOR EACH LAYER AS PER IWDFLG
      DO K=1,NLAY
        LAYWET(K) = IWDFLG
      ENDDO
C
C3B-----PRINT VALUES
      IF(ISS.EQ.0) WRITE(IOUT,3)
    3 FORMAT(1X,'TRANSIENT SIMULATION')
      IF(ISS.NE.0) WRITE(IOUT,4)
    4 FORMAT(1X,'STEADY-STATE SIMULATION')
      IF(IBCFCB.LT.0) WRITE(IOUT,8)
    8 FORMAT(1X,'CONSTANT-HEAD CELL-BY-CELL FLOWS WILL BE PRINTED',
     1     ' WHEN ICBCFL IS NOT 0')
      IF(IBCFCB.GT.0) WRITE(IOUT,9) IBCFCB
    9 FORMAT(1X,'CELL-BY-CELL FLOWS WILL BE SAVED ON UNIT',I3)
      WRITE(IOUT,11) HDRY
   11 FORMAT(1X,'HEAD AT CELLS THAT CONVERT TO DRY=',G13.5)
      IF(IWDFLG.NE.0) GO TO 35
      WRITE(IOUT,12)
   12 FORMAT(1X,'WETTING CAPABILITY IS NOT ACTIVE')
      GO TO 50
C
   35 WRITE(IOUT,36)
   36 FORMAT(1X,'WETTING CAPABILITY IS ACTIVE')
      IF(IWETIT.LE.0) IWETIT=1
      WRITE(IOUT,37)WETFCT,IWETIT
   37 FORMAT(1X,'WETTING FACTOR=',F10.5,
     1     '     WETTING ITERATION INTERVAL=',I4)
      WRITE(IOUT,38)IHDWET
   38 FORMAT(1X,'FLAG THAT SPECIFIES THE EQUATION TO USE FOR HEAD',
     1    ' AT WETTED CELLS=',I4)
   50 CONTINUE
C
      IF(IKVFLAG.NE.0)WRITE(IOUT,42)
   42 FORMAT(1X,'IKVFLAG=1, NODAL INPUT OF KV')
      IF(IKVFLAG.EQ.0)WRITE(IOUT,43)
   43 FORMAT(1X,'IKVFLAG=0, NODAL INPUT OF CV')
C
      IF(IUNSTR.NE.0)THEN
        IF(IKCFLAG.EQ.0)WRITE(IOUT,39)
   39   FORMAT(1X,'IKCFLAG=0, NODAL INPUT OF HY AND CV')
        IF(IKCFLAG.EQ.1)WRITE(IOUT,41)
   41   FORMAT(1X,'IKCFLAG=1, CONNECTIVITY INPUT OF HY ',1X,
     1         '(OR TRAN FOR CONFINED) AND CV')
        IF(IKCFLAG.EQ.-1)WRITE(IOUT,44)
   44   FORMAT(1X,'IKCFLAG=-1, CONNECTIVITY INPUT OF CONDUCTANCE')
      ENDIF
CC
C4------READ LAYCON & PRINT TITLE FOR LAYCON TABLE.
      IF(IFREFM.EQ.0) THEN
        READ(IN,51) (LAYCON(I),I=1,NLAY)
      ELSE
        READ(IN,*) (LAYCON(I),I=1,NLAY)
      ENDIF
C
   51 FORMAT(40I2)
      WRITE(IOUT,52)
   52 FORMAT(1X,5X,'LAYER  LAYER-TYPE CODE     INTERBLOCK T',
     1      /1X,5X,44('-'))
C
C5------LOOP THROUGH LAYERS CALCULATING LAYAVG, PRINTING THE LAYER-TYPE
C5------CODE, AND COUNTING LAYERS THAT NEED TOP & BOT ARRAYS.
      NBOT=0
      NTOP=0
      ILAYUNSAT = 0 ! NO LAYERS ARE UNSATURATED
C      
      DO 100 I=1,NLAY
      IF(LAYCON(I).EQ.30 .OR. LAYCON(I).EQ.32) LAYCON(I)=LAYCON(I)-10
      INAM=LAYCON(I)/10
      LAYAVG(I)=INAM*10
      IF(LAYAVG(I).LT.0 .OR. LAYAVG(I).GT.30) THEN
         WRITE(IOUT,53) LAYAVG(I)
   53    FORMAT(1X,'INVALID INTERBLOCK T CODE:',I4)
         STOP
      END IF
      LAYCON(I)=LAYCON(I)-LAYAVG(I)
      L=LAYCON(I)
      INAM=INAM+1
      WRITE(IOUT,55) I,L,LAYAVG(I),AVGNAM(INAM)
   55 FORMAT(1X,I9,I13,I11,' -- ',A)
C5A-----CONVERT LAYAVG TO CONFORM TO THAT IN LPF PACKAGE
      LAYAVG(I) = LAYAVG(I)/10
      IF(LAYAVG(I).EQ.1)THEN
          LAYAVG(I)=3
      ELSEIF(LAYAVG(I).EQ.2)THEN
          LAYAVG(I)=1
      ELSEIF(LAYAVG(I).EQ.3)THEN
          LAYAVG(I)=2
      ENDIF
C5A2 ---ALSO CONVERT LAVAVGV TO 1 (WEIGHTED ARITH AVERAGING FOR LPF)
C5A2----WHEN LAYAVG IS 4 FOR FINITE ELEMENT AVERAGING
      IF(LAYAVG(I).EQ.4)THEN
        LAYAVGV(I) = 1
      ENDIF
C5B-----NEW LAYCON OF 4 IS FOR UNCONFINED T, CONVERTABLE S
C5B-----WITH T COMPUTED USING UPSTREAM WATER TABLE DEPTH
      IF(LAYCON(I).LT.0 .OR. LAYCON(I).GT.5) THEN
         WRITE(IOUT,56) LAYCON(I)
   56    FORMAT(1X,'INVALID LAYER TYPE:',I4)
         STOP
      END IF
C
      IF(LAYCON(I).EQ.5) ILAYUNSAT = 1 ! AT LEAST 1 LAYER USES RICHARDS EQUATION    
C
C5C-----SET GLOBAL HEAD-DEPENDENT THICKNESS FLAGS.
      IF (L.EQ.0) THEN
        LAYHDT(I)=0
        LAYHDS(I)=0
      ELSEIF (L.EQ.1) THEN
        LAYHDT(I)=1
        LAYHDS(I)=0
      ELSEIF (L.EQ.2) THEN
        LAYHDT(I)=0
        LAYHDS(I)=1
      ELSE
        LAYHDT(I)=1
        LAYHDS(I)=1
      ENDIF
C
C5D-----ONLY THE TOP LAYER CAN BE UNCONFINED(LAYCON=1).
      IF(L.NE.1 .OR. I.EQ.1) GO TO 70
      WRITE(IOUT,57)
   57 FORMAT(1X,/1X,'LAYER TYPE 1 IS ONLY ALLOWED IN TOP LAYER')
      STOP
C
C5E-----LAYER TYPES 1, 4, AND 3 NEED A BOTTOM. ADD 1 TO KB.
   70 IF(L.EQ.1 .OR. L.EQ.3 .OR. L.EQ.4) NBOT=NBOT+1
C
C5F-----LAYER TYPES 2, 4, AND 3 NEED A TOP. ADD 1 TO KT.
      IF(L.EQ.2 .OR. L.EQ.3 .OR. L.EQ.4) NTOP=NTOP+1
  100 CONTINUE
C ------------------------------------------------------------------------
C6------ALLOCATE SPACE FOR ARRAYS.
      IF(ISS.EQ.0) THEN
         ALLOCATE(SC1(NODES))
      ELSE
         ALLOCATE(SC1(1))
      END IF
      IF(NTOP.GT.0 .AND. ISS.EQ.0 .or. itrnsp.ne.0) THEN
         ALLOCATE(SC2(NODES))
         SC2 = 0.0
      ELSE
         ALLOCATE(SC2(1))
      END IF
      ALLOCATE(CHANI(NLAY))
      ALLOCATE(HK(NODES))
      ALLOCATE(CV(NODES))
      IF(IWDFLG.NE.0 .AND. NBOT.GT.0) THEN
         ALLOCATE(WETDRY(NODES))
      ELSE
         ALLOCATE(WETDRY(1))
      END IF
C      
      if(iunsat.eq.1.OR.ILAYUNSAT.EQ.1)then
        if(itabrich.eq.0)then  
          ALLOCATE(alpha(NODES),beta(nodes),sr(nodes),brook(nodes))
          DO N=1,NODES
            ALPHA(N) =0.0
            BETA(N) = 1.0
            SR(N) = 0.999
            BROOK(N) = 1.0
          ENDDO
        endif  
        IF(IUNSAT.EQ.1)THEN
          DO K=1,NLAY
            IF(LAYCON(K).NE.0) LAYCON(K)=5
          ENDDO 
        ENDIF  
      ENDIF
C------------------------------------------------------------------------------------
C7------READ ANISOTROPIES FOR EACH LAYER IF NODAL K-VALUES ARE READ
      IF(IKCFLAG.EQ.0)THEN
        K = 0
        CALL U1DREL(CHANI(1),ANAME(1),NLAY,K,IN,IOUT)
C7A-------SET HORIZONTAL ANISOTROPY FLAG
        IHANISO = 0
        DO K=1,NLAY
          IF(ABS(CHANI(K) - 1.0).GT.1.0E-6) IHANISO = 1
        ENDDO
      ENDIF
      IF(IHANISO.EQ.1)THEN 
       ALLOCATE(IOUTNORMAL)
       IOUTNORMAL = 0
       ALLOCATE(ARAD(NJAS))
       ENDIF
C---------------------------------------------------------------------------------
C8------READ PARAMETERS AND CONVERT FOR UNSTRUCTURED AND STRUCTURED GRIDS
      IF(IUNSTR.EQ.0) THEN
        CALL SGWF2BCFU1S(IN)
      ELSE
        CALL SGWF2BCFU1G(IN)
      ENDIF
C8A-----FOR LAYERS WITH LAYCON=0 RESET TRANSMISSIVITIES TO HK VALUES
      DO K=1,NLAY
        IF(LAYCON(K).EQ.0.OR.LAYCON(K).EQ.2)THEN
          NNDLAY = NODLAY(K)
          NSTRT = NODLAY(K-1)+1
          DO N=NSTRT,NNDLAY
            IF(IBOUND(N).NE.0) THEN
CB-------------CALCULATE HK FROM TRANSMISSIVITY.
              HD=HNEW(N)
              BBOT=BOT(N)
              TTOP=TOP(N)
              TOTTHICK = TTOP - BBOT
              HK(N) = HK(N) / TOTTHICK
            ENDIF
          ENDDO
        ENDIF
      ENDDO
C--------------------------------------------------------------------------------
C9------SET CONSTANT TERMS IN PGF ARRAY IF IT IS NOT READ DIRECTLY
C-------FOR NODAL INPUT OF CONDUCTIVITIES
      IF(IKCFLAG.EQ.0) THEN
C
C9A------FILL VERTICAL DIRECTION IN PGF
        IF(NLAY.GT.1)THEN
          DO K=1,NLAY
            CALL SGWF2BCFU1VCONDV(K)
          ENDDO
        ENDIF
C
C9B-------FILL HORIZONTAL DIRECTION IN PGF  - HY FOR LAYCON 4 AND T FOR LAYCONS 0 OR 2
        CALL FILLPGFH
      ENDIF
C--------------------------------------------------------------------------------
C10------SET INITIAL SATURATIONS.
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
C
C--------------------------------------------------------------------------------
C11------IF TRANSIENT, CALCULATE STORAGE CAPACITY.
      IF(ISS.NE.0) GO TO 200
C
C11A-----MULTIPLY PRIMARY STORAGE COEFFICIENT BY AREA TO GET
C11A-----PRIMARY STORAGE CAPACITY.
      DO 80 N=1,NODES
      SC1(N)=SC1(N)*AREA(N)
   80 CONTINUE
C
C11B-----IF LAYER IS CONF/UNCONF MULTIPLY SECONDARY STORAGE COEFFICIENT
C11B-----BY AREA TO GET SECONDARY STORAGE CAPACITY(SC2).
      DO K=1,NLAY
        IF(LAYCON(K).NE.3.AND.LAYCON(K).NE.2.AND.LAYCON(K).NE.4
     1   .AND.LAYCON(K).NE.5)GO TO 90
        NNDLAY = NODLAY(K)
        NSTRT = NODLAY(K-1)+1
        DO 85 N=NSTRT,NNDLAY
          SC2(N)=SC2(N)*AREA(N)
   85   CONTINUE
   90   CONTINUE
      ENDDO
C
  200 CONTINUE
C
C--------------------------------------------------------------------------------
C12------DEALLOCATE UNWANTED ARRAYS
csp need cv for merging Kv with boundary leakance for RIV      DEALLOCATE(CV)
      IDEALLOC_HY = 0
C
      ILAYCON13=0
      DO I=1,NLAY
        IF(LAYCON(I).EQ.1.OR.LAYCON(I).EQ.3)ILAYCON13=1
      ENDDO
      IF(ILAYCON13.EQ.0) THEN
        IDEALLOC_HY = 1
        IF(INCLN.GT.0) IDEALLOC_HY = 2  !DEALLOCATE ONLY AFTER CLN1AR IS DONE
        IF(INSGB.GT.0) IDEALLOC_HY = 0  !DO NOT DEALLOCATE IF SGB IS ACTIVATED
      ENDIF
      IF(INLAK.NE.0) IDEALLOC_HY = 0    ! LAKE PACKAGE USES HY SO DO NOT DEALLOCATE
      IF(INTVM.NE.0) IDEALLOC_HY = 0    ! TVM PACKAGE USES HY SO DO NOT DEALLOCATE
C
      IF(ILAYCON13.EQ.0)THEN
CC        DEALLOCATE(CHANI)
CC--CHECK THIS, IT MAY BE NEEDED FOR TRANSPORT
CC        IF(IUNSTR.EQ.1)THEN
CC          DEALLOCATE(CL1,CL2)
CC        ENDIF
      ENDIF
C
      IF(IDEALLOC_HY.EQ.1)  DEALLOCATE(HK)
C
C13-----RETURN
      RETURN
      END
C
      SUBROUTINE SGWF2BCFU1S(IN)
C     ******************************************************************
C-----READ PARAMETERS FOR STRUCTURED GRID AND CONVERT TO UNSTRUCTURED
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,  ONLY:IOUT,NCOL,NROW,NLAY,ITRSS,LAYHDT,LAYHDS,ARAD,
     1  IPRCONN,IFREFM,IUNSTR,NODES,NJA,NJAS,IA,JA,JAS,ISYM,iunsat
      USE GWFBCFMODULE,ONLY:IBCFCB,IWDFLG,IWETIT,IHDWET,WETFCT,HDRY,CV,
     1    LAYCON,LAYAVG,HK,SC1,SC2,WETDRY,CHANI,
     1    IKVFLAG,IHANISO,laywet,alpha,beta,sr,brook,BP,IBPN,ITABRICH
      USE GWTBCTMODULE, ONLY: CBCF
C
      CHARACTER*24 ANAME(11)
      REAL, DIMENSION(:,:),ALLOCATABLE  ::TEMP
      REAL, DIMENSION (:), ALLOCATABLE :: TEMPPL
      DATA ANAME(1) /'    PRIMARY STORAGE COEF'/
      DATA ANAME(2) /'    TRANSMIS. ALONG ROWS'/
      DATA ANAME(3) /'   HYD. COND. ALONG ROWS'/
      DATA ANAME(4) /'VERT HYD COND /THICKNESS'/
      DATA ANAME(5) /'  SECONDARY STORAGE COEF'/
      DATA ANAME(6) /'        WETDRY PARAMETER'/
      DATA ANAME(7) /'                   alpha'/
      DATA ANAME(8) /'                    beta'/
      DATA ANAME(9) /'                      sr'/
      DATA ANAME(10) /'                   brook'/
      DATA ANAME(11) /'     BUBBLING POINT HEAD'/
      REAL PI
C     ------------------------------------------------------------------
C1------ALLOCATE TEMP ARRAY FOR STORING 3-D INFORMATION
      ALLOCATE(TEMP(NCOL,NROW))
C--------------------------------------------------------------------
C
C2------SET ANGLE INTO ARAD IF ANISOTROPIC
      IF(IHANISO.EQ.1)THEN
        PI = 3.1415926536
C2A----SET FACE ANGLES IN ARAD
        DO N=1,NODES
          DO II = IA(N)+1,IA(N+1)-1
            JJ = JA(II)
            IF(JJ.GE.N) CYCLE
            IIS = JAS(II)
            IF((N-JJ).EQ.1) THEN
              ARAD(IIS) = 0.0
            ELSE
              ARAD(IIS) = PI/2.0
            ENDIF
          ENDDO
        ENDDO
C
C2B-------WRITE FACE ANGLES ARRAY
        IF(IPRCONN.NE.0)THEN
          WRITE(IOUT,*)'FACE ANGLE IS BELOW, 22G15.6, UNSYMMETRIC'
          ALLOCATE(TEMPPL(NJA))
          DO N=1,NODES
          DO II = IA(N)+1,IA(N+1)-1
            JJ = JA(II)
            IF(JJ.GE.N)THEN
              IIS = JAS(II)
              TEMPPL(II) = ARAD(IIS)
              TEMPPL(ISYM(II)) = ARAD(IIS)
            ENDIF
          ENDDO
          ENDDO
          WRITE(IOUT,55)(TEMPPL(J),J=1,NJA)
55      FORMAT(1P,22G15.6)
CSP          WRITE(IOUT,55)(ARAD(J),J=1,NJAS) !COMMENTED OUT SYMMETRIC WRITE
          DEALLOCATE (TEMPPL)
        ENDIF
      ENDIF
C------------------------------------------------------------------------
C
C3------READ ARRAYS FOR EACH LAYER.
      DO 300 K=1,NLAY
      KK=K
C---------------------------------------------------------
C3A-----READ PRIMARY STORAGE COEFFICIENT INTO ARRAY SC1 IF TRANSIENT.
      IF(ITRSS.NE.0)THEN
        CALL U2DREL(TEMP(1,1),ANAME(1),NROW,NCOL,KK,IN,IOUT)
        DO I=1,NROW
        DO J=1,NCOL
          N=J+(I-1)*NCOL+(K-1)*NROW*NCOL
          SC1(N) = TEMP(J,I)
        ENDDO
        ENDDO
      ENDIF
C---------------------------------------------------------
C3B-----READ TRANSMISSIVITY INTO ARRAY CC IF LAYER TYPE IS 0 OR 2.
      IF(LAYCON(K).GE.3 .OR. LAYCON(K).EQ.1) GO TO 105
      CALL U2DREL(TEMP(1,1),ANAME(2),NROW,NCOL,KK,IN,IOUT)
      GO TO 110
C
C3C-----READ HYDRAULIC CONDUCTIVITY(HY) IF LAYER TYPE IS 1 OR 3.
  105 CALL U2DREL(TEMP(1,1),ANAME(3),NROW,NCOL,KK,IN,IOUT)
C3D-----PUT TRANSMISSIVITY OR HYDRAULIC CONDUCTIVITY INTO HY ARRAY
C3D-----HY CONTAINS BOTH T AND K VALUES DEPENDING ON LAYCON OF THE LAYER
110   CONTINUE
      DO I=1,NROW
      DO J=1,NCOL
        N=J+(I-1)*NCOL+(K-1)*NROW*NCOL
        HK(N) = TEMP(J,I)
      ENDDO
      ENDDO
C-----------------------------------------------------------
C3E-----READ VCONT OR KV DEPENDING ON IKVFLAG
      IF(IKVFLAG.EQ.0)THEN
C
C3E1-----READ VERTICAL HYCOND/THICK INTO ARRAY CV IF NOT BOTTOM LAYER;
        IF(K.LT.NLAY) THEN
          CALL U2DREL(TEMP(1,1),ANAME(4),NROW,NCOL,KK,IN,IOUT)
          DO I=1,NROW
          DO J=1,NCOL
            N=J+(I-1)*NCOL+(K-1)*NROW*NCOL
            CV(N) = TEMP(J,I)
          ENDDO
          ENDDO
        ENDIF
      ELSE IF(NLAY.GT.1)THEN
C
C3E2-----READ VERTICAL HYCOND INTO ARRAY CV FOR ALL NODES
        CALL U2DREL(TEMP(1,1),ANAME(4),NROW,NCOL,KK,IN,IOUT)
        DO I=1,NROW
        DO J=1,NCOL
          N=J+(I-1)*NCOL+(K-1)*NROW*NCOL
          CV(N) = TEMP(J,I)
        ENDDO
        ENDDO
      ENDIF
C---------------------------------------------------------
C3F-----READ SECONDARY STORAGE COEFFICIENT INTO ARRAY SC2 IF TRANSIENT
C3F-----AND LAYER TYPE IS 2 OR 3 OR 4.
  120 IF(LAYCON(K).NE.3.AND.LAYCON(K).NE.2.AND.LAYCON(K).NE.4
     1    .AND.LAYCON(K).NE.5)GOTO 130
      IF(ITRSS.NE.0)THEN
        CALL U2DREL(TEMP(1,1),ANAME(5),NROW,NCOL,KK,IN,IOUT)
        DO I=1,NROW
        DO J=1,NCOL
          N=J+(I-1)*NCOL+(K-1)*NROW*NCOL
          SC2(N) = TEMP(J,I)
        ENDDO
        ENDDO
      ENDIF
C---------------------------------------------------------
C3G-----READ WETDRY CODES IF LAYER TYPE IS 1 OR 3 AND WETTING
C3G-----CAPABILITY HAS BEEN INVOKED (IWDFLG NOT 0).
  130 IF(LAYCON(K).NE.3.AND.LAYCON(K).NE.1)GO TO 200
      IF(laywet(k).EQ.0)GO TO 200
      CALL U2DREL(TEMP(1,1),ANAME(6),NROW,NCOL,KK,IN,IOUT)
      DO I=1,NROW
      DO J=1,NCOL
        N=J+(I-1)*NCOL+(K-1)*NROW*NCOL
        WETDRY(N) = TEMP(J,I)
      ENDDO
      ENDDO
c
  200 CONTINUE
C4----read parameters for Richards Equation if solving unsaturated zone flow
      IF(LAYCON(K).NE.5)GOTO 300
      IF(ITABRICH.EQ.0) THEN
C5-----READ alpha, beta, sr, brook
      CALL U2DREL(TEMP(1,1),ANAME(7),NROW,NCOL,KK,IN,IOUT)
      DO I=1,NROW
      DO J=1,NCOL
        N=J+(I-1)*NCOL+(K-1)*NROW*NCOL
        alpha(N) = TEMP(J,I)
      ENDDO
      ENDDO
      CALL U2DREL(TEMP(1,1),ANAME(8),NROW,NCOL,KK,IN,IOUT)
      DO I=1,NROW
      DO J=1,NCOL
        N=J+(I-1)*NCOL+(K-1)*NROW*NCOL
        beta(N) = TEMP(J,I)
      ENDDO
      ENDDO
      CALL U2DREL(TEMP(1,1),ANAME(9),NROW,NCOL,KK,IN,IOUT)
      DO I=1,NROW
      DO J=1,NCOL
        N=J+(I-1)*NCOL+(K-1)*NROW*NCOL
        sr(N) = TEMP(J,I)
      ENDDO
      ENDDO
      CALL U2DREL(TEMP(1,1),ANAME(10),NROW,NCOL,KK,IN,IOUT)
      DO I=1,NROW
      DO J=1,NCOL
        N=J+(I-1)*NCOL+(K-1)*NROW*NCOL
        brook(N) = TEMP(J,I)
      ENDDO
      ENDDO
      IF(IBPN.GT.0)THEN
        CALL U2DREL(TEMP(1,1),ANAME(11),NROW,NCOL,KK,IN,IOUT)
        DO I=1,NROW
        DO J=1,NCOL
          N=J+(I-1)*NCOL+(K-1)*NROW*NCOL
          bP(N) = TEMP(J,I)
        ENDDO
        ENDDO
      ENDIF
C6------read tabular input for retention and relative permeability curves      
      ELSE


      ENDIF
  300 CONTINUE
C
C4------DEALLOCATE TEMPORARY ARRAY
      DEALLOCATE(TEMP)
C
C6------RETURN
      RETURN
      END
C----------------------------------------------------------------------------
      SUBROUTINE SGWF2BCFU1G(IN)
C     ******************************************************************
C-----READ PARAMETERS AND CONVERT FOR UNSTRUCTURED (GENERALIZED) GRID
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,   ONLY:IOUT,NCOL,NROW,NLAY,ITRSS,LAYHDT,LAYHDS,IFREFM,
     1              IUNSTR,NODES,NJA,NJAS,NJAG,IA,PGF,FAHL,ARAD,JA,JAS,
     2              NODLAY,IDSYMRD,IATMP,NJATMP,TOP,BOT,CL1,CL2,iunsat,
     3              IVC  
      USE GWFBCFMODULE,ONLY:IBCFCB,IWDFLG,IWETIT,IHDWET,WETFCT,HDRY,CV,
     1    LAYCON,LAYAVG,HK,SC1,SC2,WETDRY,IHANISO,IKCFLAG,IKVFLAG,
     1    laywet,alpha,beta,sr,brook,BP,IBPN,ITABRICH
C
      CHARACTER*24 ANAME(14)
      REAL, DIMENSION(:),ALLOCATABLE  ::TEMP
      DATA ANAME(1) /'    PRIMARY STORAGE COEF'/
      DATA ANAME(2) /'   TRANSMIS. OR HYD COND'/
      DATA ANAME(3) /'   CONNECTION HYD. COND.'/
      DATA ANAME(4) /'VERT HYD COND /THICKNESS'/
      DATA ANAME(5) /'  SECONDARY STORAGE COEF'/
      DATA ANAME(6) /'COLUMN TO ROW ANISOTROPY'/
      DATA ANAME(7) /'        WETDRY PARAMETER'/
      DATA ANAME(8) /'           VERT HYD COND'/
      DATA ANAME(9) /'              FACE ANGLE'/
      DATA ANAME(10) /'                   alpha'/
      DATA ANAME(11) /'                    beta'/
      DATA ANAME(12) /'                      sr'/
      DATA ANAME(13) /'                   brook'/  
      DATA ANAME(14) /'     BUBBLING POINT HEAD'/
C     ------------------------------------------------------------------
C
C1------CHECK CONSISTENCY OF FLAGS
C-------FLAG TO READ NODAL CONDUCTANCE OR CONNECTIVITY CONDUCTANCE
C     IKCFLAG = 0 IF CONDUCTANCE IS NODAL, =1 IF CONDUCTANCE IS ON CONNECTIVITY
      IF(IKCFLAG.NE.0.OR.IWDFLG.NE.0)THEN
        ILAYCON13=0
        DO I=1,NLAY
          IF(LAYCON(I).EQ.1.OR.LAYCON(I).EQ.3)ILAYCON13=1
        ENDDO
        IF(IKCFLAG.NE.0.AND.ILAYCON13.EQ.1)THEN
          WRITE(IOUT,10)
          STOP
        ENDIF
      ENDIF
10    FORMAT(5X,'LAYCON=1 OR 3 NOT ALLOWED WITH IKCFLAG=1 OR -1 SINCE ',
     1 'K OF CONNECTION IS NOT STORED SEPARATELY')
C
C--------------------------------------------------------------------
C2------READ FACE ANGLES IF ANISOTROPIC
      IF(IHANISO.EQ.1)THEN
        CALL U1DRELNJA(ARAD,IATMP,ANAME(9),NJATMP,IN,IOUT,IDSYMRD)
      ENDIF
C---------------------------------------------------------------------------
C
C3------LOOP OVER ALL LAYERS AND READ ARRAYS FOR NODES
      DO 200 K = 1,NLAY
      KK = K
      NNDLAY = NODLAY(K)
      NSTRT = NODLAY(K-1)+1
      NDSLAY = NNDLAY - NODLAY(K-1)
C---------------------------------------------------------
C
C3A-----READ PRIMARY STORAGE COEFFICIENT INTO ARRAY SC1 IF TRANSIENT.
      IF(ITRSS.NE.0)THEN
        CALL U1DREL(SC1(NSTRT),ANAME(1),NDSLAY,K,IN,IOUT)
      ENDIF
C-------------------------------------------------
      IF(IKCFLAG.NE.0)GO TO 120
C-------------------------------------------------
C
C3B-----READ TRANSMISSIVITY INTO ARRAY HY IF LAYER TYPE IS 0 OR 2,
C3B-----OR READ HYDRAULIC CONDUCTIVITY INTO ARRAY HY IF LAYER TYPE IS 1 OR 3.
      CALL U1DREL(HK(NSTRT),ANAME(2),NDSLAY,K,IN,IOUT)
C
        IF(IKVFLAG.EQ.0)THEN
C3C-------READ VERTICAL HYCOND/THICK INTO ARRAY CV FOR ALL EXCEPT LAST LAYER;
C3C-------MULTIPLIED BY CELL AREA TO CONVERT TO CONDUCTANCE LATER.
          IF(K.NE.NLAY) CALL U1DREL(CV(NSTRT),ANAME(4),NDSLAY,K,IN,IOUT)
        ELSE
C3D-------READ KV INTO ARRAY CV FOR FOR ALL LAYERS;
          IF(NLAY.GT.1) CALL U1DREL(CV(NSTRT),ANAME(8),NDSLAY,K,IN,IOUT)
        ENDIF
C-------------------------------------------------
120   CONTINUE
C-------------------------------------------------
C
C3E-----READ SECONDARY STORAGE COEFFICIENT INTO ARRAY SC2 IF TRANSIENT
C3E-----AND LAYER TYPE IS 2 OR 3 OR 4.
      IF((LAYCON(K).EQ.2.OR.LAYCON(K).EQ.3.OR.LAYCON(K).EQ.4) 
     1  .OR.LAYCON(K).EQ.5)THEN
        IF(ITRSS.NE.0)THEN
          CALL U1DREL(SC2(NSTRT),ANAME(5),NDSLAY,K,IN,IOUT)
        ENDIF
      ENDIF
C---------------------------------------------------------
C3F-----READ WETDRY CODES IF LAYER TYPE IS 1 OR 3 AND WETTING
C3F-----CAPABILITY HAS BEEN INVOKED (IWDFLG NOT 0).
  130 CONTINUE
      IF(laywet(k).EQ.0)GO TO 140
      IF(LAYCON(K).EQ.1.OR.LAYCON(K).EQ.3)
     * CALL U1DREL(WETDRY(NSTRT),ANAME(7),NDSLAY,K,IN,IOUT)
140   CONTINUE
C---------------------------------------------------------
C      if(iunsat.eq.0) go to 300
      IF(LAYCON(K).NE.5)GOTO 300
      IF(ITABRICH.EQ.0) THEN       
C3H-----READ alpha, beta, sr, brook
      CALL U1DREL(alpha(NSTRT),ANAME(10),NDSLAY,K,IN,IOUT)
      CALL U1DREL(beta(NSTRT),ANAME(11),NDSLAY,K,IN,IOUT)
      CALL U1DREL(sr(NSTRT),ANAME(12),NDSLAY,K,IN,IOUT)
      CALL U1DREL(brook(NSTRT),ANAME(13),NDSLAY,K,IN,IOUT)
      IF(IBPN.GT.0)THEN
        CALL U1DREL(bP(NSTRT),ANAME(14),NDSLAY,K,IN,IOUT)  
      ENDIF
C
      ENDIF      
  300 CONTINUE
C--------------------------------------------------------
200   CONTINUE   
C--------------------------------------------------------
C3G-----READ EFFECTIVE SATURATED CONDUCTIVITY OF CONNECTION
      IF(IKCFLAG.NE.0)THEN
        ALLOCATE(TEMP(NJAS))
        CALL U1DRELNJA(TEMP(1),IATMP,ANAME(3),NJATMP,IN,IOUT,IDSYMRD)
        DO N=1,NODES
          THICK1 = TOP(N) - BOT(N)
          DO II = IA(N)+1,IA(N+1)-1
            JJ = JA(II)
            IF(JJ.LE.N.OR.JJ.GT.NODES) CYCLE
            IIS = JAS(II)
            IF(IKCFLAG.EQ.1)THEN
              THICK2 = TOP(JJ) - BOT(JJ)
              THICK = 0.5 * (THICK1 + THICK2)
              PGF(IIS) = PGF(IIS) * TEMP(IIS) * THICK
            ELSE
              PGF(IIS) = TEMP(IIS)
            ENDIF
          ENDDO
        ENDDO 
C-------SET HK FOR THEIM SOLUTION CONNECTION
          DO N=1,NODES
            THICK = TOP(N) - BOT(N)
            AKN = 0.0
            IKN = 0
C-----------GO OVER CONNECTIONS OF NODE N AND FILL FOR UPPER SYMMETRIC PART
            DO II = IA(N)+1,IA(N+1)-1
              JJ = JA(II)
              IF(JJ.LE.NODES)THEN
                IIS = JAS(II)  
                IF(IVC(IIS).EQ.0)THEN
                  IF(IKCFLAG.EQ.1) THEN
                    AKN = AKN + PGF(IIS) / THICK 
                  ELSE
                    AKN = AKN + PGF(IIS) / THICK * CL1(IIS)  
                  ENDIF
                ENDIF     
                IKN = IKN + 1
              ENDIF
            ENDDO
            IF(IKN.GT.0) THEN 
              HK(N) = AKN / IKN
            ENDIF
          ENDDO
C-----------------------------------------------------                
        DEALLOCATE(TEMP)
      ENDIF
C
C4------RETURN
      RETURN
      END
C -----------------------------------------------------------------------
      SUBROUTINE FILLPGFH
C     ******************************************************************
C     FILL CONSTANT TERMS INTO PGF - HK FOR LAYCON 4 AND T FOR LAYCONS 0 OR 2
C     AND CV FOR VERTICAL CONNECTION FOR AN UNSTRUCTURED GRID INPUT.
C     ******************************************************************
      USE GLOBAL, ONLY:CL1,CL2,TOP,BOT,NODES,NLAY,NROW,NCOL,NJA,IA,PGF,
     1            NODLAY,IBOUND,JA,JAS,IVC,ISYM,ICONCV,AREA
      USE GWFBCFMODULE,ONLY:HK,CV,LAYCON,LAYAVG,IKVFLAG
      DOUBLE PRECISION ANUM,THIK1,THIK2,BNUM
C----------------------------------------------------------
C1------loop over all nodes
C1A-----loop over all layers
      DO K=1,NLAY
        NNDLAY = NODLAY(K)
        NSTRT = NODLAY(K-1)+1
C1B-------loop over all nodes within each layer
        DO N=NSTRT,NNDLAY
C2----------loop over all connections of node N and fill upper triangle with PGF term
          DO II = IA(N)+1,IA(N+1)-1
            JJ = JA(II)
C3------------fill only for upper triangle of porous medium nodes
            IF(JJ.LE.N.OR.JJ.GT.NODES) CYCLE
            IIS = JAS(II)
C
C5--------------for horizontal direction connection
              IF(IVC(IIS).NE.1)THEN !HORIZONTAL DIRECTION CONNECTION
              IF(LAYCON(K).EQ.1.OR.LAYCON(K).EQ.3)CYCLE !K portion of PGF included in FM
              EL1 = CL1(IIS)
              EL2 = CL2(IIS)
              THIK1 = 1.0
              THIK2 = 1.0
              IF(LAYCON(K).EQ.4.OR.LAYCON(K).EQ.5)THEN
                IF(IVC(IIS).EQ.2)THEN
                  BNUM = MIN( (TOP(N) - BOT(N)) , (TOP(JJ) - BOT(JJ)))
                ELSE
                  BNUM = 0.5*( (TOP(N) - BOT(N)) + (TOP(JJ) - BOT(JJ)))
                ENDIF
              ELSE
                BNUM = 1.0
                THIK1 = TOP(N) - BOT(N)
                THIK2 = TOP(JJ) - BOT(JJ)
              ENDIF
C5A-------------compute intercell conductance
              CALL CBCK12(II,N,JJ,K,EL1,EL2,ANUM,HK,
     *          THIK1,THIK2)
C5B-------------fill pgf with conductance term
              PGF(IIS) = PGF(IIS) * ANUM*BNUM
            ENDIF
          ENDDO
C
        ENDDO
      ENDDO
C6------return
      RETURN
      END
C------------------------------------------------------------------------
      SUBROUTINE GWF2BCFU1FM(KITER,KSTP,KPER)
C     ******************************************************************
C     ADD LEAKAGE CORRECTION AND STORAGE TO AMAT AND RHS, AND CALCULATE
C     CONDUCTANCE AS REQUIRED
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,  ONLY:IOUT,NCOL,NROW,NLAY,IBOUND,TOP,NODLAY,IVC,ISYM,
     1    HNEW,AMAT,RHS,HOLD,ISSFLG,IA,JA,JAS,PGF,NJA,NJAS,NODES,NOVFC,
     2    ICONCV,So,Sn,BOT,TOP,NEQS,INCLN,INGNC,INGNC2,INGNCn,
     3    AKR,AKRC,ILAYCON4,IWADI, iunsat
      USE CLN1MODULE, ONLY: NCLNNDS,ACLNNDS
      USE GWFBASMODULE,ONLY:DELT
      USE SMSMODULE, ONLY: IBFLAG
      USE GWFBCFMODULE,ONLY:LAYCON,SC1,SC2,CV,HWADIGW,IALTSTO
      DOUBLE PRECISION HTMP,TLED,RHO,RHO1,RHO2,TP,SOLD,SNEW,THCK,DS,HD,
     *          EX,FTERM,DFTERM,TOTTHICK,BBOT,TTOP,SATN,SATO,X,Y,EPS
C     ------------------------------------------------------------------
C1------SET STEADY-STATE FLAG
      ISS=ISSFLG(KPER)
      KB=0
      KT=0
C
C-----------------------------------------------------------------------
C2A-------FOR EACH LAYER COMPUTE NONLINEAR NODAL TRANSMISSIVITY TERMS
      DO K=1,NLAY
        KK = K
        IF(LAYCON(K).EQ.3. OR.LAYCON(K).EQ.1) THEN
C2A-------SORT OUT WET AND DRY NODES
C2A-------AND COMPUTE NODAL TRANSMISSIVITY IN BUFF
          CALL SGWF2BCFU1WDS(KK,KITER,KSTP,KPER)
        ELSEIF(LAYCON(K).EQ.4. OR. LAYCON(K).EQ.5)THEN
C-----------------------------------------------------------------------
C2B------FOR EACH LAYER, CALCULATE Sn IF LAYCON IS 4 OR 5
          CALL SGWF2BCFU1WDS4(K)
        ENDIF
      ENDDO
C2C------ALLOCATE SPACE FOR KR IF LAYCON IS 4 OR IF CLN PACKAGE IS ON
      IF(ILAYCON4.EQ.1.OR.INCLN.GT.0)THEN
        IF(ALLOCATED(AKRC)) DEALLOCATE (AKRC)  
        IF(ALLOCATED(AKR)) DEALLOCATE (AKR) 
        ALLOCATE(AKRC(NJAS),AKR(NEQS))
        AKRC = 1.0
        AKR = 1.0
        CALL SGWF2FILLKRS
      ENDIF
C
C-----------------------------------------------------------------------
C3------FOR EACH LAYER: IF T VARIES CALCULATE HORIZONTAL CONDUCTANCES
      DO 100 K=1,NLAY
        KK=K
C
        IF(LAYCON(K).EQ.3. OR.LAYCON(K).EQ.1) THEN
C
C3A-------FOR LAYER TYPES 1, & 3 CALL SGWF2BCFU1H TO CALCULATE
C3A-------HORIZONTAL AND VERTICAL CONDUCTANCES IN AMAT.
          CALL SGWF2BCFU1H(KK,KITER)
        ELSEIF(LAYCON(K).EQ.4.OR.LAYCON(K).EQ.5)THEN
C
C3B-------FOR LAYER TYPE 4 CALL SGWF2BCFU1H TO CALCULATE
C3B-------HORIZONTAL AND VERTICAL CONDUCTANCES IN AMAT.
          CALL SGWF2BCFU1H4(KK)
        ELSE
C3C-------TRANSFER PGF INTO AMAT FOR LAYCON OF 0 AND 2
          NNDLAY = NODLAY(K)
          NSTRT = NODLAY(K-1)+1
          DO N=NSTRT,NNDLAY
            IF(IBOUND(N).NE.0)THEN
C-------------GO OVER CONNECTIONS OF NODE N AND FILL FOR UPPER AND SYMMETRIC PART
              DO II = IA(N)+1,IA(N+1)-1
                JJ = JA(II)
                IF(JJ.GE.N.AND.IBOUND(JJ).NE.0.AND.JJ.LE.NODES)THEN
                  IIS = JAS(II)
                  AMAT(II) = PGF(IIS)
                  AMAT(ISYM(II)) = PGF(IIS)
                ENDIF
              ENDDO
            ENDIF
          ENDDO
        ENDIF
C
  100 CONTINUE
C
C-----------------------------------------------------------------------
C4------ACCOMMODATE LPF's ICONCV=0 AND NOCVCO=0 - ICONCV IS 0 ONLY FOR LPF
      IF(NLAY.GT.1.AND.ICONCV.EQ.0)THEN
        DO 101 K=1,NLAY
          IF(K.NE.NLAY)THEN
            IF(LAYCON(K).EQ.3. OR.LAYCON(K).EQ.1. OR.
     *       LAYCON(K+1).EQ.3. OR.LAYCON(K+1).EQ.1)
     *       CALL SGWF2LPFU1VCONDV(K)
          ENDIF
  101   CONTINUE
      ENDIF
C
C5------CALCULATE PIVOT OF FLOW TERM IN AMAT
      DO N=1,NEQS
        ID = IA(N)
c        AMAT(ID) = 0.0
        DO II = IA(N)+1,IA(N+1)-1
          AMAT(ID) = AMAT(ID) - AMAT(II)
        ENDDO
      ENDDO
C
C----------------------------------------------------------------------
C6------IF THE SIMULATION IS TRANSIENT ADD STORAGE TO DIAGONAL AND RHS
      IF(ISS.NE.0) GO TO 201
      TLED=1.0D0/DELT
      DO 200 K=1,NLAY
      NNDLAY = NODLAY(K)
      NSTRT = NODLAY(K-1)+1
C
C6A------SEE IF THIS LAYER IS CONVERTIBLE OR NON-CONVERTIBLE.
      IF(LAYCON(K).EQ.3 .OR. LAYCON(K).EQ.2) GO TO 150
      IF(LAYCON(K).EQ.4.OR.LAYCON(K).EQ.5) GO TO 160
C6B------NON-CONVERTIBLE LAYER, SO USE PRIMARY STORAGE
      DO 140 N=NSTRT,NNDLAY
        IF(IBOUND(N).EQ.0) GO TO 140
        RHO=SC1(N)*TLED
        AMAT(IA(N))=AMAT(IA(N))-RHO
        RHS(N)=RHS(N)-RHO*HOLD(N)
  140 CONTINUE
      GO TO 200
C
C6C------A CONVERTIBLE LAYER, SO CHECK OLD AND NEW HEADS TO DETERMINE
C6C------WHEN TO USE PRIMARY AND SECONDARY STORAGE
  150 CONTINUE
      DO 180 N=NSTRT,NNDLAY
C
C--------IF THE CELL IS EXTERNAL THEN SKIP IT.
      IF(IBOUND(N).EQ.0) GO TO 180
      TP=TOP(N)
      RHO2=SC2(N)*TLED
      RHO1=SC1(N)*TLED
C
C6C1-----FIND STORAGE FACTOR AT START OF TIME STEP.
      SOLD=RHO2
      IF(HOLD(N).GT.TP) SOLD=RHO1
C
C6C2-----FIND STORAGE FACTOR AT END OF TIME STEP.
      HTMP=HNEW(N)
      SNEW=RHO2
      IF(HTMP.GT.TP) SNEW=RHO1
C
C6C3-----ADD STORAGE TERMS TO RHS AND TO DIAGONAL OF AMAT.
      AMAT(IA(N))=AMAT(IA(N))-SNEW
      RHS(N)=RHS(N) - SOLD*(HOLD(N)-TP) - SNEW*TP
C
  180 CONTINUE
      GO TO 200
C
C6D------A CONVERTIBLE LAYER WITH CONTINUOUS FUNCTION,
C6D------SO USE NEWTON EXPANSION OF STORGE TERM
  160 CONTINUE
      DO 190 N=NSTRT,NNDLAY
C
C--------IF THE CELL IS EXTERNAL THEN SKIP IT.
      IF(IBOUND(N).EQ.0) GO TO 190
C6D1-----COMPUTE PORE STORAGE TERM AS PER NEWTON RAPHSON
      SATO = So(N)
      BBOT=BOT(N)
      TTOP=TOP(N)
      TOTTHICK = TTOP - BBOT
      SATN = Sn(N)
      EPS = 1.0E-4
      HD=HNEW(N)+ EPS
      CALL SAT_THIK(N,HD,TOTTHICK,BBOT,THCK,K,TTOP)
      DS = (THCK - SATN)/EPS
      IF(DS.LT.1.0E-20) DS = 1.0E-20
      RHO2  = SC2(N) * DS * TOTTHICK * TLED
      AMAT(IA(N)) = AMAT(IA(N)) - RHO2
      RHS(N) = RHS(N) - RHO2*HNEW(N) + SC2(N)*TOTTHICK*TLED*(SATN-SATO)
C-----COMPUTE COMPRESSIBLE STORAGE TERM
c      RHO1 = SATN * SC1(N) * TLED
c      AMAT(IA(N)) = AMAT(IA(N)) - RHO1
c      RHS(N) = RHS(N) - RHO1 * HOLD(N)
C7C --------CHANGE SATN FOR ALTERNATE COMPRESSIBLE STORAGE FORMULATION
      IF(IALTSTO.EQ.1) THEN  ! COMPRESSIBLE STORAGE KICKS ON ONLY NEAR TOP OF CELL
        EPS = TOTTHICK * 0.001     !SMOOTH AT 0.1 PERCENT OF CELL THICKNESS
        EX = HNEW(N) - (TTOP - EPS)
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
C6D2----COMPUTE COMPRESSIBLE STORAGE TERM VIA NEWTON RAPHSON      
      RHO1 = SATN * SC1(N) * TLED
      FTERM = RHO1 * (HNEW(N) - HOLD(N))
      DFTERM = RHO1 + SC1(N) * TLED * (HNEW(N) - HOLD(N)) * DS
      AMAT(IA(N)) = AMAT(IA(N)) - DFTERM
      RHS(N) = RHS(N) - DFTERM * HNEW(N) + FTERM  
  190 CONTINUE
C
  200 CONTINUE
C
  201 CONTINUE
C--------------------------------------------------------------------------
C7-------FOR EACH LAYER DETERMINE IF CORRECTION TERMS ARE NEEDED FOR
C7-------FLOW DOWN INTO PARTIALLY SATURATED LAYERS.
C
C7A-----skip leakage correction if NOVFC flag is on
      IF(NOVFC.EQ.1) GO TO 301
      IF(.NOT. ALLOCATED(HWADIGW)) 
     1    ALLOCATE( HWADIGW(NODES))
      HWADIGW = 0.0
C7B-----FILL HWADI TERM FOR GW DOMAIN NODES 
      DO N=1,NODES
        X = HNEW(N) - TOP(N)
        CALL WADIFN(X,Y)
        HWADIGW(N) = Y + TOP(N)
      ENDDO
C7C-----loop over all layers for leakage correction
      DO 300 K=1,NLAY
C7A-----skip leakage correction for UPS formulation
      IF(laycon(k).eq.4.OR.LAYCON(K).EQ.5) GO TO 300
      NNDLAY = NODLAY(K)
      NSTRT = NODLAY(K-1)+1
C
C7D-----SEE IF THIS LAYER MAY NEED CORRECTION FOR LEAKAGE FROM ABOVE
      IF(K.EQ.1) GO TO 250
      IF(LAYCON(K).NE.3 .AND. LAYCON(K).NE.2) GO TO 250
C
C7D1----SEE IF HEAD IS BELOW TOP AND ADD CORRECTION TERMS TO RHS
      DO 220 N=NSTRT,NNDLAY
C
C7D2----IF CELL IS EXTERNAL (IBOUND<=0) THEN SKIP IT.
      IF(IBOUND(N).EQ.0) GO TO 220
      HTMP=HNEW(N)
C
C7D3----WITH HEAD BELOW TOP FIND PGF LOCATION AND CORRECT
C7D3----GO OVER CONNECTIONS OF NODE N FOR ALL CONNECTIONS ABOVE IT
      DO II = IA(N)+1,IA(N+1)-1
        JJ = JA(II)
        IF(JJ.GE.N) CYCLE !UPPER LAYER HAS LESSER NODE NUMBER
        IIS = JAS(II)
        IF(IVC(IIS).EQ.1)THEN !VERTICAL DIRECTION CONNECTION UP SO N IS DOWNSTREAM AND JJ IS UPSTREAM
        IF(IBOUND(JJ).EQ.0) CYCLE
C-----------PUT CORRECTION FOR SYMMETRIC MATRIX ON RHS OF N AND JJ
          RHS(N)=RHS(N) + AMAT(II)*(HWADIGW(N)-HTMP)
          RHS(JJ)=RHS(JJ) - AMAT(II)*(HWADIGW(N)-HTMP)
        ENDIF
      ENDDO
  220 CONTINUE
  250 CONTINUE
  300 CONTINUE
  301 CONTINUE
C
C8------RETURN
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE WADIFN(X,Y)
C     ******************************************************************
C     COMPUTE SMOOTH FUNCTION FOR WADI HEAD
C     ******************************************************************
C
C      SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,ONLY:WADIEPS
      DOUBLE PRECISION X,Y
C     ------------------------------------------------------------------
      WADIEPS = 1.0E-4
      IF(X.GT.WADIEPS)THEN
          Y = X
      ELSEIF(X.GT.-WADIEPS)THEN
          Y = 0.25*X**2/WADIEPS + 0.5*X + 0.25*WADIEPS
      ELSE
          Y=0.0
      ENDIF
C
C8------RETURN
      RETURN
      END
      SUBROUTINE DWADIFN(X,Y)
C     ******************************************************************
C     COMPUTE SMOOTH FUNCTION'S DERIVATIVE FOR WADI DH
C     ******************************************************************
C
C      SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,ONLY:WADIEPS
      DOUBLE PRECISION X,Y
C     ------------------------------------------------------------------
      WADIEPS = 1.0E-4
      IF(X.GT.WADIEPS)THEN
          Y = 1.0
      ELSEIF(X.GT.-WADIEPS)THEN
          Y = 0.5*X/WADIEPS + 0.5
      ELSE
          Y=0.0
      ENDIF
C
C8------RETURN
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE SGWF2FILLKRS
C     ******************************************************************
C     COMPUTE AKR AND AKRC OF FLOW TERM FOR POROUS MATRIX NODES
C     ******************************************************************
C
C      SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,     ONLY:NCOL,NROW,NLAY,IBOUND,HNEW,BUFF,BOT,TOP,ISYM,
     1                IOUT,NODLAY,AMAT,RHS,IA,JA,JAS,ICONCV,AKRC,AKR,Sn,
     2                NODES,NEQS,iunsat
      USE GWFBCFMODULE,ONLY:LAYCON
      USE SMSMODULE, ONLY: DKDH,NONMETH,DKDHC
C
      DOUBLE PRECISION HD,BBOT,TTOP,THCK,ZERO,CONSTERM,FLOWTERM,
     *  TERM,TOTTHICK,EKR,SW
C     ------------------------------------------------------------------
      ZERO=0.
C-----------------------------------------------------------------------------
C2------GET AKR FROM SATURATION FOR RICHARDS EQUATION FOR ALL NODES
      DO 100 K=1,NLAY
        IF(LAYCON(K).EQ.4.OR.LAYCON(K).EQ.5) THEN
          NNDLAY = NODLAY(K)
          NSTRT = NODLAY(K-1)+1
C
          DO 150 N=NSTRT,NNDLAY
          IF(IBOUND(N).NE.0) THEN
C1A--------------CALCULATE AKR FOR EACH NODE.
            SW=Sn(N)
            HD = HNEW(N)
            CALL KR_CAL(N,SW,EKR,K,HD)
            AKR(N) = EKR
          ENDIF
  150   CONTINUE
        ENDIF
  100 CONTINUE
C-----------------------------------------------------------------------------
C2------FILL AKRC WITH UPSTREAM AKR OF THE CONNECTION FOR ALL CONNECTIONS
C-----------------------------------------------------------------------------
      DO 200 K=1,NLAY
        NNDLAY = NODLAY(K)
        NSTRT = NODLAY(K-1)+1
        DO N=NSTRT,NNDLAY
          IF(IBOUND(N).NE.0)THEN
C2A---------GO OVER CONNECTIONS OF NODE N AND FILL FOR UPPER SYMMETRIC PART
            DO II = IA(N)+1,IA(N+1)-1
              JJ = JA(II)
              IF(JJ.GE.N.AND.IBOUND(JJ).NE.0.AND.JJ.LE.NODES)THEN
                IIS = JAS(II)
C2B---------FIND UPSTREAM NODE AND HIGHER BOT NODE
                IUPS = N
                IF(HNEW(JJ).GT.HNEW(N)) IUPS = JJ
                IHBOT = N
                IF(BOT(JJ).GT.BOT(N)) IHBOT = JJ
C2C---------FILL AKRC FOR CONNECTION BETWEEN N AND JJ
                INDK = 0
                IF(IUPS.EQ.IHBOT) INDK = 1
                IF(ABS(BOT(JJ)-BOT(N)).LT.0.01) INDK = 1
                IF(INDK.EQ.1)THEN
                  AKRC(IIS) = AKR(IUPS)
                ELSEIF(IBOUND(IUPS).NE.0) THEN
                  BBOT=BOT(IHBOT)
                  TTOP=TOP(IHBOT)
                  TOTTHICK = TTOP - BBOT
                  HD=HNEW(IUPS)
                  CALL SAT_THIK(N,HD,TOTTHICK,BBOT,SW,K,TTOP)
                  CALL KR_CAL(N,SW,EKR,K,HD)
                  AKRC(IIS) = EKR
                ENDIF
              ENDIF
            ENDDO
          ENDIF
        ENDDO
  200 CONTINUE
C4------RETURN
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE SGWF2BCFU1WDS(K,KITER,KSTP,KPER)
C     ******************************************************************
C     FOR LAYER K, ADJUST WET/DRY CELLS AND COMPUTE CONDUCTANCE  FROM
C     SATURATED THICKNESS AND HYDRAULIC CONDUCTIVITY FOR LAYCON OF 1 OR 3
C     ******************************************************************
C
C      SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,     ONLY:NCOL,NROW,NLAY,IBOUND,HNEW,BUFF,BOT,TOP,ISYM,
     1                IOUT,NODLAY,AMAT,IA,JA,JAS,PGF,IUNSTR,IVC
      USE GWFBCFMODULE,ONLY:IWDFLG,WETFCT,IHDWET,IWETIT,LAYCON,HDRY,
     1                      HK,WETDRY,LAYAVG,laywet
C
      DOUBLE PRECISION HD,BBOT,TTOP,THCK,HTMP,ANUM
      CHARACTER*3 ACNVRT
      DIMENSION ICNVRT(5),JCNVRT(5),ACNVRT(5)
C     ------------------------------------------------------------------
C
C1------LOOP THROUGH EACH CELL IN LAYER AND CALCULATE TRANSMISSIVITY AT
C1------EACH ACTIVE CELL.
      ZERO=0.
      NCNVRT=0
      IHDCNV=0
      ITFLG=1
      IF(laywet(k).NE.0) ITFLG=MOD(KITER,IWETIT)
C
      NNDLAY = NODLAY(K)
      NSTRT = NODLAY(K-1)+1
      DO 200 N=NSTRT,NNDLAY
C
C2------IF CELL IS ACTIVE, THEN SKIP TO CODE THAT CALCULATES SATURATED
C2------THICKNESS.
      IF(IBOUND(N).NE.0) GO TO 20
C
C3------DETERMINE IF THE CELL CAN CONVERT BETWEEN CONFINED AND
C3------UNCONFINED.  IF NOT, SKIP TO CODE THAT SETS TRANSMISSIVITY TO 0.
      IF(ITFLG.NE.0) GO TO 6
      IF(WETDRY(N).EQ.ZERO)GO TO 6
      WD=WETDRY(N)
      IF(WD.LT.ZERO) WD=-WD
      TURNON=BOT(N)+WD
C
C3A-----CHECK HEAD IN CELLS BELOW TO SEE IF WETTING THRESHOLD HAS BEEN
C3A-----REACHED.
      DO II = IA(N)+1,IA(N+1)-1
        JJ = JA(II)
        IIS = JAS(II)
        IF(IVC(IIS).EQ.1.AND.JJ.GT.N)THEN !CELL IS BELOW DIRECTION
          HTMP=HNEW(JJ)
          IF(IBOUND(JJ).GT.0.AND.HTMP.GE.TURNON)GO TO 9
        ENDIF
        IF(IVC(IIS).EQ.1.AND.JJ.LT.N)GO TO 16 !CELL IS ABOVE DIRECTION
C
C3B-----CHECK HEAD IN ADJACENT HORIZONTAL CELLS TO SEE IF WETTING
C3B-----THRESHOLD HAS BEEN REACHED.
        IF(WETDRY(N).LT.ZERO) GO TO 16
        HTMP=HNEW(JJ)
        IF(IBOUND(JJ).GT.0.AND.IBOUND(JJ).NE.30000.AND.
     1                            HTMP.GE.TURNON)GO TO 9
   16   CONTINUE
      ENDDO
C
C3C-----CELL IS DRY AND STAYS DRY.  SET TRANSMISSIVITY TO 0, SET
C3C-----SATURATED THICKNESS (BUFF) TO 0, AND SKIP TO THE NEXT CELL.
    6 CONTINUE
      IF(LAYAVG(K).EQ.3.OR.LAYAVG(K).EQ.4) BUFF(N)=ZERO
      GO TO 200
C
C4------CELL BECOMES WET.  SET INITIAL HEAD AND VERTICAL CONDUCTANCE.
    9 IF(IHDWET.NE.0) HNEW(N)=BOT(N)+WETFCT*WD
      IF(IHDWET.EQ.0) HNEW(N)=BOT(N)+WETFCT*(HTMP-BOT(N))
      DO II = IA(N)+1,IA(N+1)-1
        JJ = JA(II)
        IF(JJ.LE.N) CYCLE
        IIS = JAS(II)
        IF(IVC(IIS).EQ.1)THEN
          IF(IBOUND(JJ).NE.0) AMAT(II)= PGF(IIS)
          IS = ISYM(II)
          IF(IBOUND(JJ).NE.0) AMAT(IS)= PGF(IIS)

        ENDIF
      ENDDO
   14 IBOUND(N)=30000
C---------------------------------------------------------------------------
C
C4A-----PRINT MESSAGE SAYING CELL HAS BEEN CONVERTED TO WET.
      NCNVRT=NCNVRT+1
      IF(IUNSTR.EQ.0)THEN !GET LAYER, ROW AND COLUMN FOR STRUCTURED GRID
        KLAYER = (N) / (NCOL*NROW) + 1
        IJ = N - (KLAYER-1)*NCOL*NROW
        IROW = (IJ-1)/NCOL + 1
        JCOLMN = IJ - (IROW-1)*NCOL
        ICNVRT(NCNVRT)=IROW
        JCNVRT(NCNVRT)=JCOLMN
      ELSE
        ICNVRT(NCNVRT)=N
        JCNVRT(NCNVRT)=0
      ENDIF
        ACNVRT(NCNVRT)='WET'
      IF(NCNVRT.LT.5) GO TO 20
      IF(IUNSTR.EQ.0)THEN
         IF(IHDCNV.EQ.0) WRITE(IOUT,17) KITER,K,KSTP,KPER
   17    FORMAT(1X,/1X,'CELL CONVERSIONS FOR ITER.=',I3,'  LAYER=',
     1    I3,'  STEP=',I7,'  PERIOD=',I7,'   (ROW,COL)')
         IHDCNV=1
         WRITE(IOUT,18) (ACNVRT(L),ICNVRT(L),JCNVRT(L),L=1,NCNVRT)
   18    FORMAT(1X,3X,5(A,'(',I7,',',I7,')   '))
         NCNVRT=0
       ELSE
         IF(IHDCNV.EQ.0) WRITE(IOUT,37) KITER,K,KSTP,KPER
   37    FORMAT(1X,/1X,'CELL CONVERSIONS FOR ITER.=',I3,'  LAYER=',
     1    I7,'  STEP=',I7,'  PERIOD=',I7,'   (NODE)')
         IHDCNV=1
         WRITE(IOUT,38) (ACNVRT(L),ICNVRT(L),L=1,NCNVRT)
   38    FORMAT(1X,3X,5(A,'(',I7,')   '))
         NCNVRT=0
       ENDIF
C----------------------------------------------------------------------------
C
C5------CALCULATE SATURATED THICKNESS.
   20 HD=HNEW(N)
      BBOT=BOT(N)
      IF(LAYCON(K).EQ.1) GO TO 50
      TTOP=TOP(N)
      IF(BBOT.GT.TTOP) THEN
         IF(IUNSTR.EQ.0)THEN
         WRITE(IOUT,35) K,I,J
   35    FORMAT(1X,'Negative cell thickness at (Layer,row,col)',
     1   I4,',',I4,',',I4)
         STOP
         ELSE
         WRITE(IOUT,36) N
36       FORMAT(1X,'Negative cell thickness at node',i6)
         ENDIF
      END IF
      IF(HD.GT.TTOP) HD=TTOP
   50 THCK=HD-BBOT
C
C6------CHECK TO SEE IF SATURATED THICKNESS IS GREATER THAN ZERO.
      IF(THCK.LE.ZERO) GO TO 100
C
C6A-----IF SATURATED THICKNESS>0 THEN CALCULATE TRANSMISSIVITY
C6A-----IN BUFF AND MOVE TO NEXT NODE.
         BUFF(N)=THCK * HK(N)
      GO TO 200
C--------------------------------------------------------------------------
C
C6B-----WHEN SATURATED THICKNESS < 0, PRINT A MESSAGE AND SET
C6B-----ROW IN AMAT, IBOUND, AND VERTICAL CONDUCTANCE =0
  100 NCNVRT=NCNVRT+1
       IF(IUNSTR.EQ.0)THEN ! GET LAYER, ROW AND COLUMN FOR STRUCTURED GRID
        KLAYER = (N) / (NCOL*NROW) + 1
        IJ = N - (KLAYER-1)*NCOL*NROW
        IROW = (IJ-1)/NCOL + 1
        JCOLMN = IJ - (IROW-1)*NCOL
        ICNVRT(NCNVRT)=IROW
        JCNVRT(NCNVRT)=JCOLMN
      ELSE
        ICNVRT(NCNVRT)=N
        JCNVRT(NCNVRT)=0
      ENDIF
      ACNVRT(NCNVRT)='DRY'
      IF(NCNVRT.LT.5) GO TO 150
      IF(IUNSTR.EQ.0)THEN
         IF(IHDCNV.EQ.0) WRITE(IOUT,17) KITER,K,KSTP,KPER
         IHDCNV=1
         WRITE(IOUT,18) (ACNVRT(L),ICNVRT(L),JCNVRT(L),L=1,NCNVRT)
         NCNVRT=0
       ELSE
         IF(IHDCNV.EQ.0) WRITE(IOUT,37) KITER,K,KSTP,KPER
         IHDCNV=1
         WRITE(IOUT,38) (ACNVRT(L),ICNVRT(L),L=1,NCNVRT)
         NCNVRT=0
       ENDIF
C-------------------------------------------------------------------------
C
  150 HNEW(N)=HDRY
C
      DO II = IA(N)+1,IA(N+1)-1
        AMAT(II) = ZERO
      ENDDO
C
      IF(IBOUND(N).GE.0) GO TO 160
         WRITE(IOUT,151)
  151    FORMAT(1X,/1X,'CONSTANT-HEAD CELL WENT DRY',
     1          ' -- SIMULATION ABORTED')
      IF(IUNSTR.EQ.0)THEN !GET LAYER, ROW AND COLUMN FOR STRUCTURED GRID
        K = N / (NCOL*NROW) + 1
        IJ = N - (K-1)*NCOL*NROW
        I = (IJ-1)/NCOL + 1
        J = IJ - (I-1)*NCOL
        WRITE(IOUT,152) K,I,J,KITER,KSTP,KPER
  152   FORMAT(1X,'LAYER=',I7,'   ROW=',I7,'   COLUMN=',I7,
     1    '   ITERATION=',I8,'   TIME STEP=',I8,'   STRESS PERIOD=',I8)
      ELSE
        WRITE(IOUT,1152) N,KITER,KSTP,KPER
 1152   FORMAT(1X,'NODE =',I9,
     1    '   ITERATION=',I8,'   TIME STEP=',I8,'   STRESS PERIOD=',I8)
      ENDIF
         WRITE(IOUT,*) BBOT,HD
         STOP
  160 IBOUND(N)=0
  200 CONTINUE
C
C7------PRINT ANY REMAINING CELL CONVERSIONS NOT YET PRINTED
      IF(NCNVRT.EQ.0) GO TO 203
      IF(IUNSTR.EQ.0)THEN
         IF(IHDCNV.EQ.0) WRITE(IOUT,17) KITER,K,KSTP,KPER
         IHDCNV=1
         WRITE(IOUT,18) (ACNVRT(L),ICNVRT(L),JCNVRT(L),L=1,NCNVRT)
         NCNVRT=0
       ELSE
         IF(IHDCNV.EQ.0) WRITE(IOUT,37) KITER,K,KSTP,KPER
         IHDCNV=1
         WRITE(IOUT,38) (ACNVRT(L),ICNVRT(L),L=1,NCNVRT)
         NCNVRT=0
       ENDIF
C
C8------CHANGE IBOUND VALUE FOR CELLS THAT CONVERTED TO WET THIS
C8------ITERATION FROM 30000 to 1.
  203 IF(laywet(k).EQ.0) GO TO 210
      DO 205 N=NSTRT,NNDLAY
      IF(IBOUND(N).EQ.30000) IBOUND(N)=1
  205 CONTINUE
  210 CONTINUE
C9------return
      RETURN
      END
C------------------------------------------------------------------------
      SUBROUTINE SGWF2BCFU1H(K,KITER)
C     ******************************************************************
C     FOR LAYER K, FILL MATRIX AMAT WITH FLOW TERMS FOR LAYCON OF 1 OR 3
C     ******************************************************************
C
C      SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,     ONLY:NCOL,NROW,NLAY,IBOUND,HNEW,BUFF,BOT,TOP,ISYM,
     1                    IOUT,NODLAY,AMAT,IA,CL1,CL2,JA,JAS,PGF,IUNSTR,
     1                    IVC,NODES
      USE GWFBCFMODULE,ONLY:IWDFLG,WETFCT,IHDWET,IWETIT,LAYCON,HDRY,
     1                      HK,WETDRY,LAYAVG,laywet
C
      DOUBLE PRECISION HD,BBOT,TTOP,THCK,HTMP,ANUM,THIK1,THIK2
      CHARACTER*3 ACNVRT
      DIMENSION ICNVRT(5),JCNVRT(5),ACNVRT(5)
C     ------------------------------------------------------------------
C
C1------Set Flags
      ZERO=0.
      NCNVRT=0
      IHDCNV=0
      ITFLG=1
      IF(laywet(k).NE.0) ITFLG=MOD(KITER,IWETIT)
C2------Get start and end node number for layer
      NNDLAY = NODLAY(K)
      NSTRT = NODLAY(K-1)+1
C
C3------COMPUTE HORIZONTAL BRANCH CONDUCTANCES FROM TRANSMISSIVITY.
C3------AND FILL VERTICAL CONDUCTANCE INTO AMAT FROM PGF
C
C3A--------loop over all nodes of the layer
      DO N=NSTRT,NNDLAY
        IF(IBOUND(N).NE.0)THEN
C3B---------loop over all connections of node
          DO II = IA(N)+1,IA(N+1)-1
            JJ = JA(II)
C3C-------------skip computations for backward connection or inactive node
            IF(JJ.LE.N.OR.IBOUND(JJ).EQ.0.OR.JJ.GT.NODES) GO TO 301
            IIS = JAS(II)
            ISLOC = ISYM(II)
C3D-------------compute conductance term for horizontal connecction
            IF(IVC(IIS).EQ.0.OR.IVC(IIS).EQ.2)THEN !HORIZONTAL DIRECTION
              EL1 = CL1(IIS)
              EL2 = CL2(IIS)
              THIK1 = 1.0
              THIK2 = 1.0
              CALL CBCK12(II,N,JJ,K,EL1,EL2,ANUM,BUFF,
     *              THIK1,THIK2)
C3D1---------------fill term in AMAT for forward and backward connection
              AMAT(II) = PGF(IIS)*ANUM
              AMAT(ISLOC) = PGF(IIS)*ANUM
C3F---------------copy vertical conductance from PGF to AMAT
            ELSE
              AMAT(II) = PGF(IIS)
              AMAT(ISLOC) = PGF(IIS)
            ENDIF
  301       CONTINUE
          ENDDO
C
        ENDIF
      ENDDO
C
C4------RETURN.
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE SGWF2BCFU1H4(K)
C     ******************************************************************
C     COMPUTE CONDUCTANCE FOR ONE LAYER FROM PGF AND SATURATED THICKNESS
C     FOR LAYCON OF 4 - ALSO FILL NEWTON TERM IF NR IS ACTIVE WITH LAYCON=4
C     ******************************************************************
C
C      SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,    ONLY:NCOL,NROW,NLAY,NODES,IBOUND,HNEW,BUFF,BOT,TOP,
     1 AKRC,ISYM,IOUT,NODLAY,AMAT,RHS,IA,JA,JAS,PGF,IVC,ICONCV,Sn,iunsat
      USE GWFBCFMODULE,ONLY:IWDFLG,WETFCT,IHDWET,IWETIT,LAYCON,HDRY,
     1                      WETDRY,LAYAVG
c      USE SMSMODULE, ONLY: DKDH,NONMETH
C
      DOUBLE PRECISION HD,BBOT,TTOP,THCK,TOTTHICK,ZERO
C     ------------------------------------------------------------------
C1------Initialize Flags
      ZERO=0.
      if(LAYCON(K).EQ.5) iconcv = 0
C
C----------------------------------------------------------------------------
C2----COMPUTE HORIZONTAL BRANCH CONDUCTANCES FROM UPSTREAM SATURATED
C2----THICKNESS AND PGF
      NNDLAY = NODLAY(K)
      NSTRT = NODLAY(K-1)+1
C2A--------loop over all nodes of the layer
      DO N=NSTRT,NNDLAY
        IF(IBOUND(N).NE.0)THEN
C2B---------loop over all connections of node
          DO II = IA(N)+1,IA(N+1)-1
            JJ = JA(II)
C2C-------------skip computations for backward connection or inactive node
            IF(JJ.GT.N.AND.IBOUND(JJ).NE.0.AND.JJ.LE.NODES)THEN
              IIS = JAS(II)
C2D-------------if vertical direction has constant CV then fill directly from PGF
              IF(IVC(IIS).EQ.1.AND.ICONCV.EQ.1)THEN !VERTICAL DIRECTION HAS CONSTANT CV
                AMAT(II) = PGF(IIS)
                AMAT(ISYM(II)) = PGF(IIS)
C2E-------------fill AMAT with non-constant term times PGF
              ELSE !HORIZONTAL DIRECTION OR NON-CONSTANT CV
                AMAT(II) = PGF(IIS)*AKRC(IIS)
                AMAT(ISYM(II)) = PGF(IIS)*AKRC(IIS)
              ENDIF
            ENDIF
          ENDDO
        ENDIF
      ENDDO
C
C3-----RETURN.
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE SGWF2BCFU1WDS4(K)
C     ******************************************************************
C     COMPUTE SATURATED THICKNESS FOR LAYCON OF 4 FOR EACH LAYER AND STORE IN Sn.
C     ******************************************************************
C
C      SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,     ONLY:NCOL,NROW,NLAY,IBOUND,HNEW,BUFF,BOT,TOP,ISYM,
     1               IOUT,NODLAY,AMAT,RHS,IA,JA,PGF,ICONCV,Sn,iunsat
      USE GWFBCFMODULE,ONLY:IWDFLG,WETFCT,IHDWET,IWETIT,LAYCON,HDRY,
     1                      WETDRY,LAYAVG
c      USE SMSMODULE, ONLY: DKDH,NONMETH
C
      DOUBLE PRECISION HD,BBOT,TTOP,THCK,TOTTHICK,ZERO
C     ------------------------------------------------------------------
      ZERO=0.
      if(LAYCON(K).EQ.5) iconcv = 0
C
C1------LOOP THROUGH EACH CELL IN LAYER
C
      NNDLAY = NODLAY(K)
      NSTRT = NODLAY(K-1)+1
      DO 200 N=NSTRT,NNDLAY
        IF(IBOUND(N).NE.0) THEN
C
C2--------CALCULATE SATURATED THICKNESS.
          HD=HNEW(N)
          BBOT=BOT(N)
          TTOP=TOP(N)
          TOTTHICK = TTOP - BBOT
          CALL SAT_THIK(N,HD,TOTTHICK,BBOT,THCK,K,TTOP)
C
C3--------STORE IN Sn AND MOVE TO NEXT NODE.
          Sn(N)=THCK
        ENDIF
  200 CONTINUE
C
C4------RETURN.
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE KR_CAL(N,SW,EKR,K,HD)
C     ******************************************************************
C     COMPUTE KR FOR RICHARDS EQUATION FORMULATION
C     ******************************************************************
C
C      SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,     ONLY:ISSFLG,iunsat,TOP,BOT,NODES
      USE GWFBCFMODULE,ONLY:alpha,beta,sr,brook,LAYCON,
     *  retcrvs,nutabrows,nuzones,iuzontab,ITABRICH,INTRICH  
      DOUBLE PRECISION EKR,SW,seff,one, tabS,pc, topS, botS,cellS, 
     *  BBOT,TTOP,HD,bott,mul,div
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
C      
C      IF(LAYCON(K).EQ.4)THEN
      IF(METHOD.EQ.3) THEN  
        EKR = SW
cc        IF(EKR.LT.1.0E-7) EKR = 1.0E-7
C      ELSEIF(LAYCON(K).EQ.5)THEN
      ELSEIF(METHOD.EQ.4) THEN   
C---------Unsaturated zone functions
        IF(SW.LT.SR(N)) THEN
          EKR = 0.0
        ELSE  
          one = 1.0
          SEFF = (SW - Sr(N))/(one - Sr(N))
cc ......use van genuchten kr when brookn is less than zero         
          IF(BROOK(N). LT. 0.0) THEN  
           IF(Seff. GT. one) Seff = one   
cc ........van genuchten is below
            EKR = (Seff**0.5)*(one-
     1        (one-Seff**(one/(one-one/beta(n))))**(one-one/beta(n)))**2
          ELSE  
Cc ........brooks corey is below               
            EKR = SEFF ** BROOK(N)
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
          izon = iuzontab(n)
          iv = 3  !variable is Kr
          EKR = tabS(pc,retcrvs,izon,iv,nutabrows,nuzones,ttop)
        endif
      ELSEIF(METHOD.EQ.6) THEN
C---------FOR USE WITH TABULAR INPUT OF INTERPOLATED TABLE
        iv = 3  !variable is Kr
        izon = iuzontab(n)
        BBOT = BOT(N)
        bott = bbot   ! bott changes to water table if below water table in tabS        
        TTOP = TOP(N)
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
C---------------------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE SAT_THIK(N,HD,TOTTHICK,BBOT,THCK,K,TTOP)
C     ******************************************************************
C     COMPUTE Kr FOR LAYCON OF 4 AND 5
C     ******************************************************************
C
C      SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,     ONLY:ISSFLG,NODES,iunsat,HNEW,INDDF 
      USE DDFMODULE, ONLY: ZETASWI,ISHARP
cs      USE DDFMODULE, ONLY: ZETASWI,ISHARP
      USE GWFBCFMODULE,ONLY:alpha,beta,brook,sr,BP,IBPN,LAYCON,IDRY,
     *  retcrvs,nutabrows,nuzones,iuzontab,ITABRICH,INTRICH 
      DOUBLE PRECISION THCK,HD,BBOT,TTOP,X,S,V,COF1,COF2,FACTOR1,FACTOR2
     *  ,eps,acof,y,TOTTHICK,seff,gamma, tabS,pc, topS, botS,cellS,bott,
     *  mul,div,THICKKP,ZETA  
C     ------------------------------------------------------------------
C
C ----DEPENDING ON THE METHOD, FIND THE SATURATED THICKNESS
C
C--------------------------------------------------------
      ZERO=0.
      METHOD = 3  !  METHOD 3 IS THE UNCONFINED SATURATED THICKNESS CURVE OF NWT
C1-------USE METHODS 4, 5 AND 6 FOR RICHARDS EQUATION      
      IF(N.LE.NODES)THEN
        IF(LAYCON(K).EQ.5) THEN ! METHOD 4 USES VAN GENUCHTEN AND BROOKS-COREY FUNCTIONS
          METHOD = 4  
          IF(ITABRICH.NE.0)THEN !METHOD 5 USES TABULAR INPUT
            METHOD = 5
            IF(INTRICH.NE.0)THEN  !METHOD 6 USES INTEGRATED FUNCTIONS (DISCARD THIS)
              METHOD = 6
            ENDIF
          ENDIF
        ENDIF  
      ENDIF  
C-------------------------------------------------------------------------      
C2--------ADJUST DIMENSIONS FOR FRESHWATER AND SALTWATER EQUATIONS OF SWI      
      IF(INDDF. NE. 0) THEN 
cs      IF(INDDF. NE. 0) THEN 
C2A-----RECOMPUTE ZETA          
        IF(ISHARP.NE.0) CALL ZETANODE (N, ZETA,HD)  
cs        IF(ISHARP.NE.0) CALL ZETANODE (N, ZETA,HD)  
C2B--------SELECT EQUATION        
        IEQSWI = 1 ! RIGHT NOW HARDWIRE FOR ONLY FRESHWATER ****************
        IF(ISHARP.NE.0 .AND. IEQSWI .EQ. 1) THEN 
          THICKKP = TOTTHICK
cs        IEQSWI = 1 ! RIGHT NOW HARDWIRE FOR ONLY FRESHWATER ****************
cs        IF(ISHARP.NE.0 .AND. IEQSWI .EQ. 1) THEN 
cs          THICKKP = TOTTHICK
C2C--------- EQUATION IS FRESHWATER EQUATION
          IF(ZETA. GT. BBOT) THEN 
           BBOT = ZETA    
           TOTTHICK = TTOP - BBOT
           IF(TOTTHICK. LE. 0.0) THEN
              THCK = 0.0
              TOTTHICK = THICKKP
              RETURN
            ENDIF 
          ENDIF  
        ELSEIF(ISHARP. NE.0. AND. IEQSWI. EQ. 2) THEN  
cs          IF(ZETA. GT. BBOT) THEN 
cs           BBOT = ZETA    
cs           TOTTHICK = TTOP - BBOT
cs           IF(TOTTHICK. LE. 0.0) THEN
cs              THCK = 0.0
cs              TOTTHICK = THICKKP
cs              RETURN
cs            ENDIF 
cs          ENDIF  
cs        ELSEIF(ISHARP. NE.0. AND. IEQSWI. EQ. 2) THEN  
C2D--------- EQUATION IS SALTWATER EQUATION
          IF(HNEW(N). LT. TTOP) THEN         ! THIS SHOULD BE HNEW OF EQUATION 1 FRESHWATER *********
            TTOP = HNEW(N)                   ! THIS SHOULD BE HNEW OF EQUATION 1 FRESHWATER *********
            TOTTHICK = TTOP - BBOT
            HD = ZETA
            IF(TOTTHICK. LE. 0.0) THEN
              THCK = 0.0
              TOTTHICK = THICKKP
              RETURN
            ENDIF 
          ENDIF 
        ENDIF  
      ENDIF    
cs          IF(HNEW(N). LT. TTOP) THEN         ! THIS SHOULD BE HNEW OF EQUATION 1 FRESHWATER *********
cs            TTOP = HNEW(N)                   ! THIS SHOULD BE HNEW OF EQUATION 1 FRESHWATER *********
cs            TOTTHICK = TTOP - BBOT
cs            HD = ZETA
cs            IF(TOTTHICK. LE. 0.0) THEN
cs              THCK = 0.0
cs              TOTTHICK = THICKKP
cs              RETURN
cs            ENDIF 
cs          ENDIF 
cs        ENDIF  
cs      ENDIF    
C      if(IUNSAT.EQ.1.AND.N.LE.NODES) method = 4
C------------------------------------------------------------------------------      
      IF(METHOD.EQ.1)THEN
C3-------STRAIGHT LINE, NO SMOOTHING
        TTOP = BBOT + TOTTHICK
        IF(HD.GT.TTOP) HD=TTOP
        THCK = (HD - BBOT) / TOTTHICK
        IF(THCK.LT.ZERO) THCK=0.0
      ELSEIF(METHOD.EQ.2)THEN
C4-------STRAIGHT LINE WITH CUBIC SMOOTHING
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
C5-------STRAIGHT LINE WITH PARABOLIC SMOOTHING
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
C6---------vanG FUNCTION WITH MODIFICATIONS FOR BUBBLEPT AND FULLYDRY
        bpn = 0.0 
        if(ibpn.eq.1)then
          bpn = bp(N)
        endif
        TTOP = BBOT + TOTTHICK
        pc = 0.5*(ttop+bbot) - (hd-bpn)
        if(pc.le.0)then
          thck = 1.0
        else
          gamma = 1.-1./beta(n)
          Seff = (1. + (alpha(n)*pc)**beta(n))**gamma
          Seff = 1.0 / Seff
          if(idry.eq.0) then
            thck = seff * (1-sr(n)) + sr(n)
          else
            thck = seff
          endif  
        endif
C------------------------------------
      ELSEIF(METHOD.EQ.5) THEN
C7-------FOR TABULAR INPUT FIND SATURATION FROM TABLE          
        TTOP = BBOT + TOTTHICK
        pc = 0.5*(ttop+bbot) - hd
        if(pc.le.0)then
          thck = 1.0  
        else  
          izon = iuzontab(n)  
          iv = 2  !variable is saturation
          thck = tabS(pc,retcrvs,izon,iv,nutabrows,nuzones,ttop)
        endif    
      ELSEIF(METHOD.EQ.6) THEN
C8---------FOR USE WITH TABULAR INPUT OF INTERPOLATED TABLE  
        iv = 2  !variable is saturation
        izon = iuzontab(n) 
        TTOP = BBOT + TOTTHICK
        bott = bbot   ! bott changes to water table if below water table in tabS
        topS = tabS((ttop-hd),retcrvs,izon,iv,nutabrows,nuzones,ttop)
        botS = tabS((bbot-hd),retcrvs,izon,iv,nutabrows,nuzones,bott)
        if(hd.lt. bbot) then 
c9----------water level is below bottom so integrate curve from top to bottom
          if((ttop-bott).gt.1.0e-5)then
            cellS=-(topS - botS)/(ttop - bott)
            if(cellS.lt.1.0e-20) cellS = 1.0e-20
          else ! past bottom of table so use last two entries to compute value
           mul = retcrvs(iv,nutabrows,izon)-retcrvs(iv,nutabrows-1,izon)
           div = retcrvs(1,nutabrows,izon)-retcrvs(1,nutabrows-1,izon)
           cellS = mul/div   
          endif    
        elseif(hd.gt. ttop) then
c10----------water level is above top so integrated saturation of cell is 1
          cellS = 1.0
        else
c11----------water level is within cell so integrate 
          cellS = (botS - topS + 1.0 * (hd -bbot))/(ttop - bbot)   
        endif
        thck = cellS 
C-------------------------------------------------------------------
      ENDIF
C-------------------------------------------------------------------------      
C2--------GET BACK TOTTHICK AND MAKE THCK FRACTION OF THAT      
      IF(INDDF. NE. 0) THEN
        IF(ISHARP.NE.0 ) THEN 
cs      IF(INDDF. NE. 0) THEN
cs        IF(ISHARP.NE.0 ) THEN 
C2A--------- EQUATION IS FRESHWATER OR SALTWATER
          THCK = THCK * TOTTHICK/THICKKP
          TOTTHICK = THICKKP
        ENDIF  
      ENDIF    
cs          THCK = THCK * TOTTHICK/THICKKP
cs          TOTTHICK = THICKKP
cs        ENDIF  
cs      ENDIF    
C      if(IUNSAT.EQ.1.AND.N.LE.NODES) method = 4
C------------------------------------------------------------------------------    
      RETURN
      END
C-----------------------------------------------------------------------
      double precision function tabS 
     *       (pc,retcrvs,izon,iv,nutabrows,nuzones,celel)
C     ******************************************************************
C     compute value of S (iv=2) or Kr (iv=3) from table for given pc
C     ******************************************************************      
C           
C       specifications
C     ------------------------------------------------------------------
      real*8 pc, hd, a, b, hir, fn ,Fir, celel
      dimension retcrvs(3,nutabrows,nuzones)
C     ------------------------------------------------------------------
c       ! set at first location if head is less than at first location
      hd = retcrvs(1,1,izon)
      if(pc. le. hd) then 
        tabS =  retcrvs(iv,1,izon) 
        celel = pc + hd
        return
      endif    
c      ! set at last location if head is greater than at last location 
      hd = retcrvs(1,nutabrows,izon)
      if(pc. ge. hd) then 
        tabS =  retcrvs(iv,nutabrows,izon) 
        celel = pc + hd
        return
      endif          
c      ! find location in table and interpolate for anything in between 
cc      ifound = 0
      do irow = 1,nutabrows
         hir = retcrvs(1,irow,izon)
         if (pc. gt. hir) go to 10 
cc         ifound = 1
         ir = irow
         irm1 = irow - 1
         a = retcrvs(iv,irow,izon) - retcrvs(iv,irm1,izon)
         b = hir - retcrvs(1,irm1,izon)
         fn = a / b * (pc - hir)
         fir = retcrvs(iv,irow,izon)
         tabS = Fir + fn
         go to 20
10       continue         
      enddo    
20    continue       
cc      ! pick end value when overshooting
cc      if(ifound.eq.0)then
cc        tabS = retcrvs(iv,nutabrows,izon) 
cc      endif
c      
      return
      end function tabS
C---------------------------------------------------------------------------------
      SUBROUTINE GWF2BCFU1AD(KPER)
C     ******************************************************************
C     SET HOLD TO BOT WHENEVER A WETTABLE CELL IS DRY
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,      ONLY:NCOL,NROW,NLAY,ISSFLG,IBOUND,HOLD,BOT,NODLAY
      USE GWFBCFMODULE,ONLY:IWDFLG,WETDRY,LAYCON,CV
C     ------------------------------------------------------------------
C
      ISS=ISSFLG(KPER)
C
C1------RETURN IF STEADY STATE OR IF NOT USING WETTING CAPABILITY
      IF(IWDFLG.EQ.0 .OR. ISS.NE.0) RETURN
C
C2------LOOP THROUGH ALL LAYERS TO SET HOLD=BOT IF A WETTABLE CELL IS DRY
      ZERO=0.
      KB=0
      DO 100 K=1,NLAY
C
C2A-----SKIP LAYERS THAT CANNOT CONVERT BETWEEN WET AND DRY
      IF(LAYCON(K).NE.3 .AND. LAYCON(K).NE.1) GO TO 100
      KB=KB+1
      NNDLAY = NODLAY(K)
      NSTRT = NODLAY(K-1)+1
      DO 90 N=NSTRT,NNDLAY
C
C2B-----SKIP CELLS THAT ARE CURRENTLY WET OR ARE NOT WETTABLE
      IF(IBOUND(N).NE.0) GO TO 90
      IF(WETDRY(N).EQ.ZERO) GO TO 90
C
C2C-----SET HOLD=BOT
      HOLD(N)=BOT(N)
   90 CONTINUE
  100 CONTINUE
C
C3-----RETURN
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE GWF2BCFU1BDS(KSTP,KPER,ITVMUNIT,iuddf)
C     ******************************************************************
C     COMPUTE STORAGE BUDGET FLOW TERM FOR BCF.
C     ******************************************************************
C
C     SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL, ONLY:NCOL,NROW,NLAY,ISSFLG,IBOUND,HNEW,HOLD,FLOWJA,IA,
     1      BUFF,TOP,IOUT,NODES,NODLAY,IUNSTR,Sn,So,TOP,BOT,iunsat
      USE GWFBASMODULE,ONLY:MSUM,ICBCFL,VBVL,VBNM,DELT,PERTIM,TOTIM
      USE GWFBCFMODULE,ONLY:IBCFCB,LAYCON,SC1,SC2,IALTSTO
C
      CHARACTER*16 TEXT
      DOUBLE PRECISION STOIN,STOUT,SSTRG,SIN,SOUT,TLED,HSING,STRG,
     *  RHO,RHO1,RHO2,SNEW,SOLD,ONE,BBOT,TTOP,TOTTHICK,TP
C
      DATA TEXT /'         STORAGE'/
C     ------------------------------------------------------------------
C  
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
      IF(IBCFCB.GT.0) IBD=ICBCFL
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
      IF(LC.EQ.3 .OR. LC.EQ.2) KT=KT+1
      NNDLAY = NODLAY(K)
      NSTRT = NODLAY(K-1)+1
      DO 300 N=NSTRT,NNDLAY
C
C6------SKIP NO-FLOW AND CONSTANT-HEAD CELLS.
      IF(IBOUND(N).EQ.0) GO TO 300
      HSING=HNEW(N)
C-------COMBINED COMPUTATION FOR LAYCON=4 OR 5
      IF(LC.EQ.4.OR.LC.EQ.5) GO TO 333
C
C7-----CHECK LAYER TYPE TO SEE IF ONE STORAGE CAPACITY OR TWO.
      IF(LC.NE.3 .AND. LC.NE.2) GO TO 285
C
C7A----TWO STORAGE CAPACITIES.
      TP=TOP(N)
      RHO2=SC2(N)*TLED
      RHO1=SC1(N)*TLED
      SOLD=RHO2
      IF(HOLD(N).GT.TP) SOLD=RHO1
      SNEW=RHO2
      IF(HSING.GT.TP) SNEW=RHO1
      STRG=SOLD*(HOLD(N)-TP) + SNEW*TP - SNEW*HSING
      GO TO 288
C
C7B----ONE STORAGE CAPACITY.
  285 RHO=SC1(N)*TLED
      STRG=RHO*HOLD(N) - RHO*HSING
      GO TO 288
333   CONTINUE
C7C-----COMBINED COMPUTATION FOR LAYCON=4 OR 5
      HSING=HNEW(N)
      SOLD = So(N)
      BBOT=BOT(N)
      TTOP=TOP(N)
      TOTTHICK = TTOP - BBOT
      CALL SAT_THIK(N,HSING,TOTTHICK,BBOT,SNEW,K,TTOP)
      RHO2 = SC2(N) * TLED * (SOLD - SNEW) * TOTTHICK
C7C --------CHANGE SATN FOR ALTERNATE COMPRESSIBLE STORAGE FORMULATION       
      IF(IALTSTO.EQ.1) THEN  ! COMPRESSIBLE STORAGE KICKS ON ONLY NEAR TOP OF CELL
        EPS = TOTTHICK * 0.001     !SMOOTH AT 0.1 PERCENT OF CELL THICKNESS
        EX = HNEW(N) - (TTOP - EPS)
        IF(EX .LT. 0) THEN ! BELOW (TOP - SMOOTH) SO NO COMPRESSIBLE STORAGE
          SNEW = 0.0  
        ELSEIF(EX. LT. EPS) THEN ! WITHIN SMOOTHENING FUNCTION
          SNEW = -2.0/EPS**3 * EX**3 + 3.0/EPS**2 * EX**2
        ELSEIF(EX. GE. EPS) THEN  ! APPLY FULL COMPRESSIBLE STORAGE TERM    
          SNEW = 1.0
        ENDIF    
      ENDIF      
C      
      RHO1 = SC1(N) * TLED * SNEW * (HOLD(N) - HSING)
      STRG = RHO1 + RHO2
C
C8-----STORE CELL-BY-CELL FLOW IN BUFFER AND ADD TO ACCUMULATORS.
  288 BUFF(N)=STRG
      FLOWJA(IA(N)) = -STRG
      SSTRG=STRG
      IF(STRG.LT.ZERO) THEN
        STOUT=STOUT-SSTRG
      ELSE
        STOIN=STOIN+SSTRG
      END IF
C
  300 CONTINUE
C-------If TVM is used, call it to update BUFF, STOIN and STOUT for any
C-------transient storage changes
      IF(ITVMUNIT.GT.0) THEN
          CALL TVMU2BDS(KSTP,KPER,STOIN,STOUT)
      END IF
C7C4C------IF DENSITY DRIVEN FLOW THEN CORRECT FOR DENSITY TERMS HERE
          IF(IUDDF. GT.0) CALL DDF1BD(KSTP,KPER)       
C9------record contents of buffer for structured and unstructured grids
      IF(IUNSTR.EQ.0)THEN
C
C9A-----IF IBD FLAG IS SET RECORD THE CONTENTS OF THE BUFFER.
        IF(IBD.EQ.1) CALL UBUDSV(KSTP,KPER,TEXT,
     1                       IBCFCB,BUFF,NCOL,NROW,NLAY,IOUT)
        IF(IBD.EQ.2) CALL UBDSV1(KSTP,KPER,TEXT,IBCFCB,
     1            BUFF,NCOL,NROW,NLAY,IOUT,DELT,PERTIM,TOTIM,IBOUND)
      ELSE
C
C9B-----IF IBD FLAG IS SET RECORD THE CONTENTS OF THE BUFFER.
        IF(IBD.EQ.1) CALL UBUDSVU(KSTP,KPER,TEXT,IBCFCB,BUFF(1),NODES,
     1         IOUT,PERTIM,TOTIM)
        IF(IBD.EQ.2) CALL UBDSV1U(KSTP,KPER,TEXT,IBCFCB,BUFF(1),NODES,
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
C--------------------------------------------------------------------------------------
      SUBROUTINE GWF2BCFU1BDCH(KSTP,KPER)
C     ******************************************************************
C     COMPUTE FLOW FROM CONSTANT-HEAD CELLS
C     ******************************************************************
C
C     SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL, ONLY:NCOL,NROW,NLAY,IBOUND,HNEW,BUFF,ITRNSP,NOVFC,
     1 TOP,IOUT,NODES,NEQS,NODLAY,IA,JA,JAS,IUNSTR,IVC,ISYM,INCLN
      USE CLN1MODULE, ONLY: NCLNNDS,ICLNCB
      USE GWFBASMODULE,ONLY:MSUM,VBVL,VBNM,DELT,PERTIM,TOTIM,ICBCFL,
     1                      ICHFLG
      USE GWFBCFMODULE,ONLY:IBCFCB,LAYCON
      USE GWTBCTMODULE, ONLY: CBCH
      USE SMSMODULE, ONLY: AMATFL
C
      DOUBLE PRECISION HD,CHIN,CHOUT,XX1,TMP,RATE,CHCH1,HDIFF,
     *  X1,CIN,COUT
C     ------------------------------------------------------------------
C
C1------SET IBD TO INDICATE IF CELL-BY-CELL BUDGET VALUES WILL BE SAVED.
      IBD=0
      IF(IBCFCB.LT.0 .AND. ICBCFL.NE.0) IBD=-1
      IF(IBCFCB.GT.0) IBD=ICBCFL
      IF(ITRNSP.GT.0.AND.IBD.EQ.0) IBD = 999
C      IF(IBD.EQ.0) RETURN
C
C2------CLEAR BUDGET ACCUMULATORS.
      ZERO=0.
      CHIN=ZERO
      CHOUT=ZERO
      IBDLBL=0
C
C3------CLEAR BUFFER.
      DO 5 N=1,NEQS
        BUFF(N)=ZERO
5     CONTINUE
C
C4------LOOP THROUGH EACH CELL AND CALCULATE FLOW INTO MODEL FROM EACH
C4------CONSTANT-HEAD CELL.
      DO 200 N=1,NEQS
C
C5------IF CELL IS NOT CONSTANT HEAD SKIP IT & GO ON TO NEXT CELL.
        IF (IBOUND(N).GE.0)GO TO 200
C
C5A-------FIND LAYER NUMBER AND LAYCON FOR SUBSURFACE NODES
        IF(N.LE.NODES)THEN
          DO K=1,NLAY
            NNDLAY = NODLAY(K)
            NSTRT = NODLAY(K-1)+1
            LC=LAYCON(K)
            IF(N.GE.NSTRT.AND.N.LE.NNDLAY) GO TO 22
          ENDDO
22        CONTINUE
        ENDIF
C
C6------CLEAR VALUES FOR FLOW RATE THROUGH EACH FACE OF CELL.
        X1=ZERO
        CHCH1=ZERO
        RATE = ZERO
C
C7------CALCULATE FLOW THROUGH CONNECTING FACES.
        DO 30 II = IA(N)+1,IA(N+1)-1
          JJ = JA(II)
C
C7A-------IF THERE IS NO FLOW TO CALCULATE THROUGH THIS FACE, THEN GO ON
C7A-------TO NEXT FACE.  NO FLOW OCCURS TO AN
C7A-------ADJACENT NO-FLOW CELL, OR TO AN ADJACENT CONSTANT-HEAD CELL.
          IF(IBOUND(JJ).EQ.0) GO TO 30
          IF(IBOUND(JJ).LT.0 .AND. ICHFLG.EQ.0) GO TO 30
          IIS = JAS(II)
          IF(JJ.LT.N.AND.IVC(IIS).EQ.1)THEN  !VERTICAL DIRECTION UP
            HD = HNEW(N)
            IF(LC.NE.3 .AND. LC.NE.2) GO TO 122
            IF(NOVFC.EQ.1) GO TO 122
            TMP=HD
            IF(TMP.LT.TOP(N)) HD=TOP(N)
 122        HDIFF = HD - HNEW(JJ)
          ELSEIF(IVC(IIS).EQ.1)THEN !VERTICAL DIRECTION DOWN
            KK = K+1
            IF(NLAY.EQ.1) KK = K  
            HD=HNEW(JJ)
            IF(LAYCON(KK).NE.3 .AND. LAYCON(KK).NE.2) GO TO 123 !CSP 2 SHOULD BE 1?
            IF(NOVFC.EQ.1) GO TO 123
            TMP=HD
            IF(TMP.LT.TOP(JJ)) HD=TOP(JJ)
123         HDIFF=HNEW(N)-HD
          ELSE ! HORIZONTAL DIRECTION OR CONDUIT-NODE CONNECTION
            HDIFF = HNEW(N) - HNEW(JJ)
          ENDIF
C
C7B-------CALCULATE FLOW THROUGH THIS FACE INTO THE ADJACENT CELL.
          CHCH1= HDIFF*AMATFL(II)
          RATE = RATE + CHCH1
          IF(IBOUND(JJ).LT.0.AND. ICHFLG.EQ.0) GO TO 30
          X1 = CHCH1
          XX1=X1
C
C7C--------ACCUMULATE POSITIVE AND NEGATIVE FLOW.
          IF (X1.LT.ZERO) THEN
            CHOUT=CHOUT-XX1
          ELSE
            CHIN=CHIN+XX1
          END IF
   30   CONTINUE
C
C8--------STORE SUM IN BUFFER.
        BUFF(N)=RATE
        IF(ITRNSP.GT.0) CBCH(N) = RATE
C
  200 CONTINUE
C
C9------RETURN.
      RETURN
      END 
C    
       SUBROUTINE GWF2BCFU1BDCHWR(KSTP,KPER)
C     ******************************************************************
C     SAVE FLOW FROM CONSTANT-HEAD CELLS
C     ******************************************************************
C
C     SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL, ONLY:NCOL,NROW,NLAY,IBOUND,HNEW,BUFF,ITRNSP,NOVFC,
     1   TOP,IOUT,NODES,NEQS,NODLAY,IA,JA,JAS,IUNSTR,IVC,ISYM,INCLN,
     2   FLOWJA,FMBE
      USE GWFBASMODULE,ONLY:MSUM,VBVL,VBNM,DELT,PERTIM,TOTIM,ICBCFL,
     1                      ICHFLG
      USE GWFBCFMODULE,ONLY:IBCFCB,LAYCON
      USE GWTBCTMODULE, ONLY: CBCH
      USE SMSMODULE, ONLY: AMATFL
C
      CHARACTER*16 TEXT(1)
      DOUBLE PRECISION HD,CHIN,CHOUT,XX1,TMP,RATE,CHCH1,HDIFF,
     *  X1,CIN,COUT,ZERO
C
      DATA TEXT(1) /'   CONSTANT HEAD'/
C     ------------------------------------------------------------------
C
C1------SET IBD TO INDICATE IF CELL-BY-CELL BUDGET VALUES WILL BE SAVED.
      IBD=0
      IF(IBCFCB.LT.0 .AND. ICBCFL.NE.0) IBD=-1
      IF(IBCFCB.GT.0) IBD=ICBCFL
      IF(ITRNSP.GT.0.AND.IBD.EQ.0) IBD = 999
C      IF(IBD.EQ.0) RETURN
      CHIN = 0.0
      CHOUT = 0.0
C
C1A-----CLEAR BUFF
      DO N=1,NODES
          BUFF(N)=0.
      ENDDO
C
C--------------------------------------------------------------------------
C2------GWF DOMAIN
C---------------------------------------------------------------------------
      IF(IBD.EQ.2) THEN
C2A-----IF SAVING CELL-BY-CELL FLOW IN A LIST, COUNT CONSTANT-HEAD
C2A-----CELLS AND WRITE HEADER RECORDS.
        NCH=0
        DO 7 N=1,NODES
          IF(IBOUND(N).LT.0) NCH=NCH+1
7       CONTINUE
C2B-------WRITE HEADER FOR THE LIST
        IF(IUNSTR.EQ.0)THEN
          CALL UBDSV2(KSTP,KPER,TEXT(1),IBCFCB,NCOL,NROW,NLAY,
     1         NCH,IOUT,DELT,PERTIM,TOTIM,IBOUND)
        ELSE
          CALL UBDSV2U(KSTP,KPER,TEXT(1),IBCFCB,NODES,
     1         NCH,IOUT,DELT,PERTIM,TOTIM,IBOUND)
        ENDIF
      END IF
C
C3------LOOP THROUGH EACH CELL AND WRITE FLOW FROM EACH
C3------CONSTANT-HEAD CELL.
      IBDLBL = 0
      ZERO = 0.0
      DO 200 N=1,NODES
C
C4------IF CELL IS NOT CONSTANT HEAD SKIP IT & GO ON TO NEXT CELL.
        IF (IBOUND(N).GE.0)GO TO 200
C
C5--------ACCUMULATE INDIVIDUAL CHD FLOWS FROM FLOWJA ARRAY
          RATE=0.0
          DO II = IA(N),IA(N+1)-1
            X1=FLOWJA(II)
            XX1=X1
            IF(X1.LT.ZERO) THEN
              CHOUT=CHOUT-XX1
            ELSE
              CHIN=CHIN+XX1
            ENDIF
            RATE=RATE+X1
          ENDDO
C
C6--------PRINT THE FLOW FOR THE CELL IF REQUESTED.
        IF(IBD.LT.0) THEN
          IF(IBDLBL.EQ.0) WRITE(IOUT,899) TEXT(1),KPER,KSTP
  899     FORMAT(1X,/1X,A,'   PERIOD',I3,'   STEP',I3)
          IF(IUNSTR.EQ.0)THEN
            K = N / (NCOL*NROW) + 1
            IJ = N - (K-1)*NCOL*NROW
            I = (IJ-1)/NCOL + 1
            J = IJ - (I-1)*NCOL
            WRITE(IOUT,900) K,I,J,RATE
  900       FORMAT(1X,'LAYER',I3,'   ROW',I4,'   COL',I4,
     1       '   RATE',1PG15.6)
          ELSE
            WRITE(IOUT,910) N,RATE
  910       FORMAT(1X,'NODE',I8,'   RATE',1PG15.6)
          ENDIF
          IBDLBL=1
        END IF
C
C7------IF SAVING CELL-BY-CELL FLOW IN LIST, WRITE FLOW FOR CELL.
        IF(IBD.EQ.2)THEN
          SRATE = RATE
          IF(IUNSTR.EQ.0)THEN
            K = N / (NCOL*NROW) + 1
            IJ = N - (K-1)*NCOL*NROW
            I = (IJ-1)/NCOL + 1
            J = IJ - (I-1)*NCOL
            CALL UBDSVA(IBCFCB,NCOL,NROW,J,I,K,SRATE,IBOUND,NLAY)
          ELSE
            CALL UBDSVAU(IBCFCB,NODES,N,SRATE,IBOUND)
          ENDIF
C--------------------------------------------------------------
        ENDIF
C
C8--------STORE SUM IN BUFFER.
        BUFF(N)=RATE
        FMBE(N) = FMBE(N) + RATE
        IF(ITRNSP.GT.0) CBCH(N) = RATE
C
  200 CONTINUE
C--------------------------------------------------------------------------
C
C8-----IF SAVING CELL-BY-CELL FLOW IN 3-D ARRAY, WRITE THE ARRAY.
      IF(IBD.EQ.1)THEN
        IF(IUNSTR.EQ.0)THEN
          CALL UBUDSV(KSTP,KPER,TEXT(1),
     1                IBCFCB,BUFF,NCOL,NROW,NLAY,IOUT)
        ELSE
          CALL UBUDSVU(KSTP,KPER,TEXT(1),IBCFCB,BUFF,NODES,
     1                 IOUT,PERTIM,TOTIM)
        ENDIF
      ENDIF
C
C9------SAVE TOTAL CONSTANT HEAD FLOWS AND VOLUMES IN VBVL TABLE
C9------FOR INCLUSION IN BUDGET. PUT LABELS IN VBNM TABLE.
      CIN=CHIN
      COUT=CHOUT
      VBVL(1,MSUM)=VBVL(1,MSUM)+CIN*DELT
      VBVL(2,MSUM)=VBVL(2,MSUM)+COUT*DELT
      VBVL(3,MSUM)=CIN
      VBVL(4,MSUM)=COUT
      VBNM(MSUM)=TEXT(1)
      MSUM=MSUM+1
C
C10------RETURN.
      RETURN
      END
C
      SUBROUTINE GWF2BCFU1BDADJ(KSTP,KPER)
C     ******************************************************************
C     COMPUTE FLOW BETWEEN ADJACENT CELLS IN A SUBREGION OF THE GRID
C     ******************************************************************
C
C     SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL, ONLY:NCOL,NROW,NLAY,IBOUND,HNEW,BUFF,NODLAY,NEQS,
     1  TOP,IOUT,NODES,NJA,IA,JA,JAS,IUNSTR,IVC,ISYM,ITRNSP,issflg,
     1  Sn,So,INGNC,INGNC2,INGNCn,FLOWJA,NOVFC,iunsat
      USE GWFBASMODULE,ONLY:ICBCFL,DELT,PERTIM,TOTIM,ICHFLG
      USE GWFBCFMODULE,ONLY:IBCFCB,LAYCON
      USE GWTBCTMODULE, ONLY: CBCF
       USE SMSMODULE, ONLY: AMATFL
C
      DOUBLE PRECISION HD,TMP,HDIFF,HNEWN
C     ------------------------------------------------------------------
C
C1------IF CELL-BY-CELL FLOWS WILL BE SAVED IN A FILE, SET FLAG IBD.
C1------RETURN IF FLOWS ARE NOT BEING SAVED OR RETURNED.
      ZERO=0.
      IBD=0
      IF(IBCFCB.GT.0) IBD=ICBCFL
      IF(ITRNSP.GT.0.AND.IBD.EQ.0) IBD = 999
!      IF(IBD.EQ.0) RETURN
C
C3------FOR EACH CELL CALCULATE FLOW THRU ADJACENT FACE & STORE IN BUFFER.
      DO 200 K=1,NLAY
        LC=LAYCON(K)
        NNDLAY = NODLAY(K)
        NSTRT = NODLAY(K-1)+1
        DO 200 N=NSTRT,NNDLAY
C
C4------IF CELL IS NOT ACTIVE GO ON TO NEXT CELL.
        IF (IBOUND(N).EQ.0)GO TO 200
C
C5------CALCULATE FLOW THROUGH CONNECTING FACES.
        DO 30 II = IA(N)+1,IA(N+1)-1
          JJ = JA(II)
          IIS = JAS(II)
C
C5A-------IF THERE IS NO FLOW TO CALCULATE THROUGH THIS FACE, THEN GO ON
C5A-------TO NEXT FACE.
          IF(ICHFLG.EQ.0) THEN
            IF((IBOUND(N).LE.0) .AND. (IBOUND(JJ).LE.0)) GO TO 30
          ELSE
            IF(IBOUND(JJ).EQ.0) GO TO 30
          END IF
          HD=HNEW(JJ)
          HNEWN = HNEW(N)
          IF(IVC(IIS).EQ.1.AND.JJ.GT.N)THEN !VERTICAL DIRECTION DOWN
            KK = K+1
            IF(NLAY.EQ.1) KK = K  
            IF(LAYCON(KK).NE.3 .AND. LAYCON(KK).NE.2) GO TO 122 !CSP 2 SHOULD BE 1?
            TMP=HD
            IF(NOVFC.EQ.0)THEN
              IF(TMP.LT.TOP(JJ)) HD=TOP(JJ)
            ENDIF
          ENDIF
C ---------WHEN NOVFC IS ZERO AND VERTICAL DIRECTION IS UP
          IF(IVC(IIS).EQ.1.AND.JJ.LT.N)THEN !VERTICAL DIRECTION UP
            KK = K-1
            IF(NLAY.EQ.1) KK = K  
            IF(LAYCON(KK).NE.3 .AND. LAYCON(KK).NE.2) GO TO 122 !CSP 2 SHOULD BE 1?
            TMP=HNEWN
            IF(NOVFC.EQ.0)THEN
              IF(TMP.LT.TOP(N)) HNEWN=TOP(N)  
            ENDIF
          ENDIF           
C
C5B-------CALCULATE FLOW THROUGH THIS FACE INTO THE ADJACENT CELL.
  122     HDIFF=HNEWN-HD
          FLOWJA(II)= HDIFF*AMATFL(II)
   30   CONTINUE
C
  200 CONTINUE
c6----return
      RETURN
      END
      SUBROUTINE GWF2BCFU1BDADJWR(KSTP,KPER)
C     ******************************************************************
C     COMPUTE FLOW BETWEEN ADJACENT CELLS IN A SUBREGION OF THE GRID
C     FLOWS ARE WRITTEN SO THAT FLOW IS POSITIVE INTO A CELL.
C     ******************************************************************
C
C     SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL, ONLY:NCOL,NROW,NLAY,IBOUND,HNEW,BUFF,NODLAY,NEQS,
     1  TOP,IOUT,NODES,NJA,IA,JA,JAS,IUNSTR,IVC,ISYM,ITRNSP,issflg,
     1  Sn,So,INGNC,INGNC2,INGNCn,FLOWJA,JAFL,NJAG,iunsat
      USE GWFBASMODULE,ONLY:ICBCFL,DELT,PERTIM,TOTIM,ICHFLG
      USE GWFBCFMODULE,ONLY:IBCFCB,LAYCON
      USE GWTBCTMODULE, ONLY: CBCF
       USE SMSMODULE, ONLY: AMATFL
C
      CHARACTER*16 TEXT(4)
      DOUBLE PRECISION HD,TMP,HDIFF
      REAL, DIMENSION(:),ALLOCATABLE  :: FLOWGWS(:,:),FLOWJAG(:)
C
      DATA TEXT(1),TEXT(2),TEXT(3),TEXT(4)
     1 /'FLOW RIGHT FACE ','FLOW FRONT FACE ','FLOW LOWER FACE ',
     1  '   FLOW JA FACE '/
C     ------------------------------------------------------------------
C1-----RETURN IF BUDGETS ARE NOT REQUIRED
      IBD=0
      IF(IBCFCB.GT.0) IBD=ICBCFL
      IF(ITRNSP.GT.0.AND.IBD.EQ.0) IBD = 999
      IF(IBD.EQ.0) RETURN
      ZERO = 0.
C
C---------------------------------------------------------------------
C2-----FILL FLOW TERM INTO SYMMETRIC CBCF ARRAY IF TRANSPORT IS ACTIVE  
      IF(ITRNSP.GT.0)THEN
        DO N=1,NODES
          DO II = IA(N)+1,IA(N+1)-1
            JJ = JA(II)
            IIS = JAS(II)
            IF(JJ.LE.N) CYCLE ! FILL  UPPER TRIANGLE
            IF(IIS.GT.0) CBCF(IIS) = FLOWJA(II)
          ENDDO
        ENDDO
      ENDIF
C---------------------------------------------------------------------
      ITESTCBC = 0
      IF(ITESTCBC.EQ.1)THEN
        write(iout,*)' fluxes are below for all faces from a node'
        do i=1,nodes
          write(iout,66) (FLOWJA(ii),ii=ia(i),ia(i+1)-1)
        enddo
  66    format(10e15.6)
      ENDIF
C
C3------RECORD CONTENTS OF ARRAY FLOWGWS AND RETURN.
      IF(IUNSTR.EQ.0)THEN
        ALLOCATE(FLOWGWS(NODES,3))
        FLOWGWS = 0.0
        DO 300 K=1,NLAY
          LC=LAYCON(K)
          NNDLAY = NODLAY(K)
          NSTRT = NODLAY(K-1)+1
          DO 300 N=NSTRT,NNDLAY
C
C4---------CALCULATE FLOW THROUGH CONNECTING FACES FOR STRUCTURED GRID.
          DO 40 II = IA(N)+1,IA(N+1)-1
            JJ = JA(II)
            IF(JJ.LE.N) CYCLE
            IIS = JAS(II)
            IF(IVC(IIS).EQ.1)THEN !LOWER LAYER
              FLOWGWS(N,3) = FLOWJA(II)
            ELSEIF(JJ.EQ.N+NCOL)THEN !FORWARD FACE
              FLOWGWS(N,2) = FLOWJA(II)
            ELSEIF(JJ.EQ.N+1)THEN !RIGHT FACE
              FLOWGWS(N,1) = FLOWJA(II)
            ELSEIF(JJ.GT.NODES)THEN ! MATRIX-TO-CONDUIT CONNECTION
c             THIS CAN BE WRITTEN TO A SEPARATE FILE, BUT DIFFERENCE IN
c             FLOW FROM ONE CONDUIT-NODE TO THE NEXT IS FLOW TO/FROM MATRIX.
            ENDIF
 40       CONTINUE
 300    CONTINUE
C
        DO IDIR = 1,3
          IF(IBD.EQ.1) CALL UBUDSV(KSTP,KPER,TEXT(IDIR),IBCFCB,
     1     FLOWGWS(1,IDIR),NCOL,NROW,NLAY,IOUT)
          IF(IBD.EQ.2) CALL UBDSV1(KSTP,KPER,TEXT(IDIR),IBCFCB,
     1     FLOWGWS(1,IDIR),NCOL,NROW,NLAY,IOUT,DELT,PERTIM,TOTIM,IBOUND)
        ENDDO
        DEALLOCATE(FLOWGWS)
C
C5--------SAVE FLOW THROUGH CONNECTING FACES FOR UNSTRUCTURED GRID.
      ELSE
        IF(IBD.NE.0)THEN
C6------SAVE ONLY FOR GROUNDWATER NODES, SO COMPRESS OUT OTHER DOMAINS FROM FLOWJA
          ALLOCATE(FLOWJAG(NJAG))
            FLOWJAG = ZERO
C
          IJAG = 1
          DO N=1,NODES
            DO II = IA(N),IA(N+1)-1
              JJ = JA(II)
              JJG = JAFL(IJAG)
              IF(JJ.NE.JJG) CYCLE
              FLOWJAG(IJAG) = -FLOWJA(II)
              IJAG = IJAG + 1
            ENDDO
          ENDDO
C6B2------STORE UNSYMMETRIC CBC ARRAY IN FULL IA AND JA STRUCTURE FOR GW NODES
        IF(IBD.EQ.1)
     1   CALL UBUDSVU(KSTP,KPER,TEXT(4),IBCFCB,FLOWJAG(1),NJAG,IOUT,
     1         PERTIM,TOTIM)
        IF(IBD.EQ.2) CALL UBDSV1U(KSTP,KPER,TEXT(4),IBCFCB,FLOWJAG(1),
     1     NJAG,IOUT,DELT,PERTIM,TOTIM,IBOUND,NODES)
        ENDIF
        DEALLOCATE(FLOWJAG)
      ENDIF
c-----------------------------------------------------------------------------
c-----finally, for transport, if iunsat=1 and steady-state, get So from Sn
CC      iss = issflg(kper)
CCCCSP      if(itrnsp.gt.0.and.iunsat.eq.1.and.iss.eq.1)then
CC      if(itrnsp.gt.0.and.iss.eq.1)then
CC        do n=1,nodes
CC         So(N) = Sn(N)
CC        enddo
CC      endif
C
c7----return
      RETURN
      END
C------------------------------------------------------------------------
      SUBROUTINE GWF2BCFU1BDCLNWR(KSTP,KPER)
C     ******************************************************************
C     WRITE THE CLN-GW FLOWS TO THE GW BUDGET FILE
C     ******************************************************************
C
C      SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,        ONLY:NODES,NCOL,NROW,NLAY,IBOUND,BUFF,ITRNSP,
     1                        IOUT,FLOWJA,IA,JA,IUNSTR
      USE GWFBASMODULE,  ONLY:ICBCFL,DELT,PERTIM,TOTIM,ICHFLG
      USE GWFBCFMODULE,  ONLY:IBCFCB
      USE CLN1MODULE,    ONLY:NCLNGWC,ACLNNDS,ACLNGWC
      CHARACTER*16 TEXT
      LOGICAL FOUND
      DATA TEXT /'             CLN'/
C     ------------------------------------------------------------------
C
C1-----RETURN IF BUDGETS ARE NOT REQUIRED
      IBD=0
      IF(IBCFCB.GT.0) IBD=ICBCFL
      IF(ITRNSP.GT.0.AND.IBD.EQ.0) IBD = 999
      IF(IBD.EQ.0) RETURN    
C
C2-----INITIALIZE BUFF
      ZERO = 0.
      DO I=1,NODES
        BUFF(I) = ZERO
      ENDDO
C
C3-----WRITE HEADER FOR COMPACT BUDGET
      IF (IBD.EQ.2) THEN
        IF(IUNSTR.EQ.0)THEN
          CALL UBDSV2(KSTP,KPER,TEXT,IBCFCB,NCOL,NROW,NLAY,
     1          NCLNGWC,IOUT,DELT,PERTIM,TOTIM,IBOUND)
        ELSE 
          CALL UBDSV2U(KSTP,KPER,TEXT,IBCFCB,NODES,
     1          NCLNGWC,IOUT,DELT,PERTIM,TOTIM,IBOUND)
        ENDIF
      ENDIF
C
C4-----LOOP THROUGH EACH GW-CLN CONNECTION, RETRIEVE FLOW FROM
C4-----FLOWJA AND ACCUMULATE IN BUFF OR WRITE COMPACT BUDGET RECORD
      DO NN = 1,NCLNGWC
        IH = ACLNGWC(NN,1)
        ND1 = ACLNNDS(IH,1)
        N = ACLNGWC(NN,2)
        FOUND = .FALSE.
        DO II = IA(ND1)+1,IA(ND1+1)-1
          JJ = JA(II)
          IF(JJ.EQ.N) THEN
            FOUND = .TRUE.
            EXIT
          ENDIF
        ENDDO
        IF(.NOT.FOUND) CALL USTOP('error in GWF2BCFU1BDCLNWR')
        RATE=FLOWJA(II)
        IF(ICHFLG.EQ.0) THEN
          IF((IBOUND(ND1).LT.0) .AND. (IBOUND(JJ).LT.0)) RATE=ZERO
        END IF
        IF(IBOUND(ND1).EQ.0) RATE=ZERO
        IF(IBOUND(JJ).EQ.0) RATE=ZERO

        IF(IBD.EQ.1) THEN
          BUFF(N) = BUFF(N) + RATE
        ELSEIF(IBD.EQ.2) THEN
          IF(IUNSTR.EQ.0) THEN
            K = N / (NCOL*NROW) + 1
            IJ = N - (K-1)*NCOL*NROW
            I = (IJ-1)/NCOL + 1
            J = IJ - (I-1)*NCOL
            CALL UBDSVA(IBCFCB,NCOL,NROW,J,I,K,RATE,IBOUND,NLAY)
          ELSE
            CALL UBDSVAU(IBCFCB,NODES,N,RATE,IBOUND)
          ENDIF
        ENDIF
      ENDDO
C
C5-----WRITE CLN-GW FLOWS FOR NONCOMPACT BUDGET
      IF(IBD.EQ.1) THEN
        IF(IUNSTR.EQ.0) THEN
           CALL UBUDSV(KSTP,KPER,TEXT,
     1                 IBCFCB,BUFF,NCOL,NROW,NLAY,IOUT)
        ELSE
          CALL UBUDSVU(KSTP,KPER,TEXT,IBCFCB,BUFF,NODES,IOUT,
     1           PERTIM,TOTIM)
        ENDIF
      ENDIF
C
C6------RETURN.
      RETURN
      END
C------------------------------------------------------------------------
      SUBROUTINE GWF2BCFU1FMBE(KSTP,KPER)
C     ******************************************************************
C     COMPUTE FLOW MASS BALANCE FOR CBC AND CONSTANT HEAD FLOWS, AND STORAGE.
C     ******************************************************************
C
C     SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL, ONLY:NEQS,IOUT,NODES,NJA,IA,JA,JAS,IUNSTR,ISYM,ITRNSP,
     1  FLOWJA,JAFL,FMBE,NEQS
C     ------------------------------------------------------------------
C      
C1-----FILL FLOW MASS BALANCE ERROR (FMBE) ARRAY BY ADDING UP FLOWJAs.
C1----- STORAGE IS ALREADY FILLED IN FLOWJA DIAGONAL
        DO N=1,NEQS
          FMBE(N) = 0.0 ! INITIALIZE FMBE  
          DO II = IA(N),IA(N+1)-1
            FMBE(N) = FMBE(N) - FLOWJA(II) ! outflow is positive in flowja so flip sign
          ENDDO          
        ENDDO
C---------------------------------------------------------------------      
c7----return
      RETURN
      END      
C------------------------------------------------------------------------
      SUBROUTINE SGWF2BCFU1VCONDV(K)
C     ******************************************************************
C     STORE VERTICAL CONDUCTANCE TERM IN PGF. ALSO, IF KV IS READ THEN
C     COMPUTE LEAKANCE BETWEEN A LAYER AND THE NEXT LOWER LAYER.
C     ******************************************************************
C
C      SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,        ONLY:NCOL,NROW,NLAY,IBOUND,HNEW,DELR,DELC,
     1                        TOP,BOT,LAYCBD,IOUT,STRT,AREA,PGF,FAHL,
     1                        IUNSTR,NODLAY,IA,JA,JAS,IVC,ISYM,NOCVCO
      USE GWFBCFMODULE,ONLY:IBCFCB,IWDFLG,IWETIT,IHDWET,WETFCT,HDRY,CV,
     1                      LAYCON,LAYAVG,IKVFLAG,LAYAVGV
C
      DOUBLE PRECISION BBOT,TTOP,CONAREA,HHD,BEE
C     ------------------------------------------------------------------
C
      ZERO=0.
      HALF=0.5
C
C1------LOOP THROUGH ALL CELLS IN THE LAYER.
      NNDLAY = NODLAY(K)
      NSTRT = NODLAY(K-1)+1
      DO 100 N=NSTRT,NNDLAY
C2--------GO OVER CONNECTIONS OF NODE N
          DO II = IA(N)+1,IA(N+1)-1
            JJ = JA(II)
            IF(JJ.LE.N) CYCLE
            IIS = JAS(II)
C3------------FILL ONLY FOR VERTICAL CONNECTION
            IF(IVC(IIS).EQ.1)THEN !VERTICAL DIRECTION CONNECTION
              CONAREA = FAHL(IIS)
              IF(IKVFLAG.EQ.1)THEN
C4--------------COMPUTE LEAKANCE FROM K AND FIILL PGF
                  KK = K+1
C
C4A----------------CALCULATE VERTICAL HYDRAULIC CONDUCTIVITY FOR CELL.
                  HYC1=CV(N)
                  IF(HYC1.GT.ZERO) THEN
C4B------------------CALCULATE VERTICAL HYDRAULIC CONDUCTIVITY FOR ADJACENT CELL.
                    HYC2=CV(JJ)
                    IF(HYC2.GT.ZERO) THEN
C-------------------------------------------------------------------------------------                        
C4C-------------------IF LAYAVGV IS ZERO THEN COMPUTE HARMONIC AVERAGING AS IN MODFLOW
                      IF(LAYAVGV(K).EQ.0)THEN
C
C4D---------------------CALCULATE INVERSE LEAKANCE FOR CELL.
                        BBOT=BOT(N)
                        TTOP=TOP(N)
                        BOVK1=(TTOP-BBOT)*HALF/HYC1
C
C4D----------------------CALCULATE INVERSE LEAKANCE FOR ADJACENT CELL.
                        BBOT=BOT(JJ)
                        TTOP=TOP(JJ)
                        BOVK2= (TTOP-BBOT)*HALF/HYC2
C
C4E----------------------CALCULATE LEAKANCE
CSP                      CV(N) = 1.0 / (BOVK1 + BOVK2)
                        PGF(IIS) = CONAREA / (BOVK1 + BOVK2)
                      ELSE
C-------------------------------------------------------------------------------------                          
C4F---------------------COMPUTE ARITHMETIC AVERAGING AS IN FINITE ELEMENT CODES                          
C
C4G---------------------CALCULATE  LEAKANCE FOR CELL.
                        BBOT=BOT(N)
                        TTOP=TOP(N)
                        BOVK1 = HYC1 / ((TTOP-BBOT)*HALF)
C
C4H----------------------CALCULATE  LEAKANCE FOR ADJACENT CELL.
                        BBOT=BOT(JJ)
                        TTOP=TOP(JJ)
                        BOVK2 = HYC2 / ((TTOP-BBOT)*HALF)
C
C4E----------------------CALCULATE LEAKANCE
                        PGF(IIS) = CONAREA * (BOVK1 + BOVK2)         
                      ENDIF    
                    END IF !(HYC2.GT.ZERO)
                  END IF !(HYC1.GT.ZERO)
              ELSE ! IF(IKVFLAG.EQ.0)THEN
C5--------------CV IS LEAKANCE - FILL PGF DIRECTLY
                PGF(IIS) = CV(N) * CONAREA
              ENDIF !(IKVFLAG.EQ.0)
            END IF !(IVC(IIS).EQ.1)
          ENDDO
  100 CONTINUE
C
C6------RETURN.
      RETURN
      END
C---------------------------------------------------------------------------------
      SUBROUTINE READLPF(IN,INLAK)
C     ******************************************************************
C     ALLOCATE AND READ DATA FOR LAYER PROPERTY FLOW PACKAGE
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL, ONLY:NCOL,NROW,NLAY,ITRSS,LAYHDT,LAYHDS,LAYCBD,
     1            NCNFBD,IBOUND,BUFF,NBOTM,DELR,DELC,IOUT,NOVFC,itrnsp,
     2            NODES,IFREFM,IUNSTR,ICONCV,NOCVCO,NJA,NJAS,IWADI,
     3            IDEALLOC_LPF,HNEW,Sn,So,NODLAY,BOT,TOP,ARAD,iunsat,
     4            IOUTNORMAL
      USE GWFBCFMODULE,ONLY:IBCFCB,IWDFLG,IWETIT,IHDWET,WETFCT,HDRY,CV,
     1                 LAYCON,LAYAVG,SC1,SC2,WETDRY,CHANI,IHANISO,
     2                 IKCFLAG,LAYWET,ISFAC,ITHFLG,LAYAVGV,
     2                 LAYTYP,LAYVKA,LAYSTRT,alpha,beta,sr,brook,
     3                 LAYFLG,VKA,VKCB,HANI,HK,IBPN,BP,IDRY,IALTSTO,
     4   ITABRICH,INTRICH,IUZONTAB,RETCRVS,NUTABROWS,NUZONES 
C
      DOUBLE PRECISION HD,THCK,TTOP,BBOT,TOTTHICK
      CHARACTER*14 LAYPRN(5),AVGNAM(5),TYPNAM(3),VKANAM(2),WETNAM(2),
     1            HANNAM 
      DATA AVGNAM/'      HARMONIC','   LOGARITHMIC','   LOG+ARITHM ',
     1   '   ARITHMETIC ','FINITE ELEMENT'/
      DATA TYPNAM/'     CONFINED ','  CONVERTIBLE ','     UPSTREAM '/
      DATA VKANAM/'    VERTICAL K','    ANISOTROPY'/
      DATA WETNAM/'  NON-WETTABLE','      WETTABLE'/
      DATA HANNAM/'      VARIABLE'/
      CHARACTER*400 LINE
      CHARACTER*24 ANAME(10),STOTXT
      CHARACTER*4 PTYP
C
      DATA ANAME(1) /'   HYD. COND. ALONG ROWS'/
      DATA ANAME(2) /'  HORIZ. ANI. (COL./ROW)'/
      DATA ANAME(3) /'     VERTICAL HYD. COND.'/
      DATA ANAME(4) /' HORIZ. TO VERTICAL ANI.'/
      DATA ANAME(5) /'QUASI3D VERT. HYD. COND.'/
      DATA ANAME(6) /'        SPECIFIC STORAGE'/
      DATA ANAME(7) /'          SPECIFIC YIELD'/
      DATA ANAME(8) /'        WETDRY PARAMETER'/
      DATA ANAME(9) /'     STORAGE COEFFICIENT'/
      DATA ANAME(10) /'UNSAT PARAMETER ZONE MAP'/      
C
C     ------------------------------------------------------------------
C1------Allocate scalar data.
      ALLOCATE(ISFAC,ITHFLG)
      ZERO=0.
C
C2------IDENTIFY PACKAGE
      WRITE(IOUT,1) IN
    1 FORMAT(1X,/1X,'LPF -- LAYER-PROPERTY FLOW PACKAGE, VERSION 7',
     1', 5/2/2005',/,9X,'INPUT READ FROM UNIT ',I4)
C
C3------READ COMMENTS AND ITEM 1.
ccsp commented out to be backward compatible with LPF package of MF2K5.
ccsp      IF(IFREFM.EQ.0) THEN
ccsp        IF(IUNSTR.NE.0)THEN
ccsp          READ(IN,2)IBCFCB,HDRY,NPLPF,IKCFLAG
ccsp        ELSE
ccsp          READ(IN,2)IBCFCB,HDRY,NPLPF
ccsp        ENDIF
ccsp      ELSE

        CALL URDCOM(IN,IOUT,LINE)
        LLOC=1
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IBCFCB,R,IOUT,IN)
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,I,HDRY,IOUT,IN)
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NPLPF,R,IOUT,IN)
        IF(IUNSTR.NE.0)
     1    CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IKCFLAG,R,IOUT,IN)
ccsp      ENDIF
2     FORMAT(I10,F10.3,2I10)
C
C3A-----WRITE ITEM 1
      IF(IBCFCB.LT.0) WRITE(IOUT,8)
    8 FORMAT(1X,'CONSTANT-HEAD CELL-BY-CELL FLOWS WILL BE PRINTED',
     1  ' WHEN ICBCFL IS NOT 0')
      IF(IBCFCB.GT.0) WRITE(IOUT,9) IBCFCB
    9 FORMAT(1X,'CELL-BY-CELL FLOWS WILL BE SAVED ON UNIT ',I4)
      WRITE(IOUT,11) HDRY
   11 FORMAT(1X,'HEAD AT CELLS THAT CONVERT TO DRY=',1PG13.5)
      IF(NPLPF.GT.0) THEN
         WRITE(IOUT,15) NPLPF
   15    FORMAT(1X,I5,' Named Parameters     ')
      ELSE
         NPLPF=0
         WRITE(IOUT,'(A)') ' No named parameters'
      END IF
C
      IF(IUNSTR.NE.0)THEN
        IF(IKCFLAG.EQ.0)WRITE(IOUT,39)
   39   FORMAT(1X,'IKCFLAG=0, NODAL INPUT OF HY AND CV')
        IF(IKCFLAG.EQ.1)WRITE(IOUT,41)
   41   FORMAT(1X,'IKCFLAG=1, CONNECTIVITY INPUT OF HY',1X,
     1        '(OR TRAN FOR CONFINED) AND CV')
        IF(IKCFLAG.EQ.-1)WRITE(IOUT,44)
   44   FORMAT(1X,'IKCFLAG=-1, CONNECTIVITY INPUT OF CONDUCTANCE')
      ENDIF
C
C3B-----GET OPTIONS.
      ISFAC=0
      ICONCV=0
      ITHFLG=0
      NOCVCO=0
      NOPCHK=0
      STOTXT=ANAME(6)
   20 CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,I,R,IOUT,IN)
      IF(LINE(ISTART:ISTOP).EQ.'STORAGECOEFFICIENT') THEN
         ISFAC=1
         STOTXT=ANAME(9)
         WRITE(IOUT,21)
   21    FORMAT(1X,'STORAGECOEFFICIENT OPTION:',/,
     1     1X,'Read storage coefficient rather than specific storage')
      ELSE IF(LINE(ISTART:ISTOP).EQ.'CONSTANTCV') THEN
         ICONCV=1
         WRITE(IOUT,23)
   23    FORMAT(1X,'CONSTANTCV OPTION:',/,1X,'Constant vertical',
     1         ' conductance for convertible layers')
      ELSE IF(LINE(ISTART:ISTOP).EQ.'THICKSTRT') THEN
         ITHFLG=1
         WRITE(IOUT,25)
   25    FORMAT(1X,'THICKSTRT OPTION:',/,1X,'Negative LAYTYP indicates',
     1 ' confined layer with thickness computed from STRT-BOT')
      ELSE IF(LINE(ISTART:ISTOP).EQ.'NOCVCORRECTION') THEN
         NOCVCO=1
         WRITE(IOUT,27)
   27    FORMAT(1X,'NOCVCORRECTION OPTION:',/,1X,
     1    'Do not adjust vertical conductance when applying',
     2              ' the vertical flow correction')
      ELSE IF(LINE(ISTART:ISTOP).EQ.'NOVFC') THEN
         NOVFC=1
         IWADI = 0
         NOCVCO=1
         WRITE(IOUT,29)
   29    FORMAT(1X,'NOVFC OPTION:',/,1X,
     1    'Do not apply the vertical flow correction')
      ELSE IF(LINE(ISTART:ISTOP).EQ.'NOPARCHECK') THEN
         NOPCHK=1
         WRITE(IOUT,30)
   30    FORMAT(1X,'NOPARCHECK  OPTION:',/,1X,
     1    'For data defined by parameters, do not check to see if ',
     2        'parameters define data at all cells')
      ELSE IF(LINE(ISTART:ISTOP).EQ.'BUBBLEPT') THEN
         IBPN=1
         WRITE(IOUT,31)
   31    FORMAT(1X,'BUBBLEPT  OPTION:',/,1X,
     1    'For Richards equation the bubble point head is also ',
     2        'input at all cells')
         allocate(bp(nodes))
         DO N=1,NODES
           BP(NODES) = 0.0
         ENDDO
      ELSE IF(LINE(ISTART:ISTOP).EQ.'FULLYDRY') THEN
         IDRY=1
         WRITE(IOUT,32)
   32    FORMAT(1X,'FULLYDRY  OPTION:',/,1X,
     1    'For Richards equation, residual saturation is only applied ',
     2    'to relative permeability and soil saturation can vary from ',
     3    '0 to 1')   
       ELSE IF(LINE(ISTART:ISTOP).EQ.'TABRICH') THEN
         ITABRICH=1
         ALLOCATE(NUZONES,NUTABROWS)         
         WRITE(IOUT,33)
   33    FORMAT(1X,'TABULAR INPUT OPTION:',/,1X,
     1    'For Richards equation the retention ',
     2    'and relative permeability curves are provided ',
     3    'as tabular input for different soil types') 
         CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NUZONES,R,IOUT,IN)
         WRITE(IOUT,36) NUZONES
36       FORMAT(1X,'NUMBER OF SOIL TYPE ZONES = ', I10)
         CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NUTABROWS,R,IOUT,IN)
         WRITE(IOUT,37) NUTABROWS
37       FORMAT(1X,'NUMBER OF ROWS OF TABULAR INPUT = ', I10)          
      ELSE IF(LINE(ISTART:ISTOP).EQ.'INTRICH') THEN
         INTRICH=1
         WRITE(IOUT,34)
   34    FORMAT(1X,'INTEGRATED CURVE OPTION:',/,1X,
     1    'For Richards equation with tabular input, the retention ',
     2    'and relative permeability curves are integrated over the ',
     3    'grid block thickness')          
      ELSE IF(LINE(ISTART:ISTOP).EQ.'ALTSTO') THEN
         IALTSTO=1
         WRITE(IOUT,35)
   35    FORMAT(1X,'ALTERNATIVE STORAGE FORMULATION (IALTSTO = 1)')    
      END IF
C
      IF(LLOC.LT.200) GO TO 20
C
C-------READ STUFF TABULAR INPUT OF RICHARDS EQUATION   
      IF(ITABRICH.NE.0) THEN
C-------DIMENSION AND READ ZONE MAP          
        ALLOCATE(IUZONTAB(NODES))
        CALL U1DINT(IUZONTAB,ANAME(10),NODES,0,IN,IOUT)  
C-------DIMENSION AND READ TABLES        
        ALLOCATE(RETCRVS(3,NUTABROWS,NUZONES))
C      
        DO IZON = 1,NUZONES
          WRITE(IOUT,65) IZON  
65        FORMAT(/10X,'RETENTION TABLE NUMBER',I5/10X,27('-')/
     *  10X,'CAPILLARY HEAD',7X,'SATURATION',3X,'RELATIVE PERMEABILITY')
          DO ITROWS = 1,NUTABROWS  
            READ (IN,*) (RETCRVS(I,ITROWS,IZON),I=1,3)  
            WRITE(IOUT,66)(RETCRVS(I,ITROWS,IZON),I=1,3) 
66          FORMAT(10X,E14.6,3X,E14.6,7X,E14.6) 
          ENDDO  
        ENDDO    
      ENDIF   
C
C4------ALLOCATE AND READ LAYTYP, LAYAVG, CHANI, LAYVKA, LAYWET, LAYSTRT.
      ALLOCATE(LAYTYP(NLAY))
      ALLOCATE(CHANI(NLAY))
      ALLOCATE(LAYVKA(NLAY))
      ALLOCATE(LAYSTRT(NLAY))
      READ(IN,*) (LAYTYP(K),K=1,NLAY)
      READ(IN,*) (LAYAVG(K),K=1,NLAY)
      READ(IN,*) (CHANI(K),K=1,NLAY)
      READ(IN,*) (LAYVKA(K),K=1,NLAY)
      READ(IN,*) (LAYWET(K),K=1,NLAY)
C
C4A-----PRINT A TABLE OF VALUES FOR LAYTYP, LAYAVG, CHANI, LAYVKA, LAYWET.
      WRITE(IOUT,47)
   47 FORMAT(1X,/3X,'LAYER FLAGS:',/1X,
     1 'LAYER       LAYTYP        LAYAVG         CHANI ',
     2 '       LAYVKA        LAYWET',/1X,75('-'))
      DO 50 K=1,NLAY
      WRITE(IOUT,48) K,LAYTYP(K),LAYAVG(K),CHANI(K),LAYVKA(K),LAYWET(K)
   48 FORMAT(1X,I4,2I14,1PE14.3,2I14)
C-----------------------------------------------------------------------------------
C4A2 ---ALSO CONVERT LAVAVGV TO 1 (WEIGHTED ARITH AVERAGING FOR LPF)
C4A2----WHEN LAYAVG IS 4 FOR FINITE ELEMENT AVERAGING
      IF(LAYAVG(K).EQ.4)THEN
        LAYAVGV(K) = 1
      ENDIF      
C
C4B-----SET OPTIONS FOR BCF
      IF(LAYTYP(K).EQ.0)THEN
        LAYCON(K) = 0
      ELSEIF(LAYTYP(K).GT.0.AND.LAYTYP(K).LT.4)THEN
        LAYCON(K) = 3
      ELSEIF(LAYTYP(K).EQ.4)THEN
        LAYCON(K) = 4
      ELSEIF(LAYTYP(K).EQ.5)THEN
        LAYCON(K) = 5        
      ELSEIF(LAYTYP(K).LT.0)THEN
        IF(ITHFLG.EQ.1)THEN
          LAYCON(K) = 0
        ELSE
          LAYCON(K) = 3
        ENDIF
      ENDIF
      IDEALLOC_LPF = 1
      IF(INLAK.GT.0) IDEALLOC_LPF = 2  !DEALLOCATE ONLY AFTER LAK7U1RP IS DONE
      IF(ICONCV.EQ.0) IDEALLOC_LPF = 0 ! NEED LPF ARRAYS FOR VARIABLE CV OPTION
      IF(CHANI(1).LE.0) IDEALLOC_LPF = 0 ! NEED LPF ARRAYS FOR VARIABLE ANISOTROPY
C
C4C-----SET GLOBAL HEAD-DEPENDENT TRANSMISSIVITY AND STORAGE FLAGS.
      IF (LAYTYP(K).NE.0) THEN
        LAYHDT(K)=1
        LAYHDS(K)=1
      ELSE
        LAYHDT(K)=0
        LAYHDS(K)=0
      ENDIF
   50 CONTINUE
      
C
C4D-----SET LAYSTRT AND RESET LAYTYP IF THICKSTRT OPTION IS ACTIVE.
      DO 60 K=1,NLAY
      LAYSTRT(K)=0
      IF(LAYTYP(K).LT.0 .AND. ITHFLG.NE.0) THEN
         LAYSTRT(K)=1
         LAYTYP(K)=0
         LAYHDT(K)=0
         LAYHDS(K)=0
         WRITE(IOUT,57) K
   57    FORMAT(1X,'Layer',I5,
     1  ' is confined because LAYTYP<0 and THICKSTRT option is active')
      END IF
   60 CONTINUE
C
C4E-----SET HORIZONTAL ANISOTROPY FLAG
      IHANISO = 0
      DO K=1,NLAY
        IF(ABS(CHANI(K) - 1.0).GT.1.0E-6) IHANISO = 1
      ENDDO
      IF(IHANISO.EQ.1) THEN 
        ALLOCATE(IOUTNORMAL)  
        IOUTNORMAL = 0
        ALLOCATE(ARAD(NJAS))
      ENDIF  
C
C4F-----BASED ON LAYTYP, LAYAVG, CHANI, LAYWET, COUNT THE NUMBER OF EACH
C4F-----TYPE OF 2-D ARRAY; CHECK VALUES FOR CONSISTENCY; AND SETUP
C4F-----POINTERS IN LAYTYP, CHANI, AND LAYWET FOR CONVENIENT ACCESS
C4F-----TO SC2, HANI, and WETDRY.  PRINT INTERPRETED VALUES OF FLAGS.
      NCNVRT=0
      NHANI=0
      NWETD=0
      WRITE(IOUT,67)
   67 FORMAT(1X,/3X,'INTERPRETATION OF LAYER FLAGS:',/1X,
     1  '                       INTERBLOCK     HORIZONTAL',
     2  '    DATA IN',/1X,
     3  '        LAYER TYPE   TRANSMISSIVITY   ANISOTROPY',
     4  '   ARRAY VKA   WETTABILITY',/1X,
     5  'LAYER      (LAYTYP)      (LAYAVG)       (CHANI)',
     6  '      (LAYVKA)      (LAYWET)',/1X,75('-'))
      ILAYUNSAT = 0
      DO 100 K=1,NLAY
C
      IF(LAYTYP(K).EQ.5) ILAYUNSAT = 1 ! AT LEAST 1 LAYER USES RICHARDS EQUATION    
C          
      IF(LAYTYP(K).NE.0) THEN
         NCNVRT=NCNVRT+1
         LAYTYP(K)=NCNVRT
      END IF
      IF(CHANI(K).LE.ZERO) THEN
         NHANI=NHANI+1
         CHANI(K)=-NHANI
      END IF
      IF(LAYWET(K).NE.0) THEN
         IF(LAYTYP(K).EQ.0) THEN
            WRITE(IOUT,*)
     1          ' LAYWET is not 0 and LAYTYP is 0 for layer:',K
            WRITE(IOUT,*) ' LAYWET must be 0 if LAYTYP is 0'
            CALL USTOP(' ')
         ELSE
            NWETD=NWETD+1
            LAYWET(K)=NWETD
         END IF
      END IF
      IF(LAYAVG(K).LT.0 .OR. LAYAVG(K).GT.4) THEN
         WRITE(IOUT,74) LAYAVG(K)
   74    FORMAT(1X,I8,
     1    ' IS AN INVALID LAYAVG VALUE -- MUST BE 0, 1, 2, 3 or 4')
         CALL USTOP(' ')
      END IF
      LAYPRN(1)=TYPNAM(1)
      IF(LAYTYP(K).NE.0) LAYPRN(1)=TYPNAM(2)
      IF(LAYCON(K).EQ.4.OR.LAYCON(K).EQ.5) LAYPRN(1)=TYPNAM(3)
      LAYPRN(2)=AVGNAM(LAYAVG(K)+1)
      IF(LAYAVGV(K).NE.0) LAYPRN(2) = AVGNAM(5)
      IF(CHANI(K).LE.0) THEN
         LAYPRN(3)=HANNAM
      ELSE
         WRITE(LAYPRN(3),'(1PE14.3)') CHANI(K)
      END IF
      LAYPRN(4)=VKANAM(1)
      IF(LAYVKA(K).NE.0) LAYPRN(4)=VKANAM(2)
      LAYPRN(5)=WETNAM(1)
      IF(LAYWET(K).NE.0) LAYPRN(5)=WETNAM(2)
      WRITE(IOUT,78) K,(LAYPRN(I),I=1,5)
   78 FORMAT(1X,I4,5A)
  100 CONTINUE
C
C4G-----PRINT WETTING INFORMATION.
      IF(NWETD.EQ.0) THEN
         WRITE(IOUT,13)
   13    FORMAT(1X,/,1X,'WETTING CAPABILITY IS NOT ACTIVE IN ANY LAYER')
         IWDFLG=0
      ELSE
         WRITE(IOUT,12) NWETD
   12    FORMAT(1X,/,1X,'WETTING CAPABILITY IS ACTIVE IN',I4,' LAYERS')
         IWDFLG=1
         READ(IN,*) WETFCT,IWETIT,IHDWET
         IF(IWETIT.LE.0) IWETIT=1
         WRITE(IOUT,*) ' WETTING FACTOR=',WETFCT
         WRITE(IOUT,*) ' WETTING ITERATION INTERVAL=',IWETIT
         WRITE(IOUT,*) ' IHDWET=',IHDWET
      END IF
C
C5------ALLOCATE MEMORY FOR ARRAYS.
      ALLOCATE(CV(NODES))
C
      ALLOCATE(LAYFLG(6,NLAY))
      ALLOCATE(HK(NODES))
      ALLOCATE(VKA(NODES))
      IF(NCNFBD.GT.0) THEN
         ALLOCATE(VKCB(NODES))
      ELSE
         ALLOCATE(VKCB(1))
      END IF
      IF(ITRSS.NE.0) THEN
         ALLOCATE(SC1(NODES))
      ELSE
         ALLOCATE(SC1(1))
      END IF
      IF(ITRSS.NE.0 .AND. NCNVRT.GT.0 .or. itrnsp.ne.0) THEN
         ALLOCATE(SC2(NODES))
         SC2 = 0.0
      ELSE
         ALLOCATE(SC2(1))
      END IF
      IF(NHANI.GT.0) THEN
         ALLOCATE(HANI(NODES))
      ELSE
         ALLOCATE(HANI(1))
      END IF
      IF(NWETD.GT.0) THEN
         ALLOCATE(WETDRY(NODES))
      ELSE
         ALLOCATE(WETDRY(1))
      END IF
C      
      if(iunsat.eq.1.OR.ILAYUNSAT.EQ.1)then
        if(itabrich.eq.0)then  
          ALLOCATE(alpha(NODES),beta(nodes),sr(nodes),brook(nodes))
          DO N=1,NODES
            ALPHA(N) =0.0
            BETA(N) = 1.0
            SR(N) = 0.999
            BROOK(N) = 1.0
          ENDDO
        endif
        IF(IUNSAT.EQ.1)THEN
          DO K=1,NLAY
            IF(LAYCON(K).NE.0) LAYCON(K)=5
          ENDDO 
        ENDIF  
      ENDIF     
C
C6------READ PARAMETER DEFINITIONS
      NPHK=0
      NPVKCB=0
      NPVK=0
      NPVANI=0
      NPSS=0
      NPSY=0
      NPHANI=0
      IF(NPLPF.GT.0) THEN
         WRITE(IOUT,115)
  115    FORMAT(/,' PARAMETERS DEFINED IN THE LPF PACKAGE')
         DO 120 K=1,NPLPF
         CALL UPARARRRP(IN,IOUT,N,1,PTYP,1,0,-1)
C   Note that NPHK and the other NP variables in
C   this group are used only as flags, not counts
         IF(PTYP.EQ.'HK') THEN
            NPHK=1
         ELSE IF(PTYP.EQ.'HANI') THEN
C6A-----WHEN A HANI PARAMETER IS USED, THEN ALL HORIZONTAL ANISOTROPY
C6A-----MUST BE DEFINED USING PARAMETERS.  ENSURE THAT ALL CHANI <= 0
            DO 118 I = 1, NLAY
              IF (CHANI(I).GT.0.0) THEN
                WRITE(IOUT,117)
  117           FORMAT(/,
     &' ERROR: WHEN A HANI PARAMETER IS USED, CHANI FOR ALL LAYERS',/,
     &' MUST BE LESS THAN OR EQUAL TO 0.0 -- STOP EXECUTION',
     &' (GWF2LPFU1AR)')
                CALL USTOP(' ')
              ENDIF
  118       CONTINUE
            NPHANI=1
         ELSE IF(PTYP.EQ.'VKCB') THEN
            NPVKCB=1
         ELSE IF(PTYP.EQ.'VK') THEN
            NPVK=1
            CALL SGWF2LPFU1CK(IOUT,N,'VK  ')
         ELSE IF(PTYP.EQ.'VANI') THEN
            NPVANI=1
            CALL SGWF2LPFU1CK(IOUT,N,'VANI')
         ELSE IF(PTYP.EQ.'SS') THEN
            NPSS=1
         ELSE IF(PTYP.EQ.'SY') THEN
            NPSY=1
         ELSE
            WRITE(IOUT,*) ' Invalid parameter type for LPF Package'
            CALL USTOP(' ')
         END IF
  120    CONTINUE
      END IF
C
C7------READ PARAMETERS AND CONVERT FOR UNSTRUCTURED AND STRUCTURED GRIDS
      IF(IUNSTR.EQ.0) THEN
        CALL SGWF2LPFU1S(IN,NPHK,NPHANI,NPVK,NPVANI,NPSS,NPSY,NPVKCB,
     *    STOTXT,NOPCHK)
      ELSE
        CALL SGWF2LPFU1G(IN,NPHK,NPHANI,NPVK,NPVANI,NPSS,NPSY,NPVKCB,
     *    STOTXT,NOPCHK)
      ENDIF
C
C--------------------------------------------------------------------------------
C8------SET INITIAL  GRID-BLOCK SATURATED THICKNESS FRACTIONS AND TRANSMISSIVITY WHEN NEEDED
      DO K=1,NLAY
        IF(LAYCON(K).EQ.4.OR.LAYCON(K).EQ.5) THEN
C8A-------SET INITIAL SATURATED GRID-BLOCK FRACTIONS FOR LAYCON=4
          NNDLAY = NODLAY(K)
          NSTRT = NODLAY(K-1)+1
          DO N=NSTRT,NNDLAY
            IF(IBOUND(N).NE.0) THEN
C-------------CALCULATE SATURATED THICKNESS/TOTAL THICKNESS.
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
C--------------------------------------------------------------------
C9------SET CONSTANT TERMS IN PGF ARRAY IF IT IS NOT READ DIRECTLY
      IF(IKCFLAG.EQ.0)THEN
C9A--------CHECK CV CONSISTENCY
        CALL SGWF2LPFU1N
C
C10--------FILL PGF ARRAY
C
C10A--------FILL VERTICAL TERMS INTO PGF
        IF(NLAY.GT.1) CALL SGWF2LPFU1VCOND
C
C10B------FILL HORIZONTAL TERMS INTO PGF - HY FOR LAYCON 4 AND T FOR LAYCONS 0 OR 2
        CALL FILLPGFH
      ENDIF
C
C-----------------------------------------------------------------------------------
C11------SET UP STORAGE CAPACITIES FROM COEFFICIENTS
      IF(ITRSS.NE.0)THEN
        IF(ISFAC.EQ.0) THEN
          CALL SGWF2LPFU1SC(SC1(1),1)
        ELSE
          CALL SGWF2LPFU1SC(SC1(1),0)
        END IF
        IF(NCNVRT.GT.0) THEN
          CALL SGWF2LPFU1SC(SC2(1),0)
        ENDIF
      END IF
C
C--------------------------------------------------------------------------------
C12-----DEALLOCATE UNWANTED ARRAYS
csp need cv for merging Kv with boundary leakance for RIV      DEALLOCATE(CV)
C------NEED HK FOR CONDUIT CELLS SO KEEP
c      ILAYCON13=0
c      DO I=1,NLAY
c        IF(LAYCON(I).EQ.1.OR.LAYCON(I).EQ.3)ILAYCON13=1
c      ENDDO
c      IF(ILAYCON13.EQ.0)THEN
c        DEALLOCATE(HK)
c      ENDIF
C13-----RETURN
      RETURN
      END
C-------------------------------------------------------------------------
      SUBROUTINE SGWF2LPFU1S(IN,NPHK,NPHANI,NPVK,NPVANI,NPSS,NPSY,
     *  NPVKCB,STOTXT,NOPCHK)
C     ******************************************************************
C     ALLOCATE AND READ DATA FOR LAYER PROPERTY FLOW PACKAGE FOR STRUCTURED GRID
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL, ONLY:NCOL,NROW,NLAY,ITRSS,LAYHDT,LAYHDS,LAYCBD,ISYM,
     1            NCNFBD,IBOUND,BUFF,NBOTM,DELR,DELC,IOUT,NODES,IVC,
     2            IFREFM,IUNSTR,IA,JA,JAS,NJA,NJAS,ARAD,IPRCONN,iunsat
      USE GWFBCFMODULE,ONLY:IBCFCB,IWDFLG,IWETIT,IHDWET,WETFCT,HDRY,CV,
     1                      LAYCON,LAYAVG,HK,SC1,SC2,WETDRY,
     2                      IKCFLAG,laywet,ISFAC,ITHFLG,
     2                      LAYTYP,CHANI,LAYVKA,LAYSTRT,
     3                      LAYFLG,VKA,VKCB,HANI,IHANISO,
     4                      alpha,beta,sr,brook,BP,IBPN,itabrich
C
      REAL, DIMENSION(:,:),ALLOCATABLE  ::TEMP
      REAL, DIMENSION (:), ALLOCATABLE :: TEMPPL
C
      CHARACTER*24 ANAME(15),STOTXT
      CHARACTER*4 PTYP
C
      DATA ANAME(1) /'   HYD. COND. ALONG ROWS'/
      DATA ANAME(2) /'  HORIZ. ANI. (COL./ROW)'/
      DATA ANAME(3) /'     VERTICAL HYD. COND.'/
      DATA ANAME(4) /' HORIZ. TO VERTICAL ANI.'/
      DATA ANAME(5) /'QUASI3D VERT. HYD. COND.'/
      DATA ANAME(6) /'        SPECIFIC STORAGE'/
      DATA ANAME(7) /'          SPECIFIC YIELD'/
      DATA ANAME(8) /'        WETDRY PARAMETER'/
      DATA ANAME(9) /'     STORAGE COEFFICIENT'/
      DATA ANAME(10) /'        WETDRY PARAMETER'/
      DATA ANAME(11) /'                   alpha'/
      DATA ANAME(12) /'                    beta'/
      DATA ANAME(13) /'                      sr'/
      DATA ANAME(14) /'                   brook'/  
      DATA ANAME(15) /'     BUBBLING POINT HEAD'/
      REAL PI
C     ------------------------------------------------------------------
C1-------ALLOCATE TEMP ARRAY FOR STORING 3-D INFORMATION
      ALLOCATE(TEMP(NCOL,NROW))
      ZERO = 0.0
C2------SET ANGLE INTO ARAD WHEN THERE IS HORIZONTAL ANISOTROPY
      IF(IHANISO.EQ.1)THEN
        PI = 3.1415926536
C2A----SET FACE ANGLES IN ARAD
        DO N=1,NODES
          DO II = IA(N)+1,IA(N+1)-1
            JJ = JA(II)
            IF(JJ.GE.N) CYCLE
            IIS = JAS(II)
            IF(IVC(IIS).EQ.1) CYCLE 
            IF((N-JJ).EQ.1) THEN
              ARAD(IIS) = pi
            ELSEIF((JJ-N).EQ.1) THEN 
              ARAD = 0  
            ELSEIF(JJ .LT. N) THEN
              ARAD(IIS) =  pi/2.
            ELSE
              ARAD(IIS) = -pi/2. 
            ENDIF
          ENDDO
        ENDDO
C
C2B-------WRITE FACE ANGLES ARRAY
        IF(IPRCONN.NE.0)THEN
          WRITE(IOUT,*)'FACE ANGLE IS BELOW, 22G15.6, UNSYMMETRIC'
          ALLOCATE(TEMPPL(NJA))
          DO N=1,NODES
          DO II = IA(N)+1,IA(N+1)-1
            JJ = JA(II)
            IF(JJ.GE.N)THEN
              IIS = JAS(II)
              TEMPPL(II) = ARAD(IIS)
              TEMPPL(ISYM(II)) = ARAD(IIS)
            ENDIF
          ENDDO
          ENDDO
          WRITE(IOUT,55)(TEMPPL(J),J=1,NJA)
55      FORMAT(1P,22G15.6)
CSP          WRITE(IOUT,55)(ARAD(J),J=1,NJAS) !COMMENTED OUT SYMMETRIC WRITE
          DEALLOCATE (TEMPPL)
        ENDIF
      ENDIF
C
C3------DEFINE DATA FOR EACH LAYER -- VIA READING OR NAMED PARAMETERS.
      DO 200 K=1,NLAY
      KK=K
C
C3A-----DEFINE HORIZONTAL HYDRAULIC CONDUCTIVITY (HK)
      IF(NPHK.EQ.0) THEN
         CALL U2DREL(TEMP(1,1),ANAME(1),NROW,NCOL,KK,IN,IOUT)
      ELSE
         READ(IN,*) LAYFLG(1,K)
         WRITE(IOUT,121) ANAME(1),K,LAYFLG(1,K)
  121    FORMAT(1X,/1X,A,' FOR LAYER',I4,
     1   ' WILL BE DEFINED BY PARAMETERS',/1X,'(PRINT FLAG=',I4,')')
         CALL UPARARRSUB1(TEMP(1,1),NCOL,NROW,KK,'HK',
     1      IOUT,ANAME(1),LAYFLG(1,KK))
         IF(NOPCHK.EQ.0) CALL UPARARRCK(BUFF(1),IBOUND,IOUT,K,NCOL,
     1    NLAY,NROW,IUNSTR,'HK  ')
      END IF
      DO I=1,NROW
      DO J=1,NCOL
        N=J+(I-1)*NCOL+(K-1)*NROW*NCOL
        HK(N) = TEMP(J,I)
      ENDDO
      ENDDO
C
C3B-----READ HORIZONTAL ANISOTROPY IF CHANI IS NON-ZERO
      IF(CHANI(K).LE.ZERO) THEN
        KHANI=-CHANI(K)
        IF(NPHANI.EQ.0) THEN
           CALL U2DREL(TEMP(1,1),ANAME(2),NROW,NCOL,KK,IN,IOUT)
        ELSE
           READ(IN,*) LAYFLG(6,K)
           WRITE(IOUT,121) ANAME(2),K,LAYFLG(6,K)
           CALL UPARARRSUB1(TEMP(1,1),NCOL,NROW,KK,'HANI',
     &      IOUT,ANAME(2),LAYFLG(6,KK))
           IF(NOPCHK.EQ.0)CALL UPARARRCK(BUFF(1),IBOUND,IOUT,K,NCOL,
     &      NLAY,NROW,IUNSTR,'HANI')
        END IF
        DO I=1,NROW
        DO J=1,NCOL
          N=J+(I-1)*NCOL+(K-1)*NROW*NCOL
          HANI(N) = TEMP(J,I)
        ENDDO
        ENDDO
      END IF
C
C3C-----DEFINE VERTICAL HYDRAULIC CONDUCTIVITY OR HORIZONTAL TO VERTICAL
C3C-----ANISOTROPY (VKA).
      IANAME=3
      PTYP='VK'
      IF(LAYVKA(K).NE.0) THEN
         IANAME=4
         PTYP='VANI'
      END IF
      IF(NPVK.EQ.0 .AND. NPVANI.EQ.0) THEN
         CALL U2DREL(TEMP(1,1),ANAME(IANAME),NROW,NCOL,KK,IN,IOUT)
      ELSE
         READ(IN,*) LAYFLG(2,K)
         WRITE(IOUT,121) ANAME(IANAME),K,LAYFLG(2,K)
         CALL UPARARRSUB1(TEMP(1,1),NCOL,NROW,KK,PTYP,IOUT,
     &                       ANAME(IANAME),LAYFLG(2,KK))
         IF(NOPCHK.EQ.0)CALL UPARARRCK(BUFF(1),IBOUND,IOUT,K,NCOL,
     &     NLAY,NROW,IUNSTR,PTYP)
      END IF
      DO I=1,NROW
      DO J=1,NCOL
        N=J+(I-1)*NCOL+(K-1)*NROW*NCOL
        VKA(N) = TEMP(J,I)
      ENDDO
      ENDDO
C
C3D-----DEFINE SPECIFIC STORAGE OR STORAGE COEFFICIENT IN ARRAY SC1 IF TRANSIENT.
      IF(ITRSS.NE.0) THEN
         IF(NPSS.EQ.0) THEN
            CALL U2DREL(TEMP(1,1),STOTXT,NROW,NCOL,KK,IN,IOUT)
         ELSE
            READ(IN,*) LAYFLG(3,K)
            WRITE(IOUT,121) STOTXT,K,LAYFLG(3,K)
            CALL UPARARRSUB1(TEMP(1,1),NCOL,NROW,KK,'SS',
     1           IOUT,STOTXT,LAYFLG(3,KK))
            IF(NOPCHK.EQ.0) CALL UPARARRCK(BUFF(1),IBOUND(1),IOUT,K,
     1       NCOL,NLAY,NROW,IUNSTR,'SS  ')
         END IF
         DO I=1,NROW
         DO J=1,NCOL
          N=J+(I-1)*NCOL+(K-1)*NROW*NCOL
          SC1(N) = TEMP(J,I)
         ENDDO
         ENDDO
      END IF
C
C3E-----DEFINE SPECIFIC YIELD IN ARRAY SC2 IF TRANSIENT AND LAYER IS
C3E-----IS CONVERTIBLE.
      IF(LAYTYP(K).NE.0) THEN
         IF(ITRSS.NE.0) THEN
            IF(NPSY.EQ.0) THEN
               CALL U2DREL(TEMP(1,1),ANAME(7),NROW,NCOL,KK,IN,
     1                 IOUT)
            ELSE
               READ(IN,*) LAYFLG(4,K)
               WRITE(IOUT,121) ANAME(7),K,LAYFLG(4,K)
               CALL UPARARRSUB1(TEMP(1,1),NCOL,
     1         NROW,KK,'SY',IOUT,ANAME(7),LAYFLG(4,KK))
               IF(NOPCHK.EQ.0) CALL UPARARRCK(BUFF(1),IBOUND(1),IOUT,K,
     1          NCOL,NLAY,NROW,IUNSTR,'SY  ')
            END IF
            DO I=1,NROW
            DO J=1,NCOL
              N=J+(I-1)*NCOL+(K-1)*NROW*NCOL
              SC2(N) = TEMP(J,I)
            ENDDO
            ENDDO
         END IF
      END IF
C
C3F-----READ CONFINING BED VERTICAL HYDRAULIC CONDUCTIVITY (VKCB)
      IF(LAYCBD(K).NE.0) THEN
         IF(NPVKCB.EQ.0) THEN
            CALL U2DREL(TEMP(1,1),ANAME(5),NROW,NCOL,KK,IN,
     1             IOUT)
         ELSE
            READ(IN,*) LAYFLG(5,K)
            WRITE(IOUT,121) ANAME(5),K,LAYFLG(5,K)
            CALL UPARARRSUB1(TEMP(1,1),NCOL,NROW,KK,
     1         'VKCB',IOUT,ANAME(5),LAYFLG(5,KK))
            IF(NOPCHK.EQ.0) CALL UPARARRCK(BUFF(1),IBOUND(1),IOUT,K,
     1       NCOL,NLAY,NROW,IUNSTR,'VKCB')
         END IF
         DO I=1,NROW
         DO J=1,NCOL
           N=J+(I-1)*NCOL+(K-1)*NROW*NCOL
           VKCB(N) = TEMP(J,I)
         ENDDO
         ENDDO
      END IF
C
C3G-----READ WETDRY CODES IF WETTING CAPABILITY HAS BEEN INVOKED
C3G-----(LAYWET NOT 0).
      IF(LAYWET(K).NE.0) THEN
         CALL U2DREL(TEMP(1,1),ANAME(8),NROW,NCOL,KK,IN,IOUT)
         DO I=1,NROW
         DO J=1,NCOL
           N=J+(I-1)*NCOL+(K-1)*NROW*NCOL
           WETDRY(N) = TEMP(J,I)
         ENDDO
         ENDDO
      END IF
C4----read parameters for Richards Equation if solving unsaturated zone flow
      IF(LAYCON(K).NE.5)GOTO 300
      IF(ITABRICH.EQ.0) THEN
C3H-----READ alpha, beta, sr, brook
      CALL U2DREL(TEMP(1,1),ANAME(11),NROW,NCOL,KK,IN,IOUT)
      DO I=1,NROW
      DO J=1,NCOL
        N=J+(I-1)*NCOL+(K-1)*NROW*NCOL
        alpha(N) = TEMP(J,I)
      ENDDO
      ENDDO
      CALL U2DREL(TEMP(1,1),ANAME(12),NROW,NCOL,KK,IN,IOUT)
      DO I=1,NROW
      DO J=1,NCOL
        N=J+(I-1)*NCOL+(K-1)*NROW*NCOL
        beta(N) = TEMP(J,I)
      ENDDO
      ENDDO
      CALL U2DREL(TEMP(1,1),ANAME(13),NROW,NCOL,KK,IN,IOUT)
      DO I=1,NROW
      DO J=1,NCOL
        N=J+(I-1)*NCOL+(K-1)*NROW*NCOL
        sr(N) = TEMP(J,I)
      ENDDO
      ENDDO
      CALL U2DREL(TEMP(1,1),ANAME(14),NROW,NCOL,KK,IN,IOUT)
      DO I=1,NROW
      DO J=1,NCOL
        N=J+(I-1)*NCOL+(K-1)*NROW*NCOL
        brook(N) = TEMP(J,I)
      ENDDO
      ENDDO
      IF(IBPN.GT.0)THEN
        CALL U2DREL(TEMP(1,1),ANAME(15),NROW,NCOL,KK,IN,IOUT)
        DO I=1,NROW
        DO J=1,NCOL
          N=J+(I-1)*NCOL+(K-1)*NROW*NCOL
          bP(N) = TEMP(J,I)
        ENDDO
        ENDDO
      ENDIF
C6------read tabular input for retention and relative permeability curves      
      ELSE


      ENDIF      
  300 CONTINUE
C-------------------------------------------------------------------------      
  200 CONTINUE
C-------------------------------------------------------------------------
      DEALLOCATE(TEMP)
C
C4------RETURN
      RETURN
      END
C-------------------------------------------------------------------------
      SUBROUTINE SGWF2LPFU1G(IN,NPHK,NPHANI,NPVK,NPVANI,NPSS,NPSY,
     *  NPVKCB,STOTXT,NOPCHK)
C     ******************************************************************
C     ALLOCATE AND READ DATA FOR LAYER PROPERTY FLOW PACKAGE FOR UNSTRUCTURED (GENERAL) GRID
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,      ONLY:NCOL,NROW,NLAY,ITRSS,LAYHDT,LAYHDS,LAYCBD,
     1                 NCNFBD,IBOUND,BUFF,NBOTM,DELR,DELC,IOUT,ARAD,JAS,
     2                 NODES,IFREFM,IUNSTR,PGF,NJA,NJAS,NJAG,TOP,BOT,
     3                 NODLAY,IA,JA,IDSYMRD,IATMP,NJATMP,IVC,iunsat,cl1,
     4                 IOUTNORMAL 
      USE GWFBCFMODULE,ONLY:IBCFCB,IWDFLG,IWETIT,IHDWET,WETFCT,HDRY,CV,
     1                      LAYCON,LAYAVG,SC1,SC2,WETDRY,
     2                      IKCFLAG,laywet,ISFAC,ITHFLG,
     2                      LAYTYP,CHANI,LAYVKA,LAYSTRT,
     3                      LAYFLG,VKA,VKCB,HANI,HK,IHANISO,
     4                      alpha,beta,sr,brook,BP,IBPN,ITABRICH
C
      REAL, DIMENSION(:),ALLOCATABLE  ::TEMP
C
      CHARACTER*24 ANAME(16),STOTXT
      CHARACTER*4 PTYP
C
      DATA ANAME(1) /'   HYD. COND. ALONG ROWS'/
      DATA ANAME(2) /'  HORIZ. ANI. (COL./ROW)'/
      DATA ANAME(3) /'     VERTICAL HYD. COND.'/
      DATA ANAME(4) /' HORIZ. TO VERTICAL ANI.'/
      DATA ANAME(5) /'QUASI3D VERT. HYD. COND.'/
      DATA ANAME(6) /'        SPECIFIC STORAGE'/
      DATA ANAME(7) /'          SPECIFIC YIELD'/
      DATA ANAME(8) /'        WETDRY PARAMETER'/
      DATA ANAME(9) /'     STORAGE COEFFICIENT'/
      DATA ANAME(10) /' CONNECTION CONDUCTIVITY'/
      DATA ANAME(11) /'              FACE ANGLE'/
      DATA ANAME(12) /'                   alpha'/
      DATA ANAME(13) /'                    beta'/
      DATA ANAME(14) /'                      sr'/
      DATA ANAME(15) /'                   brook'/ 
      DATA ANAME(16) /'     BUBBLING POINT HEAD'/
C     ------------------------------------------------------------------
C
        ZERO=0.
C1------READ FACE ANGLES IF ANISOTROPIC
      IF(IHANISO.EQ.1)THEN
        CALL U1DRELNJA(ARAD(1),IATMP,ANAME(11),NJATMP,IN,IOUT,IDSYMRD)
C---------
C8A.......CHANGE OUTWARD NORMAL FOR COMPATIBIITY WITH OLDER CROSS-DISPERSION FORMULATION
C --------OLDER FORMULATION DID DIFFERENTLY AND Y WAS + UPWARD IN MANY APPLICATIONS         
c        IF(IOUTNORMAL. EQ.0) THEN 
c          DO N=1,NODES
c            IF(IBOUND(N).EQ.0) CYCLE
cC3-------GO OVER UPPER CONNECTIONS OF NODE N AND FILL Vx, Vy, Vz IN BOTH
c            DO II = IA(N)+1,IA(N+1)-1
c              JJ = JA(II)
c              IF(JJ.GT.N.AND.JJ.LE.NODES)THEN
c                IIS = JAS(II)
c                IF(IBOUND(JJ).NE.0)THEN   
c                   ANGLE = ARAD(IIS) 
c                   IF(ANGLE .GT. 1.57079. AND. ANGLE. LT. 1.57080)THEN
c                     ARAD(IIS) = -1.570796 ! CONVERT + 90 DEGREES TO - 90 DEGREES  
c                   ELSEIF(ANGLE .GT. 4.712. AND. ANGLE. LT. 4.713)THEN 
c                     ARAD(IIS) = 1.570796  ! CONVERT - 90 DEGREES TO + 90 DEGREES
c                   ENDIF  
c                ENDIF    
c              ENDIF     
c            ENDDO       
c          ENDDO         
c        ENDIF            
      ENDIF
C2------LOOP OVER ALL LAYERS TO DEFINE ARRAYS
      DO 200 K = 1,NLAY
      KK = K
      NNDLAY = NODLAY(K)
      NSTRT = NODLAY(K-1)+1
      NDSLAY = NNDLAY - NODLAY(K-1)
C
C2A-------Perform checks for unstructured grid formulations
      IF(IKCFLAG.NE.0)THEN
        IF(LAYCON(K).EQ.1.OR.LAYCON(K).EQ.3)WRITE(IOUT,10)K
      ENDIF
10    FORMAT(5X,'**LAYTYP=1 IS NOT ALLOWED WITH IKCFLAG = 1 OR -1,',
     *  1X,'SINCE K OF CONNECTIVITY IS READ. CHECK LAYER',I8)
C-------------------------------------------------
      IF(IKCFLAG.NE.0)GO TO 120
C-------------------------------------------------
C3------DEFINE ARRAYS FOR EACH LAYER
C3A-----DEFINE HORIZONTAL HYDRAULIC CONDUCTIVITY (HK)
      IF(NPHK.EQ.0) THEN
         CALL U1DREL(HK(NSTRT),ANAME(1),NDSLAY,K,IN,IOUT)
      ELSE
         READ(IN,*) LAYFLG(1,K)
         WRITE(IOUT,121) ANAME(1),K,LAYFLG(1,K)
  121    FORMAT(1X,/1X,A,' FOR LAYER',I4,
     1   ' WILL BE DEFINED BY PARAMETERS',/1X,'(PRINT FLAG=',I4,')')
         CALL UPARARRSUB1(HK(NSTRT),NDSLAY,1,KK,'HK',
     1      IOUT,ANAME(1),LAYFLG(1,KK))
         IF(NOPCHK.EQ.0) CALL UPARARRCK(BUFF,IBOUND(NSTRT),IOUT,K,
     1    NDSLAY,1,1,IUNSTR,'HK  ')
      END IF
C
C3B-----READ HORIZONTAL ANISOTROPY IF CHANI IS NON-ZERO
      IF(CHANI(K).LE.ZERO) THEN
        KHANI=-CHANI(K)
        IF(NPHANI.EQ.0) THEN
           CALL U1DREL(HANI(NSTRT),ANAME(2),NDSLAY,K,IN,IOUT)
        ELSE
           READ(IN,*) LAYFLG(6,K)
           WRITE(IOUT,121) ANAME(2),K,LAYFLG(6,K)
           CALL UPARARRSUB1(HANI(NSTRT),NDSLAY,1,KK,'HANI',
     &      IOUT,ANAME(2),LAYFLG(6,KK))
           IF(NOPCHK.EQ.0) CALL UPARARRCK(BUFF,IBOUND(NSTRT),IOUT,K,
     &      NDSLAY,1,1,IUNSTR,'HANI')
        END IF
      END IF
C
C3C-----DEFINE VERTICAL HYDRAULIC CONDUCTIVITY OR HORIZONTAL TO VERTICAL
C3C-----ANISOTROPY (VKA).
      IANAME=3
      PTYP='VK'
      IF(LAYVKA(K).NE.0) THEN
         IANAME=4
         PTYP='VANI'
      END IF
      IF(NPVK.EQ.0 .AND. NPVANI.EQ.0) THEN
         CALL U1DREL(VKA(NSTRT),ANAME(IANAME),NDSLAY,K,IN,IOUT)
      ELSE
         READ(IN,*) LAYFLG(2,K)
         WRITE(IOUT,121) ANAME(IANAME),K,LAYFLG(2,K)
         CALL UPARARRSUB1(VKA(NSTRT),NDSLAY,1,KK,PTYP,IOUT,
     &                       ANAME(IANAME),LAYFLG(2,KK))
         IF(NOPCHK.EQ.0) CALL UPARARRCK(BUFF,IBOUND(NSTRT),IOUT,K,
     &    NDSLAY,1,1,IUNSTR,PTYP)
      END IF
C-------------------------------------------------
120   CONTINUE
C-------------------------------------------------
C
C3D-----DEFINE SPECIFIC STORAGE OR STORAGE COEFFICIENT IN ARRAY SC1 IF TRANSIENT.
      IF(ITRSS.NE.0) THEN
         IF(NPSS.EQ.0) THEN
            CALL U1DREL(SC1(NSTRT),STOTXT,NDSLAY,K,IN,IOUT)
         ELSE
            READ(IN,*) LAYFLG(3,K)
            WRITE(IOUT,121) STOTXT,K,LAYFLG(3,K)
            CALL UPARARRSUB1(SC1(NSTRT),NDSLAY,1,KK,'SS',
     1           IOUT,STOTXT,LAYFLG(3,KK))
            IF(NOPCHK.EQ.0)CALL UPARARRCK(BUFF,IBOUND(NSTRT),IOUT,K,
     1       NDSLAY,1,1,IUNSTR,'SS  ')
         END IF
      END IF
C
C3E-----DEFINE SPECIFIC YIELD IN ARRAY SC2 IF TRANSIENT AND LAYER IS
C3E-----IS CONVERTIBLE.
      IF(LAYTYP(K).NE.0) THEN
         IF(ITRSS.NE.0) THEN
            IF(NPSY.EQ.0) THEN
               CALL U1DREL(SC2(NSTRT),ANAME(7),NDSLAY,K,IN,IOUT)
            ELSE
               READ(IN,*) LAYFLG(4,K)
               WRITE(IOUT,121) ANAME(7),K,LAYFLG(4,K)
               CALL UPARARRSUB1(SC2(NSTRT),NDSLAY,
     1         1,KK,'SY',IOUT,ANAME(7),LAYFLG(4,KK))
             IF(NOPCHK.EQ.0)CALL UPARARRCK(BUFF,IBOUND(NSTRT),IOUT,K,
     1        NDSLAY,1,1,IUNSTR,'SY  ')
            END IF
         END IF
      END IF
C
C3F-----READ CONFINING BED VERTICAL HYDRAULIC CONDUCTIVITY (VKCB) IF NODAL INPUT
      IF(IKCFLAG.EQ.0.AND.LAYCBD(K).NE.0) THEN
         IF(NPVKCB.EQ.0) THEN
            CALL U1DREL(VKCB(NSTRT),ANAME(5),NDSLAY,K,IN,
     1             IOUT)
         ELSE
            READ(IN,*) LAYFLG(5,K)
            WRITE(IOUT,121) ANAME(5),K,LAYFLG(5,K)
            CALL UPARARRSUB1(VKCB(NSTRT),NDSLAY,1,KK,
     1         'VKCB',IOUT,ANAME(5),LAYFLG(5,KK))
            IF(NOPCHK.EQ.0) CALL UPARARRCK(BUFF,IBOUND(NSTRT),IOUT,K,
     1        NDSLAY,1,1,IUNSTR,'VKCB')
         END IF
      END IF
C
C3G-----READ WETDRY CODES IF WETTING CAPABILITY HAS BEEN INVOKED
C3G-----(LAYWET NOT 0).
      IF(LAYWET(K).NE.0) THEN
         CALL U1DREL(WETDRY(NSTRT),ANAME(8),NDSLAY,K,IN,IOUT)
      END IF
c
C---------------------------------------------------------
      if(LAYCON(k).NE.5) go to 300
      IF(ITABRICH.EQ.0) THEN      
C3H-----READ alpha, beta, brook
      CALL U1DREL(alpha(NSTRT),ANAME(12),NDSLAY,K,IN,IOUT)
      CALL U1DREL(beta(NSTRT),ANAME(13),NDSLAY,K,IN,IOUT)
      CALL U1DREL(sr(NSTRT),ANAME(14),NDSLAY,K,IN,IOUT)
      CALL U1DREL(brook(NSTRT),ANAME(15),NDSLAY,K,IN,IOUT)
      IF(IBPN.GT.0)THEN
        CALL U1DREL(bP(NSTRT),ANAME(15),NDSLAY,K,IN,IOUT)  
      ENDIF
C6------read tabular input for retention and relative permeability curves      
      ELSE


      ENDIF      
  300 CONTINUE
C---------------------------------------------------------------
  200 CONTINUE
C---------------------------------------------------------------      
      IF(IKCFLAG.NE.0)THEN
C4--------READ EFFECTIVE SATURATED K OF CONNECTION
        ALLOCATE(TEMP(NJAS))
        CALL U1DRELNJA(TEMP(1),IATMP,ANAME(10),NJATMP,IN,IOUT,IDSYMRD)
        IF(IKCFLAG.EQ.1)THEN
          DO IIS=1,NJAS
            PGF(IIS) = PGF(IIS) * TEMP(IIS)
          ENDDO
C-----------INCLUDE THICKNESS TERM          
          DO N=1,NODES
            THICK1 = TOP(N) - BOT(N)
C-----------GO OVER CONNECTIONS OF NODE N AND FILL FOR UPPER SYMMETRIC PART
            DO II = IA(N)+1,IA(N+1)-1
              JJ = JA(II)
              IF(JJ.GE.N.AND.JJ.LE.NODES)THEN
                IIS = JAS(II)                  
C                IF(IVC(IIS).NE.0) CYCLE ! DO ONLY FOR HORIZONTAL CONNECTION
                THICK2 = TOP(JJ) - BOT(JJ)
                THICK = 0.5 * (THICK1 + THICK2)
                PGF(IIS) = PGF(IIS) * THICK
              ENDIF
            ENDDO
          ENDDO
        ELSE
          DO IIS=1,NJAS
            PGF(IIS) = TEMP(IIS)
          ENDDO
        ENDIF
C-------SET HK FOR THEIM SOLUTION CONNECTION
          DO N=1,NODES
            THICK = TOP(N) - BOT(N)
            AKN = 0.0
            IKN = 0
C-----------GO OVER CONNECTIONS OF NODE N AND FILL FOR UPPER SYMMETRIC PART
            DO II = IA(N)+1,IA(N+1)-1
              JJ = JA(II)
              IF(JJ.LE.NODES)THEN
                IIS = JAS(II)
                IF(IVC(IIS).EQ.0)THEN
                  IF(IKCFLAG.EQ.1) THEN
                   AKN = AKN + PGF(IIS) / THICK 
                  ELSE
                    AKN = AKN + PGF(IIS) / THICK * CL1(IIS)  
                  ENDIF
                ENDIF
                IKN = IKN + 1
              ENDIF
            ENDDO
            IF(IKN.GT.0) THEN
              HK(N) = AKN / IKN
            ENDIF
          ENDDO        
        DEALLOCATE(TEMP)
      ENDIF
C
C5------RETURN
      RETURN
      END
C -----------------------------------------------------------------------
      SUBROUTINE CBCK12(II,N1,N2,K,EL1,EL2,ANUM,HK,TH1,TH2)
C     ******************************************************************
C-------COMPUTE MEAN T OR K BETWEEN NODE 1 AND NODE 2
C     ******************************************************************
C
C      SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,      ONLY:NODES,ISYM,ARAD,JAS
      USE GWFBCFMODULE,ONLY:LAYAVG,CHANI,HANI,IHANISO
      DIMENSION HK(NODES)
      DOUBLE PRECISION ANUM,RATIO,FRAC1,FRAC2,T,HY1,HY2,TH1,TH2,HYMEAN,
     *  AVT
C     ------------------------------------------------------------------
C1------GET HY OF CONNECTION FOR ISOTRIPIC MEDIUM
C1------THIS HY IS THE X-DIRECTION COMPONENT FOR ANISOTROPIC MEDIUM
      ZERO=0.
      IF(LAYAVG(K).NE.2)THEN
        HY1 = HK(N1)*TH1
        HY2 = HK(N2)*TH2
      ELSE
        HY1 = HK(N1)
        HY2 = HK(N2)
      ENDIF
      CALL CBCKISO(HY1,HY2,K,EL1,EL2,HYMEAN)
      ANUM = HYMEAN
C2------FOR ANISOTROPIC SYSTEM ABOVE IS X-DIRECTION COMPONENT;
C2------GET Y-DIRECTION COMPONENT AND TAKE RESULTANT
      IF(IHANISO.EQ.1)THEN
        IF(CHANI(K).LE.ZERO) THEN
          KHANI=-CHANI(K)
          HY1=HY1*HANI(N1)
          HY2=HY2*HANI(N2)
        ELSE
          HY1 = HY1* CHANI(K)
          HY2 = HY2* CHANI(K)
        ENDIF
C2A-------Y-DIRECTION COMPONENT OF HY BETWEEN NODES N1 AND N2
        CALL CBCKISO(HY1,HY2,K,EL1,EL2,HYMEAN)
C2B-------HY OF CONNECTION IS RESULTANT OF X- AND Y-DIRECTION COMPONENTS
        IIS = JAS(II)
        ANGLE = ARAD(IIS)
        ACA = COS(ANGLE) * COS(ANGLE)         ! ABS(COS(ANGLE))
        ASA = SIN(ANGLE) * SIN(ANGLE)         ! ABS(SIN(ANGLE))
        IF(ANUM. LT. 1.0E-20. AND. HYMEAN. LT. 1.0E-20) THEN
          ANUM = 0.0
        ELSE  
          ANUM = 1.0 / (ACA/ANUM + ASA/HYMEAN)  ! ANUM * ACA + HYMEAN*ASA
        ENDIF  
      ENDIF
C3------FOR LAYAVG OF 2, MULTIPLY LOG MEAN K BY AVERAGE SAT THICKNESS
      IF(LAYAVG(K).EQ.2)THEN
        AVT = (TH1 + TH2) * 0.5
CSP        AVT = (TH1 * EL1 + TH2 * EL2) / (EL1 + EL2)    !WITH WEIGHTED AVERAGE
        ANUM = ANUM *AVT
      ENDIF
C
C4------RETURN
      RETURN
      END
C -----------------------------------------------------------------------
      SUBROUTINE CBCKISO(HY1,HY2,K,EL1,EL2,HYMEAN)
C     ******************************************************************
C-------COMPUTE MEAN OF HY1 AND HY2 IN HYMEAN DEPENDING ON LAYAVG
C     ******************************************************************
C
C      SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,      ONLY:NODES
      USE GWFBCFMODULE,ONLY:LAYAVG,CHANI,HANI
      DOUBLE PRECISION HYMEAN,RATIO,FRAC1,FRAC2,T,HY1,HY2,DDENOM
C     ------------------------------------------------------------------
C1------COMPUTE WEIGHTED HARMONIC AVERAGE
      IF(LAYAVG(K).EQ.0)THEN
        DDENOM = HY1*EL2 + HY2*EL1
        IF(DDENOM.LT.1.0E-20) DDENOM=1.0E-20
        HYMEAN = HY1*HY2*(EL1+EL2)/DDENOM
C2------COMPUTE ARITHMETIC AVERAGE
      ELSEIF(LAYAVG(K).EQ.3)THEN
        HYMEAN = (HY1 + HY2)*0.5
      ELSEIF(LAYAVG(K).EQ.4)THEN
C3-------FINITE ELEMENT CONNECTIVITY IS WEIGHTED ARITHMETIC AVERAGE
        HYMEAN = (HY1*EL1 + HY2*EL2)/(EL1+EL2)
C4------COMPUTE LOGARITHMIC AVERAGE
      ELSEIF(LAYAVG(K).EQ.1.OR.LAYAVG(K).EQ.2)THEN
        FRAC1=1.005
        FRAC2=0.995
        RATIO=HY2/HY1
        IF(RATIO.GT.FRAC1.OR.RATIO.LT.FRAC2) THEN
           T=(HY2-HY1)/LOG(RATIO)
        ELSE
           T=(HY1 + HY2)*0.5
        END IF
        HYMEAN = T
      ENDIF
C4------RETURN
      RETURN
      END
C-------------------------------------------------------------------
      SUBROUTINE SGWF2LPFU1N
C     ******************************************************************
C     INITIALIZE AND CHECK LPF DATA
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,      ONLY:NCOL,NROW,NLAY,IBOUND,HNEW,LAYCBD,IUNSTR,
     1                      DELR,DELC,IOUT,NODLAY,IA,JA,JAS,IVC,ISYM
      USE GWFBCFMODULE,ONLY:IBCFCB,IWDFLG,IWETIT,IHDWET,WETFCT,HDRY,CV,
     1                      LAYCON,LAYAVG,SC1,SC2,WETDRY,
     2                      IKCFLAG,laywet,HK,VKCB,LAYTYP,VKA
C     ------------------------------------------------------------------
C
C1------DEFINE CONSTANTS.
      ZERO=0.
      HCNV=888.88
C
C2-------INSURE THAT EACH ACTIVE CELL HAS AT LEAST ONE NON-ZERO
C2-------TRANSMISSIVE PARAMETER.
      DO 60 K=1,NLAY
        NNDLAY = NODLAY(K)
        NSTRT = NODLAY(K-1)+1
      IF(LAYWET(K).NE.0) THEN
C
C3------WETTING IS ACTIVE.
        DO 40 N=NSTRT,NNDLAY
         IF(IBOUND(N).EQ.0 .AND. WETDRY(N).EQ.ZERO)
     1                GO TO 40
C
C3A-----CHECK HORIZONTAL HYDRAULIC CONDUCTIVITY (HK).
         IF(HK(N).NE.ZERO) GO TO 40
C
         JPREV=0
         JNEXT=0
         DO II=IA(N)+1,IA(N+1)-1
          JJ=JA(II)
          IIS = JAS(II)
          IF(JJ.LT.N.AND.IVC(IIS).EQ.1)THEN
            JPREV=JJ
          ELSEIF(IVC(IIS).EQ.1)THEN
            JNEXT=JJ
          ENDIF
        ENDDO
C
C3B-----CHECK VERTICAL HYDRAULIC CONDUCTIVITY AND CONFINING BED
C3B-----VERTICAL HYDRAULIC CONDUCTIVITY.
         IF(NLAY.GT.1) THEN
            IF(VKA(N).NE.ZERO) THEN
               IF(K.NE.NLAY) THEN
                  IF (VKA(JNEXT).NE.ZERO) THEN
                     IF(LAYCBD(K).NE.0) THEN
                        IF(VKCB(N).NE.ZERO) GO TO 40
                     ELSE
                        GO TO 40
                     END IF
                  END IF
               END IF
               IF(K.NE.1) THEN
                  IF (VKA(JPREV).NE.ZERO) THEN
                     IF (LAYCBD(K-1).NE.0) THEN
                        IF(VKCB(JPREV).NE.ZERO) GO TO 40
                     ELSE
                        GO TO 40
                     END IF
                  ENDIF
               END IF
            END IF
         END IF
C
C3C-----ALL TRANSMISSIVE TERMS ARE ALL 0, SO CONVERT CELL TO NO FLOW.
         IBOUND(N)=0
         HNEW(N)=HCNV
         WETDRY(N)=ZERO
         IF(IUNSTR.EQ.0)THEN
           KK = (N-1) / (NCOL*NROW) + 1
           IJ = N - (KK-1)*NCOL*NROW
           I = (IJ-1)/NCOL + 1
           J = IJ - (I-1)*NCOL
           WRITE(IOUT,43) KK,I,J
         ELSE
           WRITE(IOUT,43) N
         ENDIF
   40    CONTINUE
C
      ELSE
C
C4------WETTING IS INACTIVE
         DO 50 N=NSTRT,NNDLAY
         IF(IBOUND(N).EQ.0) GO TO 50
C
C4A-----CHECK HORIZONTAL HYDRAULIC CONDUCTIVITY (HK).
         IF(HK(N).NE.ZERO) GO TO 50
C
         JPREV=0
         JNEXT=0
         DO II=IA(N)+1,IA(N+1)-1
          JJ=JA(II)
          IIS = JAS(II)
          IF(JJ.LT.N.AND.IVC(IIS).EQ.1)THEN
            JPREV=JJ
          ELSEIF(IVC(IIS).EQ.1)THEN
            JNEXT=JJ
          ENDIF
        ENDDO
C
C4B-----CHECK VERTICAL HYDRAULIC CONDUCTIVITY AND CONFINING BED
C4B-----VERTICAL HYDRAULIC CONDUCTIVITY.
         IF(NLAY.GT.1) THEN
            IF(VKA(N).NE.ZERO) THEN
               IF(K.NE.NLAY) THEN
                  IF (VKA(JNEXT).NE.ZERO) THEN
                     IF(LAYCBD(K).NE.0) THEN
                        IF(VKCB(N).NE.ZERO) GO TO 50
                     ELSE
                        GO TO 50
                     END IF
                  END IF
               END IF
               IF(K.NE.1) THEN
                  IF (VKA(JPREV).NE.ZERO) THEN
                     IF (LAYCBD(K-1).NE.0) THEN
                        IF(VKCB(JPREV).NE.ZERO) GO TO 50
                     ELSE
                        GO TO 50
                     END IF
                  ENDIF
               END IF
            END IF
         END IF
C
C4C-----ALL TRANSMISSIVE TERMS ARE 0, SO CONVERT CELL TO NO FLOW.
         IBOUND(N)=0
         HNEW(N)=HCNV
         WRITE(IOUT,43) N
   43    FORMAT(1X,'NODE (LAYER,ROW,COL) ',I8,
     1 ' ELIMINATED BECAUSE ALL HYDRAULIC',/,
     2 ' CONDUCTIVITIES TO NODE ARE 0')
   50    CONTINUE
      END IF
   60 CONTINUE
C
C5------RETURN.
      RETURN
      END
      SUBROUTINE SGWF2LPFU1SC(SC,ISPST)
C     ******************************************************************
C     COMPUTE STORAGE CAPACITY
C     ******************************************************************
C
C     SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,        ONLY:NCOL,NROW,DELR,DELC,LAYCBD,NODES,TOP,BOT,
     1                   AREA
C
      DIMENSION SC(NODES)
C     ------------------------------------------------------------------
C
C1------MULTIPLY SPECIFIC STORAGE BY THICKNESS, DELR, AND DELC TO GET
C1------CONFINED STORAGE CAPACITY.
      IF(ISPST.NE.0) THEN
         DO 80 N=1,NODES
         THICK=TOP(N)-BOT(N)
         SC(N)=SC(N)*THICK*AREA(N)
   80    CONTINUE
      ELSE
C
C2------MULTIPLY SPECIFIC YIELD BY DELR AND DELC TO GET UNCONFINED
C2------STORAGEE CAPACITY(SC2).
         DO 85 N=1,NODES
         SC(N)=SC(N)*AREA(N)
   85    CONTINUE
      END IF
C3------return
      RETURN
      END
      SUBROUTINE SGWF2LPFU1VCOND
C     ******************************************************************
C     COMPUTE VERTICAL BRANCH CONDUCTANCE BETWEEN A LAYER AND THE NEXT
C     LOWER LAYER FROM VERTICAL HYDRAULIC CONDUCTIVITY IN PGF.
C     ******************************************************************
C
C      SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,        ONLY:NCOL,NROW,NLAY,IBOUND,HNEW,DELR,DELC,FAHL,
     1                        TOP,BOT,LAYCBD,IOUT,STRT,AREA,CL1,CL2,PGF,
     1                        IUNSTR,NODLAY,IA,JA,JAS,IVC,ISYM
      USE GWFBCFMODULE,ONLY:IBCFCB,IWDFLG,IWETIT,IHDWET,WETFCT,HDRY,CV,
     1                      LAYCON,LAYAVG,LAYAVGV,
     2                      LAYTYP,CHANI,LAYVKA,HK,VKA,VKCB,LAYSTRT
C
      DOUBLE PRECISION BBOT,TTOP,CONAREA
C     ------------------------------------------------------------------
C
      IF(NLAY.EQ.1) RETURN
      ZERO=0.
      HALF=0.5
C
C1------LOOP THROUGH ALL CELLS IN THE LAYER.
      DO 100 K=1,NLAY
      NNDLAY = NODLAY(K)
      NSTRT = NODLAY(K-1)+1
      DO 100 N=NSTRT,NNDLAY
C2----------GO OVER CONNECTIONS OF NODE N
          DO II = IA(N)+1,IA(N+1)-1
            JJ = JA(II)
            IF(JJ.LE.N) CYCLE
            IIS = JAS(II)
            IF(IUNSTR.EQ.1)THEN
C              CONAREA = PGF(IIS) * (CL1(IIS)+CL2((IIS))
              CONAREA = FAHL(IIS)
            ELSE
              CONAREA = AREA(N)
            ENDIF
C3------------COMPUTE WHEN VERTICAL DIRECTION IS FOUND
            IF(IVC(IIS).EQ.1)THEN !VERTICAL DIRECTION CONNECTION
                KK = K+1
C
C4--------------CALCULATE VERTICAL HYDRAULIC CONDUCTIVITY FOR CELL.
                hyc1 = 0.0
                IF(LAYVKA(K).EQ.0) THEN
                  HYC1=VKA(N)
                ELSE
                   if(vka(n).gt.1.0e-20) HYC1=HK(N)/VKA(N)
                END IF
                IF(HYC1.GT.ZERO) THEN
C5----------------CALCULATE VERTICAL HYDRAULIC CONDUCTIVITY FOR ADJACENT CELL.
                  hyc2 = 0.0  
                  IF(LAYVKA(KK).EQ.0) THEN
                    HYC2=VKA(JJ)
                  ELSE
                    if(vka(jj).gt.1.0e-20) HYC2=(HK(JJ)/VKA(JJ))
                  END IF
                  IF(HYC2.GT.ZERO) THEN
C
C6------------------CALCULATE INVERSE LEAKANCE FOR CELL.
                    BBOT=BOT(N)
                    TTOP=TOP(N)
                    IF(LAYSTRT(K).NE.0) TTOP=STRT(N)
                    BOVK1=(TTOP-BBOT)*HALF/HYC1
                    IF(BOVK1.LT.1.0E-20) BOVK1 = 1.0E-20
C
C7------------------CALCULATE INVERSE LEAKANCE FOR ADJACENT CELL.
                    BBOT=BOT(JJ)
                    TTOP=TOP(JJ)
                    IF(LAYSTRT(KK).NE.0) TTOP=STRT(JJ)
                    BOVK2=(TTOP-BBOT)*HALF/HYC2
                    IF(BOVK2.LT.1.0E-20) BOVK2 = 1.0E-20
C
C8------------------CALCULATE VERTICAL HYDRAULIC CONDUCTIVITY FOR CONFINING BED.
                    IUP = N
                    IDN = JJ
                    IK = K
C
                    IF(LAYCBD(IK).NE.0) THEN
                      IF(VKCB(IUP).GT.ZERO) THEN
C
C9----------------------CALCULATE INVERSE LEAKANCE FOR CONFINING BED.
                        B=BOT(IUP)-TOP(IDN)
                        IF(B.LT.ZERO) THEN
                          IF(IUNSTR.EQ.0)THEN
                            KK = (N-1) / (NCOL*NROW) + 1
                            IJ = N - (KK-1)*NCOL*NROW
                            I = (IJ-1)/NCOL + 1
                            J = IJ - (I-1)*NCOL
                            WRITE(IOUT,45) KK,I,J
                          ELSE
                            WRITE(IOUT,47) N
                          ENDIF
   45                     FORMAT(1X,/1X,
     1  'Negative confining bed thickness below cell (Layer,row,col)',
     2                  I4,',',I5,',',I5)
   47                     FORMAT(1X,/1X,
     1  'Negative confining bed thickness below cell (Node)',I9)
                          WRITE(IOUT,46) BOT(IUP),TOP(IDN)
   46            FORMAT(1X,'Top elevation, bottom elevation:',1P,2G13.5)
                          CALL USTOP(' ')
                        END IF
                        CBBOVK=B/VKCB(IUP)
                        IF(CBBOVK.LT.1.0E-20) CBBOVK = 1.0E-20
C10------------------------COMPUTE PGF TERM
C10A------------------------COMPUTE PGF TERM WITH CONFINING BED
                        IF(LAYAVGV(K).EQ.0)THEN
C10A1 ----------------------USE HARMONIC MEAN IF LAYAVGV = 0
                          PGF(IIS)=CONAREA/(BOVK1+CBBOVK+BOVK2)
                        ELSE
C10A1 ----------------------USE ARITHMETIC MEAN IF LAYAVGV = 1
                          PGF(IIS)=CONAREA*
     1                       (1.0/BOVK1+1.0/CBBOVK+1.0/BOVK2) / 3.0    
                        ENDIF
                      ELSE
                        PGF(IIS) = 0.0
                      END IF
C10B------------------------COMPUTE PGF TERM WITH CONFINING BED
                    ELSE
                      IF(LAYAVGV(K).EQ.0)THEN
C10A1 --------------------USE HARMONIC MEAN IF LAYAVGV = 0                        
                        PGF(IIS)=CONAREA/(BOVK1+BOVK2)
                      ELSE
C10A1 ------------------USE ARITHMETIC MEAN IF LAYAVGV = 1
                        PGF(IIS)=CONAREA*(1.0/BOVK1+1.0/BOVK2) / 2.0
                      ENDIF                      
                    END IF
                  END IF
                END IF
            END IF
          ENDDO
  100 CONTINUE
C
C11------RETURN.
      RETURN
      END
C------------------------------------------------------------------------
      SUBROUTINE SGWF2LPFU1VCONDV(K)
C     ******************************************************************
C     COMPUTE VERTICAL BRANCH CONDUCTANCE BETWEEN A LAYER AND THE NEXT
C     LOWER LAYER FROM VERTICAL HYDRAULIC CONDUCTIVITY IN AMAT.
C     ******************************************************************
C
C      SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,        ONLY:NCOL,NROW,NLAY,IBOUND,HNEW,DELR,DELC,FAHL,
     1                   TOP,BOT,LAYCBD,IOUT,STRT,AREA,PGF,
     1                   AMAT,IUNSTR,NODLAY,IA,JA,JAS,IVC,ISYM,NOCVCO
      USE GWFBCFMODULE,ONLY:IBCFCB,IWDFLG,IWETIT,IHDWET,WETFCT,HDRY,CV,
     1                      LAYCON,LAYAVG,LAYAVGV,
     2                      LAYVKA,HK,VKA,VKCB,LAYSTRT
C
      DOUBLE PRECISION BBOT,TTOP,CONAREA,HHD,BEE
C     ------------------------------------------------------------------
C
      ZERO=0.
      HALF=0.5
C
C1------LOOP THROUGH ALL CELLS IN THE LAYER.
      NNDLAY = NODLAY(K)
      NSTRT = NODLAY(K-1)+1
      DO 100 N=NSTRT,NNDLAY
        IF(IBOUND(N).NE.0)THEN
C2----------GO OVER CONNECTIONS OF NODE N
          DO II = IA(N)+1,IA(N+1)-1
            JJ = JA(II)
            IF(JJ.LE.N) CYCLE
            IIS = JAS(II)
            IF(IUNSTR.EQ.1)THEN
              CONAREA = FAHL(IIS)
            ELSE
              CONAREA = AREA(N)
            ENDIF
C3------------COMUTE WHEN VERTICAL DIRECTION IS FOUND
            IF(IVC(IIS).EQ.1)THEN
              ISLOC = ISYM(II)
              IF(IBOUND(JJ).NE.0) THEN
                KK = K+1
                IF(NLAY.EQ.1) KK = K                
C
C4--------------CALCULATE VERTICAL HYDRAULIC CONDUCTIVITY FOR CELL.
                IF(LAYVKA(K).EQ.0) THEN
                  HYC1=VKA(N)
                ELSE
                   HYC1=HK(N)/VKA(N)
                END IF
                IF(HYC1.GT.ZERO) THEN
C5----------------CALCULATE VERTICAL HYDRAULIC CONDUCTIVITY FOR ADJACENT CELL.
                  IF(LAYVKA(KK).EQ.0) THEN
                    HYC2=VKA(JJ)
                  ELSE
                    HYC2=(HK(JJ)/VKA(JJ))
                  END IF
                  IF(HYC2.GT.ZERO) THEN
C
C6------------------CALCULATE INVERSE LEAKANCE FOR CELL.
                    BBOT=BOT(N)
                    TTOP=TOP(N)
                    IF(LAYSTRT(K).NE.0) TTOP=STRT(N)
                    IF(LAYCON(K).NE.0.AND.LAYCON(K).NE.2)THEN
                      HHD=HNEW(N)
                      IF(HHD.LT.TTOP) TTOP=HHD
                    ENDIF
                    BOVK1=(TTOP-BBOT)*HALF/HYC1
C
C7------------------CALCULATE INVERSE LEAKANCE FOR ADJACENT CELL.
                    BBOT=BOT(JJ)
                    TTOP=TOP(JJ)
                    IF(LAYSTRT(KK).NE.0) TTOP=STRT(JJ)
                    BEE = (TTOP-BBOT)*HALF
                    IF(NOCVCO.EQ.0.AND.
     *               LAYCON(K+1).NE.0.AND.LAYCON(K+1).NE.2)THEN
                      HHD=HNEW(JJ)
                      IF(HHD.LT.TTOP) BEE = 0.0
                    ENDIF
                    BOVK2= BEE/HYC2
C
C8------------------CALCULATE VERTICAL HYDRAULIC CONDUCTIVITY FOR CONFINING BED.
                    IUP = N
                    IDN = JJ
                    IK = K
C
                    IF(LAYCBD(IK).NE.0) THEN
                      IF(VKCB(IUP).GT.ZERO) THEN
C
C9A----------------------CALCULATE INVERSE LEAKANCE FOR CONFINING BED.
                        B=BOT(IUP)-TOP(IDN)
                        CBBOVK=B/VKCB(IUP)
C10-----------------------COMPUTE AMAT TERM
C10A----------------------COMPUTE AMAT TERM WITH CONFINING BED
                        IF(LAYAVGV(K).EQ.0)THEN
C10A1---------------------HARMONIC AVERAGING FOR LAYAVGV=0                            
                          AMAT(II)=CONAREA/(BOVK1+CBBOVK+BOVK2)
                        ELSE
C10A2---------------------ARITHMETIC AVERAGING FOR LAYAVGV=1                            
                          AMAT(II)=CONAREA*(1./BOVK1+1./CBBOVK+1./BOVK2)
     1                       /3.0 
                        ENDIF
                        AMAT(ISLOC) = AMAT(II)
                      END IF
C10B--------------------COMPUTE AMAT TERM WITH CONFINING BED
                    ELSE
                      IF(LAYAVGV(K).EQ.0)THEN
C10B1-------------------HARMONIC AVERAGING FOR LAYAVGV=0                           
                        AMAT(II)= CONAREA/(BOVK1+BOVK2)
                      ELSE
C10B2---------------------ARITHMETIC AVERAGING FOR LAYAVGV=1                            
                        AMAT(II)=CONAREA*(1./BOVK1+1./BOVK2) /2.0   
                      ENDIF
                      AMAT(ISLOC)= AMAT(II)
                    END IF
                  END IF
                END IF
              END IF
            END IF
          ENDDO
        ENDIF
  100 CONTINUE
C
C11------RETURN.
      RETURN
      END
C------------------------------------------------------------------------
      SUBROUTINE SGWF2LPFU1CK(IOUT,NP,PTYP)
C     ******************************************************************
C     CHECK THAT JUST-DEFINED PARAMETER OF TYPE 'VK' OR 'VANI' IS USED
C     CONSISTENTLY WITH LAYVKA ENTRIES FOR LAYERS LISTED IN CLUSTERS FOR
C     THE PARAMETER
C     ******************************************************************
C
C      SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GWFBCFMODULE,  ONLY:LAYTYP,CHANI,LAYVKA
      USE PARAMMODULE
C
      CHARACTER*4 PTYP
C     ------------------------------------------------------------------
C
C1------LOOP THROUGH THE CLUSTERS FOR THIS PARAMETER.
      DO 10 ICL = IPLOC(1,NP),IPLOC(2,NP)
        LAY = IPCLST(1,ICL)
        LV = LAYVKA(LAY)
        IF (PTYP.EQ.'VK  ' .AND. LV.NE.0) THEN
          WRITE (IOUT,590) LAY,LV,LAY,PARNAM(NP),'VK'
  590     FORMAT(/,
     &1X,'LAYVKA entered for layer ',i3,' is: ',i3,'; however,',
     &' layer ',i3,' is',/,' listed in a cluster for parameter "',a,
     &'" of type ',a,' and')
          WRITE (IOUT,600)
  600     FORMAT(
     &1X,'parameters of type VK can apply only to layers for which',
     &/,' LAYVKA is specified as zero -- STOP EXECUTION (SGWF2LPFU1CK)')
          CALL USTOP(' ')
        ELSEIF (PTYP.EQ.'VANI' .AND. LV.EQ.0) THEN
          WRITE (IOUT,590) LAY,LV,LAY,PARNAM(NP),'VANI'
          WRITE (IOUT,610)
  610     FORMAT(
     &1X,'parameters of type VANI can apply only to layers for which',/,
     &' LAYVKA is not specified as zero -- STOP EXECUTION',
     &' (SGWF2LPFU1CK)')
          CALL USTOP(' ')
        ENDIF
   10 CONTINUE
C
C2------Return.
      RETURN
      END
C
C-----------------------------------------------------------------------
      SUBROUTINE GWF2BCFU1DA(INLPF)
      USE GWFBCFMODULE
      USE GLOBAL, ONLY:NLAY,IUNSTR,IDEALLOC_HY,iunsat
      INTEGER ALLOC_ERR
C
      IF(IDEALLOC_HY.EQ.0)  DEALLOCATE(HK)      
      DEALLOCATE(IBCFCB)
      DEALLOCATE(IWDFLG)
      DEALLOCATE(IWETIT)
      DEALLOCATE(IHDWET)
      DEALLOCATE(WETFCT)
      DEALLOCATE(HDRY)
      DEALLOCATE(LAYCON)
      DEALLOCATE(LAYAVG)
      DEALLOCATE(LAYAVGV)
      DEALLOCATE(SC1)
      DEALLOCATE(SC2)
      DEALLOCATE(WETDRY)
      DEALLOCATE(laywet)
      DEALLOCATE(ISFAC)
      DEALLOCATE(ITHFLG, STAT = ALLOC_ERR)
      DEALLOCATE(IHANISO)
      DEALLOCATE(CHANI)
      DEALLOCATE(ALPHA, BETA, SR, BROOK, STAT = ALLOC_ERR)
      IF(INLPF.NE.0)THEN
        DEALLOCATE(LAYTYP)          
        DEALLOCATE(LAYVKA)
        DEALLOCATE(LAYSTRT)
        DEALLOCATE(LAYFLG)
        DEALLOCATE(VKA)
        DEALLOCATE(VKCB)
        DEALLOCATE(HANI)
      ENDIF
      DEALLOCATE (CV) 
C
      RETURN
      END

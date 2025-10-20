C-----------------------------------------------------------------------
      SUBROUTINE GWT2DPTU1AR(IN,NPER)
C     ******************************************************************
C     INITIALIZE VARIABLES AND READ DATA FOR DUAL POROSITY TRANSPORT
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL, ONLY:IFREFM,IOUT,NODES,NLAY,IUNSAT,NODLAY,IDPT,IDPF,
     1    IBOUND,TOP,BOT,AREA,ISSFLG,PGF,IA,JAS,JA,Sn,ISSFLG,IDPIN
      USE GWTDPTMODULE
      USE GWTBCTMODULE, ONLY: MCOMP,ICBUND,ZODRW,ZODRS,FODRW,FODRS,
     1  PRSITY,ADSORB,IZOD,IFOD,IADSORB,CINACT,MCOMPT,NTCOMP,IHEAT,
     1  ICHAIN,ISPRCT,NPARENT,RHOW,HTCAPW,ZODRAW,FODRAW
      USE GWFBCFMODULE, ONLY: SC1,SC2,HK,LAYCON
      USE GWFDPFMODULE, ONLY: PHIF,SnIM,SoIM,SC1IM,SC2IM,HNEWIM
      USE AW_ADSORBIMMODULE, ONLY: IAREA_FNIM,IKAWI_FNIM,ITAB_AWIIM
      CHARACTER*400 LINE
      DOUBLE PRECISION MASLOCW,MASLOCS,MASLOCC,REFHD
      INTEGER IFRAHK,IMSAT
C
      REAL, DIMENSION(:),ALLOCATABLE  ::TEMPC
      REAL*8, DIMENSION(:),ALLOCATABLE  ::TEMPC8
      CHARACTER*24 ANAME(18),CNAME
      DATA ANAME(1) /'      DPT BOUNDARY ARRAY'/
      DATA ANAME(2) /'       FRACTURE POROSITY'/
      DATA ANAME(3) /'IMMOBILE DOMAIN POROSITY'/
      DATA ANAME(4) /'         DD BULK DENSITY'/
      DATA ANAME(5) /'    DD DISPERSION COEFF.'/
      DATA ANAME(6) /'   DD MASS TRANSFER RATE'/
      DATA ANAME(7) /'      DD ADSORPTION COEF'/
      DATA ANAME(8) /'  DD FREUNDLICH EXPONENT'/
      DATA ANAME(9) /'DD ZERO-ORDER DECAY WAT.'/
      DATA ANAME(10) /'DD ZERO-ORDER DECAY SOIL'/
      DATA ANAME(11) /'DD 1-ST ORDER DECAY WAT.'/
      DATA ANAME(12) /'DD 1-ST ORDER DECAY SOIL'/
      DATA ANAME(13) /'        DPT INITIAL CONC'/
      DATA ANAME(14) /'     IMMOBILE DOMAIN SAT'/
      DATA ANAME(15) /'      SOIL HEAT CAPACITY'/
      DATA ANAME(16) /'             TEMPERATURE'/      
      DATA ANAME(17) /'    ZERO-ORDER DECAY A-W'/       
      DATA ANAME(18) /'   FIRST-ORDER DECAY A-W'/              
C     ------------------------------------------------------------------
C
C1------IDENTIFY PACKAGE
      IDPT = 1 ! FLAG FOR DUAL POROSITY TRANSPORT IS ON
      WRITE(IOUT,1) IN
    1 FORMAT(1X,/1X,'DPT -- DUAL POROSITY TRANSPORT PACKAGE, VERSION 1',
     1',12/26/2012',/,9X,'INPUT READ FROM UNIT',I3)
C
C-----------------------------------------------------------------------
C2------ALLOCATE VARIABLES
      ALLOCATE(IDPTCB,IDPTCON,IC_IBNDIN_FLG,IADSORBIM,
     1  IDISPIM,IZODIM,IFODIM,IAW_ADSORBIM)
      IAW_ADSORBIM = 0
C
C-----------------------------------------------------------------------
C3------READ GENERAL DPT INFORMATION AND FLAGS
      IF(IFREFM.EQ.0) THEN
        CALL URDCOM(IN,IOUT,LINE)
        READ(LINE,2)IDPTCB,IDPTCON,IC_BNDIM_FLG,IADSORBIM,IDISPIM,
     1    IZODIM, IFODIM
        LLOC=71
      ELSE
        CALL URDCOM(IN,IOUT,LINE)
        LLOC=1
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IDPTCB,R,IOUT,IN)
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IDPTCON,R,IOUT,IN)
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IC_BNDIM_FLG,R,IOUT,IN)
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IADSORBIM,R,IOUT,IN)
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IDISPIM,R,IOUT,IN)
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IZODIM,R,IOUT,IN)
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IFODIM,R,IOUT,IN)
      ENDIF
2     FORMAT(7I10)
C      
C3B-----PRINT VALUES
      WRITE(IOUT,11) IDPTCB,IDPTCON,IC_BNDIM_FLG,IADSORBIM,
     *  IDISPIM,IZODIM,IFODIM
   11 FORMAT(1X,'C-B-C FLUX FLAG OR UNIT NUMBER (IDPTCB)',27X,' =',I3
     *  /1X,'IMMOBILE DOMAIN CONC SAVE FLAG OR UNIT NUMBER (IDPTCON)',
     *  11X,' =',I3,
     *  /1X,'ACTIVE IMMOBILE DOMAIN SAME AS FRACTURE DOMAIN FLAG',1X,
     *  '(IC_BNDIM_FLG) =',I3
     *  /1X,'IMMOBILE DOMAIN ADSORPTION FLAG (IADSORBIM)',23X,' =',I3
     *  /1X,'IMMOBILE DOMAIN DISPERSION FLAG (IDISPIM)',25X,' =',I3,
     */1X,'IMMOBILE DOMAIN ZERO-ORDER DECAY FLAG (IZODIM)',20X,' =',I3,
     */1X,'IMMOBILE DOMAIN FIRST-ORDER DECAY FLAG (IFODIM)',19X,' =',I3)
C
      IF(IDPTCB.GT.0) WRITE(IOUT,9) IDPTCB
    9 FORMAT(1X,'DUAL POROSITY CELL-BY-CELL MASS FLUX WILL BE SAVED ON',
     1 1X,'UNIT IDPTCB  =',I3)
      IF(IDPTCON.GT.0) WRITE(IOUT,10) IDPTCON
10    FORMAT(1X,'IMMOBILE DOMAIN CONC WILL BE SAVED ON UNIT IDPTCON',
     1  14X,'=',I3)
C-------------------------------------------------------------------------------
C3C--------GET OPTIONS FOR WHEN FLOW IS NOT DUAL POROSITY
      IF(IDPF.EQ.0)THEN
        IFRAHK = 0
        IMSAT = 1
100     CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,I,R,IOUT,IN)
        IF(LINE(ISTART:ISTOP).EQ.'FRAHK') THEN
C3D-------SET FLAG FOR HK TO BE FRACTURE-DOMAIN VALUE.   
          IFRAHK = 1
        ENDIF 
C3E-------SET FLAG FOR IMMOBILE DOMAIN SATURATION TREATMENT.          
        IF(LINE(ISTART:ISTOP).EQ.'MOBILESAT') THEN
C3F-------SET FLAG.   
          IMSAT = 2
        ELSEIF(LINE(ISTART:ISTOP).EQ.'INPUTSAT') THEN
C3G-------SET FLAG.   
          IMSAT = 3
        ENDIF      
        IF(LLOC.LT.200) GO TO 100 
C3H----------PRINT OPTION FLAG VALUES
        IF(IFRAHK.EQ.1) WRITE(IOUT,12) IFRAHK
12      FORMAT(1X,'CONDUCTANCE ARE FOR FRACTURE VOLUME; (IFRAHK)',12X,
     1   '=',I4)
        IF(IFRAHK.EQ.0) WRITE(IOUT,13) IFRAHK
13      FORMAT(1X,'CONDUCTANCE ARE FOR TOTAL (FRACTURE + MATRIX)', 
     1            1X,'DOMAIN; (IFRAHK)  =',I4)
        IF(IMSAT.EQ.1) WRITE(IOUT,14) IMSAT
14      FORMAT(1X,'IMMOBILE DOMAIN SATURATIONS SET TO UNITY)', 
     1            12X, '(IMSAT)  =',I4)
        IF(IMSAT.EQ.2) WRITE(IOUT,15) IMSAT
15      FORMAT(1X,'IMMOBILE DOMAIN SATURATIONS SET SAME AS FOR MOBILE', 
     1            1X, 'DOMAIN (IMSAT)  =',I4)        
        IF(IMSAT.EQ.3) WRITE(IOUT,16) IMSAT
16      FORMAT(1X,'IMMOBILE DOMAIN SATURATIONS ARE INPUT)', 
     1            12X, '(IMSAT)  =',I4)        
        ENDIF
C3I ----READ OTHER OPTIONS (A-W ADSORPTION)         
      LLOC = 1
105   CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,I,R,IOUT,IN)
      IF(LINE(ISTART:ISTOP).EQ.'A-W_ADSORBIM') THEN
C3D-------SET FLAG FOR HK TO BE FRACTURE-DOMAIN VALUE.   
        IAW_ADSORBIM = 1
        ALLOCATE(IAREA_FNIM,IKAWI_FNIM,ITAB_AWIIM)
        ITAB_AWIIM = 0
        WRITE(IOUT,33)
33      FORMAT(1X,'AIR-WATER INTERFACE ADSORPTION OPTION IN IMMOBILE',
     1    1X,'DOMAIN')
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IAREA_FNIM,R,IOUT,IN)
        WRITE(IOUT,46) IAREA_FNIM
46      FORMAT(1X,'AREA FUNCTION INDEX (IAREA_FNIM) = ', I10)
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IKAWI_FNIM,R,IOUT,IN)
        WRITE(IOUT,47) IKAWI_FNIM
47      FORMAT(1X,'K_WAI FUNCTION INDEX (IKAWI_FNIM) = ', I10)
      ELSEIF(LINE(ISTART:ISTOP).EQ.'SOMETHINGELSE') THEN
            
      ENDIF         
      IF(LLOC.LT.200) GO TO 105
C ---------------------------------------------------------------------
C-------ALLOCATE SPACE FOR AIR-WATER ADSORPTION ARRAYS
        IF(IAW_ADSORBIM.EQ.1) THEN
          CALL AW_ADSORBIM1AL (IN)
        ENDIF
C-----------------------------------------------------------------------
C4-----ALLOCATE ARRAYS AND INITIALIZE
      ALLOCATE(DDTTR(NODES),CONCIM(NODES,NTCOMP),CONCOIM(NODES,NTCOMP),
     1 ICBUNDIM(NODES),PRSITYIM(NODES),BULKDIM(NODES),MASSBCTIM(NTCOMP))
      ALLOCATE(DIADDT(NODES),RDDT(NODES),OFFDDTM(NODES),OFFDDTIM(NODES))
C      
      IF(ICHAIN.NE.0.AND.ISPRCT.NE.0) THEN 
        ALLOCATE (SPTLRCTIM(NODES,NTCOMP,NTCOMP))  
      ENDIF    
C      
      IF(IADSORBIM.NE.0) ALLOCATE (ADSORBIM(NODES,MCOMPT))
      IF(IADSORBIM.EQ.2.OR.IADSORBIM.EQ.3) 
     1  ALLOCATE (FLICHIM(NODES,MCOMP))
      IF(IDISPIM.GT.0)THEN
        ALLOCATE(DLIM(NODES))
      ENDIF
      IF(IZODIM.EQ.1)THEN
        ALLOCATE(ZODRWIM(NODES,MCOMP))
        ZODRW = 0.
      ELSE
        IF(IADSORBIM.GT.0)THEN
          IF(IZODIM.EQ.2)THEN
             ALLOCATE(ZODRSIM(NODES,MCOMP))
             ZODRSIM = 0.
          ELSEIF(IZODIM.EQ.3)THEN
            ALLOCATE(ZODRWIM(NODES,MCOMP))
            ZODRWIM = 0.
            ALLOCATE(ZODRSIM(NODES,MCOMP))
            ZODRSIM = 0.
          ENDIF
        ENDIF  
        IF(IAW_ADSORBIM.GT.0) THEN
          IF(IZODIM.EQ.4)THEN
            ALLOCATE(ZODRWIM(NEQS,MCOMP))
            ZODRW = 0.
            ALLOCATE(ZODRSIM(NEQS,MCOMP))
            ZODRS = 0.
            ALLOCATE(ZODRAWIM(NEQS,MCOMP))
            ZODRAWIM = 0.           
          ENDIF            
        ENDIF  
      ENDIF
      IF(IFODIM.EQ.1)THEN
        ALLOCATE(FODRWIM(NODES,MCOMP))
        FODRWIM = 0.
      ELSE
        IF(IADSORBIM.GT.0)THEN
          IF(IFODIM.EQ.2)THEN
             ALLOCATE(FODRSIM(NODES,MCOMP))
             FODRSIM = 0.
          ELSEIF(IFODIM.EQ.3)THEN
            ALLOCATE(FODRWIM(NODES,MCOMP))
            FODRWIM = 0.
            ALLOCATE(FODRSIM(NODES,MCOMP))
            FODRSIM = 0.
          ENDIF
        ENDIF  
        IF(IAW_ADSORBIM.GT.0) THEN
          IF(IFODIM.EQ.4)THEN
            ALLOCATE(FODRWIM(NEQS,MCOMP))
            FODRW = 0.
            ALLOCATE(FODRSIM(NEQS,MCOMP))
            FODRS = 0.
            ALLOCATE(FODRAWIM(NEQS,MCOMP))
            FODRAW = 0.           
          ENDIF            
        ENDIF              
      ENDIF
C------------------------------------------------------------------
C5------FILL BOUNDARY ARRAY(ICBUNDIM).
      IF(IC_BNDIM_FLG.EQ.0)THEN
C5A-----READ ICBUNDIM ARRAY IF FLAG IS SET
        DO K = 1,NLAY
          KK = K
          NNDLAY = NODLAY(K)
          NSTRT = NODLAY(K-1)+1
          NDSLAY = NNDLAY - NODLAY(K-1)
          CALL U1DINT(ICBUNDIM(NSTRT),ANAME(1),NDSLAY,K,IN,IOUT)
        ENDDO
C
        DO N=1,NODES
          IF(ICBUND(N).EQ.0) ICBUNDIM(N) = 0
        ENDDO
      ELSE
C5B-------SET ICBUNDIM ARRAY TO ICBUND IF FLAG IS NOT SET
        DO N=1,NODES
          ICBUNDIM(N) = ICBUND(N)
        ENDDO
      ENDIF
C
C-----------------------------------------------------------------------
C6------READ MATERIAL PROPERTIES FOR IMMOBILE DOMAIN.
C-----------------------------------------------------------------------
C6A-----READ FRACTURE DOMAIN FRACTION INTO ARRAY PHIF IF FLOW WAS SINGLE DOMAIN
      IF(IDPF.EQ.0)THEN
        ALLOCATE(PHIF(NODES))
        DO K = 1,NLAY
          NNDLAY = NODLAY(K)
          NSTRT = NODLAY(K-1)+1
          NDSLAY = NNDLAY - NODLAY(K-1)
          CALL U1DREL(PHIF(NSTRT),ANAME(2),NDSLAY,K,IN,IOUT)
        ENDDO
      ENDIF
c
C6B-----READ IMMOBILE DOMAIN POROSITY INTO ARRAY PRSITYIM
      DO K = 1,NLAY
        NNDLAY = NODLAY(K)
        NSTRT = NODLAY(K-1)+1
        NDSLAY = NNDLAY - NODLAY(K-1)
        CALL U1DREL(PRSITYIM(NSTRT),ANAME(3),NDSLAY,K,IN,IOUT)
      ENDDO
C
C6C-----READ BULK DENSITY INTO ARRAY TEMP IF ADSORPTION.
      IF(IADSORBIM.NE.0.OR.IHEAT.EQ.1)THEN
        DO K = 1,NLAY
          NNDLAY = NODLAY(K)
          NSTRT = NODLAY(K-1)+1
          NDSLAY = NNDLAY - NODLAY(K-1)
          CALL U1DREL(BULKDIM(NSTRT),ANAME(4),NDSLAY,K,IN,IOUT)
        ENDDO
      ENDIF
C
C6D-----READ DISPERSION COEFFICIENT IF DP FLOW WAS ON
      IF(IDPF.NE.0.AND.IDISPIM.NE.0)THEN
        DO K = 1,NLAY
          NNDLAY = NODLAY(K)
          NSTRT = NODLAY(K-1)+1
          NDSLAY = NNDLAY - NODLAY(K-1)
          CALL U1DREL(DLIM(NSTRT),ANAME(5),NDSLAY,K,IN,IOUT)
        ENDDO
      ENDIF
C6E-----READ IMMOBILE DOMAIN MASS TRANSFER RATE INTO ARRAY DDTTR
      DO K = 1,NLAY
        NNDLAY = NODLAY(K)
        NSTRT = NODLAY(K-1)+1
        NDSLAY = NNDLAY - NODLAY(K-1)
        CALL U1DREL(DDTTR(NSTRT),ANAME(6),NDSLAY,K,IN,IOUT)
      ENDDO
C6F-----SET OR READ IMMOBILE DOMAIN SATURATIONS IF DUAL DOMAIN FLOW IS NOT SIMULATED      
      IF(IDPF.EQ.0)THEN
C6F1-----ALLOCATE ARRAY
        ALLOCATE(SnIM(NODES))
        SoIM => SnIM
        IF(IMSAT.EQ.1)THEN
C6F2---------IMMOBILE DOMAIN SATURATIONS SET TO UNITY            
          DO N=1,NODES
            SNIM(N) = 1.0
          ENDDO
        ELSEIF(IMSAT.EQ.2)THEN
C6F3---------IMMOBILE DOMAIN SATURATIONS SET SAME AS FOR MOBILE DOMAIN            
          DO N=1,NODES
            SNIM(N) = Sn(N)
          ENDDO  
        ELSEIF(IMSAT.EQ.3)THEN
C6F4---------IMMOBILE DOMAIN SATURATIONS ARE INPUT
          DO K=1,NLAY
            NNDLAY = NODLAY(K)
            NSTRT = NODLAY(K-1)+1
            NDSLAY = NNDLAY - NODLAY(K-1)
            CALL U1DREL(TEMPC(NSTRT),ANAME(14),NDSLAY,K,IN,IOUT)
          ENDDO
          DO N = 1, NODES
            SNIM(N) = TEMPC(N)
          ENDDO
        ENDIF
      ENDIF   
C
C-----------------------------------------------------------------------
C6G----READ HEAT PARAMETERS IF IHEAT IS ON
      IF(IHEAT.EQ.1) THEN
C6G1 -----READ AND PREPARE THE SOIL HEAT CAPACITY TERM
        DO K=1,NLAY
          NNDLAY = NODLAY(K)
          NSTRT = NODLAY(K-1)+1
          NDSLAY = NNDLAY - NODLAY(K-1)
          CALL U1DREL(TEMPC(NSTRT),ANAME(15),NDSLAY,K,IN,IOUT)
        ENDDO
C6G2-------PREPARE CONSTANT PART OF HEAT CAPACITY TERM
        DO N = 1, NODES
          ADSORBIM(N,MCOMPT) = TEMPC(N)
          ADSORBIM(N,MCOMPT) = ADSORBIM(N,MCOMPT)/RHOW/HTCAPW*BULKDIM(N)
        ENDDO
      ENDIF
C-----------------------------------------------------------------------
C6H -----READ AIR-WATER INTERFACE ADSORPTION PARAMETERS IF AW_ADSORB OPTION
      IF(IAW_ADSORBIM.EQ.1) THEN
        CALL AW_ADSORBIM1RP1  (IN)
      ENDIF
C ---------------------------------------------------------------------
C7------READ SPECIES DEPENDENT ARRAYS FOR ALL LAYERS.
C-----------------------------------------------------------------------
      ALLOCATE(TEMPC(Nodes))
      ALLOCATE(TEMPC8(Nodes))
      DO ICOMP=1,NTCOMP
        WRITE(IOUT,20) ICOMP
20      FORMAT(80('-')/1X,'THE FOLLOWING ARRAYS ARE READ FOR SPECIES ',
     *    1X,'NUMBER ',I3/80('-')) 
C8 -----------READ STOCHIOMETRY OF SPATIALLY VARIABLE
        IF (ICHAIN.NE.0.AND.ISPRCT.NE.0) THEN      
          NPAREN = NPARENT(ICOMP)
          DO JCOMP = 1,NPAREN
            WRITE(IOUT,21) ICOMP, JCOMP
21          FORMAT(1X,'*** GENERATION FACTOR OF DAUGHTER SPECIES ',
     1      I5,'FROM PARENT SPECIES ',I5,' FOR EACH LAYER ***')
            DO K=1,NLAY
              NNDLAY = NODLAY(K)
              NSTRT = NODLAY(K-1)+1
              NDSLAY = NNDLAY - NODLAY(K-1)
              CALL U1DREL(TEMPC(NSTRT),ANAME(4),NDSLAY,K,IN,IOUT)
            ENDDO              
            DO N = 1, NODES
              SPTLRCTIM(N,ICOMP,JCOMP) = TEMPC(N)
            ENDDO              
          ENDDO 
        ENDIF      
C --------------------------------------------------------
C9 --------SKIP ADSORPTION AND DECAY FOR IMMOBILE COMPONENTS
        IF(IHEAT.EQ.1.AND.ICOMP.EQ.MCOMPT) GO TO 101 !TEMPERATURE DOES NOT HAVE ADSORPTION OR DECAY OR PARENT
        IF(ICOMP .GT. MCOMPT) GO TO 101     !IMMOBILE COMPONENTS DO NOT HAVE ADSORPTION OR DECAY OR PARENT        
C---------------------------------------------------------
C10 -------IF AW_ADSORB OPTION IS ON THEN READ AW-PARTITIONING DATA
        IF(IAW_ADSORBIM.NE.0)THEN
          CALL AW_ADSORBIM1RP2(IN,ICOMP)
        ENDIF
C----------------------------------------------------------
C11------READ ADSORPTION COEFFICIENT FOR EACH SPECIES
      IF(IADSORBIM.NE.0)THEN
        DO K=1,NLAY
          NNDLAY = NODLAY(K)
          NSTRT = NODLAY(K-1)+1
          NDSLAY = NNDLAY - NODLAY(K-1)
          CALL U1DREL(TEMPC(NSTRT),ANAME(7),NDSLAY,K,IN,IOUT)
        ENDDO
C12-------PREPARE CONSTANT PART OF ADSORPTION TERM
        DO N = 1, NODES
          ADSORBIM(N,ICOMP) = TEMPC(N)
          ADSORBIM(N,ICOMP) = ADSORBIM(N,ICOMP) * BULKDIM(N)
        ENDDO
C---------------------------------------------------------
C13 -----READ ADSORPTION FREUNDLICH EXPONENT OR LANGMUIR COEFFICIENT FOR EACH SPECIES
        IF(IADSORBIM.EQ.2.OR.IADSORBIM.EQ.3)THEN
          DO K=1,NLAY
            NNDLAY = NODLAY(K)
            NSTRT = NODLAY(K-1)+1
            NDSLAY = NNDLAY - NODLAY(K-1)
            CALL U1DREL(TEMPC(NSTRT),ANAME(8),NDSLAY,K,IN,IOUT)
          ENDDO
          DO N = 1, NODES
            FLICHIM(N,ICOMP) = TEMPC(N)
          ENDDO
        ENDIF
      ENDIF
C---------------------------------------------------------
C10-----READ ZERO ORDER DECAY COEFFICIENTS
      IF(IZODIM.EQ.1)THEN
C10A---ZERO ORDER DECAY IN WATER IF FLAG IZODIM=1
        DO K=1,NLAY
          NNDLAY = NODLAY(K)
          NSTRT = NODLAY(K-1)+1
          NDSLAY = NNDLAY - NODLAY(K-1)
          CALL U1DREL(TEMPC(NSTRT),ANAME(9),NDSLAY,K,IN,IOUT)
        ENDDO
        DO N = 1, NODES
          ZODRWIM(N,ICOMP) = TEMPC(N)
        ENDDO
      ELSE
        IF(IADSORBIM.GT.0)THEN
          IF(IZODIM.EQ.2)THEN
C10B------ZERO ORDER DECAY ON SOIL IF FLAG IZODIM=2
            DO K=1,NLAY
              NNDLAY = NODLAY(K)
              NSTRT = NODLAY(K-1)+1
              NDSLAY = NNDLAY - NODLAY(K-1)
              CALL U1DREL(TEMPC(NSTRT),ANAME(10),NDSLAY,K,IN,IOUT)
            ENDDO
            DO N = 1, NODES
              ZODRSIM(N,ICOMP) = TEMPC(N)
            ENDDO
          ELSEIF(IZODIM.EQ.3)THEN
C10C-------ZERO ORDER DECAY IN WATER AND ON SOIL IF FLAG IZODIM=3
            DO K=1,NLAY
              NNDLAY = NODLAY(K)
              NSTRT = NODLAY(K-1)+1
              NDSLAY = NNDLAY - NODLAY(K-1)
              CALL U1DREL(TEMPC(NSTRT),ANAME(9),NDSLAY,K,IN,IOUT)
            ENDDO
            DO N = 1, NODES
              ZODRWIM(N,ICOMP) = TEMPC(N)
            ENDDO
            DO K=1,NLAY
              NNDLAY = NODLAY(K)
              NSTRT = NODLAY(K-1)+1
              NDSLAY = NNDLAY - NODLAY(K-1)
              CALL U1DREL(TEMPC(NSTRT),ANAME(10),NDSLAY,K,IN,IOUT)
            ENDDO
            DO N = 1, NODES
              ZODRSIM(N,ICOMP) = TEMPC(N)
            ENDDO
          ENDIF
        ENDIF
        IF(IAW_ADSORBIM.EQ.1)THEN 
          IF(IZODIM.EQ.4) THEN
C10C-------ZERO ORDER DECAY IN WATER, ON SOIL AND ON AIR-WATER INTERFACE IF FLAG IZODIM=4
            DO K=1,NLAY
              NNDLAY = NODLAY(K)
              NSTRT = NODLAY(K-1)+1
              NDSLAY = NNDLAY - NODLAY(K-1)
              CALL U1DREL(TEMPC(NSTRT),ANAME(9),NDSLAY,K,IN,IOUT)
            ENDDO
            DO N = 1, NODES
              ZODRWIM(N,ICOMP) = TEMPC(N)
            ENDDO
            DO K=1,NLAY
              NNDLAY = NODLAY(K)
              NSTRT = NODLAY(K-1)+1
              NDSLAY = NNDLAY - NODLAY(K-1)
              CALL U1DREL(TEMPC(NSTRT),ANAME(10),NDSLAY,K,IN,IOUT)
            ENDDO
            DO N = 1, NODES
              ZODRSIM(N,ICOMP) = TEMPC(N)
            ENDDO 
            DO K=1,NLAY
              NNDLAY = NODLAY(K)
              NSTRT = NODLAY(K-1)+1
              NDSLAY = NNDLAY - NODLAY(K-1)
              CALL U1DREL(TEMPC(NSTRT),ANAME(17),NDSLAY,K,IN,IOUT)
            ENDDO
            DO N = 1, NODES
              ZODRAWIM(N,ICOMP) = TEMPC(N)
            ENDDO 
          ENDIF    
        ENDIF    
      ENDIF
C---------------------------------------------------------
C11-----READ FIRST ORDER DECAY COEFFICIENTS
      IF(IFODIM.EQ.1)THEN
C11A----FIRST ORDER DECAY IN WATER IF FLAG IFODIM=1
        DO K=1,NLAY
          NNDLAY = NODLAY(K)
          NSTRT = NODLAY(K-1)+1
          NDSLAY = NNDLAY - NODLAY(K-1)
          CALL U1DREL(TEMPC(NSTRT),ANAME(11),NDSLAY,K,IN,IOUT)
        ENDDO
        DO N = 1, NODES
          FODRWIM(N,ICOMP) = TEMPC(N)
        ENDDO
      ELSE
        IF(IADSORBIM.GT.0)THEN
          IF(IFODIM.EQ.2)THEN
C11B------FIRST ORDER DECAY ON SOIL IF FLAG IFODIM=2
            DO K=1,NLAY
              NNDLAY = NODLAY(K)
              NSTRT = NODLAY(K-1)+1
              NDSLAY = NNDLAY - NODLAY(K-1)
              CALL U1DREL(TEMPC(NSTRT),ANAME(12),NDSLAY,K,IN,IOUT)
            ENDDO
            DO N = 1, NODES
              FODRSIM(N,ICOMP) = TEMPC(N)
            ENDDO
          ELSEIF(IFODIM.EQ.3)THEN
C11C-----FIRST ORDER DECAY IN WATER AND ON SOIL IF FLAG IFODIM=3
            DO K=1,NLAY
              NNDLAY = NODLAY(K)
              NSTRT = NODLAY(K-1)+1
              NDSLAY = NNDLAY - NODLAY(K-1)
              CALL U1DREL(TEMPC(NSTRT),ANAME(11),NDSLAY,K,IN,IOUT)
            ENDDO
            DO N = 1, NODES
              FODRWIM(N,ICOMP) = TEMPC(N)
            ENDDO
            DO K=1,NLAY
              NNDLAY = NODLAY(K)
              NSTRT = NODLAY(K-1)+1
              NDSLAY = NNDLAY - NODLAY(K-1)
              CALL U1DREL(TEMPC(NSTRT),ANAME(12),NDSLAY,K,IN,IOUT)
            ENDDO
            DO N = 1, NODES
              FODRSIM(N,ICOMP) = TEMPC(N)
            ENDDO
          ENDIF
          ENDIF
        IF(IAW_ADSORBIM.EQ.1)THEN 
          IF(IFODIM.EQ.4) THEN
C10C-------ZERO ORDER DECAY IN WATER, ON SOIL AND ON AIR-WATER INTERFACE IF FLAG IFODIM=4
            DO K=1,NLAY
              NNDLAY = NODLAY(K)
              NSTRT = NODLAY(K-1)+1
              NDSLAY = NNDLAY - NODLAY(K-1)
              CALL U1DREL(TEMPC(NSTRT),ANAME(9),NDSLAY,K,IN,IOUT)
            ENDDO
            DO N = 1, NODES
              FODRWIM(N,ICOMP) = TEMPC(N)
            ENDDO
            DO K=1,NLAY
              NNDLAY = NODLAY(K)
              NSTRT = NODLAY(K-1)+1
              NDSLAY = NNDLAY - NODLAY(K-1)
              CALL U1DREL(TEMPC(NSTRT),ANAME(10),NDSLAY,K,IN,IOUT)
            ENDDO
            DO N = 1, NODES
              FODRSIM(N,ICOMP) = TEMPC(N)
            ENDDO 
            DO K=1,NLAY
              NNDLAY = NODLAY(K)
              NSTRT = NODLAY(K-1)+1
              NDSLAY = NNDLAY - NODLAY(K-1)
              CALL U1DREL(TEMPC(NSTRT),ANAME(18),NDSLAY,K,IN,IOUT)
            ENDDO
            DO N = 1, NODES
              FODRAWIM(N,ICOMP) = TEMPC(N)
            ENDDO 
          ENDIF    
        ENDIF              
      ENDIF
C
C---------------------------------------------------------
C12------CONCENTRATION OF EACH SPECIES
101   CONTINUE
      CNAME = ANAME(13)
      IF(IHEAT.EQ.1.AND.ICOMP.EQ.MCOMPT) CNAME = ANAME(16)      
      IF(IDPIN.EQ.0) THEN !--------------------------------SINGLE PRECISION READ
      DO K=1,NLAY
        NNDLAY = NODLAY(K)
        NSTRT = NODLAY(K-1)+1
        NDSLAY = NNDLAY - NODLAY(K-1)
        CALL U1DREL(TEMPC(NSTRT),CNAME,NDSLAY,K,IN,IOUT)
      ENDDO
      DO N = 1,NODES
        CONCIM(N,ICOMP) = TEMPC(N)
        IF(ICBUNDIM(N).EQ.0) CONCIM(N,ICOMP)=CINACT
        CONCOIM(N,ICOMP) = CONCIM(N,ICOMP)
      ENDDO
      ELSE !--------------------------------DOUBLE PRECISION READ
      DO K=1,NLAY
        NNDLAY = NODLAY(K)
        NSTRT = NODLAY(K-1)+1
        NDSLAY = NNDLAY - NODLAY(K-1)
        CALL U1DREL8(TEMPC8(NSTRT),CNAME,NDSLAY,K,IN,IOUT)
      ENDDO
      DO N = 1,NODES
        CONCIM(N,ICOMP) = TEMPC8(N)
        IF(ICBUNDIM(N).EQ.0) CONCIM(N,ICOMP)=CINACT
        CONCOIM(N,ICOMP) = CONCIM(N,ICOMP)
      ENDDO
      ENDIF
C
C---------------------------------------------------------
      ENDDO  ! END SPECIES DO LOOP ICOMP
C---------------------------------------------------------
C13   DEALLOCATE UNWANTED ARRAYS
      IF(IADSORBIM.NE.0.OR.IHEAT.EQ.1) DEALLOCATE(BULKDIM)
      DEALLOCATE (TEMPC)
      DEALLOCATE (TEMPC8)
C-----------------------------------------------------------------------
C14-----PREPARE ARRAYS OF PARAMETERS
C-----------------------------------------------------------------------
C15-----MULTIPLY DUAL DOMAIN RATE TERM BY VOLUME OF GRID BLOCK TO GIVE TOTAL TRANSFER FOR CELL
      CALL SGWF2LPFU1SC(DDTTR(1),1)
C
C-----------------------------------------------------------------------
C16-----IF FLOW IS DUAL POROSITY THEN CONDUCTANCE AND STORAGE TERMS ARE ADJUSTED IN DPT. OTHERWISE, 
C ------ADJUST CONDUCTANCE AND STORAGE HERE, IF ALSO IFRAHK=1 (PROPERTIES WERE INPUT FOR FRACTURE DOMAIN)
      IF(IDPF.EQ.0)THEN
C16B-----SCALE STORAGE COEFFICIENTS BY MOBILE FRACTION IF TRANSIENT FLOW
        ITISTR = 0
        DO NPE=1,NPER
          IF(ISSFLG(NPE).EQ.0) THEN 
            ITISTR = 1  
            GO TO 123
          ENDIF  
        ENDDO    
 123    CONTINUE
C          
        IF(ITISTR.EQ.0) GO TO 211  
C --------SCALE SC1 AND SC2 ONLY IF IFRAHK = 1 AND GIVEN FOR FRACTURE WHERE FLOW OCCURS        
        IF(IFRAHK.EQ.1)THEN 
        DO N=1,NODES
          SC1(N) = SC1(N) * PHIF(N)
        ENDDO
        NCNVRT=0
        DO K = 1,NLAY
          IF(LAYCON(K).NE.0)THEN
            NCNVRT=1
          ENDIF
        ENDDO  
        IF(NCNVRT.GT.0)THEN
          DO N=1,NODES
            SC2(N) = SC2(N) * PHIF(N)
          ENDDO
          IF(IDPF.NE.0)THEN
            DO N=1,NODES
              SC2IM(N) = SC2IM(N) * (1.0 - PHIF(N))
            ENDDO
C          ELSE
C            ALLOCATE(SC2IM(NODES))
C            DO N=1,NODES
C              SC2IM(N) = 0.0
C            ENDDO
          ENDIF  
        ENDIF
        ELSE
C4B ----SPECIFIC YIELD IS FOR TOTAL DOMAIN AND EFFECTIVE POROSITY IS FOR FRACTURE 
C4B-----CHECK FOR ERRORS, POROSITY CANNOT BE LESS THAN SPECIFIC YIELD IN FRACTURE 
        IEFFPOR = 0
        DO K = 1,NLAY
        NNDLAY = NODLAY(K)
        NSTRT = NODLAY(K-1)+1
        NDSLAY = NNDLAY - NODLAY(K-1)
        DO N = NSTRT,NNDLAY 
          IF(LAYCON(K).NE.3.AND.LAYCON(K).NE.2.AND.LAYCON(K).NE.4
     1       .AND.LAYCON(K).NE.5)GO TO 90
          IF(ICBUND(N).EQ.0.OR.IBOUND(N).EQ.0)GO TO 90
          SYIF = SC2(N)/AREA(N) / PHIF(N)
          IF(PRSITY(N)+1.e-5.LT.SYIF) THEN
              PRSITY(N) = SYIF+1.e-5 
              IEFFPOR = 1
c            WRITE(IOUT,56) N, PRSITY(N),SC2(N)/AREA(N),PHIF(N) 
c56          FORMAT(5X,'*** Porosity * phif < Sy FOR NODE: ',I10,' ***'/
c     1      9X,'POROSITY IS ',G10.4,' SC2 IS',G10.4,' PHIF IS',G10.4)  
c            STOP
          ENDIF
90        CONTINUE
        ENDDO  
        ENDDO
        IF(IEFFPOR.EQ.1) THEN 
          WRITE(IOUT,57) 
57        FORMAT(5X,'*** EFFECTIVE POROSITY INPUT IS ASSUMED TO BE FOR',
     1    1X,'TOTAL DOMAIN AND IS SCALED FOR FRACTURE DOMAIN ***' /
     2    9X,'POROSITY OF FRACURE HAS TO BE LARGER THAN ITS SPECIFIC',
     2    1X,'YIELD WHICH WAS INPUT FOR TOTAL DOMAIN (IFRACHK = 0)')
        ENDIF    
        ENDIF
  211   CONTINUE   
C16C----SCALE CONDUCTIVITY IF PROVIDED FOR FRACTURES ONLY 
        IF(IFRAHK.EQ.1) THEN
        DO K=1,NLAY
          NNDLAY = NODLAY(K)
          NSTRT = NODLAY(K-1)+1
C1B---------loop over all nodes within each layer
          DO N=NSTRT,NNDLAY
C2------------loop over all connections of node N and fill upper triangle with PGF term
            DO II = IA(N)+1,IA(N+1)-1
              JJ = JA(II)
C3--------------only for upper triangle of porous medium nodes
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
      ENDIF  
C-----------------------------------------------------------------------      
C17----FOR ZERO ORDER DECAY IN WATER
      IF(IZOD.EQ.1.OR.IZOD.EQ.3.OR.IZOD.EQ.4)THEN
        DO ICOMP = 1,MCOMP
        DO N = 1,NODES
          ZODRW(N,ICOMP) = ZODRW(N,ICOMP) * PHIF(N)
        ENDDO
        ENDDO
      ENDIF
      IF(IZODIM.EQ.1.OR.IZODIM.EQ.3.OR.IZODIM.EQ.4)THEN
        DO ICOMP = 1,MCOMP
        DO N = 1,NODES
          ZODRWIM(N,ICOMP)=ZODRWIM(N,ICOMP)*PRSITYIM(N)*(1.0 - PHIF(N))
        ENDDO
        ENDDO
      ENDIF      
C18----FOR ZERO ORDER DECAY IN SOIL ** DOES NOT NEED PHIF SINCE ADSORB TERM HAS IT FOR M AND IM
C      IF((IZOD.EQ.2.OR.IZOD.EQ.3).AND.IADSORB.NE.0)THEN
C        DO ICOMP = 1,MCOMP
C        DO N = 1,NODES
C          ZODRS(N,ICOMP) = ZODRS(N,ICOMP) * PHIF(N)
C        ENDDO
C        ENDDO
C      ENDIF
C      IF((IZODIM.EQ.2.OR.IZODIM.EQ.3).AND.IADSORBIM.NE.0)THEN
C        DO ICOMP = 1,MCOMP
C        DO N = 1,NODES
C          ZODRSIM(N,ICOMP) = ZODRSIM(N,ICOMP) * (1.0 - PRSITYIM(N))
C     1          *(1.0 - PHIF(N))
C        ENDDO
C        ENDDO
C      ENDIF

C19----FOR FIRST ORDER DECAY IN WATER
      IF(IFOD.EQ.1.OR.IFOD.EQ.3.OR.IFOD.EQ.4)THEN
        DO ICOMP = 1,MCOMP
        DO N = 1,NODES
          FODRW(N,ICOMP) = FODRW(N,ICOMP) * PHIF(N)
        ENDDO
        ENDDO
      ENDIF
      IF(IFODIM.EQ.1.OR.IFODIM.EQ.3.OR.IFOD.EQ.4)THEN
        DO ICOMP = 1,MCOMP
        DO N = 1,NODES
          FODRWIM(N,ICOMP)=FODRWIM(N,ICOMP)*PRSITYIM(N)*(1.0 - PHIF(N))
        ENDDO
        ENDDO
      ENDIF      
C20----FOR FIRST ORDER DECAY ON SOIL ** DOES NOT NEED PHIF SINCE ADSORB TERM HAS IT FOR M AND IM 
C      IF((IFOD.EQ.2.OR.IFOD.EQ.3).AND.IADSORB.NE.0)THEN
C        DO ICOMP = 1,MCOMP
C        DO N = 1,NODES
C          FODRS(N,ICOMP) = FODRS(N,ICOMP) * PHIF(N)
C          FODRSIM(N,ICOMP) = FODRSIM(N,ICOMP) * (1.0 - PHIF(N))
C        ENDDO
C        ENDDO
C      ENDIF
C21----FOR POROSITY FRACTIONED INTO MOBILE AND IMMOBILE DOMAINS
      DO N = 1,NODES
        PRSITY(N) = PRSITY(N) * PHIF(N)
        PRSITYIM(N) = PRSITYIM(N) * (1.0 - PHIF(N))
      ENDDO
C22----FOR KD VALUES ON SOIL
      IF(IADSORB.NE.0)THEN
        DO ICOMP = 1,MCOMP
        DO N = 1,NODES
          ADSORB(N,ICOMP) = ADSORB(N,ICOMP) * PHIF(N)
        ENDDO
        ENDDO
      ENDIF
      IF(IADSORBIM.NE.0)THEN
        DO ICOMP = 1,MCOMP
        DO N = 1,NODES
          ADSORBIM(N,ICOMP) = ADSORBIM(N,ICOMP) * (1.0 - PHIF(N))
        ENDDO
        ENDDO
      ENDIF      
C  ---------------------------------------------------------------------
C10------RETURN
      RETURN
      END
C
C-----------------------------------------------------------------------
      SUBROUTINE INITMASSIM 
C     ******************************************************************
C      COMPUTE THE INITIAL MASS OF IMMOBILE COMPONENTS
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL, ONLY:IOUT,NODES,IDPT,IDPF,TOP,BOT,AREA,ISSFLG,IDPF
      USE GWTDPTMODULE
      USE GWTBCTMODULE, ONLY: MCOMP,ICBUND,CINACT,
     1    ISOLUBILITY,SOLLIM,SOLSLOPE
      USE GWFDPFMODULE, ONLY: PHIF,SoIM,SC1IM,SC2IM,HOLDIM
      USE AW_ADSORBIMMODULE, ONLY: AK_AWIIM, AREA_AWIIM
      DOUBLE PRECISION MASLOCW,MASLOCS,MASLOCC,REFHD,MASLOPR,MASLOAWIIM
C     --------------------------------------------------------------------
C1-----COMPUTE MASS FOR MOBILE COMPONENTS OF IMMOBILE DOMAIN
      WRITE(IOUT,12)
12    FORMAT(/5X,'IN IMMOBILE DOMAIN'/
     *        5X,'------------------'/
     *       5X,'COMPONENT INDEX',5X,'MASS'/
     *       5X,'---------------',5X,'----')
      DO ICOMP = 1,MCOMP
        MASSBCTIM (ICOMP) = 0.0
      ENDDO
      MASLOCW = 0.0
      MASLOCS = 0.0
      MASLOCC = 0.0
      MASLOPR = 0.0
      MASLOAWIIM = 0.0
      ISS=ISSFLG(1)
C1A-------INITIALIZE A-W SURFACE AREA AND K_AW       
      IF(IAW_ADSORBIM.EQ.1)THEN       
        CALL AWI_AREAIM
        DO ICOMP = 1,MCOMP
          CALL AWI_GETKAWIIM(ICOMP)      
        ENDDO  
      ENDIF        
C2-----COMPUTE MASS FOR EACH ACTIVE NODE
      DO N = 1,NODES
        IF(ICBUNDIM(N).EQ.0) CYCLE
C3-------COMPUTE VOLUME OF NODE
        ALENG = TOP(N) - BOT(N)
        VOLU = AREA(N) * ALENG
C4-------FOR EACH COMPONENT
        DO ICOMP = 1,MCOMP
C5-------COMPUTE MASS IN WATER
          IF(ISS.EQ.1.OR.IDPF.EQ.0)THEN
            MASLOCW = PRSITYIM(N) * SoIM(N) * CONCIM(N,ICOMP) * VOLU
          ELSE
C6---------INITIAL MASS IN PORE STORAGE
            MASLOCW=(PRSITYIM(N)*VOLU-SC2IM(N)*ALENG+
     *        SC2IM(N)*ALENG*SoIM(N))*CONCIM(N,ICOMP)
C7-----------INITIAL MASS IN COMPRESSIBLE STORAGE ASSUMED FROM TOP OF CELL
            REFHD = HOLDIM(N) - TOP(N)
            IF(REFHD.LT.0.0) REFHD = 0.0
            MASLOCC = SC1IM(N) * REFHD * CONCIM(N,ICOMP)
          ENDIF
C8---------COMPUTE MASS ON SOIL
          IF(IADSORBIM.EQ.1)THEN  ! LINEAR
            MASLOCC = ADSORBIM(N,ICOMP) * CONCIM(N,ICOMP) * VOLU  
          ELSEIF(IADSORBIM.EQ.2)THEN ! FREUNDLICH    
            ETA = FLICHIM(N,ICOMP)
            MASLOCC = ADSORBIM(N,ICOMP) * CONCIM(N,ICOMP)**ETA * VOLU
          ELSEIF(IADSORBIM.EQ.3)THEN ! FREUNDLICH    
            ETA = FLICHIM(N,ICOMP)
            MASLOCC = ADSORBIM(N,ICOMP) * VOLU * * CONCIM(N,ICOMP)
            MASLOCC = MASLOCC / (1.0 + ETA * CONCIM(N,ICOMP))
          ENDIF
C9-----------COMPUTE MASS ON AIR-WATER INTERFACE
          IF(IAW_ADSORBIM.EQ.1)THEN   
            MASLOAWIIM= AREA_AWIIM(N)* AK_AWIIM(N,ICOMP)
     1        *CONCIM(N,ICOMP)* VOLU
          ENDIF
C10-----------COMPUTE MASS PRECIPITATED ABOVE SOLUBILITY LIMIT
          IF(ISOLUBILITY.NE.0)THEN 
            MASLOPR = 0.0
            CW = CONCIM(N,ICOMP)
            SLIM = SOLLIM(ICOMP)
            SSLOPE = SOLSLOPE(ICOMP)
            IF(CW.GT.SLIM) MASLOPR = (CW-SLIM) * SSLOPE              
            MASLOPR =  MASLOPR * VOLU
          ENDIF                      
C 11----------ADD UP ALL TERMS 
          MASSBCTIM(ICOMP) = MASSBCTIM(ICOMP)+MASLOCW+MASLOCS+MASLOCC 
     *     + MASLOPR  + MASLOAWIIM   
        ENDDO
      ENDDO
C11-----WRITE MASS TO OUTPUT FILE
      DO ICOMP = 1,MCOMP
        WRITE(IOUT,13) ICOMP, MASSBCTIM(ICOMP)
      ENDDO
13    FORMAT(I15,5X,E15.6)
C  ---------------------------------------------------------------------
C10------RETURN
      RETURN
      END
C
C-----------------------------------------------------------------------
      SUBROUTINE GWT2DPTU1FM(KPER,ICOMP,ISS)
C     ******************************************************************
C     FORMULATE DUAL POROSITY TERMS FOR TRANSPORT
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,  ONLY:NODES,NLAY,IUNSAT,NODLAY,AMAT,IA,HNEW,AKR,
     1  TOP,BOT,IBOUND,ISSFLG,Sn,IDPF
      USE GWTDPTMODULE
      USE GWFDPFMODULE, ONLY: HNEWIM,SnIM,CBCFIM
      USE GWTBCTMODULE, ONLY: ICBUND,IHEAT,MCOMPT,ICHAIN,ISOLUBILITY
      USE GWFBASMODULE, ONLY: DELT
      USE GWFBCFMODULE, ONLY: ISFAC,SC1,SC2,LAYCON
      USE SMSMODULE, ONLY: NONMETH,DKDH
      DOUBLE PRECISION DDFLOW,RHO,TLED,TTOP,BBOT, TOTTHICK,EPS,HD,
     *  DS,RHO2,RHO1,THCK,SATN,SATO,FTERM,DFTERM,SW,EKR,AKRCIMN,
     *  DIAGS,RTERMS,DIAGD,RTERMD,QIJ,UPSAT,DTERMS
C
C     ------------------------------------------------------------------
C1--------INITALIZE ARRAYS
      DO N=1,NODES
        DIADDT(N) = 0.0
        RDDT(N) = 0.0
        OFFDDTM(N) = 0.0
        OFFDDTIM(N) = 0.0
      ENDDO
C-----------------------------------------------------------------------------
C2------FILL STORAGE TERM ON DIAGONAL AND RHS OF IMMOBILE DOMAIN EQUATION
C-----------------------------------------------------------------------------
      TLED=1.0D0/DELT
      DO 200 K=1,NLAY
        NNDLAY = NODLAY(K)
        NSTRT = NODLAY(K-1)+1
C-----------------------------------------------------------------------------
        DO 140 N=NSTRT,NNDLAY
          IF(ICBUNDIM(N).EQ.0) CYCLE
          CALL GWT2DPT1STOIM(N,ICOMP,DIAGS,RTERMS,ISS)
          DIADDT(N)=DIADDT(N) + DIAGS
          RDDT(N)=RDDT(N) + RTERMS 
  140   CONTINUE
200   CONTINUE
C-----------------------------------------------------------------------------
C2------FILL DECAY TERM ON DIAGONAL AND RHS OF IMMOBILE DOMAIN EQUATION
C-----------------------------------------------------------------------------
      IF(IHEAT.EQ.1. AND. ICOMP.GE.MCOMPT) GO TO 52 ! NO DECAY FOR HEAT OR IMMOBILE COMPONENTS     
      DO 201 K=1,NLAY
        NNDLAY = NODLAY(K)
        NSTRT = NODLAY(K-1)+1
C-----------------------------------------------------------------------------
        DO 141 N=NSTRT,NNDLAY
          IF(ICBUNDIM(N).EQ.0) CYCLE
          CALL GWT2DPT1DCYIM(N,ICOMP,DIAGD,RTERMD)
          DIADDT(N)=DIADDT(N) + DIAGD
          RDDT(N)=RDDT(N) + RTERMD
  141   CONTINUE
201   CONTINUE
52    CONTINUE      
C-----------------------------------------------------------------------
C5I---------FILL PARENT GENERATION TERMS ON RHS
      IF(ICHAIN.NE.0) THEN
        IF(IHEAT.EQ.1. AND. ICOMP.EQ.MCOMPT) GO TO 53 ! NO GENERATION FOR HEAT FROM COMPONENT DECAY         
        DO N=1,NODES
          IF(ICBUNDIM(N).NE.0)THEN
            CALL GWT2BCT1GENIM(N,ICOMP)
          ENDIF
        ENDDO
53      CONTINUE        
      ENDIF     
C-----------------------------------------------------------------------
C5J---------FILL SOLUBILITY PRECIPITATION TERMS 
      IF(ISOLUBILITY.NE.0) THEN
        IF(IHEAT.EQ.1. AND. ICOMP.EQ.MCOMPT) GO TO 54
        DO N=1,NODES
          IF(ICBUNDIM(N).NE.0)THEN
            CALL SOLUBILITYIM(N,ICOMP,DTERMS,RTERMS,ISS)
            DIADDT(N) = DIADDT(N) + DTERMS
            RDDT(N) = RDDT(N) + RTERMS            
          ENDIF
        ENDDO
54      CONTINUE
      ENDIF           
C----------------------------------------------------------------------
C5K-------FILL PARTITION TERM FOR AIR-WATER INTERFACE IF ADSORPTION ON INTERFACE
        IF(IAW_ADSORBIM.EQ.1)THEN
          IF(IHEAT.EQ.1. AND. ICOMP.GE.MCOMPT) GO TO 55
          CALL AWI_KAWIIM (ICOMP)
55        CONTINUE
        ENDIF
C
C-----------------------------------------------------------------------------
C5------FILL ADVECTIVE/DISPERSIVE MASS TRANSFER TERM ON OFF-DIAGONALS AND DIAGONALS
C-----------------------------------------------------------------------------
      IF(IDPF.NE.0) THEN  !DUAL POROSITY FLOW
        DO 220 N = 1, NODES
          IF(ICBUNDIM(N).EQ.0) CYCLE
C---------FILL UPSTREAM TERM
          QIJ = CBCFIM(N)
          IF(QIJ.LE.0)THEN ! MOBILE DOMAIN IS UPSTREAM
            AMAT(IA(N)) = AMAT(IA(N)) + QIJ
            OFFDDTIM(N) = OFFDDTIM(N) - QIJ
          ELSE ! IMMOBILE DOMAIN IS UPSTREAM
            DIADDT(N) = DIADDT(N) + QIJ
            OFFDDTM(N) = OFFDDTM(N) - QIJ
          ENDIF
220     CONTINUE
      ENDIF
C-----------------------------------------------------------------------------
C5------FILL DIFFUSIVE MASS TRANSFER TERM ON OFF-DIAGONALS AND DIAGONALS
C-----------------------------------------------------------------------------
      DO 300 K=1,NLAY
        NNDLAY = NODLAY(K)
        NSTRT = NODLAY(K-1)+1
        IF(LAYCON(K).EQ.0)THEN
C-----------------------------------------------------------------------------
C6-------NON-CONVERTIBLE LAYER, SO COPY CONSTANT TERM TO OFF-DIAGONALS AND ADD TO DIAGONALS
          DO 240 N=NSTRT,NNDLAY
            IF(ICBUNDIM(N).NE.0.AND.ICBUND(N).NE.0) THEN
              DDFLOW = DDTTR(N)
              OFFDDTM(N) = OFFDDTM(N) + DDFLOW
              OFFDDTIM(N) = OFFDDTIM(N) + DDFLOW
              AMAT(IA(N)) =  AMAT(IA(N)) - DDFLOW
              DIADDT(N) = DIADDT(N) - DDFLOW
            ENDIF
240       CONTINUE
        ELSE
C-----------------------------------------------------------------------------
C7-------CONVERTIBLE LAYER IS NONLINEAR, COMPUTE NONLINEAR TERMS
          DO 250 N=NSTRT,NNDLAY
            IF(ICBUNDIM(N).EQ.0.OR.ICBUND(N).EQ.0) CYCLE
C
C8------------FIND UPSTREAM SATURATION BETWEEN M AND IM DOMAINS
            IF(Sn(N).GT.SnIM(N))THEN
C8A-------------MOBILE DOMAIN IS UPSTREAM
              UPSAT = Sn(N)
            ELSE
C8B-------------IMMOBILE DOMAIN IS UPSTREAM
              UPSAT = SnIM(N)
            ENDIF
C
C9--------------FILL NONLINEAR TERMS IN OFF-DIAGONALS AND ADD TO DIAGONALS
            DDFLOW = DDTTR(N) * UPSAT
            OFFDDTM(N) = OFFDDTM(N) + DDFLOW
            OFFDDTIM(N) = OFFDDTIM(N) + DDFLOW
            AMAT(IA(N)) =  AMAT(IA(N)) - DDFLOW
            DIADDT(N) = DIADDT(N) - DDFLOW
250       CONTINUE
        ENDIF
300   CONTINUE
      itestmat = 0
      if(itestmat.eq.1)then
c        open(99,file='mat_USGs.TXT')
        WRITE(99,*)'NODE, offddtm(n),offddtim(n),amat(ia(n)),diaddt(n):'
        DO N=1,Nodes
          if(ibound(n). eq.0) cycle
          WRITE(99,66)N,offddtm(n),offddtim(n),amat(ia(n)),diaddt(n)
        ENDDO
c66      FORMAT(I9,1X,G15.3,2X,100G15.3)
66      FORMAT(I9,1X,G15.4,2X,100G15.4)
        CLOSE (99)
        stop
      endif      
C
C11------RETURN
      RETURN
        END
C----------------------------------------------------------------------
      SUBROUTINE GWT2DPT1STOIM(N,ICOMP,DTERMS,RTERMS,ISS)
C     ******************************************************************
C     FORMULATE STORAGE TERM FOR EACH NODE AND COMPONENT IN IMMOBILE DOMAIN
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,   ONLY:IOUT,NCOL,NROW,NLAY,IBOUND,TOP,NODLAY,ISYM,
     1 HNEW,AMAT,RHS,HOLD,ISSFLG,IA,JA,NJA,NODES,AREA,BOT,TOP,Sn,So
      USE GWTBCTMODULE, ONLY: ICT,IHEAT,MCOMPT
      USE GWFDPFMODULE, ONLY: PHIF,SnIM,SoIM
      USE CLN1MODULE, ONLY: ACLNNDS
      USE GWFBASMODULE,ONLY:DELT
      USE GWTDPTMODULE, ONLY: CONCIM,CONCOIM,IADSORBIM,PRSITYIM,
     * ADSORBIM,FLICHIM,IZODIM,IFODIM,ZODRWIM,FODRWIM,
     1  ZODRSIM,FODRSIM
      DOUBLE PRECISION VODT,ADSTERM,CW,FL,CWO,DT,RT,ALENG,VOLU,
     1  X,Y,EPSS,EPS,CT,CEPS,QA,QEPS,DQ,OMP,DTERMS,RTERMS
C     ------------------------------------------------------------------
      DTERMS = 0.0
      RTERMS = 0.0
      IF(N.LE.NODES)THEN
        ALENG = TOP(N) - BOT(N)
      ELSE
        ALENG = ACLNNDS(N-NODES,5)
      ENDIF
      VOLU = AREA(N) * ALENG
      VODT = VOLU / DELT
C----------------------------------------------------------------------
      IF(ICT.EQ.0)THEN  !----------WATER PHASE CONCENTRATION FORMULATION
C-------ADD STORAGE TERMS ON SOIL
        IFILLADS = 0
        IF(IHEAT.EQ.1.AND.ICOMP.EQ.MCOMPT) IFILLADS = 1
        IF(N.LE.NODES.AND.IADSORBIM.EQ.1) IFILLADS = 1
        IF(N.LE.NODES.AND.IADSORBIM.EQ.2) IFILLADS = 2
        IF(N.LE.NODES.AND.IADSORBIM.EQ.3) IFILLADS = 3
        IF(N.GT.NODES) IFILLADS = 0
C          
        IF(IFILLADS.EQ.1)THEN
C---------LINEAR ISOTHERM
          ADSTERM = ADSORBIM(N,ICOMP) * VODT
          DTERMS = DTERMS - ADSTERM
          RTERMS = RTERMS - ADSTERM * CONCOIM(N,ICOMP)
C-----------------------------------------------------------------------
        ELSEIF(IFILLADS.EQ.2.OR.IFILLADS.EQ.3)THEN
C---------NONLINEAR ISOTHERM FILLED AS MODIFIED PICARD
          ADSTERM = ADSORBIM(N,ICOMP) * VODT
          FL = FLICHIM(N,ICOMP)
          CW = CONCIM(N,ICOMP)
          IF(CW.LT.0.0) CW = 0.0
          CWO = CONCOIM(N,ICOMP)
          IF(CWO.LT.0.0) CWO = 0.0
          DT = 0.0
          RT = 0.0
          IF(IFILLADS.EQ.2)THEN
            RT = CW**FL - CWO**FL 
            RT = RT * ADSTERM 
            if( fl.lt.1.0 .and. cw. lt. 1.e-7) cw = 1e-7        !No NWT term if/when slope approaches infinity
            DT = ADSTERM*FL * CW**(FL - 1.0)
          ELSEIF(IFILLADS.EQ.3) THEN
            RT = CW/(1.0 + FL*CW) - CWO/(1.0 + FL*CWO)
            RT = RT *ADSTERM
            DT = (1+FL*CW)*ADSTERM - ADSTERM*CW*FL ! ANALYTIC DERIV OF U/V = (VdU - UdV)/V^2
            DT = DT / (1.0 + FL*CW)**2  
          ENDIF    
          DTERMS = DTERMS - DT
          RTERMS = RTERMS - DT * CW + RT
        ENDIF
C-----------------------------------------------------------------------
C-------ADD STORAGE TERM IN WATER
        CALL GWT2BCT1STOIMW (N,ICOMP,DTERMS,RTERMS,VODT,VOLU,ALENG,ISS)
      ELSE       !-----------------------TOTAL CONCENTRATION FORMULATION
C-------NET STORAGE TERM FOR TOTAL CONCENTRATION FORMULATION
        DTERMS = DTERMS - VODT
        RTERMS = RTERMS - VODT * CONCOIM(N,ICOMP)
      ENDIF
C
C9------RETURN
      RETURN
      END
C---------------------------------------------------------------------------
      SUBROUTINE GWT2BCT1STOIMW(N,ICOMP,DTERMS,RTERMS,VODT,VOLU,ALENG,
     *  ISS)
C     ******************************************************************
C     FORMULATE STORAGE TERM FOR EACH NODE AND COMPONENT IN WATER
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,   ONLY:IOUT,NCOL,NROW,NLAY,IBOUND,TOP,NODLAY,ISYM,
     1 ISSFLG,IA,JA,NJA,NODES,AREA,BOT,TOP,LAYNOD,IDPF
      USE CLN1MODULE, ONLY: ACLNNDS
      USE GWFBASMODULE,ONLY:DELT
      USE GWFBCFMODULE,ONLY:LAYCON
      USE GWTBCTMODULE, ONLY: MCOMP,CONC,CONCO,ICT,IADSORB,PRSITY,
     * ADSORB
      USE GWFDPFMODULE, ONLY: PHIF,SnIM,SoIM,SC1IM,SC2IM,HNEWIM,HOLDIM
      USE GWTDPTMODULE, ONLY: CONCIM,CONCOIM,IADSORBIM,PRSITYIM,
     * ADSORBIM
      DOUBLE PRECISION VODT,ADSTERM,CW,FL,CWO,DT,RT,ALENG,VOLU,
     1  X,Y,EPSS,EPS,CT,CEPS,QA,QEPS,DQ,OMP,DTERMS,RTERMS,SATN,SATO
C-------------------------------------------------------------------------------------
      SATN = SnIM(N)
      SATO = SoIM(N)
      if(satn.lt.1.0e-3) satn = 1.0e-3
      if(sato.lt.1.0e-3) sato = 1.0e-3      
      IF(ISS.EQ.1.OR.N.GT.NODES.OR.IDPF.EQ.0) THEN
C-------USE TRADITIONAL STORAGE TERM FOR CLN DOMAIN OR SS FLOW
        IF(N.LE.NODES) VODT = VODT * PRSITYIM(N)
        DTERMS = DTERMS - VODT * SATN
        RTERMS = RTERMS - VODT * SATO * CONCOIM(N,ICOMP)
      ELSE
        K = LAYNOD(N)
        LC=LAYCON(K)
        IF(LC.EQ.4.OR.LC.EQ.5)THEN !
C---------TWO STORAGES FOR UPSTREAM WEIGHTING FORMULATION
          DTERMS = DTERMS -
     *      ((PRSITYIM(N)*VOLU-SC2IM(N)*ALENG) + ALENG*SC2IM(N)*SATN
     *      +SATN*SC1IM(N)*(HNEWIM(N)-HOLDIM(N))) / DELT
          RTERMS = RTERMS -
     *      ((PRSITYIM(N)*VOLU-SC2IM(N)*ALENG) + ALENG*SC2IM(N)*SATO)
     *       * CONCOIM(N,ICOMP) / DELT
        ELSEIF(LC.EQ.0.OR.LC.EQ.2) THEN
C---------ONE CONFINED STORAGE CAPACITY FOR FLOW
            DTERMS = DTERMS - (PRSITYIM(N)*VOLU+SC1IM(N)*
     *       (HNEWIM(N)-HOLDIM(N)))/DELT
            RTERMS = RTERMS - (PRSITYIM(N)*VOLU)/DELT
     *       * CONCOIM(N,ICOMP)
        ELSEIF(LC.EQ.1) THEN ! IN BCF, IF LAYCON = 1, THE SPECIRFIC YIELD IS READ IN SC2 ?
C---------ONE UNCONFINED STORAGE CAPACITY FOR FLOW
          DTERMS = DTERMS -
     *   ((PRSITYIM(N)*VOLU-SC2IM(N)*ALENG)+SC2IM(N)*ALENG*SATN)/DELT
          RTERMS = RTERMS -
     *   ((PRSITYIM(N)*VOLU-SC2IM(N)*ALENG)+SC2IM(N)*ALENG*SATO)/DELT
     *     * CONCOIM(N,ICOMP)
        ELSEIF(LC.EQ.3) THEN
C---------TWO STORAGE CAPACITIES FOR FLOW
          TP=TOP(N)
          IF(HOLDIM(N).GT.TP.AND.HNEWIM(N).GT.TP) THEN
C-----------CELL IS CONFINED
            DTERMS = DTERMS - (PRSITYIM(N)*VOLU+SC1IM(N)*
     *       (HNEWIM(N)-HOLDIM(N)))/DELT
            RTERMS = RTERMS - (PRSITYIM(N)*VOLU)/DELT
     *       * CONCOIM(N,ICOMP)
          ELSEIF (HOLDIM(N).LT.TP.AND.HNEWIM(N).LT.TP) THEN
C-----------CELL IS UNCONFINED
            DTERMS = DTERMS -
     *   ((PRSITYIM(N)*VOLU-SC2IM(N)*ALENG)+SC2IM(N)*ALENG*SATN)/DELT
            RTERMS = RTERMS -
     *   ((PRSITYIM(N)*VOLU-SC2IM(N)*ALENG)+SC2IM(N)*ALENG*SATO)/DELT
     *     * CONCOIM(N,ICOMP)
          ELSE
C-----------CELL HAS CONVERTED THIS TIME STEP
            IF(HNEWIM(N).GT.HOLDIM(N)) THEN
C-------------WLEs RAISED - BECOMES CONFINED
              H1 = HNEWIM(N)
              H2 = TP
            ELSE
C-------------WLEs LOWERED - BECOMES UNCONFINED
              H1 = TP
              H2 = HOLDIM(N)
            ENDIF
            DTERMS = DTERMS - ((PRSITYIM(N)*VOLU-SC2IM(N)*ALENG) +
     *       SC2IM(N)*ALENG*SATN + SC1IM(N)*(H1-H2)) / DELT
            RTERMS = RTERMS - ((PRSITYIM(N)*VOLU-SC2IM(N)*ALENG) +
     *       SC2IM(N)*ALENG*SATO) / DELT
     *       * CONCOIM(N,ICOMP)
          ENDIF
        ENDIF
      ENDIF
C
C9------RETURN
      RETURN
      END
C----------------------------------------------------------------------
      SUBROUTINE GWT2DPT1DCYIM(N,ICOMP,DTERMD,RTERMD)
C     ******************************************************************
C     FORMULATE DECAY TERM FOR EACH NODE AND COMPONENT IN IMMOBILE DOMAIN
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,   ONLY:IOUT,NCOL,NROW,NLAY,IBOUND,TOP,NODLAY,ISYM,
     1 HNEW,AMAT,RHS,HOLD,ISSFLG,IA,JA,NJA,NODES,AREA,BOT,TOP,Sn,So
      USE GWTBCTMODULE, ONLY: ICT
      USE GWFDPFMODULE, ONLY: PHIF,SnIM,SoIM
      USE CLN1MODULE, ONLY: ACLNNDS
      USE GWFBASMODULE,ONLY:DELT
      USE GWTDPTMODULE, ONLY: CONCIM,CONCOIM,IADSORBIM,PRSITYIM,
     * ADSORBIM,FLICHIM,IZODIM,IFODIM,ZODRWIM,FODRWIM,
     1  ZODRSIM,FODRSIM, IAW_ADSORBIM
      DOUBLE PRECISION ADSTERM,CW,FL,CWO,DT,RT,ALENG,VOLU,
     1  X,Y,EPSS,EPS,CT,CEPS,QA,QEPS,DQ,OMP,DTERMD,RTERMD
C     ------------------------------------------------------------------
      DTERMD = 0.0
      RTERMD = 0.0
      IF(N.LE.NODES)THEN
        ALENG = TOP(N) - BOT(N)
      ELSE
        ALENG = ACLNNDS(N-NODES,5)
      ENDIF
      VOLU = AREA(N) * ALENG
C----------------------------------------------------------------------
      IF(ICT.EQ.0)THEN  !----------WATER PHASE CONCENTRATION FORMULATION
C------- ADD DECAY TERMS ON SOIL (NO ADSORPTION ON CLN)
        IF(N.LE.NODES)THEN
C---------ZERO ORDER DECAY ON SOIL - APPLY NEWTON EXPANSION OF SUPPLY-DEMAND FUNCTION
          IF(IZODIM.GE.2.AND.IADSORBIM.GT.0)THEN
            CT = - VOLU * ZODRSIM(N,ICOMP)
            EPS = 0.01
            CEPS = MAX(0.0,CONCIM(N,ICOMP))
            X = CEPS /EPS
            CALL SMOOTH(X,Y)
            QA =  CT * Y
C-----------CALCULATE DQ/DH
            EPSS = 0.001 * EPS
            CEPS = MAX(0.0,CONCIM(N,ICOMP)+EPSS)
            X = (CEPS)/EPS
            CALL SMOOTH(X,Y)
            QEPS = CT * Y
            DQ = (QEPS - QA) / EPSS
            DTERMD = DTERMD + DQ
            RTERMD = RTERMD - QA + DQ*CONCIM(N,ICOMP)
          ENDIF
C
C---------FIRST ORDER DECAY ON SOIL
          IF(IFODIM.GE.2.AND.IADSORBIM.GT.0)THEN
            CT = -ADSORBIM(N,ICOMP) * VOLU * FODRSIM(N,ICOMP)
            IF(IADSORBIM.EQ.1)THEN
C--------------FOR LINEAR ADSORPTION
              DTERMD = DTERMD + CT
            ELSE
C--------------FOR NON-LINEAR ADSORPTION FILL AS NEWTON
              ETA = FLICHIM(N,ICOMP)
              IF(IADSORBIM.EQ.2) THEN
                QA =  CT * CONCIM(N,ICOMP) ** ETA
              ELSEIF(IADSORBIM.EQ.3) THEN
                QA =  CT * CONCIM(N,ICOMP) / (1.0+ETA*CONCIM(N,ICOMP))  
              ENDIF    
C--------------CALCULATE DQ/DH
              EPSS = 0.00001
              CEPS = CONCIM(N,ICOMP)+EPSS
              IF(IADSORBIM.EQ.2) THEN
                QEPS = CT * CEPS ** ETA
              ELSEIF(IADSORBIM.EQ.3) THEN
                QEPS =  CT * CEPS / (1.0 + ETA*CEPS)  
              ENDIF   
              DQ = (QEPS - QA) / EPSS
              DTERMD = DTERMD + DQ
              RTERMD = RTERMD - QA + DQ*CONCIM(N,ICOMP)
            ENDIF
          ENDIF
C --------DECAY ON AIR-WATER INTERFACE ADSORBED SOLUTES          
          IF(IAW_ADSORBIM.NE.0) 
     *      CALL GWT2DPT1DCYAWIM(N,ICOMP,DTERMD,RTERMD,VOLU)
        ENDIF
C-----------------------------------------------------------------------
C---------ZERO ORDER DECAY IN WATER - APPLY NEWTON EXPANSION OF SUPPLY-DEMAND FUNCTION
        IF(IZODIM.EQ.1.OR.IZODIM.EQ.3)THEN
          CT = -SnIM(N)* VOLU * ZODRWIM(N,ICOMP)
          EPS = 0.01
          CEPS = MAX(0.0,CONCIM(N,ICOMP))
          X = CEPS /EPS
          CALL SMOOTH(X,Y)
          QA =  CT * Y
C-----------CALCULATE DQ/DH
          EPSS = 0.001 * EPS
          CEPS = MAX(0.0,CONCIM(N,ICOMP)+EPSS)
          X = (CEPS)/EPS
          CALL SMOOTH(X,Y)
          QEPS = CT * Y
          DQ = (QEPS - QA) / EPSS
          DTERMD = DTERMD + DQ
          RTERMD = RTERMD - QA + DQ*CONCIM(N,ICOMP)
        ENDIF
C
C---------FIRST ORDER DECAY IN WATER
        IF(IFODIM.EQ.1.OR.IFODIM.EQ.3)THEN
          CT =  -SnIM(N)* VOLU * FODRWIM(N,ICOMP)
          DTERMD = DTERMD + CT
        ENDIF
      ELSE       !-----------------------TOTAL CONCENTRATION FORMULATION
C
C-------NET DECAY TERM FOR TOTAL CONCENTRATION FORMULATION
CSP TO DO
      ENDIF
C
C9------RETURN
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE GWT2BCT1GENIM(N,ICOMP)
C     ******************************************************************
C     FORMULATE GENERATION TERM FOR EACH NODE AND COMPONENT
C     THIS SUBROUTINE IS SIMILAR TO GWT2BCT1DCY BUT LOOPING OVER
C     PARENT COMPONENTS. THEN, FORMULATE RHS AS = RTERM - DTERM * CONC
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,   ONLY:IOUT,NCOL,NROW,NLAY,IBOUND,TOP,NODLAY,ISYM,
     1 HNEW,AMAT,RHS,HOLD,ISSFLG,IA,JA,NJA,NODES,AREA,BOT,TOP,Sn,So
      USE CLN1MODULE, ONLY: ACLNNDS
      USE GWFBASMODULE,ONLY:DELT
      USE GWTBCTMODULE, ONLY: MCOMP,CONC,CONCO,ICT,IADSORB,PRSITY,
     * ADSORB,FLICH,IZOD,IFOD,ZODRW,FODRW,
     1  ZODRS,FODRS,ICHAIN,NPARENT,JPARENT,STOTIO,ISPRCT
      USE GWTDPTMODULE, ONLY: CONCIM,CONCOIM,IADSORBIM,PRSITYIM,
     * ADSORBIM,FLICHIM,IZODIM,IFODIM,ZODRWIM,FODRWIM,
     1  ZODRSIM,FODRSIM, SPTLRCTIM, IAW_ADSORBIM     
      USE GWFDPFMODULE, ONLY: PHIF,SnIM,SoIM      
      DOUBLE PRECISION ADSTERM,CW,FL,CWO,DT,RT,ALENG,VOLU,
     1  X,Y,EPSS,EPS,CT,CEPS,QA,QEPS,DQ,OMP,DTERMD,RTERMD
C     ------------------------------------------------------------------
C1------INITIALIZE TERMS
      IF(N.LE.NODES)THEN
        ALENG = TOP(N) - BOT(N)
      ELSE
        ALENG = ACLNNDS(N-NODES,4)
      ENDIF
      VOLU = AREA(N) * ALENG
C----------------------------------------------------------------------
C ------ DO FOR ALL PARENTS OF COMPONENT ICOMP IF IT HAS PARENTS
      IF(NPARENT(ICOMP).EQ.0) RETURN
      DO NPAREN = 1, NPARENT(ICOMP)
        DTERMD = 0.0
        RTERMD = 0.0
        IPAREN = JPARENT(ICOMP,NPAREN)
        IF(ISPRCT.EQ.0)THEN 
          STOIT  = STOTIO(ICOMP,NPAREN)
        ELSE
          STOIT = SPTLRCTIM(N,ICOMP,NPAREN)
        ENDIF        
C----------------------------------------------------------------------
      IF(ICT.EQ.0)THEN  !----------WATER PHASE CONCENTRATION FORMULATION
C2-------GENERATION TERMS ON SOIL (NO ADSORPTION ON CLN)
        IF(N.LE.NODES)THEN
C3---------ZERO ORDER DECAY ON SOIL - APPLY NEWTON EXPANSION OF SUPPLY-DEMAND FUNCTION
          IF(IZODIM.GE.2.AND.IADSORBIM.GT.0)THEN
            CT = VOLU * ZODRSIM(N,IPAREN) * STOIT
            EPS = 0.01
            CEPS = MAX(0.0,CONCIM(N,IPAREN))
            X = CEPS /EPS
            CALL SMOOTH(X,Y)
            QA =  CT * Y
C3A-----------CALCULATE DQ/DH
            EPSS = 0.001 * EPS
            CEPS = MAX(0.0,CONCIM(N,IPAREN)+EPSS)
            X = (CEPS)/EPS
            CALL SMOOTH(X,Y)
            QEPS = CT * Y
            DQ = (QEPS - QA) / EPSS
            DTERMD = DTERMD + DQ
            RTERMD = RTERMD - QA + DQ*CONCIM(N,IPAREN)
          ENDIF
C
C4---------FIRST ORDER DECAY ON SOIL
          IF(IFODIM.GE.2.AND.IADSORBIM.GT.0)THEN
            CT = ADSORBIM(N,IPAREN) * VOLU * FODRSIM(N,IPAREN) * STOIT
            IF(IADSORBIM.EQ.1)THEN
C5--------------FOR LINEAR ADSORPTION
              DTERMD = DTERMD + CT
            ELSE
C--------------FOR NON-LINEAR ADSORPTION FILL AS NEWTON
              ETA = FLICHIM(N,IPAREN)
              IF(IADSORBIM.EQ.2) THEN
                QA =  CT * CONCIM(N,IPAREN) ** ETA
              ELSEIF(IADSORBIM.EQ.3) THEN
                QA =  CT * CONCIM(N,IPAREN) / (1.0+ETA*CONCIM(N,IPAREN))
              ENDIF    
C--------------CALCULATE DQ/DH
              EPSS = 0.00001
              CEPS = CONCIM(N,IPAREN)+EPSS
              IF(IADSORBIM.EQ.2) THEN
                QEPS = CT * CEPS ** ETA
              ELSEIF(IADSORBIM.EQ.3) THEN
                QEPS =  CT * CEPS / (1.0 + ETA*CEPS)  
              ENDIF   
              DQ = (QEPS - QA) / EPSS
              DTERMD = DTERMD + DQ
              RTERMD = RTERMD - QA + DQ*CONCIM(N,IPAREN)              
            ENDIF
          ENDIF
C --------GENERATION ON AIR-WATER INTERFACE ADSORBED SOLUTES            
          IF(IAW_ADSORBIM.NE.0)          
     *      CALL GWT2DPT1GENAWIM(N,ICOMP,DTERMD,RTERMD,VOLU,NPAREN) 
        ENDIF
C-----------------------------------------------------------------------
C7---------ZERO ORDER DECAY IN WATER - APPLY NEWTON EXPANSION OF SUPPLY-DEMAND FUNCTION
        IF(IZODIM.EQ.1.OR.IZODIM.EQ.3)THEN
          CT = SnIM(N)* VOLU * ZODRWIM(N,IPAREN)*STOIT
          EPS = 0.01
          CEPS = MAX(0.0,CONCIM(N,IPAREN))
          X = CEPS /EPS
          CALL SMOOTH(X,Y)
          QA =  CT * Y
C7A-----------CALCULATE DQ/DH
          EPSS = 0.001 * EPS
          CEPS = MAX(0.0,CONCIM(N,IPAREN)+EPSS)
          X = (CEPS)/EPS
          CALL SMOOTH(X,Y)
          QEPS = CT * Y
          DQ = (QEPS - QA) / EPSS
          DTERMD = DTERMD + DQ
          RTERMD = RTERMD - QA + DQ*CONCIM(N,IPAREN)
        ENDIF
C
C8---------FIRST ORDER DECAY IN WATER
        IF(IFODIM.EQ.1.OR.IFODIM.EQ.3)THEN
          CT = SnIM(N)* VOLU * FODRWIM(N,IPAREN) * STOIT
          DTERMD = DTERMD + CT
        ENDIF
      ELSE       !-----------------------TOTAL CONCENTRATION FORMULATION
C
C-------NET DECAY TERM FOR TOTAL CONCENTRATION FORMULATION
CSP TO DO
      ENDIF
      RHS(N) = RHS(N) + RTERMD - DTERMD * CONCIM(N,IPAREN)
      ENDDO   !---------------------------END DO LOOP OF PARENTS
C
C9------RETURN
      RETURN
      END
C-------------------------------------------------------------------------
      SUBROUTINE GWT2STOIM1BDW(KSTP,KPER,ICOMP)
C     ******************************************************************
C     CALCULATE MASS BUDGET TERMS FOR IMMOBILE DOMAIN FOR ALL TRANSPORT CELLS
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,   ONLY:IOUT,NCOL,NROW,NLAY,IBOUND,BUFF,NODES,IUNSTR,
     1  AMAT,IA,JA,TOP,BOT,AREA,Sn,So,ISSFLG
      USE CLN1MODULE, ONLY: ACLNNDS,NCLNNDS
      USE GWFBASMODULE,ONLY:MSUM,ISPCFL,IAUXSV,DELT,PERTIM,TOTIM
      USE GWTBCTMODULE, ONLY: MSUMT,VBVLT,VBNMT,ICT
      USE GWTDPTMODULE, ONLY: ICBUNDIM,CONCIM,CONCOIM,IDPTCB,
     1 IADSORBIM,ADSORBIM,FLICHIM,PRSITYIM,CONCOIM
      USE GWFDPFMODULE, ONLY: PHIF,SnIM,SoIM
C
      CHARACTER*16 TEXT
      DOUBLE PRECISION RATIN,RATOUT,QQ,VODT,ADSTERM,FL,CW,CWO,ALENG,
     *  DTERMS,RTERMS,VOLU
      DATA TEXT /' IM MASS STORAGE'/
C     ------------------------------------------------------------------
C
C1------CLEAR RATIN AND RATOUT ACCUMULATORS, AND SET CELL-BY-CELL
C1------BUDGET FLAG.
      ZERO=0.
      RATIN=ZERO
      RATOUT=ZERO
      IBD=0
      ISS = ISSFLG(KPER)
      IF(IDPTCB.LT.0 .AND. ISPCFL.NE.0) IBD=-1
      IF(IDPTCB.GT.0) IBD=ISPCFL
      IBDLBL=0
C
C3------CLEAR THE BUFFER.
      DO 50 N=1,NODES
      BUFF(N)=ZERO
50    CONTINUE
C
C5------LOOP THROUGH EACH NODE AND CALCULATE STORAGE
      DO 100 N=1,NODES
C
C5B-----IF THE CELL IS NOT PCB OR WRONG COMPONENT SPECIES, IGNORE IT.
      IF(ICBUNDIM(N).EQ.0)GO TO 99
C
      ALENG = TOP(N) - BOT(N)
      VOLU = AREA(N) * ALENG
      VODT = VOLU / DELT
      QQ = 0.0
      IF(ICT.EQ.0)THEN
C-------STORAGE TERM IN WATER
CSP        IF(N.LE.NODES) VODT = VODT * PRSITYIM(N)
CSP        QQ = QQ + VODT * SnIM(N) * CONCIM(N,ICOMP)
CSP     *          - VODT * SoIM(N) * CONCOIM(N,ICOMP)
        DTERMS = 0.0
        RTERMS = 0.0
        CALL GWT2BCT1STOIMW(N,ICOMP,DTERMS,RTERMS,VODT,VOLU,ALENG,ISS)
        QQ = QQ - DTERMS * CONCIM(N,ICOMP) + RTERMS
      ELSE
        QQ = QQ + VODT * CONCIM(N,ICOMP)
     *          - VODT * CONCOIM(N,ICOMP)
      ENDIF
      QQ = - QQ  ! STORAGE TERM NEGATIVE IS INFLOW AS PER MODFLOW CONVENTION
      Q = QQ
C
C5D-----PRINT FLOW RATE IF REQUESTED.
      IF(IBD.LT.0) THEN
         IF(IBDLBL.EQ.0) WRITE(IOUT,61) TEXT,KPER,KSTP
   61    FORMAT(1X,/1X,A,'   PERIOD ',I4,'   STEP ',I3)
        IF(IUNSTR.EQ.0.AND.N.LE.NODES)THEN
          IL = (N-1) / (NCOL*NROW) + 1
          IJ = N - (IL-1)*NCOL*NROW
          IR = (IJ-1)/NCOL + 1
          IC = IJ - (IR-1)*NCOL
           WRITE(IOUT,62) L,IL,IR,IC,Q
   62    FORMAT(1X,'CBC  ',I6,'   LAYER ',I3,'   ROW ',I5,'   COL ',I5,
     1       '   FLUX ',1PG15.6)
        ELSE
           WRITE(IOUT,63) L,N,Q
   63    FORMAT(1X,'CBC  ',I6,'    NODE ',I8,'   FLUX ',1PG15.6)
        ENDIF
        IBDLBL=1
      END IF
C
C5E-----ADD FLOW RATE TO BUFFER.
      BUFF(N)=BUFF(N)+Q
C
C5F-----SEE IF FLUX IS POSITIVE OR NEGATIVE.
      IF(QQ.GE.ZERO) THEN
C
C5G-----POSITIVE FLOW RATE. ADD IT TO RATIN
        RATIN=RATIN+QQ
      ELSE
C
C5H-----NEGATIVE FLOW RATE. ADD IT TO RATOUT
        RATOUT=RATOUT-QQ
      END IF
   99 CONTINUE
C
100   CONTINUE
C
C6------IF CELL-BY-CELL FLOWS WILL BE SAVED AS A 3-D ARRAY,
C6------CALL UBUDSV TO SAVE THEM.
      IF(IBD.EQ.1)THEN
        IF(IUNSTR.EQ.0)THEN
          CALL UBUDSV(KSTP,KPER,TEXT,IDPTCB,BUFF,NCOL,NROW,
     1                          NLAY,IOUT)
        ELSE
          CALL UBUDSVU(KSTP,KPER,TEXT,IDPTCB,BUFF,NODES,
     1                          IOUT,PERTIM,TOTIM)
        ENDIF
      ENDIF
C
C7------MOVE RATES, VOLUMES & LABELS INTO ARRAYS FOR PRINTING.
  200 RIN=RATIN
      ROUT=RATOUT
      VBVLT(3,MSUMT,ICOMP)=RIN
      VBVLT(4,MSUMT,ICOMP)=ROUT
      VBVLT(1,MSUMT,ICOMP)=VBVLT(1,MSUMT,ICOMP)+RATIN*DELT
      VBVLT(2,MSUMT,ICOMP)=VBVLT(2,MSUMT,ICOMP)+RATOUT*DELT
      VBNMT(MSUMT,ICOMP)=TEXT
C
C8------INCREMENT BUDGET TERM COUNTER(MSUM).
      MSUMT=MSUMT+1
C
C9------RETURN
      RETURN
      END
C-------------------------------------------------------------------------
      SUBROUTINE GWT2STOIM1BDS(KSTP,KPER,ICOMP)
C     ******************************************************************
C     CALCULATE MASS BUDGET TERMS FOR IMMOBILE DOMAIN FOR ALL TRANSPORT CELLS
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,   ONLY:IOUT,NCOL,NROW,NLAY,IBOUND,BUFF,NODES,IUNSTR,
     1  AMAT,IA,JA,TOP,BOT,AREA,Sn,So,ISSFLG
      USE CLN1MODULE, ONLY: ACLNNDS,NCLNNDS
      USE GWFBASMODULE,ONLY:MSUM,ISPCFL,IAUXSV,DELT,PERTIM,TOTIM
      USE GWTBCTMODULE, ONLY: MSUMT,VBVLT,VBNMT,ICT
      USE GWTDPTMODULE, ONLY: ICBUNDIM,CONCIM,CONCOIM,IDPTCB,
     1 IADSORBIM,ADSORBIM,FLICHIM,PRSITYIM,CONCOIM
      USE GWFDPFMODULE, ONLY: PHIF,SnIM,SoIM
C
      CHARACTER*16 TEXT
      DOUBLE PRECISION RATIN,RATOUT,QQ,VODT,ADSTERM,FL,CW,CWO,ALENG,
     *  DTERMS,RTERMS,VOLU
      DATA TEXT /'IM ADSORBED MASS'/ 
C     ------------------------------------------------------------------
C
C1------CLEAR RATIN AND RATOUT ACCUMULATORS, AND SET CELL-BY-CELL
C1------BUDGET FLAG.
      ZERO=0.
      RATIN=ZERO
      RATOUT=ZERO
      IBD=0
      ISS = ISSFLG(KPER)
      IF(IDPTCB.LT.0 .AND. ISPCFL.NE.0) IBD=-1
      IF(IDPTCB.GT.0) IBD=ISPCFL
      IBDLBL=0
C
C3------CLEAR THE BUFFER.
      DO 50 N=1,NODES
      BUFF(N)=ZERO
50    CONTINUE
C
C5------LOOP THROUGH EACH NODE AND CALCULATE STORAGE
      DO 100 N=1,NODES
C
C5B-----IF THE CELL IS NOT PCB OR WRONG COMPONENT SPECIES, IGNORE IT.
      IF(ICBUNDIM(N).EQ.0)GO TO 99
C
      ALENG = TOP(N) - BOT(N)
      VOLU = AREA(N) * ALENG
      VODT = VOLU / DELT
      QQ = 0.0
      IF(ICT.EQ.0)THEN
C-------STORAGE TERM ON SOIL
        IF(IADSORBIM.EQ.1)THEN
          ADSTERM = ADSORBIM(N,ICOMP) * VODT
          QQ = ADSTERM * (CONCIM(N,ICOMP) - CONCOIM(N,ICOMP))
        ELSEIF(IADSORBIM.EQ.2)THEN
          ADSTERM = ADSORBIM(N,ICOMP) * VODT
          FL = FLICHIM(N,ICOMP)
          CW = 0.0
          CWO = 0.0
          IF(CONCIM(N,ICOMP).GT.0.0) CW = CONCIM(N,ICOMP)
          IF(CONCOIM(N,ICOMP).GT.0.0) CWO = CONCOIM(N,ICOMP)
          QQ = ADSTERM * (CW**FL - CWO**FL)
        ELSEIF(IADSORBIM.EQ.3)THEN
          ADSTERM = ADSORBIM(N,ICOMP) * VODT
          FL = FLICHIM(N,ICOMP)
          CW = 0.0
          CWO = 0.0
          IF(CONCIM(N,ICOMP).GT.0.0) CW = CONCIM(N,ICOMP)
          IF(CONCOIM(N,ICOMP).GT.0.0) CWO = CONCOIM(N,ICOMP)
          QQ = ADSTERM * (CW/(1.0+FL*CW) - CWO/(1.0+FL*CWO))          
        ENDIF
      ELSE
        QQ = QQ + VODT * CONCIM(N,ICOMP)
     *          - VODT * CONCOIM(N,ICOMP)
      ENDIF
      QQ = - QQ  ! STORAGE TERM NEGATIVE IS INFLOW AS PER MODFLOW CONVENTION
      Q = QQ
C
C5D-----PRINT FLOW RATE IF REQUESTED.
      IF(IBD.LT.0) THEN
         IF(IBDLBL.EQ.0) WRITE(IOUT,61) TEXT,KPER,KSTP
   61    FORMAT(1X,/1X,A,'   PERIOD ',I4,'   STEP ',I3)
        IF(IUNSTR.EQ.0.AND.N.LE.NODES)THEN
          IL = (N-1) / (NCOL*NROW) + 1
          IJ = N - (IL-1)*NCOL*NROW
          IR = (IJ-1)/NCOL + 1
          IC = IJ - (IR-1)*NCOL
           WRITE(IOUT,62) L,IL,IR,IC,Q
   62    FORMAT(1X,'CBC  ',I6,'   LAYER ',I3,'   ROW ',I5,'   COL ',I5,
     1       '   FLUX ',1PG15.6)
        ELSE
           WRITE(IOUT,63) L,N,Q
   63    FORMAT(1X,'CBC  ',I6,'    NODE ',I8,'   FLUX ',1PG15.6)
        ENDIF
        IBDLBL=1
      END IF
C
C5E-----ADD FLOW RATE TO BUFFER.
      BUFF(N)=BUFF(N)+Q
C
C5F-----SEE IF FLUX IS POSITIVE OR NEGATIVE.
      IF(QQ.GE.ZERO) THEN
C
C5G-----POSITIVE FLOW RATE. ADD IT TO RATIN
        RATIN=RATIN+QQ
      ELSE
C
C5H-----NEGATIVE FLOW RATE. ADD IT TO RATOUT
        RATOUT=RATOUT-QQ
      END IF
   99 CONTINUE
C
100   CONTINUE
C
C6------IF CELL-BY-CELL FLOWS WILL BE SAVED AS A 3-D ARRAY,
C6------CALL UBUDSV TO SAVE THEM.
      IF(IBD.EQ.1)THEN
        IF(IUNSTR.EQ.0)THEN
          CALL UBUDSV(KSTP,KPER,TEXT,IDPTCB,BUFF,NCOL,NROW,
     1                          NLAY,IOUT)
        ELSE
          CALL UBUDSVU(KSTP,KPER,TEXT,IDPTCB,BUFF,NODES,
     1                          IOUT,PERTIM,TOTIM)
        ENDIF
      ENDIF
C
C7------MOVE RATES, VOLUMES & LABELS INTO ARRAYS FOR PRINTING.
  200 RIN=RATIN
      ROUT=RATOUT
      VBVLT(3,MSUMT,ICOMP)=RIN
      VBVLT(4,MSUMT,ICOMP)=ROUT
      VBVLT(1,MSUMT,ICOMP)=VBVLT(1,MSUMT,ICOMP)+RATIN*DELT
      VBVLT(2,MSUMT,ICOMP)=VBVLT(2,MSUMT,ICOMP)+RATOUT*DELT
      VBNMT(MSUMT,ICOMP)=TEXT
C
C8------INCREMENT BUDGET TERM COUNTER(MSUM).
      MSUMT=MSUMT+1
C
C9------RETURN
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE GWT2DCYIM1BD(KSTP,KPER,ICOMP)
C     ******************************************************************
C     CALCULATE MASS DECAY TERMS FOR ALL TRANSPORT CELLS
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,   ONLY:IOUT,NCOL,NROW,NLAY,IBOUND,BUFF,NODES,IUNSTR,
     1  AMAT,IA,JA,TOP,BOT,AREA,Sn,So
      USE CLN1MODULE, ONLY: ACLNNDS,NCLNNDS
      USE GWFBASMODULE,ONLY:MSUM,ISPCFL,IAUXSV,DELT,PERTIM,TOTIM,DELT
      USE GWTBCTMODULE, ONLY: MSUMT,VBVLT,VBNMT,ICT
      USE GWTDPTMODULE, ONLY: ICBUNDIM,CONCIM,IDPTCB,
     1 IADSORBIM,ADSORBIM,FLICHIM,PRSITYIM,CONCOIM,IZODIM,IFODIM,
     1 ZODRWIM,FODRWIM,ZODRSIM,FODRSIM, IAW_ADSORBIM
      USE GWFDPFMODULE, ONLY: PHIF,SnIM
C
      CHARACTER*16 TEXT(2)
      DOUBLE PRECISION RATIN,RATOUT,QQ,VOLU,ADSTERM,FL,CW,CWO,ALENG,X,Y,
     1  CEPS,EPS,CT
      DATA TEXT(1) /'   IM MASS DECAY'/
      DATA TEXT(2) /'IMCLN MASS DECAY'/
C     ------------------------------------------------------------------
C1------RETURN IF NO DECAY IN SIMULATION
      IF(IZODIM.EQ.0.AND.IFODIM.EQ.0) RETURN
C2------CLEAR RATIN AND RATOUT ACCUMULATORS, AND SET CELL-BY-CELL
C2------BUDGET FLAG.
      ZERO=0.
      RATIN=ZERO
      RATOUT=ZERO
      IBD=0
      IF(IDPTCB.LT.0 .AND. ISPCFL.NE.0) IBD=-1
      IF(IDPTCB.GT.0) IBD=ISPCFL
      IBDLBL=0
C
C3------CLEAR THE BUFFER.
      DO 50 N=1,NODES
      BUFF(N)=ZERO
50    CONTINUE
C
C5------LOOP THROUGH EACH NODE AND CALCULATE STORAGE
      DO 100 N=1,NODES
C
C5B-----IF THE CELL IS NOT PCB OR WRONG COMPONENT SPECIES, IGNORE IT.
      IF(ICBUNDIM(N).EQ.0)GO TO 99
C
      ALENG = TOP(N) - BOT(N)
      VOLU = AREA(N) * ALENG
      QQ = 0.0
C-----------------------------------------------------------------------------
      IF(ICT.EQ.0)THEN  !----------WATER PHASE CONCENTRATION FORMULATION
C-----------------------------------------------------------------------------
C-------DECAY TERM ON SOIL (NO ADSORPTION ON CLN)
        IF(N.LE.NODES)THEN
C
C---------ZERO ORDER DECAY ON SOIL
            IF(IZODIM.GE.2.AND.IADSORBIM.GT.0)THEN
            CT = -(1.0-PRSITYIM(N)) * VOLU * ZODRSIM(N,ICOMP)
            EPS = 0.01
            CEPS = MAX(0.0,CONCIM(N,ICOMP))
            X = CEPS /EPS
            CALL SMOOTH(X,Y)
            QQ =  CT * Y
          ENDIF
C
C---------FIRST ORDER DECAY ON SOIL
          IF(IFODIM.GE.2.AND.IADSORBIM.GT.0)THEN
            CT = -ADSORBIM(N,ICOMP) * VOLU * FODRSIM(N,ICOMP) ! * (1.0-PRSITYIM(N)) ???
            IF(IADSORBIM.EQ.1)THEN
C--------------FOR LINEAR ADSORPTION
              QQ = QQ + CT * CONCIM(N,ICOMP)
            ELSEIF(IADSORBIM.EQ.2)THEN
C--------------FOR NON-LINEAR FREUNDLICH ADSORPTION 
              ETA = FLICHIM(N,ICOMP)
              QQ =  CT * CONCIM(N,ICOMP) ** ETA
            ELSEIF(IADSORBIM.EQ.3)THEN
C--------------FOR NON-LINEAR LANGMUIR ADSORPTION 
              ETA = FLICHIM(N,ICOMP)
              QQ =  CT * CONCIM(N,ICOMP) / (1.0 + ETA*CONCIM(N,ICOMP)) 
            ENDIF
          ENDIF
C --------DECAY ON AIR-WATER INTERFACE ADSORBED SOLUTES              
          IF(IAW_ADSORBIM.NE.0)
     *     CALL GWT2DPT1DCYAWIM1BD(N,ICOMP,QQ,VOLU)
        ENDIF
C-----------------------------------------------------------------------------          
C-------DECAY TERM IN WATER
C-----------------------------------------------------------------------------
C-------ZERO ORDER DECAY IN WATER 
        IF(IZODIM.EQ.1.OR.IZODIM.EQ.3)THEN
          CT = -SnIM(N)* VOLU * ZODRWIM(N,ICOMP)
          EPS = 0.01
          CEPS = MAX(0.0,CONCIM(N,ICOMP))
          X = CEPS /EPS
          CALL SMOOTH(X,Y)
          QQ =  QQ + CT * Y
        ENDIF
C
C-------FIRST ORDER DECAY IN WATER
        IF(IFODIM.EQ.1.OR.IFODIM.EQ.3)THEN
          CT =  -SnIM(N)* VOLU * FODRWIM(N,ICOMP)
          QQ = QQ + CT * CONCIM(N,ICOMP)
        ENDIF
      ELSE
CSP FINISH TOTAL CONCENTRATION FORMULATION
      ENDIF
      Q = QQ
C
C5D-----PRINT FLOW RATE IF REQUESTED.
      IF(IBD.LT.0) THEN
         IF(IBDLBL.EQ.0) WRITE(IOUT,61) TEXT(1),KPER,KSTP
   61    FORMAT(1X,/1X,A,'   PERIOD ',I4,'   STEP ',I3)
        IF(IUNSTR.EQ.0.AND.N.LE.NODES)THEN
          IL = (N-1) / (NCOL*NROW) + 1
          IJ = N - (IL-1)*NCOL*NROW
          IR = (IJ-1)/NCOL + 1
          IC = IJ - (IR-1)*NCOL
           WRITE(IOUT,62) L,IL,IR,IC,Q
   62    FORMAT(1X,'CBC  ',I6,'   LAYER ',I3,'   ROW ',I5,'   COL ',I5,
     1       '  DECAY ',1PG15.6)
        ELSE
           WRITE(IOUT,63) L,N,Q
   63    FORMAT(1X,'CBC  ',I6,'    NODE ',I8,'  DECAY ',1PG15.6)
        ENDIF
        IBDLBL=1
      END IF
C
C5E-----ADD FLOW RATE TO BUFFER.
      BUFF(N)=BUFF(N)+Q
C
C5F-----SEE IF FLUX IS POSITIVE OR NEGATIVE.
      IF(QQ.GE.ZERO) THEN
C
C5G-----POSITIVE FLOW RATE. ADD IT TO RATIN
        RATIN=RATIN+QQ
      ELSE
C
C5H-----NEGATIVE FLOW RATE. ADD IT TO RATOUT
        RATOUT=RATOUT-QQ
      END IF
   99 CONTINUE
C
100   CONTINUE
C
C6------IF CELL-BY-CELL FLOWS WILL BE SAVED AS A 3-D ARRAY,
C6------CALL UBUDSV TO SAVE THEM.
      IF(IBD.EQ.1)THEN
        IF(IUNSTR.EQ.0)THEN
          CALL UBUDSV(KSTP,KPER,TEXT(1),IDPTCB,BUFF,NCOL,NROW,
     1                          NLAY,IOUT)
        ELSE
          CALL UBUDSVU(KSTP,KPER,TEXT(1),IDPTCB,BUFF,NODES,
     1                          IOUT,PERTIM,TOTIM)
        ENDIF
      ENDIF
C
C7------MOVE RATES, VOLUMES & LABELS INTO ARRAYS FOR PRINTING.
  200 RIN=RATIN
      ROUT=RATOUT
      VBVLT(3,MSUMT,ICOMP)=RIN
      VBVLT(4,MSUMT,ICOMP)=ROUT
      VBVLT(1,MSUMT,ICOMP)=VBVLT(1,MSUMT,ICOMP)+RATIN*DELT
      VBVLT(2,MSUMT,ICOMP)=VBVLT(2,MSUMT,ICOMP)+RATOUT*DELT
      VBNMT(MSUMT,ICOMP)=TEXT(1)
C
C8------INCREMENT BUDGET TERM COUNTER(MSUM).
      MSUMT=MSUMT+1
C
C9------RETURN
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE GWT2GENIM1BD(KSTP,KPER,ICOMP)
C     ******************************************************************
C     CALCULATE MASS GENERATION TERMS FOR ALL TRANSPORT CELLS
C     THIS SUBROUTINE IS SIMILAR TO GWT2GEN1BD BUT LOOPING OVER
C     PARENT COMPONENTS.
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,   ONLY:IOUT,NCOL,NROW,NLAY,IBOUND,BUFF,NODES,IUNSTR,
     1  AMAT,IA,JA,TOP,BOT,AREA,Sn,So,NEQS,INCLN
      USE CLN1MODULE, ONLY: ACLNNDS,NCLNNDS
      USE GWFBASMODULE,ONLY:MSUM,ICBCFL,IAUXSV,DELT,PERTIM,TOTIM,DELT
      USE GWTBCTMODULE, ONLY: ICBUND,CONC,MSUMT,VBVLT,VBNMT,IBCTCB,
     1 IADSORB,ADSORB,FLICH,PRSITY,CONCO,ICT,IZOD,IFOD,ZODRW,FODRW,
     1  ZODRS,FODRS,ICHAIN,NPARENT,JPARENT,STOTIO,ISPRCT
      USE GWTDPTMODULE, ONLY: ICBUNDIM,CONCIM,IDPTCB,
     1 IADSORBIM,ADSORBIM,FLICHIM,PRSITYIM,CONCOIM,IZODIM,IFODIM,
     1 ZODRWIM,FODRWIM,ZODRSIM,FODRSIM,SPTLRCTIM, IAW_ADSORBIM
      USE GWFDPFMODULE, ONLY: PHIF,SnIM      
C
      CHARACTER*16 TEXT(2)
      DOUBLE PRECISION RATIN,RATOUT,QQ,VOLU,ADSTERM,FL,CW,CWO,ALENG,X,Y,
     1  CEPS,EPS,CT
      DATA TEXT(1) /'     IM MASS GEN'/
      DATA TEXT(2) /' IM CLN MASS GEN'/
C     ------------------------------------------------------------------
C1------RETURN IF NO DECAY IN SIMULATION
      IF(IZODIM.EQ.0.AND.IFODIM.EQ.0) RETURN
C2------CLEAR RATIN AND RATOUT ACCUMULATORS, AND SET CELL-BY-CELL
C2------BUDGET FLAG.
      ZERO=0.
      RATIN=ZERO
      RATOUT=ZERO
      IBD=0
      IF(IBCTCB.LT.0 .AND. ICBCFL.NE.0) IBD=-1
      IF(IBCTCB.GT.0) IBD=ICBCFL
      IBDLBL=0
C
C3------CLEAR THE BUFFER.
      DO 50 N=1,NEQS
      BUFF(N)=ZERO
50    CONTINUE
C
C4------LOOP THROUGH EACH NODE AND CALCULATE STORAGE
      DO 100 N=1,NODES
C
C5-----IF THE CELL IS NOT PCB OR WRONG COMPONENT SPECIES, IGNORE IT.
      IF(ICBUNDIM(N).EQ.0)GO TO 99
C
      IF(N.LE.NODES)THEN
        ALENG = TOP(N) - BOT(N)
      ELSE
        ALENG = ACLNNDS(N-NODES,4)
      ENDIF
      VOLU = AREA(N) * ALENG
      QQ = 0.0
C----------------------------------------------------------------------
C ------ DO FOR ALL PARENTS OF COMPONENT ICOMP
      IF(NPARENT(ICOMP).EQ.0) GO TO 111
      DO NPAREN = 1, NPARENT(ICOMP)
        DTERMD = 0.0
        RTERMD = 0.0
        IPAREN = JPARENT(ICOMP,NPAREN)
        IF(ISPRCT.EQ.0)THEN 
          STOIT  = STOTIO(ICOMP,NPAREN)
        ELSE
          STOIT = SPTLRCTIM(N,ICOMP,NPAREN)
        ENDIF
C-----------------------------------------------------------------------------
      IF(ICT.EQ.0)THEN  !----------WATER PHASE CONCENTRATION FORMULATION
C-----------------------------------------------------------------------------
C6-------GENERATION TERM ON SOIL (NO ADSORPTION ON CLN)
        IF(N.LE.NODES)THEN
C
C7---------ZERO ORDER GENERATION ON SOIL
          IF(IZODIM.GE.2.AND.IADSORBIM.GT.0)THEN
            CT = VOLU * ZODRSIM(N,IPAREN) * STOIT
            EPS = 0.01
            CEPS = MAX(0.0,CONCIM(N,IPAREN))
            X = CEPS /EPS
            CALL SMOOTH(X,Y)
            QQ =  QQ + CT * Y
          ENDIF
C
C8---------FIRST ORDER GENERATION ON SOIL
          IF(IFODIM.GE.2.AND.IADSORBIM.GT.0)THEN
            CT = ADSORBIM(N,IPAREN) * VOLU * FODRSIM(N,IPAREN) * STOIT
            IF(IADSORBIM.EQ.1)THEN
C9--------------FOR LINEAR ADSORPTION
              QQ = QQ + CT * CONCIM(N,IPAREN)
            ELSEIF(IADSORBIM.EQ.2) THEN 
C10--------------FOR NON-LINEAR FREUNDLICH ADSORPTION 
              ETA = FLICHIM(N,IPAREN)
              QQ =  QQ + CT * CONCIM(N,IPAREN) ** ETA
C10--------------FOR NON-LINEAR LANGMUIR ADSORPTION 
              ETA = FLICHIM(N,IPAREN)
              QQ =  QQ + CT*CONCIM(N,IPAREN)/ (1.0+ETA*CONCIM(N,IPAREN))
            ENDIF
          ENDIF
C --------GENERATION ON AIR-WATER INTERFACE ADSORBED SOLUTES              
          IF(IAW_ADSORBIM.NE.0)
     *     CALL GWT2DPT1GENAWIM1BD(N,ICOMP,QQ,VOLU,NPAREN)             
               
        ENDIF
C-----------------------------------------------------------------------------
C11-------GENERATION TERM IN WATER
C-----------------------------------------------------------------------------
C12-------ZERO ORDER GENERATION IN WATER
        IF(IZODIM.EQ.1.OR.IZODIM.EQ.3)THEN
          CT = SnIM(N)* VOLU * ZODRWIM(N,IPAREN) * STOIT
          EPS = 0.01
          CEPS = MAX(0.0,CONCIM(N,IPAREN))
          X = CEPS /EPS
          CALL SMOOTH(X,Y)
          QQ =  QQ + CT * Y
        ENDIF
C
C13-------FIRST ORDER GENERATION IN WATER
        IF(IFODIM.EQ.1.OR.IFODIM.EQ.3)THEN
          CT =  SnIM(N)* VOLU * FODRWIM(N,IPAREN) * STOIT
          QQ = QQ + CT * CONCIM(N,IPAREN)
        ENDIF
      ELSE
CSP FINISH TOTAL CONCENTRATION FORMULATION
      ENDIF
      ENDDO    !---------------------------END DO LOOP OF PARENTS
111   CONTINUE !---------------------------SKIP HERE IF NO PARENTS
      Q = QQ
C
C14-----PRINT FLOW RATE IF REQUESTED.
      IF(IBD.LT.0) THEN
         IF(IBDLBL.EQ.0) WRITE(IOUT,61) TEXT(1),KPER,KSTP
   61    FORMAT(1X,/1X,A,'   PERIOD ',I4,'   STEP ',I3)
        IF(IUNSTR.EQ.0.AND.N.LE.NODES)THEN
          IL = (N-1) / (NCOL*NROW) + 1
          IJ = N - (IL-1)*NCOL*NROW
          IR = (IJ-1)/NCOL + 1
          IC = IJ - (IR-1)*NCOL
           WRITE(IOUT,62) L,IL,IR,IC,Q
   62    FORMAT(1X,'CBC  ',I6,'   LAYER ',I3,'   ROW ',I5,'   COL ',I5,
     1       'GENERATION',1PG15.6)
        ELSE
           WRITE(IOUT,63) L,N,Q
   63    FORMAT(1X,'CBC  ',I6,'    NODE ',I8,'GENERATION',1PG15.6)
        ENDIF
        IBDLBL=1
      END IF
C
C15-----ADD FLOW RATE TO BUFFER.
      BUFF(N)=BUFF(N)+Q
C
C16-----SEE IF FLUX IS POSITIVE OR NEGATIVE.
      IF(QQ.GE.ZERO) THEN
C
C17-----POSITIVE FLOW RATE. ADD IT TO RATIN
        RATIN=RATIN+QQ
      ELSE
C
C18-----NEGATIVE FLOW RATE. ADD IT TO RATOUT
        RATOUT=RATOUT-QQ
      END IF
   99 CONTINUE
C
100   CONTINUE
C
C19------IF CELL-BY-CELL FLOWS WILL BE SAVED AS A 3-D ARRAY,
C19------CALL UBUDSV TO SAVE THEM.
      IF(IBD.GE.1)THEN
        IF(IUNSTR.EQ.0)THEN
          CALL UBUDSV(KSTP,KPER,TEXT(1),IBCTCB,BUFF,NCOL,NROW,
     1                          NLAY,IOUT)
        ELSE
          CALL UBUDSVU(KSTP,KPER,TEXT(1),IBCTCB,BUFF,NODES,
     1                          IOUT,PERTIM,TOTIM)
        ENDIF
        IF(INCLN.GT.0)THEN
           CALL UBUDSVU(KSTP,KPER,TEXT(2),IBCTCB,BUFF(NODES+1:NEQS),
     1                 NCLNNDS,IOUT,PERTIM,TOTIM)
        ENDIF
      ENDIF
C
C20------MOVE RATES, VOLUMES & LABELS INTO ARRAYS FOR PRINTING.
  200 RIN=RATIN
      ROUT=RATOUT
      VBVLT(3,MSUMT,ICOMP)=RIN
      VBVLT(4,MSUMT,ICOMP)=ROUT
      VBVLT(1,MSUMT,ICOMP)=VBVLT(1,MSUMT,ICOMP)+RATIN*DELT
      VBVLT(2,MSUMT,ICOMP)=VBVLT(2,MSUMT,ICOMP)+RATOUT*DELT
      VBNMT(MSUMT,ICOMP)=TEXT(1)
C
C21------INCREMENT BUDGET TERM COUNTER(MSUM).
      MSUMT=MSUMT+1
C
C22------RETURN
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE SSMS2DPT1RED(KPER)
C     ******************************************************************
C     REDUCE THE IMMOBILE DOMAIN EQUATION INTO THE MOBILE DOMAIN EQUATION
C     AT EACH DUAL POROSITY CELL FOR TRANSPORT
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,       ONLY: NODES,AMAT,RHS,IA,ISSFLG
      USE GWTDPTMODULE, ONLY: DIADDT, RDDT, OFFDDTM,OFFDDTIM,ICBUNDIM
C
C     ------------------------------------------------------------------
C-----------------------------------------------------------------------------
C1------LOOP OVER NODES TO REDUCE EQUATIONS
C-----------------------------------------------------------------------------
      DO N=1,NODES
        IF(ICBUNDIM(N).EQ.0.OR.ABS(DIADDT(N)).LT.1.0E-15) CYCLE
        IPIV = IA(N)
        AMAT(IPIV) = AMAT(IPIV) - OFFDDTIM(N)/DIADDT(N)*OFFDDTM(N)
        RHS(N) = RHS(N) - RDDT(N)/DIADDT(N)*OFFDDTM(N)
      ENDDO
C
C11----RETURN.
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE SSMS2DPT1BKS(KPER,ICOMP)
C     ******************************************************************
C     BACK-SUBSTITUTE INTO IMMOBILE DOMAIN EQUATION TO GET ITS CONC
C     AT EACH DUAL POROSITY CELL
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL, ONLY:NODES,AMAT,RHS,IA,HNEW,ISSFLG
      USE GWTBCTMODULE, ONLY: CONC
      USE GWTDPTMODULE,ONLY:DIADDT,RDDT,OFFDDTM,OFFDDTIM,CONCIM,ICBUNDIM
C
C     ------------------------------------------------------------------
C-----------------------------------------------------------------------------
C1------LOOP OVER NODES TO FILL CONCIM
C-----------------------------------------------------------------------------
      DO N=1,NODES
        IF(ICBUNDIM(N).EQ.0.OR.ABS(DIADDT(N)).LT.1.0E-15) CYCLE
C1B-----COMPUTE CONCIM FOR TRANSIENT STRESS PERIOD
        CONCIM(N,ICOMP) = RDDT(N)/DIADDT(N) -
     1        OFFDDTIM(N)*CONC(N,ICOMP)/DIADDT(N)
      ENDDO
C
C2----RETURN.
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE GWT2DPT1AD
C     ******************************************************************
C     COPY NEW INTO OLD VARIABLES FOR DUAL POROSITY TRANSPORT
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,       ONLY:NODES
      USE GWTDPTMODULE, ONLY:CONCOIM,CONCIM
      USE GWTBCTMODULE, ONLY: MCOMP
C     ------------------------------------------------------------------
C
C1----COPY NEW TO OLD VARIABLES
      DO ICOMP=1,MCOMP
      DO N=1,NODES
        CONCOIM(N,ICOMP) = CONCIM(N,ICOMP)
      ENDDO
      ENDDO
C
C4------RETURN
      RETURN
      END
C -----------------------------------------------------------------------
      SUBROUTINE GWT2DPT1OT(KSTP,KPER,ISA)
C     ******************************************************************
C     OUTPUT CONC OF IMMOBILE DOMAIN
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,      ONLY:ITMUNI,IOUT,IUNSTR
      USE GWFBASMODULE,ONLY:DELT,PERTIM,TOTIM,ISPCFL,IBUDFL
      USE GWTBCTMODULE, ONLY: MCOMP,VBVLT,VBNMT,MSUMT
C     ------------------------------------------------------------------
C
C
C1------CLEAR PRINTOUT FLAG (IPFLG)
      IPFLG=0
      DO ICOMP = 1,MCOMP
        WRITE(IOUT,1)ICOMP
1       FORMAT(5X,'TRANSPORT SOLUTION COMPLETE FOR COMPONENT SPECIES',
     *    1X,'NUMBER',I5,' IN IMMOBILE DOMAIN'/5X,80('-'))
C
C3------IF CONCENTRATION PRINT FLAG (ISPCFL) IS SET WRITE
C3------CONCENTRATION, AND ICBUND IN ACCORDANCE WITH FLAGS IN IOFLG.
        IF(ISPCFL.EQ.0) GO TO 100
C
        IF(IUNSTR.EQ.0)THEN ! WRITE M2K5 STYLE FOR STRUCTURED GRID
          CALL SGWT2DPT1C(KSTP,KPER,IPFLG,ISA,ICOMP)
        ELSE
          CALL SGWT2DPT1CU(KSTP,KPER,IPFLG,ISA,ICOMP)
        ENDIF
        IPFLG = 1
C
  100   CONTINUE

C4------PRINT TOTAL BUDGET IF REQUESTED
CPS      IF(IBUDFL.EQ.0) GO TO 120
CPS      ITF = 1   ! TRANSPORT EQUATION WAS SOLVED
CPS        CALL SGWF2BAS7V(MSUMT,VBNMT(1,ICOMP),VBVLT(1,1,ICOMP),
CPS     *  KSTP,KPER,IOUT,ITF)
CPS        IPFLG=1
CPS
C5------END PRINTOUT WITH TIME SUMMARY AND FORM FEED IF ANY PRINTOUT
C5------WILL BE PRODUCED.
  120   IF(IPFLG.EQ.0) GO TO 99
        CALL SGWF2BAS7T(KSTP,KPER,DELT,PERTIM,TOTIM,ITMUNI,IOUT)
        WRITE(IOUT,102)
  102   FORMAT('1')

      ENDDO
C
C6------RETURN
99    CONTINUE
      RETURN
      END
C----------------------------------------------------------------------
      SUBROUTINE SGWT2DPT1C(KSTP,KPER,IPFLG,ISA,ICOMP)
C     ******************************************************************
C     PRINT AND RECORD IMMOBILE DOMAIN CONCS FOR STRUCTURED GWF GRID
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,      ONLY:NCOL,NROW,NLAY,IXSEC,HNEW,NODLAY,
     1                      IBOUND,IOUT,IDPOUT
      USE GWFBASMODULE,ONLY:PERTIM,TOTIM,IOFLG,ISPCFM,LBHDSV,
     2                      CSPCFM,IOFLG,ISPCUN
      USE GWTDPTMODULE, ONLY: CONCIM,IDPTCON,ICBUNDIM
      USE GWTBCTMODULE,ONLY:IMULTI,CROOTNAME      
C
      REAL,          SAVE,    DIMENSION(:,:,:),    ALLOCATABLE ::BUFF
      REAL*8,        SAVE,    DIMENSION(:,:,:),    ALLOCATABLE ::BUFF8  
      CHARACTER*400 CFILENAME
      CHARACTER*16 TEXT, ISTATUS, FORMC
C      DATA TEXT /'   IMMOBILE CONC'/
      INCLUDE 'openspec.inc'
C     ------------------------------------------------------------------
      ALLOCATE(BUFF(NCOL,NROW,NLAY))
      ALLOCATE(BUFF8(NCOL,NROW,NLAY))
      CALL GET_TEXTCIM(ICOMP,TEXT)
C
C1------FOR EACH LAYER MOVE CONCIM TO BUFF IF PRINT OR SAVE IS REQUESTED.
      DO 59 K=1,NLAY
C
C2------IS CONC NEEDED FOR THIS LAYER?
      KL=K
      IF(IXSEC.NE.0) KL=1
      IF(IOFLG(KL,5).EQ.0 .AND. IOFLG(KL,6).EQ.0) GO TO 59
C
C3------MOVE CONCIM TO BUFF FOR THE LAYER.
      DO 58 I=1,NROW
      DO 58 J=1,NCOL
      N = (K-1)*NROW*NCOL + (I-1)*NCOL + J
      BUFF8(J,I,K)=CONCIM(N,ICOMP)
      BUFF(J,I,K) = BUFF8(J,I,K)
   58 CONTINUE
   59 CONTINUE
C
C4------FOR EACH LAYER: DETERMINE IF CONC SHOULD BE PRINTED.
C4------IF SO THEN CALL ULAPRS OR ULAPRW TO PRINT CONC.
      IF(ISA.NE.0) THEN
         IF(IXSEC.EQ.0) THEN
           DO 69 K=1,NLAY
           KK=K
           IF(IOFLG(K,5).EQ.0) GO TO 69
           IF(ISPCFM.LT.0) CALL ULAPRS(BUFF(1,1,K),TEXT,KSTP,KPER,
     1               NCOL,NROW,KK,-ISPCFM,IOUT)
           IF(ISPCFM.GE.0) CALL ULAPRW(BUFF(1,1,K),TEXT,KSTP,KPER,
     1               NCOL,NROW,KK,ISPCFM,IOUT)
           IPFLG=1
   69      CONTINUE
C
C4A-----PRINT CONC FOR CROSS SECTION.
         ELSE
           IF(IOFLG(1,5).NE.0) THEN
             IF(ISPCFM.LT.0) CALL ULAPRS(BUFF,TEXT,KSTP,KPER,
     1                 NCOL,NLAY,-1,-ISPCFM,IOUT)
             IF(ISPCFM.GE.0) CALL ULAPRW(BUFF,TEXT,KSTP,KPER,
     1                 NCOL,NLAY,-1,ISPCFM,IOUT)
             IPFLG=1
           END IF
         END IF
      END IF
C
C5------FOR EACH LAYER: DETERMINE IF CONC SHOULD BE SAVED ON DISK.
C5------IF SO THEN CALL ULASAV OR ULASV2 TO SAVE CONC.
      IFIRST=1
      IBUN = IDPTCON ! UNIT NUMBER FOR OUTPUT IS IDPTCON IF NOT MULTISPECIES
      IF(IDPTCON.LE.0) GO TO 80
C--------------------------------------------------------------------------------------------------------------------
C5A ----FOR MULTIFILE OUTPUT OPTION, ADD COMPONENT SPECIES TO ROOT NAME AND OPEN FILE ON ISPCFM AND GO TO END OF FILE
      IF(IMULTI.NE.0)THEN
        CALL GET_TEXTC(ICOMP,TEXT)  
        CFILENAME = trim(CROOTNAME) // TRIM(TEXT) //'_DPT.CON'
        ISTATUS = 'OLD'
        FORMC = 'FORMATTED'
        IF(CSPCFM.EQ.' ') FORMC = FORM      ! - use FORM from openspec.inc rather than hardwiring BINARY or UNFORMATTED  
CC        IF(KSTP.EQ.1.AND.KPER.EQ.1) THEN
C---------MAKE STATUS OF FILE AS NEW WHEN FIRST ENTERING FILE
CC          ISTATUS = 'REPLACE'
CC        ENDIF
C
!        OPEN(IDPTCON, FILE = CFILENAME, FORM = FORMC,   
        OPEN(ISPCUN, FILE = CFILENAME, FORM = FORMC,
     1  STATUS = ISTATUS, POSITION = 'APPEND')
        IBUN = ISPCUN 
      ENDIF
C-------------------------------------------------------------------------------------------------------------            
      IF(IXSEC.EQ.0) THEN
        DO 79 K=1,NLAY
        NSTRT = NODLAY(K-1)+1
        KK=K
        IF(IOFLG(K,6).EQ.0) GO TO 79
        IF(IFIRST.EQ.1) WRITE(IOUT,74) IBUN,KSTP,KPER
   74   FORMAT(1X,/1X,'IMMOBILE DOMAIN CONC WILL BE SAVED ON UNIT ',I4,
     1      ' AT END OF TIME STEP ',I3,', STRESS PERIOD ',I4)
        IFIRST=0
        IF(IDPOUT.EQ.1) THEN
           WRITE(IOUT,*)'  CONC AND TIME VARIABLES ARE SAVED AS REAL*8' 
           CALL ULASAV8(BUFF8(1,1,K),TEXT,KSTP,KPER,PERTIM,TOTIM,NCOL,
     1                NROW,KK,IBUN)        
        ELSEIF(CSPCFM.EQ.' ') THEN
           CALL ULASAV(BUFF(1,1,K),TEXT,KSTP,KPER,PERTIM,TOTIM,NCOL,
     1                NROW,KK,IBUN)
        ELSE
           CALL ULASV2(BUFF(1,1,K),TEXT,KSTP,KPER,PERTIM,TOTIM,NCOL,
     1                NROW,KK,IBUN,CSPCFM,LBHDSV,ICBUNDIM(NSTRT))
        END IF
   79   CONTINUE
C
C5A-----SAVE CONC FOR CROSS SECTION.
      ELSE
        IF(IOFLG(1,6).NE.0) THEN
          WRITE(IOUT,74) IBUN,KSTP,KPER
          IF(IDPOUT.EQ.1)THEN
           WRITE(IOUT,*)'  CONC AND TIME VARIABLES ARE SAVED AS REAL*8' 
             CALL ULASAV8(BUFF8,TEXT,KSTP,KPER,PERTIM,TOTIM,NCOL,
     1                NLAY,-1,IBUN)          
          ELSEIF(CSPCFM.EQ.' ') THEN
             CALL ULASAV(BUFF,TEXT,KSTP,KPER,PERTIM,TOTIM,NCOL,
     1                NLAY,-1,IBUN)
          ELSE
             CALL ULASV2(BUFF,TEXT,KSTP,KPER,PERTIM,TOTIM,NCOL,
     1                  NLAY,-1,IBUN,CSPCFM,LBHDSV,ICBUNDIM)
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
C --------------------------------------------------------------------------
      SUBROUTINE GET_TEXTCIM(ICOMP,TEXT)
C     ******************************************************************
C     SELECT TEXT FOR OUTPUT FILE DEPENDING ON SPECIES NUMBER
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      
      USE GWTBCTMODULE, ONLY: MCOMP
      CHARACTER*16 TEXT
C     ------------------------------------------------------------------
C
      IF(MCOMP.EQ.1) THEN
        TEXT = '   IMMOBILE CONC'  
      ELSE
        IF(ICOMP.EQ.1) THEN  
          TEXT = ' IMMOBILE CONC01'
        ELSEIF(ICOMP.EQ.2) THEN 
          TEXT = ' IMMOBILE CONC02'
        ELSEIF(ICOMP.EQ.3) THEN             
          TEXT = ' IMMOBILE CONC03'
        ELSEIF(ICOMP.EQ.4) THEN             
          TEXT = ' IMMOBILE CONC04'
        ELSEIF(ICOMP.EQ.5) THEN             
          TEXT = ' IMMOBILE CONC05'
        ELSEIF(ICOMP.EQ.6) THEN             
          TEXT = ' IMMOBILE CONC06'
        ELSEIF(ICOMP.EQ.7) THEN             
          TEXT = ' IMMOBILE CONC07'
        ELSEIF(ICOMP.EQ.8) THEN             
          TEXT = ' IMMOBILE CONC08'
        ELSEIF(ICOMP.EQ.9) THEN             
          TEXT = ' IMMOBILE CONC09'
        ELSEIF(ICOMP.EQ.10) THEN             
          TEXT = ' IMMOBILE CONC10'
        ELSEIF(ICOMP.EQ.11) THEN             
          TEXT = ' IMMOBILE CONC11'
        ELSEIF(ICOMP.EQ.12) THEN               
          TEXT = ' IMMOBILE CONC12'
        ELSEIF(ICOMP.EQ.13) THEN   
          TEXT = ' IMMOBILE CONC13'
        ELSEIF(ICOMP.EQ.14) THEN             
          TEXT = ' IMMOBILE CONC14'
        ELSEIF(ICOMP.EQ.15) THEN               
          TEXT = ' IMMOBILE CONC15'
        ELSEIF(ICOMP.EQ.16) THEN             
          TEXT = ' IMMOBILE CONC16'
        ELSEIF(ICOMP.EQ.17) THEN             
          TEXT = ' IMMOBILE CONC17'
        ELSEIF(ICOMP.EQ.18) THEN               
          TEXT = ' IMMOBILE CONC18'
        ELSEIF(ICOMP.EQ.19) THEN             
          TEXT = ' IMMOBILE CONC19'
        ELSEIF(ICOMP.EQ.20) THEN             
          TEXT = ' IMMOBILE CONC20'
        ENDIF             
      ENDIF            
C
C6------RETURN.
      RETURN
      END
C----------------------------------------------------------------------
      SUBROUTINE SGWT2DPT1CU(KSTP,KPER,IPFLG,ISA,ICOMP)
C     ******************************************************************
C     PRINT AND RECORD IMMOBILE DOMAIN CONCS FOR UNSTRUCTURED GWF GRID
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,      ONLY:NCOL,NROW,NLAY,IXSEC,HNEW,NODLAY,
     1                      IBOUND,IOUT,NODES,BUFF,IDPOUT,NEQS
      USE GWFBASMODULE,ONLY:PERTIM,TOTIM,IOFLG,ISPCFM,LBHDSV,
     2                      CSPCFM,IOFLG, ISPCUN
      USE GWTDPTMODULE, ONLY: CONCIM,IDPTCON,ICBUNDIM
      USE GWTBCTMODULE,ONLY:IMULTI,CROOTNAME
C
      REAL*8,        SAVE,    DIMENSION(:),    ALLOCATABLE ::BUFF8
      CHARACTER*400 CFILENAME
      CHARACTER*16 TEXT, ISTATUS, FORMC
      DATA TEXT /'  IMMOBILE CONCU'/ 
      INCLUDE 'openspec.inc'
C     ------------------------------------------------------------------
C
      ALLOCATE(BUFF8(NEQS))
      CALL GET_TEXTCIM(ICOMP,TEXT)
C1------FOR EACH LAYER MOVE CONC TO BUFF IF PRINT OR SAVE IS REQUESTED.
      DO 59 K=1,NLAY
C
C2------IS CONC NEEDED FOR THIS LAYER?
      KL=K
      IF(IXSEC.NE.0) KL=1
      IF(IOFLG(KL,5).EQ.0 .AND. IOFLG(KL,6).EQ.0) GO TO 59
C
C3------MOVE CONC TO BUFF FOR THE LAYER.
      NNDLAY = NODLAY(K)
      NSTRT = NODLAY(K-1)+1
      DO 58 N=NSTRT,NNDLAY
      BUFF8(N)=CONCIM(N,ICOMP)
      BUFF(N)=BUFF8(N)
   58 CONTINUE
   59 CONTINUE
C
C4------FOR EACH LAYER: DETERMINE IF CONC SHOULD BE PRINTED.
C4------IF SO THEN CALL ULAPRU TO PRINT CONC.
      IF(ISA.NE.0) THEN
         IF(IXSEC.EQ.0) THEN
           DO 69 K=1,NLAY
           KK=K
           IF(IOFLG(K,5).EQ.0) GO TO 69
           NNDLAY = NODLAY(K)
           NSTRT = NODLAY(K-1)+1
           CALL ULAPRU(BUFF,TEXT,KSTP,KPER,
     1           NSTRT,NNDLAY,KK,IABS(ISPCFM),IOUT,PERTIM,TOTIM,NODES)
           IPFLG=1
   69      CONTINUE
C
C4A-----PRINT HEAD FOR CROSS SECTION.
         ELSE
           IF(IOFLG(1,5).NE.0) THEN
           IF(ISPCFM.NE.0) CALL ULAPRU(BUFF,TEXT,KSTP,KPER,
     1           NSTRT,NNDLAY,-1,IABS(ICONFM),IOUT,PERTIM,TOTIM,NODES)
             IPFLG=1
C
           END IF
         END IF
      END IF
C
C5------FOR EACH LAYER: DETERMINE IF CONC SHOULD BE SAVED ON DISK.
C5------IF SO THEN CALL ULASAV OR ULASV2 TO SAVE HEAD.
      IFIRST=1
      IBUN = IDPTCON ! UNIT NUMBER FOR OUTPUT IS IDPTCON IF NOT MULTISPECIES      
      IF(IDPTCON.LE.0) GO TO 80
C--------------------------------------------------------------------------------------------------------------------
C5A ----FOR MULTIFILE OUTPUT OPTION, ADD COMPONENT SPECIES TO ROOT NAME AND OPEN FILE ON ISPCFM AND GO TO END OF FILE
      IF(IMULTI.NE.0)THEN
        CALL GET_TEXTC(ICOMP,TEXT)  
        CFILENAME = trim(CROOTNAME) // TRIM(TEXT) //'_DPT.CON'
        ISTATUS = 'OLD'
        FORMC = 'FORMATTED'
        IF(CSPCFM.EQ.' ') FORMC = FORM        ! - use FORM from openspec.inc rather than hardwiring BINARY or UNFORMATTED
CC        IF(KSTP.EQ.1.AND.KPER.EQ.1) THEN
C---------MAKE STATUS OF FILE AS NEW WHEN FIRST ENTERING FILE
CC          ISTATUS = 'REPLACE'
CC        ENDIF
C
!        OPEN(IDPTCON, FILE = CFILENAME, FORM = FORMC,   
        OPEN(ISPCUN, FILE = CFILENAME, FORM = FORMC,
     1  STATUS = ISTATUS, POSITION = 'APPEND')
        IBUN = ISPCUN
      ENDIF
C-------------------------------------------------------------------------------------------------------------      
      IF(IXSEC.EQ.0) THEN
        DO 79 K=1,NLAY
        KK=K
        IF(IOFLG(K,6).EQ.0) GO TO 79
        NNDLAY = NODLAY(K)
        NSTRT = NODLAY(K-1)+1
        IF(IFIRST.EQ.1) WRITE(IOUT,74) IBUN,KSTP,KPER
   74   FORMAT(1X,/1X,'IMMOBILE CONC WILL BE SAVED ON UNIT ',I8,
     1      ' AT END OF TIME STEP ',I8,', STRESS PERIOD ',I8)
        IFIRST=0
        IF(IDPOUT.EQ.1)THEN 
          WRITE(IOUT,*)'  CONC AND TIME VARIABLES ARE SAVED AS REAL*8'
          CALL ULASAVU8(BUFF8,TEXT,KSTP,KPER,PERTIM,TOTIM,NSTRT,
     1                NNDLAY,KK,IBUN,NODES)
        ELSEIF(CSPCFM.EQ.' ') THEN
           CALL ULASAVU(BUFF,TEXT,KSTP,KPER,PERTIM,TOTIM,NSTRT,
     1                NNDLAY,KK,IBUN,NODES)
        ELSE
           CALL ULASV2U(BUFF,TEXT,KSTP,KPER,PERTIM,TOTIM,NSTRT,
     1            NNDLAY,KK,IBUN,CSPCFM,LBHDSV,ICBUNDIM(NSTRT),NODES)
        END IF
        IPFLG=1
   79   CONTINUE
C
C5A-----SAVE CONC FOR CROSS SECTION.
      ELSE
        IF(IOFLG(1,6).NE.0) THEN
          WRITE(IOUT,74) IBUN,KSTP,KPER
          IF(IDPOUT.EQ.1) THEN 
            WRITE(IOUT,*)'  CONC AND TIME VARIABLES ARE SAVED AS REAL*8'
            CALL ULASAVU8(BUFF8,TEXT,KSTP,KPER,PERTIM,TOTIM,NSTRT,
     1                NNDLAY,-1,IBUN,NODES)
          ELSEIF(CSPCFM.EQ.' ') THEN
             CALL ULASAVU(BUFF,TEXT,KSTP,KPER,PERTIM,TOTIM,NSTRT,
     1                NNDLAY,-1,IBUN,NODES)
          ELSE
             CALL ULASV2U(BUFF,TEXT,KSTP,KPER,PERTIM,TOTIM,NSTRT,
     1                  NNDLAY,-1,IBUN,CSPCFM,LBHDSV,ICBUNDIM,NODES)
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
C----------------------------------------------------------------------
      SUBROUTINE SOLUBILITYIM(N,ICOMP,DTERMS,RTERMS,ISS)

C     ******************************************************************
C     FORMULATE SOLUBILITY TERM FOR EACH NODE AND COMPONENT IN IMMOBILE DOMAIN
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,   ONLY:IOUT,NCOL,NROW,NLAY,IBOUND,TOP,NODLAY,ISYM,
     1 HNEW,AMAT,RHS,HOLD,ISSFLG,IA,JA,NJA,NODES,AREA,BOT,TOP,Sn,So
      USE CLN1MODULE, ONLY: ACLNNDS
      USE GWFBASMODULE,ONLY:DELT
      USE GWTDPTMODULE, ONLY: CONCIM, CONCOIM
      USE GWTBCTMODULE, ONLY: MCOMP,CONC,CONCO,ICT,IADSORB,PRSITY,
     *  ADSORB,IZOD,IFOD,ZODRW,FODRW,ZODRS,FODRS,MCOMPT,IHEAT,
     *  ISOLUBILITY,SOLLIM,SOLSLOPE 
      DOUBLE PRECISION VODT,CW,CWO,DT,RT,ALENG,VOLU,SLIM,SSLOPE,FN,DF,
     1  X,Y,EPSS,EPS,CT,CEPS,QA,QEPS,DQ,OMP,DTERMS,RTERMS,EMOLD,EMNEW
C     ------------------------------------------------------------------
C1------INITIALIZE TERMS
      DTERMS = 0.0
      RTERMS = 0.0
      IF(N.LE.NODES)THEN
        ALENG = TOP(N) - BOT(N)
      ELSE
        ALENG = ACLNNDS(N-NODES,4)
      ENDIF
      VOLU = AREA(N) * ALENG
      VODT = VOLU / DELT
C----------------------------------------------------------------------
      IF(ICT.EQ.0)THEN  !----------WATER PHASE CONCENTRATION FORMULATION
C2---------NONLINEAR SO FILLED AS NEWTON (KINK NOT REMOVED)
C          SSLOPE = SOLSLOPE(ICOMP) * VODT
C          SLIM = SOLLIM(ICOMP)
C          CW = CONC(N,ICOMP)
C          CWO = CONCO(N,ICOMP)
CC         
C          DT = 0.0
C          RT = 0.0
C          IF(CWO.GT.SLIM)THEN
C            RT = -CWO
C          ENDIF
C          IF(CW.GT.SLIM)THEN
C            DT = SSLOPE 
C            RT = RT + CW
C          ENDIF
C          RT = RT * SSLOPE
C          DTERMS = DTERMS - DT
C          RTERMS = RTERMS - DT * CW + RT          
C ----------------------------------------------------------------          
          SSLOPE = SOLSLOPE(ICOMP) 
          SLIM = SOLLIM(ICOMP)
          CW = CONCIM(N,ICOMP)
          CWO = CONCOIM(N,ICOMP)          
          EMOLD = 0.0
          IF(CWO.GT.SLIM) EMOLD = (CWO-SLIM) * SSLOPE
          EMNEW = 0.0
          IF(CW.GT.SLIM) EMNEW = (CW-SLIM) * SSLOPE
          FN = (EMNEW - EMOLD) * VODT
          DF = EMNEW
          CW = CW + 1.0E-4
          EMNEW = 0.0
          IF(CW.GT.SLIM) EMNEW = (CW-SLIM) * SSLOPE
          DF = (EMNEW - DF) / 1.0E-4 
          DF = DF * VODT
          RTERMS = FN - DF * CW
          DTERMS = - DF 
C-----------------------------------------------------------------------
      ELSE       !-----------------------TOTAL CONCENTRATION FORMULATION
C5-------NET STORAGE TERM FOR TOTAL CONCENTRATION FORMULATION
        DTERMS = DTERMS - VODT
        RTERMS = RTERMS - VODT * CONCOIM(N,ICOMP)
      ENDIF
C
C9------RETURN
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE SOLUBILITYMBIM(KSTP,KPER,ICOMP,ISS)
C     ******************************************************************
C     CALCULATE MASS TERMS OF PRECIPITATED SOLUTES FOR ALL TRANSPORT CELLS
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,   ONLY:IOUT,NCOL,NROW,NLAY,IBOUND,BUFF,NODES,IUNSTR,
     1  AMAT,IA,JA,TOP,BOT,AREA,Sn,So,NEQS,INCLN,IUNIT
      USE CLN1MODULE, ONLY: ACLNNDS,NCLNNDS
      USE GWFBASMODULE,ONLY:MSUM,ICBCFL,IAUXSV,DELT,PERTIM,TOTIM,DELT
      USE GWTDPTMODULE, ONLY: CONCIM, CONCOIM
      USE GWTBCTMODULE, ONLY: ICBUND,CONC,MSUMT,VBVLT,VBNMT,IBCTCB,
     1 IADSORB,ADSORB,PRSITY,CONCO,ICT,MCOMPT,
     *  ISOLUBILITY,SOLLIM,SOLSLOPE  
C
      CHARACTER*16 TEXT(1)
      DOUBLE PRECISION RATIN,RATOUT,QQ,VODT,ADSTERM,FL,CW,CWO,ALENG,
     *  DTERMS,RTERMS,VOLU,RATINTVM,RATOUTTVM,
     *  EMNEW,EMOLD,FN,SLIM,SSLOPE
      DATA TEXT(1) /'  IM MASS PRECIP'/
C     ------------------------------------------------------------------
C
C1------CLEAR RATIN AND RATOUT ACCUMULATORS, AND SET CELL-BY-CELL
C1------BUDGET FLAG.
      ZERO=0.
      RATIN=ZERO
      RATOUT=ZERO
      RATINTVM=ZERO
      RATOUTTVM=ZERO
      IBD=0
      IF(IBCTCB.LT.0 .AND. ICBCFL.NE.0) IBD=-1
      IF(IBCTCB.GT.0) IBD=ICBCFL
      IBDLBL=0
C
C2------CLEAR THE BUFFER.
      DO 50 N=1,NODES
        BUFF(N)=ZERO
50    CONTINUE
C
C3------LOOP THROUGH EACH NODE AND CALCULATE STORAGE
      DO 100 N=1,NODES
C
C4-----IF THE CELL IS NOT PCB OR WRONG COMPONENT SPECIES, IGNORE IT.
        IF(ICBUND(N).EQ.0)GO TO 99
C
        IF(N.LE.NODES)THEN
          ALENG = TOP(N) - BOT(N)
        ELSE
          ALENG = ACLNNDS(N-NODES,4)
        ENDIF
        VOLU = AREA(N) * ALENG
        VODT = VOLU / DELT
        QQ = 0.0
        IF(ICT.EQ.0)THEN
C5-------PRECIPITATED MASS - SKIP IF IMMOBILE COMPONENTS
          IF(ICOMP.GT.MCOMPT) GO TO 252  
CC        
          SSLOPE = SOLSLOPE(ICOMP) 
          SLIM = SOLLIM(ICOMP)
          CW = CONCIM(N,ICOMP)
          CWO = CONCOIM(N,ICOMP)          
          EMOLD = 0.0
          IF(CWO.GT.SLIM) EMOLD = (CWO-SLIM) * SSLOPE
          EMNEW = 0.0
          IF(CW.GT.SLIM) EMNEW = (CW-SLIM) * SSLOPE
          FN = (EMNEW - EMOLD) * VODT
          DF = EMNEW
          CW = CW + 1.0E-4
          EMNEW = 0.0
          IF(CW.GT.SLIM) EMNEW = (CW-SLIM) * SSLOPE
          DF = (EMNEW - DF) / 1.0E-4 
          DF = DF * VODT
          RTERMS = FN - DF * CW
          DTERMS = - DF 
          qq = rterms - dterms * concIM(n,icomp) 
252       CONTINUE        
C-----------------------------------------------------
        ELSE   !-----------------------TOTAL CONCENTRATION FORMULATION
C7-------NET STORAGE TERM FOR TOTAL CONCENTRATION FORMULATION
          QQ = QQ + VODT * CONCIM(N,ICOMP)
     *          - VODT * CONCOIM(N,ICOMP)
        ENDIF
        QQ = - QQ  ! STORAGE TERM NEGATIVE IS INFLOW AS PER MODFLOW CONVENTION
        Q = QQ
C
C8------PRINT FLOW RATE IF REQUESTED.
        IF(IBD.LT.0) THEN
          IF(IBDLBL.EQ.0) WRITE(IOUT,61) TEXT(1),KPER,KSTP
   61     FORMAT(1X,/1X,A,'   PERIOD ',I4,'   STEP ',I3)
          IF(IUNSTR.EQ.0)THEN
            IL = (N-1) / (NCOL*NROW) + 1
            IJ = N - (IL-1)*NCOL*NROW
            IR = (IJ-1)/NCOL + 1
            IC = IJ - (IR-1)*NCOL
            WRITE(IOUT,62) IL,IR,IC,Q
   62       FORMAT(1X,'   LAYER ',I5,'   ROW ',I6,'   COL ',I6,
     1       '   FLUX ',1PG15.6)
          ELSE
            WRITE(IOUT,63) N,Q
   63      FORMAT(1X,'    NODE ',I8,'   FLUX ',1PG15.6)
          ENDIF
          IBDLBL=1
        END IF
C
C9------ADD FLOW RATE TO BUFFER.
        BUFF(N)=BUFF(N)+QQ
C
C10-----SEE IF FLUX IS POSITIVE OR NEGATIVE.
        IF(QQ.GE.ZERO) THEN
C
C11-----POSITIVE FLOW RATE. ADD IT TO RATIN
          RATIN=RATIN+QQ
        ELSE
C
C12-----NEGATIVE FLOW RATE. ADD IT TO RATOUT
          RATOUT=RATOUT-QQ
        END IF
   99   CONTINUE
C
100   CONTINUE
C
C13------IF CELL-BY-CELL FLOWS WILL BE SAVED AS A 3-D ARRAY,
C13------CALL UBUDSV TO SAVE THEM.
      IF(IBD.GE.1)THEN
        IF(IUNSTR.EQ.0)THEN
          CALL UBUDSV(KSTP,KPER,TEXT(1),IBCTCB,BUFF,NCOL,NROW,
     1                          NLAY,IOUT)
        ELSE
          CALL UBUDSVU(KSTP,KPER,TEXT(1),IBCTCB,BUFF,NODES,
     1                          IOUT,PERTIM,TOTIM)
        ENDIF
      ENDIF  
C
C14------MOVE RATES, VOLUMES & LABELS INTO ARRAYS FOR PRINTING.
  200 RIN=RATIN
      ROUT=RATOUT
      VBVLT(3,MSUMT,ICOMP)=RIN
      VBVLT(4,MSUMT,ICOMP)=ROUT
      VBVLT(1,MSUMT,ICOMP)=VBVLT(1,MSUMT,ICOMP)+RATIN*DELT
      VBVLT(2,MSUMT,ICOMP)=VBVLT(2,MSUMT,ICOMP)+RATOUT*DELT
      VBNMT(MSUMT,ICOMP)=TEXT(1)
C
C15------INCREMENT BUDGET TERM COUNTER(MSUM).
      MSUMT=MSUMT+1
C
C16------RETURN
      RETURN
      END
C ------------------------------------------------------------
      SUBROUTINE GWT2DPT1DA
      USE GWTDPTMODULE
      USE GWFDPFMODULE, ONLY: CBCFIM
      USE AW_ADSORBIMMODULE
      INTEGER ALLOC_ERR
C
      DEALLOCATE(ICBUNDIM, STAT = ALLOC_ERR)
      DEALLOCATE(PRSITYIM, STAT = ALLOC_ERR)
      DEALLOCATE(CONCIM, STAT = ALLOC_ERR)
      DEALLOCATE(CONCOIM, STAT = ALLOC_ERR)
      DEALLOCATE(MASSBCTIM, STAT = ALLOC_ERR)
      DEALLOCATE(CBCFIM, STAT = ALLOC_ERR)

      DEALLOCATE(ADSORBIM, STAT = ALLOC_ERR)
      DEALLOCATE(FLICHIM, STAT = ALLOC_ERR)
      DEALLOCATE(ZODRWIM, STAT = ALLOC_ERR)
      DEALLOCATE(ZODRSIM, STAT = ALLOC_ERR)
      DEALLOCATE(FODRWIM, STAT = ALLOC_ERR)
      DEALLOCATE(FODRSIM, STAT = ALLOC_ERR)
      DEALLOCATE(IADSORBIM, STAT = ALLOC_ERR)
      DEALLOCATE(IDISPIM, STAT = ALLOC_ERR)
      DEALLOCATE(IFODIM, IZODIM, STAT = ALLOC_ERR)

      IF(IAW_ADSORBIM.EQ.1) THEN 
         DEALLOCATE (IAREA_FNIM, STAT = ALLOC_ERR)
         DEALLOCATE (IKAWI_FNIM, STAT = ALLOC_ERR)
         DEALLOCATE (ITAB_AWIIM, STAT = ALLOC_ERR)
         DEALLOCATE (ZODRAWIM, STAT = ALLOC_ERR)
         DEALLOCATE (FODRAWIM, STAT = ALLOC_ERR)

         CALL DPT2AWI_ADSORB1DA
      ENDIF
       DEALLOCATE (IAW_ADSORBIM, STAT = ALLOC_ERR)
C
      RETURN
      END

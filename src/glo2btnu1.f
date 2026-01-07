c -------------------------------------------------------------------------------------
      SUBROUTINE GWT2BCT1AR(IN,INDPT)
C     ******************************************************************
C     ALLOCATE ARRAYS AND READ DATA FOR BLOCK-CENTERED TRANSPORT PACKAGE
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL, ONLY:IOUT,NCOL,NROW,NLAY,ITRSS,IFREFM,IUNSTR,
     * NODES,NEQS,NJA,NJAS,ITRNSP,Sn,NIUNIT,INCLN,NODLAY,AREA,IUNIT,
     *   IFMBC,MBEGWUNF,MBEGWUNT,MBECLNUNF,MBECLNUNT,NPER,ISSFLG,
     *   IOUTNORMAL,IDPT
      USE GWFBASMODULE,ONLY:ISPCUN
      USE GWFBCFMODULE, ONLY: LAYCON, SC2,IHANISO
      USE GWTBCTMODULE
      USE CLN1MODULE, ONLY: NCLNGWC
      USE AW_ADSORBMODULE, ONLY : IAREA_FN,IKAWI_FN,ITAB_AWI 
C
      CHARACTER*16 TEXT, ISTATUS
      CHARACTER*400 CFILENAME
      CHARACTER*300 LINE
      INTEGER IC_IBOUND_FLG
C     kkz - include openspec.inc
      INCLUDE 'openspec.inc'
C
C     ------------------------------------------------------------------
C1------ALLOCATE SCALAR VARIABLES IN FORTRAN MODULE.
      ALLOCATE(IBCTCB,MCOMP,ITVD,IADSORB,ICT,NTITER,IDISP,IXDISP,
     1  IZOD,IFOD,IHEAT,MCOMPT,NTCOMP,IDISPCLN,IMULTI,NSEQITR)
      ALLOCATE(CINACT)
      ALLOCATE(CICLOSE)
      ALLOCATE(DIFFNC)
      ALLOCATE(TIMEWEIGHT)
      ALLOCATE (ICHAIN,ISPRCT)
      ALLOCATE (ISOLUBILITY)
      ALLOCATE (ISATADSORB)      
      ALLOCATE (IAW_ADSORB)
      ALLOCATE (IMASSWR)
      IF(IHANISO.EQ.0) THEN 
        ALLOCATE (IOUTNORMAL)  ! NOT ALLOCATED YET SO ALLOCATE
        IOUTNORMAL = 0
      ENDIF    
      IHEAT = 0
      NTCOMP = 0
      IDISPCLN = 0
      NSEQITR = 0
      TIMEWEIGHT = 1.0
      ICHAIN = 0
      ISPRCT = 0
      ISOLUBILITY = 0
      ISATADSORB = 0
      IAW_ADSORB = 0
      IMASSWR = 0
C
C1A------IDENTIFY PACKAGE
      WRITE(IOUT,1) IN
    1 FORMAT(1X,/1X,'BCT -- BLOCK-CENTERED TRANSPORT PACKAGE, VER. 8',
     1', 9/12/2009',/,9X,'INPUT READ FROM UNIT',I3)
C
C2-------CHECK FOR OPTIONS
      CALL URDCOM(IN,IOUT,LINE)
      IMULTI = 0  !INDEX FOR MULTIPLE OUTPUT FILES
      LLOC = 1
   10 CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,N,R,IOUT,IN)
      IF(LINE(ISTART:ISTOP).EQ.'MULTIFILE') THEN
C3B------READ KEYWORD OPTION TO OUTPUT BINARY CONC IN SEPARATE FILES.
        IMULTI = 1
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,N,R,IOUT,IN)
        CROOTNAME = LINE(ISTART:ISTOP)
        WRITE(IOUT,36) CROOTNAME
36      FORMAT(1X,'ROOT FILE NAME FOR CONCENTRATION OUTPUT = ', A20)
      ELSEIF(LINE(ISTART:ISTOP).EQ.'TIMEWEIGHT') THEN
C3C-----READ KEYWORD OPTION TO USE ALTERNATE TIME WEIGHTING
        TIMEWEIGHT = 1.0
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,I,TW,IOUT,IN)
        TIMEWEIGHT = TW
        WRITE(IOUT,37) TW
37      FORMAT(1X,'TIME WEIGHTING FACTOR (TIMEWEIGHT) = ', F10.3)
      ELSEIF(LINE(ISTART:ISTOP).EQ.'CHAINDECAY') THEN
        ICHAIN = 1
        WRITE(IOUT,38)
38      FORMAT(1X,'CHAINDECAY OPTION IS TURNED ON (ICHAIN = 1)')
      ELSEIF(LINE(ISTART:ISTOP).EQ.'SPATIALREACT') THEN
        ISPRCT = 1
        WRITE(IOUT,39)
39      FORMAT(1X,'REACTION STOCHIOMETRY IS SPATIALLY DISTRIBUTED')
      ELSEIF(LINE(ISTART:ISTOP).EQ.'SOLUBILITY') THEN
        ISOLUBILITY = 1
        WRITE(IOUT,40)
40      FORMAT(1X,'SOLUBILITY LIMITS ARE TURNED ON FOR THE SOLUTES')
      ELSEIF(LINE(ISTART:ISTOP).EQ.'VARIABLE_SOLUBILITY') THEN
        ISOLUBILITY = 2
        WRITE(IOUT,42)
42      FORMAT(1X,'SOLUBILITY LIMITS ARE TURNED ON FOR THE SOLUTES',1X,
     1   'WITH SOLUBILITY PARAMETERS VARYING OVER SPACE') 
csp      ELSEIF(LINE(ISTART:ISTOP).EQ. 'OUTWARDNORMAL') THEN 
csp        IOUTNORMAL = 1
csp        WRITE(IOUT,41)
csp41      FORMAT(1X,'ANGLE USED FOR EACH FACE IS OUTWARD NORMAL TO CELL')
      ELSEIF(LINE(ISTART:ISTOP).EQ.'ONLY_SATADSORB') THEN
        ISATADSORB = 1
        WRITE(IOUT,41)
41      FORMAT(1X,'ADSORPTION WILL OCCUR ONLY IN SATURATED PORTION OF', 
     1   ' GRID BLOCKS (ONLY_SATADSORB OPTION) FOR STEADY-STATE FLOW',
     2   ' & SP1-TS1 SAT THICKNESS FOR TRANSIENT FLOW TO CONSERVE MASS') 
c        
      ELSE IF(LINE(ISTART:ISTOP).EQ.'A-W_ADSORB') THEN
         IAW_ADSORB=1      
         ALLOCATE(IAREA_FN,IKAWI_FN,ITAB_AWI)   
         ITAB_AWI = 0
         WRITE(IOUT,33)
33       FORMAT(1X,'AIR-WATER INTERFACE ADSORPTION OPTION') 
         CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IAREA_FN,R,IOUT,IN)
         WRITE(IOUT,46) IAREA_FN
46       FORMAT(1X,'AREA FUNCTION INDEX (IAREA_FN) = ', I10)
         if(iarea_fn.eq.3)then
          write(iout,*)('iarea_fn = 3 not available, stopping ...')
          stop
         endif  
         CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IKAWI_FN,R,IOUT,IN)
         WRITE(IOUT,47) IKAWI_FN
47       FORMAT(1X,'K_WAI FUNCTION INDEX (IKAWI_FN) = ', I10) 
C         
       ELSEIF(LINE(ISTART:ISTOP).EQ.'WRITE_GWMASS') THEN
         CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IMASSWR,R,IOUT,IN)
         WRITE(IOUT,48) IMASSWR
48       FORMAT(1X,'UNIT NUMBER ON WHICH TOTAL MASS IN EACH CELL IS',1X,
     1    'SAVED ', I10)          
c             
       ELSEIF(LINE(ISTART:ISTOP).EQ.'SOMETHINGELSE') THEN

      END IF
      IF(LLOC.LT.200) GO TO 10
C
C3------READ AND PRINT IBCTCB (FLAG FOR PRINTING
C3------OR UNIT# FOR RECORDING CELL-BY-CELL FLUX TERMS), MCOMP
C3------(NUMBER OF MOBIME COMPONENTS), AND CINACT (C AT INACTIVE CELLS)
      IF(IFREFM.EQ.0) THEN
          READ(LINE,2)ITRNSP,IBCTCB,MCOMP,IC_IBOUND_FLG,ITVD,IADSORB,
     *     ICT,CINACT,CICLOSE,IDISP,IXDISP,DIFFNC,IZOD,IFOD,IFMBC,IHEAT,
     *     IMCOMP,IDISPCLN,NSEQITR
      ELSE
        LLOC = 1
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,ITRNSP,R,IOUT,IN)
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IBCTCB,R,IOUT,IN)
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,MCOMP,R,IOUT,IN)
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IC_IBOUND_FLG,R,IOUT,IN)
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,ITVD,R,IOUT,IN)
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IADSORB,R,IOUT,IN)
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,ICT,R,IOUT,IN)
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,I,CINACT,IOUT,IN)
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,I,CICLOSES,IOUT,IN)
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IDISP,R,IOUT,IN)
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IXDISP,R,IOUT,IN)
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,I,DIFFNC,IOUT,IN)
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IZOD,R,IOUT,IN)
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IFOD,R,IOUT,IN)
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IFMBC,R,IOUT,IN)
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IHEAT,R,IOUT,IN)
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IMCOMP,R,IOUT,IN)
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IDISPCLN,R,IOUT,IN)
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NSEQITR,R,IOUT,IN)
        CICLOSE = CICLOSES
        NTCOMP = IMCOMP + MCOMP
      ENDIF
2     FORMAT(7I10,2F10.3,2I10,F10.3,7I10)
      MCOMPT = MCOMP ! MCOMPT is mobile components with temperature
      IF(NTCOMP.LT.MCOMP) NTCOMP = MCOMP
      IF(IHEAT.EQ.1) THEN
        MCOMPT = MCOMPT + 1
        NTCOMP = NTCOMP + 1
      ENDIF
      NTITER = 1
      IF(ITVD.NE.0)NTITER = ABS(ITVD)
      IF(CICLOSE.LT.1.0E-20) CICLOSE = 1.0E-6
      IF(NSEQITR.LT.1) NSEQITR = 1
C
C3B-----PRINT VALUES
      WRITE(IOUT,11) ITRNSP,IBCTCB,MCOMP,IC_IBOUND_FLG,ITVD,IADSORB,ICT,
     *  CINACT,CICLOSE,IDISP,IXDISP,DIFFNC,IZOD,IFOD,IFMBC,IHEAT,IMCOMP,
     *  IDISPCLN,NSEQITR
   11 FORMAT(1X,'TRANSPORT FLAG (ITRNSP)                             =',
     *  I3/1X,'C-B-C FLUX FLAG OR UNIT NUMBER (IBCTCB)             =',I3
     *  /1X,'NUMBER OF MOBILE COMPONENTS (MCOMP)                 =',I3
     *  /1X,'ACTIVE DOMAIN SAME AS FOR FLOW FLAG (IC_IBOUND_FLG) =',I3
     *  /1X,'TVD SCHEME FLAG AND FCT COUNT (ITVD)                =',I3
     *  /1X,'ADSORPTION FLAG (IADSORB)                           =',I3
     *  /1X,'FLAG TO SOLVE FOR TOTAL CONCENTRATION (ICT)         =',I3
     */1X,'CONCENTRATION AT INACTIVE CELLS (CINACT)            =',G13.5,
     */1X,'SOLVER TOLERANCE FOR TRANSPORT (CICLOSE)            =',G13.5,
     *  /1X,'DISPERSION FLAG (IDISP)                             =',I3,
     *  /1X,'CROSS-DISPERSION FLAG (IXDISP)                      =',I3,
     */1X,'DIFFUSION COEFFICIENT (DIFFNC)                      =',G13.5,
     *  /1X,'ZERO-ORDER DECAY FLAG (IZOD)                        =',I3,
     *  /1X,'FIRST-ORDER DECAY FLAG (IFOD)                       =',I3,
     *  /1X,'FLOW ERROR CORRECTION FLAG OR UNIT NUMBER (IFMBC)   =',I3
     *  /1X,'INDEX FOR SIMULATION OF HEAT TRANSPORT (IHEAT)      =',I3,
     *  /1X,' NUMBER OF IMMOBILE COMPONENTS  (IMCOMP)             =',I3,
     *  /1x,'INDEX FOR CLN-GW DISPERSION EQUATION (IDISPCLN)     =',I3,
     *  /1x,'NUMBER OF SEQUENTIAL ITERATIONS (NSEQITR)           =',I3)
C
C3C-----IF MULTIPLE .CON OUTPUT FILES THEN OPEN THEM HERE AS NEW
      IF(IMULTI.NE.0)THEN
        DO ICOMP = 1,NTCOMP
          CALL GET_TEXTC(ICOMP,TEXT)
          CFILENAME = trim(CROOTNAME) // TRIM(TEXT) //'.CON'
!          OPEN(ISPCUN, FILE = CFILENAME, FORM = 'BINARY',  
          OPEN(ISPCUN, FILE = CFILENAME, FORM = FORM,     ! - use FORM from openspec.inc rather than hardwiring BINARY or UNFORMATTED
     1     STATUS = 'REPLACE')
          CLOSE (ISPCUN)
        ENDDO
C3D-----ALSO FOR CLN IF PRESENT
        IF(INCLN.GT.0) THEN
          DO ICOMP = 1,NTCOMP
            CALL GET_TEXTC(ICOMP,TEXT)
            CFILENAME = trim(CROOTNAME) // TRIM(TEXT) //'_CLN.CON'
!            OPEN(ISPCUN, FILE = CFILENAME, FORM = 'BINARY',  
           OPEN(ISPCUN, FILE = CFILENAME, FORM = FORM,      ! - use FORM from openspec.inc rather than hardwiring BINARY or UNFORMATTED
     1     STATUS = 'REPLACE')
            CLOSE (ISPCUN)
          ENDDO
        ENDIF
C3E-----ALSO FOR DPT IF PRESENT
        IF(INDPT.GT.0) THEN
          DO ICOMP = 1,NTCOMP
            CALL GET_TEXTC(ICOMP,TEXT)
            CFILENAME = trim(CROOTNAME) // TRIM(TEXT) //'_DPT.CON'
!            OPEN(ISPCUN, FILE = CFILENAME, FORM = 'BINARY',  
            OPEN(ISPCUN, FILE = CFILENAME, FORM = FORM,  ! - use FORM from openspec.inc rather than hardwiring BINARY or UNFORMATTED
     1     STATUS = 'REPLACE')
            CLOSE (ISPCUN)
          ENDDO 
        ENDIF
C3F-----ALSO FOR md IF PRESENT
C         IF(ImdT.GT.0) THEN

C        ENDIF
C3G-------OTHER DOMAINS IF PRESENT HERE
C
        IF(IMASSWR.GE.1)THEN        
        DO ICOMP = 1,NTCOMP
          CALL GET_TEXTM(ICOMP,TEXT)
          CFILENAME = trim(CROOTNAME) // TRIM(TEXT) //'.MAS'
!          OPEN(IMASSWR, FILE = CFILENAME, FORM = 'BINARY',  
          OPEN(IMASSWR, FILE = CFILENAME, FORM = FORM,     ! - use FORM from openspec.inc rather than hardwiring BINARY or UNFORMATTED
     1     STATUS = 'REPLACE')
          CLOSE (IMASSWR)
        ENDDO
        ENDIF
      ENDIF
C -----------------------------------------------------------------------
      IF(ITRNSP.EQ.1)THEN
        WRITE(IOUT,*)' TRANSPORT SIMULATION FOLLOWS FLOW',
     *    ' (ITRNSP = 1)'
      ELSEIF(ITRNSP.EQ.2)THEN
        WRITE(IOUT,*)' RUN TRANSPORT ON PREVIOUSLY RUN FLOW-FIELD',
     1  ' (ITRNSP = 2)'
c        WRITE(IOUT,*)'OPTION CURRENTLY NOT AVAILABLE; (ITRNSP = 2)'
c        STOP
      ELSEIF(ITRNSP.EQ.3)THEN
        WRITE(IOUT,*)' ONLY FIRST SP AND TIME OF FLOW IS SIMULATED; ',
     *    'REMAINING SP AND TIMES ARE FOR TRANSPORT ONLY (ITRNSP = 3)'
      ELSEIF(ITRNSP.EQ.4)THEN
        WRITE(IOUT,*)' RUN TRANSPORT ON PREVIOUSLY RUN FLOW-FIELD',
     1    'OF ONLY FIRST SP AND TIME   (ITRNSP = 4)'
      ELSEIF(ITRNSP.EQ.5)THEN
        WRITE(IOUT,*)' BYPASS FLOW AND RUN TRANSPORT ONLY FOR NO-FLOW',
     1    ' WITH INITIAL CONCS TREATED AS PCB (ITRNSP = 5)'
        IF(IUNIT(14). EQ. 0) IUNIT(14) = 9999999
      ENDIF
      IF(ICT.EQ.1) THEN
        WRITE(IOUT,*)' SOLVING FOR TOTAL CONCENTRATION'
        WRITE(IOUT,*)'OPTION CURRENTLY NOT AVAILABLE'
        STOP
      ENDIF
C
C4-----ALLOCATE SPACE FOR ARRAYS.
      ALLOCATE(ICBUND(NEQS))
      ALLOCATE(PRSITY(NODES),
     *  CBCF(NJAS),CBCH(NEQS),ADMAT(NJA))
      ALLOCATE (IPCBFLAG(NEQS))
      ALLOCATE (ADSORB(NODES,MCOMPT))
      ALLOCATE (FLICH(NODES,MCOMP))
      FLICH = 1.0
      ALLOCATE (Cncg(NTITER),LrcC(3,NTITER))
C      IF(INCLN.GT.0)
      ALLOCATE (FCncg(NTITER),FLrcC(NTITER))
      IF(IDISP.GT.0)THEN
        ALLOCATE(VELNOD(NEQS,3))
        ALLOCATE(VXL(NJA),VYL(NJA))        
        ALLOCATE(DLX(NEQS), DLY(NODES), DLZ(NODES))
        ALLOCATE(ATXY(NEQS), ATYZ(NODES), ATXZ(NODES))
        IF(IXDISP.NE.0)THEN
          ALLOCATE(ALXY(NEQS), ALYZ(NODES), ALXZ(NODES))
        ENDIF
      ENDIF
      IF(IZOD.EQ.1)THEN
        ALLOCATE(ZODRW(NEQS,MCOMP))
        ZODRW = 0.
      ELSE
        IF(IADSORB.GT.0)THEN
          IF(IZOD.EQ.2)THEN
             ALLOCATE(ZODRS(NEQS,MCOMP))
             ZODRS = 0.
          ELSEIF(IZOD.EQ.3)THEN
            ALLOCATE(ZODRW(NEQS,MCOMP))
            ZODRW = 0.
            ALLOCATE(ZODRS(NEQS,MCOMP))
            ZODRS = 0.
          ENDIF
        ENDIF
        IF(IAW_ADSORB.GT.0) THEN
          IF(IZOD.EQ.4)THEN
            ALLOCATE(ZODRW(NEQS,MCOMP))
            ZODRW = 0.
            ALLOCATE(ZODRS(NEQS,MCOMP))
            ZODRS = 0.
            ALLOCATE(ZODRAW(NEQS,MCOMP))
            ZODRAW = 0.           
          ENDIF  
        ENDIF    
      ENDIF
      IF(IFOD.EQ.1)THEN
        ALLOCATE(FODRW(NEQS,MCOMP))
        FODRW = 0.
      ELSE
        IF(IADSORB.GT.0)THEN
          IF(IFOD.EQ.2)THEN
             ALLOCATE(FODRS(NEQS,MCOMP))
             FODRS = 0.
          ELSEIF(IFOD.EQ.3)THEN
            ALLOCATE(FODRW(NEQS,MCOMP))
            FODRW = 0.
            ALLOCATE(FODRS(NEQS,MCOMP))
            FODRS = 0.
          ENDIF
        ENDIF
        IF(IAW_ADSORB.GT.0) THEN
          IF(IFOD.EQ.4)THEN
            ALLOCATE(FODRW(NEQS,MCOMP))
            FODRW = 0.
            ALLOCATE(FODRS(NEQS,MCOMP))
            FODRS = 0.
            ALLOCATE(FODRAW(NEQS,MCOMP))
            FODRAW = 0.           
          ENDIF  
        ENDIF    
      ENDIF
      IF(IFMBC.NE.0)THEN
        ALLOCATE(MBEGWUNF,MBEGWUNT,MBECLNUNF,MBECLNUNT)
        MBEGWUNF = 0
        MBEGWUNT = 0
        MBECLNUNF = 0
        MBECLNUNT = 0
C4A1-------...AND READ UNIT NUMBERS FOR FMBE OUTPUT FILES IF FMBC IS ON
        IF(IFREFM.EQ.0) THEN
          READ(IN,2)MBEGWUNF,MBEGWUNT,MBECLNUNF,MBECLNUNT
        ELSE
          CALL URDCOM(IN,IOUT,LINE)
          LLOC=1
          CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,MBEGWUNF,R,IOUT,IN)
          CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,MBEGWUNT,R,IOUT,IN)
          CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,MBECLNUNF,R,IOUT,IN)
          CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,MBECLNUNT,R,IOUT,IN)
        ENDIF
        WRITE(IOUT,12) MBEGWUNF,MBEGWUNT,MBECLNUNF,MBECLNUNT
12    FORMAT(1X,'UNIT NUMBER FOR OUTPUT OF FLOW IMBALANCE (MBEGWUNF) =',
     * /1X,'UNIT NUMBER FOR OUTPUT OF TRANSPORT IMBAL(MBEGWUNT) =',I3
     *  I3/1X,'UNNIT NUMBER FOR OUTPUT OF CLN FLOW IMBAL(MBECLNUNF)=',I3
     *  /1X,'UNIT NUMBER FOR OUTPUT OF TRANSP IMBAL  (MBECLNUNT) =',I3)
      ENDIF
C4A2-------IF HEAT EQUATION IS SOLVED THEN ALLOCATE SPACE FOR RELEVANT ARRAYS AND VARIABLES
      IF(IHEAT.EQ.1) THEN
        ALLOCATE(HTCONDW, RHOW, HTCAPW) ! HEAT CONDUCTIVITY, DENSITY AND HEAT CAPACITY OF WATER ARE VARIABLES
        ALLOCATE(HTCONDS(NEQS), HTCONDM(NEQS)) !HEAT CONDUCTIVITY OF SOIL AND MEDIUM ARE ARRAYS (HEAT CAPACITY IS READ IN ARRAY ADSORB)
        ALLOCATE(HTCONDT(NJA)) !INTERCELL HEAT CONDUCTANCE ARRAY
        HTCONDT = 0.0
        IF(IFREFM.EQ.0) THEN
!          READ(IN,2) HTCONDW, RHOW, HTCAPW
          READ(IN,*) HTCONDW, RHOW, HTCAPW  !kkz - format 2 not applicable to this read per JCH
        ELSE
          CALL URDCOM(IN,IOUT,LINE)
          LLOC=1
          CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,I,HTCONDW,IOUT,IN)
          CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,I,RHOW,IOUT,IN)
          CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,I,HTCAPW,IOUT,IN)
        ENDIF
        WRITE(IOUT,14) HTCONDW, RHOW, HTCAPW
14     FORMAT(1X,'HEAT CONDUCTIVITY OF WATER          (HTCONDW) =',G15.6
     * /1X,     'DENSITY OF WATER                       (RHOW) =',G15.6
     * /1X,     'HEAT CAPACITY OF WATER                (HTCAPW)=',G15.6)
C------PREPARE HTCONDW TO BE DIVIDED BY RHOW AND BY HTCAPW
        HTCONDW = HTCONDW / RHOW / HTCAPW
        ENDIF
C ---------------------------------------------------------------------
C-------ALLOCATE SPACE FOR SOLUBILITY ARRAYS
        IF(ISOLUBILITY.EQ.1) THEN
          ALLOCATE(SOLLIM(MCOMP),SOLSLOPE(MCOMP))
        ELSEIF(ISOLUBILITY.EQ.2) THEN
          ALLOCATE(SOLLIMVAR(NEQS,MCOMP),SOLSLOPEVAR(NEQS,MCOMP))
        ENDIF
C ---------------------------------------------------------------------
C-------ALLOCATE SPACE FOR PARENT DECAY
      IF(ICHAIN.NE.0)THEN
        ALLOCATE (NPARENT(NTCOMP))
        ALLOCATE (JPARENT(NTCOMP,NTCOMP))
        IF (ISPRCT.EQ.0)THEN
          ALLOCATE (STOTIO(NTCOMP,NTCOMP))
        ELSE
          ALLOCATE (SPTLRCT(NEQS,NTCOMP,NTCOMP))
        ENDIF
        NPARENT = 0
        JPARENT = 0.0
        STOITIO = 0.0
      ENDIF
C ---------------------------------------------------------------------
C-------ALLOCATE SPACE FOR AIR-WATER ADSORPTION ARRAYS
        IF(IAW_ADSORB.EQ.1) THEN
          CALL AW_ADSORB1AL (IN) 
        ENDIF
C ---------------------------------------------------------------------      
C ----------------------------------------------------------------------
C ----ALLOCATE ARRAYS FOR TOTAL COMPONENT EQUATIONS INCLUDING HEAT
      ALLOCATE(CONC(NEQS,NTCOMP),CONCO(NEQS,NTCOMP), MASSBCT(NTCOMP)) ! USE LAST SPECIES OF CONCETRATION ARRAY FOR TEMPERATURE IF HEAT IS ON
      CONC = 0.0
C4A3-----ALLOCATE SPACE FOR VBVLT, AND VBNMT ARRAYS.
      ALLOCATE (VBVLT(4,NIUNIT,NTCOMP))
      ALLOCATE (VBNMT(NIUNIT,NTCOMP))
      ALLOCATE (MSUMT)
C4B------INITIALIZE ALLOCATED ARRAYS
      ADMAT = 0.0
      CBCF = 0.
      CBCH = 0.
      LrcC = 0
      Cncg = 0.0D0
      FLrcC = 0
      FCncg = 0.0D0
      IPCBFLAG = 0
      DO 600 ICOMP = 1,NTCOMP
      DO 600 I=1,NIUNIT
      DO 600 J=1,4
      VBVLT(J,I,ICOMP)=0.0
  600 CONTINUE
C
C5----READ PARAMETERS AND CONVERT FOR UNSTRUCTURED AND STRUCTURED GRIDS
      IF(IUNSTR.EQ.0) THEN
        CALL SGWT2BCT1S(IN,IC_IBOUND_FLG)
      ELSE
        CALL SGWT2BCT1G(IN,IC_IBOUND_FLG)
      ENDIF
C
C6----READ PARAMETERS FOR CONDUIT DOMAIN NODES
      IF(INCLN.GT.0) THEN
        CALL SGWT2BCT1CLN(INCLN,IC_IBOUND_FLG)
      ENDIF
C---------------------------------------------------------------------------
C7-----PREPARE PARAMETERS
C---------------------------------------------------------------------------
C---------------------------------------------------------------------------
C7A----SET CONCO TO CONC
      DO ICOMP = 1,NTCOMP
      DO N = 1,NEQS
        CONCO(N,ICOMP) = CONC(N,ICOMP)
      ENDDO
      ENDDO
C7A----FOR ZERO ORDER DECAY IN WATER
      IF(IZOD.EQ.1.OR.IZOD.EQ.3)THEN
        DO ICOMP = 1,MCOMP
        DO N = 1,NODES
          ZODRW(N,ICOMP) = ZODRW(N,ICOMP) * PRSITY(N)
        ENDDO
        ENDDO
      ENDIF
C7B----FOR ZERO ORDER DECAY IN SOIL
      IF((IZOD.EQ.2.OR.IZOD.EQ.3).AND.IADSORB.NE.0)THEN
        DO ICOMP = 1,MCOMP
        DO N = 1,NODES
          ZODRS(N,ICOMP) = ZODRS(N,ICOMP) * (1.0 - PRSITY(N))
        ENDDO
        ENDDO
      ENDIF
C7C----FOR FIRST ORDER DECAY IN WATER
      IF(IFOD.EQ.1.OR.IFOD.EQ.3)THEN
        DO ICOMP = 1,MCOMP
        DO N = 1,NODES
          FODRW(N,ICOMP) = FODRW(N,ICOMP) * PRSITY(N)
        ENDDO
        ENDDO
      ENDIF
C7D-----FOR CROSS DISPERSION LONGITUDINAL DISPERSIVITIES
      IF(IDISP.GT.0.AND.IXDISP.NE.0)THEN
        DO N = 1,NODES
          ALXY(N) = 2.0 * DLX(N) * DLY(N)/( DLX(N) + DLY(N))
          ALYZ(N) = 2.0 * DLY(N) * DLZ(N)/( DLY(N) + DLZ(N))
          ALXZ(N) = 2.0 * DLX(N) * DLZ(N)/( DLX(N) + DLZ(N))
        ENDDO
      ENDIF
C7E-----SET SECONDARY STORAGE TO POROSITY IF IT IS NOT READ
      ILAYCON13=0
      DO K=1,NLAY
        IF(LAYCON(K).EQ.0.OR.LAYCON(K).EQ.2)THEN
          NNDLAY = NODLAY(K)
          NSTRT = NODLAY(K-1)+1
          DO N=NSTRT,NNDLAY
            SC2(N) = PRSITY(N)*AREA(N)
          ENDDO
        ENDIF
      ENDDO
C7F-----FOR ITRNSP = 5, INITIALIZE SN AND SET FLOW TO STEADY-STATE
      IF(ITRNSP.EQ.5)THEN
        DO KPER = 1,NPER
          ISSFLG(KPER) = 1
        ENDDO
      ENDIF
C
C--------------------------------------------------------------------------------
C8----RETURN
      RETURN
      END
C-------------------------------------------------------------------------
      SUBROUTINE INITMASS (KPER,KSTP,ISTARTMASS)
C     ******************************************************************
C     COMPUTE THE INITIAL MASS OF MOBILE COMPONENTS
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL, ONLY:IOUT,ISSFLG,NODES,NEQS,So,TOP,BOT,AREA,IDPT,HOLD,
     *            Sn,BUFF,IUNSTR,NROW,NCOL,NLAY,IDPOUT,NODLAY
      USE GWTBCTMODULE
      USE GWFBASMODULE, ONLY: PERTIM,TOTIM,CSPCFM,ISPCFL,IOFLG,CSPCFM,
     1  ISPCFM,ISPCUN,LBHDSV 
      USE GWFDPFMODULE, ONLY: PHIF
      USE GWFBCFMODULE, ONLY: SC1,SC2
      USE AW_ADSORBMODULE, ONLY: AREA_AWI, AK_AWI 
      DOUBLE PRECISION,  SAVE,    DIMENSION(:),    ALLOCATABLE ::BUFF8  
      DOUBLE PRECISION MASLOCW,MASLOCS,MASLOCC,REFHD,MASLOPR,MASLOAWI,
     1  ETA,VOLU,ALENG,CW,SLIM,SSLOPE 
      CHARACTER*16 TEXT,ISTATUS,FORMM  
      CHARACTER*400 CFILENAME          
      INCLUDE 'openspec.inc'  
C      DATA TEXT /'    SOLUTE MASS'/
C     ------------------------------------------------------------------     
      ALLOCATE(BUFF8(NODES)) 
C1-----INITIALIZE MASS FOR EACH COMPONENT
      DO ICOMP = 1,NTCOMP
        MASSBCT (ICOMP) = 0.0
      ENDDO
      MASLOCW = 0.0
      MASLOCS = 0.0
      MASLOCC = 0.0
      MASLOPR = 0.0
      MASLOAWI = 0.0
      ISS=ISSFLG(1)
C1A-------INITIALIZE A-W SURFACE AREA AND K_AW       
      IF(IAW_ADSORB.EQ.1)THEN       
        CALL AWI_AREA
         DO ICOMP = 1,NTCOMP
          CALL AWI_GETKAWI(ICOMP)      
         ENDDO
      ENDIF  
C4---------FOR EACH MOBILE COMPONENT PLUS TEMPERATURE
      DO ICOMP = 1,MCOMPT      
       WRITE(IOUT,76) ICOMP
76     FORMAT (10X,'TOTAL MASS COMPUTATION FOR COMPONENT NUMBER',I6)    
       CALL GET_TEXTM(ICOMP,TEXT)
C2-----COMPUTE MASS FOR EACH ACTIVE NODE
        DO N = 1,NODES
          IF(ICBUND(N).EQ.0) THEN
            CONC(N,ICOMP) = CINACT
            BUFF8(N) = 0.0
            BUFF(N) = 0.0
          ELSE
C3---------COMPUTE VOLUME OF NODE
            ALENG = TOP(N) - BOT(N)
            VOLU = AREA(N) * ALENG
C5---------COMPUTE MASS IN WATER
            IF(ISS.EQ.1)THEN
              MASLOCW = PRSITY(N) * Sn(N) * CONC(N,ICOMP) * VOLU
            ELSE
C6-----------INITIAL MASS IN PORE STORAGE
              MASLOCW=(PRSITY(N)*VOLU-SC2(N)*ALENG+SC2(N)*ALENG*So(N))
     *         *CONC(N,ICOMP)
C7-------------INITIAL MASS IN COMPRESSIBLE STORAGE ASSUMED FROM TOP OF CELL
              REFHD = HOLD(N) - TOP(N)
              IF(REFHD.LT.0.0) REFHD = 0.0
              MASLOCC = SC1(N) * REFHD * CONC(N,ICOMP)
            ENDIF
C8-----------COMPUTE MASS ON SOIL
            IF(IADSORB.EQ.1)THEN      !LINEAR 
              MASLOCS =  ADSORB(N,ICOMP)*CONC(N,ICOMP) * VOLU
            ELSEIF(IADSORB.EQ.2) THEN !FREUNDLICH   
              ETA = FLICH(N,ICOMP)
              MASLOCS =  ADSORB(N,ICOMP)*CONC(N,ICOMP)**ETA * VOLU
            ELSEIF(IADSORB.EQ.3) THEN !LANGMUIR  
              ETA = FLICH(N,ICOMP)
              MASLOCS =  ADSORB(N,ICOMP)*CONC(N,ICOMP) * VOLU /
     1           (1.0 + ETA * CONC(N,ICOMP))
            ENDIF
C9-----------COMPUTE MASS ON AIR-WATER INTERFACE
            IF(IAW_ADSORB.EQ.1)THEN       
              MASLOAWI= AREA_AWI(N)* AK_AWI(N,ICOMP)*CONC(N,ICOMP)* VOLU
            ENDIF
C10-----------COMPUTE MASS PRECIPITATED ABOVE SOLUBILITY LIMIT
            MASLOPR = 0.0
            IF(ISOLUBILITY.NE.0)THEN
              CW = CONC(N,ICOMP)
              IF(ISOLUBILITY.EQ.1)THEN
                SLIM = SOLLIM(ICOMP)
                SSLOPE = SOLSLOPE(ICOMP)
              ELSEIF(ISOLUBILITY.EQ.2)THEN
                SLIM = SOLLIMVAR(N,ICOMP)
                SSLOPE = SOLSLOPEVAR(N,ICOMP)
              ENDIF    
              IF(CW.GT.SLIM) MASLOPR = (CW-SLIM) * SSLOPE
              MASLOPR =  MASLOPR * VOLU
            ENDIF
C11------------ADD UP ALL THE TERMS
            MASSBCT(ICOMP) = MASSBCT(ICOMP)+MASLOCW +MASLOCS+MASLOCC
     1        + MASLOPR + MASLOAWI 
            BUFF8(N) = MASLOCW +MASLOCS+MASLOCC
     1        + MASLOPR + MASLOAWI  
          ENDIF
C
C11B------PRINT MASS IF REQUESTED.
          IF(IMASSWR.LT.0. AND. ISPCFL. NE. 0 ) THEN
            IF(N.EQ.1) WRITE(IOUT,61) TEXT,KPER,KSTP
   61       FORMAT(1X,/1X,A,'   PERIOD ',I4,'   STEP ',I3)
            IF(IUNSTR.EQ.0.AND.N.LE.NODES)THEN
              IL = (N-1) / (NCOL*NROW) + 1
              IJ = N - (IL-1)*NCOL*NROW
              IR = (IJ-1)/NCOL + 1
              IC = IJ - (IR-1)*NCOL
               WRITE(IOUT,62) IL,IR,IC,BUFF8(N)
   62        FORMAT(1X,'   LAYER ',I5,'   ROW ',I6,'   COL ',I6,
     1       '   FLUX ',1PG15.6)
            ELSE
             WRITE(IOUT,63) N,BUFF8(N)
   63        FORMAT(1X,'    NODE ',I8,'   FLUX ',1PG15.6)
            ENDIF
          END IF
        ENDDO
C12------IF CELL-BY-CELL MASS TERMS WILL BE SAVED AS A 3-D ARRAY,
C12------CALL UBUDSV TO SAVE THEM.
        IF(IMASSWR.GE.1. AND. ISPCFL. NE. 0 )THEN           
C--------------------------------------------------------------------------------------------------------------------
C5A ----FOR MULTIFILE OUTPUT OPTION, ADD COMPONENT SPECIES TO ROOT NAME AND OPEN FILE ON ISPCFM AND GO TO END OF FILE
      IFIRST = 1
      IF(IMULTI.NE.0)THEN
        CFILENAME = trim(CROOTNAME) // TRIM(TEXT) //'.MAS'
        ISTATUS = 'OLD'
        FORMM = 'FORMATTED'
        IF(CSPCFM.EQ.' ') FORMM = FORM            ! - use FORM from openspec.inc rather than hardwiring BINARY or UNFORMATTED
CC        IF(KSTP.EQ.1.AND.KPER.EQ.1) THEN
C---------MAKE STATUS OF FILE AS NEW WHEN FIRST ENTERING FILE
CC          ISTATUS = 'REPLACE'
CC        ENDIF

!        OPEN(ISPCUN, FILE = CFILENAME, FORM = 'BINARY',   
        OPEN(IMASSWR, FILE = CFILENAME, FORM = FORMM,
     1  STATUS = ISTATUS, POSITION = 'APPEND')
      ENDIF
C--------------------------------------------------------------------------------------------------------------------              
          IF(IUNSTR.EQ.0)THEN
C---------FOR STRUCTURED GRID *** NOTE, STILL USING THE BUDGET FORMATS INSTEAD OF CONC FORMATS FOR BINARY OUTPUT FOR STRUCTURED GRID 
            DO N = 1,NODES 
               BUFF(N) = BUFF8(N) 
            ENDDO  
            CALL UBUDSV(KSTP,KPER,TEXT,IMASSWR,BUFF,NCOL,NROW,
     1                          NLAY,IOUT)
          ELSE
C---------FOR UNSTRUCTURED GRID              
            DO 79 K=1,NLAY
            KK=K
            IF(IOFLG(K,6).EQ.0) GO TO 79
            NNDLAY = NODLAY(K)
            NSTRT = NODLAY(K-1)+1
            IF(IFIRST.EQ.1) WRITE(IOUT,74) IMASSWR,KSTP,KPER
   74       FORMAT(1X,/1X,'MASS WILL BE SAVED ON UNIT ',I4,
     1      ' AT END OF TIME STEP ',I8,', STRESS PERIOD ',I8)
            IFIRST=0
            IF(IDPOUT.EQ.1)THEN
              WRITE(IOUT,*)'MASS AND TIME VARIABLES ARE SAVED AS REAL*8'
              CALL ULASAVU8(BUFF8,TEXT,KSTP,KPER,PERTIM,TOTIM,NSTRT,
     1                NNDLAY,KK,IMASSWR,NODES)
            ELSEIF(CSPCFM.EQ.' ') THEN
               DO N = 1,NODES 
               BUFF(N) = BUFF8(N) 
               ENDDO
               CALL ULASAVU(BUFF,TEXT,KSTP,KPER,PERTIM,TOTIM,NSTRT,
     1                NNDLAY,KK,IMASSWR,NODES)
            ELSE
               DO N = 1,NODES 
               BUFF(N) = BUFF8(N) 
               ENDDO
               CALL ULASV2U(BUFF,TEXT,KSTP,KPER,PERTIM,TOTIM,NSTRT,
     1             NNDLAY,KK,IMASSWR,CSPCFM,LBHDSV,ICBUND(NSTRT),NODES)
            END IF
   79       CONTINUE
          ENDIF
        ENDIF      
      ENDDO       
C13--------FOR IMMOBILE COMPONENTS COMPUTE MASS AS ONLY BEING IN THE TOTAL VOLUME
        IF(NTCOMP.GT.MCOMPT)THEN
          DO ICOMP = MCOMPT,NTCOMP
            DO N = 1,NODES
              MASSBCT(ICOMP) = MASSBCT(ICOMP)+CONC(N,ICOMP) * VOLU
            ENDDO
          ENDDO
      ENDIF          
C ---------------------------------------------------------------------      
C14----WRITE HEADER 
      IF(ISTARTMASS.EQ.1) THEN 
      WRITE(IOUT,12)
12    FORMAT(/5X,'INITIAL MASS OF COMPONENTS:', 
     *       /5X,'--------------------------:',
     *       /5X,'IN MOBILE DOMAIN'/,
     *        5X,'------------------'/,
     *       5X,'COMPONENT INDEX',5X,'MASS'/,
     *       5X,'---------------',5X,'----')
      ELSE 
        WRITE(IOUT,14) KPER,KSTP
14      FORMAT(/5X,'MASS OF COMPONENTS AT THE ',
     *           'END OF STRESS PERIOD',I10, ' AND TIME STEP',I10/
     *        5X,'------------------------------------'
     *          ,'-----------------','----------','--------------'
     *          ,'----------'/
     *       /5X,'IN MOBILE DOMAIN'/
     *       /5X,'------------------'/        
     *       5X,'COMPONENT INDEX',5X,'MASS'/
     *       5X,'---------------',5X,'----')          
      ENDIF    
C --------------------------------------------------------------------------------------   
C15-----WRITE MASS TO OUTPUT FILE
      DO ICOMP = 1,NTCOMP
        WRITE(IOUT,13) ICOMP, MASSBCT(ICOMP)
      ENDDO
13    FORMAT(I15,5X,E15.6)
C
C--------------------------------------------------------------------------------
      DEALLOCATE(BUFF8) 
C10----RETURN
      RETURN
      END
C
C-------------------------------------------------------------------------
      SUBROUTINE INITMASSCLN 
C     ******************************************************************
C     COMPUTE THE INITIAL MASS OF MOBILE COMPONENTS IN THE CLN DOMAIN
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL, ONLY:IOUT,ISSFLG,NODES,NEQS,Sn,So,AREA,HOLD
      USE GWTBCTMODULE
      USE CLN1MODULE, ONLY: NCLNNDS,ACLNNDS
      DOUBLE PRECISION MASLOCW,SAT,CW,SLIM,SSLOPE,ALENG,VOLU
C     ------------------------------------------------------------------
C1----COMPUTE INITIAL MASS FOR MOBILE COMPONENTS
      WRITE(IOUT,12)
12    FORMAT(/5X,'IN CLN DOMAIN'/
     *        5X,'-------------'/
     *       5X,'COMPONENT INDEX',5X,'MASS'/
     *       5X,'---------------',5X,'----')
C1-----INITIALIZE MASS FOR EACH COMPONENT
      DO ICOMP = 1,NTCOMP
        MASSBCT (ICOMP) = 0.0
      ENDDO
      MASLOCW = 0.0
      MASLOCS = 0.0
      MASLOCC = 0.0
      ISS=ISSFLG(1)
C2-----COMPUTE MASS FOR EACH ACTIVE CLN NODE
      DO N = NODES+1, NEQS
        IF(ICBUND(N).EQ.0) THEN
            DO ICOMP = 1,NTCOMP
              CONC(N,ICOMP) = CINACT
            ENDDO
        ELSE
C3---------COMPUTE VOLUME OF NODE
          ALENG = ACLNNDS(N-NODES,4)
          VOLU = AREA(N) * ALENG
          SAT = So(N)
          IF(ISS.EQ.1) SAT = Sn(N)
C4---------FOR EACH COMPONENT
          DO ICOMP = 1,MCOMPT
C5---------COMPUTE MASS IN WATER (NO COMPRESSIBLE STORAGE IN CLN
            MASLOCW = SAT * CONC(N,ICOMP) * VOLU
            MASSBCT(ICOMP) = MASSBCT(ICOMP)+MASLOCW
C6------------COMPUTE MASS PRECIPITATED ABOVE SOLUBILITY LIMIT
            MASLOPR = 0.0
            IF(ISOLUBILITY.NE.0)THEN
              CW = CONC(N,ICOMP)
              IF(ISOLUBILITY.EQ.1)THEN
                SLIM = SOLLIM(ICOMP)
                SSLOPE = SOLSLOPE(ICOMP)
              ELSEIF(ISOLUBILITY.EQ.2)THEN
                SLIM = SOLLIMVAR(N,ICOMP)
                SSLOPE = SOLSLOPEVAR(N,ICOMP)
              ENDIF
              IF(CW.GT.SLIM) MASLOPR = (CW-SLIM) * SSLOPE
              MASLOPR =  MASLOPR * VOLU
              MASSBCT(ICOMP) = MASSBCT(ICOMP)+MASLOPR
            ENDIF            
          ENDDO
        ENDIF
      ENDDO
C9--------FOR IMMOBILE COMPONENTS COMPUTE MASS AS ONLY BEING IN THE TOTAL VOLUME
      IF(NTCOMP.GT.MCOMPT)THEN
        DO ICOMP = MCOMPT+1,NTCOMP
          DO N = NODES+1,NEQS
            MASSBCT(ICOMP) = MASSBCT(ICOMP)+CONC(N,ICOMP) * VOLU
          ENDDO
        ENDDO
      ENDIF
C9-----WRITE MASS TO OUTPUT FILE
      DO ICOMP = 1,NTCOMP
        WRITE(IOUT,13) ICOMP, MASSBCT(ICOMP)
      ENDDO
13    FORMAT(I15,5X,E15.6)
C
C--------------------------------------------------------------------------------
C10----RETURN
      RETURN
      END
C
C--------------------------------------------------------------------------------
      SUBROUTINE SGWT2BCT1S(IN,IC_IBOUND_FLG)
C     ******************************************************************
C-----READ PARAMETERS AND CONVERT FOR STRUCTURED GRID
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,      ONLY:IOUT,NCOL,NROW,NLAY,ITRSS,LAYHDT,LAYHDS,
     1                      IFREFM,IBOUND,NODES,NJA,NJAS,IA,JA,JAS,ARAD,
     1                      AREA,IDPIN,IVC
      USE GWTBCTMODULE,ONLY: PRSITY,CONC,ICBUND,CINACT,MCOMP,IADSORB,
     1  ADSORB,FLICH,IDISP,DLX,ATXY,DLY,ATYZ,DLZ,ATXZ,IZOD,IFOD,
     1  ZODRW,FODRW,ZODRS,FODRS,HTCONDS,IHEAT,MCOMPT,NTCOMP,RHOW,HTCAPW,
     1  ICHAIN,NPARENT,JPARENT,STOTIO,ISPRCT,SPTLRCT,SOLLIMVAR,
     1  ISOLUBILITY,SOLLIM,SOLSLOPE,IAW_ADSORB,ZODRAW,FODRAW,SOLSLOPEVAR
      USE GWFBCFMODULE, ONLY: IHANISO,LAYCON,SC2
C
      INTEGER, DIMENSION(:,:),    ALLOCATABLE ::ITEMP
      REAL,    DIMENSION(:,:),    ALLOCATABLE ::TEMP
      REAL*8,  DIMENSION(:,:),    ALLOCATABLE ::TEMP8
      REAL,    DIMENSION(:),    ALLOCATABLE ::BULKD
      CHARACTER*24 ANAME(28),CNAME
      DATA ANAME(1) /'                  ICBUND'/
      DATA ANAME(2) /'                POROSITY'/
      DATA ANAME(3) /'            BULK DENSITY'/
      DATA ANAME(4) /'  ADSORPTION COEFFICIENT'/
      DATA ANAME(5) /'           CONCENTRATION'/
      DATA ANAME(6) /'       LONG DISPERSIVITY'/
      DATA ANAME(7) /'       TRNS DISPERSIVITY'/
      DATA ANAME(8) /' X-DIR LONG DISPERSIVITY'/
      DATA ANAME(9) /' Y-DIR LONG DISPERSIVITY'/
      DATA ANAME(10) /' Z-DIR LONG DISPERSIVITY'/
      DATA ANAME(11) /'XY-DIR TRNS DISPERSIVITY'/
      DATA ANAME(12) /'YZ-DIR TRNS DISPERSIVITY'/
      DATA ANAME(13) /'XZ-DIR TRNS DISPERSIVITY'/
      DATA ANAME(14) /'     FREUNDLICH EXPONENT'/
      DATA ANAME(15) /'  ZERO-ORDER DECAY WATER'/
      DATA ANAME(16) /'   ZERO-ORDER DECAY SOIL'/
      DATA ANAME(17) /' FIRST-ORDER DECAY WATER'/
      DATA ANAME(18) /'  FIRST-ORDER DECAY SOIL'/
      DATA ANAME(19) /'     FREUNDLICH EXPONENT'/
      DATA ANAME(20) /'      SOIL HEAT CAPACITY'/
      DATA ANAME(21) /'  SOIL HEAT CONDUCTIVITY'/
      DATA ANAME(22) /'             TEMPERATURE'/
      DATA ANAME(23) /'PARENTS OF SPECIES ICOMP'/
      DATA ANAME(24) /'STOTIOMETRIC COEFFICENTS'/
      DATA ANAME(25) /'        SOLUBILITY LIMIT'/
      DATA ANAME(26) /'        SOLUBILITY SLOPE'/
      DATA ANAME(27) /'    ZERO-ORDER DECAY A-W'/       
      DATA ANAME(28) /'   FIRST-ORDER DECAY A-W'/              
      REAL PI
C     ------------------------------------------------------------------
C1-------ALLOCATE TEMP ARRAY FOR STORING 3-D INFORMATION
      ALLOCATE(TEMP(NCOL,NROW))
      ALLOCATE(TEMP8(NCOL,NROW))
      IF(IADSORB.NE.0.OR.IHEAT.EQ.1) ALLOCATE(BULKD(NODES))
      PI = 3.1415926536
C-----------------------------------------------------------------------
C2------READ/SET SPECIES INDEPENDENT ARRAYS FOR ALL LAYERS.
C---------------------------------------------------------
C3-----READ/SET ICBUND ARRAY
      IF(IC_IBOUND_FLG.EQ.0)THEN
C3A-------READ ICBUND
        ALLOCATE(ITEMP(NCOL,NROW))
        DO K=1,NLAY
        KK=K
        CALL U2DINT(ITEMP(1,1),ANAME(1),NROW,NCOL,KK,IN,IOUT)
        DO I=1,NROW
        DO J=1,NCOL
          N=J+(I-1)*NCOL+(K-1)*NROW*NCOL
          ICBUND(N) = ITEMP(J,I)
        ENDDO
        ENDDO
        ENDDO
        DEALLOCATE (ITEMP)
C3B-------CHECK FOR ERRORS
         DO N=1,NODES
           IF(ICBUND(N).NE.0.AND.IBOUND(N).EQ.0)THEN
             WRITE(IOUT,55) N
55    FORMAT(5X,'*** ACTIVE TRANSPORT NODE WHERE FLOW WAS INACTIVE AT ',
     *  'NODE: ',I10,' ***')
             STOP
           ENDIF
         ENDDO
      ELSE
C3C-------SET ICBUND
        DO N=1,NODES
          ICBUND(N) = IBOUND(N)
        ENDDO
      ENDIF
C
C-----------------------------------------------------------------------------
C4-----READ POROSITY INTO ARRAY PRSITY.
      DO K=1,NLAY
        KK=K
        CALL U2DREL(TEMP(1,1),ANAME(2),NROW,NCOL,KK,IN,IOUT)
        DO I=1,NROW
        DO J=1,NCOL
          N=J+(I-1)*NCOL+(K-1)*NROW*NCOL
          PRSITY(N) = TEMP(J,I)
C4B-------CHECK FOR ERRORS, POROSITY CANNOT BE LESS THAN SPECIFIC YIELD
          IF(LAYCON(K).NE.3.AND.LAYCON(K).NE.2.AND.LAYCON(K).NE.4
     1       .AND.LAYCON(K).NE.5)GO TO 90
          IF(ICBUND(N).EQ.0.OR.IBOUND(N).EQ.0)GO TO 90
          IF(PRSITY(N)+1.e-5.LT.SC2(N)/AREA(N)) THEN
            WRITE(IOUT,56) N
56          FORMAT(5X,'*** Porosity < Sy FOR NODE: ',I10,' ***')
            STOP
          ENDIF
90        CONTINUE
        ENDDO
        ENDDO
      ENDDO
C-----------------------------------------------------------------------------
C5-----READ BULK DENSITY INTO ARRAY BULKD.
      IF(IADSORB.NE.0)THEN
        DO K=1,NLAY
          KK=K
          CALL U2DREL(TEMP(1,1),ANAME(3),NROW,NCOL,KK,IN,IOUT)
          DO I=1,NROW
          DO J=1,NCOL
            N=J+(I-1)*NCOL+(K-1)*NROW*NCOL
            BULKD(N) = TEMP(J,I)
          ENDDO
          ENDDO
        ENDDO
      ENDIF
C----------------------------------------------------------------------------
C6-----READ DISPERSION COEFFICIENTS AND FACE ANGLES IF DISPERSION IS ON.
      IF(IDISP.NE.0)THEN
        IF(IHANISO.EQ.0)THEN
C6A-------ARAD HAS NOT BEEN ALLOCATED YET SO ALLOCATE IT
          ALLOCATE(ARAD(NJAS))
        ENDIF
C6B------SET FACE ANGLES IN ARAD
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
C7-----READ DISPERSION COEFFICIENTS
       IF(IDISP.EQ.1)THEN
C8--------READ ISOTROPIC DISPERSION COEFFICIENTS
         DO K=1,NLAY
           KK=K
           CALL U2DREL(TEMP(1,1),ANAME(6),NROW,NCOL,KK,IN,IOUT)
           DO I=1,NROW
           DO J=1,NCOL
             N=J+(I-1)*NCOL+(K-1)*NROW*NCOL
             DLX(N) = TEMP(J,I)
             DLY(N) = TEMP(J,I)
             DLZ(N) = TEMP(J,I)
           ENDDO
           ENDDO
         ENDDO
C
         DO K=1,NLAY
           KK=K
           CALL U2DREL(TEMP(1,1),ANAME(7),NROW,NCOL,KK,IN,IOUT)
           DO I=1,NROW
           DO J=1,NCOL
             N=J+(I-1)*NCOL+(K-1)*NROW*NCOL
             ATXY(N) = TEMP(J,I)
             ATYZ(N) = TEMP(J,I)
             ATXZ(N) = TEMP(J,I)
           ENDDO
           ENDDO
         ENDDO
       ELSEIF(IDISP.EQ.2)THEN
C9--------READ ANISOTROPIC DISPERSION COEFFICIENTS
         DO K=1,NLAY
           KK=K
           CALL U2DREL(TEMP(:,:),ANAME(8),NROW,NCOL,KK,IN,IOUT)
           DO I=1,NROW
           DO J=1,NCOL
             N=J+(I-1)*NCOL+(K-1)*NROW*NCOL
             DLX(N) = TEMP(J,I)
           ENDDO
           ENDDO
         ENDDO
C
         DO K=1,NLAY
           KK=K
           CALL U2DREL(TEMP(1,1),ANAME(9),NROW,NCOL,KK,IN,IOUT)
           DO I=1,NROW
           DO J=1,NCOL
             N=J+(I-1)*NCOL+(K-1)*NROW*NCOL
             DLY(N) = TEMP(J,I)
           ENDDO
           ENDDO
         ENDDO
C
         DO K=1,NLAY
           KK=K
           CALL U2DREL(TEMP(1,1),ANAME(10),NROW,NCOL,KK,IN,IOUT)
           DO I=1,NROW
           DO J=1,NCOL
             N=J+(I-1)*NCOL+(K-1)*NROW*NCOL
             DLZ(N) = TEMP(J,I)
           ENDDO
           ENDDO
         ENDDO
C
         DO K=1,NLAY
           KK=K
           CALL U2DREL(TEMP(1,1),ANAME(11),NROW,NCOL,KK,IN,IOUT)
           DO I=1,NROW
           DO J=1,NCOL
             N=J+(I-1)*NCOL+(K-1)*NROW*NCOL
             ATXY(N) = TEMP(J,I)
           ENDDO
           ENDDO
         ENDDO
C
         DO K=1,NLAY
           KK=K
           CALL U2DREL(TEMP(1,1),ANAME(12),NROW,NCOL,KK,IN,IOUT)
           DO I=1,NROW
           DO J=1,NCOL
             N=J+(I-1)*NCOL+(K-1)*NROW*NCOL
             ATYZ(N) = TEMP(J,I)
           ENDDO
           ENDDO
         ENDDO
C
         DO K=1,NLAY
           KK=K
           CALL U2DREL(TEMP(1,1),ANAME(13),NROW,NCOL,KK,IN,IOUT)
           DO I=1,NROW
           DO J=1,NCOL
             N=J+(I-1)*NCOL+(K-1)*NROW*NCOL
             ATXZ(N) = TEMP(J,I)
           ENDDO
           ENDDO
         ENDDO
       ENDIF
      ENDIF
C-----------------------------------------------------------------------
C10-----READ HEAT PARAMETERS IF IHEAT IS ON
      IF(IHEAT.EQ.1) THEN
C10A----READ AND PREPARE THE SOIL HEAT CAPACITY TERM
        DO K=1,NLAY
          KK=K
          CALL U2DREL(TEMP(1,1),ANAME(20),NROW,NCOL,KK,IN,IOUT)
          DO I=1,NROW
          DO J=1,NCOL
            N=J+(I-1)*NCOL+(K-1)*NROW*NCOL
            ADSORB(N,MCOMPT) = TEMP(J,I)
          ENDDO
          ENDDO
        ENDDO
C10B-------PREPARE CONSTANT PART OF HEAT CAPACITY TERM
        DO N = 1, NODES
          ADSORB(N,MCOMPT) = ADSORB(N,MCOMPT) /RHOW / HTCAPW * BULKD(N)
        ENDDO
C10C----READ AND PREPARE THE SOIL HEAT CONDUCTIVITY TERM
        DO K=1,NLAY
          KK=K
          CALL U2DREL(TEMP(1,1),ANAME(21),NROW,NCOL,KK,IN,IOUT)
          DO I=1,NROW
          DO J=1,NCOL
            N=J+(I-1)*NCOL+(K-1)*NROW*NCOL
            HTCONDS(N) = TEMP(J,I)
          ENDDO
          ENDDO
        ENDDO
C10D-------PREPARE CONSTANT PART OF SOIL HEAT CONDUCTIVITY TERM
        DO N = 1, NODES
          HTCONDS(N) = HTCONDS(N) * (1.0 - PRSITY(N)) /RHOW / HTCAPW
        ENDDO
      ENDIF
C-----------------------------------------------------------------------
C10E-----READ AIR-WATER INTERFACE ADSORPTION PARAMETERS IF AW_ADSORB OPTION
      IF(IAW_ADSORB.EQ.1) THEN 
        CALL AW_ADSORB1RP1  (IN)
      ENDIF    
C---------------------------------------------------------------------          
C10F----READ AND PREPARE SOLUBILITY PARAMETERS
      IF(ISOLUBILITY .EQ.1) THEN
        CALL U2DREL(SOLLIM(1),ANAME(25),1,MCOMP,1,IN,IOUT)
        CALL U2DREL(SOLSLOPE(1),ANAME(26),1,MCOMP,1,IN,IOUT)
      ENDIF
C-----------------------------------------------------------------------
C11------READ SPECIES DEPENDENT ARRAYS FOR ALL LAYERS.
C
      DO ICOMP=1,NTCOMP
        WRITE(IOUT,10) ICOMP
10      FORMAT(80('-')/1X,'THE FOLLOWING ARRAYS ARE READ FOR COMPONENT',
     *    1X,'NUMBER ',I3)
        IF(IHEAT.EQ.1.AND.ICOMP.EQ.MCOMPT) GO TO 100 !TEMPERATURE DOES NOT HAVE ADSORPTION, DECAY OR PARENT
        IF(ICOMP .GT. MCOMPT) GO TO 100     !IMMOBILE COMPONENTS DO NOT HAVE ADSORPTION, DECAY OR PARENT
C
C---------------------------------------------------------
C11AFIRST---IF CHAIN DECAY THEN READ PARENT DETAILS FOR ALL MOBILE SPECIES
        IF(ICHAIN.NE.0)THEN
          READ(IN,*) NPARENT(ICOMP)
          NPAREN = NPARENT(ICOMP)
          WRITE(IOUT,12) NPARENT(ICOMP)
12        FORMAT(1X,'FOR CHAIN DECAY, COMPONENT NUMBER',I5,' HAS',I5,
     1     ' PARENT COMPONENTS')
          K=1
          IF(NPAREN.GT.0)THEN
            CALL U1DINT(ITEMP(1,1),ANAME(23),NPAREN,K,IN,IOUT)
            DO IPAREN = 1,NPAREN
              JPARENT(ICOMP,IPAREN) = ITEMP(IPAREN,1)
            ENDDO             
C
            IF (ISPRCT.EQ.0) THEN
              WRITE(IOUT,13) ICOMP, NPAREN
13            FORMAT(1X,'*** GENERATION FACTOR OF DAUGHTER SPECIES ',I5,
     1           ' FOR EACH OF THE ',I5,' PARENT SPECIES ***')
              CALL U1DREL(TEMP(1,1),ANAME(24),NPAREN,K,IN,IOUT)
              DO IPAREN = 1,NPAREN
                STOTIO(ICOMP,IPAREN) = TEMP(NPAREN,1)  
              ENDDO    
            ELSE
              DO JCOMP = 1,NPAREN
                WRITE(IOUT,14)  ICOMP, JPARENT(ICOMP,JCOMP)
14              FORMAT(1X,'*** GENERATION FACTOR OF DAUGHTER SPECIES ',
     1          I5, 'FROM PARENT SPECIES ',I5,' FOR EACH LAYER ***')
                DO K=1,NLAY
                  KK=K
                  CALL U2DREL(TEMP(1,1),ANAME(4),NROW,NCOL,KK,IN,IOUT)
                  DO I=1,NROW
                  DO J=1,NCOL
                    N=J+(I-1)*NCOL+(K-1)*NROW*NCOL
                    SPTLRCT(N,ICOMP,JCOMP) = TEMP(J,I)
                  ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ENDIF
          ENDIF
        ENDIF
C
C---------------------------------------------------------
C11ASECOND---IF AW_ADSORB OPTION IS ON THEN READ AW-PARTITIONING DATA
        IF(IAW_ADSORB.NE.0)THEN
          CALL AW_ADSORB1RP2(IN,ICOMP)   
        ENDIF    
C----------------------------------------------------------
C11A---READ ADSORPTION COEFFICIENT OF EACH SPECIES
      IF(IADSORB.NE.0)THEN
        DO K=1,NLAY
          KK=K
          CALL U2DREL(TEMP(1,1),ANAME(4),NROW,NCOL,KK,IN,IOUT)
          DO I=1,NROW
          DO J=1,NCOL
            N=J+(I-1)*NCOL+(K-1)*NROW*NCOL
            ADSORB(N,ICOMP) = TEMP(J,I)
          ENDDO
          ENDDO
        ENDDO
C11B-------PREPARE CONSTANT PART OF ADSORPTION TERM
        DO N = 1, NODES
          ADSORB(N,ICOMP) = ADSORB(N,ICOMP) * BULKD(N)
        ENDDO
C---------------------------------------------------------
C11C-----READ FREUNDLICH EXPONENT OR LANGMUIR COEFFICIENT OF EACH SPECIES 
        IF(IADSORB.EQ.2. OR.IADSORB.EQ.3)THEN
          DO K=1,NLAY
            KK=K
            CALL U2DREL(TEMP(1,1),ANAME(14),NROW,NCOL,KK,IN,IOUT)
            DO I=1,NROW
            DO J=1,NCOL
              N=J+(I-1)*NCOL+(K-1)*NROW*NCOL
              FLICH(N,ICOMP) = TEMP(J,I)
            ENDDO
            ENDDO
          ENDDO
        ENDIF
      ENDIF
C---------------------------------------------------------
C12------READ SOLUBILITY PARAMETERS FOR VARIABLE SOLUBILITY IN SPACE
      IF(ISOLUBILITY .EQ.2) THEN
        DO K=1,NLAY
          KK=K
          CALL U2DREL(TEMP(1,1),ANAME(25),NROW,NCOL,KK,IN,IOUT)
          DO I=1,NROW
          DO J=1,NCOL
            N=J+(I-1)*NCOL+(K-1)*NROW*NCOL
            SOLLIMVAR(N,ICOMP) = TEMP(J,I)
          ENDDO
          ENDDO
        ENDDO          
C         
        DO K=1,NLAY
          KK=K
          CALL U2DREL(TEMP(1,1),ANAME(26),NROW,NCOL,KK,IN,IOUT)
          DO I=1,NROW
          DO J=1,NCOL
            N=J+(I-1)*NCOL+(K-1)*NROW*NCOL
            SOLSLOPEVAR(N,ICOMP) = TEMP(J,I)
          ENDDO
          ENDDO
        ENDDO          
      ENDIF      
C---------------------------------------------------------
C13-----READ ZERO ORDER DECAY COEFFICIENTS
      IF(IZOD.EQ.1)THEN
          DO K=1,NLAY
            KK=K
            CALL U2DREL(TEMP(1,1),ANAME(15),NROW,NCOL,KK,IN,IOUT)
            DO I=1,NROW
            DO J=1,NCOL
              N=J+(I-1)*NCOL+(K-1)*NROW*NCOL
              ZODRW(N,ICOMP) = TEMP(J,I)
            ENDDO
            ENDDO
          ENDDO
      ELSE
        IF(IADSORB.GT.0)THEN
          IF(IZOD.EQ.2)THEN
            DO K=1,NLAY
              KK=K
              CALL U2DREL(TEMP(1,1),ANAME(16),NROW,NCOL,KK,IN,IOUT)
              DO I=1,NROW
              DO J=1,NCOL
                N=J+(I-1)*NCOL+(K-1)*NROW*NCOL
                ZODRS(N,ICOMP) = TEMP(J,I)
              ENDDO
              ENDDO
            ENDDO
          ELSEIF(IZOD.EQ.3)THEN
            DO K=1,NLAY
              KK=K
              CALL U2DREL(TEMP(1,1),ANAME(15),NROW,NCOL,KK,IN,IOUT)
              DO I=1,NROW
              DO J=1,NCOL
                N=J+(I-1)*NCOL+(K-1)*NROW*NCOL
                ZODRW(N,ICOMP) = TEMP(J,I)
              ENDDO
              ENDDO
            ENDDO
            DO K=1,NLAY
              KK=K
              CALL U2DREL(TEMP(1,1),ANAME(16),NROW,NCOL,KK,IN,IOUT)
              DO I=1,NROW
              DO J=1,NCOL
                N=J+(I-1)*NCOL+(K-1)*NROW*NCOL
                ZODRS(N,ICOMP) = TEMP(J,I)
              ENDDO
              ENDDO
            ENDDO
          ENDIF
        ENDIF 
C13A-----READ ZERO ORDER DECAY COEFFICIENTS FOR A-W INTERFACE         
        IF(IAW_ADSORB.NE.0)THEN
          IF(IZOD.EQ.4)THEN
            DO K=1,NLAY
              KK=K
              CALL U2DREL(TEMP(1,1),ANAME(15),NROW,NCOL,KK,IN,IOUT)
              DO I=1,NROW
              DO J=1,NCOL
                N=J+(I-1)*NCOL+(K-1)*NROW*NCOL
                ZODRW(N,ICOMP) = TEMP(J,I)
              ENDDO
              ENDDO
            ENDDO
            DO K=1,NLAY
              KK=K
              CALL U2DREL(TEMP(1,1),ANAME(16),NROW,NCOL,KK,IN,IOUT)
              DO I=1,NROW
              DO J=1,NCOL
                N=J+(I-1)*NCOL+(K-1)*NROW*NCOL
                ZODRS(N,ICOMP) = TEMP(J,I)
              ENDDO
              ENDDO
            ENDDO
            DO K=1,NLAY
              KK=K
              CALL U2DREL(TEMP(1,1),ANAME(27),NROW,NCOL,KK,IN,IOUT)
              DO I=1,NROW
              DO J=1,NCOL
                N=J+(I-1)*NCOL+(K-1)*NROW*NCOL
                ZODRAW(N,ICOMP) = TEMP(J,I)
              ENDDO
              ENDDO
            ENDDO
          ENDIF  
        ENDIF    
      ENDIF
C---------------------------------------------------------
C14-----READ FIRST ORDER DECAY COEFFICIENTS
      IF(IFOD.EQ.1)THEN
          DO K=1,NLAY
            KK=K
            CALL U2DREL(TEMP(1,1),ANAME(17),NROW,NCOL,KK,IN,IOUT)
            DO I=1,NROW
            DO J=1,NCOL
              N=J+(I-1)*NCOL+(K-1)*NROW*NCOL
              FODRW(N,ICOMP) = TEMP(J,I)
            ENDDO
            ENDDO
          ENDDO
      ELSE
        IF(IADSORB.GT.0)THEN
          IF(IFOD.EQ.2)THEN
            DO K=1,NLAY
              KK=K
              CALL U2DREL(TEMP(1,1),ANAME(18),NROW,NCOL,KK,IN,IOUT)
              DO I=1,NROW
              DO J=1,NCOL
                N=J+(I-1)*NCOL+(K-1)*NROW*NCOL
                FODRS(N,ICOMP) = TEMP(J,I)
              ENDDO
              ENDDO
            ENDDO
          ELSEIF(IFOD.EQ.3)THEN
            DO K=1,NLAY
              KK=K
              CALL U2DREL(TEMP(1,1),ANAME(17),NROW,NCOL,KK,IN,IOUT)
              DO I=1,NROW
              DO J=1,NCOL
                N=J+(I-1)*NCOL+(K-1)*NROW*NCOL
                FODRW(N,ICOMP) = TEMP(J,I)
            ENDDO
            ENDDO
            ENDDO
            DO K=1,NLAY
              KK=K
              CALL U2DREL(TEMP(1,1),ANAME(18),NROW,NCOL,KK,IN,IOUT)
              DO I=1,NROW
              DO J=1,NCOL
                N=J+(I-1)*NCOL+(K-1)*NROW*NCOL
                FODRS(N,ICOMP) = TEMP(J,I)
              ENDDO
              ENDDO
            ENDDO
          ENDIF
        ENDIF 
C14A-----READ FIRST ORDER DECAY COEFFICIENTS FOR A-W INTERFACE         
        IF(IAW_ADSORB.NE.0)THEN
          IF(IFOD.EQ.4)THEN
            DO K=1,NLAY
              KK=K
              CALL U2DREL(TEMP(1,1),ANAME(17),NROW,NCOL,KK,IN,IOUT)
              DO I=1,NROW
              DO J=1,NCOL
                N=J+(I-1)*NCOL+(K-1)*NROW*NCOL
                FODRW(N,ICOMP) = TEMP(J,I)
              ENDDO
              ENDDO
            ENDDO
            DO K=1,NLAY
              KK=K
              CALL U2DREL(TEMP(1,1),ANAME(18),NROW,NCOL,KK,IN,IOUT)
              DO I=1,NROW
              DO J=1,NCOL
                N=J+(I-1)*NCOL+(K-1)*NROW*NCOL
                FODRS(N,ICOMP) = TEMP(J,I)
              ENDDO
              ENDDO
            ENDDO
            DO K=1,NLAY
              KK=K
              CALL U2DREL(TEMP(1,1),ANAME(28),NROW,NCOL,KK,IN,IOUT)
              DO I=1,NROW
              DO J=1,NCOL
                N=J+(I-1)*NCOL+(K-1)*NROW*NCOL
                FODRAW(N,ICOMP) = TEMP(J,I)
              ENDDO
              ENDDO
            ENDDO
          ENDIF  
        ENDIF 
      ENDIF
C
C---------------------------------------------------------
C15-----READ CONCENTRATION OF EACH SPECIES
100   CONTINUE
      CNAME = ANAME(5)
      IF(IHEAT.EQ.1.AND. ICOMP.EQ.MCOMPT) CNAME = ANAME(22)
      IF(IDPIN.EQ.0) THEN !----------------------------SINGLE PRECISION READ
      DO K=1,NLAY
        KK=K
        CALL U2DREL(TEMP(1,1),CNAME,NROW,NCOL,KK,IN,IOUT)
        DO I=1,NROW
        DO J=1,NCOL
          N=J+(I-1)*NCOL+(K-1)*NROW*NCOL
          CONC(N,ICOMP) = TEMP(J,I)
        ENDDO
        ENDDO
      ENDDO
      ELSE    !----------------------------DOUBLE PRECISION READ
      DO K=1,NLAY
        KK=K
        CALL U2DREL8(TEMP8(1,1),CNAME,NROW,NCOL,KK,IN,IOUT)
        DO I=1,NROW
        DO J=1,NCOL
          N=J+(I-1)*NCOL+(K-1)*NROW*NCOL
          CONC(N,ICOMP) = TEMP8(J,I)
        ENDDO
        ENDDO
      ENDDO          
      ENDIF
C
C---------------------------------------------------------
      ENDDO  ! SPECIES DO LOOP
C---------------------------------------------------------
      IF(IADSORB.NE.0.OR.IHEAT.EQ.1) DEALLOCATE(BULKD)
      DEALLOCATE(TEMP)
      DEALLOCATE(TEMP8)
C---------------------------------------------------------
C
      RETURN
      END
C----------------------------------------------------------------------------
      SUBROUTINE SGWT2BCT1G(IN,IC_IBOUND_FLG)
C     ******************************************************************
C-----READ PARAMETERS AND CONVERT FOR UNSTRUCTURED (GENERALIZED) GRID
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,     ONLY:IOUT,NCOL,NROW,NLAY,ITRSS,IA,JA,JAS,NJA,NJAG,
     1                     ARAD,IFREFM,IBOUND,NODES,NODLAY,IDSYMRD,NJAS,
     2                     IATMP,NJATMP,AREA,IDPIN,IOUTNORMAL
      USE GWTBCTMODULE,ONLY: PRSITY,CONC,ICBUND,CINACT,MCOMP,IADSORB,
     1  ADSORB,FLICH,IDISP,DLX,ATXY,DLY,ATYZ,DLZ,ATXZ,IZOD,IFOD,
     1  ZODRW,FODRW,ZODRS,FODRS,MCOMPT,HTCONDS,IHEAT,NTCOMP,RHOW,HTCAPW,
     1  ICHAIN,NPARENT,JPARENT,STOTIO,ISPRCT,SPTLRCT,SOLSLOPEVAR,
     1   ISOLUBILITY,SOLLIM,SOLSLOPE,IAW_ADSORB,ZODRAW,FODRAW,SOLLIMVAR
      USE GWFBCFMODULE, ONLY: IHANISO,LAYCON,SC2
C
      REAL,  DIMENSION(:),    ALLOCATABLE ::BULKD
      REAL,  DIMENSION(:),    ALLOCATABLE ::TEMPC
      REAL*8,  DIMENSION(:),    ALLOCATABLE ::TEMPC8
      INTEGER, DIMENSION(:), ALLOCATABLE:: ITEMPC
      CHARACTER*24 ANAME(28),CNAME
      DATA ANAME(1) /'                  ICBUND'/
      DATA ANAME(2) /'                POROSITY'/
      DATA ANAME(3) /'            BULK DENSITY'/
      DATA ANAME(4) /'  ADSORPTION COEFFICIENT'/
      DATA ANAME(5) /'           CONCENTRATION'/
      DATA ANAME(6) /'              FACE ANGLE'/
      DATA ANAME(7) /'       LONG DISPERSIVITY'/
      DATA ANAME(8) /'       TRNS DISPERSIVITY'/
      DATA ANAME(9) /' X-DIR LONG DISPERSIVITY'/
      DATA ANAME(10) /' Y-DIR LONG DISPERSIVITY'/
      DATA ANAME(11) /' Z-DIR LONG DISPERSIVITY'/
      DATA ANAME(12) /'XY-DIR TRNS DISPERSIVITY'/
      DATA ANAME(13) /'YZ-DIR TRNS DISPERSIVITY'/
      DATA ANAME(14) /'XZ-DIR TRNS DISPERSIVITY'/
      DATA ANAME(15) /'  ZERO-ORDER DECAY WATER'/
      DATA ANAME(16) /'   ZERO-ORDER DECAY SOIL'/
      DATA ANAME(17) /' FIRST-ORDER DECAY WATER'/
      DATA ANAME(18) /'  FIRST-ORDER DECAY SOIL'/
      DATA ANAME(19) /'     FREUNDLICH EXPONENT'/
      DATA ANAME(20) /'      SOIL HEAT CAPACITY'/
      DATA ANAME(21) /'  SOIL HEAT CONDUCTIVITY'/
      DATA ANAME(22) /'             TEMPERATURE'/
      DATA ANAME(23) /'PARENTS OF SPECIES ICOMP'/
      DATA ANAME(24) /'STOTIOMETRIC COEFFICENTS'/
      DATA ANAME(25) /'        SOLUBILITY LIMIT'/
      DATA ANAME(26) /'        SOLUBILITY SLOPE'/
      DATA ANAME(28) /'   FIRST-ORDER DECAY A-W'/      
      DATA ANAME(27) /'    ZERO-ORDER DECAY A-W'/       
      DATA ANAME(28) /'   FIRST-ORDER DECAY A-W'/        
C     ------------------------------------------------------------------
      ALLOCATE (TEMPC(NODES))
       ALLOCATE (TEMPC8(NODES))
       ALLOCATE(ITEMPC(NODES))
C1-------ALLOCATE TEMP ARRAY FOR BULK DENSITY
      IF(IADSORB.NE.0.OR.IHEAT.EQ.1) ALLOCATE(BULKD(NODES))
C-----------------------------------------------------------------------
C2------READ/SET SPECIES INDEPENDENT ARRAYS FOR ALL LAYERS.
C-----------------------------------------------------------------------
C3-----READ/SET ICBUND ARRAY
      IF(IC_IBOUND_FLG.EQ.0)THEN
C3A-------READ ICBUND
        DO K=1,NLAY
          NNDLAY = NODLAY(K)
          NSTRT = NODLAY(K-1)+1
          NDSLAY = NNDLAY - NODLAY(K-1)
          CALL U1DINT(ICBUND(NSTRT),ANAME(1),NDSLAY,K,IN,IOUT)
        ENDDO
C3B-------CHECK FOR ERRORS
         DO N=1,NODES
           IF(ICBUND(N).NE.0.AND.IBOUND(N).EQ.0)THEN
             WRITE(IOUT,55) N
55    FORMAT(5X,'*** ACTIVE TRANSPORT NODE WHERE FLOW WAS INACTIVE AT ',
     *  'NODE: ',I10,' ***')
             STOP
           ENDIF
         ENDDO
      ELSE
C3C-------SET ICBUND
        DO N=1,NODES
          ICBUND(N) = IBOUND(N)
        ENDDO
      ENDIF
C
C---------------------------------------------------------
C4-----READ POROSITY INTO ARRAY PRSITY.
      DO K = 1,NLAY
        NNDLAY = NODLAY(K)
        NSTRT = NODLAY(K-1)+1
        NDSLAY = NNDLAY - NODLAY(K-1)
        CALL U1DREL(PRSITY(NSTRT),ANAME(2),NDSLAY,K,IN,IOUT)
C4B-----CHECK FOR ERRORS, POROSITY CANNOT BE LESS THAN SPECIFIC YIELD
        DO N = NSTRT,NNDLAY
          IF(LAYCON(K).NE.3.AND.LAYCON(K).NE.2.AND.LAYCON(K).NE.4
     1       .AND.LAYCON(K).NE.5)GO TO 90
          IF(ICBUND(N).EQ.0.OR.IBOUND(N).EQ.0)GO TO 90
          IF(PRSITY(N)+1.e-5.LT.SC2(N)/AREA(N)) THEN
            WRITE(IOUT,56) N
56          FORMAT(5X,'*** Porosity < Sy FOR NODE: ',I10,' ***')
            STOP
          ENDIF
90        CONTINUE
        ENDDO
      ENDDO
C---------------------------------------------------------
C5-----READ BULK DENSITY INTO ARRAY TEMP IF ADSORPTION.
      IF(IADSORB.NE.0.OR.IHEAT.EQ.1)THEN
        DO K = 1,NLAY
          NNDLAY = NODLAY(K)
          NSTRT = NODLAY(K-1)+1
          NDSLAY = NNDLAY - NODLAY(K-1)
          CALL U1DREL(BULKD(NSTRT),ANAME(3),NDSLAY,K,IN,IOUT)
        ENDDO
      ENDIF
C-----------------------------------------------------------------------
C6-----READ DISPERSION COEFFICIENTS AND FACE ANGLES IF DISPERSION IS ON.
      IF(IDISP.NE.0)THEN
C7------READ FACE ANGLES AND STORE IN ARAD
        IF(IHANISO.EQ.0)THEN
C8--------ARAD HAS NOT BEEN ALLOCATED YET SO ALLOCATE IT
          ALLOCATE(ARAD(NJAS))
        ENDIF
        CALL U1DRELNJA(ARAD,IATMP,ANAME(6),NJATMP,IN,IOUT,IDSYMRD)
C---------
C8A.......CHANGE OUTWARD NORMAL FOR COMPATIBIITY WITH OLDER CROSS-DISPERSION FORMULATION
C --------OLDER FORMULATION DID DIFFERENTLY AND Y WAS + UPWARD IN MANY APPLICATIONS         
c        IF(IOUTNORMAL. EQ.0) THEN 
c          DO N=1,NODES
c            IF(ICBUND(N).EQ.0) CYCLE
cC3-------GO OVER UPPER CONNECTIONS OF NODE N AND FILL Vx, Vy, Vz IN BOTH
c            DO II = IA(N)+1,IA(N+1)-1
c              JJ = JA(II)
c              IF(JJ.GT.N.AND.JJ.LE.NODES)THEN
c                IIS = JAS(II)
c                IF(ICBUND(JJ).NE.0)THEN   
c                   ANGLE = ARAD(IIS) 
c                   IF(ANGLE .GT. 1.5707. AND. ANGLE. LT. 1.5708)THEN
c                     ARAD(IIS) = -1.570796 ! CONVERT + 90 DEGREES TO - 90 DEGREES  
c                   ELSEIF(ANGLE .GT. 4.712. AND. ANGLE. LT. 4.713)THEN 
c                     ARAD(IIS) = 1.570796  ! CONVERT - 90 DEGREES TO + 90 DEGREES
c                   ENDIF  
c                ENDIF    
c              ENDIF     
c            ENDDO       
c          ENDDO         
c        ENDIF     
C--------        
C9------READ DISPERSION COEFFICIENTS
        IF(IDISP.EQ.1)THEN
C9A---------READ ISOTROPIC COEFFICIENTS
          DO K = 1,NLAY
            NNDLAY = NODLAY(K)
            NSTRT = NODLAY(K-1)+1
            NDSLAY = NNDLAY - NODLAY(K-1)
            CALL U1DREL(DLX(NSTRT),ANAME(7),NDSLAY,K,IN,IOUT)
            DO NL = NSTRT,NNDLAY
              DLY(NL) = DLX(NL)
              DLZ(NL) = DLX(NL)
            ENDDO
          ENDDO
C
          DO K = 1,NLAY
            NNDLAY = NODLAY(K)
            NSTRT = NODLAY(K-1)+1
            NDSLAY = NNDLAY - NODLAY(K-1)
            CALL U1DREL(ATXY(NSTRT),ANAME(8),NDSLAY,K,IN,IOUT)
            DO NL = NSTRT,NNDLAY
              ATYZ(NL) = ATXY(NL)
              ATXZ(NL) = ATXY(NL)
            ENDDO
          ENDDO
        ELSEIF(IDISP.EQ.2)THEN
C9B---------READ ANISOTROPIC COEFFICIENTS
          DO K = 1,NLAY
            NNDLAY = NODLAY(K)
            NSTRT = NODLAY(K-1)+1
            NDSLAY = NNDLAY - NODLAY(K-1)
            CALL U1DREL(DLX(NSTRT),ANAME(9),NDSLAY,K,IN,IOUT)
          ENDDO
C
          DO K = 1,NLAY
            NNDLAY = NODLAY(K)
            NSTRT = NODLAY(K-1)+1
            NDSLAY = NNDLAY - NODLAY(K-1)
            CALL U1DREL(DLY(NSTRT),ANAME(10),NDSLAY,K,IN,IOUT)
          ENDDO
C
          DO K = 1,NLAY
            NNDLAY = NODLAY(K)
            NSTRT = NODLAY(K-1)+1
            NDSLAY = NNDLAY - NODLAY(K-1)
            CALL U1DREL(DLZ(NSTRT),ANAME(11),NDSLAY,K,IN,IOUT)
          ENDDO
C
          DO K = 1,NLAY
            NNDLAY = NODLAY(K)
            NSTRT = NODLAY(K-1)+1
            NDSLAY = NNDLAY - NODLAY(K-1)
            CALL U1DREL(ATXY(NSTRT),ANAME(12),NDSLAY,K,IN,IOUT)
          ENDDO
C
          DO K = 1,NLAY
            NNDLAY = NODLAY(K)
            NSTRT = NODLAY(K-1)+1
            NDSLAY = NNDLAY - NODLAY(K-1)
            CALL U1DREL(ATYZ(NSTRT),ANAME(13),NDSLAY,K,IN,IOUT)
          ENDDO
C
          DO K = 1,NLAY
            NNDLAY = NODLAY(K)
            NSTRT = NODLAY(K-1)+1
            NDSLAY = NNDLAY - NODLAY(K-1)
            CALL U1DREL(ATXZ(NSTRT),ANAME(14),NDSLAY,K,IN,IOUT)
          ENDDO
C
        ENDIF
      ENDIF
C------------------------------------------------------------------
C10-----READ HEAT PARAMETERS IF IHEAT IS ON
      IF(IHEAT.EQ.1) THEN
C10A----READ AND PREPARE THE SOIL HEAT CAPACITY TERM
        DO K=1,NLAY
          NNDLAY = NODLAY(K)
          NSTRT = NODLAY(K-1)+1
          NDSLAY = NNDLAY - NODLAY(K-1)
          CALL U1DREL(TEMPC(NSTRT),ANAME(20),NDSLAY,K,IN,IOUT)
        ENDDO
C10B-------PREPARE CONSTANT PART OF HEAT CAPACITY TERM
        DO N = 1, NODES
          ADSORB(N,MCOMPT) = TEMPC(N)
          ADSORB(N,MCOMPT) = ADSORB(N,MCOMPT) /RHOW / HTCAPW * BULKD(N)
        ENDDO
C10C----READ AND PREPARE THE SOIL HEAT CONDUCTIVITY TERM
        DO K=1,NLAY
          NNDLAY = NODLAY(K)
          NSTRT = NODLAY(K-1)+1
          NDSLAY = NNDLAY - NODLAY(K-1)
          CALL U1DREL(HTCONDS(NSTRT),ANAME(21),NDSLAY,K,IN,IOUT)
      ENDDO
C10D-------PREPARE CONSTANT PART OF SOIL HEAT CONDUCTIVITY TERM
        DO N = 1, NODES
          HTCONDS(N) = HTCONDS(N) * (1.0 - PRSITY(N)) / RHOW / HTCAPW
        ENDDO
      ENDIF
C-----------------------------------------------------------------------
C10E-----READ AIR-WATER INTERFACE ADSORPTION PARAMETERS IF AW_ADSORB OPTION
      IF(IAW_ADSORB.EQ.1) THEN 
        CALL AW_ADSORB1RP1 (IN) 
      ENDIF    
C
C10F----READ AND PREPARE SOLUBILITY PARAMETERS
      IF(ISOLUBILITY .EQ.1) THEN
        CALL U2DREL(SOLLIM(1),ANAME(25),1,MCOMP,1,IN,IOUT)
        CALL U2DREL(SOLSLOPE(1),ANAME(26),1,MCOMP,1,IN,IOUT)
      ENDIF
C
C-----------------------------------------------------------------------
C11------READ SPECIES DEPENDENT ARRAYS FOR ALL SPECIES AS NEEDED, LAYER BY LAYER
C---COMPONENT MCOMPT IS TEMPERATURE IF IHEAT IS ON AND COMPONENTS ABOVE THAT ARE IMMOBILE SOLUTES
C
      DO ICOMP=1,NTCOMP
        WRITE(IOUT,10) ICOMP
10      FORMAT(80('-')/1X,'THE FOLLOWING ARRAYS ARE READ FOR COMPONENT',
     *    1X,'NUMBER ',I3)
        IF(IHEAT.EQ.1.AND.ICOMP.EQ.MCOMPT) GO TO 100 !TEMPERATURE DOES NOT HAVE ADSORPTION OR DECAY OR PARENT
        IF(ICOMP .GT. MCOMPT) GO TO 100     !IMMOBILE COMPONENTS DO NOT HAVE ADSORPTION OR DECAY OR PARENT
C
C---------------------------------------------------------
C11A-FIRST---IF CHAIN DECAY THEN READ PARENT DETAILS FOR ALL MOBILE SPECIES
        IF(ICHAIN.NE.0)THEN
          READ(IN,*) NPARENT(ICOMP)
          NPAREN = NPARENT(ICOMP)
          WRITE(IOUT,12) ICOMP, NPARENT(ICOMP)
12        FORMAT(1X,'FOR CHAIN DECAY, COMPONENT NUMBER',I5,' HAS',I5,
     1     ' PARENT COMPONENTS')
          K=1
          IF(NPAREN.GT.0) THEN
            CALL U1DINT(ITEMPC(1),ANAME(23),NPAREN,K,IN,IOUT)
            DO IPAREN = 1,NPAREN
              JPARENT(ICOMP,IPAREN) = ITEMPC(IPAREN)
            ENDDO  
C
            IF (ISPRCT.EQ.0) THEN
              WRITE(IOUT,13) ICOMP, NPAREN
13            FORMAT(1X,'*** GENERATION FACTOR OF DAUGHTER SPECIES ',I5,
     1           ' FOR EACH OF THE ',I5,' PARENT SPECIES ***')
              CALL U1DREL(TEMPC(1),ANAME(24),NPAREN,K,IN,IOUT)
              DO IPAREN = 1,NPAREN
                STOTIO(ICOMP,IPAREN) = TEMPC(IPAREN)  
              ENDDO    
            ELSE
              DO JCOMP = 1,NPAREN
                WRITE(IOUT,14) ICOMP, JPARENT(ICOMP,JCOMP)
14              FORMAT(1X,'*** GENERATION FACTOR OF DAUGHTER SPECIES ',
     1          I5,'FROM PARENT SPECIES ',I5,' FOR EACH LAYER ***')
                DO K=1,NLAY
                  NNDLAY = NODLAY(K)
                  NSTRT = NODLAY(K-1)+1
                  NDSLAY = NNDLAY - NODLAY(K-1)
                  CALL U1DREL(TEMPC(NSTRT),ANAME(4),NDSLAY,K,IN,IOUT)
                ENDDO
                DO N = 1, NODES
                  SPTLRCT(N,ICOMP,JCOMP) = TEMPC(N)
                ENDDO
              ENDDO
            ENDIF
          ENDIF
        ENDIF
C---------------------------------------------------------
C11A-SECOND---IF AW_ADSORB OPTION IS ON THEN READ AW-PARTITIONING DATA
        IF(IAW_ADSORB.NE.0)THEN
          CALL AW_ADSORB1RP2(IN,ICOMP)   
      ENDIF    
C----------------------------------------------------------
C11B-----READ ADSORPTION COEFFICIENT FOR EACH SPECIES
      IF(IADSORB.NE.0)THEN
        DO K=1,NLAY
          NNDLAY = NODLAY(K)
          NSTRT = NODLAY(K-1)+1
          NDSLAY = NNDLAY - NODLAY(K-1)
          CALL U1DREL(TEMPC(NSTRT),ANAME(4),NDSLAY,K,IN,IOUT)
        ENDDO
C11C-------PREPARE CONSTANT PART OF ADSORPTION TERM
        DO N = 1, NODES
          ADSORB(N,ICOMP) = TEMPC(N)
          ADSORB(N,ICOMP) = ADSORB(N,ICOMP) * BULKD(N)
        ENDDO
C---------------------------------------------------------
C11D-----READ FREUNDLICH EXPONENT OR LANGMUIR COEFFICIENT FOR EACH SPECIES
        IF(IADSORB.EQ.2.OR.IADSORB.EQ.3)THEN
          DO K=1,NLAY
            NNDLAY = NODLAY(K)
            NSTRT = NODLAY(K-1)+1
            NDSLAY = NNDLAY - NODLAY(K-1)
            CALL U1DREL(TEMPC(NSTRT),ANAME(19),NDSLAY,K,IN,IOUT)
          ENDDO
          DO N = 1, NODES
            FLICH(N,ICOMP) = TEMPC(N)
          ENDDO
        ENDIF
      ENDIF
C---------------------------------------------------------
C11D-----READ VARIABLE SOLUBILITY      
      IF(ISOLUBILITY .EQ.2) THEN
        DO K=1,NLAY
          NNDLAY = NODLAY(K)
          NSTRT = NODLAY(K-1)+1
          NDSLAY = NNDLAY - NODLAY(K-1)
          CALL U1DREL(TEMPC(NSTRT),ANAME(25),NDSLAY,K,IN,IOUT)
        ENDDO
        DO N = 1, NODES
          SOLLIMVAR(N,ICOMP) = TEMPC(N)
        ENDDO
C
        DO K=1,NLAY
          NNDLAY = NODLAY(K)
          NSTRT = NODLAY(K-1)+1
          NDSLAY = NNDLAY - NODLAY(K-1)
          CALL U1DREL(TEMPC(NSTRT),ANAME(26),NDSLAY,K,IN,IOUT)
        ENDDO
        DO N = 1, NODES
          SOLSLOPEVAR(N,ICOMP) = TEMPC(N)
        ENDDO          
      ENDIF     
C---------------------------------------------------------
C13-----READ ZERO ORDER DECAY COEFFICIENTS
      IF(IZOD.EQ.1)THEN
        DO K=1,NLAY
          NNDLAY = NODLAY(K)
          NSTRT = NODLAY(K-1)+1
          NDSLAY = NNDLAY - NODLAY(K-1)
          CALL U1DREL(TEMPC(NSTRT),ANAME(15),NDSLAY,K,IN,IOUT)
        ENDDO
        DO N = 1, NODES
          ZODRW(N,ICOMP) = TEMPC(N)
        ENDDO
      ELSE
        IF(IADSORB.GT.0)THEN
          IF(IZOD.EQ.2)THEN
            DO K=1,NLAY
              NNDLAY = NODLAY(K)
              NSTRT = NODLAY(K-1)+1
              NDSLAY = NNDLAY - NODLAY(K-1)
              CALL U1DREL(TEMPC(NSTRT),ANAME(16),NDSLAY,K,IN,IOUT)
            ENDDO
            DO N = 1, NODES
              ZODRS(N,ICOMP) = TEMPC(N)
            ENDDO
          ELSEIF(IZOD.EQ.3)THEN
            DO K=1,NLAY
              NNDLAY = NODLAY(K)
              NSTRT = NODLAY(K-1)+1
              NDSLAY = NNDLAY - NODLAY(K-1)
              CALL U1DREL(TEMPC(NSTRT),ANAME(15),NDSLAY,K,IN,IOUT)
            ENDDO
            DO N = 1, NODES
              ZODRW(N,ICOMP) = TEMPC(N)
            ENDDO
            DO K=1,NLAY
              NNDLAY = NODLAY(K)
              NSTRT = NODLAY(K-1)+1
              NDSLAY = NNDLAY - NODLAY(K-1)
              CALL U1DREL(TEMPC(NSTRT),ANAME(16),NDSLAY,K,IN,IOUT)
            ENDDO
            DO N = 1, NODES
              ZODRS(N,ICOMP) = TEMPC(N)
            ENDDO
          ENDIF
        ENDIF  
C13A-----READ ZERO ORDER DECAY COEFFICIENTS FOR A-W INTERFACE         
        IF(IAW_ADSORB.NE.0)THEN
          IF(IZOD.EQ.4)THEN
            DO K=1,NLAY
              NNDLAY = NODLAY(K)
              NSTRT = NODLAY(K-1)+1
              NDSLAY = NNDLAY - NODLAY(K-1)
              CALL U1DREL(TEMPC(NSTRT),ANAME(15),NDSLAY,K,IN,IOUT)
            ENDDO
            DO N = 1, NODES
              ZODRW(N,ICOMP) = TEMPC(N)
            ENDDO
            DO K=1,NLAY
              NNDLAY = NODLAY(K)
              NSTRT = NODLAY(K-1)+1
              NDSLAY = NNDLAY - NODLAY(K-1)
              CALL U1DREL(TEMPC(NSTRT),ANAME(16),NDSLAY,K,IN,IOUT)
            ENDDO
            DO N = 1, NODES
              ZODRS(N,ICOMP) = TEMPC(N)
            ENDDO
            DO K=1,NLAY
              NNDLAY = NODLAY(K)
              NSTRT = NODLAY(K-1)+1
              NDSLAY = NNDLAY - NODLAY(K-1)
              CALL U1DREL(TEMPC(NSTRT),ANAME(27),NDSLAY,K,IN,IOUT)
            ENDDO
            DO N = 1, NODES
              ZODRAW(N,ICOMP) = TEMPC(N)
            ENDDO 
          ENDIF  
        ENDIF    
      ENDIF
C---------------------------------------------------------
C14-----READ FIRST ORDER DECAY COEFFICIENTS
      IF(IFOD.EQ.1)THEN
        DO K=1,NLAY
          NNDLAY = NODLAY(K)
          NSTRT = NODLAY(K-1)+1
          NDSLAY = NNDLAY - NODLAY(K-1)
          CALL U1DREL(TEMPC(NSTRT),ANAME(17),NDSLAY,K,IN,IOUT)
        ENDDO
        DO N = 1, NODES
          FODRW(N,ICOMP) = TEMPC(N)
        ENDDO
      ELSE
        IF(IADSORB.GT.0)THEN
          IF(IFOD.EQ.2)THEN
            DO K=1,NLAY
              NNDLAY = NODLAY(K)
              NSTRT = NODLAY(K-1)+1
              NDSLAY = NNDLAY - NODLAY(K-1)
              CALL U1DREL(TEMPC(NSTRT),ANAME(18),NDSLAY,K,IN,IOUT)
            ENDDO
            DO N = 1, NODES
              FODRS(N,ICOMP) = TEMPC(N)
            ENDDO
          ELSEIF(IFOD.EQ.3)THEN
            DO K=1,NLAY
              NNDLAY = NODLAY(K)
              NSTRT = NODLAY(K-1)+1
              NDSLAY = NNDLAY - NODLAY(K-1)
              CALL U1DREL(TEMPC(NSTRT),ANAME(17),NDSLAY,K,IN,IOUT)
            ENDDO
            DO N = 1, NODES
              FODRW(N,ICOMP) = TEMPC(N)
            ENDDO
            DO K=1,NLAY
              NNDLAY = NODLAY(K)
              NSTRT = NODLAY(K-1)+1
              NDSLAY = NNDLAY - NODLAY(K-1)
              CALL U1DREL(TEMPC(NSTRT),ANAME(18),NDSLAY,K,IN,IOUT)
            ENDDO
            DO N = 1, NODES
              FODRS(N,ICOMP) = TEMPC(N)
            ENDDO
          ENDIF
        ENDIF  
C14A-----READ FIRST ORDER DECAY COEFFICIENTS FOR A-W INTERFACE        
        IF(IAW_ADSORB.NE.0)THEN
          IF(IFOD.EQ.4)THEN
            DO K=1,NLAY
              NNDLAY = NODLAY(K)
              NSTRT = NODLAY(K-1)+1
              NDSLAY = NNDLAY - NODLAY(K-1)
              CALL U1DREL(TEMPC(NSTRT),ANAME(17),NDSLAY,K,IN,IOUT)
            ENDDO
            DO N = 1, NODES
              FODRW(N,ICOMP) = TEMPC(N)
            ENDDO
            DO K=1,NLAY
              NNDLAY = NODLAY(K)
              NSTRT = NODLAY(K-1)+1
              NDSLAY = NNDLAY - NODLAY(K-1)
              CALL U1DREL(TEMPC(NSTRT),ANAME(18),NDSLAY,K,IN,IOUT)
            ENDDO
            DO N = 1, NODES
              FODRS(N,ICOMP) = TEMPC(N)
            ENDDO
            DO K=1,NLAY
              NNDLAY = NODLAY(K)
              NSTRT = NODLAY(K-1)+1
              NDSLAY = NNDLAY - NODLAY(K-1)
              CALL U1DREL(TEMPC(NSTRT),ANAME(28),NDSLAY,K,IN,IOUT)
            ENDDO
            DO N = 1, NODES
              FODRAW(N,ICOMP) = TEMPC(N)
            ENDDO 
          ENDIF  
        ENDIF         
      ENDIF
C
C---------------------------------------------------------
C15-----READ CONCENTRATION OF EACH SPECIES
100   CONTINUE
      CNAME = ANAME(5)
      IF(IHEAT.EQ.1.AND.ICOMP.EQ.MCOMPT) CNAME = ANAME(22)
      IF(IDPIN.EQ.0) THEN !-----------------------------SINGLE PRECISION READ
      DO K=1,NLAY
        NNDLAY = NODLAY(K)
        NSTRT = NODLAY(K-1)+1
        NDSLAY = NNDLAY - NODLAY(K-1)
        CALL U1DREL(TEMPC(NSTRT),CNAME,NDSLAY,K,IN,IOUT)
      ENDDO
      DO N = 1,NODES
        CONC(N,ICOMP) = TEMPC(N)
      ENDDO
      ELSE !-----------------------------DOUBLE PRECISION READ
      DO K=1,NLAY
        NNDLAY = NODLAY(K)
        NSTRT = NODLAY(K-1)+1
        NDSLAY = NNDLAY - NODLAY(K-1)
        CALL U1DREL8(TEMPC8(NSTRT),CNAME,NDSLAY,K,IN,IOUT)
      ENDDO
      DO N = 1,NODES
        CONC(N,ICOMP) = TEMPC8(N)
      ENDDO          
      ENDIF 
C
C---------------------------------------------------------
      ENDDO  ! SPECIES DO LOOP
C---------------------------------------------------------
      IF(IADSORB.NE.0.OR.IHEAT.EQ.1) DEALLOCATE(BULKD)
      DEALLOCATE(TEMPC)
      DEALLOCATE(TEMPC8)
      DEALLOCATE(ITEMPC)
C---------------------------------------------------------
C16---RETURN
      RETURN
      END
C----------------------------------------------------------------------------
      SUBROUTINE SGWT2BCT1CLN(IN,IC_IBOUND_FLG)
C     ******************************************************************
C-----READ PARAMETERS FOR CLN NODES AND INTERACTION WITH GRID-BLOCKS
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,      ONLY:IOUT,NCOL,NROW,NLAY,ITRSS,IA,JA,NJA,
     1                      IFREFM,IBOUND,NODES,NODLAY,IDPIN
      USE CLN1MODULE, ONLY: NCLNNDS,NCLNGWC
      USE GWTBCTMODULE,ONLY: PRSITY,CONC,ICBUND,CINACT,MCOMP,IADSORB,
     1  ADSORB,FLICH,IDISP,DLX,ATXY,IZOD,IFOD,ZODRW,FODRW,MCOMPT,IHEAT,
     1  NTCOMP,ICHAIN,ISPRCT,SPTLRCT
C
      REAL,  DIMENSION(:),    ALLOCATABLE ::TEMPC
      REAL*8,  DIMENSION(:),    ALLOCATABLE ::TEMPC8
      CHARACTER*24 ANAME(8),CNAME
      DATA ANAME(1) /'                  ICBUND'/
      DATA ANAME(2) /'DISP COEF ALONG CLN TUBE'/
      DATA ANAME(3) /'DISP COEF FOR CLN-MATRIX'/
      DATA ANAME(4) /'  ZERO-ORDER DECAY WATER'/
      DATA ANAME(5) /' FIRST-ORDER DECAY WATER'/
      DATA ANAME(6) /'     INITIAL CONC IN CLN'/
      DATA ANAME(7) /'     INITIAL TEMP IN CLN'/
      DATA ANAME(8) /'STOTIOMETRIC COEFFICENTS'/

C     ------------------------------------------------------------------
C
      WRITE(IOUT,1)
1     FORMAT(/1X,85('-')/1X,'READ TRANSPORT PARAMETERS AND INITIAL',
     1 ' CONDITIONS FOR CLN NODES'/1X,85('-'))
C0----SET INDICES
      NSTRT = NODES+1
      NDSLAY = NCLNNDS
C1------READ/SET ICBUND ARRAY
      IF(IC_IBOUND_FLG.EQ.0)THEN
C1A-------READ ICBUND
        K = 1
        CALL U1DINT(ICBUND(NSTRT),ANAME(1),NDSLAY,K,IN,IOUT)
C1B-------CHECK FOR ERRORS
         DO N=NODES+1,NODES + NDSLAY
           IF(ICBUND(N).NE.0.AND.IBOUND(N).EQ.0)THEN
             WRITE(IOUT,55) N
55    FORMAT(5X,'*** ACTIVE TRANSPORT NODE WHERE FLOW WAS INACTIVE AT ',
     *  'NODE: ',I10,' ***')
             STOP
           ENDIF
         ENDDO
      ELSE
C1C-------SET ICBUND
        DO N=NODES+1, NODES + NCLNNDS
          ICBUND(N) = IBOUND(N)
        ENDDO
      ENDIF
C
C2-----READ DISPERSION COEFFICIENTS ALONG CONDUIT AND FOR CONDUIT-MATRIX INTERACTION
      K = 0
      IF(IDISP.NE.0)THEN
        NSTRT = NODES+1
        NDSLAY = NCLNNDS
        CALL U1DREL(DLX(NSTRT),ANAME(2),NDSLAY,K,IN,IOUT)
        CALL U1DREL(ATXY(NSTRT),ANAME(3),NDSLAY,K,IN,IOUT)
      ENDIF
C
C------------------------------------------------------------------
C3-------READ COMPONENT DEPENDENT PARAMETERS AND INITIAL CONCENTRATIONS
      ALLOCATE(TEMPC(NCLNNDS))
      DO ICOMP=1,NTCOMP
        WRITE(IOUT,10) ICOMP
10      FORMAT(80('-')/1X,'THE FOLLOWING ARRAYS ARE READ FOR COMPONENT',
     *    1X,'NUMBER ',I3)
C------------------------------------------------------------------
C3B------READ SPATIALLY VARIABLE STOCHIOMETRY
        IF(ISPRCT.NE.0)THEN
          DO JCOMP = 1,MCOMP
            WRITE(IOUT,14) ICOMP, JCOMP
14          FORMAT(1X,'*** GENERATION FACTOR OF DAUGHTER SPECIES ',
     1      I5, 'FROM PARENT SPECIES ',I5,' FOR EACH LAYER ***')
            NSTRT = NODES+1
            NDSLAY = NCLNNDS
            CALL U1DREL(TEMPC,ANAME(8),NDSLAY,K,IN,IOUT)
            DO N = 1,NCLNNDS
              SPTLRCT(N+NODES,ICOMP,JCOMP) = TEMPC(N)
            ENDDO
          ENDDO
        ENDIF
C------------------------------------------------------------------
C4-------ZERO-ORDER DECAY IN WATER IF MOBILE SOLUTE
        IF(IHEAT.EQ.1.AND.ICOMP.EQ.MCOMPT) GO TO 100 !TEMPERATURE DOES NOT HAVE DECAY TERMS
        IF(ICOMP .GT. MCOMPT) GO TO 100     !IMMOBILE COMPONENTS DO NOT HAVE DECAY TERMS
        IF(IZOD.EQ.1.OR.IZOD.EQ.3)THEN
          NSTRT = NODES+1
          NDSLAY = NCLNNDS
          CALL U1DREL(TEMPC,ANAME(4),NDSLAY,K,IN,IOUT)
          DO N = 1,NCLNNDS
            ZODRW(N+NODES,ICOMP) = TEMPC(N)
          ENDDO
        ENDIF
C------------------------------------------------------------------
C5-------FIRST-ORDER DECAY IN WATER
        IF(IFOD.EQ.1.OR.IFOD.EQ.3)THEN
          NSTRT = NODES+1
          NDSLAY = NCLNNDS
          CALL U1DREL(TEMPC,ANAME(5),NDSLAY,K,IN,IOUT)
          DO N = 1,NCLNNDS
            FODRW(N+NODES,ICOMP) = TEMPC(N)
          ENDDO
        ENDIF
C------------------------------------------------------------------
C6-------INITIAL CONCENTRATIONS
100     CONTINUE
        CNAME = ANAME(6)
        IF(IHEAT.EQ.1.AND.ICOMP.EQ.MCOMPT) CNAME = ANAME(7)
        IF(IDPIN.EQ.0) THEN !----------------------SINGLE PRECISION READ
        CALL U1DREL(TEMPC,CNAME,NDSLAY,K,IN,IOUT)
        DO N = 1,NCLNNDS
        CONC(N+NODES,ICOMP) = TEMPC(N)
        ENDDO
        ELSE !----------------------DOUBLE PRECISION READ
        ALLOCATE(TEMPC8(NCLNNDS))    
        CALL U1DREL8(TEMPC8,CNAME,NDSLAY,K,IN,IOUT)
        DO N = 1,NCLNNDS
        CONC(N+NODES,ICOMP) = TEMPC(N)
        ENDDO         
        DEALLOCATE(TEMPC8)    
      ENDIF   
C---------------END SPECIES LOOP      
      ENDDO 
C---------------END SPECIES LOOP        
      DEALLOCATE(TEMPC)
C
C------------------------------------------------------------------
C7----RETURN
      RETURN
      END
C----------------------------------------------------------------------
      SUBROUTINE GWT2BCT1SOLVE(KITER,KSTP,KPER)
C     ******************************************************************
C     SOLVE THE TRANSPORT EQUATION FOR ALL COMPONENTS
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,   ONLY:IOUT,NCOL,NROW,NLAY,IBOUND,NODLAY,ISYM,
     1  AMAT,RHS,IA,JA,PGF,NJA,NODES,IUNIT,IUNSTR,NEQS,INCLN,IDPT,
     2  ISSFLG,HNEW,HOLD,TOP,BOT,Sn,So,LAYNOD,FMBE,IFMBC,ImdT
      USE CLN1MODULE, ONLY: NCLNNDS
      USE GWFBCFMODULE, ONLY: LAYCON
      USE GWFBASMODULE,ONLY:DELT,TOTIM,ISPFASTC,ITSFASTC
      USE TVMU2MODULE
      USE GWTBCTMODULE, ONLY: MCOMP,CONC,CONCO,NTITER,ICBUND,CICLOSE,
     1    CINACT,CNCG,LRCC,MSUMT,IDISP,FCNCG,FLRCC,ADMAT,IXDISP,ITVD,
     1    DXCS,DYCS,DZCS,NTCOMP,MCOMPT,IHEAT,HTCONDT,NSEQITR,DCXL,DCYL,
     1    NODETIBC,NICB1,CTIB,TIMEWEIGHT,ICHAIN,ISOLUBILITY,VXL,VYL,
     1    IAW_ADSORB,ISATADSORB,ADSORB,IADSORB,DIFFNC
      USE GWTDPTMODULE, ONLY: IAW_ADSORBIM,ADSORBIM,IADSORBIM
      USE GWFDPFMODULE, ONLY: SnIM 
      USE GWTmdTMODULE, ONLY: TSHIFTMD
      USE GWFLAKMODULE,ONLY: ILKTRNSPT
       DOUBLE PRECISION, SAVE, DIMENSION(:), ALLOCATABLE :: CNEW
       DOUBLE PRECISION DELTO,DTERMS,RTERMS,DTERMD,RTERMD,VODT,SAT,CTMP,
     1    TWM1,TOTIMMD
      save itp,DELTO
C     ------------------------------------------------------------------
C1 --------INITIALIZE TRANSPORT FOR ALL COMPONENTS
C     ------------------------------------------------------------------
C1A--------INDICATE THAT SOLUTION IS FOR TRANSPORT
          WRITE(IOUT,38)
38        FORMAT(1X,108('-'))
          CALL UMESPR('SOLVING FOR TRANSPORT',' ',IOUT)
          WRITE(*,25)KPER,KSTP
   25     FORMAT(' Solving:  Stress period: ',i5,4x,
     &       'Time step: ',i5,4x,'Groundwater Transport Eqn.')
C1B------SET FLAGS
      ICNVG = 0
      ISS=ISSFLG(KPER)
      IF(KSTP.EQ.1.AND.KPER.EQ.1) DELTO = 0.0
C-------------------------------------------------------------------------
C1C----ADVANCE SOLUTION
      CALL GWT2BCT1AD(KITER,KSTP,KPER,ISS)
C1D--------advance the solution on MD terms
C----------------------------------------------------------------------
      IF(ImdT.GT.0) THEN
        TOTIMMD = TOTIM + TSHIFTMD
        CALL matrixdsr2(TOTIMMD)
      ENDIF
C1E ----SET ADSORBTION FOR SATURATED PORTION OF DOMAIN OF SATADSORB OPTION IS ON
C ------NOTE THAT ADSORPTION CAN NEVER OCCUR DIFFERENTLY FROM SATURATED PORTION FOR SS CASE
C ----- (OF FIRST TIME STEP IF SIMULATION IS TRANSIENT)      
      IF(KSTP.EQ.1. AND. KPER.EQ.1) THEN 
        IF(ISATADSORB.EQ.1) THEN
          DO ICOMP = 1, MCOMP
            DO N = 1, NODES
              ADSORB(N,ICOMP) = ADSORB(N,ICOMP) * SN(N)
            ENDDO  
          ENDDO
C1F  -------ALSO FOR DUAL POROSITY DOMAIN IF PRESENT          
          IF(IDPT.GT.0. AND. IADSORBIM. GT. 0) THEN 
             DO ICOMP = 1, MCOMP
              DO N = 1, NODES
                ADSORBIM(N,ICOMP) = ADSORBIM(N,ICOMP) * SNIM(N)
              ENDDO  
            ENDDO  
          ENDIF    
        ENDIF
      ENDIF  
C-----------------------------------------------------------------------
C2----FILL TERMS CONSTANT TO COMPONENTS AND TRANSPORT ITERATIONS IN ADMAT
C-----------------------------------------------------------------------
      IADMATF = 1
      IF(KSTP.GT.1.AND.ISS.EQ.1)IADMATF = 0
C
      IF(ISPFASTC.GT.0.AND.ITSFASTC.GT.0)THEN
C---------FILL MATRIX IF COMING IN FROM A RESTART
        IF(ISPFASTC.EQ.KPER.AND.ITSFASTC+1.EQ.KSTP) IADMATF = 1
      ENDIF
C2A-----NO NEED TO REFILL ADMAT IF IADMATF IS ZERO
      IF(IADMATF.EQ.1)THEN
C-------ZERO OUT ADMAT ARRAY
        DO N=1,NJA
          ADMAT(N) = 0.0
        ENDDO
        IF(IHEAT.EQ.1) THEN 
          htcondt = 0.0
        ENDIF             
C-----------------------------------------------------------------------
C2B-------COMPUTE DISPERSION TERM
        IF(IDISP.NE.0) THEN
C2C---------CALCULATE VELOCITY COMPONENTS FOR DISPERSION TERM 
          CALL GWT2BCT1VELCALC
          DO N=1,NEQS
            IF(ICBUND(N).NE.0)THEN
C2D-------------FILL DISPERSIVE TERM IN ROWS OF ACTIVE TRANSPORT NODES
              CALL GWT2BCT1DSP(N)
            ENDIF
          ENDDO
        ENDIF
C-----------------------------------------------------------------------
C2E-------COMPUTE ADVECTION TERM
        DO N=1,NEQS
          IF(ICBUND(N).NE.0)THEN
C2F-----------FILL UPSTREAM ADVECTIVE TERM FOR ACTIVE TRANSPORT NODES
            CALL GWT2BCT1ADV(N)
          ENDIF
        ENDDO
C2G-----------SUBTRACT FLOW BALANCE ERROR FROM DIAGONAL
        IF(IFMBC.NE.0)THEN
          amaxerr = 0.0
          nmax = 0
          DO N=1,NEQS
            IF(ICBUND(N).NE.0)THEN
              aerr = abs(fmbe(n)  )
              if(aerr.gt.amaxerr) then
                  amaxerr = aerr
                  nmax = n
              endif
              ADMAT(IA(N)) = ADMAT(IA(N)) - FMBE(N)
            ENDIF
          ENDDO
          write(iout,*)'max flow fmbe and n are',amaxerr,nmax
        ENDIF
C----------------------------------------------------------------------
C2H-------FILL HCONDT MATRIX FOR HEAT CONDUCTANCE IF HEAT IS SOLVED
        IF(IHEAT.EQ.1)THEN
C2I-------------FILL HEAT CONDUCTIVITY TERM IN ROWS OF ACTIVE TRANSPORT NODES
          CALL GWT2HEATCOND
        ENDIF
C----------------------------------------------------------------------
C2J-------FILL AREA TERM FOR AIR-WATER INTERFACE IF ADSORPTION ON INTERFACE
        IF(IAW_ADSORB.EQ.1)THEN
          CALL AWI_AREA
        ENDIF
C----------------------------------------------------------------------
C2K-------FILL AREA TERM FOR AIR-WATER INTERFACE IN DUAL POROSITY DOMAIN IF ADSORPTION ON INTERFACE
        IF(IDPT.GT.0) THEN  
          IF(IAW_ADSORBIM.EQ.1)THEN
            CALL AWI_AREAIM
          ENDIF
        ENDIF
C-----------------------------------------------------------------------
      ENDIF   ! ENDIF FOR IF(IADMATF.EQ.1)
C-----------------------------------------------------------------------
C3----LOOP OVER ALL MOBILE TRANSPORT COMPONENTS AND SEQUENTIAL ITERATIONS
C-----------------------------------------------------------------------
      DO ISEQITR = 1,NSEQITR !SEQUENTIAL ITERATION LOOP
      WRITE(IOUT,34)ISEQITR
34    FORMAT(/10X,'SEQUENTIAL TRANSPORT ITERATION LOOP NUMBER',I5/
     1        10X,'-----------------------------------------------')

      DO ICOMP = 1,MCOMPT ! COMPONENT LOOP
C4i-------FOR HEAT EQUATION AFTER SOLUTE, REFILL ADMAT WITH ADVECTIVE TERM ONLY TO REMOVE DIFFUSION TERM
        IF(MCOMPT.GT.1. AND. ICOMP.EQ.MCOMPT) THEN 
          ADMAT = 0.0
          DIFFNCKP = DIFFNC  ! DIFFUSION FOR SOLUTE SHOULD NOT ADD TO HEAT CONDUCTION
          DIFFNC = 0.0       ! DIFFUSION FOR SOLUTE SHOULD NOT ADD TO HEAT CONDUCTION   
C4ii-------REFILL ADVECTION TERM
          DO N=1,NEQS
            IF(ICBUND(N).NE.0)THEN
C4iii-----------FILL UPSTREAM ADVECTIVE TERM FOR ACTIVE TRANSPORT NODES
              CALL GWT2BCT1ADV(N)
C4iv-----------FILL DISPERSIVE TERM IN ROWS OF ACTIVE TRANSPORT NODES
              CALL GWT2BCT1DSP(N)              
            ENDIF
          ENDDO
          DIFFNC = DIFFNCKP  ! DIFFUSION FOR SOLUTE SHOULD NOT ADD TO HEAT CONDUCTION
          IADMATF = 1   ! WILL NEED TO REFILL MATRIX AFTER THIS
        ENDIF
C ----------------------------------------------------------------------------
C4------LOOP OVER NONLINEAR TRANSPORT ITERATIONS (TVD)
        DO ITITER = 1,NTITER
C4A--------FORMULATE CELL-BY-CELL TRANSPORT TERMS IN LHS MATRIX AND RHS
          RHS = 0.0
C5---------FILL TERMS CONSTANT TO ITERATIONS AND COMPONENTS ON LHS MATRIX
          AMAT = ADMAT
C-----------------------------------------------------------------------
C5A--------ADD HEAT CONDUCTANCE TERM TO AMAT IF IHEAT IS ON AND COMPONENT IS HEAT (MCOMPT)
          IF(IHEAT.EQ.1. AND. ICOMP.EQ.MCOMPT)THEN
             AMAT = AMAT + HTCONDT
          ENDIF
C
C-----------------------------------------------------------------------
C5B---------FILL CROSS DISPERSION TERM ON RHS - ONLY FOR MATRIX BLOCKS
          IF(IDISP.GT.0.AND.IXDISP.NE.0)THEN
C5C-----------COMPUTE DCDX, DCDY, DCDZ FOR CROSS DISPERSION TERM ON RHS
            ALLOCATE(DXCS(NODES), DCXL(NJA))
            ALLOCATE(DYCS(NODES), DCYL(NJA))
            ALLOCATE(DZCS(NODES))
C
            CALL GWT2BCT1DC(ICOMP)
C5D-----------FILL CROSS DISPERSION TERM INTO RHS
            DO N=1,NODES
              IF(ICBUND(N).NE.0)THEN
                 CALL GWT2BCT1DSPX(N)
              ENDIF
            ENDDO
            DEALLOCATE(DXCS,DCXL)
            DEALLOCATE(DYCS,DCYL)
            DEALLOCATE(DZCS)
          ENDIF
C-----------------------------------------------------------------------
C5E---------FILL TVD TERM ON RHS
          IF(ITVD.GT.0)THEN
            DO N=1,NEQS
              IF(ICBUND(N).NE.0)THEN
                CALL GWT2BCT1TVD(N,ICOMP)
              ENDIF
            ENDDO
          ENDIF
C-----------------------------------------------------------------------
C5F---------ADJUST FOR TIME WEIGHTING ON CONCENTRATION TERM
C  ---------MOVE C5H (DECAY MODULE) ABOVE THIS TO INCLUDE DECAY
C-----------WITH CRANK-NICOLSON TERM AND MODIFY DECAY MASS BALANCE ALSO
          IF(TIMEWEIGHT.LT.0.99999)THEN
            TWM1 = 1.0 - TIMEWEIGHT
C5F1----------ADJUST RHS WITH TIME WEIGHTING FACTOR
            DO N=1,NEQS
              IF(ICBUND(N).NE.0)THEN
                CALL GWT2BCT1TWA(N,ICOMP,TWM1)
              ENDIF
            ENDDO
C5F2----------ADJUST MATRIX WITH TIME WEIGHTING FACTOR
            DO I=1,NJA
             AMAT(I) = AMAT(I) * TIMEWEIGHT
            ENDDO
          ENDIF
C-----------------------------------------------------------------------
C5G---------FILL STORAGE  TERMS ON RHS AND DIAGONAL
          DO N=1,NEQS
            IF(ICBUND(N).NE.0)THEN
              CALL GWT2BCT1STO(N,ICOMP,DTERMS,RTERMS,ISS)
              IPIV = IA(N)
              AMAT(IPIV) = AMAT(IPIV) + DTERMS
              RHS(N) = RHS(N) + RTERMS
            ENDIF
          ENDDO
C-----------------------------------------------------------------------
C5H---------FILL DECAY TERMS ON RHS AND DIAGONAL (no decay for heat equation OR IMMOBILE COMPONENTS)
          IF(IHEAT.EQ.1. AND. ICOMP.GE.MCOMPT) GO TO 52
          DO N=1,NEQS
            IF(ICBUND(N).NE.0)THEN
              CALL GWT2BCT1DCY(N,ICOMP,DTERMD,RTERMD)
              IPIV = IA(N)
              AMAT(IPIV) = AMAT(IPIV) +  DTERMD
              RHS(N) = RHS(N) +  RTERMD
            ENDIF
          ENDDO
52    CONTINUE
C-----------------------------------------------------------------------
C5I---------FILL PARENT GENERATION TERMS ON RHS
      IF(ICHAIN.NE.0) THEN
        IF(IHEAT.EQ.1. AND. ICOMP.EQ.MCOMPT) GO TO 53
        DO N=1,NEQS
          IF(ICBUND(N).NE.0)THEN
            CALL GWT2BCT1GEN(N,ICOMP)
          ENDIF
        ENDDO
53      CONTINUE
      ENDIF
C-----------------------------------------------------------------------
C5J---------FILL SOLUBILITY PRECIPITATION TERMS
      IF(ISOLUBILITY.NE.0) THEN
        IF(IHEAT.EQ.1. AND. ICOMP.EQ.MCOMPT) GO TO 54
        DO N=1,NEQS
          IF(ICBUND(N).NE.0)THEN
            CALL SOLUBILITY(N,ICOMP,DTERMS,RTERMS,ISS)
            IPIV = IA(N)
            AMAT(IPIV) = AMAT(IPIV) + DTERMS
            RHS(N) = RHS(N) + RTERMS
          ENDIF
        ENDDO
54      CONTINUE
      ENDIF
C----------------------------------------------------------------------
C5K-------FILL PARTITION TERM FOR AIR-WATER INTERFACE IF ADSORPTION ON INTERFACE
        IF(IAW_ADSORB.EQ.1)THEN
          IF(IHEAT.EQ.1. AND. ICOMP.GE.MCOMPT) GO TO 55  
          CALL AWI_KAWI (ICOMP)
55        CONTINUE          
        ENDIF      
C-----------------------------------------------------------------------
C5L---------ADJUST STORAGE TERMS FOR TVM PACKAGE
          IF(IUNIT(64).GT.0) THEN
          DO INODE = 1, NTVMPOR
            N=ITVMPOR(INODE)
            IF(ICBUND(N).NE.0)THEN
              CALL GWT2BCT1STOTVM(N,ICOMP,DTERMS,RTERMS,ISS)
C              IPIV = IA(N)
C              AMAT(IPIV) = AMAT(IPIV) + DTERMS + DTERMD
              RHS(N) = RHS(N) + RTERMS
            ENDIF
          ENDDO
          ENDIF
C----------------------------------------------------------------------
C5M--------FILL MATRIX DIFFUSION TRANSPORT TERMS
          IF(ImdT.GT.0) THEN
            DO N = 1,NODES
              ICOMPM1 = ICOMP-1
              IF(ICOMP.EQ.1) ICOMPM1 = ICOMP
              CALL matrixdsr3 (N,ICOMP,DELT,CONC(N,ICOMP),
     1          CONC(N,ICOMPM1),CONCO(N,ICOMP))
              CALL matrixdsr4rhs(n,ICOMP,RTERMD)
              CALL matrixdsr5diag(n,ICOMP,DTERMD)
              IPIV = IA(N)
              AMAT(IPIV) = AMAT(IPIV) - DTERMD
              RHS(N) = RHS(N) - RTERMD
            ENDDO
          ENDIF
C----------------------------------------------------------------------
C6--------FILL TRANSPORT BOUNDARY CONDITIONS
C----------------------------------------------------------------------
C6A-------PRESCRIBED CONCENTRATION BOUNDARY
          IF(IUNIT(14).GT.0) CALL GWT2PCB1FM(ICOMP)
C6B-------PRESCRIBED HEAD BOUNDARY
          IF (IUNIT(20).GT.0) CALL GWT2PHB1FM(ICOMP)
C6C-------WELL INFLOW / OUTFLOW BOUNDARY
          IF(IUNIT(2).GT.0) CALL GWT2WEL1FM(ICOMP)
C6D-------GHB INFLOW / OUTFLOW BOUNDARY
          IF(IUNIT(7).GT.0) CALL GWT2GHB1FM(ICOMP)
C6E-------DRAIN OUTFLOW BOUNDARY
          IF(IUNIT(3).GT.0) CALL GWT2DRN1FM(ICOMP)
C6F-------DRT INFLOW / OUTFLOW BOUNDARY
          IF(IUNIT(40).GT.0) CALL GWT2DRT1FM(ICOMP)
C6F-------QRT INFLOW / OUTFLOW BOUNDARY
          IF(IUNIT(41).GT.0) CALL GWT2QRT1FM(ICOMP)
C6G-------RIVER OUTFLOW BOUNDARY
          IF(IUNIT(4).GT.0) CALL GWT2RIV1FM(ICOMP)
C6H-------EVT OUTFLOW BOUNDARY
          IF(IUNIT(5).GT.0) CALL GWT2EVT1FM(ICOMP)
C6I-------ETS OUTFLOW BOUNDARY
          IF(IUNIT(39).GT.0) CALL GWT2ETS1FM(ICOMP)
C6J-------RECH INFLOW / OUTFLOW (WITH SEEPAGE) BOUNDARY
          IF(IUNIT(8).GT.0) CALL GWT2RCH1FM(ICOMP)
C6K-------LAK INFLOW / OUTFLOW BOUNDARY
          IF(IUNIT(22).GT.0) THEN
	      IF(ILKTRNSPT.EQ.1) CALL CALCLAKECONC(ICOMP,IUNIT(44))
            CALL GWT2LAK1FM(ICOMP)
      ENDIF
C----------------------------------------------------------------------
C7--------FILL DUAL POROSITY TRANSPORT TERMS
C----------------------------------------------------------------------
      itestmat = 0
      if(itestmat.eq.1)then
        open(99,file='mat_USGs.TXT')
        WRITE(99,*)'NODE, RHS, AMAT FOLLOW'
        DO N=1,NEQS
          if(ibound(n). eq.0) cycle
          I1 = IA(N)
          I2 = IA(N+1)-1
          WRITE(99,66)N,RHS(N),(AMAT(I),I=I1,I2)
c          write(99,66)n,rhs(n),(ja(i),i=i1,i2)
        ENDDO
c66      FORMAT(I9,1X,G15.3,2X,100G15.3)
66      FORMAT(I9,1X,G15.4,2X,100G15.4)
        CLOSE (99)
        stop
      endif      
          IF(IDPT.GT.0) CALL GWT2DPTU1FM(KPER,ICOMP,ISS)
C7A------REDUCE IMMOBILE DOMAIN EQUATION INTO MATRIX AND RHS
          IF(IDPT.GT.0) CALL SSMS2DPT1RED(KPER)
C----------------------------------------------------------------------
C8 ------SOLVE THE GROUNDWATER TRANSPORT EQUATION FOR THIS SPECIES
         CALL GWT2SOLVE(ITITER,ICNVG,ICOMP,KSTP,KPER,ITP,IN_ITER)
C----------------------------------------------------------------------
C9--------BACK-SUBSTITUTE FOR IMMOBILE DOMAIN CONCENTRATIONS
C----------------------------------------------------------------------
          IF(IDPT.GT.0) CALL SSMS2DPT1BKS(KPER,ICOMP)
C----------------------------------------------------------------------
C10---------EVALUATE AND STORE MAXIMUM CONCENTRATION CHANGE AND CONVERGENCE
          CALL CONCCONV(ICOMP,ITITER,ICNVG)
C
          IF(ICNVG.EQ.1) GO TO 205
C--------------------------------------------------------------------
        ENDDO            !END TRANSPORT ITERATION LOOP (ITITER)
C--------------------------------------------------------------------
205     CONTINUE
C
C11-------WRITE TVD ITERATION SUMMARY
        ITERS = ITITER-1
        IF(ICNVG.EQ.1) ITERS = ITITER
        WRITE(IOUT,1010) ITERS,KSTP,KPER,ICOMP
 1010 FORMAT(/1X,120('-')
     &  /1X, I5,' CALLS TO SPARSE MATRIX SOLVER PACKAGE ',
     &  'IN TRANSPORT TIME STEP',I8,' STRESS PERIOD',I8,' SPECIES',I3)
C
C11A------FOR BCF NODES
       CALL SSMS2BCFU1P(Cncg,LrcC,ITERS,NTITER,IOUT,IUNSTR)
C11B------FOR CLN NODES
       CALL SSMS2CLN7P(ITERS,NTITER,IOUT,IUNSTR,FCNCG,FLRCC,INCLN)
C
        WRITE(IOUT,1)ICOMP
1       FORMAT(5X,'TRANSPORT SOLUTION COMPLETE FOR COMPONENT SPECIES',
     *    1X,'NUMBER',I5/5X,60('-'))
C----------------------------------------------------------------------
      ENDDO              ! END COMPONENT LOOP (ICOMP)
C----------------------------------------------------------------------
C12-------UPDATE MASS IN MATRIX FOR MD MODULE HERE
      IF(ImdT. GT.0) CALL matrixdsr6int (DELT)
C----------------------------------------------------------------------
C12-------CALL MULTICOMPONENT REACTION OR GEOCHEMICAL MODULE HERE




C----------------------------------------------------------------------
      ENDDO              ! END SEQUENTIAL ITERATION LOOP (NSEQITR)
C--------------------------------------------------------------------------
CC-------DETERMINE WHICH PRINTOUT IS REQUIRED.
      CALL GWT2BCT1OC(KSTP,KPER,1)
C----------------------------------------------------------------------
C13------COMPUTE MASS BALANCE TERMS FOR ALL COMPONENTS
      DO ICOMP = 1,NTCOMP
C---------------------------------------------------------------------
        MSUMT = 1
C13A-------MASS STORAGE TERM IN WATER AND SOIL
        CALL GWT2STO1BDW(KSTP,KPER,ICOMP,ISS)
        IF(IADSORB. GT. 0) CALL GWT2STO1BDS(KSTP,KPER,ICOMP,ISS)
C13B-------MASS DECAY TERMS - COMBINE ZERO AND FIRST ORDER (ALSO NOT FOR HEAT OR IMMOBILE COMPONENTS)
        IF(IHEAT.EQ.1. AND. ICOMP.GE.MCOMPT) GO TO 152
        CALL GWT2DCY1BD(KSTP,KPER,ICOMP)
152     CONTINUE
        IF(ICHAIN.NE.0)THEN
          IF(IHEAT.EQ.1. AND. ICOMP.EQ.MCOMPT) GO TO 153  ! SKIP FOR HEAT
          CALL GWT2GEN1BD(KSTP,KPER,ICOMP)
153       CONTINUE
        ENDIF
C
C13C-------MASS STORAGE, DECAY, GENERATION AND PRECIP TERMS IN IMMOBILE DOMAIN
        IF(IDPT.GT.0) THEN
          CALL GWT2STOIM1BDW(KSTP,KPER,ICOMP)
          IF(IADSORBIM.GT.0) CALL GWT2STOIM1BDS(KSTP,KPER,ICOMP)
          IF(IHEAT.EQ.1. AND. ICOMP.GE.MCOMPT) GO TO 154
          CALL GWT2DCYIM1BD(KSTP,KPER,ICOMP)
154       CONTINUE
          IF(ICHAIN.NE.0)THEN
            IF(IHEAT.EQ.1. AND. ICOMP.EQ.MCOMPT) GO TO 155  ! SKIP FOR HEAT
            CALL GWT2GENIM1BD(KSTP,KPER,ICOMP)
155         CONTINUE
          ENDIF
          IF(ISOLUBILITY.NE. 0) THEN
            CALL SOLUBILITYMBIM(KSTP,KPER,ICOMP,ISS)
          ENDIF
C13F-------MASS BALANCE DUE TO ADSORPTION ON AIR-WATER INTERFACE
          IF(IAW_ADSORBIM.EQ.1) THEN
            IF(ICOMP.GT.MCOMP) GO TO 252
            CALL AW_ADSORBIM1BD(KSTP,KPER,ICOMP,ISS)
252         CONTINUE
          ENDIF
        ENDIF
C
C13C2-------MASS balance in MD domain
      IF(ImdT.GT.0) THEN
        CALL matrixdsr8 (kstp,kper,icomp)
      ENDIF
C13D-------MASS ERROR DUE TO FLOW MASS BALANCE ERROR ADJUSTMENT
        IF(IFMBC.NE.0)THEN
          CALL GWT2FMBE1BD(KSTP,KPER,ICOMP)
        ENDIF
        IF(ISOLUBILITY.NE.0) THEN
          CALL SOLUBILITYMB(KSTP,KPER,ICOMP,ISS)
      ENDIF
C13F-------MASS BALANCE DUE TO ADSORPTION ON AIR-WATER INTERFACE      
      IF(IAW_ADSORB.EQ.1) THEN
        IF(ICOMP.GT.MCOMP) GO TO 254                  
        CALL AW_ADSORB1BD(KSTP,KPER,ICOMP,ISS)
254     CONTINUE        
      ENDIF 
CCCC-----------------------------------------------------------------------
CCCC13E---------ADJUST STORAGE TERMS FOR TVM PACKAGE
CCC        IF(IUNIT(64).GT.0) THEN
CCC          CALL GWT2STO1BDTVM(KSTP,KPER,ICOMP,ISS)
CCC        ENDIF
C
C----------------------------------------------------------------------
C14---------MASS BALANCE FOR BOUNDARY COMPONENTS
C
C14A-------PRESCRIBED CONCENTRATION BOUNDARY
        IF(IUNIT(14).GT.0) CALL GWT2PCB1BD(KSTP,KPER,ICOMP)
C14B-------PRESCRIBED HEAD BOUNDARY
        IF(IUNIT(20).GT.0) CALL GWT2PHB1BD(KSTP,KPER,ICOMP)
C14C-------WELL INFLOW / OUTFLOW BOUNDARY
        IF(IUNIT(2).GT.0) CALL GWT2WEL1BD(KSTP,KPER,ICOMP)
C14D-------GHB INFLOW / OUTFLOW BOUNDARY
        IF(IUNIT(7).GT.0) CALL GWT2GHB1BD(KSTP,KPER,ICOMP)
C14E-------DRAIN INFLOW / OUTFLOW BOUNDARY
        IF(IUNIT(3).GT.0) CALL GWT2DRN1BD(KSTP,KPER,ICOMP)
C14F-------DRT INFLOW / OUTFLOW BOUNDARY
        IF(IUNIT(40).GT.0) CALL GWT2DRT1BD(KSTP,KPER,ICOMP)
C14F-------QRT INFLOW / OUTFLOW BOUNDARY
        IF(IUNIT(41).GT.0) CALL GWT2QRT1BD(KSTP,KPER,ICOMP)
C14G-------RIVER INFLOW / OUTFLOW BOUNDARY
        IF(IUNIT(4).GT.0) CALL GWT2RIV1BD(KSTP,KPER,ICOMP)
C14H-------EVT OUTFLOW BOUNDARY
        IF(IUNIT(5).GT.0) CALL GWT2EVT1BD(KSTP,KPER,ICOMP)
C14I-------ETS OUTFLOW BOUNDARY
        IF(IUNIT(39).GT.0) CALL GWT2ETS1BD(KSTP,KPER,ICOMP)
C14J-------RCH OUTFLOW BOUNDARY
        IF(IUNIT(8).GT.0) CALL GWT2RCH1BD(KSTP,KPER,ICOMP)
C14K-------LAK INFLOW / OUTFLOW BOUNDARY
        IF(IUNIT(22).GT.0) CALL GWT2LAK1BD(KSTP,KPER,ICOMP,IUNIT(44))
      ENDDO
C--------------------------------------------------------------------------
C14L-------UPDATE SOLUTE MASS IN MD ROUTINE FROM NEW TO OLD
      IF(ImdT.GT.0) CALL matrixdsr7updt
C--------------------------------------------------------------------------
C15-----PRINT AND/OR SAVE TRANSPORT RESULTS.
      CALL GWT2BCT1OT(KSTP,KPER,1)
C15A-----PRINT AND/OR SAVE DUAL POROSITY TRANSPORT RESULTS.
      IF(IDPT.GT.0) CALL GWT2DPT1OT(KSTP,KPER,1)
      IF(ImdT.GT.0) CALL GWT2mdt1OT(KSTP,KPER) 
C       
      DELTO = DELT
C
C16-----RETURN
      RETURN
      END
C----------------------------------------------------------------------
      SUBROUTINE GWT2BCT1AD(KITER,KSTP,KPER,ISS)
C     ******************************************************************
C     ADVANCE TRANSPORT SOLUTION FOR ALL COMPONENTS
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,   ONLY:IOUT,IBOUND,NODES,IUNIT,IDPT,NEQS,INCLN,
     2  ISSFLG,HNEW,HOLD,TOP,BOT,Sn,So,LAYNOD
      USE GWFBCFMODULE, ONLY: LAYCON
      USE GWTBCTMODULE, ONLY: MCOMP,CONC,CONCO,NTITER,ICBUND,CICLOSE,
     1    CINACT,CNCG,LRCC,MSUMT,IDISP,FCNCG,FLRCC,ADMAT,IXDISP,ITVD,
     1    DXCS,DYCS,DZCS,NICB1,NODETIBC,NTCOMP,MCOMPT
C     ------------------------------------------------------------------
C1----SET/UPDATE OLD ARRAYS USING NEW VALUES (CONCO FROM CONC) FOR ALL SPECIES
      DO ICOMP=1,NTCOMP
        DO N=1,NEQS
          CONCO(N,ICOMP) = CONC(N,ICOMP)
        ENDDO
      ENDDO
C2----COPY NEW VALUES INTO OLD VALUE ARRAYS FOR IMMOBILE DOMAIN
      IF(IDPT.GT.0) CALL GWT2DPT1AD
C3-----COMPUTE Sn AND So FOR GWF DOMAIN WHEN LAYCON IS 1 OR 3 FOR TRANSIENT FLOW SITUATIONS
      IF(ISS.EQ.0)THEN
        DO N=1,NODES
           K = LAYNOD(N)
           LC = LAYCON(K)
           IF (IBOUND(N).EQ.0) THEN 
             Sn(N) = 0.0
             So(N) = 0.0
           ELSE
             IF(LAYCON(K).EQ.1.OR.LAYCON(K).EQ.3)THEN
               IF (KSTP.EQ.1.AND.KPER.EQ.1) THEN
                 SAT = (HOLD(N) - BOT(N)) / (TOP(N) - BOT(N))
                 IF(SAT.GT.1.0) SAT = 1.0
                 IF(SAT.LT.0.0) SAT = 0.0
                 So(N) = SAT
               ELSE
                 So(N) = Sn(n)
               ENDIF
               SAT = (HNEW(N) - BOT(N)) / (TOP(N) - BOT(N))
               IF(SAT.GT.1.0) SAT = 1.0
               IF(SAT.LT.0.0) SAT = 0.0
               Sn(N) = SAT
             ENDIF
           ENDIF
        ENDDO
      ENDIF
C3B -----SET So TO Sn IF STEADY STATE FLOW
      IF(ISS.EQ.1) THEN 
        DO N=1,NEQS
          SO (N) = Sn(N)
        ENDDO  
      ENDIF
C4------COMPUTE INITIAL MASS AND WRITE TO LISTING FILE AT FIRST TIME OF FIRST STRESS PERIOD
      IF(KSTP.EQ.1.AND.KPER.EQ.1)THEN
        ISTARTMASS = 1 
        CALL INITMASS (KPER,KSTP,ISTARTMASS)
        IF(IDPT.NE.0) CALL INITMASSIM 
       IF(INCLN.NE.0) CALL INITMASSCLN 
      ENDIF
C
C5------RETURN
      RETURN
      END
C----------------------------------------------------------------------
      SUBROUTINE GWT2BCT1VELCALC
C     ******************************************************************
C     COMPUTE DARCY VELOCITY COMPONENTS FOR EACH NODE FOR DISPERSION COMPUTATIONS
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,  ONLY:IOUT,NCOL,NROW,NLAY,IBOUND,TOP,NODLAY,IVC,ISYM,
     1  ARAD,AMAT,RHS,ISSFLG,IA,JA,JAS,NJA,NODES,FAHL,BOT,NEQS,INCLN
      USE CLN1MODULE, ONLY: NCLNNDS
      USE GWFBASMODULE,ONLY:DELT
      USE GWTBCTMODULE, ONLY: MCOMP,CONC,CBCF,ICBUND,VELNOD,PRSITY,
     *   VXL,VYL
      REAL, DIMENSION(:),    ALLOCATABLE  ::XCS, YCS, ZCS
      DOUBLE PRECISION VFACE,VFACE1,VFACE2,VEX,VEY,CA,SA,ACA,ASA,
     1  CA2,SA2,ACA2,ASA2 
      REAL ANGLE
C     ------------------------------------------------------------------
C1-----FIRST VELOCITIES FOR POROUS MATRIX GRID-BLOCKS
      ALLOCATE(XCS(NODES))
      ALLOCATE(YCS(NODES))
      ALLOCATE(ZCS(NODES))
      XCS = 0.0
      YCS = 0.0
      ZCS = 0.0
      VELNOD = 0.0
C2-----LOOP OVER ALL ACTIVE NODES AND COMPUTE VELOCITY AS AN AVERAGE FROM EACH FACE FLUX
      DO N=1,NODES
        IF(ICBUND(N).EQ.0) CYCLE
C3-------GO OVER UPPER CONNECTIONS OF NODE N AND FILL Vx, Vy, Vz IN BOTH
        DO II = IA(N)+1,IA(N+1)-1
          JJ = JA(II)
          IF(JJ.GT.N.AND.JJ.LE.NODES)THEN
            IIS = JAS(II)
            IF(ICBUND(JJ).NE.0)THEN
              IF(IVC(IIS).EQ.1)THEN !VERTICAL DIRECTION CONNECTION - z flux was in negative direction
                VFACE = CBCF(IIS)/FAHL(IIS)
                VELNOD(N,3) = VELNOD(N,3) + VFACE
                VELNOD(JJ,3) = VELNOD(JJ,3) + VFACE
                ZCS(N) = ZCS(N) + 1
                ZCS(JJ) = ZCS(JJ) + 1
              ELSE                 !HORIZONTAL DIRECTION CONNECTION
                VFACE = CBCF(IIS)/FAHL(IIS)
                VFACE1 = VFACE / (TOP(N) - BOT(N))
                VFACE2 = -VFACE / (TOP(JJ) - BOT(JJ))    
C4---------------X-DIRECTION COMPONENT IS V*Cos(Angle)
                ANGLE = ARAD(IIS)   
                CA = COS(ANGLE)
                SA = SIN(ANGLE)
                ACA = CA*CA
                ASA = SA*SA
                IF(ABS(CA). LT. 1.0E-6) CA = 1.0E-6
                IF(ABS(SA). LT. 1.0E-6) SA = 1.0E-6                
c                
                ANGLE2 = ANGLE + 3.14159265358979
                CA2 = COS(ANGLE2)
                SA2 = SIN(ANGLE2)      
                ACA2 = CA2*CA2
                ASA2 = SA2*SA2
                IF(ABS(CA2). LT. 1.0E-6) CA2 = 1.0E-6
                IF(ABS(SA2). LT. 1.0E-6) SA2 = 1.0E-6                 
c                
                IF(ABS(CA).GT.1.0E-06) 
     1            VELNOD(N,1) = VELNOD(N,1) + VFACE1 / CA * ACA
                IF(ABS(CA2).GT.1.0E-06)
     1            VELNOD(JJ,1) = VELNOD(JJ,1) + VFACE2 / CA2 * ACA2
C5---------------Y-DIRECTION COMPONENT IS V/Sin(Angle)
                IF(ABS(SA).GT.1.0E-06)
     1            VELNOD(N,2) = VELNOD(N,2) + VFACE1 / SA * ASA
                IF(ABS(SA2).GT.1.0E-06)
     1            VELNOD(JJ,2) = VELNOD(JJ,2) + VFACE2 / SA2 * ASA2
C6---------------ADD UP THE DIRECTIONAL COMPONENTS FOR AVERAGING
                XCS(N) = XCS(N) + ACA
                XCS(JJ) = XCS(JJ) + ACA2
                YCS(N) = YCS(N) + ASA
                YCS(JJ) = YCS(JJ) + ASA2
              ENDIF
            ENDIF
          ENDIF
        ENDDO
      ENDDO
C7-------COMPUTE THE AVERAGE VELOCITY FOR THE NODE
      DO N=1,NODES
        IF(XCS(N).GT.1.0E-06)VELNOD(N,1) = VELNOD(N,1)/ XCS(N)
        IF(YCS(N).GT.1.0E-06)VELNOD(N,2) = VELNOD(N,2)/ YCS(N)
        IF(ZCS(N).GT.1.0E-06)VELNOD(N,3) = VELNOD(N,3)/ ZCS(N)
cc        WRITE(IOUT,*) N,velnod(N,1),velnod(N,2),XCS(N),YCS(N)
      ENDDO
      DEALLOCATE(XCS)
      DEALLOCATE(YCS)
      DEALLOCATE(ZCS)
C--------------------------------------------------------------------------------
C8-------COMPUTE CELL HORIZONTAL VELOCITY COMPONENTS AT FACES
      DO N=1,NODES
        IF(ICBUND(N).EQ.0) CYCLE
C8A---------GET HORIZV AND ITS ANGLE WITH GLOBAL AXIS
          VEX = VELNOD(N,1)
          VEY = VELNOD(N,2)
C          HORZV = (VEX**2 + VEY**2)**0.5 
C          ANGLEV = VEY / VEX
C          ANGLEV = ATAN(ANGLEV)     
C8B-------GO OVER UPPER CONNECTIONS OF NODE N AND FILL Vx, Vy, Vz FOR N AND JJ
        DO II = IA(N)+1,IA(N+1)-1
          JJ = JA(II)
          IF(JJ.GT.N.AND.JJ.LE.NODES)THEN
              
            IIS = JAS(II)   
            ANGLE = ARAD(IIS) 
C8C------------GET LOCAL COORDINATE ANGLES AND Vs            
C            ANGLE = ANGLEV - ANGLE
            CA = COS(ANGLE)
            SA = SIN(ANGLE)
            VXL(II) = VEX * CA + VEY * SA
            VYL(II) = -VEX * SA + VEY * CA
C            
            VEX2 = VELNOD(JJ,1)
            VEY2 = VELNOD(JJ,2)
            ANGLE2 = ANGLE + 3.14159265358979
            CA2 = COS(ANGLE2)
            SA2 = SIN(ANGLE2)
            VXL(ISYM(II)) = VEX2 * CA2 + VEY2 * SA2
            VYL(ISYM(II)) = -VEX2 * SA2 + VEY2 * CA2
          ENDIF
        ENDDO  
cc        write(iout,*)n,(vxl(ii),II = IA(N)+1,IA(N+1)-1)
cc        write(iout,*)n,(vyl(ii),II = IA(N)+1,IA(N+1)-1)
      ENDDO 
C------------------------------------------------------------------------------------
C8-----VELOCITIES FOR CONDUIT DOMAIN NODES -
C8----- VELNOD(1) CONTAINS CONDUIT-CONDUIT AND VELNOD(2) CONTAINS CONDUIT-MATRIX
      IF(INCLN.EQ.0) GO TO 100
      ALLOCATE(XCS(NCLNNDS))
      ALLOCATE(YCS(NCLNNDS))
      XCS = 0.0
      YCS = 0.0
ccsp      VELNOD = 0.0
C9-----LOOP OVER ALL ACTIVE CLN NODES AND COMPUTE VELOCITY AS AN AVERAGE FROM EACH FACE FLUX
      DO ICLN=1,NCLNNDS
        N = NODES + ICLN
        IF(ICBUND(N).EQ.0) CYCLE
C10-------GO OVER UPPER CONNECTIONS OF NODE N AND FILL Vx, Vy, Vz IN BOTH
        DO II = IA(N)+1,IA(N+1)-1
          JJ = JA(II)
          IF(JJ.GT.N.AND.ICBUND(JJ).NE.0)THEN
            IIS = JAS(II)
            IF(IVC(IIS).EQ.3)THEN !CONDUIT-CONDUIT CONNECTION (CAN BE MULTIPLE)
              VFACE = CBCF(IIS)/FAHL(IIS)
              VELNOD(N,1) = VELNOD(N,1) + VFACE
              VELNOD(JJ,1) = VELNOD(JJ,1) + VFACE
              XCS(N-NODES) = XCS(N-NODES) + 1
              XCS(JJ-NODES) = XCS(JJ-NODES) + 1
            ELSEIF(IVC(IIS).EQ.4)THEN   !CONDUIT-MATRIX CONNECTION (also can be multiple)
              VFACE = CBCF(IIS)/FAHL(IIS)
              VELNOD(N,2) = VELNOD(N,2) + VFACE
              VELNOD(JJ,2) = VELNOD(JJ,2) + VFACE
              YCS(N-NODES) = YCS(N-NODES) + 1
              YCS(JJ-NODES) = YCS(JJ-NODES) + 1
            ENDIF
          ENDIF
        ENDDO
      ENDDO
C11-------COMPUTE THE AVERAGE VELOCITY FOR THE NODES
      DO ICLN=1,NCLNNDS
        N = NODES + ICLN
        IF(XCS(N-NODES).GT.0.0)VELNOD(N,1) = VELNOD(N,1)/ XCS(N-NODES)
        IF(YCS(N-NODES).GT.0.0)VELNOD(N,1) = VELNOD(N,1)/ YCS(N-NODES)
      ENDDO
      DEALLOCATE(XCS)
      DEALLOCATE(YCS)
C------------------------------------------------------------------------------------
100   CONTINUE
C
C12-----RETURN
      RETURN
      END
C----------------------------------------------------------------------
      SUBROUTINE GWT2BCT1DSP(N)
C     ******************************************************************
C     FORMULATE PRINCIPAL DISPERSIVE TERM FOR EACH NODE IMPLICITLY
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,  ONLY:IOUT,NCOL,NROW,NLAY,IBOUND,TOP,NODLAY,IVC,ISYM,
     1     IA,JA,JAS,ARAD,NJA,NODES,FAHL,TOP,BOT,CL1,CL2,Sn,NEQS
      USE GWFBASMODULE,ONLY:DELT
      USE GWTBCTMODULE, ONLY: MCOMP,CONC,ICBUND,ITVD,ADMAT,
     *  VELNOD,DLX,ATXY,DLY,ATYZ,DLZ,ATXZ,PRSITY,DIFFNC,IDISPCLN
      DOUBLE PRECISION VX,VY,VZ,DXX,DYY,DZZ,DISP,ARE,THIK1,THIK2,SATN,
     *  SATJJ,SATC,DIFF,DXX2,DYY2,DZZ2,DXX1,DYY1,DZZ1,RL1,RL2,DIFFUSION,
     *  PORN,PORJJ,AVPOR,VTOT,DENOM,ACAODXX,ASAODYY
      REAL ANGLE,ACA,ASA
C     ------------------------------------------------------------------
C1-----COMPUTE Dxx, Dyy AND Dzz FOR THE CELL N
      VX = VELNOD(N,1)
      VY = VELNOD(N,2)
      VZ = VELNOD(N,3)
      VTOT = SQRT (VX*VX + VY*VY + VZ*VZ)
      SATN = Sn(N)
C
      DXX1 = 0.0
      DYY1 = 0.0
      DZZ1 = 0.0
      IF(VTOT.GT.0.0)THEN
        IF(N.LT.NODES) THEN ! MATRIX BLOCKS
          DXX1 = (DLX(N)*VX*VX + ATXY(N)*VY*VY + ATXZ(N)*VZ*VZ) / VTOT
          DYY1 = (ATXY(N)*VX*VX + DLY(N)*VY*VY + ATYZ(N)*VZ*VZ) / VTOT
          DZZ1 = (ATXZ(N)*VX*VX + ATYZ(N)*VY*VY + DLZ(N)*VZ*VZ) / VTOT
        ELSE      ! CONDUIT CONNECTIONS
          DXX1 = DLX(N)*VX*VX / VTOT
          DYY1 = ATXY(N)*VY*VY / VTOT
        ENDIF
      ENDIF
C2-----LOOP OVER CONNECTIONS OF NODE N AND FILL
      DO II = IA(N)+1,IA(N+1)-1
        ICMC = 0   ! INDICATES A CONDUIT TO MATRIX CONNECTION
        JJ = JA(II)
        IF(JJ.LE.N) CYCLE !ONLY FOR UPPER TRIANGLE
        IF(ICBUND(JJ).NE.0)THEN
          IIS = JAS(II)
C3-----COMPUTE Dxx, Dyy AND Dzz FOR THE CELL JJ
          VX = VELNOD(JJ,1)
          VY = VELNOD(JJ,2)
          VZ = VELNOD(JJ,3)
          VTOT = SQRT (VX*VX + VY*VY + VZ*VZ)
C
          DXX2 = 0.0
          DYY2 = 0.0
          DZZ2 = 0.0
          IF(VTOT.GT.0.0)THEN
            IF(JJ.LT.NODES) THEN ! MATRIX BLOCKS
              DXX2 =(DLX(JJ)*VX*VX +ATXY(JJ)*VY*VY +ATXZ(JJ)*VZ*VZ)/VTOT
              DYY2 =(ATXY(JJ)*VX*VX +DLY(JJ)*VY*VY +ATYZ(JJ)*VZ*VZ)/VTOT
              DZZ2 =(ATXZ(JJ)*VX*VX +ATYZ(JJ)*VY*VY +DLZ(JJ)*VZ*VZ)/VTOT
            ELSE      ! CONDUIT CONNECTIONS
              DXX2 = DLX(JJ)*VX*VX / VTOT
              DYY2 = ATXY(JJ)*VY*VY / VTOT
            ENDIF
          ENDIF
C4-----COMPUTE Dxx, Dyy AND Dzz FOR FACE BETWEEN CELLS N AND JJ FOR GW NODES
          IF(N.LE.NODES.AND.JJ.LE.NODES)THEN
            DXX = (DXX1 + DXX2)*0.5
            DYY = (DYY1 + DYY2)*0.5
            DZZ = (DZZ1 + DZZ2)*0.5
          ENDIF
C5-----COMPUTE THE DISPERSION CONDUCTANCE TERM BETWEEN CELLS N AND JJ
          IF(IVC(IIS).EQ.1)THEN !VERTICAL DIRECTION CONNECTION
            ARE = FAHL(IIS)
            DISP = DZZ
          ELSEIF(JJ.GT.NODES.OR.N.GT.NODES)THEN !CONDUIT-MATRIX OR CONDUIT-CONDUIT CONNECTION
            ARE = FAHL(IIS)
            IF(N.GT.NODES.AND.JJ.GT.NODES)THEN !CONDUIT-CONDUIT
              DXX = (DXX1 + DXX2)*0.5
              DISP = DXX
            ELSE                               !CONDUIT-MATRIX
                DYY = DYY2  ! ONLY THE CONDUIT HAS DISPERSIVITY AND SHOULD NOT BE MIXED WITH MATRIX
                DISP = DYY
              ICMC = 1     !INDEX FOR CONDUIT-MATRIX CONNECTION IS ON
            ENDIF
          ELSE                                 !HORIZONTAL DIRECTION CONNECTION
            THIK1 = (TOP(N) - BOT(N))
            THIK2 = (TOP(JJ) - BOT(JJ))
            ANGLE = ARAD(IIS)
            ARE =  (THIK1 + THIK2)*0.5 * FAHL(IIS)
C
C6-------------FOR FACE OF LAYER WITH STACKED LAYERS AREA IS MINIMUM THICKNESS
c            ACA = ABS(COS(ANGLE))
c            ASA = ABS(SIN(ANGLE))
c            DISP = DXX * ACA + DYY*ASA
C---------------------------------------------------------------------------
            ACA = COS(ANGLE) * COS(ANGLE)
            ASA = SIN(ANGLE) * SIN(ANGLE)
CSP-------------USE WEIGHTED ARITHMITIC MEAN ROTATION
            DISP = DXX*ACA + DYY*ASA
CSP--------------OR, PERFORM HARMONIC MEAN ROTATION LIKE FOR ANISOTROPIC K TERM
cc          IF(DXX. LT. 1.0E-20. AND. DYY. LT. 1.0E-20) THEN
cc            DISP = 0.0
cc          ELSE
cc            IF(DXX.LT.1.0E-20) DXX = 1.0E-20
cc            IF(DYY.LT.1.0E-20) DYY = 1.0E-20
cc            DISP = 1.0 / (ACA/DXX + ASA/DYY)
cc          ENDIF
C-----------GET RID OF SINGULARITIES
CSP            IF(DXX. LT. 1.0E-10. AND. DYY. LT. 1.0E-10) THEN
CSP              DISP = 0.0
CSP            ELSE
CSP              IF(DXX.LT.1.0E-10) THEN
CSP                  ACAODXX = 1.0E20
CSP                  IF(ACA.LT.1.0E-10) ACAODXX = 0.0
CSP              ELSE
CSP                  ACAODXX = ACA/DXX
CSP              ENDIF
CSP              IF(DYY.LT.1.0E-10) THEN
CSP                  ASAODYY = 1.0E20
CSP                  IF(ASA.LT.1.0E-10) ASAODYY = 0.0
CSP              ELSE
CSP                  ASAODYY = ASA/DYY
CSP              ENDIF
CSP              DISP = 1.0 / (ACAODXX + ASAODYY)
CSP            ENDIF
C-------------------------------------------------------------------
          ENDIF
C
C7---------ADD DIFFUSION TERM
          DIFFUSION = DIFFNC
          PORN = 1.0
          PORJJ = 1.0
          IF(N.LE.NODES) PORN = PRSITY(N)
          IF(JJ.LE.NODES) PORJJ = PRSITY(JJ)
          AVPOR =  0.5 * (PORN + PORJJ)
          DIFFUSION = DIFFUSION * AVPOR
          DISP = DISP + DIFFUSION
C-----------USE MINIMUM OF TWO TO LIMIT DIFFUSION TO DRY NODES
          SATC = 0.0
          SATJJ = Sn(JJ)
          SATC = SATJJ
          IF(SATC. GT. SATN) THEN
            SATC = SATN
          ENDIF
C-----------NET DISPERSION SCALED BY SATURATION
          DISP = SATC * DISP
C
C8---------CONVERT DISPERSION  TO DISPERSION CONDUCTANCE
          IF(IDISPCLN.EQ.1. AND. ICMC.EQ.1) THEN
            FRAD = CL1(IIS)
            RO = CL2(IIS)
            DISP = DISP * FAHL(IIS) / FRAD /LOG(RO/FRAD)
          ELSE
            DISP = DISP * ARE / (CL1(IIS) + CL2(IIS))
          ENDIF
C
C9---------ADD TERM TO MATRIX FOR THIS NODE AND ITS EFFECT ON CONNECTING NODE JJ
          ADMAT(II) = ADMAT(II) + DISP
          ADMAT(IA(N)) = ADMAT(IA(N)) - DISP
          ADMAT(ISYM(II)) = ADMAT(ISYM(II)) + DISP
          ADMAT(IA(JJ)) = ADMAT(IA(JJ)) - DISP
        ENDIF
      ENDDO
C
C10-----RETURN
      RETURN
      END
C----------------------------------------------------------------------
C----------------------------------------------------------------------
      SUBROUTINE GWT2BCT1DC(ICOMP)
C     ******************************************************************
C     COMPUTE CONCENTRATION GRADIENTS IN X-, Y-, AND Z-DIRECTIONS
C     FOR RHS CROSS-DISPERSION COMPUTATION
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,  ONLY:IOUT,NCOL,NROW,NLAY,IBOUND,TOP,NODLAY,IVC,ISYM,
     1    ARAD,AMAT,RHS,ISSFLG,IA,JA,JAS,NJA,NODES,FAHL,BOT,NEQS,INCLN,
     1    CL1,CL2
      USE CLN1MODULE, ONLY: NCLNNDS
      USE GWFBASMODULE,ONLY:DELT
      USE GWTBCTMODULE, ONLY: MCOMP,CONC,ICBUND,PRSITY,DXCS,DYCS,DZCS,
     1    DCXL,DCYL
      REAL, DIMENSION(:),    ALLOCATABLE  ::XCS, YCS, ZCS
      REAL ANGLE
      REAL*8 DEX, DEY, DCDL,CA,SA,ACA,ASA,CA2,SA2,ACA2,ASA2
C     ------------------------------------------------------------------
C1-----VELOCITIES FOR POROUS MATRIX GRID-BLOCKS
      ALLOCATE(XCS(NODES))
      ALLOCATE(YCS(NODES))
      ALLOCATE(ZCS(NODES))
      DO N=1,NODES
        XCS(N) = 0.0
        YCS(N) = 0.0
        ZCS(N) = 0.0
        DXCS(N) = 0.0
        DYCS(N) = 0.0
        DZCS(N) = 0.0
      ENDDO
C2-----LOOP OVER ALL ACTIVE NODES AND COMPUTE VELOCITY AS AN AVERAGE FROM EACH FACE FLUX
      DO N=1,NODES
        IF(ICBUND(N).EQ.0) CYCLE
C3-------GO OVER UPPER CONNECTIONS OF NODE N AND FILL Vx, Vy, Vz IN BOTH
        DO II = IA(N)+1,IA(N+1)-1
          JJ = JA(II)
          IF(JJ.GT.N.AND.JJ.LE.NODES)THEN
            IIS = JAS(II)
            IF(IVC(IIS).EQ.1)THEN !VERTICAL DIRECTION CONNECTION
              EL1 = TOP(N) - BOT(N) 
              EL2 = TOP(JJ) - BOT(JJ)
            ELSE  
              EL1 = CL1(IIS)
              EL2 = CL2(IIS)
            ENDIF  
            DCDL = 0.0E0
            IIS = JAS(II)
            IF(ICBUND(JJ).NE.0)THEN
              DCDL = (CONC(JJ,ICOMP) - CONC(N,ICOMP))/(EL1 + EL2)     
              IF(IVC(IIS).EQ.1)THEN !VERTICAL DIRECTION CONNECTION
                DZCS(N) = DZCS(N) + DCDL
                DZCS(JJ) = DZCS(JJ) + DCDL
                ZCS(N) = ZCS(N) + 1
                ZCS(JJ) = ZCS(JJ) + 1
              ELSE                 !HORIZONTAL DIRECTION CONNECTION
C4---------------X-DIRECTION COMPONENT IS V*Cos(Angle)
                ANGLE = ARAD(IIS)
                CA = COS(ANGLE)
                SA = SIN(ANGLE)
                ACA = CA*CA
                ASA = SA*SA
                IF(ABS(CA). LT. 1.0E-6) CA = 1.0E-6
                IF(ABS(SA). LT. 1.0E-6) SA = 1.0E-6                   
c
                ANGLE2 = ANGLE + 3.14159265358979 
                CA2 = COS(ANGLE2)
                SA2 = SIN(ANGLE2)
                ACA2 = CA2*CA2
                ASA2 = SA2*SA2
                IF(ABS(CA2). LT. 1.0E-6) CA2 = 1.0E-6
                IF(ABS(SA2). LT. 1.0E-6) SA2 = 1.0E-6                   
                DCDL2 = - DCDL
c                
                IF(ABS(CA).GT.1.0E-06)
     1            DXCS(N) = DXCS(N) + DCDL / CA * ACA
                IF(ABS(CA2).GT.1.0E-06)
     1            DXCS(JJ) = DXCS(JJ) + DCDL2 / CA2 * ACA2   
C5---------------Y-DIRECTION COMPONENT IS V/Sin(Angle)
                IF(ABS(SA).GT.1.0E-06)
     1            DYCS(N) = DYCS(N) + DCDL / SA * ASA
                IF(ABS(SA2).GT.1.0E-06)
     1            DYCS(JJ) = DYCS(JJ) + DCDL2 / SA2 * ASA2    
C6---------------ADD UP THE DIRECTIONAL COMPONENTS FOR AVERAGING
                XCS(N) = XCS(N) + ACA
                XCS(JJ) = XCS(JJ) + ACA2
                YCS(N) = YCS(N) + ASA
                YCS(JJ) = YCS(JJ) + ASA2
              ENDIF
            ENDIF
          ENDIF
        ENDDO
      ENDDO
C7-------COMPUTE THE AVERAGE GRADIENTS DCDX, DCDY, DCDZ FOR THE NODE
      DO N=1,NODES
        IF(XCS(N).GT.1.0E-06)DXCS(N) = DXCS(N)/ XCS(N)
        IF(YCS(N).GT.1.0E-06)DYCS(N) = DYCS(N)/ YCS(N)
        IF(ZCS(N).GT.1.0E-06)DZCS(N) = DZCS(N)/ ZCS(N)
cc        WRITE(IOUT,*)'nodal dcdx and dcdy', N,DXCS(N),DYCS(N)
      ENDDO
      DEALLOCATE(XCS)
      DEALLOCATE(YCS)
      DEALLOCATE(ZCS)
C8A------COMPUTE CELL HORIZONTAL dC COMPONENTS AT FACES
      DO N=1,NODES
        IF(ICBUND(N).EQ.0) CYCLE     
        DEX = DXCS(N)
        DEY = DYCS(N)
C8B-------GO OVER ALL CONNECTIONS OF NODE N AND FILL Vx, Vy, Vz
        DO II = IA(N)+1,IA(N+1)-1
          JJ = JA(II)
          IF(JJ.GT.N.AND.JJ.LE.NODES)THEN
            IIS = JAS(II)   
            ANGLE = ARAD(IIS) 
C8C------------GET LOCAL COORDINATE ANGLES AND LOCAL DX AND DY           
            CA = COS(ANGLE)
            SA = SIN(ANGLE)            
            DCXL(II) = DEX * CA + DEY * SA
            DCYL(II) = -DEX * SA + DEY * CA
C            
            DEX2 = DXCS(JJ)
            DEY2 = DYCS(JJ)
            ANGLE2 = ARAD(IIS) + 3.14159265358979
            CA2 = COS(ANGLE2)
            SA2 = SIN(ANGLE2)
            DCXL(ISYM(II)) = DEX2 * CA2 + DEY2 * SA2
            DCYL(ISYM(II)) = -DEX2 * SA2 + DEY2 * CA2
          ENDIF
        ENDDO  
      ENDDO  
C     
cc      write(iout,*)'dcdx then dcdy on all faces of a cell'      
cc      DO N=1,NODES   
cc        write(iout,*)n,(DCXL(ii),II = IA(N)+1,IA(N+1)-1)
cc        write(iout,*)n,(DCYL(ii),II = IA(N)+1,IA(N+1)-1) 
cc      ENDDO        
C
C8------RETURN
      RETURN
      END
C------------------------------------------------------------------------
      SUBROUTINE GWT2BCT1DSPX(N)
C     ******************************************************************
C     FORMULATE CROSS-DISPERSIVE TERM FOR EACH NODE ON RHS
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,  ONLY:IOUT,NCOL,NROW,NLAY,IBOUND,TOP,NODLAY,IVC,ISYM,
     1     RHS,IA,JA,JAS,NJA,ARAD,NODES,FAHL,TOP,BOT,CL1,CL2,Sn
      USE GWFBASMODULE,ONLY:DELT
      USE GWTBCTMODULE, ONLY: MCOMP,CONC,CBCF,ICBUND,ITVD,
     *  VELNOD,DLX,ATXY,DLY,ATYZ,DLZ,ATXZ,PRSITY,DIFFNC,DXCS,DYCS,DZCS,
     1  ALXY,ALYZ,ALXZ,VXL,VYL,DCXL,DCYL
      DOUBLE PRECISION VX,VY,VZ,DXY1,DYZ1,DXZ1,DYX1,DZY1,DZX1,AREA,
     *  THIK1,THIK2,SATN,SATJJ,SATC,DCDX1,DCDY1,DCDZ1,DCDX2,DCDY2,DCDZ2,
     *  DZ,DL,TERM1,TERM2,SA,CA,EL,DXY2,DYZ2,DXZ2,DYX2,DZY2,DZX2,VTOT,
     *  dxy,dcdy,dyx,dcdx,dxz,dcdz,dyz,dzy,dzx,disp,disp2 
      REAL ANGLE,ACA,ASA
      DOUBLE PRECISION DCX1,DCX2,DCY1,DCY2
C
C     ------------------------------------------------------------------
C
       SATN = Sn(N) 
C1-----GO OVER ALL CONNECTIONS OF CELL N AND FILL
C4-----Dxy, Dyz AND Dxz FOR N AND JJ ALONG LOCAL DIRECTION      
      DO II = IA(N)+1,IA(N+1)-1
        JJ = JA(II)
C3A-----ONLY GO OVER UPPER TRIANGLE OF MATRIX (FORWARD FACES)
C3A-----ALSO NO CROSS DISPERSION COMPONENT ON CLN CELLS
        IF(JJ.LT.N.OR.JJ.GT.NODES) CYCLE
        IF(ICBUND(JJ).NE.0)THEN
          IIS = JAS(II)
          ANGLE = ARAD(IIS)
          CA = (COS(ANGLE))
          SA = (SIN (ANGLE))
          DISP = 0.0
CD2          DISP2 = 0.0
C------------------------------------------------------------------
C4-------COMPUTE VELOCITY AND DC FOR CELLS N AND JJ ALONG LOCAL COORD DIRECTION
          VX1 = VXL(II)
          VY1 = VYL(II)
          VZ1 = (VELNOD(N,3))
          VTOT1 = SQRT (VX1*VX1 + VY1*VY1 + VZ1*VZ1)
          VX2 = -VXL(ISYM(II))
          VY2 = -VYL(ISYM(II))
          VZ2 = (VELNOD(JJ,3))
          VTOT2 = SQRT (VX2*VX2 + VY2*VY2 + VZ2*VZ2)    
          DCX1 = DCXL(II)
          DCY1 = DCYL(II)
          DCX2 = -DCXL(ISYM(II))
          DCY2 = -DCYL(ISYM(II))
C----------CROSS D TERMS ALONG LOCAL COORD DIRECTION
          DXY1 = 0.0
          DYZ1 = 0.0
          DXZ1 = 0.0
          IF(VTOT1.GT.0.0)THEN
            DXY1 = (ALXY(N) - ATXY(N))*VX1*VY1 / VTOT1
            DYZ1 = (ALYZ(N) - ATYZ(N))*VY1*VZ1 / VTOT1
            DXZ1 = (ALXZ(N) - ATXZ(N))*VX1*VZ1 / VTOT1
          ENDIF
          DYX1 = DXY1
          DZY1 = DYZ1
          DZX1 = DXZ1
C                   
          DXY2 = 0.0
          DYZ2 = 0.0
          DXZ2 = 0.0
          IF(VTOT2.GT.0.0)THEN
            DXY2 = (ALXY(JJ) - ATXY(JJ))*VX2*VY2 / VTOT2
            DYZ2 = (ALYZ(JJ) - ATYZ(JJ))*VY2*VZ2 / VTOT2
            DXZ2 = (ALXZ(JJ) - ATXZ(JJ))*VX2*VZ2 / VTOT2
          ENDIF
          DYX2 = DXY2
          DZY2 = DYZ2
          DZX2 = DXZ2
C------------------------------------------------------------------
C5-------FILL DCDX, DCDY AND DCDZ ALONG LOCAL COORD DIRECTIONS
C          DCDX1 = DXCS(N) * CA + DYCS(N) * SA
C          DCDY1 = -DXCS(N) * SA + DYCS(N) * CA
C          DCDZ1 = DZCS(N)
C          THIK1 = TOP(N) - BOT(N)
C          DCDX2 = DXCS(JJ) * CA + DYCS(JJ) * SA
C          DCDY2 = -DXCS(JJ) * SA + DYCS(JJ) * CA
C          DCDZ2 = DZCS(JJ)
C          THIK2 = TOP(JJ) - BOT(JJ)          
C------------------------------------------------------------------
C5-------FILL DCDX, DCDY AND DCDZ ALONG LOCAL COORD DIRECTIONS
          DCDX1 = DCX1
          DCDY1 = DCY1
          DCDZ1 = DZCS(N)
          THIK1 = TOP(N) - BOT(N)
          DCDX2 = DCX2
          DCDY2 = DCY2
          DCDZ2 = DZCS(JJ)
          THIK2 = TOP(JJ) - BOT(JJ)                    
C------------------------------------------------------------------
C6------COMPUTE DISPERSION TERM FOR CONNECTION BEING IN VERTICAL DIRECTION
          IF(IVC(IIS).EQ.1)THEN !VERTICAL DIRECTION CONNECTION
C------------------------------------------------------------------
C6A------d/dz (Dzx dC/dx)
             DZX = (DZX1 + DZX2)*0.5
             DCDX = (DCDX1 + DCDX2) * 0.5
             DISP = DISP + DZX * DCDX 
CD2             DISP = DISP + DZX1 * DCDX1
CD2             DISP2 = DISP2 + DZX2 * DCDX2 
C------------------------------------------------------------------
C6B------d/dz (Dzy dC/dy)
             DZY = (DZY1 + DZY2)*0.5
             DCDY = (DCDY1 + DCDY2) * 0.5
             DISP = DISP + DZY * DCDY
CD2             DISP = DISP + DZY1 * DCDY1
CD2             DISP2 = DISP2 + DZY2 * DCDY2
C
C6C------CONVERT DISPERSION  TO DISPERSION CONDUCTANCE
            AREA = FAHL(IIS)
            DISP = DISP * AREA
CD2            DISP2 = DISP2 * AREA
C------------------------------------------------------------------
C7------COMPUTE DISPERSION TERM IN HORIZONTAL DIRECTION
          ELSE                  !HORIZONTAL DIRECTION CONNECTION
C------------------------------------------------------------------
C 
C7A------d/dx (Dxy dC/dy)
             DXY = (DXY1 + DXY2)*0.5
             DCDY = (DCDY1 + DCDY2) * 0.5
             DISP = DISP + DXY * DCDY 
CD2             DISP = DISP + DXY1 * DCDY1  
CD2             DISP2 = DISP2 + DXY2 * DCDY2 
C7C------d/dx (Dxz dC/dz)
             DXZ = (DXZ1 + DXZ2)*0.5
             DCDZ = (DCDZ1 + DCDZ2) * 0.5
             DISP = DISP + DXZ * DCDZ
CD2             DISP = DISP + DXZ1 * DCDZ1
CD2             DISP2 = DISP2 + DXZ2 * DCDZ2
C 
C7E------CONVERT DISPERSION  TO DISPERSION CONDUCTANCE
            THIK1 = (TOP(N) - BOT(N))
            THIK2 = (TOP(JJ) - BOT(JJ))
            AREA =  (THIK1 + THIK2)*0.5 * FAHL(IIS)
            DISP = DISP * AREA
CD2            DISP2 = DISP2 * AREA
          ENDIF
C---------------------------------------------------------------------------
C8-----------USE MINIMUM OF TWO SATURATIONS TO LIMIT DISPERSION TO DRY NODES
          SATC = 0.0
          SATJJ = Sn(JJ)
          SATC = SATJJ
          IF(SATC. GT. SATN) THEN
            SATC = SATN
          ENDIF
C 
C-----------NET DISPERSION SCALED
          DISP = SATC * DISP 
CD2          DISP2 = SATC * DISP2 
C-----------------------------------------------------------------------------------
C8---------ADD TERM TO RHS VECTOR FOR THIS NODE AND ITS EFFECT ON CONNECTING NODE JJ
          RHS(N) = RHS(N) - DISP
          RHS(JJ) = RHS(JJ) + DISP 
CD2          RHS(JJ) = RHS(JJ) + DISP2         
        ENDIF
      ENDDO
C
C9------RETURN
      RETURN
      END
C----------------------------------------------------------------------
      SUBROUTINE GWT2BCT1ADV(N)
C     ******************************************************************
C     FORMULATE ADVECTIVE TERM FOR EACH NODE AND COMPONENT
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,  ONLY:IOUT,NCOL,NROW,NLAY,IBOUND,TOP,NODLAY,ISYM,
     1     AMAT,RHS,ISSFLG,IA,JA,JAS,NJA,NODES
      USE GWFBASMODULE,ONLY:DELT
      USE GWTBCTMODULE, ONLY: MCOMP,CONC,CBCF,ICBUND,ITVD,ADMAT
      DOUBLE PRECISION QIJ
C     ------------------------------------------------------------------
C1-----GO OVER CONNECTIONS OF NODE N AND FILL
      DO II = IA(N)+1,IA(N+1)-1
        JJ = JA(II)
        IF(JJ.GT.N.AND.ICBUND(JJ).NE.0)THEN
        IIS = JAS(II)
C----------------------------------------------------------------------
C2---------FILL UPSTREAM TERM
          QIJ = CBCF(IIS)
          IF(QIJ.GT.0)THEN !N IS UPSTREAM OF CELLS N AND JJ
            ADMAT(IA(N)) = ADMAT(IA(N)) - QIJ
            ADMAT(ISYM(II)) = ADMAT(ISYM(II)) + QIJ
          ELSE ! JJ IS UPSTREAM OF CELLS N AND JJ
            ADMAT(II) = ADMAT(II) - QIJ
            ADMAT(IA(JJ)) = ADMAT(IA(JJ)) + QIJ
          ENDIF
C-----------------------------------------------------------------------
        ENDIF
      ENDDO
C
C3------RETURN
      RETURN
      END
C----------------------------------------------------------------------
      SUBROUTINE GWT2BCT1TVD(N,ICOMP)
C     ******************************************************************
C     FILL TVD TERM ON RHS
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,  ONLY:IOUT,NCOL,NROW,NLAY,IBOUND,TOP,NODLAY,ISYM,
     1     RHS,ISSFLG,IA,JA,JAS,NJA,NODES
      USE GWFBASMODULE,ONLY:DELT
      USE GWTBCTMODULE, ONLY: MCOMP,CONC,CBCF,ICBUND,ITVD
      DOUBLE PRECISION QIJ
C     ------------------------------------------------------------------
C1-----GO OVER CONNECTIONS OF NODE N AND FILL
      DO II = IA(N)+1,IA(N+1)-1
        JJ = JA(II)
        IF(JJ.GT.N.AND.ICBUND(JJ).NE.0)THEN
          IIS = JAS(II)
          QIJ = CBCF(IIS)
          CALL SGWT2BCT1TVDS(N,ICOMP,II,JJ,QIJ)
        ENDIF
      ENDDO
C
C2------RETURN
      RETURN
      END
C----------------------------------------------------------------------
      SUBROUTINE SGWT2BCT1TVDS(N,ICOMP,II,JJ,QIJ)
C     ******************************************************************
C     FORMULATE TVD TERM FOR EACH NODE AND COMPONENT
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,  ONLY:IOUT,NCOL,NROW,NLAY,IBOUND,TOP,NODLAY,ISYM,
     1     AMAT,RHS,ISSFLG,IA,JA,JAS,NJA,NODES,CL1,CL2
      USE GWFBASMODULE,ONLY:DELT
      USE GWTBCTMODULE, ONLY: MCOMP,CONC,CBCF,ICBUND,ITVD,
     1   CONCO,TIMEWEIGHT
      DOUBLE PRECISION QIJ,Q2,ELUPDN,ELUP2UP,ALIMITER,SMOOTH,CDIFF,
     1   CDN,CUP,C2UP,TW,TWM1
C1---------------FIND UPSTREAM NODE
      IIS = JAS(II)
      IF(QIJ.GT.0)THEN !N IS UPSTREAM OF CELLS N AND JJ
        IUP = N
        IDN = JJ
      ELSE ! JJ IS UPSTREAM OF CELLS N AND JJ
        IUP = JJ
        IDN = N
      ENDIF
      ELUPDN = (CL1(IIS) + CL2(IIS))
C2-----------FIND SECOND POINT UPSTREAM TO POINT IUP
      I2UP = 0
      Q2 = 0.0
      DO I2 = IA(IUP)+1,IA(IUP+1)-1
        J2 = JA(I2)
        IF(ICBUND(J2).NE.0)THEN
          I2S = JAS(I2)
          IF(J2.GT.IUP)THEN   ! J2 IS HIGHER NODE NUMBER - Q INTO IUP IS NEGATIVE
            IF(-CBCF(I2S).GT.Q2) THEN
              Q2 = -CBCF(I2S)
              I2UP = J2
              ELUP2UP = (CL1(I2S) + CL2(I2S))
            ENDIF
          ELSEIF(IUP.GT.J2)THEN !J2 IS LOWER NODE NUMBER
            IF(CBCF(I2S).GT.Q2) THEN
              Q2 = CBCF(I2S)
              I2UP = J2
              ELUP2UP = (CL1(I2S) + CL2(I2S))
            ENDIF
          ENDIF
        ENDIF
      ENDDO
C3---------------COMPUTE SMOOTHNESS SENSOR AND LIMITER TERM
      IF(I2UP.NE.0)THEN
        CDN = CONC(IDN,ICOMP)
        CUP = CONC(IUP,ICOMP)
        C2UP = CONC(I2UP,ICOMP)
C---------------FOR TIME WEIGHTING
        IF(TIMEWEIGHT.LT.0.99999)THEN
          TW = TIMEWEIGHT
          TWM1 = 1.0 - TW
          CDN = (TW * CDN + TWM1 * CONCO(IDN,ICOMP))
          CUP = (TW * CUP + TWM1 * CONCO(IUP,ICOMP))
          C2UP = (TW * C2UP + TWM1 * CONCO(I2UP,ICOMP))
        ENDIF
C--------------------------------------
        SMOOTH = 0.0
        CDIFF = ABS(CDN - CUP)
        IF(CDIFF.GT.1.0E-10)
     *  SMOOTH = (CUP - C2UP) /ELUP2UP * ELUPDN / (CDN - CUP)
        ALIMITER = 0.0
        IF(SMOOTH.GT.0.0) ALIMITER = 2. * SMOOTH / (1. + SMOOTH)
C
        RHS(N) = RHS(N) + 0.5*ALIMITER*QIJ * (CDN - CUP)
        RHS(JJ) = RHS(JJ) - 0.5*ALIMITER*QIJ * (CDN - CUP)
      ENDIF
C
C4------RETURN
      RETURN
      END
C----------------------------------------------------------------------
      SUBROUTINE GWT2BCT1TWA(N,ICOMP,TWM1)
C     ******************************************************************
C     TIME WEIGHTING ADJUSTED FOR USER-DEFINED WEIGHT (DEFAULT IS FULLY IMPLICIT)
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,  ONLY:IOUT,NCOL,NROW,NLAY,IBOUND,TOP,NODLAY,ISYM,
     1     AMAT,RHS,ISSFLG,IA,JA,JAS,NJA,NODES
      USE GWFBASMODULE,ONLY:DELT
      USE GWTBCTMODULE, ONLY: MCOMP,CONC,CONCO,CBCF,ICBUND,ITVD,
     1  TIMEWEIGHT
      DOUBLE PRECISION QIJ,TWM1
C     ------------------------------------------------------------------
C1-----GO OVER CONNECTIONS OF NODE N AND FILL
      DO II = IA(N),IA(N+1)-1
        JJ = JA(II)
        IF(ICBUND(JJ).NE.0)THEN
          RHS(N) = RHS(N) - AMAT(II)*TWM1 * CONCO(JJ,ICOMP)
        ENDIF
      ENDDO
C
C2------RETURN
      RETURN
      END
C----------------------------------------------------------------------
      SUBROUTINE GWT2BCT1STO(N,ICOMP,DTERMS,RTERMS,ISS)

C     ******************************************************************
C     FORMULATE STORAGE TERM FOR EACH NODE AND COMPONENT
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,   ONLY:IOUT,NCOL,NROW,NLAY,IBOUND,TOP,NODLAY,ISYM,
     1 HNEW,AMAT,RHS,HOLD,ISSFLG,IA,JA,NJA,NODES,AREA,BOT,TOP,Sn,So
      USE CLN1MODULE, ONLY: ACLNNDS
      USE GWFBASMODULE,ONLY:DELT
      USE GWTBCTMODULE, ONLY: MCOMP,CONC,CONCO,ICT,IADSORB,PRSITY,
     * ADSORB,FLICH,IZOD,IFOD,ZODRW,FODRW,ZODRS,FODRS,MCOMPT,IHEAT
      DOUBLE PRECISION VODT,ADSTERM,CW,FL,CWO,DT,RT,ALENG,VOLU,
     1  X,Y,EPSS,EPS,CT,CEPS,QA,QEPS,DQ,OMP,DTERMS,RTERMS
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
C2-------STORAGE TERMS ON SOIL
        IFILLADS = 0
        IF(IHEAT.EQ.1.AND.ICOMP.EQ.MCOMPT) IFILLADS = 1
        IF(N.LE.NODES.AND.IADSORB.EQ.1) IFILLADS = 1
        IF(N.LE.NODES.AND.IADSORB.EQ.2) IFILLADS = 2
        IF(N.LE.NODES.AND.IADSORB.EQ.3) IFILLADS = 3
        IF(N.GT.NODES) IFILLADS = 0
        IF(IFILLADS.EQ.1)THEN
C---------LINEAR ISOTHERM
          ADSTERM = ADSORB(N,ICOMP) * VODT
          DTERMS = DTERMS - ADSTERM
          RTERMS = RTERMS - ADSTERM * CONCO(N,ICOMP)
C-----------------------------------------------------------------------
        ELSEIF(IFILLADS.EQ.2.OR.IFILLADS.EQ.3)THEN
C3---------NONLINEAR FREUNDLICH ISOTHERM FILLED AS NEWTON
          ADSTERM = ADSORB(N,ICOMP) * VODT
          FL = FLICH(N,ICOMP)
          CW = CONC(N,ICOMP)
          IF(CW.LT.0.0) CW = 0.0
          CWO = CONCO(N,ICOMP)
          IF(CWO.LT.0.0) CWO = 0.0
          DT = 0.0
          RT = 0.0
          IF(IFILLADS.EQ.2) THEN ! FREUNDLICH FILLED AS NEWTON
            RT = CW**FL - CWO**FL 
            RT = RT * ADSTERM 
            if( fl.lt.1.0 .and. cw. lt. 1.e-7) cw = 1e-7        !No NWT term if/when slope approaches infinity
            DT = ADSTERM*FL * CW**(FL - 1.0)           
          ELSEIF(IFILLADS.EQ.3) THEN ! LANGMUIR FILLED AS NEWTON
            RT = CW/(1.0 + FL*CW) - CWO/(1.0 + FL*CWO)
            RT = RT *ADSTERM
            DT = (1.0+FL*CW)*ADSTERM - ADSTERM*CW*FL ! ANALYTIC DERIV OF U/V = (VdU - UdV)/V^2
            DT = DT / (1.0 + FL*CW)**2
          ENDIF    
          DTERMS = DTERMS - DT
          RTERMS = RTERMS - DT * CW + RT
        ENDIF
C-----------------------------------------------------------------------
C4-------STORAGE TERM IN WATER
        CALL GWT2BCT1STOW (N,ICOMP,DTERMS,RTERMS,VODT,VOLU,ALENG,ISS)
C-----------------------------------------------------------------------
      ELSE       !-----------------------TOTAL CONCENTRATION FORMULATION
C5-------NET STORAGE TERM FOR TOTAL CONCENTRATION FORMULATION
        DTERMS = DTERMS - VODT
        RTERMS = RTERMS - VODT * CONCO(N,ICOMP)
      ENDIF
C
C9------RETURN
      RETURN
      END
C---------------------------------------------------------------------------
      SUBROUTINE GWT2BCT1STOW(N,ICOMP,DTERMS,RTERMS,VODT,VOLU,ALENG,ISS)
C     ******************************************************************
C     FORMULATE STORAGE TERM FOR EACH NODE AND COMPONENT IN WATER
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,   ONLY:IOUT,NCOL,NROW,NLAY,IBOUND,TOP,NODLAY,ISYM,
     1 HNEW,AMAT,RHS,HOLD,ISSFLG,IA,JA,NJA,NODES,AREA,BOT,TOP,Sn,So,
     2 LAYNOD
      USE CLN1MODULE, ONLY: ACLNNDS
      USE GWFBASMODULE,ONLY:DELT
      USE GWFBCFMODULE,ONLY:LAYCON,SC1,SC2,IALTSTO
      USE GWTBCTMODULE, ONLY: MCOMP,CONC,CONCO,ICT,IADSORB,PRSITY,
     * ADSORB,FLICH,IZOD,IFOD,ZODRW,FODRW,
     1  ZODRS,FODRS
      DOUBLE PRECISION VODT,ADSTERM,CW,FL,CWO,DT,RT,ALENG,VOLU,
     1  X,Y,EPSS,EPS,CT,CEPS,QA,QEPS,DQ,OMP,DTERMS,RTERMS,
     1  H1,H2,TP,SATN,SATO
C-------------------------------------------------------------------------------------
C
      SATN = Sn(N)
      SATO = So(N)
      if(satn.lt.1.0e-3) satn = 1.0e-3
      if(sato.lt.1.0e-3) sato = 1.0e-3
      IF(ISS.EQ.1.OR.N.GT.NODES) THEN
        IF(ISS.EQ.1) SATO = SATN
C1-------USE TRADITIONAL STORAGE TERM FOR CLN DOMAIN OR SS FLOW
        IF(N.LE.NODES) VODT = VODT * PRSITY(N)
        DTERMS = DTERMS - VODT * SATN
        RTERMS = RTERMS - VODT * SATO * CONCO(N,ICOMP)
      ELSE
        K = LAYNOD(N)
        LC=LAYCON(K)
        IF(LC.EQ.4.OR.LC.EQ.5)THEN !
C2---------TWO STORAGES FOR UPSTREAM WEIGHTING FORMULATION
        IF(IALTSTO.EQ.1) THEN  ! COMPRESSIBLE STORAGE KICKS ON ONLY NEAR TOP OF CELL
          IF(SATN .LT. 0.99999) THEN ! BELOW (TOP - SMOOTH) SO NO COMPRESSIBLE STORAGE
            SATN = 0.0
          ELSE
            SATN = 1.0
          ENDIF
        ENDIF
          DTERMS = DTERMS -
     *      ((PRSITY(N)*VOLU-SC2(N)*ALENG) + ALENG*SC2(N)*SATN +
     *      SATN*SC1(N)*(HNEW(N)-HOLD(N))) / DELT
          RTERMS = RTERMS -
     *      ((PRSITY(N)*VOLU-SC2(N)*ALENG) + ALENG*SC2(N)*SATO )
     *       * CONCO(N,ICOMP) / DELT
        ELSEIF(LC.EQ.0.OR.LC.EQ.2) THEN
C3---------ONE CONFINED STORAGE CAPACITY FOR FLOW
            DTERMS = DTERMS - (PRSITY(N)*VOLU + SC1(N)*
     *        (HNEW(N)-HOLD(N))) / DELT
            RTERMS = RTERMS - (PRSITY(N)*VOLU) / DELT
     *       * CONCO(N,ICOMP)
        ELSEIF(LC.EQ.1) THEN ! IN BCF, IF LAYCON = 1, THE SPECIRFIC YIELD IS READ IN SC2 ?
C4---------ONE UNCONFINED STORAGE CAPACITY FOR FLOW
          DTERMS = DTERMS -
     *     ((PRSITY(N)*VOLU-SC2(N)*ALENG) + SC2(N)*ALENG*SATN) / DELT
          RTERMS = RTERMS -
     *     ((PRSITY(N)*VOLU-SC2(N)*ALENG) + SC2(N)*ALENG*SATO) / DELT
     *     * CONCO(N,ICOMP)
        ELSEIF(LC.EQ.3) THEN
C5---------TWO STORAGE CAPACITIES FOR FLOW
          TP=TOP(N)
          IF(HOLD(N).GT.TP.AND.HNEW(N).GT.TP) THEN
C6-----------CELL IS CONFINED
            DTERMS = DTERMS - (PRSITY(N)*VOLU + SC1(N)*
     *       (HNEW(N)-HOLD(N))) / DELT
            RTERMS = RTERMS - (PRSITY(N)*VOLU) / DELT
     *       * CONCO(N,ICOMP)
          ELSEIF (HOLD(N).LT.TP.AND.HNEW(N).LT.TP) THEN
C7-----------CELL IS UNCONFINED
            DTERMS = DTERMS -
     *       ((PRSITY(N)*VOLU-SC2(N)*ALENG) + SC2(N)*ALENG*SATN) / DELT
            RTERMS = RTERMS -
     *       ((PRSITY(N)*VOLU-SC2(N)*ALENG) + SC2(N)*ALENG*SATO) / DELT
     *     * CONCO(N,ICOMP)
          ELSE
C8-----------CELL HAS CONVERTED THIS TIME STEP
            IF(HNEW(N).GT.HOLD(N)) THEN
C8A-------------WLEs RAISED - BECOMES CONFINED
              H1 = HNEW(N)
              H2 = TP
            ELSE
C8B-------------WLEs LOWERED - BECOMES UNCONFINED
              H1 = TP
              H2 = HOLD(N)
            ENDIF
            DTERMS = DTERMS - ((PRSITY(N)*VOLU-SC2(N)*ALENG) +
     *       SC2(N)*ALENG*SATN + SC1(N)*(H1-H2)) / DELT
            RTERMS = RTERMS - ((PRSITY(N)*VOLU-SC2(N)*ALENG) +
     *       SC2(N)*ALENG*SATO) / DELT * CONCO(N,ICOMP)
          ENDIF
        ENDIF
      ENDIF
C
C9------RETURN
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE GWT2BCT1DCY(N,ICOMP,DTERMD,RTERMD)
C     ******************************************************************
C     FORMULATE DECAY TERM FOR EACH NODE AND COMPONENT
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,   ONLY:IOUT,NCOL,NROW,NLAY,IBOUND,TOP,NODLAY,ISYM,
     1 HNEW,AMAT,RHS,HOLD,ISSFLG,IA,JA,NJA,NODES,AREA,BOT,TOP,Sn,So
      USE CLN1MODULE, ONLY: ACLNNDS
      USE GWFBASMODULE,ONLY:DELT
      USE GWTBCTMODULE, ONLY: MCOMP,CONC,CONCO,ICT,IADSORB,PRSITY,
     * ADSORB,FLICH,IZOD,IFOD,ZODRW,FODRW,ZODRS,FODRS,IAW_ADSORB
      DOUBLE PRECISION ADSTERM,CW,FL,CWO,DT,RT,ALENG,VOLU,
     1  X,Y,EPSS,EPS,CT,CEPS,QA,QEPS,DQ,OMP,DTERMD,RTERMD
C     ------------------------------------------------------------------
C1------INITIALIZE TERMS
      DTERMD = 0.0
      RTERMD = 0.0
      IF(N.LE.NODES)THEN
        ALENG = TOP(N) - BOT(N)
      ELSE
        ALENG = ACLNNDS(N-NODES,4)
      ENDIF
      VOLU = AREA(N) * ALENG
C----------------------------------------------------------------------
      IF(ICT.EQ.0)THEN  !----------WATER PHASE CONCENTRATION FORMULATION
C2-------DECAY TERMS ON SOIL (NO ADSORPTION ON CLN)
        IF(N.LE.NODES)THEN
C3---------ZERO ORDER DECAY ON SOIL - APPLY NEWTON EXPANSION OF SUPPLY-DEMAND FUNCTION
          IF(IZOD.GE.2.AND.IADSORB.GT.0)THEN
            CT = - VOLU * ZODRS(N,ICOMP)
            EPS = 0.01
            CEPS = MAX(0.0,CONC(N,ICOMP))
            X = CEPS /EPS
            CALL SMOOTH(X,Y)
            QA =  CT * Y
C3A-----------CALCULATE DQ/DH
            EPSS = 0.001 * EPS
            CEPS = MAX(0.0,CONC(N,ICOMP)+EPSS)
            X = (CEPS)/EPS
            CALL SMOOTH(X,Y)
            QEPS = CT * Y
            DQ = (QEPS - QA) / EPSS
            DTERMD = DTERMD + DQ
            RTERMD = RTERMD - QA + DQ*CONC(N,ICOMP)
          ENDIF
C
C4---------FIRST ORDER DECAY ON SOIL
          IF(IFOD.GE.2.AND.IADSORB.GT.0)THEN
            CT = -ADSORB(N,ICOMP) * VOLU * FODRS(N,ICOMP)
            IF(IADSORB.EQ.1)THEN
C5--------------FOR LINEAR ADSORPTION
              DTERMD = DTERMD + CT
            ELSEIF(IADSORB.EQ.2) THEN ! FREUNDLICH
C6--------------FOR NON-LINEAR ADSORPTION FILL AS NEWTON
              ETA = FLICH(N,ICOMP)
              QA =  CT * CONC(N,ICOMP) ** ETA
C6A--------------CALCULATE DQ/DH
C              EPSS = 0.00001
C              CEPS = CONC(N,ICOMP)+EPSS
C              QEPS = CT * CEPS ** ETA
C              DQ = (QEPS - QA) / EPSS
CC              DQ = CT * ETA * CONC(N,ICOMP) ** (ETA-1.0)
CC              DTERMD = DTERMD + DQ
CC              RTERMD = RTERMD - QA + DQ*CONC(N,ICOMP)
              DTERMD = DTERMD + CT * CONC(N,ICOMP) **(ETA - 1.0) * ETA
            ELSEIF(IADSORB.EQ.3) THEN ! LANGMUIR
C6--------------FOR NON-LINEAR ADSORPTION FILL AS NEWTON
              ETA = FLICH(N,ICOMP)
              QA =  CT * CONC(N,ICOMP) / (1.0 + ETA*CONC(N,ICOMP))
C6A--------------CALCULATE DQ/DH
              EPSS = 0.00001
              CEPS = CONC(N,ICOMP)+EPSS
              QEPS = CT * CEPS / (1.0 + ETA * CEPS)
              DQ = (QEPS - QA) / EPSS
              DTERMD = DTERMD + DQ
              RTERMD = RTERMD - QA + DQ*CONC(N,ICOMP) 
            ENDIF
          ENDIF
C --------DECAY ON AIR-WATER INTERFACE ADSORBED SOLUTES          
          IF(IAW_ADSORB.NE.0) 
     *      CALL GWT2DCYAW(N,ICOMP,DTERMD,RTERMD,VOLU)            
        ENDIF
C-----------------------------------------------------------------------
C7---------ZERO ORDER DECAY IN WATER - APPLY NEWTON EXPANSION OF SUPPLY-DEMAND FUNCTION
        IF(IZOD.EQ.1.OR.IZOD.EQ.3)THEN
          CT = -Sn(N)* VOLU * ZODRW(N,ICOMP)
          EPS = 0.01
          CEPS = MAX(0.0,CONC(N,ICOMP))
          X = CEPS /EPS
          CALL SMOOTH(X,Y)
          QA =  CT * Y
C7A-----------CALCULATE DQ/DH
          EPSS = 0.001 * EPS
          CEPS = MAX(0.0,CONC(N,ICOMP)+EPSS)
          X = (CEPS)/EPS
          CALL SMOOTH(X,Y)
          QEPS = CT * Y
          DQ = (QEPS - QA) / EPSS
          DTERMD = DTERMD + DQ
          RTERMD = RTERMD - QA + DQ*CONC(N,ICOMP)
        ENDIF
C
C8---------FIRST ORDER DECAY IN WATER
        IF(IFOD.EQ.1.OR.IFOD.EQ.3)THEN
          CT = -Sn(N)* VOLU * FODRW(N,ICOMP)
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
      SUBROUTINE GWT2BCT1GEN(N,ICOMP)
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
     * ADSORB,FLICH,IZOD,IFOD,ZODRW,FODRW,IAW_ADSORB,
     1  ZODRS,FODRS,ICHAIN,NPARENT,JPARENT,STOTIO,ISPRCT,SPTLRCT
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
          STOIT = SPTLRCT(N,ICOMP,NPAREN)
        ENDIF
C----------------------------------------------------------------------
      IF(ICT.EQ.0)THEN  !----------WATER PHASE CONCENTRATION FORMULATION
C2-------GENERATION TERMS ON SOIL (NO ADSORPTION ON CLN)
        IF(N.LE.NODES)THEN
C3---------ZERO ORDER DECAY ON SOIL - APPLY NEWTON EXPANSION OF SUPPLY-DEMAND FUNCTION
          IF(IZOD.GE.2.AND.IADSORB.GT.0)THEN
            CT = VOLU * ZODRS(N,IPAREN) * STOIT
            EPS = 0.01
            CEPS = MAX(0.0,CONC(N,IPAREN))
            X = CEPS /EPS
            CALL SMOOTH(X,Y)
            QA =  CT * Y
C3A-----------CALCULATE DQ/DH
            EPSS = 0.001 * EPS
            CEPS = MAX(0.0,CONC(N,IPAREN)+EPSS)
            X = (CEPS)/EPS
            CALL SMOOTH(X,Y)
            QEPS = CT * Y
            DQ = (QEPS - QA) / EPSS
            DTERMD = DTERMD + DQ
            RTERMD = RTERMD - QA + DQ*CONC(N,IPAREN)
          ENDIF
C
C4---------FIRST ORDER DECAY ON SOIL
          IF(IFOD.GE.2.AND.IADSORB.GT.0)THEN
            CT = ADSORB(N,IPAREN) * VOLU * FODRS(N,IPAREN) * STOIT
            IF(IADSORB.EQ.1)THEN
C5--------------FOR LINEAR ADSORPTION
              DTERMD = DTERMD + CT
            ELSEIF(IADSORB.EQ.2) THEN
C6--------------FOR NON-LINEAR ADSORPTION FILL AS NEWTON
              ETA = FLICH(N,IPAREN)
              QA =  CT * CONC(N,IPAREN) ** ETA
C6A--------------CALCULATE DQ/DH
C              EPSS = 0.00001
C              CEPS = CONC(N,IPAREN)+EPSS
C              QEPS = CT * CEPS ** ETA
C              DQ = (QEPS - QA) / EPSS
CC              DQ = CT * ETA * CONC(N,IPAREN) ** (ETA-1.0)
CC              DTERMD = DTERMD + DQ
CC              RTERMD = RTERMD - QA + DQ*CONC(N,IPAREN)
              DTERMD = DTERMD + CT * CONC(N,IPAREN) **(ETA - 1.0)
            ELSEIF(IADSORB.EQ.3) THEN ! LANGMUIR
C6--------------FOR NON-LINEAR ADSORPTION FILL AS NEWTON
              ETA = FLICH(N,ICOMP)
              QA =  CT * CONC(N,ICOMP) / (1.0 + ETA*CONC(N,ICOMP))
C6A--------------CALCULATE DQ/DH
              EPSS = 0.00001
              CEPS = CONC(N,ICOMP)+EPSS
              QEPS = CT * CEPS / (1.0 + ETA * CEPS)
              DQ = (QEPS - QA) / EPSS
              DTERMD = DTERMD + DQ
              RTERMD = RTERMD - QA + DQ*CONC(N,ICOMP) 
            ENDIF
          ENDIF
C --------GENERATION ON AIR-WATER INTERFACE ADSORBED SOLUTES  
            IF(IAW_ADSORB.NE.0)          
     *      CALL GWT2GENAW(N,ICOMP,DTERMD,RTERMD,VOLU,NPAREN)    
        ENDIF
C-----------------------------------------------------------------------
C7---------ZERO ORDER DECAY IN WATER - APPLY NEWTON EXPANSION OF SUPPLY-DEMAND FUNCTION
        IF(IZOD.EQ.1.OR.IZOD.EQ.3)THEN
          CT = Sn(N)* VOLU * ZODRW(N,IPAREN)*STOIT
          EPS = 0.01
          CEPS = MAX(0.0,CONC(N,IPAREN))
          X = CEPS /EPS
          CALL SMOOTH(X,Y)
          QA =  CT * Y
C7A-----------CALCULATE DQ/DH
          EPSS = 0.001 * EPS
          CEPS = MAX(0.0,CONC(N,IPAREN)+EPSS)
          X = (CEPS)/EPS
          CALL SMOOTH(X,Y)
          QEPS = CT * Y
          DQ = (QEPS - QA) / EPSS
          DTERMD = DTERMD + DQ
          RTERMD = RTERMD - QA + DQ*CONC(N,IPAREN)
        ENDIF
C
C8---------FIRST ORDER DECAY IN WATER
        IF(IFOD.EQ.1.OR.IFOD.EQ.3)THEN
          CT = Sn(N)* VOLU * FODRW(N,IPAREN) * STOIT
          DTERMD = DTERMD + CT
        ENDIF
      ELSE       !-----------------------TOTAL CONCENTRATION FORMULATION
C
C-------NET DECAY TERM FOR TOTAL CONCENTRATION FORMULATION
CSP TO DO
      ENDIF
      RHS(N) = RHS(N) + RTERMD - DTERMD * CONC(N,IPAREN)
      ENDDO   !---------------------------END DO LOOP OF PARENTS
C
C9------RETURN
      RETURN
      END
C----------------------------------------------------------------------
      SUBROUTINE GWT2BCT1STOTVM(N,ICOMP,DTERMS,RTERMS,ISS)

C     ******************************************************************
C     FORMULATE STORAGE TERM FOR EACH TRANSIENT PROPERTY NODE AS A CORRECTION
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
     1  ZODRS,FODRS
      DOUBLE PRECISION VODT,ADSTERM,CW,FL,CWO,DT,RT,ALENG,VOLU,
     1  X,Y,EPSS,EPS,CT,CEPS,QA,QEPS,DQ,OMP,DTERMS,RTERMS
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
C4-------STORAGE TERM IN WATER
        CALL GWT2BCT1STOWTVM (N,ICOMP,DTERMS,RTERMS,VODT,VOLU,ALENG,ISS)
      ENDIF
C
C9------RETURN
      RETURN
      END
      SUBROUTINE GWT2BCT1STOWTVM(N,ICOMP,DTERMS,RTERMS,VODT,VOLU,ALENG,
     1  ISS)
C     ******************************************************************
C     FORMULATE STORAGE TERM IN WATER FOR EACH TRANSIENT PROPERTY NODE AS A CORRECTION
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,   ONLY:IOUT,NCOL,NROW,NLAY,IBOUND,TOP,NODLAY,ISYM,
     1 HNEW,AMAT,RHS,HOLD,ISSFLG,IA,JA,NJA,NODES,AREA,BOT,TOP,Sn,So,
     2 LAYNOD
      USE CLN1MODULE, ONLY: ACLNNDS
      USE GWFBASMODULE,ONLY:DELT
      USE GWFBCFMODULE,ONLY:LAYCON,SC1,SC2
      USE GWTBCTMODULE, ONLY: MCOMP,CONC,CONCO,ICT,IADSORB,PRSITY,
     * ADSORB,FLICH,IZOD,IFOD,ZODRW,FODRW,
     1  ZODRS,FODRS
      USE TVMU2MODULE
      DOUBLE PRECISION VODT,ADSTERM,CW,FL,CWO,DT,RT,ALENG,VOLU,
     1  X,Y,EPSS,EPS,CT,CEPS,QA,QEPS,DQ,OMP,DTERMS,RTERMS,
     1  H1,H2,TP,SATN,SATO,VODTO,VODTN
C-------------------------------------------------------------------------------------
C
      SATN = Sn(N)
      SATO = So(N)
      IF(ISS.EQ.1.OR.N.GT.NODES) THEN
        IF(ISS.EQ.1) SATO = SATN
C1-------USE TRADITIONAL STORAGE TERM FOR CLN DOMAIN OR SS FLOW
        IF(N.LE.NODES) THEN
          VODTO = VODT * POROLD(N)
          VODTN = VODT * PRSITY(N)
          RTERMS = RTERMS + (VODTN*SATN - VODTO*SATO) * CONCO(N,ICOMP)
        ENDIF
      ELSE
        K = LAYNOD(N)
        LC=LAYCON(K)
        IF(LC.EQ.4.OR.LC.EQ.5)THEN !
C2---------TWO STORAGES FOR UPSTREAM WEIGHTING FORMULATION
          RTERMS = RTERMS +
     *      ( ((PRSITY(N)*VOLU-SC2(N)*ALENG) + ALENG*SC2(N)*SATO +
     *        SATN*SC1(N)*HOLD(N)) -
     *        ((POROLD(N)*VOLU-SYOLD(N)*ALENG) + ALENG*SYOLD(N)*SATO +
     *        SATN*SSOLD(N)*HOLD(N))
     *      ) * CONCO(N,ICOMP) / DELT
        ELSEIF(LC.EQ.0.OR.LC.EQ.2) THEN
C3---------ONE CONFINED STORAGE CAPACITY FOR FLOW
            RTERMS = RTERMS +
     *        ( (PRSITY(N)*VOLU + SC1(N)*HOLD(N)) -
     *          (POROLD(N)*VOLU + SSOLD(N)*HOLD(N))
     *        ) * CONCO(N,ICOMP) / DELT
        ELSEIF(LC.EQ.1) THEN ! IN BCF, IF LAYCON = 1, THE SPECIRFIC YIELD IS READ IN SC2 ?
C4---------ONE UNCONFINED STORAGE CAPACITY FOR FLOW
          RTERMS = RTERMS +
     *      ( ((PRSITY(N)*VOLU-SC2(N)*ALENG) + SC2(N)*ALENG*SATO) -
     *        ((POROLD(N)*VOLU-SYOLD(N)*ALENG) + SYOLD(N)*ALENG*SATO)
     *      ) * CONCO(N,ICOMP) / DELT
        ELSEIF(LC.EQ.3) THEN
C5---------TWO STORAGE CAPACITIES FOR FLOW
          TP=TOP(N)
          IF(HOLD(N).GT.TP.AND.HNEW(N).GT.TP) THEN
C6-----------CELL IS CONFINED
            RTERMS = RTERMS +
     *        ( (PRSITY(N)*VOLU + SC1(N)*HOLD(N)) -
     *          (POROLD(N)*VOLU + SSOLD(N)*HOLD(N))
     *        ) * CONCO(N,ICOMP) / DELT
          ELSEIF (HOLD(N).LT.TP.AND.HNEW(N).LT.TP) THEN
C7-----------CELL IS UNCONFINED
            RTERMS = RTERMS +
     *        ( ((PRSITY(N)*VOLU-SC2(N)*ALENG) + SC2(N)*ALENG*SATO) -
     *          ((POROLD(N)*VOLU-SYOLD(N)*ALENG) + SYOLD(N)*ALENG*SATO)
     *        ) * CONCO(N,ICOMP) / DELT
          ELSE
C8-----------CELL HAS CONVERTED THIS TIME STEP
            IF(HNEW(N).GT.HOLD(N)) THEN
C8A-------------WLEs RAISED - BECOMES CONFINED
              H1 = HNEW(N)
              H2 = TP
            ELSE
C8B-------------WLEs LOWERED - BECOMES UNCONFINED
              H1 = TP
              H2 = HOLD(N)
            ENDIF
            RTERMS = RTERMS +
     *        ( ((PRSITY(N)*VOLU-SC2(N)*ALENG) + SC2(N)*ALENG*SATO +
     *          SC1(N)*H2) -
     *          ((POROLD(N)*VOLU-SYOLD(N)*ALENG) + SYOLD(N)*ALENG*SATO +
     *          SSOLD(N)*H2)
     *        ) * CONCO(N,ICOMP)/ DELT
          ENDIF
        ENDIF
      ENDIF
C
C9------RETURN
      RETURN
      END
C--------------------------------------------------------------------------
      SUBROUTINE GWT2SOLVE(ITITER,ICNVG,ICOMP,KSTP,KPER,ITP,IN_ITER)
C     ******************************************************************
C     SOLVE THE GROUNDWATER TRANSPORT EQUATION FOR THIS COMPONENT
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,  ONLY: IOUT,NEQS,IBOUND,AMAT,RHS,IA,NJA
      USE SMSMODULE,ONLY:HTEMP,LINMETH
      USE XMDMODULE,ONLY:IACL
      USE PCGUMODULE, ONLY: ILINMETH
      USE GWTBCTMODULE, ONLY: MCOMP,CONC,ICBUND,CICLOSE,
     1 CINACT,CNCG,LRCC,MSUMT,FCNCG,FLRCC,CONCO
      DOUBLE PRECISION BIG,ADIAG,ROD
C     ------------------------------------------------------------------
C1-----TAKE CARE OF LOOSE ENDS FOR ALL NODES BEFORE CALL TO SOLVER
          BIG = 1.0E20
          DO N=1,NEQS
C2---------TAKE CARE OF ZERO ROW DIAGONAL
            ADIAG = ABS(AMAT(IA(N)))
            IF(ADIAG.LT.1.0E-06)THEN
              AMAT(IA(N)) = 1.0E08
              RHS(N) = CONC(N,ICOMP)*1.0E08       !   + RHS(N)
            ENDIF
C2B--------TAKE CARE OF ZERO OFFDIAGONALS (NO INFLOW) BUT LARGE RHS/DIAGONAL
c            IOFFDIAG = 0
c            DO JJ = IA(N)+1,IA(N+1)-1
c              IF (ABS(AMAT(JJ)).GT.1.0E-6) THEN
c                IOFFDIAG = 1
c                GO TO 11
c              ENDIF
c            ENDDO
c11          CONTINUE
c            ROD = ABS(RHS(N) / AMAT(IA(N)))
c            IF(IOFFDIAG.EQ.0. AND. ROD. GT.1.0E6) THEN
c              AMAT(IA(N)) = 1.0E00
c              RHS(N) = CONC(N,ICOMP)
c            ENDIF
          ENDDO
C
C3 --------STORE CONC IN TEMP LOCATION
          DO N=1,NEQS
            HTEMP(N) = CONC(N,ICOMP)
          ENDDO
C4 ---------ENSURE SOLVER IS UNSYMMETRIC
          IF(LINMETH.EQ.1)THEN    ! XMD SOLVER
            IACLTMP = -1
            IF(IACL.EQ.0)THEN
              IACLTMP = IACL
              IACL = 1
            ENDIF
          ELSE                    ! PCGU SOLVER
            ILINMETHTMP = -1
            IF(ILINMETH.EQ.1)THEN
              ILINMETHTMP = ILINMETH
              ILINMETH = 2
            ENDIF
          ENDIF
C5 ---------EVALUATE CONDITIONS WHEN ILU NEED NOT BE PERFORMED
          ILUFLAG = 1
cccc          IF(ITITER.GT.1) ILUFLAG = 0
cccc          IF(IADMATF.EQ.1.AND.DABS(DELT - DELTO).LT.1.0E-15) ILUFLAG = 0
C6---------CALL LINEAR SOLVER ROUTINES
          CALL SOLVERS(IOUT,ITITER,ICNVG,KSTP,KPER,AMAT,CONC(1,ICOMP),
     *      RHS,ICBUND,CICLOSE,CINACT,ITP,NEQS,NJA,ILUFLAG,IN_ITER,
     *      CONCO(1,ICOMP) )
C
C7--------RESET SOLVER TO SYMMETRIC IF CHANGED
         IF(LINMETH.EQ.1)THEN    ! XMD SOLVER
            IF(IACLTMP.EQ.0)THEN
              IACL = 0
            ENDIF
          ELSE                    ! PCGU SOLVER
            IF(ILINMETHTMP.EQ.1)THEN
              ILINMETH = 1
            ENDIF
          ENDIF
C
C9------RETURN
      RETURN
      END
C----------------------------------------------------------------------
      SUBROUTINE CONCCONV(ICOMP,ITITER,ICNVG)
C     ******************************************************************
C     COMPUTE MAX CONCENTRATION CHANGES AND CHECK CONVERGENCE
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,  ONLY:IUNSTR,NODES,NEQS,NLAY,NROW,NCOL,INCLN,
     1  IBOUND
      USE CLN1MODULE, ONLY: NCLNNDS
      USE SMSMODULE,ONLY:HTEMP
      USE GWTBCTMODULE, ONLY: MCOMP,CONC,ICBUND,CICLOSE,
     1 CINACT,CNCG,LRCC,MSUMT,FCNCG,FLRCC
      DOUBLE PRECISION BIGCC,ABIGCC,CDIF,ACDIF,BIGCL,ABIGCL
C     ------------------------------------------------------------------
C
C1---------STORE LARGEST CONC CHANGE FOR ITERATION
          NB=1
          ICNVG=0
          BIGCC=0.0
          ABIGCC=0.0
          DO N=1,NODES
            IF(ICBUND(N).EQ.0) CYCLE
            CDIF=CONC(N,ICOMP)-HTEMP(N)
            ACDIF=ABS(CDIF)
            IF(ACDIF.GE.ABIGCC)THEN
              BIGCC= CDIF
              ABIGCC= ACDIF
              NB = N
            ENDIF
          ENDDO
          IF(ABIGCC.LE.CICLOSE) ICNVG=1
C
C2---------STORE MAXIMUM CHANGE VALUE AND LOCATION
          CNCG(ITITER) = BIGCC
C
          IF(IUNSTR.EQ.0)THEN !GET LAYER, ROW AND COLUMN FOR STRUCTURED GRID
            KLAYER = (NB-1) / (NCOL*NROW) + 1
            IJ = NB - (KLAYER-1)*NCOL*NROW
            IROW = (IJ-1)/NCOL + 1
            JCOLMN = IJ - (IROW-1)*NCOL
            LRCC(1,ITITER) = KLAYER
            LRCC(2,ITITER) = IROW
            LRCC(3,ITITER) = JCOLMN
          ELSE
            LRCC(1,ITITER) = NB
          ENDIF
C
C3---------CHECK OUTER ITERATION CONVERGENCE FOR CONDUIT-NODES
          IF(INCLN.EQ.0) GO TO 204
          NB=1
          ICNVGL=0
          BIGCCL=0.0
          ABIGCCL=0.0
          DO N=NODES+1,NODES+NCLNNDS
            IF(IBOUND(N).EQ.0) CYCLE
            CDIF=CONC(N,ICOMP)-HTEMP(N)
            ACDIF=ABS(CDIF)
            IF(ACDIF.GE.ABIGCCL)THEN
              BIGCCL= CDIF
              ABIGCCL= ACDIF
              NB = N
            ENDIF
          ENDDO
C
          IF(ABIGCCL.LE.CICLOSE) ICNVGL=1
C
C4---------STORE MAXIMUM CHANGE VALUE AND LOCATION
          FCNCG(ITITER) = BIGCCL
          FLRCC(ITITER) = NB - NODES
C5---------NOT CONVERGED, IF EITHER IS NOT CONVERGED
          IF(ICNVG.EQ.0.OR.ICNVGL.EQ.0) ICNVG = 0
204       CONTINUE
C
C6------RETURN
      RETURN
      END

C-----------------------------------------------------------------------
      SUBROUTINE GWT2BCT1OC(KSTP,KPER,ISA)
C     ******************************************************************
C     SET PRINT FLAGS FOR TRANSPORT
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,      ONLY:ITMUNI,IOUT,IUNSTR,ISSFLG,PERLEN
      USE GWFBASMODULE,ONLY:DELT,PERTIM,TOTIM,ISPCFL,IBUDFL,IATS,TMINAT,
     1  NPTIMES,NPSTPS,ITIMOTC,TIMOTC,ISPCFLAT
      USE GWTBCTMODULE, ONLY: MCOMP,VBVLT,VBNMT,MSUMT
C     ------------------------------------------------------------------
C
C1------DETERMINE IF TRANSPORT PRINTING IS NEEDED AT THIS TIME
      IF(ISSFLG(KPER).EQ.1 .AND. IATS.NE.0)THEN
C      IF(IATS.NE.0)THEN
        IATSPR = 0
C1B-----PRINT EVERY NPSTPS TIMES
        IF (MOD(KSTP,NPSTPS).EQ.0) IATSPR = 1
C1C-----PRINT AT END OF STRESS PERIOD
        IF(ABS(PERTIM-PERLEN(KPER)).LT.TMINAT) IATSPR = 1
C1D-----PRINT WHEN OUTPUT TIME IS REACHED
        IF(NPTIMES.GT.0)THEN
          IF(TOTIM.GE.(TIMOTC(ITIMOTC) - TMINAT))THEN
            IATSPR = 1
            ITIMOTC = ITIMOTC + 1
          ENDIF
        ENDIF
C2------SET FLAGS IF PRINTING
        ISPCFL = 0
        IF(IATSPR.EQ.1) THEN
          ISPCFL = ISPCFLAT
        ENDIF
      ENDIF
C
C3------RETURN
      RETURN
      END
C---------------------------------------------------------------------------------
      SUBROUTINE GWT2BCT1OT(KSTP,KPER,ISA)
C     ******************************************************************
C     OUTPUT CONCENTRATION FOR EACH SPECIES AND MASS BALANCE SUMMARY
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,      ONLY:ITMUNI,IOUT,IUNSTR,ISSFLG,PERLEN,INCLN,IDPT
      USE GWFBASMODULE,ONLY:DELT,PERTIM,TOTIM,ISPCFL,IBUDFL,IATS,TMINAT,
     1  NPTIMES,NPSTPS,ITIMOT,TIMOT,ISPCFLAT
      USE GWTBCTMODULE, ONLY: MCOMP,VBVLT,VBNMT,MSUMT,MCOMPT,NTCOMP
C     ------------------------------------------------------------------
C
C1------CLEAR PRINTOUT FLAG (IPFLG)
      IPFLG=0
      WRITE(IOUT,1)
1     FORMAT(/1X,'TRANSPORT SOLUTION COMPLETE FOR ALL SPECIES, SORTING',
     1  ' OUTPUT FOR CONCENTRATIONS AND MASS BALANCE'/1X,100('-'))
      DO ICOMP = 1,NTCOMP
        WRITE(IOUT,2)ICOMP
2       FORMAT(5X,'TRANSPORT OUTPUT FOR COMPONENT SPECIES',
     *    1X,'NUMBER',I5/5X,60('-'))
C
C2------IF CONCENTRATION PRINT FLAG (ISPCFL) IS SET WRITE
C2------CONCENTRATION, AND ICBUND IN ACCORDANCE WITH FLAGS IN IOFLG.
        IF(ISPCFL.EQ.0) GO TO 100
C
        IF(IUNSTR.EQ.0)THEN ! WRITE M2K5 STYLE FOR STRUCTURED GRID
          CALL SGWT2BCT1C(KSTP,KPER,IPFLG,ISA,ICOMP)
csp          CALL SGWT2BCT1IB(KSTP,KPER)
        ELSE
          CALL SGWT2BCT1CU(KSTP,KPER,IPFLG,ISA,ICOMP)
csp          CALL SGWT2BCT1IBU(KSTP,KPER)
        ENDIF
        IPFLG = 1
C
C-------------------------------------------------------------
C2B-----FOR CONDUIT NODES
C-------------------------------------------------------------
        IF(INCLN.GT.0)THEN
          CALL SCLN1C(KSTP,KPER,IPFLG,ISA,ICOMP)
        ENDIF
C-------------------------------------------------------------
C
  100   CONTINUE

C3------PRINT TOTAL BUDGET IF REQUESTED
      IF(IBUDFL.EQ.0) GO TO 120
      ITF = 1   ! TRANSPORT EQUATION WAS SOLVED      
        CALL SGWF2BAS7V(MSUMT,VBNMT(1,ICOMP),VBVLT(1,1,ICOMP),
     *  KSTP,KPER,IOUT,ITF)
      IPFLG=1
C3a------PRINT TOTAL MASS IN SYSTEM TO OUTPUT FILE
      ISTARTMASS = 0 
      CALL INITMASS (KPER,KSTP,ISTARTMASS)
      IF(IDPT.NE.0) CALL INITMASSIM 
      IF(INCLN.NE.0) CALL INITMASSCLN         
C
C4------END PRINTOUT WITH TIME SUMMARY AND FORM FEED IF ANY PRINTOUT
C4------WILL BE PRODUCED.
  120   IF(IPFLG.EQ.0) GO TO 99
        CALL SGWF2BAS7T(KSTP,KPER,DELT,PERTIM,TOTIM,ITMUNI,IOUT)
        WRITE(IOUT,101)
  101   FORMAT('1')
C
      ENDDO
C
C5------RETURN
99    CONTINUE
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE SGWT2BCT1C(KSTP,KPER,IPFLG,ISA,ICOMP)
C     ******************************************************************
C     PRINT AND RECORD CONCS
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,      ONLY:NCOL,NROW,NLAY,IXSEC,HNEW,NODLAY,
     1                      IBOUND,IOUT,IDPOUT
      USE GWFBASMODULE,ONLY:PERTIM,TOTIM,ISPCFM,ISPCUN,LBHDSV,
     2                      CSPCFM,IOFLG
      USE GWTBCTMODULE, ONLY: ICBUND,CONC,IMULTI,CROOTNAME  
C
      REAL,          SAVE,    DIMENSION(:,:,:),    ALLOCATABLE ::BUFF
      REAL*8,        SAVE,    DIMENSION(:,:,:),    ALLOCATABLE ::BUFF8
      CHARACTER*400 CFILENAME      
      CHARACTER*16 TEXT,ISTATUS,FORMC 
      INCLUDE 'openspec.inc'      
CSP      DATA TEXT /'            CONC'/
C     ------------------------------------------------------------------
      ALLOCATE(BUFF(NCOL,NROW,NLAY))
      ALLOCATE(BUFF8(NCOL,NROW,NLAY))
      CALL GET_TEXTC (ICOMP,TEXT)
C
C1------FOR EACH LAYER MOVE HNEW TO BUFF IF PRINT OR SAVE IS REQUESTED.
      DO 59 K=1,NLAY
C
C2------IS CONC NEEDED FOR THIS LAYER?
      KL=K
      IF(IXSEC.NE.0) KL=1
      IF(IOFLG(KL,5).EQ.0.AND.IOFLG(KL,6).EQ.0) GO TO 59
C
C3------MOVE CONC TO BUFF FOR THE LAYER.
      DO 58 I=1,NROW
      DO 58 J=1,NCOL
      N = (K-1)*NROW*NCOL + (I-1)*NCOL + J
      BUFF8(J,I,K)=CONC(N,ICOMP)
      BUFF(J,I,K)=BUFF8(J,I,K)
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
             IF(ISPCFM.LT.0) CALL ULAPRS(BUFF(1,1,1),TEXT,KSTP,KPER,
     1                 NCOL,NLAY,-1,-ISPCFM,IOUT)
             IF(ISPCFM.GE.0) CALL ULAPRW(BUFF(1,1,1),TEXT,KSTP,KPER,
     1                 NCOL,NLAY,-1,ISPCFM,IOUT)
             IPFLG=1
           END IF
         END IF
      END IF
C
C5------FOR EACH LAYER: DETERMINE IF CONC SHOULD BE SAVED ON DISK.
C5------IF SO THEN CALL ULASAV OR ULASV2 TO SAVE CONC.
      IFIRST=1
      IF(ISPCUN.LE.0) GO TO 80
C--------------------------------------------------------------------------------------------------------------------
C5A ----FOR MULTIFILE OUTPUT OPTION, ADD COMPONENT SPECIES TO ROOT NAME AND OPEN FILE ON ISPCFM AND GO TO END OF FILE
      IF(IMULTI.NE.0)THEN
        CFILENAME = trim(CROOTNAME) // TRIM(TEXT) //'.CON'
        ISTATUS = 'OLD'
        FORMC = 'FORMATTED'
        IF(CSPCFM.EQ.' ') FORMC = FORM            ! - use FORM from openspec.inc rather than hardwiring BINARY or UNFORMATTED
CC        IF(KSTP.EQ.1.AND.KPER.EQ.1) THEN
C---------MAKE STATUS OF FILE AS NEW WHEN FIRST ENTERING FILE
CC          ISTATUS = 'REPLACE'
CC        ENDIF

!        OPEN(ISPCUN, FILE = CFILENAME, FORM = 'BINARY',   
        OPEN(ISPCUN, FILE = CFILENAME, FORM = FORMC,
     1  STATUS = ISTATUS, POSITION = 'APPEND')
      ENDIF
C--------------------------------------------------------------------------------------------------------------------      
      IF(IXSEC.EQ.0) THEN
        DO 79 K=1,NLAY
        NSTRT = NODLAY(K-1)+1
        KK=K
        IF(IOFLG(K,6).EQ.0) GO TO 79
        IF(IFIRST.EQ.1) WRITE(IOUT,74) ISPCUN,KSTP,KPER
   74   FORMAT(1X,/1X,'CONC WILL BE SAVED ON UNIT ',I4,
     1      ' AT END OF TIME STEP ',I8,', STRESS PERIOD ',I8)
        IFIRST=0
        IF(IDPOUT.EQ.1) THEN
           WRITE(IOUT,*)'  CONC AND TIME VARIABLES ARE SAVED AS REAL*8' 
           CALL ULASAV8(BUFF8(1,1,K),TEXT,KSTP,KPER,PERTIM,TOTIM,NCOL,
     1                NROW,KK,ISPCUN)        
        ELSEIF(CSPCFM.EQ.' ') THEN
           CALL ULASAV(BUFF(1,1,K),TEXT,KSTP,KPER,PERTIM,TOTIM,NCOL,
     1                NROW,KK,ISPCUN)
        ELSE
           CALL ULASV2(BUFF(1,1,K),TEXT,KSTP,KPER,PERTIM,TOTIM,NCOL,
     1                NROW,KK,ISPCUN,CSPCFM,LBHDSV,ICBUND(NSTRT))
        END IF
   79   CONTINUE
C
C5A-----SAVE CONC FOR CROSS SECTION.
      ELSE
        IF(IOFLG(1,6).NE.0) THEN
          WRITE(IOUT,74) ISPCUN,KSTP,KPER
          IF(IDPOUT.EQ.1)THEN
           WRITE(IOUT,*)'  CONC AND TIME VARIABLES ARE SAVED AS REAL*8' 
             CALL ULASAV8(BUFF8,TEXT,KSTP,KPER,PERTIM,TOTIM,NCOL,
     1                NLAY,-1,ISPCUN)          
          ELSEIF(CSPCFM.EQ.' ') THEN
             CALL ULASAV(BUFF(1,1,1),TEXT,KSTP,KPER,PERTIM,TOTIM,NCOL,
     1                NLAY,-1,ISPCUN)
          ELSE
             CALL ULASV2(BUFF(1,1,1),TEXT,KSTP,KPER,PERTIM,TOTIM,NCOL,
     1                  NLAY,-1,ISPCUN,CSPCFM,LBHDSV,ICBUND)
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
      SUBROUTINE GET_TEXTC(ICOMP,TEXT)
C     ******************************************************************
C     SELECT TEXT FOR OUTPUT FILE DEPENDING ON SPECIES NUMBER
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------

      USE GWTBCTMODULE, ONLY: MCOMP,IHEAT,MCOMPT,IMULTI
      CHARACTER*16 TEXT
C     ------------------------------------------------------------------
C
      IF(IHEAT.EQ.1.AND.ICOMP.EQ.MCOMPT) THEN
        TEXT = ' TMPR '
        RETURN
      ENDIF
      IF(MCOMP.EQ.1) THEN
        TEXT = ' CONC'
      ELSE
        IF(ICOMP.EQ.1) THEN
          TEXT = ' CONC01'
        ELSEIF(ICOMP.EQ.2) THEN
          TEXT = ' CONC02'
        ELSEIF(ICOMP.EQ.3) THEN
          TEXT = ' CONC03'
        ELSEIF(ICOMP.EQ.4) THEN
          TEXT = ' CONC04'
        ELSEIF(ICOMP.EQ.5) THEN
          TEXT = ' CONC05'
        ELSEIF(ICOMP.EQ.6) THEN
          TEXT = ' CONC06'
        ELSEIF(ICOMP.EQ.7) THEN
          TEXT = ' CONC07'
        ELSEIF(ICOMP.EQ.8) THEN
          TEXT = ' CONC08'
        ELSEIF(ICOMP.EQ.9) THEN
          TEXT = ' CONC09'
        ELSEIF(ICOMP.EQ.10) THEN
          TEXT = ' CONC10'
        ELSEIF(ICOMP.EQ.11) THEN
          TEXT = ' CONC11'
        ELSEIF(ICOMP.EQ.12) THEN
          TEXT = ' CONC12'
        ELSEIF(ICOMP.EQ.13) THEN
          TEXT = ' CONC13'
        ELSEIF(ICOMP.EQ.14) THEN
          TEXT = ' CONC14'
        ELSEIF(ICOMP.EQ.15) THEN
          TEXT = ' CONC15'
        ELSEIF(ICOMP.EQ.16) THEN
          TEXT = ' CONC16'
        ELSEIF(ICOMP.EQ.17) THEN
          TEXT = ' CONC17'
        ELSEIF(ICOMP.EQ.18) THEN
          TEXT = ' CONC18'
        ELSEIF(ICOMP.EQ.19) THEN
          TEXT = ' CONC19'
        ELSEIF(ICOMP.EQ.20) THEN
          TEXT = ' CONC20'
        ELSEIF(ICOMP.EQ.21) THEN
          TEXT = ' CONC21'
        ELSEIF(ICOMP.EQ.22) THEN
          TEXT = ' CONC22'
        ELSEIF(ICOMP.EQ.23) THEN
          TEXT = ' CONC23'
        ELSEIF(ICOMP.EQ.24) THEN
          TEXT = ' CONC24'
        ELSEIF(ICOMP.EQ.25) THEN
          TEXT = ' CONC25'
        ELSEIF(ICOMP.EQ.26) THEN
          TEXT = ' CONC26'
        ELSEIF(ICOMP.EQ.27) THEN
          TEXT = ' CONC27'
        ELSEIF(ICOMP.EQ.28) THEN
          TEXT = ' CONC28'
        ELSEIF(ICOMP.EQ.29) THEN
          TEXT = ' CONC29'
        ELSEIF(ICOMP.EQ.30) THEN
          TEXT = ' CONC30'
        ELSEIF(ICOMP.EQ.31) THEN
          TEXT = ' CONC31'
        ELSEIF(ICOMP.EQ.32) THEN
          TEXT = ' CONC32'
        ELSEIF(ICOMP.EQ.33) THEN
          TEXT = ' CONC33'
        ELSEIF(ICOMP.EQ.34) THEN
          TEXT = ' CONC34'
        ELSEIF(ICOMP.EQ.35) THEN
          TEXT = ' CONC35'
        ELSEIF(ICOMP.EQ.36) THEN
          TEXT = ' CONC36'
        ELSEIF(ICOMP.EQ.37) THEN
          TEXT = ' CONC37'
        ELSEIF(ICOMP.EQ.38) THEN
          TEXT = ' CONC38'
        ELSEIF(ICOMP.EQ.39) THEN
          TEXT = ' CONC39'
        ELSEIF(ICOMP.EQ.40) THEN
          TEXT = ' CONC40'
        ELSEIF(ICOMP.EQ.41) THEN
          TEXT = ' CONC41'
        ELSEIF(ICOMP.EQ.42) THEN
          TEXT = ' CONC42'
        ELSEIF(ICOMP.EQ.43) THEN
          TEXT = ' CONC43'
        ELSEIF(ICOMP.EQ.44) THEN
          TEXT = ' CONC44'
        ELSEIF(ICOMP.EQ.45) THEN
          TEXT = ' CONC45'
        ELSEIF(ICOMP.EQ.46) THEN
          TEXT = ' CONC46'
        ELSEIF(ICOMP.EQ.47) THEN
          TEXT = ' CONC47'
        ELSEIF(ICOMP.EQ.48) THEN
          TEXT = ' CONC48'
        ELSEIF(ICOMP.EQ.49) THEN
          TEXT = ' CONC49'
        ELSEIF(ICOMP.EQ.50) THEN
          TEXT = ' CONC50'
        ELSEIF(ICOMP.EQ.51) THEN
          TEXT = ' CONC51'
        ELSEIF(ICOMP.EQ.52) THEN
          TEXT = ' CONC52'
        ELSEIF(ICOMP.EQ.53) THEN
          TEXT = ' CONC53'
        ELSEIF(ICOMP.EQ.54) THEN
          TEXT = ' CONC54'
        ELSEIF(ICOMP.EQ.55) THEN
          TEXT = ' CONC55'
        ELSEIF(ICOMP.EQ.56) THEN
          TEXT = ' CONC56'
        ELSEIF(ICOMP.EQ.57) THEN
          TEXT = ' CONC57'
        ELSEIF(ICOMP.EQ.58) THEN
          TEXT = ' CONC58'
        ELSEIF(ICOMP.EQ.59) THEN
          TEXT = ' CONC59'
        ELSEIF(ICOMP.EQ.60) THEN
          TEXT = ' CONC60'
        ELSEIF(ICOMP.EQ.61) THEN
          TEXT = ' CONC61'
        ELSEIF(ICOMP.EQ.62) THEN
          TEXT = ' CONC62'
        ELSEIF(ICOMP.EQ.63) THEN
          TEXT = ' CONC63'
        ELSEIF(ICOMP.EQ.64) THEN
          TEXT = ' CONC64'
        ELSEIF(ICOMP.EQ.65) THEN
          TEXT = ' CONC65'
        ELSEIF(ICOMP.EQ.66) THEN
          TEXT = ' CONC66'
        ELSEIF(ICOMP.EQ.67) THEN
          TEXT = ' CONC67'
        ELSEIF(ICOMP.EQ.68) THEN
          TEXT = ' CONC68'
        ELSEIF(ICOMP.EQ.69) THEN
          TEXT = ' CONC69'
        ELSEIF(ICOMP.EQ.70) THEN
          TEXT = ' CONC70'
        ELSEIF(ICOMP.EQ.71) THEN
          TEXT = ' CONC71'
        ELSEIF(ICOMP.EQ.72) THEN
          TEXT = ' CONC72'
        ELSEIF(ICOMP.EQ.73) THEN
          TEXT = ' CONC73'
        ELSEIF(ICOMP.EQ.74) THEN
          TEXT = ' CONC74'
        ELSEIF(ICOMP.EQ.75) THEN
          TEXT = ' CONC75'
        ELSEIF(ICOMP.EQ.76) THEN
          TEXT = ' CONC76'
        ELSEIF(ICOMP.EQ.77) THEN
          TEXT = ' CONC77'
        ELSEIF(ICOMP.EQ.78) THEN
          TEXT = ' CONC78'
        ELSEIF(ICOMP.EQ.79) THEN
          TEXT = ' CONC79'
        ELSEIF(ICOMP.EQ.80) THEN
          TEXT = ' CONC80'
       ELSEIF(ICOMP.EQ.81) THEN
          TEXT = ' CONC81'
        ELSEIF(ICOMP.EQ.82) THEN
          TEXT = ' CONC82'
        ELSEIF(ICOMP.EQ.83) THEN
          TEXT = ' CONC83'
        ELSEIF(ICOMP.EQ.84) THEN
          TEXT = ' CONC84'
        ELSEIF(ICOMP.EQ.85) THEN
          TEXT = ' CONC85'
        ELSEIF(ICOMP.EQ.86) THEN
          TEXT = ' CONC86'
        ELSEIF(ICOMP.EQ.87) THEN
          TEXT = ' CONC87'
        ELSEIF(ICOMP.EQ.88) THEN
          TEXT = ' CONC88'
        ELSEIF(ICOMP.EQ.89) THEN
          TEXT = ' CONC89'
        ELSEIF(ICOMP.EQ.90) THEN
          TEXT = ' CONC90'
        ELSEIF(ICOMP.EQ.91) THEN
          TEXT = ' CONC91'
        ELSEIF(ICOMP.EQ.92) THEN
          TEXT = ' CONC92'
        ELSEIF(ICOMP.EQ.93) THEN
          TEXT = ' CONC93'
        ELSEIF(ICOMP.EQ.94) THEN
          TEXT = ' CONC94'
        ELSEIF(ICOMP.EQ.95) THEN
          TEXT = ' CONC95'
        ELSEIF(ICOMP.EQ.96) THEN
          TEXT = ' CONC96'
        ELSEIF(ICOMP.EQ.97) THEN
          TEXT = ' CONC97'
        ELSEIF(ICOMP.EQ.98) THEN
          TEXT = ' CONC98'
        ELSEIF(ICOMP.EQ.99) THEN
          TEXT = ' CONC99'
        ENDIF
      ENDIF
C
C6------RETURN.
      RETURN
      END
C
C--------------------------------------------------------------------------
      SUBROUTINE GET_TEXTM(ICOMP,TEXT)
C     ******************************************************************
C     SELECT TEXT FOR OUTPUT FILE DEPENDING ON SPECIES NUMBER
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------

      USE GWTBCTMODULE, ONLY: MCOMP,IHEAT,MCOMPT,IMULTI
      CHARACTER*16 TEXT
C     ------------------------------------------------------------------
C
      IF(IHEAT.EQ.1.AND.ICOMP.EQ.MCOMPT) THEN
        TEXT = ' TMPR '
        RETURN
      ENDIF
      IF(MCOMP.EQ.1) THEN
        TEXT = '     SOLUTE MASS'
      ELSE
        IF(ICOMP.EQ.1) THEN
          TEXT = '  SOLUTE MASS01'
        ELSEIF(ICOMP.EQ.2) THEN
          TEXT = '  SOLUTE MASS02'
        ELSEIF(ICOMP.EQ.3) THEN
          TEXT = '  SOLUTE MASS03'
        ELSEIF(ICOMP.EQ.4) THEN
          TEXT = '  SOLUTE MASS04'
        ELSEIF(ICOMP.EQ.5) THEN
          TEXT = '  SOLUTE MASS05'
        ELSEIF(ICOMP.EQ.6) THEN
          TEXT = '  SOLUTE MASS06'
        ELSEIF(ICOMP.EQ.7) THEN
          TEXT = '  SOLUTE MASS07'
        ELSEIF(ICOMP.EQ.8) THEN
          TEXT = '  SOLUTE MASS08'
        ELSEIF(ICOMP.EQ.9) THEN
          TEXT = '  SOLUTE MASS09'
        ELSEIF(ICOMP.EQ.10) THEN
          TEXT = '  SOLUTE MASS10'
        ELSEIF(ICOMP.EQ.11) THEN
          TEXT = '  SOLUTE MASS11'
        ELSEIF(ICOMP.EQ.12) THEN
          TEXT = '  SOLUTE MASS12'
        ELSEIF(ICOMP.EQ.13) THEN
          TEXT = '  SOLUTE MASS13'
        ELSEIF(ICOMP.EQ.14) THEN
          TEXT = '  SOLUTE MASS14'
        ELSEIF(ICOMP.EQ.15) THEN
          TEXT = '  SOLUTE MASS15'
        ELSEIF(ICOMP.EQ.16) THEN
          TEXT = '  SOLUTE MASS16'
        ELSEIF(ICOMP.EQ.17) THEN
          TEXT = '  SOLUTE MASS17'
        ELSEIF(ICOMP.EQ.18) THEN
          TEXT = '  SOLUTE MASS18'
        ELSEIF(ICOMP.EQ.19) THEN
          TEXT = '  SOLUTE MASS19'
        ELSEIF(ICOMP.EQ.20) THEN
          TEXT = '  SOLUTE MASS20'
        ELSEIF(ICOMP.EQ.21) THEN
          TEXT = '  SOLUTE MASS21'
        ELSEIF(ICOMP.EQ.22) THEN
          TEXT = '  SOLUTE MASS22'
        ELSEIF(ICOMP.EQ.23) THEN
          TEXT = '  SOLUTE MASS23'
        ELSEIF(ICOMP.EQ.24) THEN
          TEXT = '  SOLUTE MASS24'
        ELSEIF(ICOMP.EQ.25) THEN
          TEXT = '  SOLUTE MASS25'
        ELSEIF(ICOMP.EQ.26) THEN
          TEXT = '  SOLUTE MASS26'
        ELSEIF(ICOMP.EQ.27) THEN
          TEXT = '  SOLUTE MASS27'
        ELSEIF(ICOMP.EQ.28) THEN
          TEXT = '  SOLUTE MASS28'
        ELSEIF(ICOMP.EQ.29) THEN
          TEXT = '  SOLUTE MASS29'
        ELSEIF(ICOMP.EQ.30) THEN
          TEXT = '  SOLUTE MASS30'
        ELSEIF(ICOMP.EQ.31) THEN
          TEXT = '  SOLUTE MASS31'
        ELSEIF(ICOMP.EQ.32) THEN
          TEXT = '  SOLUTE MASS32'
        ELSEIF(ICOMP.EQ.33) THEN
          TEXT = '  SOLUTE MASS33'
        ELSEIF(ICOMP.EQ.34) THEN
          TEXT = '  SOLUTE MASS34'
        ELSEIF(ICOMP.EQ.35) THEN
          TEXT = '  SOLUTE MASS35'
        ELSEIF(ICOMP.EQ.36) THEN
          TEXT = '  SOLUTE MASS36'
        ELSEIF(ICOMP.EQ.37) THEN
          TEXT = '  SOLUTE MASS37'
        ELSEIF(ICOMP.EQ.38) THEN
          TEXT = '  SOLUTE MASS38'
        ELSEIF(ICOMP.EQ.39) THEN
          TEXT = '  SOLUTE MASS39'
        ELSEIF(ICOMP.EQ.40) THEN
          TEXT = '  SOLUTE MASS40'
        ELSEIF(ICOMP.EQ.41) THEN
          TEXT = '  SOLUTE MASS41'
        ELSEIF(ICOMP.EQ.42) THEN
          TEXT = '  SOLUTE MASS42'
        ELSEIF(ICOMP.EQ.43) THEN
          TEXT = '  SOLUTE MASS43'
        ELSEIF(ICOMP.EQ.44) THEN
          TEXT = '  SOLUTE MASS44'
        ELSEIF(ICOMP.EQ.45) THEN
          TEXT = '  SOLUTE MASS45'
        ELSEIF(ICOMP.EQ.46) THEN
          TEXT = '  SOLUTE MASS46'
        ELSEIF(ICOMP.EQ.47) THEN
          TEXT = '  SOLUTE MASS47'
        ELSEIF(ICOMP.EQ.48) THEN
          TEXT = '  SOLUTE MASS48'
        ELSEIF(ICOMP.EQ.49) THEN
          TEXT = '  SOLUTE MASS49'
        ELSEIF(ICOMP.EQ.50) THEN
          TEXT = '  SOLUTE MASSC50'
        ELSEIF(ICOMP.EQ.51) THEN
          TEXT = '  SOLUTE MASS51'
        ELSEIF(ICOMP.EQ.52) THEN
          TEXT = '  SOLUTE MASS52'
        ELSEIF(ICOMP.EQ.53) THEN
          TEXT = '  SOLUTE MASS53'
        ELSEIF(ICOMP.EQ.54) THEN
          TEXT = '   SOLUTE MASS54'
        ELSEIF(ICOMP.EQ.55) THEN
          TEXT = '   SOLUTE MASS55'
        ELSEIF(ICOMP.EQ.56) THEN
          TEXT = '  SOLUTE MASS56'
        ELSEIF(ICOMP.EQ.57) THEN
          TEXT = '  SOLUTE MASS57'
        ELSEIF(ICOMP.EQ.58) THEN
          TEXT = '  SOLUTE MASS58'
        ELSEIF(ICOMP.EQ.59) THEN
          TEXT = '  SOLUTE MASS59'
        ELSEIF(ICOMP.EQ.60) THEN
          TEXT = '  SOLUTE MASS60'
        ELSEIF(ICOMP.EQ.61) THEN
          TEXT = '  SOLUTE MASS61'
        ELSEIF(ICOMP.EQ.62) THEN
          TEXT = '  SOLUTE MASS62'
        ELSEIF(ICOMP.EQ.63) THEN
          TEXT = '  SOLUTE MASS63'
        ELSEIF(ICOMP.EQ.64) THEN
          TEXT = '  SOLUTE MASSC64'
        ELSEIF(ICOMP.EQ.65) THEN
          TEXT = '  SOLUTE MASS65'
        ELSEIF(ICOMP.EQ.66) THEN
          TEXT = '  SOLUTE MASS66'
        ELSEIF(ICOMP.EQ.67) THEN
          TEXT = '  SOLUTE MASS67'
        ELSEIF(ICOMP.EQ.68) THEN
          TEXT = '  SOLUTE MASS68'
        ELSEIF(ICOMP.EQ.69) THEN
          TEXT = '  SOLUTE MASS69'
        ELSEIF(ICOMP.EQ.70) THEN
          TEXT = '  SOLUTE MASS70'
        ELSEIF(ICOMP.EQ.71) THEN
          TEXT = '  SOLUTE MASS71'
        ELSEIF(ICOMP.EQ.72) THEN
          TEXT = '  SOLUTE MASS72'
        ELSEIF(ICOMP.EQ.73) THEN
          TEXT = '  SOLUTE MASS73'
        ELSEIF(ICOMP.EQ.74) THEN
          TEXT = '  SOLUTE MASS74'
        ELSEIF(ICOMP.EQ.75) THEN
          TEXT = '  SOLUTE MASS75'
        ELSEIF(ICOMP.EQ.76) THEN
          TEXT = '  SOLUTE MASS76'
        ELSEIF(ICOMP.EQ.77) THEN
          TEXT = '  SOLUTE MASS77'
        ELSEIF(ICOMP.EQ.78) THEN
          TEXT = '  SOLUTE MASS78'
        ELSEIF(ICOMP.EQ.79) THEN
          TEXT = '  SOLUTE MASS79'
        ELSEIF(ICOMP.EQ.80) THEN
          TEXT = '  SOLUTE MASS80'
       ELSEIF(ICOMP.EQ.81) THEN
          TEXT = '  SOLUTE MASS81'
        ELSEIF(ICOMP.EQ.82) THEN
          TEXT = '  SOLUTE MASS82'
        ELSEIF(ICOMP.EQ.83) THEN
          TEXT = '  SOLUTE MASS83'
        ELSEIF(ICOMP.EQ.84) THEN
          TEXT = '  SOLUTE MASS84'
        ELSEIF(ICOMP.EQ.85) THEN
          TEXT = '  SOLUTE MASS85'
        ELSEIF(ICOMP.EQ.86) THEN
          TEXT = '  SOLUTE MASS86'
        ELSEIF(ICOMP.EQ.87) THEN
          TEXT = '  SOLUTE MASS87'
        ELSEIF(ICOMP.EQ.88) THEN
          TEXT = '  SOLUTE MASS88'
        ELSEIF(ICOMP.EQ.89) THEN
          TEXT = '  SOLUTE MASS89'
        ELSEIF(ICOMP.EQ.90) THEN
          TEXT = '  SOLUTE MASS90'
        ELSEIF(ICOMP.EQ.91) THEN
          TEXT = '  SOLUTE MASS91'
        ELSEIF(ICOMP.EQ.92) THEN
          TEXT = '  SOLUTE MASS92'
        ELSEIF(ICOMP.EQ.93) THEN
          TEXT = '  SOLUTE MASS93'
        ELSEIF(ICOMP.EQ.94) THEN
          TEXT = '  SOLUTE MASS94'
        ELSEIF(ICOMP.EQ.95) THEN
          TEXT = '  SOLUTE MASS95'
        ELSEIF(ICOMP.EQ.96) THEN
          TEXT = '  SOLUTE MASS96'
        ELSEIF(ICOMP.EQ.97) THEN
          TEXT = '  SOLUTE MASS97'
        ELSEIF(ICOMP.EQ.98) THEN
          TEXT = '  SOLUTE MASS98'
        ELSEIF(ICOMP.EQ.99) THEN
          TEXT = '  SOLUTE MASS99'
        ENDIF
      ENDIF
C
C6------RETURN.
      RETURN
      END
C
      SUBROUTINE SGWT2BCT1IB(KSTP,KPER)
C     ******************************************************************
C     RECORD IBOUND
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,      ONLY:NCOL,NROW,NLAY,IXSEC,IBOUND,IOUT
      USE GWFBASMODULE,ONLY:PERTIM,TOTIM,IBOUUN,LBBOSV,CBOUFM,IOFLG
      USE GWTBCTMODULE, ONLY:ICBUND
C
      INTEGER,  SAVE,    DIMENSION(:,:,:),    ALLOCATABLE ::ITEMP
      CHARACTER*16 TEXT
      DATA TEXT /'          ICBUND'/
C     ------------------------------------------------------------------
C1------TRANSFER IBOUND TO TEMP LOCATION
      IF(IBOUUN.LE.0) RETURN
      ALLOCATE (ITEMP(NCOL,NROW,NLAY))
C
      N=0
      DO 58 K=1,NLAY
      DO 58 I=1,NROW
      DO 58 J=1,NCOL
      N = N+1
      ITEMP(J,I,K)=ICBUND(N)
   58 CONTINUE
C
C2------FOR EACH LAYER: SAVE ICBUND WHEN REQUESTED.
      IFIRST=1
      IF(IXSEC.EQ.0) THEN
        DO 79 K=1,NLAY
        KK=K
        IF(IOFLG(K,7).EQ.0) GO TO 79
        IF(IFIRST.EQ.1) WRITE(IOUT,74) IBOUUN,KSTP,KPER
   74   FORMAT(1X,/1X,'ICBUND WILL BE SAVED ON UNIT ',I4,
     1      ' AT END OF TIME STEP ',I8,', STRESS PERIOD ',I8)
        IFIRST=0
        CALL ULASV3(ITEMP(1,1,K),TEXT,KSTP,KPER,PERTIM,TOTIM,NCOL,
     1                NROW,KK,IBOUUN,CBOUFM,LBBOSV)
   79   CONTINUE
C
C3-----SAVE IBOUND FOR CROSS SECTION.
      ELSE
        IF(IOFLG(1,5).NE.0) THEN
          WRITE(IOUT,74) IBOUUN,KSTP,KPER
          CALL ULASV3(ITEMP(1,1,1),TEXT,KSTP,KPER,PERTIM,TOTIM,NCOL,
     1                  NLAY,-1,IBOUUN,CBOUFM,LBBOSV)
        END IF
      END IF
C
C4------RETURN.
      DEALLOCATE(ITEMP)
      RETURN
      END
C----------------------------------------------------------------------
      SUBROUTINE SGWT2BCT1CU(KSTP,KPER,IPFLG,ISA,ICOMP)
C     ******************************************************************
C     PRINT AND RECORD CONCENTRATIONS FOR UNSTRUCTURED GRIDS
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,      ONLY:NCOL,NROW,NLAY,IXSEC,HNEW,NODLAY,
     1                      IBOUND,IOUT,NODES,BUFF,IDPOUT
      USE GWFBASMODULE,ONLY:PERTIM,TOTIM,ISPCFM,ISPCUN,LBHDSV,
     2                      CSPCFM,IOFLG
      USE GWTBCTMODULE,ONLY:ICBUND,CONC,IMULTI,CROOTNAME
      DOUBLE PRECISION,  SAVE,    DIMENSION(:),    ALLOCATABLE ::BUFF8
C
      CHARACTER*400 CFILENAME
      CHARACTER*16 TEXT, ISTATUS, FORMC
C     kkz - include openspec.inc
      INCLUDE 'openspec.inc'
csp      DATA TEXT /'            CONC'/
C     ------------------------------------------------------------------
C
      ALLOCATE(BUFF8(NODES))    
      CALL GET_TEXTC(ICOMP,TEXT)
C1------FOR EACH LAYER MOVE HCONC TO BUFF IF PRINT OR SAVE IS REQUESTED.
      DO 59 K=1,NLAY
C
C2------IS CONC NEEDED FOR THIS LAYER?
      KL=K
      IF(IXSEC.NE.0) KL=1
      IF(IOFLG(KL,5).EQ.0 .AND. IOFLG(KL,6).EQ.0) GO TO 59
C
C3------MOVE HNEW TO BUFF FOR THE LAYER.
      NNDLAY = NODLAY(K)
      NSTRT = NODLAY(K-1)+1
      DO 58 N=NSTRT,NNDLAY
       BUFF8(N)=CONC(N,ICOMP)
       BUFF(N)= BUFF8(N)
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
C4A-----PRINT CONC FOR CROSS SECTION.
         ELSE
           IF(IOFLG(1,5).NE.0) THEN
             CALL ULAPRU(BUFF,TEXT,KSTP,KPER,
     1           NSTRT,NNDLAY,-1,IABS(ISPCFM),IOUT,PERTIM,TOTIM,NODES)
             IPFLG=1
C
           END IF
         END IF
      END IF
C
C5------FOR EACH LAYER: DETERMINE IF CONC SHOULD BE SAVED ON DISK.
C5------IF SO THEN CALL ULASAV OR ULASV2 TO SAVE HEAD.
      IFIRST=1
      IF(ISPCUN.LE.0) GO TO 80
C--------------------------------------------------------------------------------------------------------------------
C5A ----FOR MULTIFILE OUTPUT OPTION, ADD COMPONENT SPECIES TO ROOT NAME AND OPEN FILE ON ISPCFM AND GO TO END OF FILE
      IF(IMULTI.NE.0)THEN
        CFILENAME = trim(CROOTNAME) // TRIM(TEXT) //'.CON'
        ISTATUS = 'OLD'
        FORMC = 'FORMATTED'
        IF(CSPCFM.EQ.' ') FORMC = FORM     ! - use FORM from openspec.inc rather than hardwiring BINARY or UNFORMATTED    
CC        IF(KSTP.EQ.1.AND.KPER.EQ.1) THEN
C---------MAKE STATUS OF FILE AS NEW WHEN FIRST ENTERING FILE
CC          ISTATUS = 'REPLACE'
CC        ENDIF

csp        OPEN(ISPCUN, FILE = CFILENAME, FORM = 'BINARY',   
        OPEN(ISPCUN, FILE = CFILENAME, FORM = FORMC,
     1  STATUS = ISTATUS, POSITION = 'APPEND')
      ENDIF
C--------------------------------------------------------------------------------------------------------------------
      IF(IXSEC.EQ.0) THEN
        DO 79 K=1,NLAY
        KK=K
        IF(IOFLG(K,6).EQ.0) GO TO 79
        NNDLAY = NODLAY(K)
        NSTRT = NODLAY(K-1)+1
        IF(IFIRST.EQ.1) WRITE(IOUT,74) ISPCUN,KSTP,KPER
   74   FORMAT(1X,/1X,'CONC WILL BE SAVED ON UNIT ',I4,
     1      ' AT END OF TIME STEP ',I8,', STRESS PERIOD ',I8)
        IFIRST=0
        IF(IDPOUT.EQ.1)THEN 
          WRITE(IOUT,*)'  CONC AND TIME VARIABLES ARE SAVED AS REAL*8'
          CALL ULASAVU8(BUFF8,TEXT,KSTP,KPER,PERTIM,TOTIM,NSTRT,
     1                NNDLAY,KK,ISPCUN,NODES)
        ELSEIF(CSPCFM.EQ.' ') THEN
           CALL ULASAVU(BUFF,TEXT,KSTP,KPER,PERTIM,TOTIM,NSTRT,
     1                NNDLAY,KK,ISPCUN,NODES)
        ELSE
           CALL ULASV2U(BUFF,TEXT,KSTP,KPER,PERTIM,TOTIM,NSTRT,
     1             NNDLAY,KK,ISPCUN,CSPCFM,LBHDSV,ICBUND(NSTRT),NODES)
        END IF
   79   CONTINUE
C
C5A-----SAVE CONC FOR CROSS SECTION.
      ELSE
        IF(IOFLG(1,6).NE.0) THEN
          WRITE(IOUT,74) ISPCUN,KSTP,KPER
          IF(IDPOUT.EQ.1) THEN 
            WRITE(IOUT,*)'  CONC AND TIME VARIABLES ARE SAVED AS REAL*8'
            CALL ULASAVU8(BUFF8,TEXT,KSTP,KPER,PERTIM,TOTIM,NSTRT,
     1                NNDLAY,-1,ISPCUN,NODES)
          ELSEIF(CSPCFM.EQ.' ') THEN
             CALL ULASAVU(BUFF,TEXT,KSTP,KPER,PERTIM,TOTIM,NSTRT,
     1                NNDLAY,-1,ISPCUN,NODES)
          ELSE
             CALL ULASV2U(BUFF,TEXT,KSTP,KPER,PERTIM,TOTIM,NSTRT,
     1                  NNDLAY,-1,ISPCUN,CSPCFM,LBHDSV,ICBUND,NODES)
          END IF
        END IF
      END IF
C
C6------RETURN.
   80 CONTINUE
      DEALLOCATE(BUFF8)
      RETURN
C
      END
      SUBROUTINE SGWT2BCT1IBU(KSTP,KPER)
C     ******************************************************************
C     RECORD IBOUND FOR UNSTRUCTURED GRIDS
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,     ONLY:NCOL,NROW,NLAY,IXSEC,IOUT,NODLAY,NODES
      USE GWFBASMODULE,ONLY:PERTIM,TOTIM,IBOUUN,LBBOSV,CBOUFM,IOFLG
      USE GWTBCTMODULE,ONLY:ICBUND
C
      CHARACTER*16 TEXT
      DATA TEXT /'          IBOUND'/
C     ------------------------------------------------------------------
      IF(IBOUUN.LE.0) RETURN
C
C1------FOR EACH LAYER: SAVE ICBUND WHEN REQUESTED.
      IFIRST=1
      IF(IXSEC.EQ.0) THEN
        DO 79 K=1,NLAY
        KK=K
        IF(IOFLG(K,7).EQ.0) GO TO 79
        NNDLAY = NODLAY(K)
        NSTRT = NODLAY(K-1)+1
C
        IF(IFIRST.EQ.1) WRITE(IOUT,74) IBOUUN,KSTP,KPER
   74   FORMAT(1X,/1X,'IBOUND WILL BE SAVED ON UNIT ',I4,
     1      ' AT END OF TIME STEP ',I8,', STRESS PERIOD ',I8)
        IFIRST=0
        CALL ULASV3U(ICBUND,TEXT,KSTP,KPER,PERTIM,TOTIM,NSTRT,
     1                NNDLAY,KK,IBOUUN,CBOUFM,LBBOSV,NODES)
   79   CONTINUE
C
C1A-----SAVE ICBUND FOR CROSS SECTION.
      ELSE
        IF(IOFLG(1,7).NE.0) THEN
          WRITE(IOUT,74) IBOUUN,KSTP,KPER
          CALL ULASV3U(ICBUND,TEXT,KSTP,KPER,PERTIM,TOTIM,NSTRT,
     1                  NNDLAY,-1,IBOUUN,CBOUFM,LBBOSV,NODES)
        END IF
      END IF
C
C6------RETURN.
      RETURN
      END
C----------------------------------------------------------------------
      SUBROUTINE SCLN1C(KSTP,KPER,IPFLG,ISA,ICOMP)
C     ******************************************************************
C     PRINT AND SAVE CONCS IN CLN CELLS
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,      ONLY:HNEW,IBOUND,IOUT,NODES,IDPOUT
      USE CLN1MODULE, ONLY:  NCLNNDS,ICLNHD,ICLNCN
      USE GWFBASMODULE,ONLY:PERTIM,TOTIM,IHEDFM,IHEDUN,LBDDSV,
     2                      CHEDFM,CDDNFM,IOFLG,ISPCFM,CSPCFM,ISPCUN
      USE GWTBCTMODULE, ONLY: ICBUND,CONC,IMULTI,CROOTNAME
C
      CHARACTER*400 CFILENAME
      CHARACTER*16 TEXT,ISTATUS, FORMC
      DOUBLE PRECISION SSTRT
      REAL,          SAVE,    DIMENSION(:),    ALLOCATABLE ::BUFF
      REAL*8,        SAVE,    DIMENSION(:),    ALLOCATABLE ::BUFF8
C
      INCLUDE 'openspec.inc'
C
C     ------------------------------------------------------------------
      CALL GET_TEXTC(ICOMP,TEXT)
      ALLOCATE(BUFF(NCLNNDS))
      ALLOCATE(BUFF8(NCLNNDS))
C
C1------FOR EACH CLN NODE PUT CONCS IN BUFF IF PRINT OR SAVE IS REQUESTED.
      DO 59 N=1,NCLNNDS
C
C2------Save CONCS in buffer array BUFF
        NG = N+NODES
        BUFF8(N)=CONC(NG,ICOMP)
        BUFF(N)=BUFF8(N)
   59 CONTINUE
C
C3------CALL ULAPRS OR ULAPRW TO PRINT CONCS.
      IF(ISA.NE.0) THEN
        IF(IOFLG(1,5).NE.0) THEN
          IF(ISPCFM.LT.0) CALL ULAPRS(BUFF(1),TEXT,KSTP,KPER,
     1                  NCLNNDS,1,1,-ISPCFM,IOUT)
          IF(ISPCFM.GE.0) CALL ULAPRW(BUFF(1),TEXT,KSTP,KPER,
     1                  NCLNNDS,1,1,ISPCFM,IOUT)
          IPFLG=1
        ENDIF
C
      END IF
C
C4------DETERMINE IF CONC SHOULD BE SAVED.
C4------IF SO THEN CALL A ULASAV OR ULASV2 TO RECORD CONC.
      IFIRST=1
      IF(ICLNCN.LE.0) GO TO 80
        NSTRT = NODES+1
        IF(IOFLG(1,6).EQ.0) GO TO 80
C--------------------------------------------------------------------------------------------------------------------
C5A ----FOR MULTIFILE OUTPUT OPTION, ADD COMPONENT SPECIES TO ROOT NAME AND OPEN FILE ON ISPCFM AND GO TO END OF FILE
      IF(IMULTI.NE.0)THEN
        CFILENAME = trim(CROOTNAME) // TRIM(TEXT) //'_CLN.CON'
        ISTATUS = 'OLD'
        FORMC = 'FORMATTED'
        IF(CSPCFM.EQ.' ') FORMC = FORM        ! - use FORM from openspec.inc rather than hardwiring BINARY or UNFORMATTED     
CC        IF(KSTP.EQ.1.AND.KPER.EQ.1) THEN
C---------MAKE STATUS OF FILE AS NEW WHEN FIRST ENTERING FILE
CC          ISTATUS = 'REPLACE'
CC        ENDIF
C
!        OPEN(ICLNCN, FILE = CFILENAME, FORM = 'BINARY',   
        OPEN(ICLNCN, FILE = CFILENAME, FORM = FORMC,
     1  STATUS = ISTATUS, POSITION = 'APPEND')
      ENDIF
C--------------------------------------------------------------------------------------------------------------------
        IF(IFIRST.EQ.1) WRITE(IOUT,74) ICLNCN,KSTP,KPER
   74   FORMAT(1X,/1X,'CLN CONC WILL BE SAVED ON UNIT ',I4,
     1      ' AT END OF TIME STEP ',I3,', STRESS PERIOD ',I4)
        IFIRST=0
        IF(IDPOUT.EQ.1) THEN 
          WRITE(IOUT,*)' CLNCONC AND TIME VARIABLES ARE SAVED AS REAL*8'
            CALL ULASAV8(BUFF8(1),TEXT,KSTP,KPER,PERTIM,TOTIM,NCLNNDS,
     1                1,1,ICLNCN)
        ELSEIF(CSPCFM.EQ.' ') THEN
           CALL ULASAV(BUFF(1),TEXT,KSTP,KPER,PERTIM,TOTIM,NCLNNDS,
     1                1,1,ICLNCN)
        ELSE
           CALL ULASV2(BUFF(1),TEXT,KSTP,KPER,PERTIM,TOTIM,NCLNNDS,
     1                1,1,ICLNCN,CSPCFM,LBHDSV,ICBUND(NSTRT))
        END IF
C
80    CONTINUE
      DEALLOCATE(BUFF)
      DEALLOCATE(BUFF8)

C
C5------RETURN.
      RETURN
      END
C------------------------------------------------------------------------
      SUBROUTINE GWT2BCT1DA
      USE GWTBCTMODULE
      USE AW_ADSORBMODULE
      USE GLOBAL, ONLY: ARAD,IOUTNORMAL
      INTEGER ALLOC_ERR
C
      DEALLOCATE(IBCTCB, STAT = ALLOC_ERR)
      DEALLOCATE(MCOMP, STAT = ALLOC_ERR)
      DEALLOCATE(ITVD, STAT = ALLOC_ERR)
      DEALLOCATE(CINACT, STAT = ALLOC_ERR)
      DEALLOCATE(CICLOSE, STAT = ALLOC_ERR)
      DEALLOCATE(ICBUND, STAT = ALLOC_ERR)
      DEALLOCATE(PRSITY, STAT = ALLOC_ERR)
      DEALLOCATE(CONC, STAT = ALLOC_ERR)
      DEALLOCATE(CONCO, STAT = ALLOC_ERR)
      DEALLOCATE(MASSBCT, STAT = ALLOC_ERR)
      DEALLOCATE(ADMAT, STAT = ALLOC_ERR)
      DEALLOCATE(IPCBFLAG, STAT = ALLOC_ERR)
      DEALLOCATE(CBCF, STAT = ALLOC_ERR)
      DEALLOCATE(CBCH, STAT = ALLOC_ERR)
      DEALLOCATE (Cncg,LrcC, STAT = ALLOC_ERR)
c      IF(IADSORB.NE.0) DEALLOCATE(ADSORB, STAT = ALLOC_ERR)
c      IF(IADSORB.EQ.2) DEALLOCATE(FLICH, STAT = ALLOC_ERR)
      DEALLOCATE(ADSORB, STAT = ALLOC_ERR)
      DEALLOCATE(FLICH, STAT = ALLOC_ERR)
c      IF(IDISP.NE.0)THEN
       DEALLOCATE(VELNOD, STAT = ALLOC_ERR)
       DEALLOCATE(VXL,VYL, STAT = ALLOC_ERR)
       DEALLOCATE(DLX, DLY, DLZ, ATXY,ATYZ, ATXZ, STAT = ALLOC_ERR)
c       IF(IXDISP.NE.0) DEALLOCATE(ALXY, ALYZ, ALXZ)
       DEALLOCATE(ALXY, ALYZ, ALXZ, STAT = ALLOC_ERR)
c      ENDIF
c      IF(IZOD.EQ.1)THEN
c        DEALLOCATE(ZODRW, STAT = ALLOC_ERR)
c      ELSEIF(IADSORB.GT.0)THEN
c        IF(IZOD.EQ.2)THEN
c           DEALLOCATE(ZODRS, STAT = ALLOC_ERR)
c        ELSEIF(IZOD.EQ.3)THEN
          DEALLOCATE(ZODRW, STAT = ALLOC_ERR)
          DEALLOCATE(ZODRS, STAT = ALLOC_ERR)
c        ENDIF
c      ENDIF
c      IF(IFOD.EQ.1)THEN
c        DEALLOCATE(FODRW, STAT = ALLOC_ERR)
c      ELSEIF(IADSORB.GT.0)THEN
c        IF(IFOD.EQ.2)THEN
c           DEALLOCATE(FODRS, STAT = ALLOC_ERR)
c        ELSEIF(IFOD.EQ.3)THEN
          DEALLOCATE(FODRW, STAT = ALLOC_ERR)
          DEALLOCATE(FODRS, STAT = ALLOC_ERR)
c        ENDIF
c      ENDIF
      DEALLOCATE(IADSORB, STAT = ALLOC_ERR)
      DEALLOCATE(IDISP, STAT = ALLOC_ERR)
      DEALLOCATE(IXDISP, STAT = ALLOC_ERR)
      DEALLOCATE(DIFFNC, STAT = ALLOC_ERR)
      DEALLOCATE(VBVLT, STAT = ALLOC_ERR)
      DEALLOCATE(VBNMT, STAT = ALLOC_ERR)
      DEALLOCATE(MSUMT, STAT = ALLOC_ERR)
      DEALLOCATE(IFOD, IZOD, STAT = ALLOC_ERR)
      DEALLOCATE(FCncg, STAT = ALLOC_ERR)
      DEALLOCATE(FLrcC, STAT = ALLOC_ERR)
      DEALLOCATE(ARAD, STAT = ALLOC_ERR)
      IF(ISOLUBILITY.EQ.1) THEN
        DEALLOCATE(SOLLIM, SOLSLOPE, STAT = ALLOC_ERR)
      ENDIF
      IF(ISOLUBILITY.EQ.2) THEN
        DEALLOCATE(SOLLIMVAR, SOLSLOPEVAR, STAT = ALLOC_ERR)
      ENDIF      
      DEALLOCATE (ISOLUBILITY, STAT = ALLOC_ERR)
      DEALLOCATE (IOUTNORMAL, STAT = ALLOC_ERR)
      IF(IAW_ADSORB.EQ.1) THEN 
         DEALLOCATE (IAREA_FN, STAT = ALLOC_ERR)
         DEALLOCATE (IKAWI_FN, STAT = ALLOC_ERR)
         DEALLOCATE (ITAB_AWI, STAT = ALLOC_ERR)
         DEALLOCATE (ZODRAW, STAT = ALLOC_ERR)
         DEALLOCATE (FODRAW, STAT = ALLOC_ERR)
         
         CALL GWT2AWI_ADSORB1DA
      ENDIF
       DEALLOCATE (IAW_ADSORB, STAT = ALLOC_ERR)
C
      RETURN
      END
C----------------------------------------------------------------------
      SUBROUTINE GWT2HEATCOND
C     ******************************************************************
C     FORMULATE HEAT CONDUCTANCE TERM FOR EACH NODE IN ARRAY HTCONDT FOR GROUNDWATER NODES
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,  ONLY:IOUT,NCOL,NROW,NLAY,IBOUND,TOP,NODLAY,IVC,ISYM,
     1  IA,JA,JAS,ARAD,NJA,NODES,FAHL,TOP,BOT,CL1,CL2,Sn,NEQS,AREA,INCLN
      USE GWFBASMODULE,ONLY:DELT
      USE GWTBCTMODULE, ONLY: MCOMP,CONC,ICBUND,ITVD,PRSITY,
     *  HTCONDS,HTCONDW,HTCONDT,IDISPCLN,HTCAPW,RHOW
      USE CLN1MODULE, ONLY: ACLNNDS,IBHETYP,BHEPROP
      DOUBLE PRECISION VOLN,VOLJJ,SATN,SATJJ,SATC,HTCN,HTCJJ,HTC,ARE,
     1                THIKN,THIKJJ,ALNTHN,ALNTHJJ,DIV
C     ------------------------------------------------------------------
      DO N=1,NEQS
        IFTYP1 = 0     ! INDEX OF CELL TYPE (=0 FOR GWF; >0 FOR CLN TYPE IDENTIFIER)
        IF(ICBUND(N).NE.0)THEN
C1-----COMPUTE HEAT CONDUCTIVITY FOR CELL N
          SATN = Sn(N)
          IF(N.LE.NODES) THEN
            HTCN = HTCONDS(N) + PRSITY(N) * SATN * HTCONDW
            THIKN = TOP(N) - BOT(N)
            VOLN = THIKN * AREA(N)
          ELSE
            NC1 = N-NODES
            IFTYP1 = ACLNNDS(NC1,2)     ! (=0 FOR GWF; >0 FOR CLN TYPE IDENTIFIER)
            IF(IBHETYP.EQ.0) THEN
              HTCN = HTCONDW * SATN   ! CONDUIT HAS ONLY WATER
            ELSE
              HTFLUID = BHEPROP(IFTYP1,3)/HTCAPW/RHOW
              IF(HTFLUID.LT. 0.0) HTFLUID = HTCONDW
              HTCN = HTFLUID * SATN   ! CONDUIT HAS ONLY WATER
            ENDIF
            VOLN = AREA(N) * ACLNNDS(NC1,4)
          ENDIF
C          IF(N.LE.NODES) THEN
C            HTCN = HTCN / PRSITY(N)    ! DIFFUSION TIMES WATER CONTENT IS EQUIVALENCED
C          ENDIF    
C
C2-----LOOP OVER CONNECTIONS OF NODE N AND FILL
          DO II = IA(N)+1,IA(N+1)-1
            ICMC = 0   ! INDICATES A CONDUIT TO MATRIX CONNECTION
            IFTYP2 = 0 ! INDEX OF CELL TYPE (=0 FOR GWF; >0 FOR CLN TYPE IDENTIFIER)
            JJ = JA(II)
            IF(JJ.LE.N) CYCLE !ONLY FOR UPPER TRIANGLE
            IF(ICBUND(JJ).NE.0)THEN
              IIS = JAS(II)
C3-----COMPUTE HEAT CONDUCTIVITY FOR CELL JJ
              SATJJ = Sn(JJ)
              ALNTHN = CL1(IIS)
              ALNTHJJ = CL2(IIS)
              IF(JJ.LE.NODES) THEN
                HTCJJ = HTCONDS(JJ) + PRSITY(JJ) * SATJJ * HTCONDW
                THIKJJ = TOP(JJ) - BOT(JJ)
                VOLJJ = THIKJJ * AREA(JJ)
              ELSE
                NC2 = JJ-NODES
                IFTYP2 = ACLNNDS(NC2,2)          ! (=0 FOR GWF; >0 FOR CLN TYPE IDENTIFIER)
                IF(IBHETYP.EQ.0) THEN
                  HTCJJ = HTCONDW *SATJJ   ! CONDUIT HAS ONLY WATER
                ELSE
                  HTFLUID = BHEPROP(IFTYP2,3)/HTCAPW/RHOW
                  IF(HTFLUID.LT. 0.0) HTFLUID = HTCONDW
                  HTCJJ = HTFLUID * SATN   ! CONDUIT HAS ONLY WATER
                ENDIF
                VOLN = AREA(JJ) * ACLNNDS(NC2,4)
              ENDIF
C              IF(JJ.LE.NODES) THEN
C                HTCJJ = HTCJJ / PRSITY(JJ)    ! DIFFUSION TIMES WATER CONTENT IS EQUIVALENCED
C              ENDIF               
C-----COMPUTE VOLUME WEIGHTED AVERAGE HEAT CONDUCTIVITY FOR CELLS N AND JJ
C              HTC = HTCN*HTCJJ * (VOLN + VOLJJ)
C              HTC = HTC / (HTCN*VOLJJ + HTCJJ*VOLN)
C4--------FLAG CONDUIT-MATRIX CONNECTION
              IF(( N.GT.NODES.AND.JJ.LE.NODES).OR.
     1          ( JJ.GT.NODES.AND. N.LE.NODES)) THEN
                ICMC = 1     !INDEX FOR CONDUIT-MATRIX CONNECTION IS ON
              ENDIF
              IF(INCLN.NE.0) THEN 
                IF(IBHETYP.EQ.0) ICMC = 0 !SHUT OFF IF NO SPECIAL BHE
              ENDIF  
C5-----COMPUTE LENGTH WEIGHTED AVERAGE HEAT CONDUCTIVITY FOR CELLS N AND JJ
              IF(ICMC.EQ.0)THEN
                HTC = HTCN*HTCJJ * (ALNTHN + ALNTHJJ)
                DIV = (HTCN*ALNTHJJ + HTCJJ*ALNTHN)
                IF(DIV. GT. 1.0E-20) HTC = HTC / DIV
              ELSE
                IFTYP = IFTYP1
                IF(IFTYP1.EQ.0) IFTYP = IFTYP2
                HTP = BHEPROP(IFTYP,2) / BHEPROP(IFTYP,1) *HTCAPW * RHOW      ! HEAT RESISTANCE OF PIPE MATERIAL
                HTC=ALNTHN/HTCN + ALNTHJJ/HTCJJ + HTP   ! ADD UP ALL HEAT RESISTANCES
                HTC = (ALNTHN + ALNTHJJ + BHEPROP(IFTYP,2))/ HTC   !CONVERT RESISTANCE TO CONDUCTANCE
                HTC = HTC + BHEPROP(IFTYP,4)/HTCAPW/RHOW              ! ADD CONVECTIVE CONDUCTANCE COEFFICIENT
              ENDIF
C6---------CONVERT HEAT CONDUCTIVITY TO HEAT CONDUCTANCE
              IF(IDISPCLN.EQ.1. AND. ICMC.EQ.1) THEN
                FRAD = CL1(IIS)
                RO = CL2(IIS)
                HTC = HTC * FAHL(IIS) / FRAD /LOG(RO/FRAD)
              ELSE
                IF(IVC(IIS).EQ.1)THEN !VERTICAL DIRECTION CONNECTION
                  ARE = FAHL(IIS)
                ELSEIF(JJ.GT.NODES.OR.N.GT.NODES)THEN !CONDUIT-MATRIX OR CONDUIT-CONDUIT CONNECTION
                  ARE = FAHL(IIS)
                ELSE
                  ARE = FAHL(IIS) * (THIKN + THIKJJ) * 0.5
                ENDIF
                HTC = HTC * ARE / (CL1(IIS) + CL2(IIS))
              ENDIF
C
C7---------ADD TERM TO MATRIX FOR THIS NODE AND ITS EFFECT ON CONNECTING NODE JJ
              HTCONDT(II) = HTCONDT(II) + HTC
              HTCONDT(IA(N)) =  HTCONDT(IA(N)) - HTC
              HTCONDT(ISYM(II)) = HTCONDT(ISYM(II)) + HTC
              HTCONDT(IA(JJ)) = HTCONDT(IA(JJ)) - HTC
            ENDIF
          ENDDO
        ENDIF
      ENDDO
C
C8-----RETURN
      RETURN
      END
C----------------------------------------------------------------------
      SUBROUTINE GWT2HEATCONDCLN
C     ******************************************************************
C     FORMULATE HEAT CONDUCTANCE TERM FOR EACH NODE IN ARRAY HTCONDT FOR CLN NODES
C     ******************************************************************
      USE GLOBAL, ONLY:NODES,NJA,IA,PGF,FAHL,TOP,BOT,CL1,CL2,NODES,NLAY,
     1            NODLAY,IBOUND,JA,JAS,IVC,ISYM,AREA,IDEALLOC_HY,Sn
      USE CLN1MODULE, ONLY: ACLNNDS,NCLNNDS,NNDCLN,ACLNGWC,NCLNGWC,
     1      ACLNREC,ACLNCOND,NCLN,CLNCON,IA_CLN,JA_CLN,IDXGLO_CLN
      USE GWFBCFMODULE,ONLY:HK,CV,LAYCON,LAYAVG,IKVFLAG
      USE GWTBCTMODULE, ONLY: HTCONDT,HTCONDW,IDISPCLN,PRSITY,HTCONDS
      DOUBLE PRECISION AREAF,EL,RADFSQ,RO,ROD,CWND,DEX,DEXN,DEZ,
     1  AREAF1,AREAF2,FPER,FRAD,FANISO,EL1,EL2,AKVAL,PERIF,SATN,SATJJ,
     1  HTCN,HTCJJ,HTC,VOLN,VOLJJ,ARE
C
C--------------------------------------------------------------------------------------
C1-----CONNECT CLN NODES TO EACH OTHER
C--------------------------------------------------------------------------------------
      PI = 3.1415926
C-----------------------------------------------------------------------------------
C1A---------loop over all CLN nodes
      DO NC1 = 1,NCLNNDS
C1-----COMPUTE HEAT CONDUCTANCE FOR CLN CELL N
        N = ACLNNDS(NC1,1)
        SATN = Sn(N)
C        HTCN =  SATN * HTCONDW
C2----------loop over all connections of node NC1
        DO II_CLN = IA_CLN(NC1)+1,IA_CLN(NC1+1)-1
          NC2 = JA_CLN(II_CLN)
          IF(NC2.GT.NC1) CYCLE
          JJ = ACLNNDS(NC2,1)
          IC1 = ACLNNDS(NC1,2)     !CLN TYPE FOR NODE 1
          IC2 = ACLNNDS(NC2,2)     !CLN TYPE FOR NODE 2
          II = IDXGLO_CLN(II_CLN)
          IIS = JAS(II)
C3-----COMPUTE HEAT CONDUCTANCE FOR CELL JJ
C          SATJJ = Sn(JJ)
C          HTCJJ = SATJJ * HTCONDW
C4-----COMPUTE HEAT CONDUCTANCE FOR FACE BETWEEN CELLS N AND JJ
C          VOLN = AREA(N) * ACLNNDS(NC1,4)
C          VOLJJ = AREA(JJ) * ACLNNDS(NC2,4)
C          HTC = HTCN*HTCJJ * (VOLN + VOLJJ)
C          HTC = HTC / (HTCN*VOLJJ + HTCJJ*VOLN)
          HTC = HTCONDW  ! NO NEED TO HARMONIC AVERAGE SINCE ONLY WATER IN CLN.
C5-----COMPUTE THE CONDUCTANCE TERM BETWEEN CELLS N AND JJ
          ARE = FAHL(IIS)
          HTC = HTC * ARE / (CL1(IIS) + CL2(IIS))
C9---------ADD TERM TO MATRIX FOR THIS NODE AND ITS EFFECT ON CONNECTING NODE JJ
          HTCONDT(II) =  HTC
          HTCONDT(IA(N)) =  - HTC
          HTCONDT(ISYM(II)) =  + HTC
          HTCONDT(IA(JJ)) =  - HTC
        ENDDO
      ENDDO
C
C--------------------------------------------------------------------------------------
C4-----CONNECT CLN NODES WITH POROUS MATRIX
C-------------------------------------------------------------------------------------
C4A-----loop over all conduit node to GW connections
      DO IFN = 1,NCLNGWC
        IH = ACLNGWC(IFN,1)
        NH = ACLNNDS(IH,1)
        NL = ACLNGWC(IFN,2)
        SATN = Sn(NL)
        HTCN = HTCONDS(NL) + PRSITY(NL) * SATN * HTCONDW
        DO II = IA(NL)+1,IA(NL+1)-1
          JJ = JA(II)
          IF(JJ.NE.NH) CYCLE
          IIS = JAS(II)
          SATJJ = Sn(JJ)
          HTCJJ = SATJJ * HTCONDW    ! CONDUIT ONLY HAS WATER IN IT
C4-----COMPUTE HEAT CONDUCTANCE FOR FACE BETWEEN CELLS N AND JJ
          VOLN = (TOP(NL) - BOT(NL)) * AREA(NL)
          VOLJJ = AREA(JJ) * ACLNNDS(IH,4)
          HTC = HTCN*HTCJJ * (VOLN + VOLJJ)
          HTC = HTC / (HTCN*VOLJJ + HTCJJ*VOLN)
C--------------------------------------------------------------------------------------
C4B---------FIRST COMPUTE EFFECTIVE CELL RADIUS FOR THIS CONNECTION
          IFNC = ACLNGWC(IFN,1)
          IFDIR = ACLNNDS(IFNC,3)
          IF(IFDIR.EQ.0)THEN
C4B1---------GET RO FOR VERTICAL WELL USING ISOTROPIC THIEM EQUATION,
            RO = 0.0
            DEXN = 0.0
            DO IJ = IA(NL)+1,IA(NL+1)-1
              JJ1 = JA(IJ)
              IF(JJ1.GT.NODES)CYCLE
              IJS = JAS(IJ)
              IF(IVC(IJS).EQ.1) CYCLE
              RO = RO + CL1(IJS)**2
              DEXN = DEXN + 1.0
            ENDDO
            if(dexn.lt.1.0e-10) dexn = 1.0
            if(ro.lt.1.0e-10) ro = 1.0
            RO = 0.28 * SQRT(2.0*RO/DEXN)
          ELSE
C4B2--------GET RO FOR HORIZONTAL WELL USING ANISOTROPIC EQUATION
            DEX = 0.0
            DEXN = 0.0
            DO IJ = IA(NL)+1,IA(NL+1)-1
              JJ2 = JA(IJ)
              IF(JJ2.GT.NODES)CYCLE
              IJS = JAS(IJ)
              IF(IVC(IJS).EQ.1) CYCLE
              DEX = DEX + CL1(IJS)
              DEXN = DEXN + 1.0
            ENDDO
            DEX = 2.0 * DEX / DEXN
            DEZ = TOP(NL) - BOT(NL)
            RO = DEX**2 + DEZ**2
            RO = 0.28 * SQRT(RO)
          ENDIF
C--------------------------------------------------------------------------------------
C5---------COMPUTE CONDUCTANCE TERM FOR THE DIFFERENT CLN-MATRIX CONNECTION TYPES
          IFTYP = ACLNNDS(IFNC,2)
          FLENG = ACLNGWC(IFN,6)
          CALL CLNR(IFTYP,FRAD)
          IF(IDISPCLN.EQ.0)THEN
C
C5A-----------CONNECTION IS ACROSS A LEAKANCE TERM LIKE CONDUIT FLOW PROCESS OF MF2K5, COMPUTE LEAKANCE
            ARE = FAHL(IIS)
            HTC =  HTC * ARE / (CL1(IIS) + CL2(IIS))
          ELSEIF(IDISPCLN.EQ.1)THEN
C
C5C-----------CONNECTION USES THIEM EQUATION LIKE MULTI-NODE WELL PACKAGE
C5C2----------COMPUTE EFFECTIVE WELL CONDUCTIVITY
            CWND = LOG(RO / FRAD)
C5C3------------COMPUTE THE CONDUCTANCE TERM

            HTC = 2.0*PI*HTC * FLENG / CWND
          ENDIF
C9---------ADD TERM TO MATRIX FOR THIS NODE AND ITS EFFECT ON CONNECTING NODE JJ
          HTCONDT(II) =  HTCONDT(II) + HTC
          HTCONDT(IA(IFN)) = HTCONDT(IA(IFN)) - HTC
          HTCONDT(ISYM(II)) = HTCONDT(ISYM(II)) + HTC
          HTCONDT(IA(JJ)) = HTCONDT(IA(JJ)) - HTC
        ENDDO
      ENDDO
C
C7------RETURN
      RETURN
      END
C----------------------------------------------------------------------
      SUBROUTINE SOLUBILITY(N,ICOMP,DTERMS,RTERMS,ISS)

C     ******************************************************************
C     FORMULATE SOLUBILITY TERM FOR EACH NODE AND COMPONENT
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,   ONLY:IOUT,NCOL,NROW,NLAY,IBOUND,TOP,NODLAY,ISYM,
     1 HNEW,AMAT,RHS,HOLD,ISSFLG,IA,JA,NJA,NODES,AREA,BOT,TOP,Sn,So
      USE CLN1MODULE, ONLY: ACLNNDS
      USE GWFBASMODULE,ONLY:DELT
      USE GWTBCTMODULE, ONLY: MCOMP,CONC,CONCO,ICT,IADSORB,PRSITY,
     *  ADSORB,FLICH,IZOD,IFOD,ZODRW,FODRW,ZODRS,FODRS,MCOMPT,IHEAT,
     *  ISOLUBILITY,SOLLIM,SOLSLOPE,SOLLIMVAR,SOLSLOPEVAR
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
          IF(ISOLUBILITY .EQ. 1) THEN
            SSLOPE = SOLSLOPE(ICOMP)
            SLIM = SOLLIM(ICOMP)
          ELSEIF(ISOLUBILITY .EQ. 2) THEN
            SSLOPE = SOLSLOPEVAR(N,ICOMP)
            SLIM = SOLLIMVAR(N,ICOMP)  
          ENDIF    
          CW = CONC(N,ICOMP)
          CWO = CONCO(N,ICOMP)
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
        RTERMS = RTERMS - VODT * CONCO(N,ICOMP)
      ENDIF
C
C9------RETURN
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE SOLUBILITYMB(KSTP,KPER,ICOMP,ISS)
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
      USE GWTBCTMODULE, ONLY: ICBUND,CONC,MSUMT,VBVLT,VBNMT,IBCTCB,
     1 IADSORB,ADSORB,FLICH,PRSITY,CONCO,ICT,MCOMPT,
     *  ISOLUBILITY,SOLLIM,SOLSLOPE,SOLLIMVAR,SOLSLOPEVAR
C
      CHARACTER*16 TEXT(2)
      DOUBLE PRECISION RATIN,RATOUT,QQ,VODT,ADSTERM,FL,CW,CWO,ALENG,
     *  DTERMS,RTERMS,VOLU,RATINTVM,RATOUTTVM,
     *  EMNEW,EMOLD,FN,SLIM,SSLOPE
      DATA TEXT(1) /'     MASS PRECIP'/
      DATA TEXT(2) /' CLN MASS PRECIP'/
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
      DO 50 N=1,NEQS
      BUFF(N)=ZERO
50    CONTINUE
C
C3------LOOP THROUGH EACH NODE AND CALCULATE STORAGE
      DO 100 N=1,NEQS
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
C          SSLOPE = SOLSLOPE(ICOMP) * VODT
C          SLIM = SOLLIM(ICOMP)
C          CW = CONC(N,ICOMP)
C          CWO = CONCO(N,ICOMP)
C          IF(CONC(N,ICOMP).LT.SLIM) CW = SLIM
C          IF(CONCO(N,ICOMP).LT.SLIM) CWO = SLIM
C          QQ = SSLOPE * (CW - CWO)
CC
csp          SSLOPE = SOLSLOPE(ICOMP)
csp          SLIM = SOLLIM(ICOMP)
csp          CW = CONC(N,ICOMP)
csp          CWO = CONCO(N,ICOMP)
csp          EMOLD = 0.0
csp          IF(CWO.GT.SLIM) EMOLD = (CWO-SLIM) * SSLOPE
csp          EMNEW = 0.0
csp          IF(CW.GT.SLIM) EMNEW = (CW-SLIM) * SSLOPE
csp          FN = (EMNEW - EMOLD) * VODT
csp          QQ = FN
cc
          IF(ISOLUBILITY .EQ. 1) THEN
            SSLOPE = SOLSLOPE(ICOMP)
            SLIM = SOLLIM(ICOMP)
          ELSEIF(ISOLUBILITY .EQ. 2) THEN
            SSLOPE = SOLSLOPEVAR(N,ICOMP)
            SLIM = SOLLIMVAR(N,ICOMP)  
          ENDIF
          CW = CONC(N,ICOMP)
          CWO = CONCO(N,ICOMP)
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
          qq = rterms - dterms * conc(n,icomp)
252     CONTINUE
C-----------------------------------------------------
      ELSE   !-----------------------TOTAL CONCENTRATION FORMULATION
C7-------NET STORAGE TERM FOR TOTAL CONCENTRATION FORMULATION
        QQ = QQ + VODT * CONC(N,ICOMP)
     *          - VODT * CONCO(N,ICOMP)
      ENDIF
      QQ = - QQ  ! STORAGE TERM NEGATIVE IS INFLOW AS PER MODFLOW CONVENTION
      Q = QQ
C
C8------PRINT FLOW RATE IF REQUESTED.
      IF(IBD.LT.0) THEN
         IF(IBDLBL.EQ.0) WRITE(IOUT,61) TEXT(1),KPER,KSTP
   61    FORMAT(1X,/1X,A,'   PERIOD ',I4,'   STEP ',I3)
        IF(IUNSTR.EQ.0.AND.N.LE.NODES)THEN
          IL = (N-1) / (NCOL*NROW) + 1
          IJ = N - (IL-1)*NCOL*NROW
          IR = (IJ-1)/NCOL + 1
          IC = IJ - (IR-1)*NCOL
           WRITE(IOUT,62) IL,IR,IC,Q
   62    FORMAT(1X,'   LAYER ',I5,'   ROW ',I6,'   COL ',I6,
     1       '   FLUX ',1PG15.6)
        ELSE
           WRITE(IOUT,63) N,Q
   63    FORMAT(1X,'    NODE ',I8,'   FLUX ',1PG15.6)
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
   99 CONTINUE
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
        IF(INCLN.GT.0)THEN
           CALL UBUDSVU(KSTP,KPER,TEXT(2),IBCTCB,BUFF(NODES+1:NEQS),
     1                 NCLNNDS,IOUT,PERTIM,TOTIM)
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

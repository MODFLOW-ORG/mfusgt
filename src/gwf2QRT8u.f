      MODULE GWFQRTMODULE
        INTEGER,SAVE,POINTER :: NQRTCL,MXQRT,NQRTVL,NQRTNP,IQRTCB
        INTEGER, SAVE, POINTER ::      MXRTCELLS, IQRTQV,IAFRQRT
        INTEGER,SAVE,POINTER   :: NPQRT,IQRTPB,IQRTFL,NOPRQT,IQCHANGEC
        INTEGER, SAVE, DIMENSION (:), ALLOCATABLE :: IQCHNGTYP        
        INTEGER,      SAVE, DIMENSION(:), ALLOCATABLE :: NodQRT
        REAL,         SAVE, DIMENSION(:,:), ALLOCATABLE ::QRTF
        REAL,         SAVE, DIMENSION(:), ALLOCATABLE ::QRTFLOW,RTAREA
        CHARACTER*16, SAVE, DIMENSION(:),   ALLOCATABLE ::QRTAUX
        REAL,    SAVE, ALLOCATABLE, DIMENSION(:,:)  ::BDQV        
        DOUBLE PRECISION,    SAVE, ALLOCATABLE, DIMENSION(:)   ::BDQTIM 
        INTEGER, SAVE, DIMENSION(:), ALLOCATABLE :: IQRTN 
        INTEGER, SAVE, POINTER :: NBDQTIM,ISTEPQ 
        DOUBLE PRECISION, SAVE, POINTER :: TIMQRT        
      END MODULE GWFQRTMODULE



      SUBROUTINE GWF2QRT8U1AR(IN)
C     ******************************************************************
C     ALLOCATE ARRAY STORAGE AND READ PARAMETERS FOR SINKS AND
C     RETURN FLOWS
C     ******************************************************************
C     SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL, ONLY:IOUT,NCOL,NROW,NLAY,IFREFM,NODES,IUNSTR,NEQS,
     1  ITRNSP
      USE GWFQRTMODULE, ONLY:NQRTCL,MXQRT,NQRTVL,NQRTNP,IQRTCB,NPQRT,
     1  IQRTPB,IQRTFL,NOPRQT,QRTF,QRTAUX,MXRTCELLS,IQRTQV,IAFRQRT,
     2  NodQRT,QRTFLOW,RTAREA,IQCHANGEC,IQCHNGTYP,NBDQTIM,BDQTIM,BDQV,
     3  ISTEPQ,IQRTN,TIMQRT 
      CHARACTER*400 LINE
C     ------------------------------------------------------------------
      ALLOCATE(NQRTCL,MXQRT,NQRTVL,NQRTNP,IQRTCB,MXRTCELLS,IQRTQV,
     1  IAFRQRT)
      ALLOCATE(NPQRT,IQRTPB,IQRTFL,NOPRQT,IQCHANGEC,NBDQTIM,ISTEPQ)
C
C1------IDENTIFY PACKAGE AND INITIALIZE NQRTCL.
      WRITE(IOUT,1)IN
    1 FORMAT(1X,/
     &1X,'QRT8 -- SINK RETURN PACKAGE, VERSION 8, 8/31/2016',/,
     &' INPUT READ FROM UNIT ',I4)
      NQRTCL=0
      NQRTNP=0
      IQRTFL=0
      IQRTQV=0 ! autoflow reduce flag (0=no; 1=yes)
      NBDQTIM = 0
      ISTEPQ = 0
C
C2------READ MAXIMUM NUMBER OF SINKS AND UNIT OR FLAG FOR
C2------CELL-BY-CELL FLOW TERMS.
C     READ COMMENTS (ITEM 0)
      CALL URDCOM(IN,IOUT,LINE)
C     READ ITEM 1
      IF (IFREFM.EQ.0) THEN
        READ(LINE,'(5I10)') MXAQRT,MXRTCELLS,IQRTCB,NPQRT,MXL
        LLOC=51
      ELSE
        LLOC=1
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,MXAQRT,R,IOUT,IN)
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,MXRTCELLS,R,IOUT,IN)
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IQRTCB,R,IOUT,IN)
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NPQRT,R,IOUT,IN)
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,MXL,R,IOUT,IN)
      ENDIF
      WRITE(IOUT,3) MXAQRT
    3 FORMAT(1X,'MAXIMUM OF ',I6,
     &' ACTIVE SINKS WITH RETURN FLOW AT ONE TIME (MXAQRT)')
      WRITE(IOUT,4) MXRTCELLS
    4 FORMAT(1X,'MAXIMUM OF ',I6,
     &' RETURN FLOW CELLS IN SIMULATION (MXRTCELLS)')
      IF (IQRTCB.LT.0) WRITE(IOUT,7)
    7 FORMAT(1X,'CELL-BY-CELL FLOWS WILL BE PRINTED WHEN ICBCFL NOT 0')
      IF (IQRTCB.GT.0) WRITE(IOUT,8) IQRTCB
    8 FORMAT(1X,'CELL-BY-CELL FLOWS WILL BE SAVED ON UNIT ',I4)
C
      IF (NPQRT.GT.0) THEN
        WRITE(IOUT,9) NPQRT,MXL
    9   FORMAT(1X,I5,' Named Parameters     ',I5,' List entries')
      ELSE
        WRITE(IOUT,'(A)') ' No named parameters'
      END IF

C3------READ AUXILIARY VARIABLES AND CBC ALLOCATION OPTION.
      ALLOCATE (QRTAUX (20))
      NAUX=0
      NOPRQT=0
      IQCHANGEC=0      
   10 CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,N,R,IOUT,IN)
      IF(LINE(ISTART:ISTOP).EQ.'AUXILIARY' .OR.
     &       LINE(ISTART:ISTOP).EQ.'AUX') THEN
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,N,R,IOUT,IN)
        IF (NAUX.LT.20) THEN
          NAUX=NAUX+1
          QRTAUX(NAUX)=LINE(ISTART:ISTOP)
          WRITE(IOUT,12) QRTAUX(NAUX)
   12     FORMAT(1X,'AUXILIARY SINK-RETURN FLOW VARIABLE: ',A)
        ENDIF
        GOTO 10
      ELSEIF (LINE(ISTART:ISTOP).EQ.'RETURNFLOW') THEN
        IQRTFL=2                               !************************************
        WRITE(IOUT,13)
        GOTO 10
   13   FORMAT(1X,'RETURN FLOW OPTION IS SELECTED')
      ELSE IF(LINE(ISTART:ISTOP).EQ.'AUTOFLOWREDUCE') THEN
         WRITE(IOUT,16)
   16    FORMAT(1X,'EXTRACTION FLUX WILL BE REDUCED WHEN SATURATED ',
     1       'THICKNESS IS LESS THAN 1 PERCENT OF CELL THICKNESS')
         IQRTQV = 1
         GO TO 10
      ELSE IF(LINE(ISTART:ISTOP).EQ.'IUNIT_AFR_QRT') THEN
         CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IAFR,R,IOUT,IN)
         WRITE(IOUT,25) IAFRQRT
   25    FORMAT(1X,'EXTRACTION REDUCTION INFO WILL BE WRITTEN TO UNIT: ',
     1       I5)
         GO TO 10         
      ELSE IF(LINE(ISTART:ISTOP).EQ.'NOPRINT') THEN
         WRITE(IOUT,14)
   14    FORMAT(1X,'LISTS OF SINK-RETURN CELLS WILL NOT BE PRINTED')
         NOPRQT = 1
         GO TO 10
      ELSEIF (LINE(ISTART:ISTOP).EQ.'CHANGEC') THEN
        IF(ITRNSP.NE.0) IQCHANGEC=1 
        ALLOCATE (IQCHNGTYP(MXAQRT))
        WRITE(IOUT,15)
        GOTO 10
   15   FORMAT(1X,'RETURN FLOW OCCURS WITH CHANGED CONC OR TEMP', 
     *   1X,'(FOR TRANSPORT)')     
      ELSEIF(LINE(ISTART:ISTOP).EQ.'TRANSIENTQ') THEN
C3C------IS KEWORD OPTION FOR TRANSIENT Q LIKE IN qrt PACKAGE
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NBDQTIM,R,IOUT,INOC)
        WRITE(IOUT,18) ABS(NBDQTIM)
        IF(NBDQTIM.LT.0) THEN
          ISTEPQ = 1
          NBDQTIM = -NBDQTIM
          WRITE (IOUT,17)
17        FORMAT(1X,'STAIRCASING WILL BE USED INSTEAD OF INTERPOLATION',
     1      '  BETWEEN DATA POINTS - OPTION STAIRCASE IS ON')
          ALLOCATE(TIMQRT)
          TIMQRT=0          
        ENDIF  
18      FORMAT(/1X,'TRANSIENT Q: TOTAL OF',I5,
     &   ' TIMES WILL BE USED TO DEFINE VARIATIONS IN FLOW.')
        ALLOCATE (BDQTIM(NBDQTIM), BDQV(NBDQTIM,MXAQRT))
        ALLOCATE (IQRTN(MXAQRT))
      ENDIF
C3A-----THERE ARE SIX INPUT VALUES(L, R, C, Q, NUMRT,Rfprop)
C3A-----PLUS ONE LOCATION FOR CELL-BY-CELL FLOW.
      NQRTVL=5+NAUX+IQRTFL                   !************************************
C
C4------ALLOCATE SPACE FOR THE QRTF ARRAY.
      IQRTPB=MXAQRT+1
      MXQRT=MXAQRT+MXL
      ALLOCATE (QRTF(NQRTVL,MXQRT))
      IF(IQRTFL. GT.0) THEN
        ALLOCATE (NodQRT(MXRTCELLS))
        ALLOCATE (QRTFLOW(MXRTCELLS))
        ALLOCATE (RTAREA(MXQRT))
      ENDIF
C ----INITIALIZE  
      DO J=1,MXQRT
        DO N=1,NQRTVL
           QRTF(N,J) = 0.0
        ENDDO
        RTAREA(J) = 0.0
      ENDDO
      DO N=1,MXRTCELLS
        QRTFLOW(N) = 0.0
        NodQRT(N) = 0
      ENDDO
C
C5------READ NAMED PARAMETERS.
      WRITE(IOUT,500) NPQRT
  500 FORMAT(1X,//1X,I5,' Sink-return parameters')
      IF (NPQRT.GT.0) THEN
        NAUX=NQRTVL-5-IQRTFL                 !************************************ 
        LSTSUM=IQRTPB
        ITERPU = 1
        IF (NOPRQT .EQ.1) ITERPU = 99
        DO 100 K=1,NPQRT
          LSTBEG=LSTSUM
C5A-----READ ITEM 2
          CALL UPARLSTRP(LSTSUM,MXQRT,IN,IOUT,IP,'QRT','QRT',1,
     &                   NUMINST)
          NLST=LSTSUM-LSTBEG
          IF (NUMINST.GT.1) NLST = NLST/NUMINST
C5B-----ASSIGN STARTING INDEX FOR READING INSTANCES
          IF (NUMINST.EQ.0) THEN
            IB=0
          ELSE
            IB=1
          ENDIF
C         READ LIST(S) OF CELLS, PRECEDED BY INSTANCE NAME IF NUMINST>0
          LB=LSTBEG
          DO 50 I=IB,NUMINST
            IF (I.GT.0) THEN
              CALL UINSRP(I,IN,IOUT,IP,1)
            ENDIF
C5C-----READ ITEM 3
            CALL SGWF2QRT8LR(NLST,QRTF,LB,NQRTVL,MXQRT,IN,IOUT,
     &           QRTAUX,5,NAUX,IFREFM,NCOL,NROW,NLAY,ITERPU,IQRTFL,
     &           IUNSTR,NEQS,NodQRT,MXRTCELLS,IQCHANGEC,IQCHNGTYP)
            LB = LB+NLST
   50     CONTINUE
  100   CONTINUE
      ENDIF
C6 -----FOR TRANSIENT Q, READ THE TIME VALUES AND Q VALUES      
      IF(NBDQTIM.NE.0)THEN 
C6A------READ TIMES AT WHICH SPECIFIED FLOW AND HEAD VALUES WILL BE READ
        READ(IN,*) IQRTUN,CNSTM
        WRITE(IOUT,310) IQRTUN,CNSTM
310     FORMAT(1X,'TIMES FOR SPECIFIED-FLOW VALUES WILL BE READ',
     & ' ON UNIT',I4,' AND',/,
     &' MULTIPLIED BY',G12.4,'.')
        READ(IQRTUN,*) (BDQTIM(L),L=1,NBDQTIM)
        DO L=1,NBDQTIM
          BDQTIM(L)=BDQTIM(L)*CNSTM
        ENDDO      
C6A------READ SPECIFIED FLOW VALUES AT EACH Q CELL OF QRT PACKAGE
        READ(IN,*) IQRTUN,CNSTM
        WRITE(IOUT,311) IQRTUN,CNSTM
311     FORMAT(1X,'SPECIFIED-FLOW VALUES WILL BE READ',
     & ' ON UNIT',I4,' AND',/,
     &' MULTIPLIED BY',G12.4,'.')
        DO N = 1,MXAQRT        
          READ(IQRTUN,*) IQRTN(N),(BDQV(L,N),L=1,NBDQTIM)
          DO  L=1,NBDQTIM
            BDQV(L,N)=BDQV(L,N)*CNSTM
          ENDDO        
        ENDDO  
C6B -------------PRINT INTERPOLATION TIMES AND NODAL FLUXES IF REQUESTED
        IF(NOPRQT.EQ.0) THEN 
         WRITE(IOUT,320) NBDQTIM
320      FORMAT(1X,I5,' TIMES FOR SPECIFYING FLOWS:')
         WRITE(IOUT,322) (L,L=1,NBDQTIM)
322      FORMAT(16X,I8,99I12)
         WRITE(IOUT,326) (BDQTIM(L),L=1,NBDQTIM)
326      FORMAT(17X,100G12.4)    
         WRITE(IOUT,327)
327      FORMAT(6X,' NODE NUMBER AND SPECIFIED FLOW VALUES:')
         DO N=1,MXAQRT
           WRITE(IOUT,328) IQRTN(N),(BDQV(L,N),L=1,NBDQTIM) 
328        FORMAT(1X, I8,5X, 100G12.4)          
         ENDDO
        ENDIF
      ENDIF
C
C6------RETURN
      RETURN
      END
      SUBROUTINE GWF2QRT8U1RP(IN)
C     ******************************************************************
C     READ SINK FLUX.  IF THE RETURNFLOW OPTION IS SELECTED,
C     READ NUMBER OF RECIPIENT CELLS AND PROPORTION FOLLOWED
C     BY THE RECIPIENT CELL NUMBERS.
C     ******************************************************************
C
C     SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL, ONLY:IOUT,NCOL,NROW,NLAY,IFREFM,IUNSTR,NODES,NEQS,
     1    AREA,IBOUND,ITRNSP  
      USE GWFQRTMODULE, ONLY:NQRTCL,MXQRT,NQRTVL,NQRTNP,NPQRT,NodQRT,
     1    IQRTPB,IQRTFL,NOPRQT,QRTF,QRTAUX,MXRTCELLS,QRTFLOW,RTAREA,
     1    IQCHANGEC,IQCHNGTYP
      DOUBLE PRECISION FR
C     ------------------------------------------------------------------
C
C1------IDENTIFY PACKAGE AND INITIALIZE NQRTCL.
      WRITE(IOUT,1)IN
    1 FORMAT(1X,/
     &1X,'QRT8 -- SINK RETURN PACKAGE, VERSION 8, 8/31/2016',/,
     &' INPUT READ FROM UNIT ',I4)
C
C1------READ ITMP (NUMBER OF SINKS OR FLAG TO REUSE DATA) AND
C1------NUMBER OF PARAMETERS.
      IF (NPQRT.GT.0) THEN
        IF (IFREFM.EQ.0) THEN
          READ(IN,'(2I10)') ITMP,NP
        ELSE
          READ(IN,*) ITMP,NP
        ENDIF
      ELSE
        NP=0
        IF (IFREFM.EQ.0) THEN
          READ(IN,'(I10)') ITMP
        ELSE
          READ(IN,*) ITMP
        ENDIF
      ENDIF
C
C------CALCULATE SOME CONSTANTS
      NAUX=NQRTVL-5-IQRTFL             !************************************
      ITERPU = 1
      IOUTU = IOUT
      IF (NOPRQT.EQ.1) THEN
        ITERPU = 99
        IOUTU = -IOUT
      ENDIF
C
C2------DETERMINE THE NUMBER OF NON-PARAMETER SINK-RETURN CELLS.
      IF (ITMP.LT.0) THEN
        WRITE(IOUT,7)
    7   FORMAT(1X,/,
     &' REUSING NON-PARAMETER SINK-RETURN CELLS FROM',
     &' LAST STRESS PERIOD')
      ELSE
        NQRTNP=ITMP
      ENDIF
C
C3------IF THERE ARE NEW NON-PARAMETER SINK-RETURN CELLS, READ THEM.
      MXAQRT=IQRTPB-1
      IF (ITMP.GT.0) THEN
        IF (NQRTNP.GT.MXAQRT) THEN
          WRITE(IOUT,500) NQRTNP,MXAQRT
  500     FORMAT(1X,/1X,'THE NUMBER OF ACTIVE QRT SINKS (',I6,
     &           ') IS GREATER THAN MXAQRT(',I6,')')
          CALL USTOP(' ')
        ENDIF
        CALL SGWF2QRT8LR(NQRTNP,QRTF,1,NQRTVL,MXQRT,IN,IOUT,
     &       QRTAUX,5,NAUX,IFREFM,NCOL,NROW,NLAY,ITERPU,IQRTFL,IUNSTR,
     &       NEQS,NodQRT,MXRTCELLS,IQCHANGEC,IQCHNGTYP)
      ENDIF
      NQRTCL=NQRTNP
C
C1C-----IF THERE ARE ACTIVE QRT PARAMETERS, READ THEM AND SUBSTITUTE
      CALL PRESET('QRT')
      IF (NP.GT.0) THEN
        NREAD=NQRTVL -1
        DO 30 N=1,NP
          CALL SGWF2QRT8LS(IN,IOUTU,QRTF,NQRTVL,MXQRT,NREAD,MXAQRT,
     &         NQRTCL,QRTAUX,20,NAUX,IQRTFL,IUNSTR,NodQRT,MXRTCELLS)
   30   CONTINUE
      ENDIF
C --------------------------------------------------------------------------
C     PREPARE ARRAYS: 
C ----------------------------------------------------------------------
      IF (IQRTFL.GT.0) THEN
        IRT = 0 ! SET POINTER FOR RETURN FLOW NODES ARRAY
        DO  L=1,NQRTCL
C -----------COMPUTE RETURN FLOW AREAS FOR EACH SINK IN RTAREA                 
          RTAREA(L) = 0
          NumRT = QRTF(5,L)
          IF (NumRT.EQ.0) CYCLE
          DO JJ = 1,NumRT
            IRT = IRT + 1  
            INR = NodQRT(IRT)
            IF (IBOUND(INR) .GT. 0) THEN
              RTAREA(L) = RTAREA(L) + AREA(INR)  
            END IF
          ENDDO
        ENDDO 
C ----------------------------------------------------------------------        
        IRT = 0 ! SET POINTER FOR RETURN FLOW NODES ARRAY
        DO  L=1,NQRTCL          
C ------------PARTITION SINK ONTO RETURN FLOW NODES IN QRTFLOW   
          NumRT = QRTF(5,L)  
          IF (NumRT.EQ.0) CYCLE
          FR = QRTF(6,L) / RTAREA(L) ! FRACTURN RETURNED DIVIDED BY AREA          
          DO JJ = 1,NumRT
            IRT = IRT + 1  
            INR = NodQRT(IRT)
            IF (IBOUND(INR) .GT. 0) THEN
              QRTFLOW(IRT) = FR * QRTF(4,L)* AREA(INR)
            END IF
          ENDDO
        ENDDO 
      ENDIF
C --------------------------------------------------------------------------
C
C3------PRINT NUMBER OF SINK-RETURN CELLS IN CURRENT STRESS PERIOD.
      WRITE (IOUT,510) NQRTCL
  510 FORMAT(1X,/1X,I6,' SINK-RETURN CELLS')
C8------PRINT INDEX OF QRT-RETURN CELL CONCENTRATION CHANGE TYPE FOR TRANSPORT.
      IF(ITRNSP.NE.0) THEN
        IF(IQCHANGEC.NE.0) THEN
          WRITE (IOUT,511) (IQCHNGTYP(I),I=1,NQRTCL)
  511    FORMAT(/1X,'LIST OF CHANGE TYPE INDICES TO RETURN FLOW SOLUTES'
     &    /20(1X,I5))
        ENDIF
      ENDIF        
C
C8------RETURN.
      RETURN
      END
      SUBROUTINE GWF2QRT8U1FM
C     ******************************************************************
C     ADD SINK-RETURN FLOW TO SOURCE TERMS FOR BOTH SINK-RETURN CELLS
C     AND RECIPIENT CELLS
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,       ONLY: HNEW,AMAT,RHS,IBOUND,IA,TOP,BOT,AREA,NODES
      USE GWFQRTMODULE, ONLY: NQRTCL,QRTF,IQRTFL,NodQRT,IQRTQV,IQRTFL,
     1                        NQRTVL,QRTFLOW
      USE CLN1MODULE, ONLY: ACLNNDS
C
      DOUBLE PRECISION SinkQ,FR,HD,RFPROP
      DOUBLE PRECISION THCK,X,Y,Q,QTHIK,BOTT,QA,EPS,QEPS,DQ
C     ------------------------------------------------------------------
C
C1------IF NQRTCL<=0 THERE ARE NO SINKS. RETURN.
      IF (NQRTCL.LE.0) RETURN
C
C2------PROCESS EACH CELL IN THE SINK-RETURN CELL LIST.
      IRT = 0                  !------------------! POINTER FOR LOCATION IN NodQRT ARRAY
      DO 100 L=1,NQRTCL
C
C3------GET NODE NUMBER OF CELL CONTAINING SINK.
        ND=QRTF(1,L)
C
C4-------IF THE CELL IS EXTERNAL SKIP IT.
        IF (IBOUND(ND).LE.0) GOTO 100
C
C5-------IF THE CELL IS INTERNAL GET THE SINK DATA.
        Q = -QRTF(4,L)
C
C----------------------------------------------------------------------------
C6------SUBTRACT Q FROM RHS BUT WITH AUTOFLOWREDUCE AS NEEDED
        IF(IQRTQV.EQ.1)THEN
          IPIV = IA(ND)
C---------HONOR SUPPLY/DEMAND CONDITIONS FOR EXTRACTION WELLS (NEWTON METHOD)
          HD = HNEW(ND)
          IF(ND.GT.NODES)THEN
            ICLN = ND-NODES
            CALL CLNV(ICLN,THCK)
            BOTT = ACLNNDS(ICLN,5)
          ELSE
            THCK = (TOP(ND) - BOT(ND)) 
            BOTT = BOT(ND)
          ENDIF
          QTHIK = THCK * 0.01  ! SMOOTH OVER 1 PERCENT OF CELL THICKNESS OR 1 FT WHICHEVER IS LESS
          QTHIK = MIN (QTHIK,1.0)
          X = (HD - BOTT) /QTHIK
          CALL SMOOTH(X,Y)
          QA = Q * Y
C---------CALCULATE DQ/DH
          EPS = 0.01 * QTHIK
          X = (HD+EPS - BOTT) /QTHIK
          CALL SMOOTH(X,Y)
          QEPS = Q*Y
          DQ = (QEPS - QA) / EPS
          AMAT(IPIV) = AMAT(IPIV) + DQ
          RHS(ND) = RHS(ND) - QA + DQ*HD
        ELSE
          RHS(ND) = RHS(ND) - Q
        ENDIF
C----------------------------------------------------------------------------
C7------INCLUDE RETURN FLOW
        IF (IQRTFL.GT.0) THEN
          NumRT = QRTF(5,L)
          IF (NumRT.EQ.0) CYCLE
          DO JJ = 1,NumRT
            IRT = IRT + 1  
            INR = NodQRT(IRT)
            IF (IBOUND(INR) .GT. 0) THEN
              RHS(INR) = RHS(INR) - QRTFLOW(IRT)
            END IF
          ENDDO
        ENDIF
  100 CONTINUE
C
C8------RETURN.
      RETURN
      END
C----------------------------------------------------------------
      SUBROUTINE GWF2QRT8U1BD(KSTP,KPER)
C     ******************************************************************
C     CALCULATE VOLUMETRIC BUDGET FOR SINK-RETURN CELLS
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,  ONLY: IOUT,HNEW,IBOUND,BUFF,NCOL,NROW,NLAY,
     *                   NODES,NEQS,IUNSTR,TOP,BOT,AREA,FMBE,INCLN
      USE GWFBASMODULE, ONLY: MSUM,VBNM,VBVL,PERTIM,TOTIM,DELT,ICBCFL,
     1                        IAUXSV
      USE GWFQRTMODULE, ONLY: QRTF,NQRTCL,MXQRT,IQRTCB,NQRTVL,IQRTFL,
     1          QRTAUX,NodQRT,QRTFLOW,IQRTQV,RTAREA,MXRTCELLS,IAFRQRT
      USE CLN1MODULE, ONLY: ACLNNDS,NCLNNDS,ICLNCB
C
      DOUBLE PRECISION HD,RATIN,RATOUT,QQ,QIN,QTHIK,X,Y,THCK,BOTT,FR
      CHARACTER*16 TEXT
      DATA TEXT /'    SINKS (QRT)'/
C     ------------------------------------------------------------------
C
C1------INITIALIZE CELL-BY-CELL FLOW TERM FLAG (IBD) AND
C1------ACCUMULATORS (RATIN AND RATOUT).
      ZERO=0.
      QIN = ZERO
      RATIN=ZERO
      RATOUT=ZERO
      NWELLAFR=0      
      IBD=0
      IF (IQRTCB.LT.0 .AND. ICBCFL.NE.0) IBD=-1
      IF (IQRTCB.GT.0) IBD=ICBCFL
      IBDLBL=0
C
C2------IF CELL-BY-CELL FLOWS WILL BE SAVED AS A LIST, WRITE HEADER.
      IF (IBD.EQ.2) THEN
        NAUX=NQRTVL-5-IQRTFL                 !************************************ 
        IF (IAUXSV.EQ.0) NAUX = 0
         IICLNCB=0
         NNCLNNDS=0
         IF(INCLN.GT.0) THEN
           IICLNCB=ICLNCB
           NNCLNNDS=NCLNNDS
         ENDIF
         CALL UBDSVHDRQRT(IUNSTR,KSTP,KPER,IOUT,IQRTCB,IICLNCB,NODES,
     1    NNCLNNDS,NCOL,NROW,NLAY,NQRTCL,NQRTVL,NAUX,IBOUND,
     2    TEXT,QRTAUX,DELT,PERTIM,TOTIM,IQRTFL,NodQRT,MXRTCELLS,QRTF)
C        
C         IF(IUNSTR.EQ.0)THEN
C        CALL UBDSV4(KSTP,KPER,TEXT,NAUX,QRTAUX,IQRTCB,NCOL,NROW,NLAY,
C     &              NQRTCL,IOUT,DELT,PERTIM,TOTIM,IBOUND)
C         ELSE
C           CALL UBDSV4U(KSTP,KPER,TEXT,NAUX,QRTAUX,IQRTCB,NEQS,
C     1          NQRTCL,IOUT,DELT,PERTIM,TOTIM,IBOUND)
C         ENDIF
      ENDIF
C
C3------CLEAR THE BUFFER.
      DO 10 N=1,NEQS
            BUFF(N)=ZERO
   10     CONTINUE
C
C4------IF THERE ARE NO SINK-RETURN CELLS THEN DO NOT ACCUMULATE FLOW.
      IF (NQRTCL.LE.0) GOTO 200
C
C5------LOOP THROUGH EACH SINK-RETURN CELL, CALCULATING FLOW.
      IRT = 0   !------------------! INDEX FOR LOCATION IN NodQRT ARRAY
      DO 100 L=1,NQRTCL
C
C5A-----GET NODE NUMBER OF CELL CONTAINING SINK.
        ND=QRTF(1,L)
        QQ=ZERO
        Q = ZERO 
C
C5B-----IF CELL IS NO-FLOW OR CONSTANT-HEAD, IGNORE IT.
        IF (IBOUND(ND).LE.0) GOTO 99
C
C5C-----GET SINK PARAMETERS FROM SINK-RETURN LIST.
        QQ = -QRTF(4,L)
C-----------------------------------------------------------------------
C
        IF(IQRTQV.EQ.1)THEN
C---------HONOR SUPPLY/DEMAND CONDITIONS FOR EXTRACTION WELLS
          HD = HNEW(ND)
          IF(ND.GT.NODES)THEN
            ICLN = ND-NODES
            CALL CLNV(ICLN,THCK)
            BOTT = ACLNNDS(ICLN,5)
          ELSE
            THCK = (TOP(ND) - BOT(ND)) 
            BOTT = BOT(ND)
          ENDIF
          QTHIK = THCK * 0.01  ! SMOOTH OVER 1 PERCENT OF CELL THICKNESS OR 1 FT WHICHEVER IS LESS
          QTHIK = MIN (QTHIK,1.0)
          X = (HD - BOTT) /QTHIK
          CALL SMOOTH(X,Y)
          QQ = QQ * Y
        ENDIF
        Q=QQ
        RATOUT=RATOUT-QQ
C5F-----ADD Q TO BUFFER.
        BUFF(ND) = BUFF(ND) + QQ
        FMBE(ND) = FMBE(ND) + QQ        
        QRTF(NQRTVL,L) = -QQ
C----------------------------------------------------------------------------
C7------INCLUDE RETURN FLOW. ** note the iteration lag on qrtflow used in fm is not done here and latest qq is used
        IF (IQRTFL.GT.0) THEN
          NumRT = QRTF(5,L)
          RFPROP = QRTF(6,L)
          IF (NumRT.EQ.0) CYCLE 
          FR = QRTF(6,L) / RTAREA(L) ! FRACTURN RETURNED DIVIDED BY AREA           
          DO JJ = 1,NumRT
            IRT = IRT + 1  
            INR = NodQRT(IRT)
            IF (IBOUND(INR) .GT. 0) THEN
              QIN = QRTFLOW(IRT)
              BUFF(INR) = BUFF(INR) + QIN
              FMBE(INR) = FMBE(INR) + QIN    
              RATIN = RATIN + QIN
            END IF
C-----------COMPUTE RETURN FLOW FOR NEXT TIME STEP (QQ MAY BE REDUCED BY AUTOFLOWREDUCE)            
c............comment for not having return flow affected by autoflowreduce
c            QRTFLOW(IRT) = - FR * QQ * AREA(INR)
          ENDDO
        ENDIF
C
C5E-----PRINT THE INDIVIDUAL RATES IF REQUESTED(IQRTCB<0).
        IF (IBD.LT.0) THEN
          IF (IBDLBL.EQ.0) WRITE(IOUT,61) TEXT,KPER,KSTP
   61     FORMAT(1X,/1X,A,'   PERIOD ',I4,'   STEP ',I3)
          WRITE(IOUT,62) L,IL,IR,IC,Q
   62     FORMAT(1X,'SINK ',I6,'   LAYER ',I3,'   ROW ',I5,
     &       '   COL ',I5,'   RATE ',1PG15.6)
          IF (NumRT.NE.0) THEN
            DO I=IRT-NUMRT+1, IRT
              WRITE(IOUT,550) L,NodQRT(I),QRTFLOW(I)
            ENDDO
  550       FORMAT(1X,'SINK ',I6,
     *        ' RETURN:  NODE ',I10,'   RATE ',1PG15.6)
          ENDIF
          IBDLBL=1
        ENDIF
C
C5G-----IF SAVING CELL-BY-CELL FLOWS IN A LIST, WRITE FLOW.
   99 CONTINUE 
      IF(IBD.EQ.2)THEN
        CALL UBDSVREC(IUNSTR,ND,NODES,NNCLNNDS,IQRTCB,IICLNCB,NQRTVL,
     1    6,NAUX,Q,QRTF(:,L),IBOUND,NCOL,NROW,NLAY)
        IF (NumRT.NE.0) THEN
          DO I=IRT-NUMRT+1, IRT
            CALL UBDSVREC(IUNSTR,I,NODES,NNCLNNDS,IQRTCB,IICLNCB,NQRTVL,
     1    6,NAUX,QRTFLOW(I),QRTF(:,L),IBOUND,NCOL,NROW,NLAY)      
          ENDDO    
        ENDIF      
      ENDIF
C
C5J-----WRITE FLOW REDUCTION INFO IF REQUESTED
      IF (IQRTQV.GT.0. .AND. IAFRQRT.GT.0) THEN
        IF(Q.GT.QRTF(4,L)) THEN
            IF (NWELLAFR.EQ.0) THEN
              WRITE(IAFRQRT, *)
              WRITE(IAFRQRT,300) KPER,KSTP
            END IF
            IF(IUNSTR.EQ.0) THEN
              IF(N.LE.NODES) THEN
                IL = (N-1) / (NCOL*NROW) + 1
                IJ = N - (IL-1)*NCOL*NROW
                IR = (IJ-1)/NCOL + 1
                IC = IJ - (IR-1)*NCOL
                IF(NWELLAFR.EQ.0) WRITE(IAFRQRT,400)
                WRITE(IAFRQRT,500) L,IL,IR,IC,QRTF(4,L),Q,HD,BOTT
              ELSE
                IL = N-NODES
                IR = 0
                IC = 0
                IF(NWELLAFR.EQ.0) WRITE(IAFRQRT,401)
                WRITE(IAFRQRT,501) L,IL,QRTF(4,L),Q,HD,BOTT
              ENDIF  
            ELSE
              IF(NWELLAFR.EQ.0) WRITE(IAFRQRT,402)
              WRITE(IAFRQRT,501) L,N,QRTF(4,L),Q,HD,BOTT
            ENDIF
            NWELLAFR = NWELLAFR + 1
        ENDIF
      END IF
  300 FORMAT(' WELLS WITH REDUCED PUMPING FOR STRESS PERIOD ',I5,
     1      ' TIME STEP ',I5)
  400 FORMAT('WELL.NO   LAY   ROW   COL         APPL.Q          ACT.Q',
     1       '        GW_HEAD       CELL_BOT')
  401 FORMAT('WELL.NO  CLN NODE     APPL.Q          ACT.Q',
     1       '        GW_HEAD       CELL_BOT')
  402 FORMAT('WELL.NO  NODE         APPL.Q          ACT.Q',
     1       '        GW_HEAD       CELL_BOT')      
  500 FORMAT(I7,3I6,4(1PG15.6))
  501 FORMAT(I7,1X,I9,4(1PG15.6))      
C
C      IF (IBD.EQ.2) THEN
C          IF(IUNSTR.EQ.0)THEN
C            IL = (ND-1) / (NCOL*NROW) + 1
C            IJ = ND - (IL-1)*NCOL*NROW
C            IR = (IJ-1)/NCOL + 1
C            IC = IJ - (IR-1)*NCOL
C            CALL UBDSVB(IQRTCB,NCOL,NROW,IC,IR,IL,Q,
C     1                  QRTF(1,L),NQRTVL,NAUX,10,IBOUND,NLAY)
C            IF (NumRT.NE.0) THEN
C              DO I=IRT-NUMRT+1, IRT
C                ILR = (ND-1) / (NCOL*NROW) + 1
C                IJR = ND - (ILR-1)*NCOL*NROW
C                IRR = (IJR-1)/NCOL + 1
C                ICR = IJR - (IRR-1)*NCOL
C                CALL UBDSVB(IQRTCB,NCOL,NROW,ICR,IRR,ILR,QRTFLOW(I),
C     &             QRTF(1,L),NQRTVL,NAUX,10,IBOUND,NLAY)
C              ENDDO
C            ENDIF
C          ELSE
C            CALL UBDSVBU(IQRTCB,NEQS,ND,Q,
C     1                  QRTF(1,L),NQRTVL,NAUX,10,IBOUND)
C            IF (NumRT.NE.0) THEN
C              DO I=IRT-NUMRT+1, IRT
C                CALL UBDSVBU(IQRTCB,NEQS,NodQRT(I),QRTFLOW(I),QRTF(1,L),
C     &                  NQRTVL,NAUX,10,IBOUND)
C              ENDDO    
C            ENDIF
C          ENDIF
C        ENDIF    
  100 CONTINUE
C
C6------IF CELL-BY-CELL FLOW WILL BE SAVED AS A 3-D ARRAY,
C6------CALL UBUDSV TO SAVE THEM.
      IF(IUNSTR.EQ.0)THEN
      IF (IBD.EQ.1) CALL UBUDSV(KSTP,KPER,TEXT,IQRTCB,BUFF,NCOL,NROW,
     &                          NLAY,IOUT)
      ELSE
C      IF (IBD.EQ.1) CALL UBUDSVU(KSTP,KPER,TEXT,IQRTCB,BUFF,NEQS,IOUT,
C     *  PERTIM,TOTIM)
      IF (IBD.EQ.1) CALL UBUDSVU(KSTP,KPER,TEXT,IQRTCB,BUFF,NODES,IOUT,
     *  PERTIM,TOTIM)
      ENDIF
C
      IF(IBD.EQ.1.AND.INCLN.GT.0)THEN
        IF(ICLNCB.GT.0) CALL UBUDSVU(KSTP,KPER,TEXT,ICLNCB,
     1    BUFF(NODES+1),NCLNNDS,IOUT,PERTIM,TOTIM)
      ENDIF
C
C7------MOVE RATES,VOLUMES & LABELS INTO ARRAYS FOR PRINTING.
  200 CONTINUE
      RIN = RATIN
      ROUT=RATOUT
      VBVL(3,MSUM)=RIN
      VBVL(1,MSUM)=VBVL(1,MSUM)+RIN*DELT
      VBVL(4,MSUM)=ROUT
      VBVL(2,MSUM)=VBVL(2,MSUM)+ROUT*DELT
      VBNM(MSUM)=TEXT
C
C8------INCREMENT BUDGET TERM COUNTER.
      MSUM=MSUM+1
C
C9------RETURN.
      RETURN
      END
C----------------------------------------------------------------------
      SUBROUTINE SGWF2QRT8LR(NLIST,QRTF,LSTBEG,NQRTVL,MXQRT,
     &        INPACK,IOUT,QRTAUX,NCAUX,NAUX,IFREFM,NCOL,NROW,NLAY,ITERP,
     &        IQRTFL,IUNSTR,NEQS,NodQRT,MXRTCELLS,IQCHANGEC,IQCHNGTYP)
C     ******************************************************************
C     Read and print a list of drain and optional associated
C     return-flow recipient cells.  NAUX of the values in the list are
C     optional -- auxiliary data.
C     ******************************************************************
      CHARACTER*57 LABEL1, LABEL2, LABEL3,LABEL4,LABEL5
      CHARACTER*16 QRTAUX(NCAUX)
      DIMENSION QRTF(NQRTVL,MXQRT),NodQRT(MXRTCELLS),IQCHNGTYP(MXQRT)
      CHARACTER*400 LINE,FNAME
C
      CHARACTER*24 ANAME
C
      DATA ANAME /'       RETURN FLOW NODES'/      
      DATA NUNOPN/99/
C     ------------------------------------------------------------------
C
      IERR = 0
      ISCLOC1 = 4
      ISCLOC2 = 4
      IN=INPACK
      ICLOSE=0
      LABEL1='SINK NO.  LAYER   ROW   COL     STRESS FACTOR'
      LABEL2='          ----SINK CELL----                        RETURN'
      LABEL3='SINK NO.  LAY     ROW    COL   NO.RECIPIENT_CELLS   PROP.'
      LABEL4='SINK NO.       NODE             STRESS FACTOR'
      LABEL5='SINK NO.       NODE       NO. OF RECIPIENT CELLS   PROP.'
C
C  Check for and decode EXTERNAL and SFAC records.
      READ(IN,'(A)') LINE
      SFAC=1.
      LLOC=1
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,I,R,IOUT,IN)
      IF (LINE(ISTART:ISTOP).EQ.'EXTERNAL') THEN
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,I,R,IOUT,IN)
        IN=I
        WRITE(IOUT,510) IN
  510   FORMAT(1X,'Reading list on unit ',I4)
        READ(IN,'(A)') LINE
      ELSEIF (LINE(ISTART:ISTOP).EQ.'OPEN/CLOSE') THEN
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,0,N,R,IOUT,IN)
        FNAME=LINE(ISTART:ISTOP)
        IN=NUNOPN
        WRITE(IOUT,520) IN,FNAME
  520   FORMAT(1X,/1X,'OPENING FILE ON UNIT ',I4,':',/1X,A)
        OPEN(UNIT=IN,FILE=FNAME)
        ICLOSE=1
        READ(IN,'(A)') LINE
      ENDIF
      LLOC=1
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,I,R,IOUT,IN)
      IF (LINE(ISTART:ISTOP).EQ.'SFAC') THEN
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,I,SFAC,IOUT,IN)
        WRITE(IOUT,530) SFAC
  530   FORMAT(1X,'LIST SCALING FACTOR= ',1PG12.5)
        IF (ISCLOC1.EQ.ISCLOC2) THEN
          WRITE(IOUT,540) ISCLOC1
  540     FORMAT(1X,'(THE SCALE FACTOR WAS APPLIED TO FIELD ',I2,')')
        ELSE
          WRITE(IOUT,550) ISCLOC1,ISCLOC2
  550     FORMAT(1X,'(THE SCALE FACTOR WAS APPLIED TO FIELDS ',
     &           I2,'-',I2,')')
        ENDIF
        READ(IN,'(A)') LINE
      ENDIF
C
C  Write a label for the list.
      WRITE(IOUT,'(1X)')
      IF(IUNSTR.EQ.0)THEN
        CALL ULSTLB(IOUT,LABEL1,QRTAUX,NCAUX,NAUX)
      ELSE
        CALL ULSTLB(IOUT,LABEL4,QRTAUX,NCAUX,NAUX)
      ENDIF
C------READ LIST OF SINK-RETURN NODES FOR STRUCTURED AND UNSTRUCTURED GRIDS
C----------------------------------------------------------------------------
      IF (IUNSTR.EQ.0) THEN
C
C  Read the list
      NREAD2=NQRTVL-1
      NREAD1=NREAD2-NAUX
      N=NLIST+LSTBEG-1
      DO 100 II=LSTBEG,N
C  Read a line into the buffer.  (The first line has already been read
C  in order to scan for EXTERNAL and SFAC records.)
        IF (II.NE.LSTBEG) READ(IN,'(A)') LINE
C
C  Read the non-optional values from the line.
        IF (IQRTFL.EQ.0) THEN
          IF (IFREFM.EQ.0) THEN
            READ(LINE,'(3I10,9F10.0)') K,I,J,(QRTF(JJ,II),JJ=4,NREAD1)
            LLOC=10*NREAD1+1
          ELSE
            LLOC=1
            CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,K,R,IOUT,IN)
            CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,I,R,IOUT,IN)
            CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,J,R,IOUT,IN)
            DO 10 JJ=4,NREAD1
              CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,IDUM,QRTF(JJ,II),
     &                    IOUT,IN)
   10       CONTINUE
          ENDIF
        ELSE
          IF (IFREFM.EQ.0) THEN
            READ(LINE,'(3I10,2F10.0,3I10,9F10.0)') K,I,J,
     &          QRTF(4,II),NumRT,(QRTF(JJ,II),JJ=6,NREAD1)
            LLOC=10*NREAD1+1
            IF(IQCHANGEC.NE.0) READ(LINE,'(I10)')IQCHNGTYP(II)        
          ELSE
            LLOC=1
            CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,K,R,IOUT,IN)
            CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,I,R,IOUT,IN)
            CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,J,R,IOUT,IN)
            CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,IDUM,QRTF(4,II),IOUT,
     &                  IN)
            CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NumRT,R,IOUT,IN)
            IF (NumRT.EQ.0 .AND. NREAD1.EQ.6 .AND. NAUX.EQ.0) THEN
              QRTF(5,II) = 0.0
            ELSE
              DO 20 JJ=6,NREAD1
                CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,IDUM,QRTF(JJ,II),
     &                      IOUT,IN)
   20         CONTINUE
            ENDIF
            IF(IQCHANGEC.NE.0) 
     &     CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IQCHNGTYP(II),R,IOUT,IN)            
          ENDIF
          QRTF(5,II) = NumRT
        ENDIF
        QRTF(1,II)=K
        QRTF(2,II)=I
        QRTF(3,II)=J
        DO 50 ILOC=ISCLOC1,ISCLOC2
          QRTF(ILOC,II)=QRTF(ILOC,II)*SFAC
   50   CONTINUE
C
C  Read the optional values from the line
        IF (NAUX.GT.0) THEN
          DO 60 JJ=NREAD1+1,NREAD2
            CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,IDUM,QRTF(JJ,II),IOUT,
     &                  IN)
   60     CONTINUE
        ENDIF
C
C  Write the values that were read and that are not related to
C  return flow.
        NN=II-LSTBEG+1
        IF (IQRTFL.EQ.0) THEN
          WRITE(IOUT,570) NN,K,I,J,(QRTF(JJ,II),JJ=4,NREAD2)
  570     FORMAT(1X,I6,I7,I7,I7,15X,14G16.4)
        ELSE
          IF (NREAD2.GE.6) THEN
            WRITE(IOUT,570) NN,K,I,J,QRTF(4,II),
     &                      (QRTF(JJ,II),JJ=10,NREAD2)
          ELSE
            WRITE(IOUT,570) NN,K,I,J,QRTF(4,II)
          ENDIF
        ENDIF
C
C  Check for illegal grid location
        IF (K.LT.1 .OR. K.GT.NLAY) THEN
          WRITE(IOUT,*) ' ERROR: Layer number is outside of the grid'
          IERR = 1
        ENDIF
        IF (I.LT.1 .OR. I.GT.NROW) THEN
          WRITE(IOUT,*) ' ERROR: Row number is outside of the grid'
          IERR = 1
        ENDIF
        IF (J.LT.1 .OR. J.GT.NCOL) THEN
          WRITE(IOUT,*) ' ERROR: Column number is outside of the grid'
          IERR = 1
        ENDIF
        IF (IERR.NE.0) CALL USTOP(' ')
  100 CONTINUE
C
C     Check and write data related to return-flow recipient cells
      IF (IQRTFL.GT.0) THEN
        WRITE(IOUT,'(/,1X,A,/,1X,A)') LABEL2,LABEL3
        NN = 0
        DO 110 II=LSTBEG,N
          NN = NN + 1
          K = QRTF(1,II)
          I = QRTF(2,II)
          J = QRTF(3,II)
          NumRT = QRTF(5,II)
          RFP = QRTF(6,II)
          WRITE(IOUT,600) NN,K,I,J,NumRT,RFP
  600     FORMAT(1X,I6,4I7,2X,F8.6)
C
C  Check for incorrect stuff
          IF (NumRT.NE.0) THEN
C
C  Check for invalid return-flow proportion
            IF (RFP.LT.0.0 .OR. RFP.GT.1.0) THEN
              WRITE(IOUT,590)
  590         FORMAT(' ERROR: Proportion must be between 0.0 and 1.0')
              IERR = 1
            ENDIF
          ENDIF
C  Change lay, row, col to global node numbers
          NODNO = NROW*NCOL*(K-1) + (I-1)*NCOL + J
          QRTF(1,II) = NODNO
          QRTF(NQRTVL,II) = QRTF(4,II)
C
          IF (IERR.NE.0) CALL USTOP(' ')
  110   CONTINUE
      ENDIF
C------READ LIST OF SINK-RETURN NODES FOR UNSTRUCTURED GRIDS
C----------------------------------------------------------------------------
      ELSE
C----------------------------------------------------------------------------
C
C  Read the list
      NREAD2=NQRTVL-1
      NREAD1=NREAD2-NAUX
      N=NLIST+LSTBEG-1
      DO 1001 II=LSTBEG,N
C  Read a line into the buffer.  (The first line has already been read
C  in order to scan for EXTERNAL and SFAC records.)
        IF (II.NE.LSTBEG) READ(IN,'(A)') LINE
C
C  Read the non-optional values from the line.
        IF (IQRTFL.EQ.0) THEN
          IF (IFREFM.EQ.0) THEN
            READ(LINE,'(3I10,9F10.0)') ND,(QRTF(JJ,II),JJ=4,NREAD1)
            LLOC=10*NREAD1+1
          ELSE
            LLOC=1
            CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,ND,R,IOUT,IN)
            DO 101 JJ=4,NREAD1
              CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,IDUM,QRTF(JJ,II),
     &                    IOUT,IN)
  101       CONTINUE
          ENDIF
        ELSE
          IF (IFREFM.EQ.0) THEN
            READ(LINE,'(I10,2F10.0,I10,9F10.0)') ND,
     &          QRTF(4,II),NumRT,(QRTF(JJ,II),JJ=9,NREAD1)
            LLOC=10*NREAD1+1
            IF(IQCHANGEC.NE.0) READ(LINE,'(I10)')IQCHNGTYP(II)          
          ELSE
            LLOC=1
            CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,ND,R,IOUT,IN)
            CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,IDUM,QRTF(4,II),IOUT,
     &                  IN)
            CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NumRT,R,IOUT,IN)
            IF (NumRT.EQ.0 .AND. NREAD1.EQ.6 .AND. NAUX.EQ.0) THEN
              QRTF(5,II) = 0.0
            ELSE
              DO 201 JJ=6,NREAD1
                CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,IDUM,QRTF(JJ,II),
     &                      IOUT,IN)
  201         CONTINUE
            ENDIF
            IF(IQCHANGEC.NE.0) 
     &     CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IQCHNGTYP(II),R,IOUT,IN)
          ENDIF
          QRTF(5,II) = NumRT
        ENDIF
        QRTF(1,II)= ND
        DO 501 ILOC=ISCLOC1,ISCLOC2
          QRTF(ILOC,II)=QRTF(ILOC,II)*SFAC
  501   CONTINUE
C
C  Read the optional values from the line
        IF (NAUX.GT.0) THEN
          DO 601 JJ=NREAD1+1,NREAD2
            CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,IDUM,QRTF(JJ,II),IOUT,
     &                  IN)
  601     CONTINUE
        ENDIF
C
C  Write the values that were read and that are not related to
C  return flow.
        NN=II-LSTBEG+1
        IF (IQRTFL.EQ.0) THEN
          WRITE(IOUT,5701) NN,ND,(QRTF(JJ,II),JJ=4,NREAD2)
 5701     FORMAT(1X,I6,I7,14G16.4)
        ELSE
          IF (NREAD2.GE.6) THEN
            WRITE(IOUT,5701) NN,ND,QRTF(4,II),
     &                      (QRTF(JJ,II),JJ=10,NREAD2)
          ELSE
            WRITE(IOUT,5701) NN,ND,QRTF(4,II)
          ENDIF
        ENDIF
C
C  Check for illegal grid location
        IF (ND.LT.1 .OR. ND.GT.NEQS) THEN
          WRITE(IOUT,*) ' ERROR: Node number is outside of the grid'
          IERR = 1
        ENDIF
        IF (IERR.NE.0) CALL USTOP(' ')
 1001 CONTINUE
C
C     Check and write data related to return-flow recipient cells
      IF (IQRTFL.GT.0) THEN
        WRITE(IOUT,'(/,1X,A,/,1X,A)') LABEL2,LABEL5
        NN = 0
        DO 1101 II=LSTBEG,N
          NN = NN + 1
          ND = QRTF(1,II)
          NumRT = QRTF(5,II)
          RFP = QRTF(6,II)
          QRTF(NQRTVL,II) = QRTF(4,II)
          WRITE(IOUT,6001) NN,ND,NumRT,RFP
 6001     FORMAT(1X,I6,I7,I7,2X,F8.6)
C
C  Check for invalid return-flow proportion
          IF (NumRT.NE.0) THEN
            IF (RFP.LT.0.0 .OR. RFP.GT.1.0) THEN
              WRITE(IOUT,5901)
 5901         FORMAT(' ERROR: Proportion must be between 0.0 and 1.0')
              IERR = 1
            ENDIF
          ENDIF
          IF (IERR.NE.0) CALL USTOP(' ')
 1101   CONTINUE
      ENDIF

C----------------------------------------------------------------------------
      ENDIF
C----------------------------------------------------------------------------
C ----READ RETURN FLOW NODES IF RETURN FLOW IS ACTIVE
      IF(IQRTFL.NE.0)THEN
        IRT = 1     ! INDEX FOR LOCATION IN NodQRT ARRAY
        DO II=LSTBEG,N
          NumRT = QRTF(5,II)
          IF(NumRT.EQ.0) CYCLE
          ItotRT = NumRT + IRT - 1
          CALL U1DINT(NodQRT(IRT),ANAME,NumRT,II,IN,IOUT)
C          IRT = IRT + NumRT - 1
C --------CHECK FOR INCORRECT GRID LOCATION
          DO IT =1,NumRT
            IRT = IRT + 1
            NQRT = NodQRT(IT)
            IF (NQRT.LT.0 .OR. NQRT.GT.NEQS) THEN
              WRITE(IOUT,*) ' ERROR: Node number is outside of the',
     &                      ' grid'
              IERR = 1
            ENDIF
          ENDDO  
        ENDDO
      ENDIF
C
      IF (ICLOSE.NE.0) CLOSE(UNIT=IN)
C
      RETURN
      END
      SUBROUTINE SGWF2QRT8LS(IN,IOUTU,QRTF,NQRTVL,MXQRT,NREAD,MXAQRT,
     &       NQRTCL,QRTAUX,NCAUX,NAUX,IQRTFL,IUNSTR,NodQRT,MXRTCELLS)
C     ******************************************************************
C     Read a list parameter name, look it up in the list of parameters,
C     and substitute values into active part of package array.
C     ******************************************************************
C     Modified 11/8/2001 to support parameter instances - ERB
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE PARAMMODULE
      CHARACTER*3 PACK, PTYP
      DIMENSION QRTF(NQRTVL,MXQRT),NodQRT(MXRTCELLS)
      CHARACTER*57 LABEL1, LABEL2, LABEL3, LABEL4, LABEL5
      CHARACTER*16 QRTAUX(NCAUX)
      CHARACTER*400 LINE
      CHARACTER*10 CTMP1, CTMP2, CTMP3, CTMP4
      CHARACTER*24 ANAME
C
      DATA ANAME /'       RETURN FLOW NODES'/          
C     ------------------------------------------------------------------
  500 FORMAT(/,' Parameter:  ',A)
  510 FORMAT(1X,'Parameter type conflict:',/
     &       1X,'Named parameter:',A,' was defined as type:',A,/
     &       1X,'However, this parameter is used in the ',A,
     &       ' file, so it should be type:',A)
  512 FORMAT(/,1X,'Blank instance name in the ',A,
     &       ' file for parameter ',A)
  514 FORMAT(3X,'Instance:  ',A)
  516 FORMAT(/,1X,'The ',A,' file specifies undefined instance "',
     &       A,'" for parameter ',A)
  520 FORMAT(1X,/1X,'THE NUMBER OF ACTIVE LIST ENTRIES (',I6,
     &       ')',/1X,'IS GREATER THAN THE MAXIMUM ALLOWED (',I6,')')
  530 FORMAT(1X,I6,I7,I7,I7,14G16.4)
 5301 FORMAT(1X,I6,I7,14G16.4)
  550 FORMAT(/,1X,'*** ERROR: PARAMETER "',A,
     &'" HAS ALREADY BEEN ACTIVATED THIS STRESS PERIOD',/,
     &' -- STOP EXECUTION (SGWF2QRT7LS)')
  600 FORMAT(1X,I6,4I7,2X,F8.6)
 6001 FORMAT(1X,I6,I7,I7,2X,F8.6)
C
      PACK = 'QRT'
      PTYP = 'QRT'
      IPVL1 = 5
      IPVL2 = 5
      LABEL1='SINK NO.  LAYER   ROW   COL     SinkQ   '
      LABEL2='          ----SINK CELL----  --RECIPIENT CELL--   RETURN'
      LABEL3='SINK NO.  LAYER   ROW   COL   LAYER   ROW   COL    PROP.'
      LABEL4='SINK NO.  NODE     SinkQ      STRESS FACTOR'
      LABEL5='SINK NO.       NODE       NO. OF RECIPIENT CELLS   PROP.'
      IOUT = ABS(IOUTU)
C
      READ(IN,'(A)') LINE
      LLOC=1
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,0,IDUM,RDUM,IOUT,IN)
      WRITE(IOUT,500) LINE(ISTART:ISTOP)
      IF(LINE(ISTART:ISTOP).EQ.' ') THEN
        WRITE(IOUT,*) ' Blank parameter name in the ',PACK,' file.'
        CALL USTOP(' ')
      END IF
C
      CTMP1=LINE(ISTART:ISTOP)
      CALL UPCASE(CTMP1)
      DO 100 IP=1,MXPAR
        CTMP2=PARNAM(IP)
        CALL UPCASE(CTMP2)
        IF(CTMP1.EQ.CTMP2) THEN
          IF(PARTYP(IP).NE.PTYP) THEN
            WRITE(IOUT,510) PARNAM(IP),PARTYP(IP),PACK,PTYP
            CALL USTOP(' ')
          ENDIF
C
C         DESIGNATE CELLS CORRESPONDING TO CORRECT PARAMETER INSTANCE
          NLST=IPLOC(2,IP)-IPLOC(1,IP)+1
          NUMINST=IPLOC(3,IP)
          ILOC=IPLOC(4,IP)
          NI=1
          IF(NUMINST.GT.0) THEN
            NLST=NLST/NUMINST
            CALL URWORD(LINE,LLOC,ISTART,ISTOP,0,IDUM,RDUM,IOUT,IN)
            CTMP3=LINE(ISTART:ISTOP)
            IF(CTMP3.EQ.' ') THEN
              WRITE(IOUT,512)PACK,PARNAM(IP)
              CALL USTOP(' ')
            ENDIF
            WRITE(IOUT,514) CTMP3
            CALL UPCASE(CTMP3)
            DO 10 KI=1,NUMINST
              CTMP4=INAME(ILOC+KI-1)
              CALL UPCASE(CTMP4)
              IF(CTMP3.EQ.CTMP4) THEN
                NI=KI
                GOTO 15
              ENDIF
   10       CONTINUE
            WRITE(IOUT,516) PACK,CTMP3,PARNAM(IP)
            CALL USTOP(' ')
   15       CONTINUE
          ENDIF
C
          IF (IACTIVE(IP).GT.0) THEN
            WRITE(IOUT,550) PARNAM(IP)
            CALL USTOP(' ')
          ENDIF
C
          IACTIVE(IP)=NI
C
          NQRTCL=NQRTCL+NLST
          IF(NQRTCL.GT.MXAQRT) THEN
            WRITE(IOUT,520) NQRTCL,MXAQRT
            CALL USTOP(' ')
          ENDIF
C
C  Write label for list values
          IF(IUNSTR.EQ.0) THEN
            IF (IOUTU.GT.0) CALL ULSTLB(IOUT,LABEL1,QRTAUX,NCAUX,NAUX)
          ELSE
            IF (IOUTU.GT.0) CALL ULSTLB(IOUT,LABEL4,QRTAUX,NCAUX,NAUX)
          ENDIF
C
C  Substitute values
          DO 60 I=1,NLST
            II=NQRTCL-NLST+I
            III=I-1+IPLOC(1,IP)+(NI-1)*NLST
            DO 20 J=1,NREAD
              QRTF(J,II)=QRTF(J,III)
   20       CONTINUE
            DO 40 IPVL=IPVL1,IPVL2
              QRTF(IPVL,II)=QRTF(IPVL,II)*B(IP)
   40       CONTINUE
C  SET VALUES FOR STRUCTURED GRID
            IF(IUNSTR.EQ.0) THEN
              IL=QRTF(1,II)
              IR=QRTF(2,II)
              IC=QRTF(3,II)
              IF (IOUTU.GT.0) THEN
                IF (IQRTFL.EQ.0) THEN
                  WRITE(IOUT,530) II,IL,IR,IC,(QRTF(JJ,II),JJ=4,NREAD)
                ELSE
                  IF (NREAD.GE.10) THEN
                    WRITE(IOUT,530) II,IL,IR,IC,QRTF(4,II),
     &                              (QRTF(JJ,II),JJ=10,NREAD)
                  ELSE
                    WRITE(IOUT,530) II,IL,IR,IC,QRTF(4,II)
                  ENDIF
                ENDIF
              ENDIF
C  SET VALUES FOR UNSTRUCTURED GRID
            ELSE
              ND=QRTF(1,II)
              IF (IOUTU.GT.0) THEN
                IF (IQRTFL.EQ.0) THEN
                  WRITE(IOUT,5301) II,ND,(QRTF(JJ,II),JJ=4,NREAD)
                ELSE
                  IF (NREAD.GE.10) THEN
                    WRITE(IOUT,5301) II,ND,QRTF(4,II),
     &                              (QRTF(JJ,II),JJ=10,NREAD)
                  ELSE
                    WRITE(IOUT,5301) II,ND,QRTF(4,II)
                  ENDIF
                ENDIF
              ENDIF
C
            ENDIF
   60     CONTINUE
          GOTO 120
        ENDIF
  100 CONTINUE
C
      WRITE(IOUT,*) ' The ',PACK,
     &   ' file specifies an undefined parameter:',LINE(ISTART:ISTOP)
      CALL USTOP(' ')
C
  120 CONTINUE
C
C----------------------------------------------------------------------------
C ----READ RETURN FLOW NODES IF RETURN FLOW IS ACTIVE
C      IF(IQRTFL.NE.0)THEN
C        IRT = 1     ! INDEX FOR LOCATION IN NodQRT ARRAY
C        DO II=1,NLST
C          NumRT = QRTF(5,II)
C          IF(NumRT.EQ.0) CYCLE
C          ItotRT = NumRT + IRT - 1
C          CALL U1DINT(NodQRT(IRT),ANAME,ItotRT,II,IN,IOUT)
CC          IRT = IRT + NumRT - 1
CC --------CHECK FOR INCORRECT GRID LOCATION
C          DO IT =1,NumRT
C            IRT = IRT + 1
C            NQRT = NodQRT(IT)
C            IF (NQRT.LT.0 .OR. NQRT.GT.NEQS) THEN
C              WRITE(IOUT,*) ' ERROR: Node number is outside of the',
C     &                      ' grid'
C              IERR = 1
C            ENDIF
C          ENDDO  
C        ENDDO
C      ENDIF
CC -------------------------------------------------------------------------      
      IF (IQRTFL.GT.0 .AND. IOUTU.GT.0) THEN
C     WRITE DATA RELATED TO RETURN-FLOW RECIPIENT CELLS FOR STRUCTURED GRIDS
        IF(IUNSTR.EQ.0) THEN
          WRITE(IOUT,'(/,1X,A,/,1X,A)') LABEL2,LABEL3
          NN = 0
          DO 140 II=NQRTCL-NLST+1,NQRTCL
            NN = NN + 1
            K = QRTF(1,II)
            I = QRTF(2,II)
            J = QRTF(3,II)
            NumRT = QRTF(5,II)
            RFP = QRTF(6,II)
            WRITE(IOUT,600) NN,K,I,J,NumRT,RFP
  140     CONTINUE
C     WRITE DATA RELATED TO RETURN-FLOW RECIPIENT CELLS FOR UNSTRUCTURED GRIDS
        ELSE
          WRITE(IOUT,'(/,1X,A,/,1X,A)') LABEL2,LABEL5
          NN = 0
          DO 1401 II=NQRTCL-NLST+1,NQRTCL
            NN = NN + 1
            ND = QRTF(1,II)
            NR = QRTF(5,II)
            RFP = QRTF(6,II)
            WRITE(IOUT,6001) NN,ND,NR,RFP
 1401     CONTINUE
        ENDIF
C
      ENDIF
C
      RETURN
      END
C -------------------------------------------------------------      
      SUBROUTINE GWF2QRT8U1AD
C     ******************************************************************
C     COMPUTE SPECIFIED QRT FLOWS AT CURRENT TIME STEP
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,      ONLY:ITRSS,HNEW,HOLD,AREA,IBOUND
      USE GWFBASMODULE,ONLY:DELT,TOTIM,TMINAT
      USE GWFQRTMODULE,ONLY:NBDQTIM,ISTEPQ,IQRTN,BDQTIM,BDQV,QRTF,
     1  TIMQRT,MXQRT,QRTFLOW,RTAREA,NQRTCL,IQRTFL,NODQRT
      DOUBLE PRECISION QFACT1,QFACT2,T1,T2,TP,DDT
C     ------------------------------------------------------------------
C
C1 -----IF QRT IS NOT TRANSIENT THEN RETURN      
      IF(NBDQTIM.EQ.0) RETURN
C      
C2------IF NBDQTIM = 1 THEN SET INITIAL VALUE AND RETURN
      IF(NBDQTIM.EQ.1) THEN
         IF (MXQRT.GT.0) THEN
           DO 5 NF=1,MXQRT
             QRTF(4,NF)=BDQV(1,NF)
5          CONTINUE
         ENDIF
         RETURN
      ENDIF  
C
C3------FIND ARRAY INDICES OF TIMES AROUND TIME AT START OF CURRENT
C3------TIME STEP
      IB4=0
      IF(MXQRT.LT.1) GO TO 200
      T2=TOTIM
      T1=TOTIM-DELT
      DO 20 L=2,NBDQTIM
      IF(T1.LE.BDQTIM(L)) THEN
         IB1=L-1
         IB2=L
         GO TO 40
      ENDIF
 20   CONTINUE
      IB1=NBDQTIM-1
      IB2=NBDQTIM
C
C3------COMPUTE FACTOR FOR INTERPOLATION OR EXTRAPOLATION OF FLOW AT
C3------START OF CURRENT TIME STEP
 40   QFACT1=(T1-BDQTIM(IB1))/(BDQTIM(IB2)-BDQTIM(IB1))
      IF(ISTEPQ.EQ.1) QFACT1 = 0.
C
C4------FIND ARRAY INDICES OF TIMES AROUND TIME AT END OF CURRENT
C4------TIME STEP
      DO 60 L=IB2,NBDQTIM
      IF(T2.LE.BDQTIM(L)) THEN
         IB3=L-1
         IB4=L
         GO TO 70
      ENDIF
 60   CONTINUE
      IB4=NBDQTIM
      IB3=NBDQTIM-1
C
C5------COMPUTE FACTOR FOR INTERPOLATION OR EXTRAPOLATION OF FLOW AT
C5------END OF CURRENT TIME STEP
 70   QFACT2=(T2-BDQTIM(IB3))/(BDQTIM(IB4)-BDQTIM(IB3))
      IF(ISTEPQ.EQ.1) QFACT2 = 0.
C
C6------COMPUTE SPECIFIED FLOW RATES FOR THIS TIME STEP
      NPI=IB4-IB2
      DO 90 NF=1,MXQRT
      QA=BDQV(IB1,NF)
      QB=BDQV(IB2,NF)
      QC=BDQV(IB3,NF)
      QD=BDQV(IB4,NF)
      Q1=(QA+QFACT1*(QB-QA))
      Q2=(QC+QFACT2*(QD-QC))
      IF(NPI.EQ.0) THEN
         QRTF(4,NF)=0.5*(Q1+Q2)
      ELSE
         TP=T1
         QP=Q1
         SUM1=0.0
         DO 80 NI=IB2,IB3
         QN=BDQV(NI,NF)
         DDT=BDQTIM(NI)-TP
         SUM1=SUM1+DDT*0.5*(QN+QP)
         TP=BDQTIM(NI)
         QP=QN
 80      CONTINUE
         DDT=T2-TP
         SUM1=SUM1+DDT*0.5*(Q2+QP)
         QRTF(4,NF)=SUM1/DELT
         QRTF(NQRTVL,NF) = QRTF(4,NF)
      ENDIF
 90   CONTINUE
C
C7-----COMPUTE VALUES OF AUXILIARY VARIABLES FOR SPECIFIED-FLOW
C7-----CELLS FOR CURRENT TIME STEP

 200  CONTINUE
C ----------------------------------------------------------------------
C8 ----- RESET RETURN FLOW FLUXES IN QRTFLOW    
      IF (IQRTFL.GT.0) THEN
        IRT = 0 ! SET POINTER FOR RETURN FLOW NODES ARRAY
        DO  L=1,NQRTCL          
C ------------PARTITION SINK ONTO RETURN FLOW NODES IN QRTFLOW   
          NumRT = QRTF(5,L)  
          IF (NumRT.EQ.0) CYCLE
          FR = QRTF(6,L) / RTAREA(L) ! FRACTURN RETURNED DIVIDED BY AREA          
          DO JJ = 1,NumRT
            IRT = IRT + 1  
            INR = NodQRT(IRT)
            IF (IBOUND(INR) .GT. 0) THEN
              QRTFLOW(IRT) = FR * QRTF(4,L)* AREA(INR)
            END IF
          ENDDO
        ENDDO 
      ENDIF
C --------------------------------------------------------------------------     
C11-----FILL TIMQRT FOR ATS
      IF(ISTEPQ.EQ.1) THEN
        IF(IB4.GT.0) TIMQRT = BDQTIM(IB4)
C11B------ADVANCE TIMQRT WHEN TIME EXCEEDS PREVIOUS TIMQRT        
        IF(TOTIM.GE.TIMQRT-TMINAT)THEN
          IF(IB4.LT.NBDQTIM)THEN
            TIMQRT = BDQTIM(IB4+1)
          ELSE
            TIMQRT = 1.0E30
          ENDIF
        ENDIF
      ENDIF
C
C11------RETURN
      RETURN
      END      
C-----------------------------------------------------------------------      
      SUBROUTINE GWF2QRT8U1DA
C  Deallocate QRT MEMORY
      USE GWFQRTMODULE
C
        DEALLOCATE(NQRTCL)
        DEALLOCATE(MXQRT)
        DEALLOCATE(NQRTVL)
        DEALLOCATE(NQRTNP)
        DEALLOCATE(IQRTCB)
        DEALLOCATE(NPQRT)
        DEALLOCATE(IQRTPB)
        DEALLOCATE(IQRTFL)
        DEALLOCATE(NOPRQT)
        DEALLOCATE(QRTF)
        DEALLOCATE(QRTAUX)
        DEALLOCATE(NodQRT)
        DEALLOCATE(QRTFLOW)
        IF(NBDQTIM.GT.0) THEN
         DEALLOCATE(IQRTN,BDQTIM,BDQV)   
        ENDIF
        IF(ISTEPQ.NE.0) THEN
          DEALLOCATE (TIMQRT) 
        ENDIF  
        DEALLOCATE(NBDQTIM,ISTEPQ)
C
      RETURN
      END

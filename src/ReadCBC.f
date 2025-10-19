      SUBROUTINE READCBC(KSTP,KPER)
C     ******************************************************************
C     READ CBC FILES FOR RUNNING TRANSPORT FROM PREVIOUS FLOW RUN (ITRNSP=2)
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL, ONLY:NCOL,NROW,NLAY,IBOUND,HNEW,BUFF,NODLAY,NEQS,
     1  TOP,IOUT,NODES,NJA,IA,JA,JAS,IUNSTR,IVC,ISYM,ITRNSP,issflg,
     1  Sn,So,INGNC,INGNC2,INGNCn,FLOWJA,JAFL,NJAG,iunsat,incln,IUNIT
      USE GWFBASMODULE,ONLY:ICBCFL,DELT,PERTIM,TOTIM,ICHFLG
      USE GWFBCFMODULE,ONLY:IBCFCB,LAYCON
      USE CLN1MODULE, ONLY: ICLNCB,NCLN,NNDCLN,CLNCON,NCLNNDS,ACLNNDS,
     1    NCLNGWC,ACLNGWC,IA_CLN,JA_CLN,NJA_CLN,IDXGLO_CLN,ICLNGWCB,
     2    ICLNPCB
      USE GWTBCTMODULE, ONLY: CBCF
C
      CHARACTER*16 TEXT(3)
      REAL DT,PTIM,DTIM
      REAL, ALLOCATABLE  :: FLOWGWS(:,:),FLOWJAG(:),FLOWCLNCLN(:)
      DATA TEXT(1) /'   FLOW JA FACE '/
      DATA TEXT(2) /'   FLOW CLN FACE'/
      DATA TEXT(3) /'             GWF'/
      LOGICAL FOUND  !kkz - need to explicitly declare FOUND as a LOGICAL for GFortran
c
C-----------------------------------------------------------------------------
C1------FILL STORAGE TERM
      CALL STORFLXRD(KSTP,KPER)
C1------FILL CONSTANT HEAD BOUNDARY FLUX FIRST
      CALL CHDFLUXRD(KSTP,KPER)
C
C-----------------------------------------------------------------------------
C2------FOR GW CELLS READ FLOWJA(II) TO FILL CELL BY CELL FLOWS INTO CBCF(IIS)
C-----------------------------------------------------------------------------
      IBD=0
      IF(IBCFCB.GT.0) IBD=ICBCFL
      IF(IBD.EQ.0)THEN
        WRITE(IOUT,61)
61      FORMAT(1X,'*** IBD = 0 FLOW INFORMATION WAS NOT WRITTEN, '
     1   'STOPPING ***')
        STOP
      ENDIF
      ZERO = 0.
C
      ALLOCATE(FLOWJAG(NJAG),FLOWJA(NJA))
      FLOWJA = ZERO
      FLOWJAG = ZERO
C2A-------READ ARRAY FLOWJAG
      IF(IBD.EQ.1)
     1   CALL UBUDSVUrd(KSTP,KPER,TEXT(1),IBCFCB,FLOWJAG(1),NJAG,IOUT,
     1         PERTIM,TOTIM)
      IF(IBD.EQ.2) CALL UBDSV1Urd(KSTP,KPER,TEXT(1),IBCFCB,FLOWJAG(1),
     1     NJAG,IOUT,DELT,PERTIM,TOTIM,IBOUND,NODES)
C
C2B------STORE CBC ARRAY IN FULL IA AND JA STRUCTURE FOR GW NODES
      IJAG = 1
      DO N=1,NODES
        DO II = IA(N),IA(N+1)-1
          JJ = JA(II)
          JJG = JAFL(IJAG)
          IF(JJ.NE.JJG) CYCLE
          FLOWJA(II) = -FLOWJAG(IJAG)
          IJAG = IJAG + 1
        ENDDO
      ENDDO
      DEALLOCATE(FLOWJAG)
C
C-----------------------------------------------------------------------------
C3------FOR CLN-CLN FLOW READ FLOWJA(II) TO FILL CELL BY CELL FLOWS INTO CBCF(IIS)
C-----------------------------------------------------------------------------
      IF(INCLN.EQ.0) GO TO 22
C-----------------------------------------------------------------------------
      IBD=0
      IF(ICLNCB.GT.0) IBD=ICBCFL
C
C3A-------READ ARRAY FLOWCLNCLN
      LCLN = NJA_CLN
      ALLOCATE(FLOWCLNCLN(LCLN))
      FLOWCLNCLN = ZERO
      IF(IBD.EQ.1)
     1  CALL UBUDSVUrd(KSTP,KPER,TEXT(2),ICLNCB,FLOWCLNCLN,LCLN,IOUT,
     1         PERTIM,TOTIM)
      IF(IBD.EQ.2)
     1  CALL UBDSV1Urd(KSTP,KPER,TEXT(2),ICLNCB,FLOWCLNCLN,LCLN,IOUT,
     2         DELT,PERTIM,TOTIM,IBOUND(NODES+1),NCLNNDS)
C
C3B------STORE CBC ARRAY IN FULL IA AND JA STRUCTURE FOR CLN NODES
      DO NC1 = 1,NCLNNDS
C4B----------loop over all connections of node NC1
        DO II_CLN = IA_CLN(NC1)+1,IA_CLN(NC1+1)-1
          NC2 = JA_CLN(II_CLN)
          ND1 = ACLNNDS(NC1,1)   !  NC1 + NODES
          ND2 = ACLNNDS(NC2,1)   !  NC2 + NODES
          IF(IBOUND(ND1).EQ.0.OR.IBOUND(ND2).EQ.0) CYCLE
          II = IDXGLO_CLN(II_CLN)
            IF(ICHFLG.EQ.0) THEN
              IF((IBOUND(ND1).LE.0) .AND. (IBOUND(ND2).LE.0)) CYCLE
            END IF
            FLOWJA(II) = -FLOWCLNCLN(II_CLN)
        ENDDO
      ENDDO
      DEALLOCATE(FLOWCLNCLN)
C
C-----------------------------------------------------------------------------
C4------FOR CLN-GWF FLOW READ FLOWJA(II) TO FILL CELL BY CELL FLOWS INTO CBCF(IIS)
C-----------------------------------------------------------------------------
C4A-----STOP IF IBD IS 1 SINCE INFORMATION IS LOST ON MULTIPLE CONNECTIONS IN CBC FILE
      IF(IBD.EQ.1)THEN
        WRITE(IOUT,20)
20      FORMAT(1X,'*** IBD=1 LOSES INFORMATION UNLESS ONE-TO-ONE ',
     1   'CONNECTION ONLY BETWEEN CLN AND GWF, STOPPING')
        STOP
      ENDIF
      IF(IBD.NE.2) GO TO 22
C
C4B-----READ HEADER FOR COMPACT BUDGET
      CALL UBDSV2Urd(KSTP,KPER,TEXT(3),ICLNCB,NCLNNDS,
     1          NCLNGWC,IOUT,DELT,PERTIM,TOTIM,IBOUND(NODES+1))
C4C-----READ EACH CLN-GWF CONNECTION FLOW
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
cccc        IF(.NOT.FOUND) CALL USTOP('error in CLN1BDGWFWR')
        CALL UBDSVAUrd(ICLNCB,NCLNNDS,IH,SRATE,IBOUND(NODES+1))
        FLOWJA(II) = -SRATE
      ENDDO
C--------------------------------------------------------------------------------
22    CONTINUE
C--------------------------------------------------------------------------------
C5-----FILL FLOW TERM FROM FLOWJA ARRAY INTO SYMMETRIC CBCF ARRAY
C--------------------------------------------------------------------------------
      DO N=1,NEQS
        DO II = IA(N)+1,IA(N+1)-1
          JJ = JA(II)
          IIS = JAS(II)
          IF(JJ.LE.N) CYCLE ! FILL  UPPER TRIANGLE
          IF(IIS.GT.0) CBCF(IIS) = FLOWJA(II)
        ENDDO
      ENDDO
C
      DEALLOCATE(FLOWJA)
C--------------------------------------------------------------------------------
C6-----FILL FLOW TERM FOR VARIOUS BOUNDARY CONDITIONS IN SAME ORDER AS BD IN MFUSG
C--------------------------------------------------------------------------------
C6A------FILL WELL BOUNDARY FLUX
      IF(IUNIT(2).GT.0) CALL WELFLUXRD(KSTP,KPER)
C6B------FILL DRN BOUNDARY FLUX
      IF(IUNIT(3).GT.0) CALL DRNFLUXRD(KSTP,KPER)
C6C------FILL RIV BOUNDARY FLUX
      IF(IUNIT(4).GT.0) CALL RIVFLUXRD(KSTP,KPER)
C6D------FILL EVT BOUNDARY FLUX
      IF(IUNIT(5).GT.0) CALL EVTFLUXRD(KSTP,KPER)
C6E------FILL GHB BOUNDARY FLUX
      IF(IUNIT(7).GT.0) CALL GHBFLUXRD(KSTP,KPER)
C6F------FILL RCH BOUNDARY FLUX
      IF(IUNIT(8).GT.0) CALL RCHFLUXRD(KSTP,KPER)
C6G------FILL STR BOUNDARY FLUX
      IF(IUNIT(18).GT.0) CALL STRFLUXRD(KSTP,KPER)
C6H------FILL ETS BOUNDARY FLUX
      IF(IUNIT(39).GT.0) CALL ETSFLUXRD(KSTP,KPER)
C6I------FILL DRT BOUNDARY FLUX
      IF(IUNIT(40).GT.0) CALL DRTFLUXRD(KSTP,KPER)
C6J------FILL QRT BOUNDARY FLUX
      IF(IUNIT(41).GT.0) CALL QRTFLUXRD(KSTP,KPER)
C6K------FILL SFR BOUNDARY FLUX
      IF(IUNIT(44).GT.0) CALL SFRFLUXRD(KSTP,KPER)
C6L------FILL LAK BOUNDARY FLUX
      IF(IUNIT(22).GT.0) CALL LAKFLUXRD(KSTP,KPER)
C6M------FILL WHATEVER OTHER BOUNDARY FLUX HERE


C--------------------------------------------------------------------------------
C
C7------RETURN
      RETURN
      END
C-----------------------------------------------------------------------------------
      SUBROUTINE STORFLXRD(KSTP,KPER)
C     ******************************************************************
C     READ STORAGE FLOW TERM. NOT USED FOR TRANSPORT SO NOT STORED AFTER READING
C     ******************************************************************
C
C     SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL, ONLY:NCOL,NROW,NLAY,ISSFLG,IBOUND,HNEW,HOLD,FLOWJA,IA,
     1      BUFF,TOP,IOUT,NODES,NODLAY,IUNSTR,Sn,So,TOP,BOT,iunsat
      USE GWFBASMODULE,ONLY:MSUM,ICBCFL,VBVL,VBNM,DELT,PERTIM,TOTIM
      USE GWFBCFMODULE,ONLY:IBCFCB,LAYCON,SC1,SC2
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
C2------IF STEADY STATE, STORAGE TERM IS ZERO
      IF(ISS.NE.0) GOTO 400
C
C3------IF CELL-BY-CELL FLOWS WILL BE SAVED, SET FLAG IBD.
      IBD=0
      IF(IBCFCB.GT.0) IBD=ICBCFL
C
C9------record contents of buffer for structured and unstructured grids
      IF(IUNSTR.EQ.0)THEN
C
C9A-----IF IBD FLAG IS SET RECORD THE CONTENTS OF THE BUFFER.
        IF(IBD.EQ.1) CALL UBUDSVRD(KSTP,KPER,TEXT,
     1                       IBCFCB,BUFF,NCOL,NROW,NLAY,IOUT)
        IF(IBD.EQ.2) CALL UBDSV1RD(KSTP,KPER,TEXT,IBCFCB,
     1            BUFF,NCOL,NROW,NLAY,IOUT,DELT,PERTIM,TOTIM,IBOUND)
      ELSE
C
C9B-----IF IBD FLAG IS SET RECORD THE CONTENTS OF THE BUFFER.
        IF(IBD.EQ.1) CALL UBUDSVURD(KSTP,KPER,TEXT,IBCFCB,BUFF(1),NODES,
     1         IOUT,PERTIM,TOTIM)
        IF(IBD.EQ.2) CALL UBDSV1URD(KSTP,KPER,TEXT,IBCFCB,BUFF(1),NODES,
     1     IOUT,DELT,PERTIM,TOTIM,IBOUND,NODES)
      ENDIF
C
C10-----ADD TOTAL RATES AND VOLUMES TO VBVL & PUT TITLE IN VBNM.
  400 CONTINUE
C
C11----RETURN.
      RETURN
      END
C--------------------------------------------------------------------------------------
C--------------------------------------------------------------------------------
      SUBROUTINE CHDFLUXRD(KSTP,KPER)
C     ******************************************************************
C     READ CONSTANT HEAD FLUX FOR RUNNING TRANSPORT FROM PREVIOUS FLOW RUN (ITRNSP=2)
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL, ONLY:NCOL,NROW,NLAY,IBOUND,HNEW,BUFF,ITRNSP,NOVFC,
     1   TOP,IOUT,NODES,NEQS,NODLAY,IA,JA,JAS,IUNSTR,IVC,ISYM,INCLN,
     2   FLOWJA,FMBE
      USE GWFBASMODULE,ONLY:MSUM,VBVL,VBNM,DELT,PERTIM,TOTIM,ICBCFL,
     1                      ICHFLG
      USE GWFBCFMODULE,ONLY:IBCFCB,LAYCON
      USE GWTBCTMODULE, ONLY: CBCH
C
      CHARACTER*16 TEXT(1)
      DATA TEXT(1) /'   CONSTANT HEAD'/
C     ------------------------------------------------------------------
C
C1------SET IBD TO INDICATE IF CELL-BY-CELL BUDGET VALUES WILL BE SAVED.
      IBD=0
      ZERO = 0.0
      IF(IBCFCB.GT.0) IBD=ICBCFL
C2-----STOP IF IBD IS 1 NOT CODED AS OTHER BOUNDARIES MAY MISS INFORMATION ON IBD=1
      IF(IBD.EQ.1)THEN
        WRITE(IOUT,60)
60      FORMAT(1X,'*** IBD=1 IS NOT SUPPORTED PLEASE USE COMPACT ',
     1   'BUDGET FORMATS WITH ICBCFL = 2. STOPPING')
        STOP
      ENDIF
      IF(IBD.NE.2) GO TO 200
      CALL UBDSV2Urd(KSTP,KPER,TEXT(1),IBCFCB,NODES,
     1         NCH,IOUT,DELT,PERTIM,TOTIM,IBOUND)
CC      DO N=1,NODES
CC        IF(IBOUND(N).GE.0) CYCLE
        DO N=1,NCH
        IF(IUNSTR.EQ.0)THEN
          K = N / (NCOL*NROW) + 1
          IJ = N - (K-1)*NCOL*NROW
          I = (IJ-1)/NCOL + 1
          J = IJ - (I-1)*NCOL
cc          CALL UBDSVArd(IBCFCB,NCOL,NROW,J,I,K,SRATE,IBOUND,NLAY) ! not for structured grid
        ELSE
          CALL UBDSVAUrd(IBCFCB,NODES,ND,SRATE,IBOUND)
        ENDIF
        CBCH(ND) = SRATE
      ENDDO
C--------------------------------------------------------------------------------
200   CONTINUE
C--------------------------------------------------------------------------------
C
C7------RETURN
      RETURN
      END
C--------------------------------------------------------------------------------
      SUBROUTINE WELFLUXRD(KSTP,KPER)
C     ******************************************************************
C     READ WEL FLUX FOR RUNNING TRANSPORT FROM PREVIOUS FLOW RUN (ITRNSP=2)
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
C
      USE GLOBAL,      ONLY:IOUT,NCOL,NROW,NLAY,IBOUND,BUFF,NODES,
     1                  IUNSTR,TOP,BOT,HNEW,NEQS,INCLN,FMBE
      USE CLN1MODULE, ONLY: ACLNNDS,NCLNNDS,ICLNCB
      USE GWFBASMODULE,ONLY:MSUM,ICBCFL,IAUXSV,DELT,PERTIM,TOTIM,
     1                      VBVL,VBNM
      USE GWFWELMODULE,ONLY:NWELLS,IWELCB,WELL,NWELVL,WELAUX,IWELQV,
     1  IAFR,WELLBOT
C
      CHARACTER*16 TEXT
      DOUBLE PRECISION RATIN,RATOUT,QQ,QTHIK,X,Y,HD,THCK,BOTT
      DATA TEXT /'           WELLS'/
C-----------------------------------------------------------------------------
C1------INITIALIZE CELL-BY-CELL FLOW TERM FLAG (IBD)
      ZERO=0.
      IBD=0
      IF(IWELCB.GT.0) IBD=ICBCFL
C2-----STOP IF IBD IS 1 SINCE INFORMATION IS LOST ON MULTIPLE WELs TO A CELL
      IF(IBD.EQ.1)THEN
        WRITE(IOUT,60)
60      FORMAT(1X,'*** IBD=1 LOSES INFORMATION UNLESS ONE ',
     1   'WELL BOUNDARY IS PROVIDED PER CELL, STOPPING')
        STOP
      ENDIF
      IF(IBD.EQ.0.AND.NWELLS.GT.0)THEN
        WRITE(IOUT,61)
61      FORMAT(1X,'*** IBD = 0 SO WELL INFORMATION WAS NOT WRITTEN, '
     1   'STOPPING')
        STOP
      ENDIF
      IF(IBD.NE.2) GO TO 200
C
C3-----READ HEADER FOR COMPACT BUDGET
      NAUX=NWELVL-5
      IF(IAUXSV.EQ.0) NAUX=0
      IICLNCB=0
      NNCLNNDS=0
      IF(INCLN.GT.0) THEN
        IICLNCB=ICLNCB
        NNCLNNDS=NCLNNDS
      ENDIF
      CALL UBDSVHDRrd(IUNSTR,KSTP,KPER,IOUT,IWELCB,IICLNCB,NODES,
     1    NNCLNNDS,NCOL,NROW,NLAY,NWELLS,NWELVL,NAUX,IBOUND,
     2    TEXT,WELAUX,DELT,PERTIM,TOTIM,WELL)
C
C4-----READ EACH WEL FLOW
      IF(NWELLS.LE.0) GO TO 200
      DO L=1,NWELLS
        N=WELL(1,L)
        IF(IBOUND(N).EQ.0) CYCLE
        Q=ZERO
        CALL UBDSVRECrd(IUNSTR,N,NODES,NNCLNNDS,IWELCB,IICLNCB,NWELVL,
     1    5,NAUX,Q,WELL(:,L),IBOUND,NCOL,NROW,NLAY)
        WELL(NWELVL,L)=Q
      ENDDO
C--------------------------------------------------------------------------------
200   CONTINUE
C--------------------------------------------------------------------------------
C
C5------RETURN
      RETURN
      END
C--------------------------------------------------------------------------------
      SUBROUTINE DRNFLUXRD(KSTP,KPER)
C     ******************************************************************
C     READ DRN FLUX FOR RUNNING TRANSPORT FROM PREVIOUS FLOW RUN (ITRNSP=2)
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,      ONLY:IOUT,NCOL,NROW,NLAY,IBOUND,HNEW,BUFF,IUNSTR,
     *                 NODES,NEQS,FMBE,INCLN
      USE CLN1MODULE,  ONLY:NCLNNDS,ICLNCB
      USE GWFBASMODULE,ONLY:MSUM,ICBCFL,IAUXSV,DELT,PERTIM,TOTIM,
     1                      VBVL,VBNM
      USE GWFDRNMODULE,ONLY:NDRAIN,IDRNCB,DRAI,NDRNVL,DRNAUX
C
      CHARACTER*16 TEXT
      DATA TEXT /'          DRAINS'/
C     ------------------------------------------------------------------
C
C1------INITIALIZE CELL-BY-CELL FLOW TERM FLAG (IBD)
      ZERO=0.
      IBD=0
      IF(IDRNCB.GT.0) IBD=ICBCFL
C2-----STOP IF IBD IS 1 SINCE INFORMATION IS LOST ON MULTIPLE DRNs TO A CELL
      IF(IBD.EQ.1)THEN
        WRITE(IOUT,60)
60      FORMAT(1X,'*** IBD=1 LOSES INFORMATION UNLESS ONE ',
     1   'DRAIN BOUNDARY IS PROVIDED PER CELL, STOPPING')
        STOP
      ENDIF
      IF(IBD.EQ.0.AND.NDRAIN.GT.0)THEN
        WRITE(IOUT,61)
61      FORMAT(1X,'*** IBD = 0 SO WELL INFORMATION WAS NOT WRITTEN, '
     1   'STOPPING')
        STOP
      ENDIF
      IF(IBD.NE.2) GO TO 200
C
C3-----READ HEADER FOR COMPACT BUDGET
      NAUX=NDRNVL-6
      IF(IAUXSV.EQ.0) NAUX=0
      IICLNCB=0
      NNCLNNDS=0
      IF(INCLN.GT.0) THEN
        IICLNCB=ICLNCB
        NNCLNNDS=NCLNNDS
      ENDIF
      CALL UBDSVHDRrd(IUNSTR,KSTP,KPER,IOUT,IDRNCB,IICLNCB,NODES,
     1    NNCLNNDS,NCOL,NROW,NLAY,NDRAIN,NDRNVL,NAUX,IBOUND,
     2    TEXT,DRNAUX,DELT,PERTIM,TOTIM,DRAI)
C
C4-----READ EACH DRN FLOW
      IF(NDRAIN.LE.0) GO TO 200
      DO L=1,NDRAIN
        N=DRAI(1,L)
        IF(IBOUND(N).EQ.0) CYCLE
        Q=ZERO
        CALL UBDSVRECrd(IUNSTR,N,NODES,NNCLNNDS,IDRNCB,IICLNCB,NDRNVL,
     1    6,NAUX,Q,DRAI(:,L),IBOUND,NCOL,NROW,NLAY)
        DRAI(NDRNVL,L)=Q
      ENDDO
C--------------------------------------------------------------------------------
200   CONTINUE
C--------------------------------------------------------------------------------
C
C5------RETURN
      RETURN
      END
C--------------------------------------------------------------------------------
      SUBROUTINE RIVFLUXRD(KSTP,KPER)
C     ******************************************************************
C     READ RIV FLUX FOR RUNNING TRANSPORT FROM PREVIOUS FLOW RUN (ITRNSP=2)
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,ONLY:IOUT,NCOL,NROW,NLAY,IBOUND,HNEW,BUFF,NODES,IUNSTR,
     1             NEQS,FMBE,INCLN
      USE CLN1MODULE,  ONLY:NCLNNDS,ICLNCB
      USE GWFBASMODULE,ONLY:MSUM,ICBCFL,IAUXSV,DELT,PERTIM,TOTIM,
     1                      VBVL,VBNM
      USE GWFRIVMODULE,ONLY:NRIVER,IRIVCB,RIVR,NRIVVL,RIVAUX
C
      CHARACTER*16 TEXT
      DATA TEXT /'   RIVER LEAKAGE'/
C     ------------------------------------------------------------------
C
C1------INITIALIZE CELL-BY-CELL FLOW TERM FLAG (IBD)
      ZERO=0.
      IBD=0
      IF(IRIVCB.GT.0) IBD=ICBCFL
C2-----STOP IF IBD IS 1 SINCE INFORMATION IS LOST ON MULTIPLE RIVs TO A CELL
      IF(IBD.EQ.1)THEN
        WRITE(IOUT,60)
60      FORMAT(1X,'*** IBD=1 LOSES INFORMATION UNLESS ONE ',
     1   'RIVER BOUNDARY IS PROVIDED PER CELL, STOPPING')
        STOP
      ENDIF
      IF(IBD.EQ.0.AND.NRIVER.GT.0)THEN
        WRITE(IOUT,61)
61      FORMAT(1X,'*** IBD = 0 SO WELL INFORMATION WAS NOT WRITTEN, '
     1   'STOPPING')
        STOP
      ENDIF
      IF(IBD.NE.2) GO TO 200
C
C3-----READ HEADER FOR COMPACT BUDGET
      NAUX=NRIVVL-7
      IF(IAUXSV.EQ.0) NAUX=0
      IICLNCB=0
      NNCLNNDS=0
      IF(INCLN.GT.0) THEN
        IICLNCB=ICLNCB
        NNCLNNDS=NCLNNDS
      ENDIF
      CALL UBDSVHDRrd(IUNSTR,KSTP,KPER,IOUT,IRIVCB,IICLNCB,NODES,
     1    NNCLNNDS,NCOL,NROW,NLAY,NRIVER,NRIVVL,NAUX,IBOUND,
     2    TEXT,RIVAUX,DELT,PERTIM,TOTIM,RIVR)
C
C4-----READ EACH RIV FLOW
      IF(NRIVER.EQ.0)GO TO 200
      DO L=1,NRIVER
        N=RIVR(1,L)
        IF(IBOUND(N).EQ.0) CYCLE
        Q = ZERO
        CALL UBDSVRECrd(IUNSTR,N,NODES,NNCLNNDS,IRIVCB,IICLNCB,NRIVVL,
     1    7,NAUX,RATE,RIVR(:,L),IBOUND,NCOL,NROW,NLAY)
        RIVR(NRIVVL,L)=Q
      ENDDO
C--------------------------------------------------------------------------------
200   CONTINUE
C--------------------------------------------------------------------------------
C
C5------RETURN
      RETURN
      END
C--------------------------------------------------------------------------------
      SUBROUTINE EVTFLUXRD(KSTP,KPER)
C     ******************************************************************
C     READ EVT FLUX FOR RUNNING TRANSPORT FROM PREVIOUS FLOW RUN (ITRNSP=2)
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,      ONLY:IOUT,NCOL,NROW,NLAY,IBOUND,HNEW,BUFF,NODES,
     *                      IA,JA,NODLAY,IUNSTR,IVC,NEQS,FMBE
      USE GWFBASMODULE,ONLY:MSUM,VBVL,VBNM,ICBCFL,DELT,PERTIM,TOTIM
      USE GWFEVTMODULE,ONLY:NEVTOP,IEVTCB,EVTR,EXDP,SURF,IEVT,
     *  NIEVT,EVTF
C
      DOUBLE PRECISION RATOUT,QQ,HH,SS,DD,XX,HHCOF,RRHS
      CHARACTER*16 TEXT
      INTEGER,ALLOCATABLE,DIMENSION(:,:) :: ITEMP
      INTEGER,ALLOCATABLE,DIMENSION(:) :: IBUFF
      DATA TEXT /'              ET'/
C     ------------------------------------------------------------------
C
C1------INITIALIZE CELL-BY-CELL FLOW TERM FLAG (IBD)
      ZERO=0.
      IBD=0
      IF(IEVTCB.GT.0) IBD=ICBCFL
C
      IF(IBD.NE.2) GO TO 200
C
C2------READ EVT FLUX ARRAY
      ALLOCATE(IBUFF(NIEVT))
      CALL UBDSV3Urd(KSTP,KPER,TEXT,IEVTCB,BUFF,IBUFF,
     1        NIEVT,NEVTOP,NODES,IOUT,DELT,PERTIM,TOTIM,IBOUND)
C3------PUT FLUX IN EVTF ARRAY
      DO NN=1,NIEVT
        N = IEVT(NN)
        EVTF(NN) = BUFF(N)
      ENDDO
C
      DEALLOCATE(IBUFF)
200   CONTINUE
C
C5------RETURN
      RETURN
      END
C--------------------------------------------------------------------------------
      SUBROUTINE GHBFLUXRD(KSTP,KPER)
C     ******************************************************************
C     READ GHB FLUX FOR RUNNING TRANSPORT FROM PREVIOUS FLOW RUN (ITRNSP=2)
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,ONLY:IOUT,NCOL,NROW,NLAY,IBOUND,HNEW,BUFF,NODES,IUNSTR,
     1            NEQS,FMBE,INCLN
      USE CLN1MODULE,  ONLY:NCLNNDS,ICLNCB
      USE GWFBASMODULE,ONLY:MSUM,ICBCFL,IAUXSV,DELT,PERTIM,TOTIM,
     1                      VBVL,VBNM
      USE GWFGHBMODULE,ONLY:NBOUND,IGHBCB,BNDS,NGHBVL,GHBAUX
C
      DOUBLE PRECISION CCGHB,CHB,RATIN,RATOUT,RRATE,HB,C
      CHARACTER*16 TEXT
      DATA TEXT /' HEAD DEP BOUNDS'/
C     ------------------------------------------------------------------
C
C1------INITIALIZE CELL-BY-CELL FLOW TERM FLAG (IBD)
      ZERO=0.
      IBD=0
      IF(IGHBCB.GT.0) IBD=ICBCFL
C2-----STOP IF IBD IS 1 SINCE INFORMATION IS LOST ON MULTIPLE RIVs TO A CELL
      IF(IBD.EQ.1)THEN
        WRITE(IOUT,60)
60      FORMAT(1X,'*** IBD=1 LOSES INFORMATION UNLESS ONE ',
     1   'GHB BOUNDARY IS PROVIDED PER CELL, STOPPING')
        STOP
      ENDIF
      IF(IBD.EQ.0.AND.NBOUND.GT.0)THEN
        WRITE(IOUT,61)
61      FORMAT(1X,'*** IBD = 0 SO WELL INFORMATION WAS NOT WRITTEN, '
     1   'STOPPING')
        STOP
      ENDIF
      IF(IBD.NE.2) GO TO 200
C
C3-----READ HEADER FOR COMPACT BUDGET
      NAUX=NGHBVL-6
      IF(IAUXSV.EQ.0) NAUX=0
      IICLNCB=0
      NNCLNNDS=0
      IF(INCLN.GT.0) THEN
        IICLNCB=ICLNCB
        NNCLNNDS=NCLNNDS
      ENDIF
      CALL UBDSVHDRrd(IUNSTR,KSTP,KPER,IOUT,IGHBCB,IICLNCB,NODES,
     1    NNCLNNDS,NCOL,NROW,NLAY,NBOUND,NGHBVL,NAUX,IBOUND,
     2    TEXT,GHBAUX,DELT,PERTIM,TOTIM,BNDS)
C
C4-----READ EACH GHB FLOW
      IF(NBOUND.EQ.0) GO TO 200
      DO L=1,NBOUND
        N=BNDS(1,L)
        IF(IBOUND(N).EQ.0) CYCLE
        RATE=ZERO
        CALL UBDSVRECrd(IUNSTR,N,NODES,NNCLNNDS,IGHBCB,IICLNCB,NGHBVL,
     1    6,NAUX,RATE,BNDS(:,L),IBOUND,NCOL,NROW,NLAY)
        BNDS(NGHBVL,L)=RATE
      ENDDO
C--------------------------------------------------------------------------------
200   CONTINUE
C--------------------------------------------------------------------------------
C
C5------RETURN
      RETURN
      END
C--------------------------------------------------------------------------------
      SUBROUTINE RCHFLUXRD(KSTP,KPER)
C     ******************************************************************
C     READ RCH FLUX FOR RUNNING TRANSPORT FROM PREVIOUS FLOW RUN (ITRNSP=2)
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,  ONLY:IOUT,NCOL,NROW,NLAY,IBOUND,BUFF,IA,JA,JAS,NODES,
     1             NODLAY,IUNSTR,IVC,hnew,NEQS,FMBE
      USE GWFBASMODULE,ONLY:MSUM,VBVL,VBNM,ICBCFL,DELT,PERTIM,TOTIM
      USE GWFRCHMODULE,ONLY:NRCHOP,IRCHCB,RECH,IRCH,NIRCH,
     1                      SELEV,ISELEV,IPONDOPT,RCHF
C
      DOUBLE PRECISION RATIN,RATOUT,QQ
      DOUBLE PRECISION RECHFLUX,acoef,eps,pe,hd,rch
      CHARACTER*16 TEXT
      INTEGER,ALLOCATABLE,DIMENSION(:,:) :: ITEMP
      INTEGER,ALLOCATABLE,DIMENSION(:) :: IBUFF
      DATA TEXT /'        RECHARGE'/
C     ------------------------------------------------------------------
C
C1------INITIALIZE CELL-BY-CELL FLOW TERM FLAG (IBD)
      ZERO=0.
      IBD=0
      IF(IRCHCB.GT.0) IBD=ICBCFL
C
      IF(IBD.NE.2) GO TO 200
C
C2------READ EVT FLUX ARRAY
      ALLOCATE(IBUFF(NIRCH))
      CALL UBDSV3Urd(KSTP,KPER,TEXT,IRCHCB,BUFF,IBUFF,
     1    NIRCH,NRCHOP,NODES,IOUT,DELT,PERTIM,TOTIM,IBOUND)
C3------PUT FLUX IN EVTF ARRAY
      DO NN=1,NIRCH
        N = IRCH(NN)
        RCHF(NN) = BUFF(N)
      ENDDO
C
      DEALLOCATE(IBUFF)
200   CONTINUE
C
C5------RETURN
      RETURN
      END
C--------------------------------------------------------------------------------
      SUBROUTINE STRFLUXRD(KSTP,KPER)
C     ******************************************************************
C     READ STR FLUX FOR RUNNING TRANSPORT FROM PREVIOUS FLOW RUN (ITRNSP=2)
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,      ONLY:IOUT,NCOL,NROW,NLAY,IBOUND,BUFF,HNEW,IUNSTR,
     *                 NODES,NEQS,FMBE
      USE GWFBASMODULE,ONLY:MSUM,VBVL,VBNM,ICBCFL,DELT,
     *                 IAUXSV,PERTIM,TOTIM
      USE GWFSTRMODULE,ONLY:MXSTRM,NSTREM,NSS,NTRIB,NDIV,ICALC,ISTCB1,
     1                      ISTCB2,IPTFLG,CONST,
     2                      STRM,ARTRIB,ISTRM,ITRBAR,IDIVAR,NDFGAR,
     3                      NSTRVL,STRAUX
C
      CHARACTER*16 TEXT,STRTXT
      DATA   TEXT/'  STREAM LEAKAGE'/
      DATA STRTXT/'STREAM FLOW OUT '/
C-----------------------------------------------------------------------------
C1------INITIALIZE CELL-BY-CELL FLOW TERM FLAG (IBD)
      ZERO=0.
      IBD=0
      IF(ISTCB1.GT.0) IBD=ICBCFL
C2-----STOP IF IBD IS 1 SINCE INFORMATION IS LOST ON MULTIPLE WELs TO A CELL
      IF(IBD.EQ.1)THEN
        WRITE(IOUT,60)
60      FORMAT(1X,'*** IBD=1 LOSES INFORMATION UNLESS ONE ',
     1   'STR BOUNDARY IS PROVIDED PER CELL, STOPPING')
        STOP
      ENDIF
      IF(IBD.EQ.0.AND.NSTREM.GT.0)THEN
        WRITE(IOUT,61)
61      FORMAT(1X,'*** IBD = 0 SO WELL INFORMATION WAS NOT WRITTEN, '
     1   'STOPPING')
        STOP
      ENDIF
      IF(IBD.NE.2) GO TO 200
C
C3-----READ HEADER FOR COMPACT BUDGET
      NAUX=NSTRVL-11
      IF(IAUXSV.EQ.0) NAUX=0
      IF(IUNSTR.EQ.0) THEN
        CALL UBDSV4rd(KSTP,KPER,TEXT,NAUX,STRAUX,ISTCB1,NCOL,NROW,NLAY,
     1       NSTREM,IOUT,DELT,PERTIM,TOTIM,IBOUND)
      ELSE
        CALL UBDSV4Urd(KSTP,KPER,TEXT,NAUX,STRAUX,ISTCB1,NEQS,
     1          NSTREM,IOUT,DELT,PERTIM,TOTIM,IBOUND)
      ENDIF
C
C4-----READ EACH STR LEAKAGE FLOW
      IF(NSTREM.LE.0) GO TO 100
      DO L=1,NSTREM
        ND=ISTRM(1,L)
        IF(IBOUND(ND).EQ.0) CYCLE
        Q=ZERO
       CALL UBDSVBUrd(ISTCB1,NEQS,ND,Q,
     1                  STRM(:,L),NSTRVL,NAUX,12,IBOUND)
        STRM(11,L) = Q
      ENDDO
C--------------------------------------------------------------------------------
100   CONTINUE
C--------------------------------------------------------------------------------
C5------IF COMPACT BUDGET, READ HEADER RECORD, BUT DO NOT INCLUDE
C5------AUX VARIABLES BECAUSE THEY WOULD BE WITH STREAM LEAKAGE.
      NAUX=0
      IF(IUNSTR.EQ.0) THEN
        CALL UBDSV4rd(KSTP,KPER,STRTXT,NAUX,STRAUX,ISTCB2,NCOL,NROW,
     1       NLAY,NSTREM,IOUT,DELT,PERTIM,TOTIM,IBOUND)
      ELSE
        CALL UBDSV4Urd(KSTP,KPER,STRTXT,NAUX,STRAUX,ISTCB2,NEQS,
     1       NSTREM,IOUT,DELT,PERTIM,TOTIM,IBOUND)
      ENDIF
C
C6------READ FOR EACH STREAM OUTFLOW
      DO L=1,NSTREM
        ND=ISTRM(1,L)
        IF(IBOUND(ND).LE.0) CYCLE
           IF(IUNSTR.EQ.0) THEN
             IL = (ND-1) / (NCOL*NROW) + 1
             IJ = ND - (IL-1)*NCOL*NROW
             IR = (IJ-1)/NCOL + 1
             IC = IJ - (IR-1)*NCOL
             CALL UBDSVBrd(ISTCB2,NCOL,NROW,IC,IR,IL,STRM(9,L),
     1                  STRM(:,L),NSTRVL,NAUX,12,IBOUND,NLAY)
          ELSE
             CALL UBDSVBU(ISTCB2,NEQS,ND,STRM(9,L),
     1                  STRM(:,L),NSTRVL,NAUX,12,IBOUND)
          ENDIF
      ENDDO
C--------------------------------------------------------------------------------
200   CONTINUE
C--------------------------------------------------------------------------------
C
C5------RETURN
      RETURN
      END
C--------------------------------------------------------------------------------
      SUBROUTINE ETSFLUXRD(KSTP,KPER)
C     ******************************************************************
C     READ ETS FLUX FOR RUNNING TRANSPORT FROM PREVIOUS FLOW RUN (ITRNSP=2)
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL ,      ONLY: IOUT,HNEW,IBOUND,BUFF,NCOL,NROW,NLAY,
     *                  NODES,IA,JA,NODLAY,IUNSTR,IVC,NEQS,FMBE
      USE GWFBASMODULE, ONLY: MSUM,VBNM,VBVL,PERTIM,TOTIM,DELT,ICBCFL
      USE GWFETSMODULE, ONLY: NETSOP,IETSCB,NETSEG,IETS,ETSR,ETSX,ETSS,
     1                        PXDP,PETM,NIETS,IETS,ETSF
C
      DOUBLE PRECISION RATOUT, QQ, HH, SS, DD, XX, HHCOF, RRHS,
     &                 PXDP1, PXDP2
      CHARACTER*16 TEXT
      INTEGER,ALLOCATABLE,DIMENSION(:,:) :: ITEMP
      INTEGER,ALLOCATABLE,DIMENSION(:) :: IBUFF
      DATA TEXT /'     ET SEGMENTS'/
C     ------------------------------------------------------------------
C
C1------INITIALIZE CELL-BY-CELL FLOW TERM FLAG (IBD)
      ZERO=0.
      IBD=0
      IF(IETSCB.GT.0) IBD=ICBCFL
C
      IF(IBD.NE.2) GO TO 200
C
C2------READ EVT FLUX ARRAY
      ALLOCATE(IBUFF(NIEVT))
      CALL UBDSV3Urd(KSTP,KPER,TEXT,IEVTCB,BUFF,IBUFF,
     1        NIEVT,NEVTOP,NODES,IOUT,DELT,PERTIM,TOTIM,IBOUND)
C3------PUT FLUX IN EVTF ARRAY
      DO NN=1,NIETS
        N = IETS(NN)
        ETSF(NN) = BUFF(N)
      ENDDO
C
      DEALLOCATE(IBUFF)
200   CONTINUE
C
C5------RETURN
      RETURN
      END
C--------------------------------------------------------------------------------
      SUBROUTINE DRTFLUXRD(KSTP,KPER)
C     ******************************************************************
C     READ DRT FLUX FOR RUNNING TRANSPORT FROM PREVIOUS FLOW RUN (ITRNSP=2)
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,  ONLY: IOUT,HNEW,IBOUND,BUFF,NCOL,NROW,NLAY,
     *                   NODES,NEQS,IUNSTR,FMBE,INCLN
      USE CLN1MODULE,  ONLY:NCLNNDS,ICLNCB
C
      USE GWFBASMODULE, ONLY: MSUM,VBNM,VBVL,PERTIM,TOTIM,DELT,ICBCFL,
     1                        IAUXSV
      USE GWFDRTMODULE, ONLY: DRTF,NDRTCL,MXDRT,IDRTCB,NDRTVL,IDRTFL,
     1                        NRFLOW,DRTAUX
C
      CHARACTER*16 TEXT
      DATA TEXT /'    DRAINS (DRT)'/
C     ------------------------------------------------------------------
C
C1------INITIALIZE CELL-BY-CELL FLOW TERM FLAG (IBD)
      ZERO=0.
      IBD=0
      IF (IDRTCB.GT.0) IBD=ICBCFL
C2-----STOP IF IBD IS 1 SINCE INFORMATION IS LOST ON MULTIPLE DRTs TO A CELL
      IF(IBD.EQ.1)THEN
        WRITE(IOUT,60)
60      FORMAT(1X,'*** IBD=1 LOSES INFORMATION UNLESS ONE ',
     1   'DRT BOUNDARY IS PROVIDED PER CELL, STOPPING')
        STOP
      ENDIF
      IF(IBD.EQ.0.AND.NDRTCL.GT.0)THEN
        WRITE(IOUT,61)
61      FORMAT(1X,'*** IBD = 0 SO WELL INFORMATION WAS NOT WRITTEN, '
     1   'STOPPING')
        STOP
      ENDIF
      IF(IBD.NE.2) GO TO 200
C
C3-----READ HEADER FOR COMPACT BUDGET
        NAUX=NDRTVL-5-2-IDRTFL
        IF (IAUXSV.EQ.0) NAUX = 0
         IICLNCB=0
         NNCLNNDS=0
         IF(INCLN.GT.0) THEN
           IICLNCB=ICLNCB
           NNCLNNDS=NCLNNDS
         ENDIF
         CALL UBDSVHDRTrd(IUNSTR,KSTP,KPER,IOUT,IDRTCB,IICLNCB,NODES,
     1    NNCLNNDS,NCOL,NROW,NLAY,NDRTCL,NDRTVL,NAUX,IBOUND,
     2    TEXT,DRTAUX,DELT,PERTIM,TOTIM,IDRTFL,DRTF)
C
C4-----READ EACH DRT FLOW
      IF (NDRTCL.LE.0) GOTO 200
      DO L=1,NDRTCL
        ND=DRTF(1,L)
        INR=DRTF(6,L)
        IF(IBOUND(ND).EQ.0) CYCLE
        Q=ZERO
        CALL UBDSVRECrd(IUNSTR,ND,NODES,NNCLNNDS,IDRTCB,IICLNCB,NDRTVL,
     1    8,NAUX,Q,DRTF(:,L),IBOUND,NCOL,NROW,NLAY)
CCSP        IF (IDRTFL.NE.0 .AND. ILR.GT.0) THEN
        IF (IDRTFL.NE.0) THEN
          CALL UBDSVRECrd(IUNSTR,INR,NODES,NNCLNNDS,IDRTCB,IICLNCB,
     1    NDRTVL,8,NAUX,QIN,DRTF(:,L),IBOUND,NCOL,NROW,NLAY)
          DRTF(NDRTVL-1,L) = QIN
        ENDIF
        DRTF(NDRTVL,L) = Q
      ENDDO
C--------------------------------------------------------------------------------
200   CONTINUE
C--------------------------------------------------------------------------------
C
C5------RETURN
      RETURN
      END
C--------------------------------------------------------------------------------
      SUBROUTINE QRTFLUXRD(KSTP,KPER)
C     ******************************************************************
C     READ QRT FLUX FOR RUNNING TRANSPORT FROM PREVIOUS FLOW RUN (ITRNSP=2)
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,  ONLY: IOUT,HNEW,IBOUND,BUFF,NCOL,NROW,NLAY,
     *                   NODES,NEQS,IUNSTR,TOP,BOT,AREA,FMBE,INCLN
      USE GWFBASMODULE, ONLY: MSUM,VBNM,VBVL,PERTIM,TOTIM,DELT,ICBCFL,
     1                        IAUXSV
      USE GWFQRTMODULE, ONLY: QRTF,NQRTCL,MXQRT,IQRTCB,NQRTVL,IQRTFL,
     1                  QRTAUX,QRTFLOW,IQRTQV,RTAREA,MXRTCELLS,NodQRT
      USE CLN1MODULE, ONLY: ACLNNDS,NCLNNDS,ICLNCB
C
      DOUBLE PRECISION HD,RATIN,RATOUT,QQ,QIN,QTHIK,X,Y,THCK,BOTT,FR
      CHARACTER*16 TEXT
      DATA TEXT /'    SINKS (QRT)'/
C     ------------------------------------------------------------------
C
C1------INITIALIZE CELL-BY-CELL FLOW TERM FLAG (IBD)
      ZERO=0.
      IBD=0
      IF (IQRTCB.GT.0) IBD=ICBCFL
C2-----STOP IF IBD IS 1 SINCE INFORMATION IS LOST ON MULTIPLE DRTs TO A CELL
      IF(IBD.EQ.1)THEN
        WRITE(IOUT,60)
60      FORMAT(1X,'*** IBD=1 LOSES INFORMATION UNLESS ONE ',
     1   'QRT BOUNDARY IS PROVIDED PER CELL, STOPPING')
        STOP
      ENDIF
      IF(IBD.EQ.0.AND.NQRTCL.GT.0)THEN
        WRITE(IOUT,61)
61      FORMAT(1X,'*** IBD = 0 SO WELL INFORMATION WAS NOT WRITTEN, '
     1   'STOPPING')
        STOP
      ENDIF
      IF(IBD.NE.2) GO TO 200
C
C3-----READ HEADER FOR COMPACT BUDGET
      NAUX=NQRTVL-5-IQRTFL                 !************************************
      IF (IAUXSV.EQ.0) NAUX = 0
      IICLNCB=0
      NNCLNNDS=0
      IF(INCLN.GT.0) THEN
        IICLNCB=ICLNCB
        NNCLNNDS=NCLNNDS
      ENDIF
      CALL UBDSVHDRQRTrd(IUNSTR,KSTP,KPER,IOUT,IQRTCB,IICLNCB,NODES,
     1  NNCLNNDS,NCOL,NROW,NLAY,NQRTCL,NQRTVL,NAUX,IBOUND,
     2  TEXT,QRTAUX,DELT,PERTIM,TOTIM,IQRTFL,NodQRT,MXRTCELLS,QRTF)
C
C4-----READ EACH DRT FLOW
      IF (NQRTCL.LE.0) GOTO 200
      DO L=1,NQRTCL
        ND=QRTF(1,L)
        IF(IBOUND(ND).EQ.0) CYCLE
        Q=ZERO
        CALL UBDSVRECrd(IUNSTR,ND,NODES,NNCLNNDS,IQRTCB,IICLNCB,NQRTVL,
     1    6,NAUX,Q,QRTF(:,L),IBOUND,NCOL,NROW,NLAY)
        QRTF(NQRTVL,L) = Q
        IF (NumRT.NE.0) THEN
          DO I=IRT-NUMRT+1, IRT
            CALL UBDSVRECrd(IUNSTR,I,NODES,NNCLNNDS,IQRTCB,IICLNCB,
     1    NQRTVL,6,NAUX,QRTFLOW(I),QRTF(:,L),IBOUND,NCOL,NROW,NLAY)
          ENDDO
        ENDIF
      ENDDO
C--------------------------------------------------------------------------------
200   CONTINUE
C--------------------------------------------------------------------------------
C
C5------RETURN
      RETURN
      END
C--------------------------------------------------------------------------------
      SUBROUTINE SFRFLUXRD(KSTP,KPER)
C     ******************************************************************
C     READ SFR FLUX FOR RUNNING TRANSPORT FROM PREVIOUS FLOW RUN (ITRNSP=2)
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GWFSFRMODULE
      USE GWFLAKMODULE, ONLY: VOL, LKARR1, STGNEW, STGOLD
      USE GLOBAL,       ONLY: NCOL, NROW, NLAY, IOUT, ISSFLG, IBOUND,
     +                        HNEW, BUFF, BOT, NODES, IVSD, NODLAY,
     +                        IA, JAS, IVC, JA, IUNSTR, NEQS, FMBE
      USE GWFBASMODULE, ONLY: MSUM, ICBCFL, IBUDFL, DELT, PERTIM, TOTIM,
     +                        VBVL, VBNM
      USE GWFBCFMODULE, ONLY: HDRY
      USE GWFRCHMODULE, ONLY: RECH  !cjm
      IMPLICIT REAL (A-H,O-Z)
      INTRINSIC FLOAT, ABS, IABS, DSQRT, DLOG10, SQRT, SNGL
      CHARACTER*16 text, strtxt, txtlst
      DATA text/'  STREAM LEAKAGE'/
      DATA strtxt/'STREAMFLOW OUT  '/
      DATA txtlst/'STREAM LISTING  '/
C-----------------------------------------------------------------------------
C1------INITIALIZE CELL-BY-CELL FLOW TERM FLAG (IBD)
      ZERO=0.
      IBD=0
      ibdst = 0
      IF ( ISTCB1.EQ.-1 .AND. ICBCFL.NE.0 ) THEN
        ibd = -1
        iblst = -1
      ELSE IF ( ISTCB1.EQ.-1 .AND. IBUDFL.GT.0 ) THEN
        ibd = -1
        iblst = -1
      ELSE IF ( ISTCB1.GT.0 .AND. ICBCFL.NE.0 ) THEN
        ibd = ICBCFL
        iout1 = ISTCB1
      END IF
      IF ( ISTCB2.GT.0 .AND. ICBCFL.NE.0 ) THEN
        ibdst = -1
        iout2 = ISTCB2
      ELSE IF ( ISTCB2.LT.0 .AND. ICBCFL.NE.0 ) THEN
        ibdst = ICBCFL
        iout2 = ABS(ISTCB2)
      END IF
C2-----STOP IF IBD IS 1 SINCE INFORMATION IS LOST ON MULTIPLE WELs TO A CELL
      IF(IBD.EQ.1.or.ibdst.eq.1)THEN
        WRITE(IOUT,60)
60      FORMAT(1X,'*** IBD=1 LOSES INFORMATION UNLESS ONE ',
     1   'SFR BOUNDARY IS PROVIDED PER CELL, STOPPING')
        STOP
      ENDIF
      IF(IBD.EQ.0.AND.NSTRM.GT.0)THEN
        WRITE(IOUT,61)
61      FORMAT(1X,'*** IBD = 0 SO WELL INFORMATION WAS NOT WRITTEN, '
     1   'STOPPING')
        STOP
      ENDIF
C
C3-----READ HEADER FOR COMPACT BUDGET
      IF(IBD.EQ.2) THEN
         IF(IUNSTR.EQ.0) THEN
CSP           CALL UBDSV2rd(Kkstp, Kkper, text, iout1, NCOL, NROW,
CSP     +                            NLAY, NSTRM, IOUT, DELT, PERTIM,
CSP     +                            TOTIM, IBOUND)
         ELSE
           CALL UBDSV2Urd(Kkstp, Kkper, text, iout1, NEQS,
     +                            NSTRM, IOUT, DELT, PERTIM,
     +                            TOTIM, IBOUND)
         ENDIF
      END IF
      IF ( ibdst.EQ.2 ) THEN
        IF(IUNSTR.EQ.0) THEN
CSP          CALL UBDSV2rd (Kkstp, Kkper, strtxt, iout2, NCOL,
CSP     +                               NROW, NLAY, NSTRM, IOUT, DELT,
CSP     +                               PERTIM, TOTIM, IBOUND)
        ELSE
          CALL UBDSV2Urd(Kkstp, Kkper, STRTXT, iout2, NEQS,
     +                            NSTRM, IOUT, DELT, PERTIM,
     +                            TOTIM, IBOUND)
        ENDIF
      ENDIF
C
C4-----READ EACH STR LEAKAGE FLOW
      IF(NSTRM.LE.0) GO TO 100
      DO L=1,NSTRM
        N=ISTRM(1,L)
        IF(IBOUND(N).EQ.0) CYCLE
        Q=ZERO
        IF ( imassroute.EQ.1 .AND. ibd.EQ.2 ) THEN
          IF(IUNSTR.EQ.0)THEN
CSP            CALL UBDSVArd(iout1, NCOL, NROW, ic, ir, il,
CSP     +                   Q, IBOUND, NLAY)
          ELSE
            CALL UBDSVAUrd(iout1, neqs, NCPT,
     +                   Q, IBOUND)
          ENDIF
        END IF
        STRM(11,L) = Q
      ENDDO
C--------------------------------------------------------------------------------
100   CONTINUE
C--------------------------------------------------------------------------------
      IF(IBDST.LE.2) GO TO 200
C5------READ FOR EACH STREAM OUTFLOW
      DO L=1,NSTRM
        NN = ISTRM(1, l)
        IF(IBOUND(NN).LE.0) CYCLE
           IF(IUNSTR.EQ.0) THEN
             IL = (ND-1) / (NCOL*NROW) + 1
             IJ = ND - (IL-1)*NCOL*NROW
             IR = (IJ-1)/NCOL + 1
             IC = IJ - (IR-1)*NCOL
             CALL UBDSVA(iout2, NCOL, NROW, ic, ir,
     +           il, STRM(9, l), IBOUND, NLAY)
          ELSE
             CALL UBDSVAU(iout2, NEQS, NN,
     +           STRM(9, l), IBOUND)
          ENDIF
      ENDDO
C--------------------------------------------------------------------------------
200   CONTINUE
C--------------------------------------------------------------------------------
C
C6------RETURN
      RETURN
      END
C--------------------------------------------------------------------------------
      SUBROUTINE LAKFLUXRD(KSTP,KPER)
C     ******************************************************************
C     READ LAK FLUX FOR RUNNING TRANSPORT FROM PREVIOUS FLOW RUN (ITRNSP=2)
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GWFLAKMODULE
      USE GLOBAL,      ONLY: NCOL, NROW, NLAY, NODES, IBOUND, IOUT,
     +                       ISSFLG, DELR, DELC, HNEW, NODLAY,IA,IUNSTR,
     +                       BUFF, AREA, TOP, JA, BOT, IVSD, FMBE, INCLN
      USE CLN1MODULE,  ONLY:NCLNNDS,ICLNCB
      USE GWFBASMODULE, ONLY: MSUM, ICBCFL, IAUXSV, DELT, PERTIM, TOTIM,
     +                        HNOFLO, VBVL, VBNM
      USE GWFSFRMODULE, ONLY: STRIN, DLKSTAGE, SLKOTFLW
!      USE GWFUZFMODULE, ONLY: SURFDEP,IUZFBND,FINF,VKS
      IMPLICIT NONE
      CHARACTER*16 TEXT
C
      DOUBLE PRECISION BOTLK,BOTCL,CONDUC,H,FLOBOT,STGON,
     1RATE,RATIN,RATOUT
      DOUBLE PRECISION THET1,SURFDPTH,CONDMX,BOTLKUP,BOTLKDN,VOL2
      DOUBLE PRECISION FLOTOUZF,SILLELEV,ADJSTAGE, voltest, areatest
      DOUBLE PRECISION RAMPGW,RAMPSTGO,RAMPSTGN,RAMPSTGON,HTEMP
      DOUBLE PRECISION CLOSEZERO, FLOBO2, FLOBO3, RUNOFF, DLSTG
      DOUBLE PRECISION RUNF, RUNFD, AREA1, RAIN
      REAL zero, FACE, R, WDRAW, OLDSTAGE, AVHD, TOTARE, SUM
      REAL STGTST, SVT1, TVOLM, STO, FLSUM, TV, PPTIN, EOUT
      REAL SEEPUZF, QIN, QOUT, QSIN, QSOUT, DENOM
      REAL TOTIMREAL,Q,DRAI
      INTEGER L1, IGRID, ISS, KPER, IBD, KCNT, LDR, NAUX, KSTP, IL
      INTEGER IR, IC, LK, ITRIB, INODE, LAKE, IUNITUZF, II, L
      INTEGER IUNITGWT, LL, K, J, I, JCLS, ICL, IC4, ICM4
      INTEGER ITYPE, ITYPE2, INOFLO, IC5, ICM, IS1, IS2, ICl1, LK3, LK1
      INTEGER ICM2, IC2, IC3, ICNR, ICM3, ICNT, ICM1, L11, ICNR1
      INTEGER ILB, IRB, ICB, LDR1, NN, IIC, JIC, LIC, IUNITGAGE
      INTEGER NSOL, ILL, IL2, IUNITSFR, IL1, NL, NT, IVERT, NL1
      INTEGER N, ND, IALY, IJ, IROW, JCOLMN, NLL, NL2, ILAY
      INTEGER IICLNCB, NNCLNNDS,IGWTOUT
      DIMENSION JCLS(NCLS,ICMX)
      DIMENSION ILB(5),IRB(5),ICB(5)
      CHARACTER*16 LAKAUX(20)
      DIMENSION FACE(LKNODE)
      DATA TEXT /'   LAKE  SEEPAGE'/
      DATA LAKAUX(1)/'IFACE'/
C     ------------------------------------------------------------------
C
C1------INITIALIZE CELL-BY-CELL FLOW TERM FLAG (IBD)
      ZERO=0.
      IBD=0
      IF(ILKCB.GT.0) IBD=ICBCFL
C2-----STOP IF IBD IS 1 SINCE INFORMATION IS LOST ON MULTIPLE DRNs TO A CELL
      IF(IBD.EQ.1)THEN
        WRITE(IOUT,60)
60      FORMAT(1X,'*** IBD=1 LOSES INFORMATION UNLESS ONE ',
     1   'LAKE BOUNDARY IS PROVIDED PER CELL, STOPPING')
        STOP
      ENDIF
      IF(IBD.EQ.0.AND.NLAKES.GT.0)THEN
        WRITE(IOUT,61)
61      FORMAT(1X,'*** IBD = 0 SO WELL INFORMATION WAS NOT WRITTEN, '
     1   'STOPPING')
        STOP
      ENDIF
      IF(IBD.NE.2) GO TO 200
C
C3-----READ HEADER FOR COMPACT BUDGET
      NAUX=1
      IICLNCB=0
      NNCLNNDS=0
      IF(INCLN.GT.0) THEN
        IICLNCB=ICLNCB
        NNCLNNDS=NCLNNDS
      ENDIF
      DO N=1, LKNODE
        FACE(N)=ILAKE(6,N)
      ENDDO
      IF(IGWTOUT.EQ.0)
     1 CALL UBDSVHDRrd(IUNSTR,KSTP,KPER,IOUT,ILKCB,IICLNCB,NODES,
     1  NNCLNNDS,NCOL,NROW,NLAY,LKNODE,1,NAUX,IBOUND,
     2  TEXT,LAKAUX,DELT,PERTIM,TOTIM,FACE)
C
C4-----READ EACH LAK FLOW
      IF(NLAKES.LE.0) GO TO 200
      DO II = 1,2
      DO L=1,LKNODE
        IF(II.EQ.2)THEN
          FACE(1)=ILAKE(6,L)
          R=RATE
          IF(IGWTOUT.EQ.0)
     1    CALL UBDSVRECrd(IUNSTR,NL,NODES,NNCLNNDS,ILKCB,IICLNCB,
     1      1,1,NAUX,Q,FACE(1),IBOUND,NCOL,NROW,NLAY)
          FLOB(LK) = Q
        ENDIF
      ENDDO
      ENDDO
C--------------------------------------------------------------------------------
200   CONTINUE
C--------------------------------------------------------------------------------
C
C5------RETURN
      RETURN
      END
C--------------------------------------------------------------------------------


C------------------------------------------------------------------------
      SUBROUTINE GWT2STO1BDW(KSTP,KPER,ICOMP,ISS)
C     ******************************************************************
C     CALCULATE MASS BUDGET TERMS FOR ALL TRANSPORT CELLS
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,   ONLY:IOUT,NCOL,NROW,NLAY,IBOUND,BUFF,NODES,IUNSTR,
     1  AMAT,IA,JA,TOP,BOT,AREA,Sn,So,NEQS,INCLN,IUNIT
      USE CLN1MODULE, ONLY: ACLNNDS,NCLNNDS,ICLNMB
      USE GWFBASMODULE,ONLY:MSUM,ICBCFL,IAUXSV,DELT,PERTIM,TOTIM,DELT
      USE GWTBCTMODULE, ONLY: ICBUND,CONC,MSUMT,VBVLT,VBNMT,IBCTCB,
     1 IADSORB,ADSORB,FLICH,PRSITY,CONCO,ICT,MCOMPT
C
      CHARACTER*16 TEXT(2)
      DOUBLE PRECISION RATIN,RATOUT,QQ,VODT,ADSTERM,FL,CW,CWO,ALENG,
     *  DTERMS,RTERMS,VOLU,RATINTVM,RATOUTTVM
      DATA TEXT(1) /'    MASS STORAGE'/
      DATA TEXT(2) /'CLN MASS STORAGE'/
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
      Q  = 0.0
      IF(ICT.EQ.0)THEN
C-----------------------------------------------------
C6-------STORAGE TERM IN WATER
        DTERMS = 0.0
        RTERMS = 0.0
        CALL GWT2BCT1STOW(N,ICOMP,DTERMS,RTERMS,VODT,VOLU,ALENG,ISS)
        QQ = QQ - DTERMS * CONC(N,ICOMP) + RTERMS
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
C-----------------------------------------------------------------------
C13---------ADJUST STORAGE TERMS FOR TVM PACKAGE
      IF(IUNIT(64).GT.0) THEN
        CALL GWT2STO1BDTVM(KSTP,KPER,ICOMP,ISS,RATINTVM,RATOUTTVM)
        RATIN=RATIN - RATOUTTVM
        RATOUT=RATOUT - RATINTVM
      ENDIF
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
        IF(INCLN.GT.0.AND.ICLNMB.NE.0)THEN
           CALL UBUDSVU(KSTP,KPER,TEXT(2),ICLNMB,BUFF(NODES+1:NEQS),
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
C-----------------------------------------------------------------------
      SUBROUTINE GWT2STO1BDS(KSTP,KPER,ICOMP,ISS)
C     ******************************************************************
C     CALCULATE ADSORBED MASS BUDGET TERMS FOR ALL TRANSPORT CELLS
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,   ONLY:IOUT,NCOL,NROW,NLAY,IBOUND,BUFF,NODES,IUNSTR,
     1  AMAT,IA,JA,TOP,BOT,AREA,Sn,So,NEQS,INCLN,IUNIT
      USE CLN1MODULE, ONLY: ACLNNDS,NCLNNDS
      USE GWFBASMODULE,ONLY:MSUM,ICBCFL,IAUXSV,DELT,PERTIM,TOTIM,DELT
      USE GWTBCTMODULE, ONLY: ICBUND,CONC,MSUMT,VBVLT,VBNMT,IBCTCB,
     1 IADSORB,ADSORB,FLICH,PRSITY,CONCO,ICT,MCOMPT
C
      CHARACTER*16 TEXT
      DOUBLE PRECISION RATIN,RATOUT,QQ,VODT,ADSTERM,FL,CW,CWO,ALENG,
     *  DTERMS,RTERMS,VOLU,RATINTVM,RATOUTTVM
      DATA TEXT /'   ADSORBED MASS'/ 
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
      Q  = 0.0
      IF(ICT.EQ.0)THEN
C5-------STORAGE TERM ON SOIL - SKIP IF IMMOBILE COMPONENTS
        IF(ICOMP.GT.MCOMPT) GO TO 252
        IF(IADSORB.EQ.1)THEN
          ADSTERM = ADSORB(N,ICOMP) * VODT
          QQ = ADSTERM * (CONC(N,ICOMP) - CONCO(N,ICOMP))
        ELSEIF(IADSORB.EQ.2)THEN
          ADSTERM = ADSORB(N,ICOMP) * VODT
          FL = FLICH(N,ICOMP)
          CW = 0.0
          CWO = 0.0
          IF(CONC(N,ICOMP).GT.0.0) CW = CONC(N,ICOMP)
          IF(CONCO(N,ICOMP).GT.0.0) CWO = CONCO(N,ICOMP)
          QQ = ADSTERM * (CW**FL - CWO**FL)
        ELSEIF(IADSORB.EQ.3)THEN
          ADSTERM = ADSORB(N,ICOMP) * VODT
          FL = FLICH(N,ICOMP)
          CW = 0.0
          CWO = 0.0
          IF(CONC(N,ICOMP).GT.0.0) CW = CONC(N,ICOMP)
          IF(CONCO(N,ICOMP).GT.0.0) CWO = CONCO(N,ICOMP)
          QQ = ADSTERM * (CW/(1.0+FL*CW) - CWO/(1.0+FL*CWO))          
        ENDIF
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
         IF(IBDLBL.EQ.0) WRITE(IOUT,61) TEXT,KPER,KSTP
   61    FORMAT(1X,/1X,A,'   PERIOD ',I4,'   STEP ',I3)
        IF(IUNSTR.EQ.0)THEN
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
C-----------------------------------------------------------------------
C13------IF CELL-BY-CELL FLOWS WILL BE SAVED AS A 3-D ARRAY,
C13------CALL UBUDSV TO SAVE THEM.
      IF(IBD.GE.1)THEN
        IF(IUNSTR.EQ.0)THEN
          CALL UBUDSV(KSTP,KPER,TEXT,IBCTCB,BUFF,NCOL,NROW,
     1                          NLAY,IOUT)
        ELSE
          CALL UBUDSVU(KSTP,KPER,TEXT,IBCTCB,BUFF,NODES,
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
      VBNMT(MSUMT,ICOMP)=TEXT
C
C15------INCREMENT BUDGET TERM COUNTER(MSUM).
      MSUMT=MSUMT+1
C
C16------RETURN
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE GWT2STO1BDTVM(KSTP,KPER,ICOMP,ISS,RATIN,RATOUT)
C     ******************************************************************
C     CALCULATE MASS BUDGET TERMS FOR TRANSIENT TRANSPORT PROPERTIES
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,   ONLY:IOUT,NCOL,NROW,NLAY,IBOUND,BUFF,NODES,IUNSTR,
     1  AMAT,IA,JA,TOP,BOT,AREA,Sn,So,NEQS,INCLN
      USE CLN1MODULE, ONLY: ACLNNDS,NCLNNDS,ICLNMB
      USE GWFBASMODULE,ONLY:MSUM,ICBCFL,IAUXSV,DELT,PERTIM,TOTIM,DELT
      USE GWTBCTMODULE, ONLY: ICBUND,CONC,MSUMT,VBVLT,VBNMT,IBCTCB,
     1 IADSORB,ADSORB,FLICH,PRSITY,CONCO,ICT
      USE TVMU2MODULE
C
      CHARACTER*16 TEXT(2)
      DOUBLE PRECISION RATIN,RATOUT,QQ,VODT,ADSTERM,FL,CW,CWO,ALENG,
     *  DTERMS,RTERMS,VOLU
      DATA TEXT(1) /'TVM MASS STORAGE'/
      DATA TEXT(2) /'CLN TVM MASS STO'/
C     ------------------------------------------------------------------
C
C1------CLEAR RATIN AND RATOUT ACCUMULATORS, AND SET CELL-BY-CELL
C1------BUDGET FLAG.
      ZERO=0.
      RATIN=ZERO
      RATOUT=ZERO
      IBD=0
      IF(IBCTCB.LT.0 .AND. ICBCFL.NE.0) IBD=-1
      IF(IBCTCB.GT.0) IBD=ICBCFL
      IBDLBL=0
CCCC
CCCC2------CLEAR THE BUFFER.
CCC      DO 50 N=1,NEQS
CCC      BUFF(N)=ZERO
CCC50    CONTINUE
C
C3------LOOP THROUGH EACH NODE AND CALCULATE STORAGE
CCC      DO 100 N=1,NEQS
      DO 100 INODE = 1, NTVMPOR
      N=ITVMPOR(INODE)
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
      Q = 0.0
      IF(ICT.EQ.0)THEN
C-----------------------------------------------------
C6-------STORAGE TERM IN WATER
        DTERMS = 0.0
        RTERMS = 0.0
        CALL GWT2BCT1STOWTVM(N,ICOMP,DTERMS,RTERMS,VODT,VOLU,ALENG,ISS)
        QQ = QQ + RTERMS
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
CCCC
CCCC13------IF CELL-BY-CELL FLOWS WILL BE SAVED AS A 3-D ARRAY,
CCCC13------CALL UBUDSV TO SAVE THEM.
CCC      IF(IBD.GE.1)THEN
CCC        IF(IUNSTR.EQ.0)THEN
CCC          CALL UBUDSV(KSTP,KPER,TEXT(1),IBCTCB,BUFF,NCOL,NROW,
CCC     1                          NLAY,IOUT)
CCC        ELSE
CCC          CALL UBUDSVU(KSTP,KPER,TEXT(1),IBCTCB,BUFF,NODES,
CCC     1                          IOUT,PERTIM,TOTIM)
CCC        ENDIF
CCC        IF(INCLN.GT.0)THEN
CCC           CALL UBUDSVU(KSTP,KPER,TEXT(2),ICLNMB,BUFF(NODES+1:NEQS),
CCC     1                 NCLNNDS,IOUT,PERTIM,TOTIM)
CCC        ENDIF
CCC      ENDIF
CCCC
CCCC14------MOVE RATES, VOLUMES & LABELS INTO ARRAYS FOR PRINTING.
CCC  200 RIN=RATIN
CCC      ROUT=RATOUT
CCC      VBVLT(3,MSUMT,ICOMP)=RIN
CCC      VBVLT(4,MSUMT,ICOMP)=ROUT
CCC      VBVLT(1,MSUMT,ICOMP)=VBVLT(1,MSUMT,ICOMP)+RATIN*DELT
CCC      VBVLT(2,MSUMT,ICOMP)=VBVLT(2,MSUMT,ICOMP)+RATOUT*DELT
CCC      VBNMT(MSUMT,ICOMP)=TEXT(1)
CCCC
CCCC15------INCREMENT BUDGET TERM COUNTER(MSUM).
CCC      MSUMT=MSUMT+1
CCCC
C16------RETURN
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE GWT2DCY1BD(KSTP,KPER,ICOMP)
C     ******************************************************************
C     CALCULATE MASS DECAY TERMS FOR ALL TRANSPORT CELLS
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,   ONLY:IOUT,NCOL,NROW,NLAY,IBOUND,BUFF,NODES,IUNSTR,
     1  AMAT,IA,JA,TOP,BOT,AREA,Sn,So,NEQS,INCLN
      USE CLN1MODULE, ONLY: ACLNNDS,NCLNNDS,ICLNMB
      USE GWFBASMODULE,ONLY:MSUM,ICBCFL,IAUXSV,DELT,PERTIM,TOTIM,DELT
      USE GWTBCTMODULE, ONLY: ICBUND,CONC,MSUMT,VBVLT,VBNMT,IBCTCB,
     1 IADSORB,ADSORB,FLICH,PRSITY,CONCO,ICT,IZOD,IFOD,ZODRW,FODRW,
     1  ZODRS,FODRS,IAW_ADSORB
C
      CHARACTER*16 TEXT(2)
      DOUBLE PRECISION RATIN,RATOUT,QQ,VOLU,ADSTERM,FL,CW,CWO,ALENG,X,Y,
     1  CEPS,EPS,CT
      DATA TEXT(1) /'      MASS DECAY'/
      DATA TEXT(2) /'  CLN MASS DECAY'/
C     ------------------------------------------------------------------
C1------RETURN IF NO DECAY IN SIMULATION
      IF(IZOD.EQ.0.AND.IFOD.EQ.0) RETURN
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
      DO 100 N=1,NEQS
C
C5-----IF THE CELL IS NOT PCB OR WRONG COMPONENT SPECIES, IGNORE IT.
      IF(ICBUND(N).EQ.0)GO TO 99
C
      IF(N.LE.NODES)THEN
        ALENG = TOP(N) - BOT(N)
      ELSE
        ALENG = ACLNNDS(N-NODES,4)
      ENDIF
      VOLU = AREA(N) * ALENG
      QQ = 0.0
      Q = 0.0
C-----------------------------------------------------------------------------
      IF(ICT.EQ.0)THEN  !----------WATER PHASE CONCENTRATION FORMULATION
C-----------------------------------------------------------------------------
C6-------DECAY TERM ON SOIL (NO ADSORPTION ON CLN)
        IF(N.LE.NODES)THEN
C
C7---------ZERO ORDER DECAY ON SOIL
          IF(IZOD.GE.2.AND.IADSORB.GT.0)THEN
            CT = - VOLU * ZODRS(N,ICOMP)
            EPS = 0.01
            CEPS = MAX(0.0,CONC(N,ICOMP))
            X = CEPS /EPS
            CALL SMOOTH(X,Y)
            QQ =  CT * Y
          ENDIF
C
C8---------FIRST ORDER DECAY ON SOIL
          IF(IFOD.GE.2.AND.IADSORB.GT.0)THEN
            CT = -ADSORB(N,ICOMP) * VOLU * FODRS(N,ICOMP)
            IF(IADSORB.EQ.1)THEN
C9--------------FOR LINEAR ADSORPTION
              QQ = QQ + CT * CONC(N,ICOMP)
            ELSEIF(IADSORB.EQ.2) THEN
C10--------------FOR NON-LINEAR FREUNDLICH ADSORPTION
              ETA = FLICH(N,ICOMP)
              QQ =  CT * CONC(N,ICOMP) ** ETA
            ELSEIF(IADSORB.EQ.3) THEN
C10--------------FOR NON-LINEAR LANGMUIR ADSORPTION 
              ETA = FLICH(N,ICOMP)
              QQ =  CT * CONC(N,ICOMP) / (1.0 + ETA*CONC(N,ICOMP))      
            ENDIF
          ENDIF
C --------DECAY ON AIR-WATER INTERFACE ADSORBED SOLUTES            
            IF(IAW_ADSORB.NE.0)
     *     CALL GWT2DCYAW1BD(N,ICOMP,QQ,VOLU)            
        ENDIF
C-----------------------------------------------------------------------------
C11-------DECAY TERM IN WATER
C-----------------------------------------------------------------------------
C12-------ZERO ORDER DECAY IN WATER
        IF(IZOD.EQ.1.OR.IZOD.EQ.3)THEN
          CT = -Sn(N)* VOLU * ZODRW(N,ICOMP)
          EPS = 0.01
          CEPS = MAX(0.0,CONC(N,ICOMP))
          X = CEPS /EPS
          CALL SMOOTH(X,Y)
          QQ =  QQ + CT * Y
        ENDIF
C
C13-------FIRST ORDER DECAY IN WATER
        IF(IFOD.EQ.1.OR.IFOD.EQ.3)THEN
          CT =  -Sn(N)* VOLU * FODRW(N,ICOMP)
          QQ = QQ + CT * CONC(N,ICOMP)
        ENDIF
      ELSE
CSP FINISH TOTAL CONCENTRATION FORMULATION
      ENDIF
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
     1       '  DECAY ',1PG15.6)
        ELSE
           WRITE(IOUT,63) L,N,Q
   63    FORMAT(1X,'CBC  ',I6,'    NODE ',I8,'  DECAY ',1PG15.6)
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
           CALL UBUDSVU(KSTP,KPER,TEXT(2),ICLNMB,BUFF(NODES+1:NEQS),
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
      SUBROUTINE GWT2GEN1BD(KSTP,KPER,ICOMP)
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
      USE CLN1MODULE, ONLY: ACLNNDS,NCLNNDS,ICLNMB
      USE GWFBASMODULE,ONLY:MSUM,ICBCFL,IAUXSV,DELT,PERTIM,TOTIM,DELT
      USE GWTBCTMODULE, ONLY: ICBUND,CONC,MSUMT,VBVLT,VBNMT,IBCTCB,
     1 IADSORB,ADSORB,FLICH,PRSITY,CONCO,ICT,IZOD,IFOD,ZODRW,FODRW,
     1  ZODRS,FODRS,ICHAIN,NPARENT,JPARENT,STOTIO,ISPRCT,SPTLRCT,
     1  IAW_ADSORB
C
      CHARACTER*16 TEXT(2)
      DOUBLE PRECISION RATIN,RATOUT,QQ,VOLU,ADSTERM,FL,CW,CWO,ALENG,X,Y,
     1  CEPS,EPS,CT
      DATA TEXT(1) /'        MASS GEN'/
      DATA TEXT(2) /'    CLN MASS GEN'/
C     ------------------------------------------------------------------
C1------RETURN IF NO DECAY IN SIMULATION
      IF(IZOD.EQ.0.AND.IFOD.EQ.0) RETURN
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
      DO 100 N=1,NEQS
C
C5-----IF THE CELL IS NOT PCB OR WRONG COMPONENT SPECIES, IGNORE IT.
      IF(ICBUND(N).EQ.0)GO TO 99
C
      IF(N.LE.NODES)THEN
        ALENG = TOP(N) - BOT(N)
      ELSE
        ALENG = ACLNNDS(N-NODES,4)
      ENDIF
      VOLU = AREA(N) * ALENG
      QQ = 0.0
      Q = 0.0
C----------------------------------------------------------------------
C ------ DO FOR ALL PARENTS OF COMPONENT ICOMP
      IF(NPARENT(ICOMP).EQ.0) GO TO 111
      DO NPAREN = 1, NPARENT(ICOMP)
        DTERMD = 0.0
        RTERMD = 0.0
        IPAREN = JPARENT(ICOMP,NPAREN)
        IF(ISPRCT.EQ.0) THEN
          STOIT  = STOTIO(ICOMP,NPAREN)
        ELSE
          STOIT = SPTLRCT(N,ICOMP,NPAREN)
        ENDIF
C-----------------------------------------------------------------------------
      IF(ICT.EQ.0)THEN  !----------WATER PHASE CONCENTRATION FORMULATION
C-----------------------------------------------------------------------------
C6-------GENERATION TERM ON SOIL (NO ADSORPTION ON CLN)
        IF(N.LE.NODES)THEN
C
C7---------ZERO ORDER GENERATION ON SOIL
          IF(IZOD.GE.2.AND.IADSORB.GT.0)THEN
            CT = VOLU * ZODRS(N,IPAREN) * STOIT
            EPS = 0.01
            CEPS = MAX(0.0,CONC(N,IPAREN))
            X = CEPS /EPS
            CALL SMOOTH(X,Y)
            QQ =  QQ + CT * Y
          ENDIF
C
C8---------FIRST ORDER GENERATION ON SOIL
          IF(IFOD.GE.2.AND.IADSORB.GT.0)THEN
            CT = ADSORB(N,IPAREN) * VOLU * FODRS(N,IPAREN) * STOIT
            IF(IADSORB.EQ.1)THEN
C9--------------FOR LINEAR ADSORPTION
              QQ = QQ + CT * CONC(N,IPAREN)
            ELSEIF(IADSORB.EQ.2)THEN
C10--------------FOR NON-LINEAR FREUNDLICH ADSORPTION 
              ETA = FLICH(N,IPAREN)
              QQ =  QQ + CT * CONC(N,IPAREN) ** ETA
            ELSEIF(IADSORB.EQ.3)THEN
C10--------------FOR NON-LINEAR LANGMUIR ADSORPTION 
              ETA = FLICH(N,IPAREN)
              QQ =  QQ + CT * CONC(N,IPAREN) /(1.0 + ETA*CONC(N,IPAREN))
            ENDIF
          ENDIF
C --------GENERATION ON AIR-WATER INTERFACE ADSORBED SOLUTES            
          IF(IAW_ADSORB.NE.0)
     *     CALL GWT2GENAW1BD(N,ICOMP,QQ,VOLU,NPAREN)   
        ENDIF
C-----------------------------------------------------------------------------
C11-------GENERATION TERM IN WATER
C-----------------------------------------------------------------------------
C12-------ZERO ORDER GENERATION IN WATER
        IF(IZOD.EQ.1.OR.IZOD.EQ.3)THEN
          CT = Sn(N)* VOLU * ZODRW(N,IPAREN) * STOIT
          EPS = 0.01
          CEPS = MAX(0.0,CONC(N,IPAREN))
          X = CEPS /EPS
          CALL SMOOTH(X,Y)
          QQ =  QQ + CT * Y
        ENDIF
C
C13-------FIRST ORDER GENERATION IN WATER
        IF(IFOD.EQ.1.OR.IFOD.EQ.3)THEN
          CT =  Sn(N)* VOLU * FODRW(N,IPAREN) * STOIT
          QQ = QQ + CT * CONC(N,IPAREN)
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
           CALL UBUDSVU(KSTP,KPER,TEXT(2),ICLNMB,BUFF(NODES+1:NEQS),
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
      SUBROUTINE GWT2FMBE1BD(KSTP,KPER,ICOMP)
C     ******************************************************************
C     CALCULATE MASS BUDGET TERMS FOR ALL TRANSPORT CELLS
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,   ONLY:IOUT,NCOL,NROW,NLAY,IBOUND,BUFF,NODES,IUNSTR,
     1  NEQS,INCLN,IFMBC,FMBE
      USE CLN1MODULE, ONLY: ACLNNDS,NCLNNDS
      USE GWFBASMODULE,ONLY:MSUM,ICBCFL,IAUXSV,DELT,PERTIM,TOTIM,DELT
      USE GWTBCTMODULE, ONLY: ICBUND,CONC,MSUMT,VBVLT,VBNMT,CONCO
C
      CHARACTER*16 TEXT(2)
      DOUBLE PRECISION RATIN,RATOUT,QQ,VODT
      DATA TEXT(1) /'   TRNSP FMB ERR'/
      DATA TEXT(2) /'TRNSP CLN FMBERR'/
C     ------------------------------------------------------------------
C
C1------CLEAR RATIN AND RATOUT ACCUMULATORS, AND SET CELL-BY-CELL
C1------BUDGET FLAG.
      ZERO=0.
      RATIN=ZERO
      RATOUT=ZERO
      IBD=0
      IF(IFMBC.LT.0 .AND. ICBCFL.NE.0) IBD=-1
      IF(IFMBC.GT.0) IBD=ICBCFL
      IBDLBL=0
C
C2------CLEAR THE BUFFER.
      DO 50 N=1,NEQS
      BUFF(N)=ZERO
50    CONTINUE
C
C3------LOOP THROUGH EACH NODE AND CALCULATE ERROR
      amaxerr = 0.0
      nmax = 0
      DO 100 N=1,NEQS
C
C4-----IF THE CELL IS NOT PCB OR WRONG COMPONENT SPECIES, IGNORE IT.
      IF(ICBUND(N).EQ.0)GO TO 99
C
      QQ = -FMBE(N) * CONC(N,ICOMP)
      Q = QQ
      aerr = abs(qq )
      if(aerr.gt.amaxerr) then
          amaxerr = aerr
          nmax = n
      endif
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
     1       '  ERROR ',1PG15.6)
        ELSE
           WRITE(IOUT,63) N,Q
   63    FORMAT(1X,'    NODE ',I8,'ERROR ',1PG15.6)
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
100      CONTINUE
       write(iout,*)'max transport fmbe and n are',amaxerr,nmax
C
C13------IF CELL-BY-CELL FLOWS WILL BE SAVED AS A 3-D ARRAY,
C13------CALL UBUDSV TO SAVE THEM.
      IF(IBD.GE.1)THEN
        IF(IUNSTR.EQ.0)THEN
          CALL SGWF2BAS7TE(KSTP,KPER,IPFLG,ISA)
        ELSE
          CALL SGWF2BAS7TEU(KSTP,KPER,IPFLG,ISA)
        ENDIF
        IF(INCLN.GT.0)THEN
          CALL SCLN1TE(KSTP,KPER,IPFLG,ISA)
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
      SUBROUTINE SGWF2BAS7TEU(KSTP,KPER,IPFLG,ISA)
C     ******************************************************************
C     PRINT AND RECORD TRANSPORT MASS BALANCE ERROR FOR UNSTRUCTURED GWF GRID
C     RESULTING FROM FLOW BALANCE ERROR
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,      ONLY:NCOL,NROW,NLAY,IXSEC,FMBE,NODLAY,
     1                      IBOUND,IOUT,NODES,BUFF,
     1   IFMBC,MBEGWUNF,MBEGWUNT,MBECLNUNF,MBECLNUNT
      USE GWFBASMODULE,ONLY:PERTIM,TOTIM,IHEDFM,IHEDUN,LBHDSV,
     2                      CHEDFM,IOFLG
C
      CHARACTER*16 TEXT
      DATA TEXT /'    TRNS BAL ERR'/
C     ------------------------------------------------------------------
C
C4------FOR EACH LAYER: DETERMINE IF FMBE SHOULD BE PRINTED.
C4------IF SO THEN CALL ULAPRU TO PRINT FMBE.
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
C4A-----PRINT FMBE FOR CROSS SECTION.
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
C5------FOR EACH LAYER: DETERMINE IF FMBE SHOULD BE SAVED ON DISK.
C5------IF SO THEN CALL ULASAV OR ULASV2 TO SAVE FMBE.
      IFIRST=1
      IF(MBEGWUNT.LE.0) GO TO 80
      IF(IXSEC.EQ.0) THEN
        DO 79 K=1,NLAY
        KK=K
        IF(IOFLG(K,3).EQ.0) GO TO 79
        NNDLAY = NODLAY(K)
        NSTRT = NODLAY(K-1)+1
        IF(IFIRST.EQ.1) WRITE(IOUT,74) MBEGWUNT,KSTP,KPER
   74   FORMAT(1X,/1X,'FMBE WILL BE SAVED ON UNIT ',I8,
     1      ' AT END OF TIME STEP ',I8,', STRESS PERIOD ',I8)
        IFIRST=0
        IF(CHEDFM.EQ.' ') THEN
           CALL ULASAVU(BUFF,TEXT,KSTP,KPER,PERTIM,TOTIM,NSTRT,
     1                NNDLAY,KK,MBEGWUNT,NODES)
        ELSE
           CALL ULASV2U(BUFF,TEXT,KSTP,KPER,PERTIM,TOTIM,NSTRT,
     1             NNDLAY,KK,MBEGWUNT,CHEDFM,LBHDSV,IBOUND(NSTRT),NODES)
        END IF
c        IPFLG=1
   79   CONTINUE
C
C5A-----SAVE FMBE FOR CROSS SECTION.
      ELSE
        IF(IOFLG(1,3).NE.0) THEN
          WRITE(IOUT,74) MBEGWUNT,KSTP,KPER
          IF(CHEDFM.EQ.' ') THEN
             CALL ULASAVU(BUFF,TEXT,KSTP,KPER,PERTIM,TOTIM,NSTRT,
     1                NNDLAY,-1,MBEGWUNT,NODES)
          ELSE
             CALL ULASV2U(BUFF,TEXT,KSTP,KPER,PERTIM,TOTIM,NSTRT,
     1                  NNDLAY,-1,MBEGWUNT,CHEDFM,LBHDSV,IBOUND,NODES)
          END IF
c          IPFLG=1
        END IF
      END IF
C
C6------RETURN.
   80 CONTINUE
      RETURN
C
      END
      SUBROUTINE SGWF2BAS7TE(KSTP,KPER,IPFLG,ISA)
C     ******************************************************************
C     PRINT AND RECORD TRANSPORT MASS BALANCE ERROR FOR STRUCTURED GWF GRID
C     RESULTING FROM FLOW BALANCE ERROR
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,      ONLY:NCOL,NROW,NLAY,IXSEC,FMBE,NODLAY,
     1                      IBOUND,IOUT,
     1   IFMBC,MBEGWUNF,MBEGWUNT,MBECLNUNF,MBECLNUNT
      USE GWFBASMODULE,ONLY:PERTIM,TOTIM,IHEDFM,IHEDUN,LBHDSV,
     2                      CHEDFM,IOFLG
C
      REAL,          SAVE,    DIMENSION(:,:,:),    ALLOCATABLE ::BUFF
      CHARACTER*16 TEXT
      DATA TEXT /'    TRNS BAL ERR'/
C     ------------------------------------------------------------------
      ALLOCATE(BUFF(NCOL,NROW,NLAY))
C
C4------FOR EACH LAYER: DETERMINE IF FMBE SHOULD BE PRINTED.
C4------IF SO THEN CALL ULAPRS OR ULAPRW TO PRINT FMBE.
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
C4A-----PRINT FMBE FOR CROSS SECTION.
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
C5------FOR EACH LAYER: DETERMINE IF FMBE SHOULD BE SAVED ON DISK.
C5------IF SO THEN CALL ULASAV OR ULASV2 TO SAVE FMBE.
      IFIRST=1
      IF(MBEGWUNT.LE.0) GO TO 80
      IF(IXSEC.EQ.0) THEN
        DO 79 K=1,NLAY
        NSTRT = NODLAY(K-1)+1
        KK=K
        IF(IOFLG(K,3).EQ.0) GO TO 79
        IF(IFIRST.EQ.1) WRITE(IOUT,74) MBEGWUNT,KSTP,KPER
   74   FORMAT(1X,/1X,'FMBE WILL BE SAVED ON UNIT ',I4,
     1      ' AT END OF TIME STEP ',I4,', STRESS PERIOD ',I4)
        IFIRST=0
        IF(CHEDFM.EQ.' ') THEN
           CALL ULASAV(BUFF(1,1,K),TEXT,KSTP,KPER,PERTIM,TOTIM,NCOL,
     1                NROW,KK,MBEGWUNT)
        ELSE
           CALL ULASV2(BUFF(1,1,K),TEXT,KSTP,KPER,PERTIM,TOTIM,NCOL,
     1                NROW,KK,MBEGWUNT,CHEDFM,LBHDSV,IBOUND(NSTRT))
        END IF
   79   CONTINUE
C
C5A-----SAVE FMBE FOR CROSS SECTION.
      ELSE
        IF(IOFLG(1,3).NE.0) THEN
          WRITE(IOUT,74) MBEGWUNT,KSTP,KPER
          IF(CHEDFM.EQ.' ') THEN
             CALL ULASAV(BUFF,TEXT,KSTP,KPER,PERTIM,TOTIM,NCOL,
     1                NLAY,-1,MBEGWUNT)
          ELSE
             CALL ULASV2(BUFF,TEXT,KSTP,KPER,PERTIM,TOTIM,NCOL,
     1                  NLAY,-1,MBEGWUNT,CHEDFM,LBHDSV,IBOUND)
          END IF
        END IF
      END IF
80    CONTINUE
      DEALLOCATE(BUFF)
C
C6------RETURN.
      RETURN
      END
      SUBROUTINE SCLN1TE(KSTP,KPER,IPFLG,ISA)
C     ******************************************************************
C     PRINT AND RECORD TRANSPORT MASS BALANCE ERROR IN CLN CELLS
C     RESULTING FROM FLOW BALANCE ERROR
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,      ONLY:FMBE,IBOUND,IOUT,NODES,
     1   IFMBC,MBEGWUNF,MBEGWUNT,MBECLNUNF,MBECLNUNT
      USE CLN1MODULE, ONLY:  NCLNNDS,ICLNHD
      USE GWFBASMODULE,ONLY:PERTIM,TOTIM,IHEDFM,IHEDUN,LBDDSV,
     2                      CHEDFM,CDDNFM,IOFLG
C
      CHARACTER*16 TEXT
      DOUBLE PRECISION SSTRT
      REAL,          SAVE,    DIMENSION(:),    ALLOCATABLE ::BUFF
C
      DATA TEXT /'CLN TRNS BAL ERR'/
C     ------------------------------------------------------------------
      ALLOCATE(BUFF(NCLNNDS))
C
C1------FOR EACH CLN NODE PUT FMBE IN BUFF IF PRINT OR SAVE IS REQUESTED.
      DO 59 N=1,NCLNNDS
C
C2------Save FMBE in buffer array BUFF
        NG = N+NODES
        BUFF(N)=FMBE(NG)
   59 CONTINUE
C
C3------CALL ULAPRS OR ULAPRW TO PRINT FMBE.
      IF(ISA.NE.0) THEN
        IF(IOFLG(1,1).NE.0) THEN
          IF(IHEDFM.LT.0) CALL ULAPRS(BUFF(1),TEXT,KSTP,KPER,
     1                  NCLNNDS,1,1,-IHEDFM,IOUT)
          IF(IHEDFM.GE.0) CALL ULAPRW(BUFF(1),TEXT,KSTP,KPER,
     1                  NCLNNDS,1,1,IHEDFM,IOUT)
          IPFLG=1
        ENDIF
C
      END IF
C
C4------DETERMINE IF FMBE SHOULD BE SAVED.
C4------IF SO THEN CALL A ULASAV OR ULASV2 TO RECORD FMBE.
      IFIRST=1
      IF(MBECLNUNT.LE.0) GO TO 80
        NSTRT = NODES+1
        IF(IOFLG(1,3).EQ.0) GO TO 80
        IF(IFIRST.EQ.1) WRITE(IOUT,74) MBECLNUNT,KSTP,KPER
   74   FORMAT(1X,/1X,'CLN FMBE WILL BE SAVED ON UNIT ',I4,
     1      ' AT END OF TIME STEP ',I3,', STRESS PERIOD ',I4)
        IFIRST=0
        IF(CHEDFM.EQ.' ') THEN
           CALL ULASAV(BUFF(1),TEXT,KSTP,KPER,PERTIM,TOTIM,NCLNNDS,
     1                1,1,MBECLNUNT)
        ELSE
           CALL ULASV2(BUFF(1),TEXT,KSTP,KPER,PERTIM,TOTIM,NCLNNDS,
     1                1,1,MBECLNUNT,CHEDFM,LBDDSV,IBOUND(NSTRT))
        END IF
C
80    CONTINUE
      DEALLOCATE(BUFF)

C
C5------RETURN.
      RETURN
      END
C
C-----------------------------------------------------------------------
      MODULE GWTPCBMODULE
        INTEGER,SAVE,POINTER  ::NPCB,MXPCB,IPCBCB,NPCBVL,IPRPCB
        INTEGER,SAVE,POINTER  ::NPPCB,IPCBPB,NNPPCB
        CHARACTER(LEN=16),SAVE, DIMENSION(:),   ALLOCATABLE     ::PCBAUX
        REAL,             SAVE, DIMENSION(:,:), ALLOCATABLE     ::PCB
        DOUBLE PRECISION, SAVE, DIMENSION(:,:), ALLOCATABLE  ::
     1      AMATDIAG,RHSKPT
        DOUBLE PRECISION, SAVE, DIMENSION(:), ALLOCATABLE :: 
     1   AMATSKP,AMATTKP
      END MODULE GWTPCBMODULE
C
      SUBROUTINE GWT2PCB1AR(IN)
C     ******************************************************************
C     ALLOCATE ARRAY STORAGE FOR PRESCRIBED CONCENTRATION BOUNDARY PACKAGE
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,  ONLY:IOUT,NCOL,NROW,NLAY,IFREFM,NODES,IUNSTR,ITRNSP,
     1                  NEQS,IA,NJA
      USE GWTPCBMODULE
      USE GWTBCTMODULE, ONLY: MCOMP,MCOMPT,CONC,IHEAT
C
      CHARACTER*400 LINE
C     ------------------------------------------------------------------
C-----PCB REQUIRED ONLY IF TRANSPORT SIMULATION IS PERFORMED
      IF(ITRNSP.EQ.0)THEN
        IN = 0
        RETURN
      ENDIF
C     ------------------------------------------------------------------
      ALLOCATE(NPCB,MXPCB,IPCBCB,NPCBVL,IPRPCB)
      ALLOCATE(NPPCB,IPCBPB,NNPPCB)
C
C1------IDENTIFY PACKAGE AND INITIALIZE NPCB.
      WRITE(IOUT,1)IN
    1 FORMAT(1X,/1X,'PCB -- PRESCRIBED CONCENTRATION PACKAGE,',1X,
     1 'VERSION 7, 2/2/2010 INPUT READ FROM UNIT ',I7)
      NPCB=0
      NNPPCB=0
C
C1A-----FOR ITRNSP=5 ALLOCATE PCB ARRAYS FOR ALL GWF CELLS AND RETURN
      IF(ITRNSP.EQ.5)THEN
        NPCB = NEQS * MCOMPT
        MXPCB = NPCB
        NPCBVL = 5
        NAUX = 0
        ALLOCATE (PCB(NPCBVL,MXPCB))
        ALLOCATE (AMATDIAG(MXPCB,MCOMPT),RHSKPT(MXPCB,MCOMPT))
C---------------------------------------------------------------------------------      
C4i-------FOR HEAT EQUATION AFTER SOLUTE, NEED TO STORE SOLUTE AND HEAT AMATS
      IF(MCOMPT.GT.1. AND. IHEAT.EQ.1) THEN  
C4ii ------ ALLOCATE ARRAYS TO KEEP SOLUTE AND TEMPERATURE AMAT ARRAYS
          NJAPCB = 1  
          DO N = 1,MXPCB
            DO II = IA(N)+1,IA(N+1)-1
              NJAPCB = NJAPCB + 1 
            ENDDO               
          ENDDO
          ALLOCATE(AMATSKP(NJAPCB), AMATTKP(NJAPCB))     
C ----NOTE: Need to store only for PCB rows as only those are needed for mass balance - need to trim       
      ENDIF      
C ---------------------------------------------------------------------------------        
        RETURN
      ENDIF
C
C2------READ MAXIMUM NUMBER OF PCBS AND UNIT OR FLAG FOR
C2------CELL-BY-CELL FLOW TERMS.
      CALL URDCOM(IN,IOUT,LINE)
      CALL UPARLSTAL(IN,IOUT,LINE,NPPCB,MXPW)
      IF(IFREFM.EQ.0) THEN
         READ(LINE,'(2I10)') MXPCB,IPCBCB
         LLOC=21
      ELSE
         LLOC=1
         CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,MXPCB,R,IOUT,IN)
         CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IPCBCB,R,IOUT,IN)
      END IF
      WRITE(IOUT,3) MXPCB
    3 FORMAT(1X,'MAXIMUM OF ',I6,' ACTIVE PRESCRIBED CONCS AT ONE TIME')
      IF(IPCBCB.LT.0) WRITE(IOUT,7)
    7 FORMAT(1X,'CELL-BY-CELL FLUXES WILL BE PRINTED WHEN ICBCFL NOT 0')
      IF(IPCBCB.GT.0) WRITE(IOUT,8) IPCBCB
    8 FORMAT(1X,'CELL-BY-CELL FLOXES WILL BE SAVED ON UNIT ',I4)
C
C3------READ PRINT FLAG.
      ALLOCATE(PCBAUX(20))
      IPRPCB=1
   10 CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,N,R,IOUT,IN)
      IF(LINE(ISTART:ISTOP).EQ.'NOPRINT') THEN
         WRITE(IOUT,13)
   13    FORMAT(1X,'LIST OF PRESCRIBED CONC CELLS WILL NOT BE PRINTED')
         IPRPCB = 0
         GO TO 10
      END IF
C
C4------ALLOCATE SPACE FOR THE PCB DATA.
      IPCBPB=MXPCB+1
      MXPCB=MXPCB+MXPW
      IF(MXPCB.LT.1) THEN
         WRITE(IOUT,17)
   17    FORMAT(1X,
     1'Deactivating the PCBl Package because MXPCB=0')
         IN=0
      END IF
      NPCBVL = 5
      NAUX = 0
      ALLOCATE (PCB(NPCBVL,MXPCB))
      ALLOCATE (AMATDIAG(MXPCB,MCOMPT),RHSKPT(MXPCB,MCOMPT)) 
C---------------------------------------------------------------------------------      
C4i-------FOR HEAT EQUATION AFTER SOLUTE, NEED TO STORE SOLUTE AND HEAT AMATS
      IF(MCOMPT.GT.1. AND. IHEAT.EQ.1) THEN  
C4ii ------ ALLOCATE ARRAYS TO KEEP SOLUTE AND TEMPERATURE AMAT ARRAYS
          NJAPCB = 1  
          DO N = 1,MXPCB
            DO II = IA(N)+1,IA(N+1)-1
              NJAPCB = NJAPCB + 1 
            ENDDO               
          ENDDO
          ALLOCATE(AMATSKP(NJAPCB), AMATTKP(NJAPCB))     
C ----NOTE: Need to store only for PCB rows as only those are needed for mass balance - need to trim       
      ENDIF      
C ---------------------------------------------------------------------------------
C5------READ NAMED PARAMETERS.
      WRITE(IOUT,18) NPPCB
   18 FORMAT(1X,//1X,I5,' PCB parameters')
      IF(NPPCB.GT.0) THEN
        LSTSUM=IPCBPB
        DO 120 K=1,NPPCB
          LSTBEG=LSTSUM
          CALL UPARLSTRP(LSTSUM,MXPCB,IN,IOUT,IP,'PCB','Q',1,
     &                   NUMINST)
          NLST=LSTSUM-LSTBEG
          IF(NUMINST.EQ.0) THEN
C5A-----READ PARAMETER WITHOUT INSTANCES.
            IF(IUNSTR.EQ.0)THEN
              CALL ULSTRD(NLST,PCB,LSTBEG,NPCBVL,MXPCB,1,IN,IOUT,
     &      '  PCB NO.  LAYER   ROW   COL   SPECIES NO.  STRESS FACTOR',
     &        PCBAUX,20,NAUX,IFREFM,NCOL,NROW,NLAY,5,5,IPRPCB)
            ELSE
             CALL ULSTRDU(NLST,PCB,LSTBEG,NPCBVL,MXPCB,1,IN,IOUT,
     &      '  PCB NO.  LAYER   ROW   COL   SPECIES NO.  STRESS FACTOR',
     &        PCBAUX,20,NAUX,IFREFM,NEQS,5,5,IPRPCB)
            ENDIF
          ELSE
C5B-----READ INSTANCES.
            NINLST=NLST/NUMINST
            DO 110 I=1,NUMINST
            CALL UINSRP(I,IN,IOUT,IP,IPRPCB)
            IF(IUNSTR.EQ.0)THEN
              CALL ULSTRD(NINLST,PCB,LSTBEG,NPCBVL,MXPCB,1,IN,IOUT,
     &      '  PCB NO.  LAYER   ROW   COL   SPECIES NO.  STRESS FACTOR',
     &        PCBAUX,20,NAUX,IFREFM,NCOL,NROW,NLAY,5,5,IPRPCB)
            ELSE
             CALL ULSTRDU(NINLST,PCB,LSTBEG,NPCBVL,MXPCB,1,IN,IOUT,
     &      '  PCB NO.  LAYER   ROW   COL   SPECIES NO.  STRESS FACTOR',
     &        PCBAUX,20,NAUX,IFREFM,NEQS,5,5,IPRPCB)
            ENDIF
            LSTBEG=LSTBEG+NINLST
  110       CONTINUE
          END IF
  120   CONTINUE
      END IF
C
C6------RETURN
      RETURN
      END
      SUBROUTINE GWT2PCB1RP(IN,KPER)
C     ******************************************************************
C     READ PRESCRIBED CONCENTRATION DATA FOR A STRESS PERIOD
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,       ONLY:IOUT,NCOL,NROW,NLAY,IFREFM,NODES,IUNSTR,
     1                       NEQS,ITRNSP,IDPIN,IUNIT
      USE GWTPCBMODULE, ONLY:NPCB,MXPCB,NPCBVL,IPRPCB,NPPCB,
     1                       IPCBPB,NNPPCB,PCBAUX,PCB
      USE GWTBCTMODULE, ONLY: MCOMPT,CONC
      USE GWTDPTMODULE, ONLY: CONCIM 
      DOUBLE PRECISION PERTIMREAD,TOTIMREAD
C
      CHARACTER*6 CPCB
C     ------------------------------------------------------------------
C0------ITRNSP = 5 LOADS UP MATRIX DIFFUSION TERM FROM FRACTURE CONC
      IF(ITRNSP.EQ.5) THEN 
C0A ----READ CONC ARRAY IF NOT FIRST TIME STEP
        IF(KPER.GT.1) THEN 
          IGWFHD = IUNIT(15)
          ICLNHD = 0  ! IUNIT(29) no need for cln input as it does not affect MD
          IDPFHD = 0
          DO ICOMP=1,MCOMPT
            WRITE(IOUT,351) ICOMP
351        FORMAT(9X,'*** READING CONC OF COMPONENT',I5,' ***')
            IF(IDPIN.EQ.0)THEN
              CALL READVAR4 (IGWFHD,ICLNHD,IDPFHD,
     1          CONC(1,ICOMP),CONCIM(1,ICOMP))
            ELSE
              CALL READVAR8 (IGWFHD,ICLNHD,IDPFHD,
     1          CONC(1,ICOMP),CONCIM(1,ICOMP)) 
            ENDIF          
          ENDDO                
        ENDIF
C0B-------TRANSFER CONC ARRAY AND LOCATION INTO PCB ARRAYS        
        L = 0
        DO ICOMP = 1,MCOMPT
          DO N = 1,NEQS
            L = L + 1
            PCB(1,L) = N        ! NODE NUMBER IS IN LOCATION 1
            PCB(4,L) = ICOMP    ! NODE NUMBER IS IN LOCATION 1
            PCB(5,L) = CONC(N,ICOMP)
          ENDDO
        ENDDO          
C        
        RETURN
      ENDIF    
C----------------------------------------------------------------------
C1------IDENTIFY PACKAGE.
      WRITE(IOUT,1)IN
   1  FORMAT(1X,/1X,'PCB -- PRESCRIBED CONCENTRATION PACKAGE,',1X,
     1 'VERSION 7, 2/2/2010 INPUT READ FROM UNIT ',I4)
C
C2----READ NUMBER OF PCBS (OR FLAG SAYING REUSE PCB DATA).
C2----AND NUMBER OF PARAMETERS
      IF(NPPCB.GT.0) THEN
        IF(IFREFM.EQ.0) THEN
           READ(IN,'(2I10)') ITMP,NP
        ELSE
           READ(IN,*) ITMP,NP
        END IF
      ELSE
         NP=0
         IF(IFREFM.EQ.0) THEN
            READ(IN,'(I10)') ITMP
         ELSE
            READ(IN,*) ITMP
         END IF
      END IF
C
C3------Calculate some constants.
      NAUX=NPCBVL-5
      IOUTU = IOUT
      IF (IPRPCB.EQ.0) IOUTU=-IOUTU
C
C4-----IF ITMP LESS THAN ZERO REUSE NON-PARAMETER DATA. PRINT MESSAGE.
C4-----IF ITMP=>0, SET NUMBER OF NON-PARAMETER PCBS EQUAL TO ITMP.
      IF(ITMP.LT.0) THEN
         WRITE(IOUT,6)
    6    FORMAT(1X,/
     1    1X,'REUSING NON-PARAMETER PCBS FROM LAST STRESS PERIOD')
      ELSE
         NNPPCB=ITMP
      END IF
C
C5-----IF THERE ARE NEW NON-PARAMETER PCBs, READ THEM.
      MXPCB=IPCBPB-1
      IF(ITMP.GT.0) THEN
         IF(NNPPCB.GT.MXPCB) THEN
            WRITE(IOUT,99) NNPPCB,MXPCB
   99       FORMAT(1X,/1X,'THE NUMBER OF ACTIVE PCBs (',I6,
     1                     ') IS GREATER THAN MXPCB(',I6,')')
            CALL USTOP(' ')
         END IF
         IF(IUNSTR.EQ.0)THEN
           CALL ULSTRD(NNPPCB,PCB,1,NPCBVL,MXPCB,0,IN,IOUT,
     &      'PCB  NO.  LAYER   ROW   COL  COMPONENT NO.  STRESS RATE',
     2             PCBAUX,20,NAUX,IFREFM,NCOL,NROW,NLAY,5,5,IPRPCB)
          ELSE
             CALL ULSTRDU(NNPPCB,PCB,1,NPCBVL,MXPCB,0,IN,IOUT,
     &      'PCB  NO.  LAYER   ROW   COL  COMPONENT NO.  STRESS FACTOR',
     &        PCBAUX,20,NAUX,IFREFM,NEQS,5,5,IPRPCB)
          ENDIF
      END IF
      NPCB=NNPPCB
C
C6-----IF THERE ARE ACTIVE PCB PARAMETERS, READ THEM AND SUBSTITUTE
      CALL PRESET('Q')
      NREAD=NPCBVL-1
      IF(NP.GT.0) THEN
         DO 30 N=1,NP
         CALL UPARLSTSUB(IN,'PCB',IOUTU,'Q',PCB,NPCBVL,MXPCB,NREAD,
     1                MXPCB,NPCB,5,5,
     &      'PCB  NO.  LAYER   ROW   COL  COMPONENT NO.  STRESS RATE',
     3            PCBAUX,20,NAUX)
   30    CONTINUE
      END IF
C
C7------PRINT NUMBER OF PCBS IN CURRENT STRESS PERIOD.
      CPCB=' PCBs '
      IF(NPCB.EQ.1) CPCB=' PCBS '
      WRITE(IOUT,101) NPCB,CPCB
  101 FORMAT(1X,/1X,I6,A)
C
C8-------FOR STRUCTURED GRID, CALCULATE NODE NUMBER AND PLACE IN LAYER LOCATION
      IF(ITMP.GT.0.AND.IUNSTR.EQ.0)THEN
        DO L=1,NPCB
          IR=PCB(2,L)
          IC=PCB(3,L)
          IL=PCB(1,L)
          N = IC + NCOL*(IR-1) + (IL-1)* NROW*NCOL
          PCB(1,L) = N
        ENDDO
      ENDIF
C
C9------RETURN
      RETURN
      END
      SUBROUTINE GWT2PCB1FM(ICOMP)
C     ******************************************************************
C     PRESCRIBE CONCENTRATIONS AT PCB CELLS
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,       ONLY:IBOUND,RHS,IA,JA,AMAT
      USE GWTPCBMODULE, ONLY:NPCB,PCB,AMATDIAG,RHSKPT,AMATSKP,AMATTKP
      USE GWTBCTMODULE, ONLY: ICBUND,CONC,IHEAT,MCOMPT
      DOUBLE PRECISION BIG,Co
C     ------------------------------------------------------------------
C
C1------IF NUMBER OF PCBs <= 0 THEN RETURN.
      IF(NPCB.LE.0) RETURN
      BIG = 1.0E20
C
C2------PROCESS EACH PCB IN THE LIST.
      IIKP = 1
      DO 100 L=1,NPCB
      N=PCB(1,L)
      IC=PCB(4,L)
      Co=PCB(5,L)
C
C2A-----IF THE CELL IS INACTIVE OR WRONG COMPONENT SPECIES THEN BYPASS PROCESSING.
      IF(ICBUND(N).EQ.0.OR.IC.NE.ICOMP) GO TO 100
C
C2B-----IF THE CELL IS PRESCRIBED CONCENTRATION THEN PROCESS
C---------------------------------------------------------------------------------      
C4i-------FOR HEAT EQUATION AFTER SOLUTE, NEED TO STORE SOLUTE AND HEAT AMATS
      IF(MCOMPT.GT.1. AND. IHEAT.EQ.1) THEN  
C4ii ------ RETRIEVE ARRAYS FOR SOLUTE AND TEMPERATURE AS APPROPRIATE
        IF(ICOMP.EQ.1)THEN !RETRIEVE SOLUTE ARRAY
          DO II = IA(N)+1,IA(N+1)-1
            AMATSKP(IIKP) = AMAT(II) 
            IIKP = IIKP + 1
          ENDDO          
        ELSEIF(ICOMP.EQ.MCOMPT) THEN  ! RETRIEVE TEMPERATURE ARRAY
          DO II = IA(N)+1,IA(N+1)-1
            AMATTKP(IIKP) = AMAT(II)  
            IIKP = IIKP + 1
          ENDDO             
        ENDIF  
      ENDIF        
C --------------------------------------------------------------------------------      
        AMATDIAG(L,ICOMP) = AMAT(IA(N))
        RHSKPT(L,ICOMP) = RHS(N)
        AMAT(IA(N)) = -1.0 * BIG
        RHS(N) = -Co * BIG
  100 CONTINUE
C
C3------RETURN
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE GWT2PCB1BD(KSTP,KPER,ICOMP)
C     ******************************************************************
C     CALCULATE MASS BUDGET FOR PCB CELLS
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,   ONLY:IOUT,NCOL,NROW,NLAY,IBOUND,BUFF,NODES,IUNSTR,
     1              AMAT,IA,JA,INCLN,IDPT,ITRNSP
      USE GWFBASMODULE,ONLY:MSUM,ICBCFL,IAUXSV,DELT,PERTIM,TOTIM
      USE GWTPCBMODULE,ONLY:NPCB,IPCBCB,PCB,NPCBVL,PCBAUX,
     *  AMATDIAG,RHSKPT,AMATSKP,AMATTKP
      USE GWTBCTMODULE, ONLY: ICBUND,CONC,MSUMT,VBVLT,VBNMT,IPCBFLAG,
     1  MCOMPT,IHEAT
      USE CLN1MODULE, ONLY: ICLNMB,NCLNNDS
      USE GWTDPTMODULE, ONLY: CONCIM,DDTTR

C
      CHARACTER*16 TEXT
      DOUBLE PRECISION RATIN,RATOUT,QQ,CDIFF
      DATA TEXT /'PRESCRIBED CONCS'/
C     ------------------------------------------------------------------
C
C1------CLEAR RATIN AND RATOUT ACCUMULATORS, AND SET CELL-BY-CELL
C1------BUDGET FLAG.
      ZERO=0.
      RATIN=ZERO
      RATOUT=ZERO
      IBD=0
      IF(IPCBCB.LT.0 .AND. ICBCFL.NE.0) IBD=-1
      IF(IPCBCB.GT.0) IBD=ICBCFL
      IF(ITRNSP.EQ.5) IBD = 0 ! OVERRIDE WHEN ITRNSP = 5
      IBDLBL=0
C
C2-----IF CELL-BY-CELL FLOWS WILL BE SAVED AS A LIST, WRITE HEADER.
      IF(IBD.EQ.2) THEN
         NAUX=NPCBVL-5
         IF(IAUXSV.EQ.0) NAUX=0
         IF(IUNSTR.EQ.0) THEN
           CALL UBDSV4(KSTP,KPER,TEXT,NAUX,PCBAUX,IPCBCB,NCOL,NROW,NLAY,
     1          NPCB,IOUT,DELT,PERTIM,TOTIM,ICBUND)
         ELSE
           CALL UBDSV4U(KSTP,KPER,TEXT,NAUX,PCBAUX,IPCBCB,NODES,
     1          NPCB,IOUT,DELT,PERTIM,TOTIM,ICBUND)
         ENDIF
      END IF
C
C3------CLEAR THE BUFFER.
      DO 50 N=1,NODES
      BUFF(N)=ZERO
      IPCBFLAG(N) = 0
50    CONTINUE
C
C4------IF THERE ARE NO PCBs, DO NOT ACCUMULATE
      IF(NPCB.EQ.0) GO TO 200
      IF(ICOMP.GT.MCOMPT) GO TO 200
C
C5------LOOP THROUGH EACH PCB CALCULATING MASS FLUX
      IIKP = 1
      DO 100 L=1,NPCB
C
C5A-----GET NODE NUMBER OF CELL CONTAINING PCB.
      N=PCB(1,L)
      IC=PCB(4,L)
      QQ=ZERO
      Q = 0.0
C
C5B-----IF THE CELL IS NOT PCB OR WRONG COMPONENT SPECIES, IGNORE IT.
      IF(ICBUND(N).EQ.0.OR.IC.NE.ICOMP)GO TO 99
      IPCBFLAG(N) = 1
C---------------------------------------------------------------------------------      
C4i-------FOR HEAT EQUATION AFTER SOLUTE, NEED TO STORE SOLUTE AND HEAT AMATS
      IF(MCOMPT.GT.1. AND. IHEAT.EQ.1) THEN  
C4ii ------ RETRIEVE ARRAYS FOR SOLUTE AND TEMPERATURE AS APPROPRIATE
        IF(ICOMP.EQ.1)THEN !RETRIEVE SOLUTE ARRAY
          DO II = IA(N)+1,IA(N+1)-1
            AMAT(II) = AMATSKP(IIKP) 
            IIKP = IIKP + 1
          ENDDO          
        ELSEIF(ICOMP.EQ.MCOMPT) THEN  ! RETRIEVE TEMPERATURE ARRAY
          DO II = IA(N)+1,IA(N+1)-1
            AMAT(II) = AMATTKP(IIKP) 
            IIKP = IIKP + 1
          ENDDO             
        ENDIF  
      ENDIF        
C --------------------------------------------------------------------------------
C5C---BACK-CALCULATE MASS FLUX FOR PCB NODE.
      QQ = -AMATDIAG(L,ICOMP) * CONC(N,ICOMP)
      DO II = IA(N)+1,IA(N+1)-1
        JJ = JA(II)
        QQ = QQ - AMAT(II) * CONC(JJ,ICOMP)
      ENDDO
C--------------------------------------------------------
C ------DUAL POROSITY TRANSPORT INTO /OUT OF MATRIX TERM
      IF(IDPT. NE.0) THEN
        QQ=QQ-(DDTTR(N) * (CONCIM(N,ICOMP)-CONC(N,ICOMP)))
      ENDIF
      QQ=QQ + RHSKPT(L,ICOMP)
      Q = QQ
C
C5D-----PRINT FLOW RATE IF REQUESTED.
      IF(IBD.LT.0) THEN
         IF(IBDLBL.EQ.0) WRITE(IOUT,61) TEXT,KPER,KSTP
   61    FORMAT(1X,/1X,A,'   PERIOD ',I5,'   STEP ',I5)
        IF(IUNSTR.EQ.0)THEN
          IL = (N-1) / (NCOL*NROW) + 1
          IJ = N - (IL-1)*NCOL*NROW
          IR = (IJ-1)/NCOL + 1
          IC = IJ - (IR-1)*NCOL
           WRITE(IOUT,62) L,IL,IR,IC,Q
   62    FORMAT(1X,'PCB  ',I6,'   LAYER ',I3,'   ROW ',I5,'   COL ',I5,
     1       '   FLUX ',1PG15.6)
        ELSE
           WRITE(IOUT,63) L,N,Q
   63    FORMAT(1X,'PCB  ',I6,'    NODE ',I8,'   FLUX ',1PG15.6)
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
C5G-----FLOW RATE IS POSITIVE (RECHARGE). ADD IT TO RATIN.
        RATIN=RATIN+QQ
      ELSE
C
C5H-----FLOW RATE IS NEGATIVE (DISCHARGE). ADD IT TO RATOUT.
        RATOUT=RATOUT-QQ
      END IF
C
C5I-----IF SAVING CELL-BY-CELL FLOWS IN A LIST, WRITE FLOW.  ALSO
C5I-----COPY FLOW TO PCB LIST.

   99 CONTINUE
      IF(IUNSTR.EQ.0)THEN
        IF(IBD.EQ.2) CALL UBDSVB(IPCBCB,NCOL,NROW,IC,IR,IL,Q,
     1                  PCB(1,L),NPCBVL,NAUX,6,ICBUND,NLAY)
      ELSE
        IF(IBD.EQ.2) CALL UBDSVBU(IPCBCB,NODES,N,Q,
     1                  PCB(1,L),NPCBVL,NAUX,6,ICBUND)
      ENDIF
100   CONTINUE
C
C6------IF CELL-BY-CELL FLOWS WILL BE SAVED AS A 3-D ARRAY,
C6------CALL UBUDSV TO SAVE THEM.
      IF(IUNSTR.EQ.0)THEN
        IF(IBD.EQ.1)CALL UBUDSV(KSTP,KPER,TEXT,IPCBCB,BUFF(1),NCOL,NROW,
     1                          NLAY,IOUT)
      ELSE
        IF(IBD.EQ.1) CALL UBUDSVU(KSTP,KPER,TEXT,IPCBCB,BUFF(1),NODES,
     1                          IOUT,PERTIM,TOTIM)
      ENDIF
      IF(IBD.EQ.1.AND.INCLN.GT.0)THEN
        IICLNCB = ICLNMB
        IF(IICLNCB.GT.0) CALL UBUDSVU(KSTP,KPER,TEXT,IICLNCB,
     1    BUFF(NODES+1),NCLNNDS,IOUT,PERTIM,TOTIM)
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
C--------------------------------------------------------------------------
      SUBROUTINE GWT2PCB1DA
C  Deallocate PCB MEMORY
      USE GWTPCBMODULE
C
        DEALLOCATE(NPCB)
        DEALLOCATE(MXPCB)
        DEALLOCATE(NPCBVL)
        DEALLOCATE(IPCBCB)
        DEALLOCATE(IPRPCB)
        DEALLOCATE(NPPCB)
        DEALLOCATE(IPCBPB)
        DEALLOCATE(NNPPCB)
        DEALLOCATE(PCBAUX)
        DEALLOCATE(PCB)
        DEALLOCATE(AMATDIAG)
        DEALLOCATE(RHSKPT)
C
      RETURN
      END
C----------------------------------------------------------------------------
      SUBROUTINE GWT2PHB1FM(ICOMP)
C     ******************************************************************
C     FORMULATE BOUNDARY CONDITION FOR TRANSPORT AT PRESCRIBED HEAD NODES
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,       ONLY:IBOUND,RHS,IA,JA,AMAT,NODES,IUNIT,NEQS
      USE GWTBCTMODULE, ONLY: ICBUND,CBCH
      USE GWFCHDMODULE,ONLY:NCHDS,MXCHD,NCHDVL,IPRCHD,NPCHD,ICHDPB,
     1                      NNPCHD,CHDAUX,CHDS

      DOUBLE PRECISION Co,Q
C     ------------------------------------------------------------------
C-----INITIALIZE CONSTANTS AND FLAGS
c      BIG = 1.0E20
c      IF(IUNIT(20).GT.0)THEN  ! READ FROM chd FILE IF CHD IS USED
c        NAUX = NCHDVL- 5
c        CALL CONCIAUX(ICOMP,NAUX,CHDAUX,IAUX)
c      ENDIF
cC
cC1------PROCESS EACH NODE IN THE LIST.
c      L = 0
c      DO 100 N=1,NEQS
cC
cC2------IF THE CELL IS INACTIVE OR NOT PRESCRIBED HEAD NODE, BYPASS PROCESSING.
c      IF(ICBUND(N).EQ.0.OR.IBOUND(N).GE.0) GO TO 100
c        Q = CBCH(N)
c        Co= 0.0
cC3--------IF CHD IS ON THEN READ CONCENTRATION FROM CHD FILE SEQUENTIALLY
cC3--------NOTE THAT IF CHD IS ON THEN SHOULD NOT HAVE IBOUND=0 AS WELL DEFINING OTHER NODES
c        IF(IUNIT(20).GT.0.AND.IAUX.NE.0)THEN
c          IF(NCHDVL.GT.5)THEN
c            L = L + 1
c            Co = CHDS(5+IAUX,L)
c          ENDIF
c        ENDIF
cC
cC4------IF THE CELL IS OUTFLOW, PUT Q ON LHS DIAGONAL
c        IF(Q.LE.0.0)THEN
c          AMAT(IA(N)) = AMAT(IA(N)) + Q
c        ELSE
cC5--------IF THE CELL IS INFLOW, PUT Q*Co ON RHS
c          RHS(N) = RHS(N) -Co * Q
c        ENDIF
cC
c  100 CONTINUE
C----------------------------------------------------------------------
      IF(NCHDS.LE.0) RETURN
      BIG = 1.0E20
      NAUX=NCHDVL- 5
      CALL CONCIAUX(ICOMP,NAUX,CHDAUX,IAUX)
C
C2------PROCESS EACH GHB IN THE LIST.
      DO 100 L=1,NCHDS
      N=CHDS(1,L)
C
C2A-----IF THE CELL IS INACTIVE THEN BYPASS PROCESSING.
      IF(ICBUND(N).EQ.0) GO TO 100
        Q = CBCH(N)   !   CHDS(NCHDVL,L)
C
C2B-----IF THE CELL IS OUTFLOW, PUT Q ON LHS DIAGONAL
        IF(Q.LE.0.0)THEN
          AMAT(IA(N)) = AMAT(IA(N)) + Q
        ELSE
C2C-------IF THE CELL IS INFLOW, PUT Q*Co ON RHS
          Co = 0.0
          IF(IAUX.NE.0) Co= CHDS(5+IAUX,L)
          RHS(N) = RHS(N) -Co * Q
        ENDIF
C
  100 CONTINUE  
C
C6------RETURN
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE GWT2PHB1BD(KSTP,KPER,ICOMP)
C     ******************************************************************
C     CALCULATE MASS BUDGET FOR TRANSPORT AT PRESCRIBED HEAD NODES
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,   ONLY:IOUT,NCOL,NROW,NLAY,IBOUND,BUFF,NODES,IUNSTR,
     1              AMAT,IA,JA,IUNIT,INCLN
      USE GWFBASMODULE,ONLY:MSUM,ICBCFL,IAUXSV,DELT,PERTIM,TOTIM
      USE GWTBCTMODULE,ONLY: ICBUND,CONC,MSUMT,VBVLT,VBNMT,CBCH,
     1  IPCBFLAG,IBCTCB,MCOMPT
      USE GWFCHDMODULE,ONLY:NCHDS,NCHDVL,CHDAUX,CHDS
      USE CLN1MODULE, ONLY: ICLNMB,NCLNNDS
C
      CHARACTER*16 TEXT
      DOUBLE PRECISION RATIN,RATOUT,QQ,BIG,CDIFF,Co,QF
      DATA TEXT /'CNST H MASS FLUX'/
C     ------------------------------------------------------------------
C
C1------CLEAR RATIN AND RATOUT ACCUMULATORS, AND SET CELL-BY-CELL
C1------BUDGET FLAG.
      ZERO=0.
      RATIN=ZERO
      RATOUT=ZERO
      IBD=0
      IF(IBCTCB.LT.0 .AND. ICBCFL.NE.0) IBD=-1
      IF(IBCTCB.GT.0) IBD=ICBCFL
      IBDLBL=0
      IF(IUNIT(20).GT.0)THEN
        NAUX = NCHDVL- 5
        CALL CONCIAUX(ICOMP,NAUX,CHDAUX,IAUX)
      ENDIF
C
C2-----IF CELL-BY-CELL FLOWS WILL BE SAVED AS A LIST, WRITE HEADER.
      IICLNCB=0
      IF(IBD.EQ.2.AND.IUNIT(20).GT.0) THEN
         NAUX=NCHDVL-5
         IF(IAUXSV.EQ.0) NAUX=0
         NNCLNNDS=0
         IF(INCLN.GT.0) THEN
           IICLNCB=ICLNMB
           NNCLNNDS=NCLNNDS
         ENDIF
        CALL UBDSVHDR(IUNSTR,KSTP,KPER,IOUT,IBCTCB,IICLNCB,NODES,
     1    NNCLNNDS,NCOL,NROW,NLAY,NCHDS,NCHDVL,NAUX,IBOUND,
     2    TEXT,CHDAUX,DELT,PERTIM,TOTIM,CHDS)
      END IF
C
C3------CLEAR THE BUFFER.
      DO 50 N=1,NODES
      BUFF(N)=ZERO
50    CONTINUE
C
C4------IF THERE ARE NO GHBs, DO NOT ACCUMULATE
      IF(NCHDS.EQ.0) GO TO 200
      IF(ICOMP.GT.MCOMPT) GO TO 200
C
C5------LOOP THROUGH EACH PRESCRIBED HEAD BOUNDARY CALCULATING MASS FLUX
      L = 0
c--------------------------------------------------------------------------      
c      DO 100 N=1,NODES
c      QQ = 0.0
c      Q = 0.0    
cC
cC6------IF THE CELL IS INACTIVE OR NOT PRESCRIBED HEAD NODE, IGNORE IT.
c      IF(ICBUND(N).EQ.0.OR.IBOUND(N).GE.0.OR.IPCBFLAG(N).EQ.1)GO TO 99
cC
cC7------COMPUTE MASS FLUX AT RPESCRIBED HEAD BOUNDARY
c        QF = CBCH(N)
c        Co= 0.0
cC8--------IF CHD IS ON THEN READ CONCENTRATION FROM CHD FILE SEQUENTIALLY
cC8--------NOTE THAT IF CHD IS ON THEN SHOULD NOT HAVE IBOUND=0 AS WELL DEFINING OTHER NODES
c          IF(IUNIT(20).GT.0.AND.IAUX.NE.0)THEN
c            IF(NCHDVL.GT.5)THEN
c              L = L + 1
c              Co = CHDS(5+IAUX,L)
c            ENDIF
c          ENDIF
c          IF(L.EQ.0) L = 1
cC
cC9------IF THE CELL IS OUTFLOW, MASS IS Q * CONC
c        IF(QF.LE.0.0)THEN
c          QQ = QF * CONC(N,ICOMP)
c        ELSE
cC8--------IF THE CELL IS INFLOW, MASS IS Q*Co
c          QQ = QF * Co
c        ENDIF
c        Q = QQ
c--------------------------------------------------------------------------
      DO 100 L=1,NCHDS
C
C6-----GET NODE NUMBER OF CELL CONTAINING GHB.
      N=CHDS(1,L)
      QQ=ZERO
      Q = 0.0
C
C7A-----IF THE CELL IS INACTIVE OR PCB, IGNORE IT.
      IF(ICBUND(N).EQ.0.OR.IPCBFLAG(N).EQ.1)GO TO 99
C
C7B-----COMPUTE MASS FLUX AT GHB BOUNDARY
        QF = CBCH(N)     !     CHDS(NCHDVL,L) 
C
C8A-----IF THE CELL IS OUTFLOW, MASS IS Q * CONC
        IF(QF.LE.0.0)THEN
          QQ = QF * CONC(N,ICOMP)
        ELSE
C8B-------IF THE CELL IS INFLOW, MASS IS Q*Co
          Co = 0.0
          IF(IAUX.NE.0) Co = CHDS(5+IAUX,L)
          QQ = QF * Co
        ENDIF
        Q = QQ
C
C9------PRINT FLOW RATE IF REQUESTED.
        IF(IBD.LT.0) THEN
          IF(IBDLBL.EQ.0) WRITE(IOUT,61) TEXT,KPER,KSTP
   61     FORMAT(1X,/1X,A,'   PERIOD ',I4,'   STEP ',I3)
          IF(IUNSTR.EQ.0)THEN
            IL = (N-1) / (NCOL*NROW) + 1
            IJ = N - (IL-1)*NCOL*NROW
            IR = (IJ-1)/NCOL + 1
            IC = IJ - (IR-1)*NCOL
            WRITE(IOUT,62) L,IL,IR,IC,Q
   62       FORMAT(1X,'PHB  ',I6,'   LAYER ',I3,'   ROW ',I5,'   COL ',
     1       I5, '   FLUX ',1PG15.6)
          ELSE
            WRITE(IOUT,63) L,N,Q
   63       FORMAT(1X,'PHB  ',I6,'    NODE ',I8,'   FLUX ',1PG15.6)
          ENDIF
          IBDLBL=1
        END IF
C
C10-----ADD FLOW RATE TO BUFFER.
        BUFF(N)=BUFF(N)+Q
C
C11-----SEE IF FLUX IS POSITIVE OR NEGATIVE.
        IF(QQ.GE.ZERO) THEN
C
C12-------FLOW RATE IS POSITIVE (RECHARGE). ADD IT TO RATIN.
          RATIN=RATIN+QQ
        ELSE
C
C13-------FLOW RATE IS NEGATIVE (DISCHARGE). ADD IT TO RATOUT.
          RATOUT=RATOUT-QQ
        END IF
C
C14-----IF SAVING CELL-BY-CELL FLOWS IN A LIST, WRITE FLOW.  ALSO
C14-----COPY FLOW TO CONSANT HEAD LIST.
C   
      IF(IBD.EQ.2.AND.IUNIT(20).GT.0)THEN
        NN = N
        CALL UBDSVREC(IUNSTR,NN,NODES,NNCLNNDS,IBCTCB,IICLNCB,NCHDVL,
     1    5,NAUX,Q,CHDS(:,L),IBOUND,NCOL,NROW,NLAY)
      ENDIF
99    CONTINUE      
100   CONTINUE
C
C15------IF CELL-BY-CELL FLOWS WILL BE SAVED AS A 3-D ARRAY,
C15------CALL UBUDSV TO SAVE THEM.
      IF(IUNSTR.EQ.0)THEN
        IF(IBD.EQ.1) CALL UBUDSV(KSTP,KPER,TEXT,IBCTCB,BUFF,NCOL,NROW,
     1                          NLAY,IOUT)
      ELSE
        IF(IBD.EQ.1) CALL UBUDSVU(KSTP,KPER,TEXT,IBCTCB,BUFF,NODES,
     1                          IOUT,PERTIM,TOTIM)
      ENDIF
      IF(IBD.EQ.1.AND.INCLN.GT.0)THEN
        IF(IICLNCB.GT.0) CALL UBUDSVU(KSTP,KPER,TEXT,IICLNCB,
     1    BUFF(NODES+1),NCLNNDS,IOUT,PERTIM,TOTIM)
      ENDIF
C
C16------MOVE RATES, VOLUMES & LABELS INTO ARRAYS FOR PRINTING.
  200 RIN=RATIN
      ROUT=RATOUT
      VBVLT(3,MSUMT,ICOMP)=RIN
      VBVLT(4,MSUMT,ICOMP)=ROUT
      VBVLT(1,MSUMT,ICOMP)=VBVLT(1,MSUMT,ICOMP)+RATIN*DELT
      VBVLT(2,MSUMT,ICOMP)=VBVLT(2,MSUMT,ICOMP)+RATOUT*DELT
      VBNMT(MSUMT,ICOMP)=TEXT
C
C17------INCREMENT BUDGET TERM COUNTER(MSUM).
      MSUMT=MSUMT+1
C
C18------RETURN
      RETURN
      END
C----------------------------------------------------------------------------
      SUBROUTINE GWT2WEL1FM(ICOMP)
C     ******************************************************************
C     FORMULATE WELL TRANSPORT BOUNDARY CONDITION
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,       ONLY:IBOUND,RHS,IA,JA,AMAT
      USE GWFWELMODULE, ONLY:NWELLS,WELL,NWELVL,WELAUX
      USE GWTBCTMODULE, ONLY: ICBUND
      DOUBLE PRECISION Co,Q
C     ------------------------------------------------------------------
C
C1------IF NUMBER OF WELLS <= 0 THEN RETURN.
      IF(NWELLS.LE.0) RETURN
      BIG = 1.0E20
      NAUX =NWELVL - 5
      CALL CONCIAUX(ICOMP,NAUX,WELAUX,IAUX)
C
C2------PROCESS EACH WELL IN THE LIST.
      DO 100 L=1,NWELLS
      N=WELL(1,L)
C
C2A-----IF THE CELL IS INACTIVE THEN BYPASS PROCESSING.
      IF(ICBUND(N).EQ.0) GO TO 100
        Q = WELL(NWELVL,L)
C
C2B-----IF THE CELL IS OUTFLOW, PUT Q ON LHS DIAGONAL
        IF(Q.LE.0.0)THEN
          AMAT(IA(N)) = AMAT(IA(N)) + Q
        ELSE
C2C-------IF THE CELL IS INFLOW, PUT Q*Co ON RHS
          Co = 0.0
          IF(IAUX.NE.0) Co=WELL(4+IAUX,L)
          RHS(N) = RHS(N) -Co * Q
        ENDIF
C
  100 CONTINUE
C
C3------RETURN
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE GWT2WEL1BD(KSTP,KPER,ICOMP)
C     ******************************************************************
C     CALCULATE MASS BUDGET FOR WELL TRANSPORT BOUNDARY CONDITIONS
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,   ONLY:IOUT,NCOL,NROW,NLAY,IBOUND,BUFF,NODES,IUNSTR,
     1              AMAT,IA,JA,INCLN
      USE GWFBASMODULE,ONLY:MSUM,ICBCFL,IAUXSV,DELT,PERTIM,TOTIM
      USE GWFWELMODULE,ONLY:NWELLS,MXWELL,NWELVL,IWELCB,IPRWEL,NPWEL,
     1                       IWELPB,NNPWEL,WELAUX,WELL
      USE GWTBCTMODULE, ONLY: ICBUND,CONC,MSUMT,VBVLT,VBNMT,IPCBFLAG,
     1  IBCTCB,MCOMPT
      USE CLN1MODULE, ONLY: ICLNMB,NCLNNDS
C
      CHARACTER*16 TEXT
      DOUBLE PRECISION RATIN,RATOUT,QQ,BIG,CDIFF,QF
      DATA TEXT /'  WELL MASS FLUX'/
C     ------------------------------------------------------------------
C
C1------CLEAR RATIN AND RATOUT ACCUMULATORS, AND SET CELL-BY-CELL
C1------BUDGET FLAG.
      ZERO=0.
      RATIN=ZERO
      RATOUT=ZERO
      IBD=0
      IF(IBCTCB.LT.0 .AND. ICBCFL.NE.0) IBD=-1
      IF(IBCTCB.GT.0) IBD=ICBCFL
      IBDLBL=0
      NAUX =NWELVL - 5
      CALL CONCIAUX(ICOMP,NAUX,WELAUX,IAUX)
C
C2-----IF CELL-BY-CELL FLOWS WILL BE SAVED AS A LIST, WRITE HEADER.
      IICLNCB=0
      IF(IBD.EQ.2) THEN
         NAUX=NWELVL-5
         IF(IAUXSV.EQ.0) NAUX=0
         NNCLNNDS=0
         IF(INCLN.GT.0) THEN
           IICLNCB=ICLNMB
           NNCLNNDS=NCLNNDS
         ENDIF
        CALL UBDSVHDR(IUNSTR,KSTP,KPER,IOUT,IBCTCB,IICLNCB,NODES,
     1    NNCLNNDS,NCOL,NROW,NLAY,NWELLS,NWELVL,NAUX,IBOUND,
     2    TEXT,WELAUX,DELT,PERTIM,TOTIM,WELL)
      END IF
C
C3------CLEAR THE BUFFER.
      DO 50 N=1,NODES
      BUFF(N)=ZERO
50    CONTINUE
C
C4------IF THERE ARE NO WELLS, DO NOT ACCUMULATE (ALSO FOR IMMOBILE COMPONENTS
      IF(NWELLS.EQ.0) GO TO 200
      IF(ICOMP.GT.MCOMPT) GO TO 200
C
C5------LOOP THROUGH EACH WELL CALCULATING MASS FLUX
      DO 100 L=1,NWELLS
C
C6-----GET NODE NUMBER OF CELL CONTAINING WELL.
      N=WELL(1,L)
      QQ=ZERO
      Q = 0.0
C
C7-----IF THE CELL IS INACTIVE OR PCB, IGNORE IT.
      IF(ICBUND(N).EQ.0.OR.IPCBFLAG(N).EQ.1)GO TO 99
C
C8-----COMPUTE MASS FLUX AT WELL BOUNDARY
        QF = WELL(NWELVL,L)
C
C9-----IF THE CELL IS OUTFLOW, MASS IS Q * CONC
        IF(QF.LE.0.0)THEN
          QQ = QF * CONC(N,ICOMP)
        ELSE
C10-------IF THE CELL IS INFLOW, MASS IS Q*Co
          Co = 0.0
          IF(IAUX.NE.0) Co=WELL(4+IAUX,L)
          QQ = QF * Co
        ENDIF
        Q = QQ
C
C11-----PRINT FLOW RATE IF REQUESTED.
        IF(IBD.LT.0) THEN
          IF(IBDLBL.EQ.0) WRITE(IOUT,61) TEXT,KPER,KSTP
   61     FORMAT(1X,/1X,A,'   PERIOD ',I4,'   STEP ',I3)
          IF(IUNSTR.EQ.0)THEN
            IL = (N-1) / (NCOL*NROW) + 1
            IJ = N - (IL-1)*NCOL*NROW
            IR = (IJ-1)/NCOL + 1
            IC = IJ - (IR-1)*NCOL
            WRITE(IOUT,62) L,IL,IR,IC,Q
   62       FORMAT(1X,'WEL  ',I6,'   LAYER ',I3,'   ROW ',I5,'   COL ',
     1       I5, '   FLUX ',1PG15.6)
          ELSE
            WRITE(IOUT,63) L,N,Q
   63       FORMAT(1X,'WEL  ',I6,'    NODE ',I8,'   FLUX ',1PG15.6)
          ENDIF
          IBDLBL=1
        END IF
C
C12-----ADD FLOW RATE TO BUFFER.
        BUFF(N)=BUFF(N)+QQ
C
C13-----SEE IF FLUX IS POSITIVE OR NEGATIVE.
        IF(QQ.GE.ZERO) THEN
C
C14-------FLOW RATE IS POSITIVE (RECHARGE). ADD IT TO RATIN.
          RATIN=RATIN+QQ
        ELSE
C
C15-------FLOW RATE IS NEGATIVE (DISCHARGE). ADD IT TO RATOUT.
          RATOUT=RATOUT-QQ
        END IF
C
C16-----IF SAVING CELL-BY-CELL FLOWS IN A LIST, WRITE FLOW.  ALSO
C16-----COPY FLOW TO WELL LIST.

   99 CONTINUE
      IF(IBD.EQ.2)THEN
        CALL UBDSVREC(IUNSTR,N,NODES,NNCLNNDS,IBCTCB,IICLNCB,NWELVL,
     1    5,NAUX,Q,WELL(:,L),IBOUND,NCOL,NROW,NLAY)
      ENDIF
100   CONTINUE
C
C17------IF CELL-BY-CELL FLOWS WILL BE SAVED AS A 3-D ARRAY,
C17------CALL UBUDSV TO SAVE THEM.
      IF(IUNSTR.EQ.0)THEN
        IF(IBD.EQ.1)CALL UBUDSV(KSTP,KPER,TEXT,IBCTCB,BUFF(1),NCOL,NROW,
     1                          NLAY,IOUT)
      ELSE
        IF(IBD.EQ.1) CALL UBUDSVU(KSTP,KPER,TEXT,IBCTCB,BUFF(1),NODES,
     1                          IOUT,PERTIM,TOTIM)
      ENDIF
      IF(IBD.EQ.1.AND.INCLN.GT.0)THEN
        IF(IICLNCB.GT.0) CALL UBUDSVU(KSTP,KPER,TEXT,IICLNCB,
     1    BUFF(NODES+1),NCLNNDS,IOUT,PERTIM,TOTIM)
      ENDIF

C
C18------MOVE RATES, VOLUMES & LABELS INTO ARRAYS FOR PRINTING.
  200 RIN=RATIN
      ROUT=RATOUT
      VBVLT(3,MSUMT,ICOMP)=RIN
      VBVLT(4,MSUMT,ICOMP)=ROUT
      VBVLT(1,MSUMT,ICOMP)=VBVLT(1,MSUMT,ICOMP)+RATIN*DELT
      VBVLT(2,MSUMT,ICOMP)=VBVLT(2,MSUMT,ICOMP)+RATOUT*DELT
      VBNMT(MSUMT,ICOMP)=TEXT
C
C19------INCREMENT BUDGET TERM COUNTER(MSUM).
      MSUMT=MSUMT+1
C
C20------RETURN
      RETURN
      END
C----------------------------------------------------------------------------
      SUBROUTINE GWT2GHB1FM(ICOMP)
C     ******************************************************************
C     FORMULATE GHB TRANSPORT BOUNDARY CONDITION
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,       ONLY:IBOUND,RHS,IA,JA,AMAT
      USE GWFGHBMODULE, ONLY:NBOUND,BNDS,NGHBVL,GHBAUX
      USE GWTBCTMODULE, ONLY: ICBUND
      DOUBLE PRECISION Co,Q
C     ------------------------------------------------------------------
C
C1------IF NUMBER OF GHBS <= 0 THEN RETURN.
      IF(NBOUND.LE.0) RETURN
      BIG = 1.0E20
      NAUX=NGHBVL-6
      CALL CONCIAUX(ICOMP,NAUX,GHBAUX,IAUX)
C
C2------PROCESS EACH GHB IN THE LIST.
      DO 100 L=1,NBOUND
      N=BNDS(1,L)
C
C2A-----IF THE CELL IS INACTIVE THEN BYPASS PROCESSING.
      IF(ICBUND(N).EQ.0) GO TO 100
        Q = BNDS(NGHBVL,L)
C
C2B-----IF THE CELL IS OUTFLOW, PUT Q ON LHS DIAGONAL
        IF(Q.LE.0.0)THEN
          AMAT(IA(N)) = AMAT(IA(N)) + Q
        ELSE
C2C-------IF THE CELL IS INFLOW, PUT Q*Co ON RHS
          Co = 0.0
          IF(IAUX.NE.0) Co= BNDS(5+IAUX,L)
          RHS(N) = RHS(N) -Co * Q
        ENDIF
C
  100 CONTINUE
C
C3------RETURN
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE GWT2GHB1BD(KSTP,KPER,ICOMP)
C     ******************************************************************
C     CALCULATE MASS BUDGET FOR GHB TRANSPORT BOUNDARY CONDITIONS
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,   ONLY:IOUT,NCOL,NROW,NLAY,IBOUND,BUFF,NODES,IUNSTR,
     1              AMAT,IA,JA,INCLN
      USE GWFBASMODULE,ONLY:MSUM,ICBCFL,IAUXSV,DELT,PERTIM,TOTIM
      USE GWFGHBMODULE,ONLY:NBOUND,IGHBCB,BNDS,NGHBVL,GHBAUX
      USE GWTBCTMODULE, ONLY: ICBUND,CONC,MSUMT,VBVLT,VBNMT,IPCBFLAG,
     1  IBCTCB,MCOMPT
      USE CLN1MODULE, ONLY: ICLNMB,NCLNNDS
C
      CHARACTER*16 TEXT
      DOUBLE PRECISION RATIN,RATOUT,QQ,BIG,CDIFF
      DATA TEXT /'   GHB MASS FLUX'/
C     ------------------------------------------------------------------
C
C1------CLEAR RATIN AND RATOUT ACCUMULATORS, AND SET CELL-BY-CELL
C1------BUDGET FLAG.
      ZERO=0.
      RATIN=ZERO
      RATOUT=ZERO
      IBD=0
      IF(IBCTCB.LT.0 .AND. ICBCFL.NE.0) IBD=-1
      IF(IBCTCB.GT.0) IBD=ICBCFL
      IBDLBL=0
      NAUX=NGHBVL-6
      CALL CONCIAUX(ICOMP,NAUX,GHBAUX,IAUX)
C
C2-----IF CELL-BY-CELL FLOWS WILL BE SAVED AS A LIST, WRITE HEADER.
      IICLNCB=0
      IF(IBD.EQ.2) THEN
         NAUX=NGHBVL-6
         IF(IAUXSV.EQ.0) NAUX=0
         NNCLNNDS=0
         IF(INCLN.GT.0) THEN
           IICLNCB=ICLNMB
           NNCLNNDS=NCLNNDS
         ENDIF
         CALL UBDSVHDR(IUNSTR,KSTP,KPER,IOUT,IBCTCB,IICLNCB,NODES,
     1    NNCLNNDS,NCOL,NROW,NLAY,NBOUND,NGHBVL,NAUX,IBOUND,
     2    TEXT,GHBAUX,DELT,PERTIM,TOTIM,BNDS)
      END IF
C
C3------CLEAR THE BUFFER.
      DO 50 N=1,NODES
      BUFF(N)=ZERO
50    CONTINUE
C
C4------IF THERE ARE NO GHBs, DO NOT ACCUMULATE
      IF(NBOUND.EQ.0) GO TO 200
      IF(ICOMP.GT.MCOMPT) GO TO 200
C
C5------LOOP THROUGH EACH GHB CALCULATING MASS FLUX
      DO 100 L=1,NBOUND
C
C6-----GET NODE NUMBER OF CELL CONTAINING GHB.
      N=BNDS(1,L)
      QQ=ZERO
      Q = 0.0
C
C7-----IF THE CELL IS INACTIVE OR PCB, IGNORE IT.
      IF(ICBUND(N).EQ.0.OR.IPCBFLAG(N).EQ.1)GO TO 99
C
C8-----COMPUTE MASS FLUX AT GHB BOUNDARY
        QF = BNDS(NGHBVL,L)
C
C9-----IF THE CELL IS OUTFLOW, MASS IS Q * CONC
        IF(QF.LE.0.0)THEN
          QQ = QF * CONC(N,ICOMP)
        ELSE
C10-------IF THE CELL IS INFLOW, MASS IS Q*Co
          Co = 0.0
          IF(IAUX.NE.0) Co = BNDS(5+IAUX,L)
          QQ = QF * Co
        ENDIF
        Q = QQ
C
C11-----PRINT FLOW RATE IF REQUESTED.
        IF(IBD.LT.0) THEN
          IF(IBDLBL.EQ.0) WRITE(IOUT,61) TEXT,KPER,KSTP
   61     FORMAT(1X,/1X,A,'   PERIOD ',I4,'   STEP ',I3)
          IF(IUNSTR.EQ.0)THEN
            IL = (N-1) / (NCOL*NROW) + 1
            IJ = N - (IL-1)*NCOL*NROW
            IR = (IJ-1)/NCOL + 1
            IC = IJ - (IR-1)*NCOL
            WRITE(IOUT,62) L,IL,IR,IC,Q
   62       FORMAT(1X,'GHB  ',I6,'   LAYER ',I3,'   ROW ',I5,'   COL ',
     1       I5, '   FLUX ',1PG15.6)
          ELSE
            WRITE(IOUT,63) L,N,Q
   63       FORMAT(1X,'GHB  ',I6,'    NODE ',I8,'   FLUX ',1PG15.6)
          ENDIF
          IBDLBL=1
        END IF
C
C12-----ADD FLOW RATE TO BUFFER.
        BUFF(N)=BUFF(N)+Q
C
C13-----SEE IF FLUX IS POSITIVE OR NEGATIVE.
        IF(QQ.GE.ZERO) THEN
C
C14-------FLOW RATE IS POSITIVE (RECHARGE). ADD IT TO RATIN.
          RATIN=RATIN+QQ
        ELSE
C
C15-------FLOW RATE IS NEGATIVE (DISCHARGE). ADD IT TO RATOUT.
          RATOUT=RATOUT-QQ
        END IF
C
C16-----IF SAVING CELL-BY-CELL FLOWS IN A LIST, WRITE FLOW.  ALSO
C16-----COPY FLOW TO GHB LIST.

   99 CONTINUE
      IF(IBD.EQ.2)THEN
        CALL UBDSVREC(IUNSTR,N,NODES,NNCLNNDS,IBCTCB,IICLNCB,NGHBVL,
     1    6,NAUX,Q,BNDS(:,L),IBOUND,NCOL,NROW,NLAY)
      ENDIF

100   CONTINUE
C
C17------IF CELL-BY-CELL FLOWS WILL BE SAVED AS A 3-D ARRAY,
C17------CALL UBUDSV TO SAVE THEM.
      IF(IUNSTR.EQ.0)THEN
        IF(IBD.EQ.1)CALL UBUDSV(KSTP,KPER,TEXT,IBCTCB,BUFF(1),NCOL,NROW,
     1                          NLAY,IOUT)
      ELSE
        IF(IBD.EQ.1) CALL UBUDSVU(KSTP,KPER,TEXT,IBCTCB,BUFF(1),NODES,
     1                          IOUT,PERTIM,TOTIM)
      ENDIF
      IF(IBD.EQ.1.AND.INCLN.GT.0)THEN
        IF(IICLNCB.GT.0) CALL UBUDSVU(KSTP,KPER,TEXT,IICLNCB,
     1    BUFF(NODES+1),NCLNNDS,IOUT,PERTIM,TOTIM)
      ENDIF
C
C18------MOVE RATES, VOLUMES & LABELS INTO ARRAYS FOR PRINTING.
  200 RIN=RATIN
      ROUT=RATOUT
      VBVLT(3,MSUMT,ICOMP)=RIN
      VBVLT(4,MSUMT,ICOMP)=ROUT
      VBVLT(1,MSUMT,ICOMP)=VBVLT(1,MSUMT,ICOMP)+RATIN*DELT
      VBVLT(2,MSUMT,ICOMP)=VBVLT(2,MSUMT,ICOMP)+RATOUT*DELT
      VBNMT(MSUMT,ICOMP)=TEXT
C
C19------INCREMENT BUDGET TERM COUNTER(MSUM).
      MSUMT=MSUMT+1
C
C20------RETURN
      RETURN
      END
C----------------------------------------------------------------------------
      SUBROUTINE GWT2DRN1FM(ICOMP)
C     ******************************************************************
C     FORMULATE DRAIN TRANSPORT BOUNDARY CONDITION
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,       ONLY:IBOUND,RHS,IA,JA,AMAT
      USE GWFDRNMODULE, ONLY:NDRAIN,DRAI,NDRNVL
      USE GWTBCTMODULE, ONLY: ICBUND
      DOUBLE PRECISION Co,Q
C     ------------------------------------------------------------------
C
C1------IF NUMBER OF DRAINS <= 0 THEN RETURN.
      IF(NDRAIN.LE.0) RETURN
      BIG = 1.0E20
C
C2------PROCESS EACH DRAIN IN THE LIST.
      DO 100 L=1,NDRAIN
      N=DRAI(1,L)
C
C2A-----IF THE CELL IS INACTIVE THEN BYPASS PROCESSING.
      IF(ICBUND(N).EQ.0) GO TO 100
        Q = DRAI(NDRNVL,L)
C
C2B-----IDRAIN IS ALWAYS OUTFLOW OUTFLOW, PUT Q ON LHS DIAGONAL
        AMAT(IA(N)) = AMAT(IA(N)) + Q
C
  100 CONTINUE
C
C3------RETURN
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE GWT2DRN1BD(KSTP,KPER,ICOMP)
C     ******************************************************************
C     CALCULATE MASS BUDGET FOR DRAIN TRANSPORT BOUNDARY CONDITIONS
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,   ONLY:IOUT,NCOL,NROW,NLAY,IBOUND,BUFF,NODES,IUNSTR,
     1              AMAT,IA,JA,INCLN
      USE GWFBASMODULE,ONLY:MSUM,ICBCFL,IAUXSV,DELT,PERTIM,TOTIM
      USE GWFDRNMODULE,ONLY:NDRAIN,IDRNCB,DRAI,NDRNVL,DRNAUX
      USE GWTBCTMODULE, ONLY: ICBUND,CONC,MSUMT,VBVLT,VBNMT,IPCBFLAG,
     1  IBCTCB,MCOMPT
      USE CLN1MODULE, ONLY: ICLNMB,NCLNNDS
C
      CHARACTER*16 TEXT
      DOUBLE PRECISION RATIN,RATOUT,QQ,BIG,CDIFF
      DATA TEXT /'   DRN MASS FLUX'/
C     ------------------------------------------------------------------
C
C1------CLEAR RATIN AND RATOUT ACCUMULATORS, AND SET CELL-BY-CELL
C1------BUDGET FLAG.
      ZERO=0.
      RATIN=ZERO
      RATOUT=ZERO
      IBD=0
      IF(IBCTCB.LT.0 .AND. ICBCFL.NE.0) IBD=-1
      IF(IBCTCB.GT.0) IBD=ICBCFL
      IBDLBL=0
C
C2-----IF CELL-BY-CELL FLOWS WILL BE SAVED AS A LIST, WRITE HEADER.
      IICLNCB=0
      IF(IBD.EQ.2) THEN
        NAUX=NDRNVL-6
        IF(IAUXSV.EQ.0) NAUX=0
        NNCLNNDS=0
        IF(INCLN.GT.0) THEN
          IICLNCB=ICLNMB
          NNCLNNDS=NCLNNDS
        ENDIF
         CALL UBDSVHDR(IUNSTR,KSTP,KPER,IOUT,IBCTCB,IICLNCB,NODES,
     1    NNCLNNDS,NCOL,NROW,NLAY,NDRAIN,NDRNVL,NAUX,IBOUND,
     2    TEXT,DRNAUX,DELT,PERTIM,TOTIM,DRAI)
      END IF
C
C3------CLEAR THE BUFFER.
      DO 50 N=1,NODES
      BUFF(N)=ZERO
50    CONTINUE
C
C4------IF THERE ARE NO DRNs, DO NOT ACCUMULATE
      IF(NDRAIN.EQ.0) GO TO 200
      IF(ICOMP.GT.MCOMPT) GO TO 200
C
C5------LOOP THROUGH EACH DRN CALCULATING MASS FLUX
      DO 100 L=1,NDRAIN
C
C6-----GET NODE NUMBER OF CELL CONTAINING DRN.
      N=DRAI(1,L)
      QQ=ZERO
      Q = 0.0
C
C7-----IF THE CELL IS INACTIVE OR PCB, IGNORE IT.
      IF(ICBUND(N).EQ.0.OR.IPCBFLAG(N).EQ.1)GO TO 99
C
C8-----COMPUTE MASS FLUX AT DRAIN BOUNDARY
        QF = DRAI(NDRNVL,L)
C
C9-----DRAIN IS ALWAYS OUTFLOW, MASS IS Q * CONC
        QQ = QF * CONC(N,ICOMP)
         Q = QQ
C
C10-----PRINT FLOW RATE IF REQUESTED.
        IF(IBD.LT.0) THEN
          IF(IBDLBL.EQ.0) WRITE(IOUT,61) TEXT,KPER,KSTP
   61     FORMAT(1X,/1X,A,'   PERIOD ',I4,'   STEP ',I3)
          IF(IUNSTR.EQ.0)THEN
            IL = (N-1) / (NCOL*NROW) + 1
            IJ = N - (IL-1)*NCOL*NROW
            IR = (IJ-1)/NCOL + 1
            IC = IJ - (IR-1)*NCOL
            WRITE(IOUT,62) L,IL,IR,IC,Q
   62       FORMAT(1X,'DRN  ',I6,'   LAYER ',I3,'   ROW ',I5,'   COL ',
     1       I5, '   FLUX ',1PG15.6)
          ELSE
            WRITE(IOUT,63) L,N,Q
   63       FORMAT(1X,'DRN  ',I6,'    NODE ',I8,'   FLUX ',1PG15.6)
          ENDIF
          IBDLBL=1
        END IF
C
C11-----ADD FLOW RATE TO BUFFER.
        BUFF(N)=BUFF(N)+Q
C
C12-----SEE IF FLUX IS POSITIVE OR NEGATIVE.
        IF(QQ.GE.ZERO) THEN
C
C13-------FLOW RATE IS POSITIVE (RECHARGE). ADD IT TO RATIN.
          RATIN=RATIN+QQ
        ELSE
C
C14-------FLOW RATE IS NEGATIVE (DISCHARGE). ADD IT TO RATOUT.
          RATOUT=RATOUT-QQ
        END IF
C
C15-----IF SAVING CELL-BY-CELL FLOWS IN A LIST, WRITE FLOW.  ALSO
C15-----COPY FLOW TO DRN LIST.

   99 CONTINUE
      IF(IBD.EQ.2)THEN
        CALL UBDSVREC(IUNSTR,N,NODES,NNCLNNDS,IBCTCB,IICLNCB,NDRNVL,
     1    6,NAUX,Q,DRAI(:,L),IBOUND,NCOL,NROW,NLAY)
      ENDIF

100   CONTINUE
C
C16------IF CELL-BY-CELL FLOWS WILL BE SAVED AS A 3-D ARRAY,
C16------CALL UBUDSV TO SAVE THEM.
      IF(IUNSTR.EQ.0)THEN
        IF(IBD.EQ.1)CALL UBUDSV(KSTP,KPER,TEXT,IBCTCB,BUFF(1),NCOL,NROW,
     1                          NLAY,IOUT)
      ELSE
        IF(IBD.EQ.1) CALL UBUDSVU(KSTP,KPER,TEXT,IBCTCB,BUFF(1),NODES,
     1                          IOUT,PERTIM,TOTIM)
      ENDIF
      IF(IBD.EQ.1.AND.INCLN.GT.0)THEN
        IF(IICLNCB.GT.0) CALL UBUDSVU(KSTP,KPER,TEXT,IICLNCB,
     1    BUFF(NODES+1),NCLNNDS,IOUT,PERTIM,TOTIM)
      ENDIF

C
C17------MOVE RATES, VOLUMES & LABELS INTO ARRAYS FOR PRINTING.
  200 RIN=RATIN
      ROUT=RATOUT
      VBVLT(3,MSUMT,ICOMP)=RIN
      VBVLT(4,MSUMT,ICOMP)=ROUT
      VBVLT(1,MSUMT,ICOMP)=VBVLT(1,MSUMT,ICOMP)+RATIN*DELT
      VBVLT(2,MSUMT,ICOMP)=VBVLT(2,MSUMT,ICOMP)+RATOUT*DELT
      VBNMT(MSUMT,ICOMP)=TEXT
C
C18------INCREMENT BUDGET TERM COUNTER(MSUM).
      MSUMT=MSUMT+1
C
C19------RETURN
      RETURN
      END
C----------------------------------------------------------------------------
      SUBROUTINE GWT2DRT1FM(ICOMP)
C     ******************************************************************
C     FORMULATE DRT TRANSPORT BOUNDARY CONDITION
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,       ONLY:IBOUND,RHS,IA,JA,AMAT
      USE GWFDRTMODULE, ONLY:NDRTCL,DRTF,NDRTVL,IDRTFL,DRTAUX,
     1    IDCHANGEC,IDCHNGTYP
      USE GWTBCTMODULE, ONLY: ICBUND,CONC
      DOUBLE PRECISION Co,Q,QIN,MASSIN
C     ------------------------------------------------------------------
C
C1------IF NUMBER OF DRT CELLS <= 0 THEN RETURN.
      IF(NDRTCL.LE.0) RETURN
      BIG = 1.0E20
      NAUX=NDRTVL-5-2-IDRTFL
      CALL CONCIAUX(ICOMP,NAUX,DRTAUX,IAUX)
C
C2------PROCESS EACH DRAIN-RETURN CELL IN THE LIST.
      DO 100 L=1,NDRTCL
      ND=DRTF(1,L)
C
C3------IF THE CELL IS INACTIVE THEN BYPASS PROCESSING.
      IF(ICBUND(ND).EQ.0) GO TO 100
        Q = DRTF(NDRTVL,L)
C
C4--------DRAIN IS OUTFLOW OUTFLOW, PUT Q ON LHS DIAGONAL
        AMAT(IA(ND)) = AMAT(IA(ND)) + Q
C
C5--------TAKE CARE OF RETURN FLOW TERM
        IF (IDRTFL.GT.0) THEN
          INR = DRTF(6,L)
          IF (INR.NE.0) THEN
            IF (ICBUND(INR) .GT. 0) THEN
              Co = CONC(ND,ICOMP)
              QIN = DRTF(NDRTVL-1,L)
              MASSIN = QIN * Co
C6------------MAKE CHANGE TO SOLUTE BEFORE RETURNING FLOW
              IF(IDCHANGEC.NE.0.AND.IAUX.NE.0)THEN
                IF(IDCHNGTYP(L).EQ.1) MASSIN = MASSIN+QIN*DRTF(9+IAUX,L)
                IF(IDCHNGTYP(L).EQ.2) MASSIN = MASSIN + DRTF(9+IAUX,L)
                IF(IDCHNGTYP(L).EQ.3) MASSIN = QIN*DRTF(9+IAUX,L)
                IF(IDCHNGTYP(L).EQ.4) MASSIN = DRTF(9+IAUX,L)
              ENDIF
C7------------PUT NET RETURN FLOW MASS ON RHS OF EQUATION
              RHS(INR) = RHS(INR) - MASSIN
            ENDIF
          ENDIF
        ENDIF
C
  100 CONTINUE
C
C8------RETURN
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE GWT2DRT1BD(KSTP,KPER,ICOMP)
C     ******************************************************************
C     CALCULATE MASS BUDGET FOR DRT TRANSPORT BOUNDARY CONDITIONS
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,   ONLY:IOUT,NCOL,NROW,NLAY,IBOUND,BUFF,NODES,IUNSTR,
     1              AMAT,IA,JA,INCLN
      USE GWFBASMODULE,ONLY:MSUM,ICBCFL,IAUXSV,DELT,PERTIM,TOTIM
      USE GWFDRTMODULE,ONLY:NDRTCL,IDRTCB,DRTF,NDRTVL,DRTAUX,IDRTFL,
     1                      NRFLOW,IDCHANGEC,IDCHNGTYP
      USE GWTBCTMODULE, ONLY: ICBUND,CONC,MSUMT,VBVLT,VBNMT,IPCBFLAG,
     1  IBCTCB,MCOMPT
      USE CLN1MODULE, ONLY: ICLNMB,NCLNNDS,ICLNCB
C
      CHARACTER*16 TEXT
      DOUBLE PRECISION RATIN,RATOUT,QQ,BIG,CDIFF,QFIN,QQIN,QF
      DATA TEXT /'   DRT MASS FLUX'/
C     ------------------------------------------------------------------
C
C1------CLEAR RATIN AND RATOUT ACCUMULATORS, AND SET CELL-BY-CELL
C1------BUDGET FLAG.
      ZERO=0.
      RATIN=ZERO
      RATOUT=ZERO
      IBD=0
      IF(IBCTCB.LT.0 .AND. ICBCSPCFL.NE.0) IBD=-1
      IF(IBCTCB.GT.0) IBD=ICBCFL
      IBDLBL=0
      NAUX=NDRTVL-5-2-IDRTFL
      CALL CONCIAUX(ICOMP,NAUX,DRTAUX,IAUX)
C
C2-----IF CELL-BY-CELL FLOWS WILL BE SAVED AS A LIST, WRITE HEADER.
      IICLNCB=0
      IF(IBD.EQ.2) THEN
         NAUX=NDRTVL-5-2
         IF(IAUXSV.EQ.0) NAUX=0
         NNCLNNDS=0
         IF(INCLN.GT.0) THEN
           IICLNCB=ICLNMB
           NNCLNNDS=NCLNNDS
         ENDIF
         CALL UBDSVHDR(IUNSTR,KSTP,KPER,IOUT,IBCTCB,IICLNCB,NODES,
     1    NNCLNNDS,NCOL,NROW,NLAY,NBOUND,NDRTVL,NAUX,IBOUND,
     2    TEXT,DRTAUX,DELT,PERTIM,TOTIM,DRTF)
C         IF(IUNSTR.EQ.0) THEN
C           CALL UBDSV4(KSTP,KPER,TEXT,NAUX,DRTAUX,IBCTCB,NCOL,NROW,NLAY,
C     1          NDRTCL+NRFLOW,IOUT,DELT,PERTIM,TOTIM,ICBUND)
C         ELSE
C           CALL UBDSV4U(KSTP,KPER,TEXT,NAUX,DRTAUX,ICLNCB,NEQS,
C     1          NDRTCL+NRFLOW,IOUT,DELT,PERTIM,TOTIM,ICBUND)
C         ENDIF
      END IF
C
C3------CLEAR THE BUFFER.
      DO 50 N=1,NODES
      BUFF(N)=ZERO
50    CONTINUE
C
C4------IF THERE ARE NO DRTs, DO NOT ACCUMULATE
      IF(NDRTCL.EQ.0) GO TO 200
      IF(ICOMP.GT.MCOMPT) GO TO 200
C
C5------LOOP THROUGH EACH DRN CALCULATING MASS FLUX
      DO 100 L=1,NDRTCL
C
C6-----GET NODE NUMBER OF CELL CONTAINING DRN.
      ND=DRTF(1,L)
      QQ=ZERO
      Q = 0.0
C
C7-----IF THE CELL IS INACTIVE OR PCB, IGNORE IT.
      IF(ICBUND(ND).EQ.0.OR.IPCBFLAG(ND).EQ.1)GO TO 99
C
C8-----COMPUTE MASS FLUX AT DRAIN BOUNDARY
        QF = DRTF(NDRTVL,L)
C
C9-----DRAIN IS OUTFLOW, MASS IS Q * CONC
        Co = CONC(ND,ICOMP)
        QQ = QF * Co
        Q = QQ
        RATOUT=RATOUT-QQ
C
C10-----COMPUTE MASS FLUX AT RETURN FLOW NODE BOUNDARY
        IF (IDRTFL.GT.0) THEN
          INR = DRTF(6,L)
          IF (INR.NE.0) THEN
            IF (ICBUND(INR) .GT. 0) THEN
              QFIN = DRTF(NDRTVL-1,L)
              QQIN = QFIN * Co
C9B-----------MAKE CHANGE TO SOLUTE BEFORE RETURNING FLOW
              IF(IDCHANGEC.NE.0.AND.IAUX.NE.0)THEN
                IF(IDCHNGTYP(L).EQ.1) QQIN = QQIN + QFIN*DRTF(9+IAUX,L)
                IF(IDCHNGTYP(L).EQ.2) QQIN = QQIN + DRTF(9+IAUX,L)
                IF(IDCHNGTYP(L).EQ.3) QQIN = QFIN*DRTF(9+IAUX,L)
                IF(IDCHNGTYP(L).EQ.4) QQIN = DRTF(9+IAUX,L)
                IF(IDCHNGTYP(L).EQ.5) THEN
                   IF(Co .GT. DRTF(9+IAUX,L)) QQIN = QFIN*DRTF(9+IAUX,L)
                ENDIF   
              ENDIF
              QIN = QQIN
              RATIN=RATIN+QQIN
            ENDIF
          ENDIF
        ENDIF
C ------------------------------------------------------------------------
C11-----PRINT FLOW RATE IF REQUESTED.
        IF(IBD.LT.0) THEN
          IF(IBDLBL.EQ.0) WRITE(IOUT,61) TEXT,KPER,KSTP
   61     FORMAT(1X,/1X,A,'   PERIOD ',I4,'   STEP ',I3)
C11A--------PRINT FOR STRUCTURED GRID
          IF(IUNSTR.EQ.0)THEN
            IL = (ND-1) / (NCOL*NROW) + 1
            IJ = ND - (IL-1)*NCOL*NROW
            IR = (IJ-1)/NCOL + 1
            IC = IJ - (IR-1)*NCOL
            WRITE(IOUT,62) L,IL,IR,IC,Q
   62       FORMAT(1X,'DRT  ',I6,'   LAYER ',I3,'   ROW ',I5,'   COL ',
     1       I5, '   FLUX ',1PG15.6)
C11B-----PRINT FOR RETURN FLOW
            IF(INR.NE.0)THEN
              IL = (INR-1) / (NCOL*NROW) + 1
              IJ = INR - (IL-1)*NCOL*NROW
              IR = (IJ-1)/NCOL + 1
              IC = IJ - (IR-1)*NCOL
              WRITE(IOUT,64) L,IL,IR,IC,QIN
   64         FORMAT(1X,'RETURN ',I6,'   LAYER ',I3,'   ROW ',I5,
     1         '   COL ',I5, '   FLUX ',1PG15.6)
            ENDIF
C11C--------PRINT FOR UNSTRUCTURED GRID
          ELSE
            WRITE(IOUT,63) L,ND,Q
   63       FORMAT(1X,'DRT  ',I6,'    NODE ',I8,'   FLUX ',1PG15.6)
          ENDIF
C11D-----PRINT FOR RETURN FLOW
          IF(INR.NE.0)THEN
            WRITE(IOUT,65) L,INR,QIN
   65       FORMAT(1X,'DRT  ',I6,'    NODE ',I8,'   FLUX ',1PG15.6)
          ENDIF
          IBDLBL=1
        END IF
C ------------------------------------------------------------------------
C
C11-----ADD FLOW RATE TO BUFFER.
        BUFF(ND)=BUFF(ND)+Q
        IF(INR.NE.0) BUFF(INR)=BUFF(INR)+QIN
C
C15-----IF SAVING CELL-BY-CELL FLOWS IN A LIST, WRITE FLOW.  ALSO
C15-----COPY FLOW TO DRN LIST.

   99 CONTINUE
      IF(IBD.EQ.2)THEN
        CALL UBDSVREC(IUNSTR,N,NODES,NNCLNNDS,IBCTCB,IICLNCB,NDRTVL,
     1    6,NAUX,Q,DRTF(:,L),IBOUND,NCOL,NROW,NLAY)
        IF (INR.NE.0) THEN
          CALL UBDSVREC(IUNSTR,N,NODES,NNCLNNDS,IBCTCB,IICLNCB,NDRTVL,
     1    6,NAUX,QIN,DRTF(:,L),IBOUND,NCOL,NROW,NLAY)
        ENDIF
      ENDIF
C      IF(IBD.EQ.2) THEN
C        IF(IUNSTR.EQ.0)THEN
C          CALL UBDSVB(IBCTCB,NCOL,NROW,IC,IR,IL,Q,
C     1                  DRTF(1,L),NDRTVL,NAUX,10,ICBUND,NLAY)
C          IF(INR.NE.0)THEN
C            CALL UBDSVB(IBCTCB,NCOL,NROW,ICR,IRR,ILR,QIN,
C     1                  DRTF(1,L),NDRTVL,NAUX,10,ICBUND,NLAY)
C          ENDIF
C        ELSE
C          CALL UBDSVBU(IBCTCB,NEQS,ND,Q,
C     1                  DRTF(1,L),NDRTVL,NAUX,10,ICBUND)
C          IF(INR.NE.0)THEN
C          CALL UBDSVBU(IBCTCB,NEQS,INR,QIN,
C     1                  DRTF(1,L),NDRTVL,NAUX,10,ICBUND)
C          ENDIF
C        ENDIF
C      ENDIF
100   CONTINUE
C
C16------IF CELL-BY-CELL FLOWS WILL BE SAVED AS A 3-D ARRAY,
C16------CALL UBUDSV TO SAVE THEM.
      IF(IUNSTR.EQ.0)THEN
        IF(IBD.EQ.1)CALL UBUDSV(KSTP,KPER,TEXT,IBCTCB,BUFF(1),NCOL,NROW,
     1                          NLAY,IOUT)
      ELSE
        IF(IBD.EQ.1) CALL UBUDSVU(KSTP,KPER,TEXT,IBCTCB,BUFF(1),NODES,
     1                          IOUT,PERTIM,TOTIM)
      ENDIF
      IF(IBD.EQ.1.AND.INCLN.GT.0)THEN
        IF(IICLNCB.GT.0) CALL UBUDSVU(KSTP,KPER,TEXT,IICLNCB,
     1    BUFF(NODES+1),NCLNNDS,IOUT,PERTIM,TOTIM)
      ENDIF
C
C17------MOVE RATES, VOLUMES & LABELS INTO ARRAYS FOR PRINTING.
  200 RIN=RATIN
      ROUT=RATOUT
      VBVLT(3,MSUMT,ICOMP)=RIN
      VBVLT(4,MSUMT,ICOMP)=ROUT
      VBVLT(1,MSUMT,ICOMP)=VBVLT(1,MSUMT,ICOMP)+RATIN*DELT
      VBVLT(2,MSUMT,ICOMP)=VBVLT(2,MSUMT,ICOMP)+RATOUT*DELT
      VBNMT(MSUMT,ICOMP)=TEXT
C
C18------INCREMENT BUDGET TERM COUNTER(MSUM).
      MSUMT=MSUMT+1
C
C19------RETURN
      RETURN
      END
C----------------------------------------------------------------------------
      SUBROUTINE GWT2QRT1FM(ICOMP)
C     ******************************************************************
C     FORMULATE QRT TRANSPORT BOUNDARY CONDITION
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,       ONLY:IBOUND,RHS,IA,JA,AMAT
      USE GWFQRTMODULE, ONLY:NQRTCL,QRTF,NQRTVL,IQRTFL,NodQRT,QRTFLOW,
     1    QRTAUX,IQCHANGEC,IQCHNGTYP
      USE GWTBCTMODULE, ONLY: ICBUND,CONC
      DOUBLE PRECISION Co,Q,QIN,MASSIN
C     ------------------------------------------------------------------
C
C1------IF NUMBER OF QRT CELLS <= 0 THEN RETURN.
      IF(NQRTCL.LE.0) RETURN
      BIG = 1.0E20
      NAUX=NQRTVL-5-IQRTFL
      CALL CONCIAUX(ICOMP,NAUX,QRTAUX,IAUX)
C
C2------PROCESS EACH CELL IN THE QRT LIST.
      IRT = 0                  !------------------! POINTER FOR LOCATION IN NodQRT ARRAY
      DO 100 L=1,NQRTCL
      ND=QRTF(1,L)
C
C3------IF THE CELL IS INACTIVE THEN BYPASS PROCESSING.
      IF(ICBUND(ND).EQ.0) GO TO 100
        Q = QRTF(NQRTVL,L)
C
C4--------QRT's SINK IS OUTFLOW, REMOVE Q FROM LHS DIAGONAL
        AMAT(IA(ND)) = AMAT(IA(ND)) - Q
C
C5--------TAKE CARE OF RETURN FLOW TERM
        IF (IQRTFL.GT.0) THEN
          NumRT = QRTF(5,L)
          IF (NumRT.NE.0) THEN
            DO JJ = 1,NumRT
              IRT = IRT + 1
              INR = NodQRT(IRT)
              IF (ICBUND(INR) .GT. 0) THEN
                Co = CONC(ND,ICOMP)
                QIN = QRTFLOW(IRT)
                MASSIN = QIN * Co
C6------------MAKE CHANGE TO SOLUTE BEFORE RETURNING FLOW
                IF(IQCHANGEC.NE.0.AND.IAUX.NE.0)THEN
                  IF(IQCHNGTYP(L).EQ.1) MASSIN=MASSIN+QIN*QRTF(6+IAUX,L)
                  IF(IQCHNGTYP(L).EQ.2) MASSIN = MASSIN + QRTF(6+IAUX,L)
                  IF(IQCHNGTYP(L).EQ.3) MASSIN= QIN*QRTF(6+IAUX,L)
                  IF(IQCHNGTYP(L).EQ.4) MASSIN =  QRTF(6+IAUX,L)
                ENDIF
                RHS(INR) = RHS(INR) - MASSIN
              END IF
            ENDDO
          ENDIF
        ENDIF
C
  100 CONTINUE
C
C3------RETURN
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE GWT2QRT1BD(KSTP,KPER,ICOMP)
C     ******************************************************************
C     CALCULATE MASS BUDGET FOR QRT TRANSPORT BOUNDARY CONDITIONS
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,   ONLY:IOUT,NCOL,NROW,NLAY,IBOUND,BUFF,NODES,IUNSTR,
     1              AMAT,IA,JA,INCLN
      USE GWFBASMODULE,ONLY:MSUM,ICBCFL,IAUXSV,DELT,PERTIM,TOTIM
      USE GWFQRTMODULE,ONLY:NQRTCL,IQRTCB,QRTF,NQRTVL,QRTAUX,IQRTFL,
     1                      QRTFLOW,NodQRT,IQCHANGEC,IQCHNGTYP
      USE GWTBCTMODULE, ONLY: ICBUND,CONC,MSUMT,VBVLT,VBNMT,IPCBFLAG,
     1  IBCTCB,MCOMPT
      USE CLN1MODULE, ONLY: ICLNMB,NCLNNDS
C
      CHARACTER*16 TEXT
      DOUBLE PRECISION RATIN,RATOUT,QQ,BIG,CDIFF,QFIN,QQIN,QF
      DATA TEXT /'   QRT MASS FLUX'/
C     ------------------------------------------------------------------
C
C1------CLEAR RATIN AND RATOUT ACCUMULATORS, AND SET CELL-BY-CELL
C1------BUDGET FLAG.
      ZERO=0.
      RATIN=ZERO
      RATOUT=ZERO
      IBD=0
      IF(IBCTCB.LT.0 .AND. ICBCFL.NE.0) IBD=-1
      IF(IBCTCB.GT.0) IBD=ICBCFL
      IBDLBL=0
      NAUX=NQRTVL-5-IQRTFL
      CALL CONCIAUX(ICOMP,NAUX,QRTAUX,IAUX)
C
C2-----IF CELL-BY-CELL FLOWS WILL BE SAVED AS A LIST, WRITE HEADER.
      IICLNCB=0
      IF(IBD.EQ.2) THEN
         NAUX=NQRTVL-5
         IF(IAUXSV.EQ.0) NAUX=0
         NNCLNNDS=0
         IF(INCLN.GT.0) THEN
           IICLNCB=ICLNMB
           NNCLNNDS=NCLNNDS
         ENDIF
         CALL UBDSVHDR(IUNSTR,KSTP,KPER,IOUT,IBCTCB,IICLNCB,NODES,
     1    NNCLNNDS,NCOL,NROW,NLAY,NBOUND,NQRTVL,NAUX,IBOUND,
     2    TEXT,QRTAUX,DELT,PERTIM,TOTIM,QRTF)
C
C         IF(IUNSTR.EQ.0) THEN
C           CALL UBDSV4(KSTP,KPER,TEXT,NAUX,QRTAUX,IBCTCB,NCOL,NROW,NLAY,
C     1          NQRTCL,IOUT,DELT,PERTIM,TOTIM,ICBUND)
C         ELSE
C           CALL UBDSV4U(KSTP,KPER,TEXT,NAUX,QRTAUX,IBCTCB,NEQS,
C     1          NQRTCL,IOUT,DELT,PERTIM,TOTIM,ICBUND)
C         ENDIF
      END IF
C
C3------CLEAR THE BUFFER.
      DO 50 N=1,NODES
      BUFF(N)=ZERO
50    CONTINUE
C
C4------IF THERE ARE NO DQRTs, DO NOT ACCUMULATE
      IF(NQRTCL.EQ.0) GO TO 200
      IF(ICOMP.GT.MCOMPT) GO TO 200
C
C5------LOOP THROUGH EACH QRN CALCULATING MASS FLUX
      IRT = 0   !------------------! INDEX FOR LOCATION IN NodQRT ARRA
      DO 100 L=1,NQRTCL
C
C6-----GET NODE NUMBER OF CELL CONTAINING SINK.
      ND=QRTF(1,L)
      QQ=ZERO
      Q = 0.0
C
C7-----IF THE CELL IS INACTIVE OR PCB, IGNORE IT.
      IF(ICBUND(ND).EQ.0.OR.IPCBFLAG(ND).EQ.1)GO TO 99
C
C8-----COMPUTE MASS FLUX AT SINK OF QRT BOUNDARY
        QF = QRTF(NQRTVL,L)
C
C9-----QRT SINK IS OUTFLOW, MASS IS Q * CONC
        Co = CONC(ND,ICOMP)
        QQ = -QF * Co
        Q = QQ
        RATOUT=RATOUT-QQ
C9A--------ADD FLOW RATE TO BUFFER.
        BUFF(ND)=BUFF(ND)+Q
C
C10-----COMPUTE MASS FLUX AT RETURN FLOW NODE BOUNDARY
        IF (IQRTFL.GT.0) THEN
          NumRT = QRTF(5,L)
          IF (NUMRT.NE.0) THEN

            DO JJ = 1,NumRT
              IRT = IRT + 1
              INR = NodQRT(IRT)
              IF (ICBUND(INR) .GT. 0) THEN
                QFIN = QRTFLOW(IRT)
                QQIN = QFIN * Co
C9B-----------MAKE CHANGE TO SOLUTE BEFORE RETURNING FLOW
              IF(IQCHANGEC.NE.0.AND.IAUX.NE.0)THEN
                IF(IQCHNGTYP(L).EQ.1) QQIN = QQIN + QFIN*QRTF(6+IAUX,L)
                IF(IQCHNGTYP(L).EQ.2) QQIN = QQIN + QRTF(6+IAUX,L)
                IF(IQCHNGTYP(L).EQ.3) QQIN = QFIN*QRTF(6+IAUX,L)
                IF(IQCHNGTYP(L).EQ.4) QQIN = QRTF(6+IAUX,L)
                IF(IQCHNGTYP(L).EQ.5) THEN
                   IF(Co .GT. QRTF(9+IAUX,L)) QQIN = QFIN*QRTF(9+IAUX,L)
                ENDIF                
              ENDIF
                QIN = QQIN
                RATIN=RATIN+QQIN
C10A------------ADD FLOW RATE TO BUFFER.
                BUFF(INR)=BUFF(INR)+QIN
              ENDIF
            ENDDO
          ENDIF
        ENDIF
C ------------------------------------------------------------------------
C11-----PRINT FLOW RATE IF REQUESTED.
        IF(IBD.LT.0) THEN
          IF(IBDLBL.EQ.0) WRITE(IOUT,61) TEXT,KPER,KSTP
   61     FORMAT(1X,/1X,A,'   PERIOD ',I4,'   STEP ',I3)
C11A--------PRINT FOR STRUCTURED GRID
          IF(IUNSTR.EQ.0)THEN
            IL = (ND-1) / (NCOL*NROW) + 1
            IJ = ND - (IL-1)*NCOL*NROW
            IR = (IJ-1)/NCOL + 1
            IC = IJ - (IR-1)*NCOL
            WRITE(IOUT,62) L,IL,IR,IC,Q
   62       FORMAT(1X,'QRT  ',I6,'   LAYER ',I3,'   ROW ',I5,'   COL ',
     1       I5, '   FLUX ',1PG15.6)
C11B-----PRINT FOR RETURN FLOW
          IF (NumRT.NE.0) THEN
            DO I=IRT-NUMRT+1, IRT
              IL = (INR-1) / (NCOL*NROW) + 1
              IJ = INR - (IL-1)*NCOL*NROW
              IR = (IJ-1)/NCOL + 1
              IC = IJ - (IR-1)*NCOL
              WRITE(IOUT,64) L,IL,IR,IC,QRTFLOW(I)*Co
   64         FORMAT(1X,'RETURN ',I6,'   LAYER ',I3,'   ROW ',I5,
     1         '   COL ',I5, '   FLUX ',1PG15.6)
            ENDDO
          ENDIF
C11C--------PRINT FOR UNSTRUCTURED GRID
          ELSE
            WRITE(IOUT,63) L,ND,Q
   63       FORMAT(1X,'QRT  ',I6,'    NODE ',I8,'   FLUX ',1PG15.6)
          ENDIF
C11D-----PRINT FOR RETURN FLOW
          IF (NumRT.NE.0) THEN
            DO I=IRT-NUMRT+1, IRT
              WRITE(IOUT,550) L,NodQRT(I),QRTFLOW(I)*Co
            ENDDO
  550       FORMAT(1X,'SINK ',I6,
     *        ' RETURN:  NODE ',I10,'   FLUX ',1PG15.6)
          ENDIF
          IBDLBL=1
        END IF
C ------------------------------------------------------------------------
C
C15-----IF SAVING CELL-BY-CELL FLOWS IN A LIST, WRITE FLOW.  ALSO
C15-----COPY FLOW TO DRN LIST.

   99 CONTINUE
      IF(IBD.EQ.2)THEN
        CALL UBDSVREC(IUNSTR,N,NODES,NNCLNNDS,IBCTCB,IICLNCB,NQRTVL,
     1    6,NAUX,Q,QRTF(:,L),IBOUND,NCOL,NROW,NLAY)
        IF (NumRT.NE.0) THEN
          DO I=IRT-NUMRT+1, IRT
            AMASS = QRTFLOW(I)*Co
            CALL UBDSVREC(IUNSTR,N,NODES,NNCLNNDS,IBCTCB,IICLNCB,NQRTVL,
     1    6,NAUX,AMASS,QRTF(:,L),IBOUND,NCOL,NROW,NLAY)
          ENDDO
        ENDIF
      ENDIF
C      IF(IBD.EQ.2) THEN
C        IF(IUNSTR.EQ.0)THEN
C          IL = (ND-1) / (NCOL*NROW) + 1
C          IJ = ND - (IL-1)*NCOL*NROW
C          IR = (IJ-1)/NCOL + 1
C          IC = IJ - (IR-1)*NCOL
C          CALL UBDSVB(IBCTCB,NCOL,NROW,IC,IR,IL,Q,
C     1                  QRTF(1,L),NQRTVL,NAUX,10,ICBUND,NLAY)
C          IF (NumRT.NE.0) THEN
C            DO I=IRT-NUMRT+1, IRT
C              ILR = (ND-1) / (NCOL*NROW) + 1
C              IJR = ND - (ILR-1)*NCOL*NROW
C              IRR = (IJR-1)/NCOL + 1
C              ICR = IJR - (IRR-1)*NCOL
C              CALL UBDSVB(IBCTCB,NCOL,NROW,ICR,IRR,ILR,QRTFLOW(I),
C     &           QRTF(1,L),NQRTVL,NAUX,10,ICBUND,NLAY)
C            ENDDO
C          ENDIF
C        ELSE
C          CALL UBDSVBU(IBCTCB,NEQS,ND,Q,
C     1                  QRTF(1,L),NQRTVL,NAUX,10,ICBUND)
C          IF (NumRT.NE.0) THEN
C            DO I=IRT-NUMRT+1, IRT
C              AMASS = QRTFLOW(I)*Co
C              CALL UBDSVBU(IBCTCB,NEQS,NodQRT(I),AMASS,QRTF(1,L),
C     &                NQRTVL,NAUX,10,ICBUND)
C            ENDDO
C          ENDIF
C        ENDIF
C      ENDIF
100   CONTINUE
C
C16------IF CELL-BY-CELL FLOWS WILL BE SAVED AS A 3-D ARRAY,
C16------CALL UBUDSV TO SAVE THEM.
      IF(IUNSTR.EQ.0)THEN
        IF(IBD.EQ.1)CALL UBUDSV(KSTP,KPER,TEXT,IBCTCB,BUFF(1),NCOL,NROW,
     1                          NLAY,IOUT)
      ELSE
        IF(IBD.EQ.1) CALL UBUDSVU(KSTP,KPER,TEXT,IBCTCB,BUFF(1),NODES,
     1                          IOUT,PERTIM,TOTIM)
      ENDIF
      IF(IBD.EQ.1.AND.INCLN.GT.0)THEN
        iiclncb = ICLNMB
        IF(IICLNCB.GT.0) CALL UBUDSVU(KSTP,KPER,TEXT,IICLNCB,
     1    BUFF(NODES+1),NCLNNDS,IOUT,PERTIM,TOTIM)
      ENDIF
C
C17------MOVE RATES, VOLUMES & LABELS INTO ARRAYS FOR PRINTING.
  200 RIN=RATIN
      ROUT=RATOUT
      VBVLT(3,MSUMT,ICOMP)=RIN
      VBVLT(4,MSUMT,ICOMP)=ROUT
      VBVLT(1,MSUMT,ICOMP)=VBVLT(1,MSUMT,ICOMP)+RATIN*DELT
      VBVLT(2,MSUMT,ICOMP)=VBVLT(2,MSUMT,ICOMP)+RATOUT*DELT
      VBNMT(MSUMT,ICOMP)=TEXT
C
C18------INCREMENT BUDGET TERM COUNTER(MSUM).
      MSUMT=MSUMT+1
C
C19------RETURN
      RETURN
      END
C----------------------------------------------------------------------------
      SUBROUTINE GWT2EVT1FM(ICOMP)
C     ******************************************************************
C     FORMULATE EVT TRANSPORT BOUNDARY CONDITION
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,       ONLY:IBOUND,RHS,IA,JA,AMAT,FMBE,IFMBC
      USE GWFEVTMODULE, ONLY:EVTF,INIEVT,IEVT,IEVTCB,ETFACTOR,NEVTOP
      USE GWTBCTMODULE, ONLY: ICBUND
      DOUBLE PRECISION Co,Q
C     ------------------------------------------------------------------
C
C1------IF NUMBER OF EVT CELLS <= 0 THEN RETURN.
      IF(INIEVT.LE.0) RETURN
      BIG = 1.0E20
C
C2------PROCESS EACH EVT NODE IN THE LIST.
      DO 100 L=1,INIEVT
      N=IEVT(L)
      IF(NEVTOP.EQ.3.AND.IBOUND(N).EQ.0)THEN
        CALL FIRST_ACTIVE_BELOW(N)
      ENDIF      
C
C2A-----IF THE CELL IS INACTIVE THEN BYPASS PROCESSING.
      IF(ICBUND(N).EQ.0) GO TO 100
        Q = EVTF(L)
C
C2B-----EVT IS ALWAYS OUTFLOW OUTFLOW, PUT Q ON LHS DIAGONAL
        AMAT(IA(N)) = AMAT(IA(N)) + Q * ETFACTOR(ICOMP)
C2C-----ADJUST FOR DRY CELLS WHERE ETRACTOR IS ZERO
        IF(IFMBC.NE.0)THEN
          IF(ABS(ETFACTOR(ICOMP)).LT.1.0E-5)THEN
            AMAT(IA(N)) = AMAT(IA(N)) + fmbe(n)  ! no correction on ET nodes (keeps diagonal zero for dry cells)
          ENDIF
        ENDIF
C
  100 CONTINUE
C
C3------RETURN
      RETURN
      END
C----------------------------------------------------------------------------
      SUBROUTINE GWT2ETS1FM(ICOMP)
C     ******************************************************************
C     FORMULATE ETS TRANSPORT BOUNDARY CONDITION
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,       ONLY:IBOUND,RHS,IA,JA,AMAT
      USE GWFETSMODULE, ONLY:ETSF,INIETS,IETS,IETSCB,ESFACTOR,NETSOP
      USE GWTBCTMODULE, ONLY: ICBUND
      DOUBLE PRECISION Co,Q
C     ------------------------------------------------------------------
C
C1------IF NUMBER OF ETS CELLS <= 0 THEN RETURN.
      IF(INIETS.LE.0) RETURN
      BIG = 1.0E20
C
C2------PROCESS EACH ETS NODE IN THE LIST.
      DO 100 L=1,INIETS
      N=IETS(L)
      IF(NETSOP.EQ.3.AND.IBOUND(N).EQ.0)THEN
        CALL FIRST_ACTIVE_BELOW(N)
      ENDIF      
C
C2A-----IF THE CELL IS INACTIVE THEN BYPASS PROCESSING.
      IF(ICBUND(N).EQ.0) GO TO 100
        Q = ETSF(L)
C
C2B-----ETS IS ALWAYS OUTFLOW OUTFLOW, PUT Q ON LHS DIAGONAL
        AMAT(IA(N)) = AMAT(IA(N)) + Q * ESFACTOR(ICOMP)
C
  100 CONTINUE
C
C3------RETURN
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE GWT2EVT1BD(KSTP,KPER,ICOMP)
C     ******************************************************************
C     CALCULATE MASS BUDGET FOR EVT TRANSPORT BOUNDARY CONDITIONS
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,   ONLY:IOUT,NCOL,NROW,NLAY,IBOUND,BUFF,NODES,IUNSTR,
     1              AMAT,IA,JA
      USE GWFBASMODULE,ONLY:MSUM,ICBCFL,IAUXSV,DELT,PERTIM,TOTIM
      USE GWFEVTMODULE,ONLY:EVTF,INIEVT,IEVT,IEVTCB,ETFACTOR,NEVTOP
      USE GWTBCTMODULE, ONLY: ICBUND,CONC,MSUMT,VBVLT,VBNMT,IPCBFLAG,
     1  IBCTCB,MCOMPT
C
      CHARACTER*16 TEXT
      DOUBLE PRECISION RATIN,RATOUT,QQ,BIG,CDIFF
      DATA TEXT /'   EVT MASS FLUX'/
C     ------------------------------------------------------------------
C
C1------CLEAR RATIN AND RATOUT ACCUMULATORS, AND SET CELL-BY-CELL
C1------BUDGET FLAG.
      ZERO=0.
      RATIN=ZERO
      RATOUT=ZERO
      IBD=0
      IF(IBCTCB.LT.0 .AND. ICBCFL.NE.0) IBD=-1
      IF(IBCTCB.GT.0 .AND. ICBCFL.NE.0) IBD=1        ! should be ICBCFL but IBD =2 is not yet in EVT 
      IBDLBL=0
C
C2------CLEAR THE BUFFER.
      DO 50 N=1,NODES
      BUFF(N)=ZERO
50    CONTINUE
C
C3------IF THERE ARE NO EVTs, DO NOT ACCUMULATE
      IF(INIEVT.EQ.0) GO TO 200
      IF(ICOMP.GT.MCOMPT) GO TO 200
C
C4------LOOP THROUGH EACH EVT CALCULATING MASS FLUX
      DO 100 L=1,INIEVT
C
C5-----GET NODE NUMBER OF CELL CONTAINING EVT.
      N=IEVT(L)
      IF(NEVTOP.EQ.3.AND.IBOUND(N).EQ.0)THEN
        CALL FIRST_ACTIVE_BELOW(N)
      ENDIF      
      QQ=ZERO
      Q = 0.0
C
C6-----IF THE CELL IS INACTIVE OR PCB, IGNORE IT.
      IF(ICBUND(N).EQ.0.OR.IPCBFLAG(N).EQ.1)GO TO 99
C
C7-----COMPUTE MASS FLUX AT EVT BOUNDARY
        QF = EVTF(L)
C
C8-----EVT IS ALWAYS OUTFLOW, MASS IS Q * CONC
        QQ = QF * CONC(N,ICOMP) * ETFACTOR(ICOMP)
         Q = QQ
C
C9-----PRINT FLOW RATE IF REQUESTED.
        IF(IBD.LT.0) THEN
          IF(IBDLBL.EQ.0) WRITE(IOUT,61) TEXT,KPER,KSTP
   61     FORMAT(1X,/1X,A,'   PERIOD ',I4,'   STEP ',I3)
          IF(IUNSTR.EQ.0)THEN
            IL = (N-1) / (NCOL*NROW) + 1
            IJ = N - (IL-1)*NCOL*NROW
            IR = (IJ-1)/NCOL + 1
            IC = IJ - (IR-1)*NCOL
            WRITE(IOUT,62) L,IL,IR,IC,Q
   62       FORMAT(1X,'EVT  ',I6,'   LAYER ',I3,'   ROW ',I5,'   COL ',
     1       I5, '   FLUX ',1PG15.6)
          ELSE
            WRITE(IOUT,63) L,N,Q
   63       FORMAT(1X,'EVT  ',I6,'    NODE ',I8,'   FLUX ',1PG15.6)
          ENDIF
          IBDLBL=1
        END IF
C
C10-----ADD FLOW RATE TO BUFFER.
        BUFF(N)=BUFF(N)+Q
C
C11-----SEE IF FLUX IS POSITIVE OR NEGATIVE.
        IF(QQ.GE.ZERO) THEN
C
C12-------FLOW RATE IS POSITIVE (RECHARGE). ADD IT TO RATIN.
          RATIN=RATIN+QQ
        ELSE
C
C13-------FLOW RATE IS NEGATIVE (DISCHARGE). ADD IT TO RATOUT.
          RATOUT=RATOUT-QQ
        END IF
C
C14-----IF SAVING CELL-BY-CELL FLOWS IN A LIST, WRITE FLOW.  ALSO
C14-----COPY FLOW TO EVT LIST.

   99 CONTINUE
C
100   CONTINUE
C
C15------IF CELL-BY-CELL FLOWS WILL BE SAVED AS A 3-D ARRAY,
C15------CALL UBUDSV TO SAVE THEM.
      IF(IUNSTR.EQ.0)THEN
        IF(IBD.EQ.1) CALL UBUDSV(KSTP,KPER,TEXT,IBCTCB,BUFF,NCOL,NROW,
     1                          NLAY,IOUT)
      ELSE
        IF(IBD.EQ.1) CALL UBUDSVU(KSTP,KPER,TEXT,IBCTCB,BUFF,NODES,
     1                          IOUT,PERTIM,TOTIM)
      ENDIF
C
C16------MOVE RATES, VOLUMES & LABELS INTO ARRAYS FOR PRINTING.
  200 RIN=RATIN
      ROUT=RATOUT
      VBVLT(3,MSUMT,ICOMP)=RIN
      VBVLT(4,MSUMT,ICOMP)=ROUT
      VBVLT(1,MSUMT,ICOMP)=VBVLT(1,MSUMT,ICOMP)+RATIN*DELT
      VBVLT(2,MSUMT,ICOMP)=VBVLT(2,MSUMT,ICOMP)+RATOUT*DELT
      VBNMT(MSUMT,ICOMP)=TEXT
C
C17------INCREMENT BUDGET TERM COUNTER(MSUM).
      MSUMT=MSUMT+1
C
C18------RETURN
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE GWT2ETS1BD(KSTP,KPER,ICOMP)
C     ******************************************************************
C     CALCULATE MASS BUDGET FOR ETS TRANSPORT BOUNDARY CONDITIONS
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,   ONLY:IOUT,NCOL,NROW,NLAY,IBOUND,BUFF,NODES,IUNSTR,
     1              AMAT,IA,JA
      USE GWFBASMODULE,ONLY:MSUM,ICBCFL,IAUXSV,DELT,PERTIM,TOTIM
      USE GWFETSMODULE,ONLY:ETSF,INIETS,IETS,IETSCB,ESFACTOR,NETSOP
      USE GWTBCTMODULE, ONLY: ICBUND,CONC,MSUMT,VBVLT,VBNMT,IPCBFLAG,
     1  IBCTCB,MCOMPT
C
      CHARACTER*16 TEXT
      DOUBLE PRECISION RATIN,RATOUT,QQ,BIG,CDIFF
      DATA TEXT /'   ETS MASS FLUX'/
C     ------------------------------------------------------------------
C
C1------CLEAR RATIN AND RATOUT ACCUMULATORS, AND SET CELL-BY-CELL
C1------BUDGET FLAG.
      ZERO=0.
      RATIN=ZERO
      RATOUT=ZERO
      IBD=0
      IF(IBCTCB.LT.0 .AND. ICBCFL.NE.0) IBD=-1
      IF(IBCTCB.GT.0 .AND. ICBCFL.NE.0) IBD= 1      ! should be ICBCFL but IBD =2 is not yet in EVT
      IBDLBL=0
C
C2------CLEAR THE BUFFER.
      DO 50 N=1,NODES
      BUFF(N)=ZERO
50    CONTINUE
C
C3------IF THERE ARE NO ETSs, DO NOT ACCUMULATE
      IF(INIETS.EQ.0) GO TO 200
      IF(ICOMP.GT.MCOMPT) GO TO 200
C
C4------LOOP THROUGH EACH ETS CALCULATING MASS FLUX
      DO 100 L=1,INIETS
C
C5-----GET NODE NUMBER OF CELL CONTAINING ETS.
      N=IETS(L)
      IF(NETSOP.EQ.3.AND.IBOUND(N).EQ.0)THEN
        CALL FIRST_ACTIVE_BELOW(N)
      ENDIF      
      QQ=ZERO
      Q = 0.0
C
C6-----IF THE CELL IS INACTIVE OR PCB, IGNORE IT.
      IF(ICBUND(N).EQ.0.OR.IPCBFLAG(N).EQ.1)GO TO 99
C
C7-----COMPUTE MASS FLUX AT ETS BOUNDARY
        QF = ETSF(L)
C
C8-----ETS IS ALWAYS OUTFLOW, MASS IS Q * CONC
        QQ = QF * CONC(N,ICOMP) * ESFACTOR(ICOMP)
         Q = QQ
C
C9-----PRINT FLOW RATE IF REQUESTED.
        IF(IBD.LT.0) THEN
          IF(IBDLBL.EQ.0) WRITE(IOUT,61) TEXT,KPER,KSTP
   61     FORMAT(1X,/1X,A,'   PERIOD ',I4,'   STEP ',I3)
          IF(IUNSTR.EQ.0)THEN
            IL = (N-1) / (NCOL*NROW) + 1
            IJ = N - (IL-1)*NCOL*NROW
            IR = (IJ-1)/NCOL + 1
            IC = IJ - (IR-1)*NCOL
            WRITE(IOUT,62) L,IL,IR,IC,Q
   62       FORMAT(1X,'ETS  ',I6,'   LAYER ',I3,'   ROW ',I5,'   COL ',
     1       I5, '   FLUX ',1PG15.6)
          ELSE
            WRITE(IOUT,63) L,N,Q
   63       FORMAT(1X,'ETS  ',I6,'    NODE ',I8,'   FLUX ',1PG15.6)
          ENDIF
          IBDLBL=1
        END IF
C
C10-----ADD FLOW RATE TO BUFFER.
        BUFF(N)=BUFF(N)+Q
C
C11-----SEE IF FLUX IS POSITIVE OR NEGATIVE.
        IF(QQ.GE.ZERO) THEN
C
C12-------FLOW RATE IS POSITIVE (RECHARGE). ADD IT TO RATIN.
          RATIN=RATIN+QQ
        ELSE
C
C13-------FLOW RATE IS NEGATIVE (DISCHARGE). ADD IT TO RATOUT.
          RATOUT=RATOUT-QQ
        END IF
C
C14-----IF SAVING CELL-BY-CELL FLOWS IN A LIST, WRITE FLOW.  ALSO
C14-----COPY FLOW TO ETS LIST.

   99 CONTINUE
C
100   CONTINUE
C
C15------IF CELL-BY-CELL FLOWS WILL BE SAVED AS A 3-D ARRAY,
C15------CALL UBUDSV TO SAVE THEM.
      IF(IUNSTR.EQ.0)THEN
        IF(IBD.EQ.1) CALL UBUDSV(KSTP,KPER,TEXT,IBCTCB,BUFF,NCOL,NROW,
     1                          NLAY,IOUT)
      ELSE
        IF(IBD.EQ.1) CALL UBUDSVU(KSTP,KPER,TEXT,IBCTCB,BUFF,NODES,
     1                          IOUT,PERTIM,TOTIM)
      ENDIF
C
C16------MOVE RATES, VOLUMES & LABELS INTO ARRAYS FOR PRINTING.
  200 RIN=RATIN
      ROUT=RATOUT
      VBVLT(3,MSUMT,ICOMP)=RIN
      VBVLT(4,MSUMT,ICOMP)=ROUT
      VBVLT(1,MSUMT,ICOMP)=VBVLT(1,MSUMT,ICOMP)+RATIN*DELT
      VBVLT(2,MSUMT,ICOMP)=VBVLT(2,MSUMT,ICOMP)+RATOUT*DELT
      VBNMT(MSUMT,ICOMP)=TEXT
C
C17------INCREMENT BUDGET TERM COUNTER(MSUM).
      MSUMT=MSUMT+1
C
C18------RETURN
      RETURN
      END
C----------------------------------------------------------------------------
      SUBROUTINE GWT2RCH1FM(ICOMP)
C     ******************************************************************
C     FORMULATE RCH TRANSPORT BOUNDARY CONDITION
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,       ONLY:IBOUND,RHS,IA,JA,AMAT
      USE GWFRCHMODULE, ONLY:RCHF,INIRCH,IRCH,IRCHCB,NRCHOP,
     &  RCHCONC,IRCHCONC
      USE GWTBCTMODULE, ONLY: ICBUND
      DOUBLE PRECISION Co,Q
C     ------------------------------------------------------------------
C
C1------IF NUMBER OF RCH CELLS <= 0 THEN RETURN.
      IF(INIRCH.LE.0) RETURN
!      IF (NOT(ALLOCATED(IRCHCONC))) RETURN
      IF (.NOT.(ALLOCATED(IRCHCONC))) RETURN      !kkz - update per JCH for GFortran; alternative is to use ALLOCATED(IRCHCONC) .EQV. .FALSE.
      BIG = 1.0E20
C-------FIND COMPONENT NUMBER IF RECHARGE CONCENTRATIONS ARE READ
      IF(IRCHCONC(ICOMP).EQ.1)THEN
          ICONCRCH = 0
          DO II=1,ICOMP
            ICONCRCH = ICONCRCH + IRCHCONC(II)
          ENDDO
        ENDIF
C
C2------PROCESS EACH RCH NODE IN THE LIST.
      DO 100 L=1,INIRCH
      N=IRCH(L)
C-------FIND TOP-MOST ACTIVE NODE IF NOT N
        IF(NRCHOP.EQ.3.AND.ICBUND(N).EQ.0)THEN
          CALL FIRST_ACTIVE_BELOW(N)
        ENDIF
C
C2A-----IF THE CELL IS INACTIVE THEN BYPASS PROCESSING.
      IF(ICBUND(N).EQ.0) GO TO 100
        Q = RCHF(L)
        Co= 0.0
        IF(IRCHCONC(ICOMP).EQ.1)THEN
          Co = RCHCONC(L,ICONCRCH)
        ENDIF
C
C2B-----IF THE CELL IS OUTFLOW, PUT Q ON LHS DIAGONAL
        IF(Q.LT.0.0)THEN
          AMAT(IA(N)) = AMAT(IA(N)) + Q
        ELSE
C2C-------IF THE CELL IS INFLOW, PUT Q*Co ON RHS
          RHS(N) = RHS(N) -Co * Q
        ENDIF
C
  100 CONTINUE
C
C3------RETURN
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE GWT2RCH1BD(KSTP,KPER,ICOMP)
C     ******************************************************************
C     CALCULATE MASS BUDGET FOR RCH TRANSPORT BOUNDARY CONDITIONS
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,   ONLY:IOUT,NCOL,NROW,NLAY,IBOUND,BUFF,NODES,IUNSTR,
     1              AMAT,IA,JA
      USE GWFBASMODULE,ONLY:MSUM,ICBCFL,IAUXSV,DELT,PERTIM,TOTIM
      USE GWFRCHMODULE,ONLY:RCHF,INIRCH,IRCH,IRCHCB,NRCHOP,
     &  RCHCONC,IRCHCONC
      USE GWTBCTMODULE, ONLY: ICBUND,CONC,MSUMT,VBVLT,VBNMT,IPCBFLAG,
     1  IBCTCB,MCOMPT
C
      CHARACTER*16 TEXT
      DOUBLE PRECISION RATIN,RATOUT,QQ,BIG,CDIFF
      DATA TEXT /'   RCH MASS FLUX'/
C     ------------------------------------------------------------------
C
C1------CLEAR RATIN AND RATOUT ACCUMULATORS, AND SET CELL-BY-CELL
C1------BUDGET FLAG.
      ZERO=0.
      RATIN=ZERO
      RATOUT=ZERO
      IBD=0
      IF(IBCTCB.LT.0 .AND. ICBCFL.NE.0) IBD=-1
      IF(IBCTCB.GT.0 .AND. ICBCFL.NE.0) IBD= 1   ! should be ICBCFL but IBD =2 is not yet in RCH
      IBDLBL=0
C
C2------CLEAR THE BUFFER.
      DO 50 N=1,NODES
      BUFF(N)=ZERO
50    CONTINUE
C
C3------IF THERE ARE NO RCHs, DO NOT ACCUMULATE
      IF(INIRCH.EQ.0) GO TO 200
      IF(ICOMP.GT.MCOMPT) GO TO 200
!      IF(NOT(ALLOCATED(IRCHCONC))) GO TO 200
      IF(.NOT.(ALLOCATED(IRCHCONC))) GO TO 200          !kkz - update per JCH for GFortran; alternative is to use ALLOCATED(IRCHCONC) .EQV. .FALSE.
C-------FIND COMPONENT NUMBER IF RECHARGE CONCENTRATIONS ARE READ
      IF(IRCHCONC(ICOMP).EQ.1)THEN
          ICONCRCH = 0
          DO II=1,ICOMP
            ICONCRCH = ICONCRCH + IRCHCONC(II)
          ENDDO
        ENDIF
C
C4------LOOP THROUGH EACH RCH CALCULATING MASS FLUX
      DO 100 L=1,INIRCH
C
C5------GET NODE NUMBER OF CELL CONTAINING RCH.
      N=IRCH(L)
C-------FIND TOP-MOST ACTIVE NODE IF NOT N
        IF(NRCHOP.EQ.3.AND.ICBUND(N).EQ.0)THEN
          CALL FIRST_ACTIVE_BELOW(N)
        ENDIF
C
      QQ=ZERO
      Q = 0.0
C
C6------IF THE CELL IS INACTIVE OR PCB, IGNORE IT.
      IF(ICBUND(N).EQ.0.OR.IPCBFLAG(N).EQ.1)GO TO 99
C
C7------COMPUTE MASS FLUX AT RCH BOUNDARY
        QF = RCHF(L)
        Co = 0.0
        IF(IRCHCONC(ICOMP).EQ.1)THEN
          Co = RCHCONC(L,ICONCRCH)
        ENDIF
C
C8------IF THE CELL IS OUTFLOW, MASS IS Q * CONC
        IF(QF.LT.0.0)THEN
          QQ = QF * CONC(N,ICOMP)
        ELSE
C9--------IF THE CELL IS INFLOW, MASS IS Q*Co
          QQ = QF * Co
        ENDIF
        Q = QQ
C
C10-----PRINT FLOW RATE IF REQUESTED.
        IF(IBD.LT.0) THEN
          IF(IBDLBL.EQ.0) WRITE(IOUT,61) TEXT,KPER,KSTP
   61     FORMAT(1X,/1X,A,'   PERIOD ',I4,'   STEP ',I3)
          IF(IUNSTR.EQ.0)THEN
            IL = (N-1) / (NCOL*NROW) + 1
            IJ = N - (IL-1)*NCOL*NROW
            IR = (IJ-1)/NCOL + 1
            IC = IJ - (IR-1)*NCOL
            WRITE(IOUT,62) L,IL,IR,IC,Q
   62       FORMAT(1X,'RCH  ',I6,'   LAYER ',I3,'   ROW ',I5,'   COL ',
     1       I5, '   FLUX ',1PG15.6)
          ELSE
            WRITE(IOUT,63) L,N,Q
   63       FORMAT(1X,'RCH  ',I6,'    NODE ',I8,'   FLUX ',1PG15.6)
          ENDIF
          IBDLBL=1
        END IF
C
C11-----ADD FLOW RATE TO BUFFER.
        BUFF(N)=BUFF(N)+Q
C
C12-----SEE IF FLUX IS POSITIVE OR NEGATIVE.
        IF(QQ.GE.ZERO) THEN
C
C13-------FLOW RATE IS POSITIVE (RECHARGE). ADD IT TO RATIN.
          RATIN=RATIN+QQ
        ELSE
C
C14-------FLOW RATE IS NEGATIVE (DISCHARGE). ADD IT TO RATOUT.
          RATOUT=RATOUT-QQ
        END IF
C
C15-----IF SAVING CELL-BY-CELL FLOWS IN A LIST, WRITE FLOW.  ALSO
C15-----COPY FLOW TO RCH LIST.

   99 CONTINUE
C
100   CONTINUE
C
C16------IF CELL-BY-CELL FLOWS WILL BE SAVED AS A 3-D ARRAY,
C16------CALL UBUDSV TO SAVE THEM.
      IF(IUNSTR.EQ.0)THEN
        IF(IBD.EQ.1) CALL UBUDSV(KSTP,KPER,TEXT,IBCTCB,BUFF,NCOL,NROW,
     1                          NLAY,IOUT)
      ELSE
        IF(IBD.EQ.1) CALL UBUDSVU(KSTP,KPER,TEXT,IBCTCB,BUFF,NODES,
     1                          IOUT,PERTIM,TOTIM)
      ENDIF
C
C17------MOVE RATES, VOLUMES & LABELS INTO ARRAYS FOR PRINTING.
  200 RIN=RATIN
      ROUT=RATOUT
      VBVLT(3,MSUMT,ICOMP)=RIN
      VBVLT(4,MSUMT,ICOMP)=ROUT
      VBVLT(1,MSUMT,ICOMP)=VBVLT(1,MSUMT,ICOMP)+RATIN*DELT
      VBVLT(2,MSUMT,ICOMP)=VBVLT(2,MSUMT,ICOMP)+RATOUT*DELT
      VBNMT(MSUMT,ICOMP)=TEXT
C
C18------INCREMENT BUDGET TERM COUNTER(MSUM).
      MSUMT=MSUMT+1
C
C19------RETURN
      RETURN
      END
C----------------------------------------------------------------------------
      SUBROUTINE GWT2RIV1FM(ICOMP)
C     ******************************************************************
C     FORMULATE RIV TRANSPORT BOUNDARY CONDITION
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,       ONLY:IBOUND,RHS,IA,JA,AMAT
      USE GWFRIVMODULE, ONLY:NRIVER,RIVR,NRIVVL,RIVAUX
      USE GWTBCTMODULE, ONLY: ICBUND
      DOUBLE PRECISION Co,Q
C     ------------------------------------------------------------------
C
C1------IF NUMBER OF RIV NODES <= 0 THEN RETURN.
      IF(NRIVER.LE.0) RETURN
      BIG = 1.0E20
      NAUX=NRIVVL-7
      CALL CONCIAUX(ICOMP,NAUX,RIVAUX,IAUX)
C
C2------PROCESS EACH RIV NODE IN THE LIST.
      DO 100 L=1,NRIVER
      N=RIVR(1,L)
C
C2A-----IF THE CELL IS INACTIVE THEN BYPASS PROCESSING.
      IF(ICBUND(N).EQ.0) GO TO 100
        Q = RIVR(NRIVVL,L)
C
C2B-----IF THE CELL IS OUTFLOW, PUT Q ON LHS DIAGONAL
        IF(Q.LE.0.0)THEN
          AMAT(IA(N)) = AMAT(IA(N)) + Q
        ELSE
C2C-------IF THE CELL IS INFLOW, PUT Q*Co ON RHS
          Co = 0.0
          IF(IAUX.NE.0) Co= RIVR(6+IAUX,L)
          RHS(N) = RHS(N) -Co * Q
        ENDIF
C
  100 CONTINUE
C
C3------RETURN
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE GWT2RIV1BD(KSTP,KPER,ICOMP)
C     ******************************************************************
C     CALCULATE MASS BUDGET FOR RIV TRANSPORT BOUNDARY CONDITIONS
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,   ONLY:IOUT,NCOL,NROW,NLAY,IBOUND,BUFF,NODES,IUNSTR,
     1              AMAT,IA,JA,INCLN
      USE GWFBASMODULE,ONLY:MSUM,ICBCFL,IAUXSV,DELT,PERTIM,TOTIM
      USE GWFRIVMODULE,ONLY:NRIVER,IRIVCB,RIVR,NRIVVL,RIVAUX
      USE GWTBCTMODULE, ONLY: ICBUND,CONC,MSUMT,VBVLT,VBNMT,IPCBFLAG,
     1  IBCTCB,MCOMPT
      USE CLN1MODULE, ONLY: ICLNMB,NCLNNDS
C
      CHARACTER*16 TEXT
      DOUBLE PRECISION RATIN,RATOUT,QQ,BIG,CDIFF
      DATA TEXT /'   RIV MASS FLUX'/
C     ------------------------------------------------------------------
C
C1------CLEAR RATIN AND RATOUT ACCUMULATORS, AND SET CELL-BY-CELL
C1------BUDGET FLAG.
      ZERO=0.
      RATIN=ZERO
      RATOUT=ZERO
      IBD=0
      IF(IBCTCB.LT.0 .AND. ICBCFL.NE.0) IBD=-1
      IF(IBCTCB.GT.0) IBD=ICBCFL
      IBDLBL=0
      NAUX=NRIVVL-7
      CALL CONCIAUX(ICOMP,NAUX,RIVAUX,IAUX)
C
C2-----IF CELL-BY-CELL FLOWS WILL BE SAVED AS A LIST, WRITE HEADER.
      IICLNCB=0
      IF(IBD.EQ.2) THEN
        NAUX=NRIVVL-7
        IF(IAUXSV.EQ.0) NAUX=0
        NNCLNNDS=0
        IF(INCLN.GT.0) THEN
          IICLNCB=ICLNMB
          NNCLNNDS=NCLNNDS
        ENDIF
         CALL UBDSVHDR(IUNSTR,KSTP,KPER,IOUT,IBCTCB,IICLNCB,NODES,
     1    NNCLNNDS,NCOL,NROW,NLAY,NRIVER,NRIVVL,NAUX,IBOUND,
     2    TEXT,RIVAUX,DELT,PERTIM,TOTIM,RIVR)
      END IF
C
C3------CLEAR THE BUFFER.
      DO 50 N=1,NODES
      BUFF(N)=ZERO
50    CONTINUE
C
C4------IF THERE ARE NO RIVs, DO NOT ACCUMULATE
      IF(NRIVER.EQ.0) GO TO 200
      IF(ICOMP.GT.MCOMPT) GO TO 200
C
C5------LOOP THROUGH EACH RIV CALCULATING MASS FLUX
      DO 100 L=1,NRIVER
C
C6-----GET NODE NUMBER OF CELL CONTAINING RIV.
      N=RIVR(1,L)
      QQ=ZERO
      Q = 0.0
C
C7-----IF THE CELL IS INACTIVE OR PCB, IGNORE IT.
      IF(ICBUND(N).EQ.0.OR.IPCBFLAG(N).EQ.1)GO TO 99
C
C8-----COMPUTE MASS FLUX AT RIVER BOUNDARY
        QF = RIVR(NRIVVL,L)
C
C9-----IF THE CELL IS OUTFLOW, MASS IS Q * CONC
        IF(QF.LE.0.0)THEN
          QQ = QF * CONC(N,ICOMP)
        ELSE
C10-------IF THE CELL IS INFLOW, MASS IS Q*Co
          Co = 0.0
          IF(IAUX.NE.0) Co = RIVR(6+IAUX,L)
          QQ = QF * Co
        ENDIF
        Q = QQ
C
C11-----PRINT FLOW RATE IF REQUESTED.
        IF(IBD.LT.0) THEN
          IF(IBDLBL.EQ.0) WRITE(IOUT,61) TEXT,KPER,KSTP
   61     FORMAT(1X,/1X,A,'   PERIOD ',I4,'   STEP ',I3)
          IF(IUNSTR.EQ.0)THEN
            IL = (N-1) / (NCOL*NROW) + 1
            IJ = N - (IL-1)*NCOL*NROW
            IR = (IJ-1)/NCOL + 1
            IC = IJ - (IR-1)*NCOL
            WRITE(IOUT,62) L,IL,IR,IC,Q
   62       FORMAT(1X,'RIV  ',I6,'   LAYER ',I3,'   ROW ',I5,'   COL ',
     1       I5, '   FLUX ',1PG15.6)
          ELSE
            WRITE(IOUT,63) L,N,Q
   63       FORMAT(1X,'RIV  ',I6,'    NODE ',I8,'   FLUX ',1PG15.6)
          ENDIF
          IBDLBL=1
        END IF
C
C12-----ADD FLOW RATE TO BUFFER.
        BUFF(N)=BUFF(N)+Q
C
C13-----SEE IF FLUX IS POSITIVE OR NEGATIVE.
        IF(QQ.GE.ZERO) THEN
C
C14-------FLOW RATE IS POSITIVE (RECHARGE). ADD IT TO RATIN.
          RATIN=RATIN+QQ
        ELSE
C
C15-------FLOW RATE IS NEGATIVE (DISCHARGE). ADD IT TO RATOUT.
          RATOUT=RATOUT-QQ
        END IF
C
C16-----IF SAVING CELL-BY-CELL FLOWS IN A LIST, WRITE FLOW.  ALSO
C16-----COPY FLOW TO RIVER LIST.

   99 CONTINUE
      IF(IBD.EQ.2) THEN
        CALL UBDSVREC(IUNSTR,N,NODES,NNCLNNDS,IBCTCB,IICLNCB,NRIVVL,
     1    7,NAUX,Q,RIVR(:,L),IBOUND,NCOL,NROW,NLAY)
      ENDIF
100   CONTINUE
C
C17------IF CELL-BY-CELL FLOWS WILL BE SAVED AS A 3-D ARRAY,
C17------CALL UBUDSV TO SAVE THEM.
      IF(IUNSTR.EQ.0)THEN
        IF(IBD.EQ.1)CALL UBUDSV(KSTP,KPER,TEXT,IBCTCB,BUFF(1),NCOL,NROW,
     1                          NLAY,IOUT)
      ELSE
        IF(IBD.EQ.1) CALL UBUDSVU(KSTP,KPER,TEXT,IBCTCB,BUFF(1),NODES,
     1                          IOUT,PERTIM,TOTIM)
      ENDIF
      IF(IBD.EQ.1.AND.INCLN.GT.0)THEN
        IF(IICLNCB.GT.0) CALL UBUDSVU(KSTP,KPER,TEXT,IICLNCB,
     1    BUFF(NODES+1),NCLNNDS,IOUT,PERTIM,TOTIM)
      ENDIF
C
C18------MOVE RATES, VOLUMES & LABELS INTO ARRAYS FOR PRINTING.
  200 RIN=RATIN
      ROUT=RATOUT
      VBVLT(3,MSUMT,ICOMP)=RIN
      VBVLT(4,MSUMT,ICOMP)=ROUT
      VBVLT(1,MSUMT,ICOMP)=VBVLT(1,MSUMT,ICOMP)+RATIN*DELT
      VBVLT(2,MSUMT,ICOMP)=VBVLT(2,MSUMT,ICOMP)+RATOUT*DELT
      VBNMT(MSUMT,ICOMP)=TEXT
C
C19------INCREMENT BUDGET TERM COUNTER(MSUM).
      MSUMT=MSUMT+1
C
C20------RETURN
      RETURN
      END
C----------------------------------------------------------------------------
      SUBROUTINE CALCLAKECONC(ICOMP,IUNITSFR)
C     ******************************************************************
C     CALCULATE LAKE CONCENTRATION
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GWFLAKMODULE,ONLY: ILAKE,LKNODE,LKNODE,FLOB,RMASLAK,VOUTLAK,
     1                       CLAKE,CLKOLD,NLAKES,SURFA,PRCPLK,RNF,
     1                       WTHDRW,EVAPLK,CPPT,CRNF,CAUG,VOL,VOLOLD,
     1                       NTRB,ITRB,NDV,IDIV,IETLAK
      USE GWTBCTMODULE,ONLY: CONC
      USE GWFBASMODULE,ONLY: DELT
      USE GWFSFRMODULE, ONLY: STRIN, STROUT, FXLKOT, DLKSTAGE, SEG
      IMPLICIT NONE
      INTEGER ICOMP,IUNITSFR
      INTEGER L,N,LK,NGW,LM,ITRIB,INODE,IDV
      REAL Q,C,CO,VO
C     ------------------------------------------------------------------
C
C--ZERO OUT TERMS
      RMASLAK=0.
      VOUTLAK=0.
C
C--ALL INFLOW TERMS-------------------------------------------------
C
C--GW TO LAK FLOW
      DO L=1,LKNODE
        NGW=ILAKE(1,L) !GW NODE
        N=ILAKE(4,L) !LAKE NUMBER
        Q=0.
        Q=FLOB(L)
        C=CONC(NGW,ICOMP)
C.......CONSIDER ONLY FLOW INTO LAKE
        IF(Q.LE.0.) THEN
          RMASLAK(N)=RMASLAK(N)+C*ABS(Q)
        ENDIF
      ENDDO
C
C--INFLOW FROM BCs
      DO LM=1,NLAKES
        N=LM
C.......PRECIPITATION
        Q=PRCPLK(LM)*SURFA(LM)
        C=CPPT(LM,ICOMP)
        RMASLAK(N)=RMASLAK(N)+C*ABS(Q)
C.......RUNOFF
        Q=RNF(LM)
        C=CRNF(LM,ICOMP)
        RMASLAK(N)=RMASLAK(N)+C*ABS(Q)
C.......INJECTION
        Q=WTHDRW(LM)
        C=CAUG(LM,ICOMP)
        IF(WTHDRW(LM).LT.0.0) RMASLAK(N)=RMASLAK(N)+C*ABS(Q)
      ENDDO
C
C--INFLOW FROM SFR - ASSUME SFR INFLOW WITH ZERO CONCENTRATION
      IF (IUNITSFR.GT.0) THEN
        DO 200 LK=1,NLAKES
        DO 200 ITRIB=1,NTRB
            INODE=ITRB(LK,ITRIB)
            IF (INODE.LE.0) GO TO 200
            Q=STRIN(INODE)
            C=0.
            RMASLAK(N)=RMASLAK(N)+C*ABS(Q)
200     CONTINUE
      END IF
C
C--UZF TO LAK FLOW
C      IF(iUnitTRNOP(7).GT.0) THEN
C        DO ICON=1,NLAK2UZF
C            N=INOD1LKUZ(ICON)
C            NUZ=INOD2LKUZ(ICON)
C            CALL NODE2KIJ(NUZ,NLAY,NROW,NCOL,K,I,J)
C            Q=QLAK2UZF(ICON)
C            Q=ABS(Q)
C            IF(IUZCODELK(ICON).EQ.1) THEN
C              CONC=CNEW(J,I,K,ICOMP)
C              RMASLAK(N)=RMASLAK(N)+CONC*ABS(Q)
C            ELSEIF(IUZCODELK(ICON).EQ.2) THEN
C              CONC=CUZINF(J,I,ICOMP)
C              RMASLAK(N)=RMASLAK(N)+CONC*ABS(Q)
C            ELSEIF(IUZCODELK(ICON).EQ.3) THEN
C              CONC=CUZINF(J,I,ICOMP)
C              RMASLAK(N)=RMASLAK(N)+CONC*ABS(Q)
C            ELSE
C              WRITE(IOUT,*) 'CHECK FTL FILE - IUZCODELK(ICON) INVALID'
C              WRITE(*,*) 'CHECK FTL FILE - IUZCODELK(ICON) INVALID'
C              STOP
C            ENDIF
C        ENDDO
C      ENDIF
C
C--MULTIPLY MASS TERM WITH TRANSPORT TIME-STEP
      RMASLAK=RMASLAK*DELT
C
C--ADD OLD MASS IN THE LAKE
C      IF(iSSTrans.EQ.0) THEN
      DO N=1,NLAKES
        CO=CLKOLD(N,ICOMP)
        VO=VOLOLD(N)
        RMASLAK(N)=RMASLAK(N)+VO*CO
      ENDDO
C      ENDIF
C-------------------------------------------------------------------
C
C--CALCULATE LAKE CONCENTRATION-------------------------------------
C
C--CALCULATE TOTAL VOLUME FOR CONC CALCULATION
C--IT IS ASSUMED THAT ET IS TAKEN OUT FIRST SO THAT LAKE AND OUT CONC IS SAME
C
C--LAK TO GW FLOW
      DO L=1,LKNODE
        NGW=ILAKE(1,L) !GW NODE
        N=ILAKE(4,L) !LAKE NUMBER
        Q=0.
        Q=FLOB(L)
C.......CONSIDER ONLY FLOW OUT OF LAKE
        IF(Q.GT.0.) THEN
          VOUTLAK(N)=VOUTLAK(N)+Q
        ENDIF
      ENDDO
C
C--OUTFLOW TO BOUNDARY CONDITIONS (ONLY WITHDRAWAL)
      DO N=1,NLAKES
        LM=N
        Q=0.
        Q=WTHDRW(LM)
        IF(WTHDRW(LM).GE.0.0) VOUTLAK(N)=VOUTLAK(N)+Q
        IF(IETLAK.NE.0) THEN
          Q=EVAPLK(LM)*SURFA(LM)
          VOUTLAK(N)=VOUTLAK(N)+Q
        ENDIF
      ENDDO
C
C--FLOW TO STREAMS
C
C--INFLOW FROM SFR - ASSUME SFR INFLOW WITH ZERO CONCENTRATION
      IF (IUNITSFR.GT.0) THEN
        DO LK=1,NLAKES
        DO IDV=1,NDV
            INODE=IDIV(LK,IDV)
            IF (INODE.GT.0) THEN
              Q=STROUT(INODE)
              VOUTLAK(LK)=VOUTLAK(LK)+ABS(Q)
            ENDIF
        ENDDO
        ENDDO
      END IF
C
C-----LAKE VOLUME CHANGE
C      IF(iSSTrans.EQ.0) THEN
      DO N=1,NLAKES
        VO=VOLOLD(N)
        VOUTLAK(N)=VOUTLAK(N)*DELT+VOL(N)
      ENDDO
C      ENDIF
C
C--CALCULATE LAKE AND OUT CONCENTRATION
      DO N=1,NLAKES
        IF(VOUTLAK(N).LE.1E-6) THEN
          CLAKE(N,ICOMP)=0.
        ELSE
          CLAKE(N,ICOMP)=RMASLAK(N)/VOUTLAK(N)
        ENDIF
      ENDDO
C
C3------RETURN
      RETURN
      END
C----------------------------------------------------------------------------
      SUBROUTINE GWT2LAK1FM(ICOMP)
C     ******************************************************************
C     FORMULATE LAK TRANSPORT BOUNDARY CONDITION
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,       ONLY:IBOUND,RHS,IA,JA,AMAT
      USE GWFLAKMODULE, ONLY:LKNODE,FLOB,ILAKE,CLAKE
      USE GWTBCTMODULE, ONLY: ICBUND
      INTEGER ICOMP
      DOUBLE PRECISION Co,Q
C     ------------------------------------------------------------------
C
C1------IF NUMBER OF LKNODE <= 0 THEN RETURN.
      IF(LKNODE.LE.0) RETURN
      BIG = 1.0E20
C
C2------PROCESS EACH LAK-NODE IN THE LIST.
      DO 100 L=1,LKNODE
      N=ILAKE(1,L)
      LK=ILAKE(4,L)
C
C2A-----IF THE CELL IS INACTIVE THEN BYPASS PROCESSING.
C      IF(ICBUND(N).EQ.0) GO TO 100
      Q = FLOB(L)
C
C2B-----IF THE CELL IS OUTFLOW, PUT Q ON LHS DIAGONAL
      IF(Q.LE.0.0)THEN
        AMAT(IA(N)) = AMAT(IA(N)) + Q
      ELSE
C2C-------IF THE CELL IS INFLOW, PUT Q*Co ON RHS
        Co= CLAKE(LK,ICOMP)
        RHS(N) = RHS(N) -Co * Q
      ENDIF
C
  100 CONTINUE
C
C3------RETURN
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE GWT2LAK1BD(KSTP,KPER,ICOMP,IUNITSFR)
C     ******************************************************************
C     CALCULATE MASS BUDGET FOR LAK TRANSPORT BOUNDARY CONDITIONS
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,   ONLY:IOUT,NCOL,NROW,NLAY,IBOUND,BUFF,NODES,IUNSTR,
     1              AMAT,IA,JA,INCLN
      USE GWFBASMODULE,ONLY:MSUM,ICBCFL,IAUXSV,DELT,PERTIM,TOTIM,ISPCFL
      USE GWFLAKMODULE,ONLY:LKNODE,ILKCB,FLOB,ILAKE,CLAKE,RMASLAK,
     1                      VOUTLAK,ILKTRNSPT,NLAKES,SURFA,PRCPLK,RNF,
     1                      WTHDRW,EVAPLK,CPPT,CRNF,CAUG,VOL,VOLOLD,
     1                      NTRB,ITRB,NDV,IDIV,IETLAK,CLKOLD,CGW2LAK,
     1       CUZF2LAK,CGWFROMLAK,CSFR2LAK,CSFRFROMLAK,CPRECLK,CRUNOFLK,
     1       CWDRLLK,CETLK,CSTORINLK,CSTOROTLK
      USE GWFSFRMODULE, ONLY: STRIN, STROUT, FXLKOT, DLKSTAGE, SEG
      USE GWTBCTMODULE, ONLY: ICBUND,CONC,MSUMT,VBVLT,VBNMT,IPCBFLAG,
     1  IBCTCB,MCOMPT
      USE CLN1MODULE, ONLY: ICLNMB,NCLNNDS
C
      IMPLICIT NONE
      CHARACTER*16 TEXT
      DOUBLE PRECISION RATIN,RATOUT,QQ,BIG,CDIFF
      INTEGER KPER,KSTP,ICOMP,IUNITSFR
      INTEGER NLAKVL,NAUX,N,L,LK,IBD,IBDLBL,NNCLNNDS,IL,IJ,IK,IR,IC,
     1        IICLNCB,LM,ITRIB,INODE,IDV
      REAL C,Q,GW2LAK,GWFROMLAK,UZF2LAK,SFR2LAK,SFRFROMLAK,PRECLK,
     1     RUNOFLK,WDRLLK,WDRLINJLK,ETLK,STORINLK,STOROTLK,Q1,Q2,ZERO,
     1     QF,CO,RATE,RIN,ROUT,VO,CN,VN,TOTMASOLD,TOTMASNEW,STORDIFF,
     1     TOTINLK,TOTOUTLK,CTOTINLK,CTOTOUTLK,DIFF,PERC,CPERC,PRTOUT
      DATA TEXT /'   LAK MASS FLUX'/
C     ------------------------------------------------------------------
C
C1------CLEAR RATIN AND RATOUT ACCUMULATORS, AND SET CELL-BY-CELL
C1------BUDGET FLAG.
      ZERO=0.
      RATIN=ZERO
      RATOUT=ZERO
      IBD=0
      NAUX=1
      NLAKVL=1
C
C
C--ZERO OUT TERMS
      C=0.
      Q=0.
      RMASLAK=0.
      VOUTLAK=0.
C
      GW2LAK=0.
      GWFROMLAK=0.
      UZF2LAK=0.
      SFR2LAK=0.
      SFRFROMLAK=0.
      PRECLK=0.
      RUNOFLK=0.
      WDRLLK=0.
      WDRLINJLK=0.
      ETLK=0.
      STORINLK=0.
      STOROTLK=0.
      Q1=0.
      Q2=0.
C
      IF(IBCTCB.LT.0 .AND. ICBCFL.NE.0) IBD=-1
      IF(IBCTCB.GT.0) IBD=ICBCFL
      IBDLBL=0
C
C2-----IF CELL-BY-CELL FLOWS WILL BE SAVED AS A LIST, WRITE HEADER.
      IICLNCB=0
      IF(IBD.EQ.2) THEN
         NNCLNNDS=0
         IF(INCLN.GT.0) THEN
           IICLNCB=ICLNMB
           NNCLNNDS=NCLNNDS
         ENDIF
CVSB         NAUX=1
CVSB         NLAKVL=1
CVSB         CALL UBDSVHDR(IUNSTR,KSTP,KPER,IOUT,IBCTCB,IICLNCB,NODES,
CVSB     1    NNCLNNDS,NCOL,NROW,NLAY,NBOUND,NLAKVL,NAUX,IBOUND,
CVSB     2    TEXT,TEXT,DELT,PERTIM,TOTIM,FLOB)
      END IF
C
C3------CLEAR THE BUFFER.
      DO N=1,NODES
      BUFF(N)=ZERO
      ENDDO
C
C4------IF THERE ARE NO LAKs, DO NOT ACCUMULATE
      IF(LKNODE.EQ.0) GO TO 200
      IF(ICOMP.GT.MCOMPT) GO TO 200
C
C5------LOOP THROUGH EACH LAK CALCULATING MASS FLUX
      DO 100 L=1,LKNODE
C
C6-----GET NODE NUMBER OF CELL CONNECTED TO LAK.
      N=ILAKE(1,L)
      LK=ILAKE(4,L)
      QQ=ZERO
C
C7-----IF THE CELL IS INACTIVE OR PCB, IGNORE IT.
      IF(ICBUND(N).EQ.0.OR.IPCBFLAG(N).EQ.1)GO TO 99
C
C8-----COMPUTE MASS FLUX AT LAK BOUNDARY
        QF = FLOB(L)
C
C9-----IF THE CELL IS OUTFLOW, MASS IS Q * CONC
        IF(QF.LE.0.0)THEN
          QQ = QF * CONC(N,ICOMP)
          GW2LAK=GW2LAK+ABS(QQ)
        ELSE
C10-------IF THE CELL IS INFLOW, MASS IS Q*Co
          Co = CLAKE(LK,ICOMP)
          QQ = QF * Co
          GWFROMLAK=GWFROMLAK+ABS(QQ)
        ENDIF
        Q = QQ
C
C11-----PRINT FLOW RATE IF REQUESTED.
        IF(IBD.LT.0) THEN
          IF(IBDLBL.EQ.0) WRITE(IOUT,11) TEXT,KPER,KSTP
   11     FORMAT(1X,/1X,A,'   PERIOD ',I4,'   STEP ',I3)
          IF(IUNSTR.EQ.0)THEN
            IL = (N-1) / (NCOL*NROW) + 1
            IJ = N - (IL-1)*NCOL*NROW
            IR = (IJ-1)/NCOL + 1
            IC = IJ - (IR-1)*NCOL
            WRITE(IOUT,12) L,IL,IR,IC,Q
   12       FORMAT(1X,'LAK  ',I6,'   LAYER ',I3,'   ROW ',I5,'   COL ',
     1       I5, '   FLUX ',1PG15.6)
          ELSE
            WRITE(IOUT,13) L,N,Q
   13       FORMAT(1X,'LAK  ',I6,'    NODE ',I8,'   FLUX ',1PG15.6)
          ENDIF
          IBDLBL=1
        END IF
C
C12-----ADD FLOW RATE TO BUFFER.
        BUFF(N)=BUFF(N)+Q
C
C13-----SEE IF FLUX IS POSITIVE OR NEGATIVE.
        IF(QQ.GE.ZERO) THEN
C
C14-------FLOW RATE IS POSITIVE (RECHARGE). ADD IT TO RATIN.
          RATIN=RATIN+QQ
        ELSE
C
C15-------FLOW RATE IS NEGATIVE (DISCHARGE). ADD IT TO RATOUT.
          RATOUT=RATOUT-QQ
        END IF
C
C16-----IF SAVING CELL-BY-CELL FLOWS IN A LIST, WRITE FLOW.  ALSO
C16-----COPY FLOW TO GHB LIST.

   99 CONTINUE
      IF(IBD.EQ.2)THEN
        CALL UBDSVREC(IUNSTR,N,NODES,NNCLNNDS,IBCTCB,IICLNCB,NLAKVL,
     1    6,NAUX,RATE,FLOB,IBOUND,NCOL,NROW,NLAY)
      ENDIF

100   CONTINUE
C
C17------IF CELL-BY-CELL FLOWS WILL BE SAVED AS A 3-D ARRAY,
C17------CALL UBUDSV TO SAVE THEM.
      IF(IUNSTR.EQ.0)THEN
        IF(IBD.EQ.1)CALL UBUDSV(KSTP,KPER,TEXT,IBCTCB,BUFF(1),NCOL,NROW,
     1                          NLAY,IOUT)
      ELSE
        IF(IBD.EQ.1) CALL UBUDSVU(KSTP,KPER,TEXT,IBCTCB,BUFF(1),NODES,
     1                          IOUT,PERTIM,TOTIM)
      ENDIF
      IF(IBD.EQ.1.AND.INCLN.GT.0)THEN
        IF(IICLNCB.GT.0) CALL UBUDSVU(KSTP,KPER,TEXT,IICLNCB,
     1    BUFF(NODES+1),NCLNNDS,IOUT,PERTIM,TOTIM)
      ENDIF
C
  200 CONTINUE
C
C18------MOVE RATES, VOLUMES & LABELS INTO ARRAYS FOR PRINTING.
      RIN=RATIN
      ROUT=RATOUT
      VBVLT(3,MSUMT,ICOMP)=RIN
      VBVLT(4,MSUMT,ICOMP)=ROUT
      VBVLT(1,MSUMT,ICOMP)=VBVLT(1,MSUMT,ICOMP)+RATIN*DELT
      VBVLT(2,MSUMT,ICOMP)=VBVLT(2,MSUMT,ICOMP)+RATOUT*DELT
      VBNMT(MSUMT,ICOMP)=TEXT
C
C19------INCREMENT BUDGET TERM COUNTER(MSUM).
      MSUMT=MSUMT+1
C
C--CALCULATE LAKE BUDGET AND WRITE LAKE MASS BALANCE
      IF(ILKTRNSPT.EQ.1) THEN
C
C--INFLOW FROM BCs
      DO LM=1,NLAKES
        N=LM
C.......PRECIPITATION
        Q=PRCPLK(LM)*SURFA(LM)
        C=CPPT(LM,ICOMP)
        PRECLK=PRECLK+C*ABS(Q)
C.......RUNOFF
        Q=RNF(LM)
        C=CRNF(LM,ICOMP)
        RUNOFLK=RUNOFLK+C*ABS(Q)
C.......WITHDRAWAL/INJECTION
        Q=WTHDRW(LM)
        C=CAUG(LM,ICOMP)
        IF(WTHDRW(LM).LT.0.0) THEN
          WDRLINJLK=WDRLINJLK+C*ABS(Q)
        ELSE
          WDRLLK=WDRLLK+CLAKE(LM,ICOMP)*ABS(Q)
        ENDIF
C.......EVAPORATION
        IF(IETLAK.NE.0) THEN
          Q=EVAPLK(LM)*SURFA(LM)
          VOUTLAK(N)=VOUTLAK(N)+Q
          ETLK=ETLK+CLAKE(LM,ICOMP)*ABS(Q)
        ENDIF
      ENDDO
C
C--INFLOW FROM SFR - ASSUME SFR INFLOW WITH ZERO CONCENTRATION
      IF (IUNITSFR.GT.0) THEN
        DO LK=1,NLAKES
        DO ITRIB=1,NTRB
            INODE=ITRB(LK,ITRIB)
            IF (INODE.LE.0) GO TO 200
            Q=STRIN(INODE)
            C=0.
            SFR2LAK=SFR2LAK+C*ABS(Q)
        ENDDO
        ENDDO
      END IF
C
C--OUTFLOW TO SFR
      IF (IUNITSFR.GT.0) THEN
        DO LK=1,NLAKES
        DO IDV=1,NDV
            INODE=IDIV(LK,IDV)
            IF (INODE.GT.0) THEN
              Q=STROUT(INODE)
              SFRFROMLAK=SFRFROMLAK+CLAKE(LK,ICOMP)*ABS(Q)
            ENDIF
        ENDDO
        ENDDO
      END IF
C
C--UZF TO LAK FLOW
C      IF(iUnitTRNOP(7).GT.0) THEN
C        DO ICON=1,NLAK2UZF
C            N=INOD1LKUZ(ICON)
C            NUZ=INOD2LKUZ(ICON)
C            CALL NODE2KIJ(NUZ,NLAY,NROW,NCOL,K,I,J)
C            Q=QLAK2UZF(ICON)
C            Q=ABS(Q)
C            IF(IUZCODELK(ICON).EQ.1) THEN
C              CONC=CNEW(J,I,K,ICOMP)
C              RMASLAK(N)=RMASLAK(N)+CONC*ABS(Q)
C            ELSEIF(IUZCODELK(ICON).EQ.2) THEN
C              CONC=CUZINF(J,I,ICOMP)
C              RMASLAK(N)=RMASLAK(N)+CONC*ABS(Q)
C            ELSEIF(IUZCODELK(ICON).EQ.3) THEN
C              CONC=CUZINF(J,I,ICOMP)
C              RMASLAK(N)=RMASLAK(N)+CONC*ABS(Q)
C            ELSE
C              WRITE(IOUT,*) 'CHECK FTL FILE - IUZCODELK(ICON) INVALID'
C              WRITE(*,*) 'CHECK FTL FILE - IUZCODELK(ICON) INVALID'
C              STOP
C            ENDIF
C        ENDDO
C      ENDIF
C
C--ADD OLD MASS IN THE LAKE
C      IF(iSSTrans.EQ.0) THEN
      DO N=1,NLAKES
        CO=CLKOLD(N,ICOMP)
        VO=VOLOLD(N)
        TOTMASOLD=VO*CO
        CN=CLAKE(N,ICOMP)
        VN=VOL(N)
        TOTMASNEW=VN*CN
        STORDIFF=TOTMASNEW-TOTMASOLD
        IF(STORDIFF.LT.0) THEN
          STORINLK=STORINLK-STORDIFF/DELT
        ELSE
          STOROTLK=STOROTLK+STORDIFF/DELT
        ENDIF
      ENDDO
C      ENDIF
C
C
C--CUMULATIVE TERMS
      CGW2LAK(ICOMP)=CGW2LAK(ICOMP)+GW2LAK*DELT
      CUZF2LAK(ICOMP)=CUZF2LAK(ICOMP)+UZF2LAK*DELT
      CGWFROMLAK(ICOMP)=CGWFROMLAK(ICOMP)+GWFROMLAK*DELT
      CSFR2LAK(ICOMP)=CSFR2LAK(ICOMP)+SFR2LAK*DELT
      CSFRFROMLAK(ICOMP)=CSFRFROMLAK(ICOMP)+SFRFROMLAK*DELT
      CPRECLK(ICOMP)=CPRECLK(ICOMP)+PRECLK*DELT
      CRUNOFLK(ICOMP)=CRUNOFLK(ICOMP)+RUNOFLK*DELT
      CWDRLLK(ICOMP)=CWDRLLK(ICOMP)+WDRLLK*DELT
      CETLK(ICOMP)=CETLK(ICOMP)+ETLK*DELT
      CSTORINLK(ICOMP)=CSTORINLK(ICOMP)+STORINLK*DELT
      CSTOROTLK(ICOMP)=CSTOROTLK(ICOMP)+STOROTLK*DELT
C
C--CALCULATE TOTAL
      TOTINLK=GW2LAK+SFR2LAK+PRECLK+RUNOFLK+STORINLK+UZF2LAK
      TOTOUTLK=GWFROMLAK+SFRFROMLAK+WDRLLK+ETLK+STOROTLK
      CTOTINLK=CGW2LAK(ICOMP)+CSFR2LAK(ICOMP)+CPRECLK(ICOMP)+
     &         CRUNOFLK(ICOMP)+CSTORINLK(ICOMP)+CUZF2LAK(ICOMP)
      CTOTOUTLK=CGWFROMLAK(ICOMP)+CSFRFROMLAK(ICOMP)+
     &          CWDRLLK(ICOMP)+CETLK(ICOMP)+CSTOROTLK(ICOMP)
C
      DIFF=TOTINLK-TOTOUTLK
      CDIFF=CTOTINLK-CTOTOUTLK
      IF(ABS(TOTINLK+TOTOUTLK).LE.1.0E-10) TOTINLK=1.0E-10
      PERC=DIFF*100/((TOTINLK+TOTOUTLK)/2.0E0)
      IF(ABS(CTOTINLK+CTOTOUTLK).LE.1.0E-10) CTOTINLK=1.0E-10
      CPERC=CDIFF*100/((CTOTINLK+CTOTOUTLK)/2.0E0)
C
C--WRITE LAKE MASS BALANCE TO OUTPUT FILE
      IF(ISPCFL.NE.0) THEN
        WRITE(IOUT,10) KSTP,KPER,ICOMP
        WRITE(IOUT,20)
        WRITE(IOUT,29) CSTORINLK(ICOMP),STORINLK
        WRITE(IOUT,30) CGW2LAK(ICOMP),GW2LAK
C        IF(iUnitTRNOP(7).GT.0) WRITE(IOUT,31) CUZF2LAK(ICOMP),UZF2LAK
        IF(IUNITSFR.GT.0) WRITE(IOUT,35) CSFR2LAK(ICOMP),SFR2LAK
        WRITE(IOUT,40) CPRECLK(ICOMP),PRECLK
        WRITE(IOUT,41) CRUNOFLK(ICOMP),RUNOFLK
        WRITE(IOUT,43)
        WRITE(IOUT,45) CTOTINLK,TOTINLK
        WRITE(IOUT,49) CSTOROTLK(ICOMP),STOROTLK
        WRITE(IOUT,50) CGWFROMLAK(ICOMP),GWFROMLAK
        IF(IUNITSFR.GT.0) WRITE(IOUT,55) CSFRFROMLAK(ICOMP),
     1  SFRFROMLAK
        WRITE(IOUT,60) CWDRLLK(ICOMP),WDRLLK
        IF(IETLAK.NE.0) WRITE(IOUT,61) CETLK(ICOMP),ETLK
        WRITE(IOUT,43)
        WRITE(IOUT,65) CTOTOUTLK,TOTOUTLK
        WRITE(IOUT,70) CDIFF,DIFF
        WRITE(IOUT,75) CPERC,PERC
      ENDIF
10    FORMAT(//31X,'LAKE MASS BUDGETS AT END OF TIME STEP',
     &         I5,', STRESS PERIOD',I5,' FOR ',
     &            'COMPONENT',I4,
     &       /21X,101('-'))
20    FORMAT(/33X,7X,1X,'CUMULATIVE MASS [M]',
     &         8X,13X,15X,' MASS RATE [M/T]',
     &       /41X,19('-'),36X,16('-'))
29    FORMAT(16X,'        LAKE DEPLETION =',G15.7,
     &       16X,'        LAKE DEPLETION =',G15.7)
30    FORMAT(16X,'            GW TO LAKE =',G15.7,
     &       16X,'            GW TO LAKE =',G15.7)
31    FORMAT(16X,'           UZF TO LAKE =',G15.7,
     &       16X,'           UZF TO LAKE =',G15.7)
35    FORMAT(16X,'        STREAM TO LAKE =',G15.7,
     &       16X,'        STREAM TO LAKE =',G15.7)
40    FORMAT(16X,'         PRECIPITATION =',G15.7,
     &       16X,'         PRECIPITATION =',G15.7)
41    FORMAT(16X,'                RUNOFF =',G15.7,
     &       16X,'                RUNOFF =',G15.7)
43    FORMAT(41X,19('-'),36X,14('-'))
45    FORMAT(16X,'              TOTAL IN =',G15.7,
     &       16X,'              TOTAL IN =',G15.7)
49    FORMAT(/16X,'     LAKE ACCUMULATION =',G15.7,
     &        16X,'     LAKE ACCUMULATION =',G15.7)
50    FORMAT(16X,'            LAKE TO GW =',G15.7,
     &       16X,'            LAKE TO GW =',G15.7)
55    FORMAT(16X,'        LAKE TO STREAM =',G15.7,
     &       16X,'        LAKE TO STREAM =',G15.7)
60    FORMAT(16X,'            WITHDRAWAL =',G15.7,
     &       16X,'            WITHDRAWAL =',G15.7)
61    FORMAT(16X,'           EVAPORATION =',G15.7,
     &       16X,'           EVAPORATION =',G15.7)
65    FORMAT(16X,'             TOTAL OUT =',G15.7,
     &       16X,'             TOTAL OUT =',G15.7)
70    FORMAT(/16X,'        NET (IN - OUT) =',G15.7,
     &        16X,'        NET (IN - OUT) =',G15.7)
75    FORMAT(16X,' DISCREPANCY (PERCENT) =',G15.7,
     &       16X,' DISCREPANCY (PERCENT) =',G15.7,//)
C
      ENDIF
C
C20------RETURN
      RETURN
      END
C----------------------------------------------------------------------------
      SUBROUTINE CONCIAUX(ICOMP,NAUX,BNDAUX,IAUX)
C     ******************************************************************
C     GET INDEX FOR AUXILIARY VARIABLE FOR COMPONENT ICOMP
C     ******************************************************************
C
C        SPECIFICATIONS:
      USE GWTBCTMODULE, ONLY: MCOMP,IHEAT,MCOMPT
      CHARACTER*16 BNDAUX(20),BAUX
C     ------------------------------------------------------------------
      BAUX = BNDAUX(1)
      IAUX = 0
      DO I=1,NAUX
        IF(IHEAT.EQ.1.AND.ICOMP.EQ.MCOMPT) THEN
          IF(BNDAUX(I).EQ.'TMPR '.AND.ICOMP.EQ.MCOMPT)THEN
            IAUX = I
            RETURN
          ENDIF
        ENDIF
C
        IF(BNDAUX(I).EQ.'C01 '.AND.ICOMP.EQ.1)THEN
          IAUX = I
          RETURN
        ELSEIF(BNDAUX(I).EQ.'C02 '.AND.ICOMP.EQ.2)THEN
          IAUX = I
          RETURN
        ELSEIF(BNDAUX(I).EQ.'C03 '.AND.ICOMP.EQ.3)THEN
          IAUX = I
          RETURN
        ELSEIF(BNDAUX(I).EQ.'C04 '.AND.ICOMP.EQ.4)THEN
          IAUX = I
          RETURN
        ELSEIF(BNDAUX(I).EQ.'C05 '.AND.ICOMP.EQ.5)THEN
          IAUX = I
          RETURN
        ELSEIF(BNDAUX(I).EQ.'C06 '.AND.ICOMP.EQ.6)THEN
          IAUX = I
          RETURN
        ELSEIF(BNDAUX(I).EQ.'C07 '.AND.ICOMP.EQ.7)THEN
          IAUX = I
          RETURN
        ELSEIF(BNDAUX(I).EQ.'C08 '.AND.ICOMP.EQ.8)THEN
          IAUX = I
          RETURN
        ELSEIF(BNDAUX(I).EQ.'C09 '.AND.ICOMP.EQ.9)THEN
          IAUX = I
          RETURN
        ELSEIF(BNDAUX(I).EQ.'C10 '.AND.ICOMP.EQ.10)THEN
          IAUX = I
          RETURN
        ELSEIF (BNDAUX(I).EQ.'C11 '.AND.ICOMP.EQ.11)THEN
          IAUX = I
          RETURN
         ELSEIF (BNDAUX(I).EQ.'C12 '.AND.ICOMP.EQ.12)THEN
          IAUX = I
          RETURN
        ELSEIF (BNDAUX(I).EQ.'C13 '.AND.ICOMP.EQ.13)THEN
          IAUX = I
          RETURN
        ELSEIF (BNDAUX(I).EQ.'C14 '.AND.ICOMP.EQ.14)THEN
          IAUX = I
          RETURN
        ELSEIF (BNDAUX(I).EQ.'C15 '.AND.ICOMP.EQ.15)THEN
          IAUX = I
          RETURN
        ELSEIF (BNDAUX(I).EQ.'C16 '.AND.ICOMP.EQ.16)THEN
          IAUX = I
          RETURN
        ELSEIF (BNDAUX(I).EQ.'C17 '.AND.ICOMP.EQ.17)THEN
          IAUX = I
          RETURN
        ELSEIF (BNDAUX(I).EQ.'C18 '.AND.ICOMP.EQ.18)THEN
          IAUX = I
          RETURN
        ELSEIF (BNDAUX(I).EQ.'C19 '.AND.ICOMP.EQ.19)THEN
          IAUX = I
          RETURN
        ELSEIF (BNDAUX(I).EQ.'C20 '.AND.ICOMP.EQ.20)THEN
          IAUX = I
          RETURN
        ENDIF
      ENDDO
C
C ------RETURN
      RETURN
      END
C----------------------------------------------------------------------------

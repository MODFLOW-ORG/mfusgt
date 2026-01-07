      SUBROUTINE UBUDSVRD(KSTP,KPER,TEXT,IBDCHN,BUFF,NCOL,NROW,NLAY,
     1  IOUT)
C     ******************************************************************
C     READ CELL-BY-CELL FLOW TERMS FOR ONE COMPONENT OF FLOW.
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      CHARACTER*16 TEXT
      DIMENSION BUFF(NCOL,NROW,NLAY)
C     ------------------------------------------------------------------
C
C1------WRITE AN UNFORMATTED RECORD IDENTIFYING DATA.
      WRITE(IOUT,1) TEXT,IBDCHN,KSTP,KPER
    1 FORMAT(1X,'UBUDSVRD READING "',A16,'" ON UNIT',I4,
     1     ' AT TIME STEP',I7,', STRESS PERIOD ',I7)
      READ(IBDCHN) KSTP,KPER,TEXT,NCOL,NROW,NLAY
C
C2------WRITE AN UNFORMATTED RECORD CONTAINING VALUES FOR
C2------EACH CELL IN THE GRID.
      READ(IBDCHN) BUFF
C
C3------RETURN
      RETURN
      END
C---------------------------------------------------------------------------------- 
      SUBROUTINE UBDSV1RD(KSTP,KPER,TEXT,IBDCHN,BUFF,NCOL,NROW,NLAY,
     1          IOUT,DELT,PERTIM,TOTIM,IBOUND)
C     ******************************************************************
C     RECORD CELL-BY-CELL FLOW TERMS FOR ONE COMPONENT OF FLOW AS A 3-D
C     ARRAY WITH EXTRA RECORD TO INDICATE DELT, PERTIM, AND TOTIM.
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      CHARACTER*16 TEXT
      DOUBLE PRECISION PERTIM,TOTIM,DELT
      REAL PTIM,TTIM,DT
      DIMENSION BUFF(NCOL,NROW,NLAY),IBOUND(NCOL,NROW,NLAY)
C     ------------------------------------------------------------------
C
C0------ASSIGN SINGLE PRECISION VARIABLES
      DT=DELT
      PTIM=PERTIM
      TTIM=TOTIM
C
C1------WRITE TWO UNFORMATTED RECORDS IDENTIFYING DATA.
      IF(IOUT.GT.0) WRITE(IOUT,1) TEXT,IBDCHN,KSTP,KPER
    1 FORMAT(1X,'UBDSV1RD READING "',A16,'" ON UNIT',I4,
     1     ' AT TIME STEP',I7,', STRESS PERIOD',I7)
      READ(IBDCHN) KSTP,KPER,TEXT,NCOL,NROW,NLAY
      READ(IBDCHN) IONE,DT,PTIM,TTIM
C
C2------WRITE AN UNFORMATTED RECORD CONTAINING VALUES FOR
C2------EACH CELL IN THE GRID.
      READ(IBDCHN) BUFF
C
C3------RETURN
      RETURN
      END      
C----------------------------------------------------------------------------------      
      SUBROUTINE UBUDSVUrd(KSTP,KPER,TEXT,IBDCHN,BUFF,NJA,IOUT,
     1  PERTIM,TOTIM)
C     ******************************************************************
C     RECORD CELL-BY-CELL FLOW TERMS FOR UNSTRUCTURED FORMAT.
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      CHARACTER*16 TEXT
      DOUBLE PRECISION PERTIM,TOTIM
      DIMENSION BUFF(NJA)
C     ------------------------------------------------------------------
C
C1------WRITE AN UNFORMATTED RECORD IDENTIFYING DATA.
      WRITE(IOUT,1) TEXT,IBDCHN,KSTP,KPER
    1 FORMAT(1X,'UBUDSVUrd reading "',A16,'" ON UNIT',I4,
     1     ' AT TIME STEP',I7,', STRESS PERIOD',I7)
      read(IBDCHN) KSTPrd,KPERrd,TEXT,NJArd,i1,i2
C
C2------WRITE AN UNFORMATTED RECORD CONTAINING VALUES FOR
C2------EACH CELL IN THE GRID.
      read(IBDCHN) BUFF
C
C3------RETURN
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE UBDSV1Urd(KSTP,KPER,TEXT,IBDCHN,BUFF,NJA,IOUT,
     1          DELT,PERTIM,TOTIM,IBOUND,NODES)
C     ******************************************************************
C     RECORD CELL-BY-CELL FLOW TERMS FOR UNSTRUCTURED FORMAT,
C     WITH EXTRA RECORD TO INDICATE DELT, PERTIM, AND TOTIM.
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      CHARACTER*16 TEXT
      DOUBLE PRECISION PERTIM,TOTIM,DELT
      DIMENSION BUFF(NJA),IBOUND(NODES)
      REAL PTIM,TTIM,DT
C     ------------------------------------------------------------------
C
C1------WRITE TWO UNFORMATTED RECORDS IDENTIFYING DATA.
      IF(IOUT.GT.0) WRITE(IOUT,1) TEXT,IBDCHN,KSTP,KPER
    1 FORMAT(1X,'UBDSV1Urd READING "',A16,'" ON UNIT',I4,
     1     ' AT TIME STEP',I7,', STRESS PERIOD',I7)
      read(IBDCHN) KSTPrd,KPERrd,TEXT,NJArd,i1,i2
      read(IBDCHN) i1,DT,PTIM,TTIM
C
C2------WRITE AN UNFORMATTED RECORD CONTAINING VALUES FOR
C2------EACH CELL IN THE GRID.
      read(IBDCHN) BUFF
C
C3------RETURN
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE UBDSV2Urd(KSTP,KPER,TEXT,IBDCHN,NODES,
     1          NLIST,IOUT,DELT,PERTIM,TOTIM,IBOUND)
C     ******************************************************************
C     WRITE HEADER RECORDS FOR CELL-BY-CELL FLOW TERMS FOR ONE COMPONENT
C     OF FLOW USING A LIST STRUCTURE.  EACH ITEM IN THE LIST IS WRITTEN
C     BY MODULE UBDSVA
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      CHARACTER*16 TEXT
      DOUBLE PRECISION PERTIM,TOTIM,DELT
      DIMENSION IBOUND(NODES)
      REAL PTIM,TTIM,DT
C     ------------------------------------------------------------------
C
C1------WRITE THREE UNFORMATTED RECORDS IDENTIFYING DATA.
      IF(IOUT.GT.0) WRITE(IOUT,1) TEXT,IBDCHN,KSTP,KPER
    1 FORMAT(1X,'UBDSV2Urd READING "',A16,'" ON UNIT',I4,
     1     ' AT TIME STEP',I7,', STRESS PERIOD',I7)
      read(IBDCHN) KSTPrd,KPERrd,TEXT,NODESrd,i1,i2
      read(IBDCHN) i3,DT,PTIM,TTIM
      read(IBDCHN) NLIST
C
C2------RETURN
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE UBDSVAUrd(IBDCHN,NODES,N,Q,IBOUND)
C     ******************************************************************
C     WRITE ONE VALUE OF CELL-BY-CELL FLOW USING A LIST STRUCTURE.
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      DIMENSION IBOUND(NODES)
C     ------------------------------------------------------------------
C
C2------WRITE CELL NUMBER AND FLOW RATE
      read(IBDCHN) N,Q
C
C3------RETURN
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE UBDSV3Urd(KSTP,KPER,TEXT,IBDCHN,BUFF,IBUFF,IDIM,NOPT,
     1              NODES,IOUT,DELT,PERTIM,TOTIM,IBOUND)
C     ******************************************************************
C     READ CELL-BY-CELL FLOW TERMS FOR ONE COMPONENT OF FLOW AS A 2-D
C     ARRAY OF FLOW VALUES AND OPTIONALLY A 2-D ARRAY OF LAYER NUMBERS
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      CHARACTER*16 TEXT
      DOUBLE PRECISION PERTIM,TOTIM,DELT
      DIMENSION BUFF(NODES),IBUFF(IDIM),
     1          IBOUND(NODES)
      REAL PTIM,TTIM,DT
C     ------------------------------------------------------------------
C
C1------READ TWO UNFORMATTED RECORDS IDENTIFYING DATA.
      IF(IOUT.GT.0) WRITE(IOUT,1) TEXT,IBDCHN,KSTP,KPER
    1 FORMAT(1X,'UBDSV3Urd READING "',A16,'" ON UNIT',I4,
     1     ' AT TIME STEP',I7,', STRESS PERIOD',I7)
      read(IBDCHN) KSTP,KPER,TEXT,IDIM,i1,i2
      IMETH=3
      IF(NOPT.EQ.1) IMETH=4
      read(IBDCHN) IMETHrd,DT,PTIM,TTIM
C
C2------READ DATA AS ONE OR TWO UNFORMATTED RECORDS CONTAINING ONE
C2------VALUE PER LAYER.
      IF(NOPT.EQ.1) THEN
C2A-----READ ONE RECORD WHEN NOPT IS 1.  THE VALUES ARE FLOW VALUES
C2A-----FOR THE NODES IN BUFF.
         read(IBDCHN) (BUFF(IBUFF(IDI)),IDI = 1,IDIM)
      ELSE
C2B-----READ TWO RECORDS WHEN NOPT IS NOT 1.  FIRST RECORD CONTAINS
C2B-----NODE NUMBERS;  SECOND RECORD CONTAINS FLOW VALUES.
         read(IBDCHN) (IBUFF(IDI),IDI=1,IDIM)
         read(IBDCHN) (BUFF(IBUFF(IDI)),IDI = 1,IDIM)
      END IF
C
C3------RETURN
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE UBDSV4rd(KSTP,KPER,TEXT,NAUX,AUXTXT,IBDCHN,
     1          NCOL,NROW,NLAY,NLIST,IOUT,DELT,PERTIM,TOTIM,IBOUND)
C     ******************************************************************
C     READ HEADER RECORDS FOR CELL-BY-CELL FLOW TERMS FOR ONE COMPONENT
C     OF FLOW PLUS AUXILIARY DATA USING A LIST STRUCTURE.  EACH ITEM IN
C     THE LIST IS READ BY MODULE UBDSVB
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      CHARACTER*16 TEXT,AUXTXT(*)
      DOUBLE PRECISION PERTIM,TOTIM,DELT
      REAL DT,PTIM,DTIM
      DIMENSION IBOUND(NCOL,NROW,NLAY)
C     ------------------------------------------------------------------
C
C1------WRITE UNFORMATTED RECORDS IDENTIFYING DATA.
      IF(IOUT.GT.0) WRITE(IOUT,1) TEXT,IBDCHN,KSTP,KPER
    1 FORMAT(1X,'UBDSV4rd READING "',A16,'" ON UNIT',I4,
     1     ' AT TIME STEP',I7,', STRESS PERIOD',I7)
      read(IBDCHN) KSTPrd,KPERrd,TEXT,NCOLrd,NROWrd,NLAYrd
      read(IBDCHN) i5,DT,PTIM,TTIM
      read(IBDCHN) NAUXrd
      IF(NAUX.GT.0) read(IBDCHN) (AUXTXT(N),N=1,NAUX)
      read(IBDCHN) NLIST
C
C2------RETURN
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE UBDSV4Urd(KSTP,KPER,TEXT,NAUX,AUXTXT,IBDCHN,
     1          NODES,NLIST,IOUT,DELT,PERTIM,TOTIM,IBOUND)
C     ******************************************************************
C     READ HEADER RECORDS FOR CELL-BY-CELL FLOW TERMS FOR ONE COMPONENT
C     OF FLOW PLUS AUXILIARY DATA USING A LIST STRUCTURE FOR UNSTRUCTURED GRID.
C     EACH ITEM IN THE LIST IS READ BY MODULE UBDSVBU
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      CHARACTER*16 TEXT,AUXTXT(*)
      DOUBLE PRECISION PERTIM,TOTIM,DELT
      DIMENSION IBOUND(NODES)
      REAL PTIM,TTIM,DT
C     ------------------------------------------------------------------
C
C1------READ UNFORMATTED RECORDS IDENTIFYING DATA.
      IF(IOUT.GT.0) WRITE(IOUT,1) TEXT,IBDCHN,KSTP,KPER
    1 FORMAT(1X,'UBDSV4Urd READING "',A16,'" ON UNIT',I4,
     1     ' AT TIME STEP',I7,', STRESS PERIOD',I7)
      read(IBDCHN) KSTPrd,KPERrd,TEXT,NODESrd,i1,i2
      read(IBDCHN) i5,DT,PTIM,TTIM
      read(IBDCHN) NAUXrd
      IF(NAUX.GT.0) read(IBDCHN) (AUXTXT(N),N=1,NAUX)
      read(IBDCHN) NLIST
C
C2------RETURN
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE UBDSVBUrd(IBDCHN,NODES,NN,Q,VAL,NVL,NAUX,LAUX,IBOUND)
C     ******************************************************************
C     READ ONE VALUE OF CELL-BY-CELL FLOW PLUS AUXILIARY DATA USING
C     A LIST FOR UNSTRUCTURED GRID.
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
       DIMENSION VAL(NVL),IBOUND(NODES)
C     ------------------------------------------------------------------
C
C1------READ CELL NUMBER AND FLOW RATE
      IF(NAUX.GT.0) THEN
         N2=LAUX+NAUX-1
         READ(IBDCHN) NN,Q,(VAL(N),N=LAUX,N2)
      ELSE
         READ(IBDCHN) NN,Q
      END IF
C
C2------RETURN
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE UBDSVBrd(IBDCHN,NCOL,NROW,J,I,K,Q,VAL,NVL,NAUX,LAUX,
     1                  IBOUND,NLAY)
C     ******************************************************************
C     READ ONE VALUE OF CELL-BY-CELL FLOW PLUS AUXILIARY DATA USING
C     A LIST STRUCTURE.
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      DIMENSION IBOUND(NCOL,NROW,NLAY),VAL(NVL)
C     ------------------------------------------------------------------
C
C1------CALCULATE CELL NUMBER
      ICRL= (K-1)*NROW*NCOL + (I-1)*NCOL + J
C
C2------READ CELL NUMBER AND FLOW RATE
      IF(NAUX.GT.0) THEN
         N2=LAUX+NAUX-1
         read(IBDCHN) ICRLrd,Q,(VAL(N),N=LAUX,N2)
      ELSE
         read(IBDCHN) ICRLrd,Q
      END IF
C
C3------RETURN
      RETURN
      END
C----------------------------------------------------------------------------------
      SUBROUTINE UBDSVHDRrd(IUNSTR,KSTP,KPER,IOUT,IBNDCB,ICLNCB,NODES,
     1  NCLNNDS,NCOL,NROW,NLAY,NBND,NBNDVL,NAUX,IBOUND,
     2  TEXT,BNDAUX,DELT,PERTIM,TOTIM,BND)
C     ******************************************************************
C     READ THE COMPACT BUDGET HEADER FOR EITHER STRUCTURED OR
C     UNSTRUCTURED MODELS.
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
        IMPLICIT NONE
        INTEGER,INTENT(IN)::IUNSTR
        INTEGER,INTENT(IN)::KSTP
        INTEGER,INTENT(IN)::KPER
        INTEGER,INTENT(IN)::IOUT
        INTEGER,INTENT(IN)::IBNDCB
        INTEGER,INTENT(IN)::ICLNCB
        INTEGER,INTENT(IN)::NODES
        INTEGER,INTENT(IN)::NCLNNDS
        INTEGER,INTENT(IN)::NCOL
        INTEGER,INTENT(IN)::NROW
        INTEGER,INTENT(IN)::NLAY
        INTEGER,INTENT(IN)::NBND
        INTEGER,INTENT(IN)::NBNDVL
        INTEGER,INTENT(IN)::NAUX
        INTEGER,DIMENSION(NODES+NCLNNDS),INTENT(IN)::IBOUND
        CHARACTER(LEN=16),INTENT(IN)::TEXT
        CHARACTER(LEN=16),DIMENSION(NAUX),INTENT(IN)::BNDAUX
        DOUBLE PRECISION,INTENT(IN)::DELT,PERTIM,TOTIM
        REAL,DIMENSION(NBNDVL,NBND),INTENT(IN)::BND
        INTEGER::NBNDGW,NBNDCLN,L,N
C     ------------------------------------------------------------------
C
C1------SEPARATE GW BOUNDARIES AND CLN BOUNDARIES IF NEEDED
      IF(ICLNCB.EQ.IBNDCB) THEN
        NBNDGW = NBND
      ELSE
C1------COUNT GW BOUNDARIES AND CLN BOUNDARIES
        NBNDGW = 0
        NBNDCLN = 0
        DO L=1,NBND
          N=BND(1,L)
          IF(N.GT.NODES) THEN
            NBNDCLN = NBNDCLN + 1
          ELSE
            NBNDGW = NBNDGW + 1
          ENDIF
        ENDDO
      ENDIF
C
C3------IF STRUCTURED, WRITE STRUCTURED HEADER, OTHERWISE WRITE
C       UNSTRUTURED HEADER
        IF(IUNSTR.EQ.0)THEN
C
C3A------WRITE STRUCTURED HEADER
          CALL UBDSV4rd(KSTP,KPER,TEXT,NAUX,BNDAUX,IBNDCB,NCOL,NROW,
     1                NLAY,NBNDGW,IOUT,DELT,PERTIM,TOTIM,IBOUND)
        ELSE
C
C3A------WRITE UNSTRUCTURED HEADER
          CALL UBDSV4Urd(KSTP,KPER,TEXT,NAUX,BNDAUX,IBNDCB,NODES,
     1                 NBNDGW,IOUT,DELT,PERTIM,TOTIM,IBOUND)
        ENDIF
C
C4------WRITE CLN BUDGET FILE HEADER IF NOT SAME UNIT NUMBER AS GWF
        IF(ICLNCB.EQ.IBNDCB) RETURN
C
        IF(ICLNCB.GT.0)
     1    CALL UBDSV4Urd(KSTP,KPER,TEXT,NAUX,BNDAUX,ICLNCB,NCLNNDS,
     2                 NBNDCLN,IOUT,DELT,PERTIM,TOTIM,IBOUND)
C
C5------RETURN
      RETURN
      END
C------------------------------------------------------------------------------
      SUBROUTINE UBDSVHDRTrd(IUNSTR,KSTP,KPER,IOUT,IBNDCB,ICLNCB,NODES,
     1  NCLNNDS,NCOL,NROW,NLAY,NBND,NBNDVL,NAUX,IBOUND,
     2  TEXT,BNDAUX,DELT,PERTIM,TOTIM,IDRTFL,BND)
C     ******************************************************************
C     READ THE COMPACT BUDGET HEADER FOR EITHER STRUCTURED OR
C     UNSTRUCTURED MODELS.
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
        IMPLICIT NONE
        INTEGER,INTENT(IN)::IUNSTR
        INTEGER,INTENT(IN)::KSTP
        INTEGER,INTENT(IN)::KPER
        INTEGER,INTENT(IN)::IOUT
        INTEGER,INTENT(IN)::IBNDCB
        INTEGER,INTENT(IN)::ICLNCB
        INTEGER,INTENT(IN)::NODES
        INTEGER,INTENT(IN)::NCLNNDS
        INTEGER,INTENT(IN)::NCOL
        INTEGER,INTENT(IN)::NROW
        INTEGER,INTENT(IN)::NLAY
        INTEGER,INTENT(IN)::NBND
        INTEGER,INTENT(IN)::NBNDVL
        INTEGER,INTENT(IN)::NAUX
        INTEGER,DIMENSION(NODES+NCLNNDS),INTENT(IN)::IBOUND
        CHARACTER(LEN=16),INTENT(IN)::TEXT
        CHARACTER(LEN=16),DIMENSION(NAUX),INTENT(IN)::BNDAUX
        DOUBLE PRECISION,INTENT(IN)::DELT,PERTIM,TOTIM
        REAL,DIMENSION(NBNDVL,NBND),INTENT(IN)::BND
        INTEGER::NBNDGW,NBNDCLN,L,N,IDRTFL
C     ------------------------------------------------------------------
C
C1------SEPARATE GW BOUNDARIES AND CLN BOUNDARIES IF NEEDED
      IF(ICLNCB.EQ.IBNDCB) THEN
        NBNDGW = NBND
      ELSE
C1------COUNT GW BOUNDARIES AND CLN BOUNDARIES
        NBNDGW = 0
        NBNDCLN = 0
        DO L=1,NBND
          N=BND(1,L)
          IF(N.GT.NODES) THEN
            NBNDCLN = NBNDCLN + 1
          ELSE
            NBNDGW = NBNDGW + 1
          ENDIF
          IF(IDRTFL.GT.0) THEN
            N=BND(6,L)
            IF(N.GT.NODES) THEN
              NBNDCLN = NBNDCLN + 1
            ELSE
              NBNDGW = NBNDGW + 1
            ENDIF
          ENDIF
        ENDDO
      ENDIF
C
C3------IF STRUCTURED, WRITE STRUCTURED HEADER, OTHERWISE WRITE
C       UNSTRUTURED HEADER
        IF(IUNSTR.EQ.0)THEN
C
C3A------READ STRUCTURED HEADER
          CALL UBDSV4rd(KSTP,KPER,TEXT,NAUX,BNDAUX,IBNDCB,NCOL,NROW,
     1                NLAY,NBNDGW,IOUT,DELT,PERTIM,TOTIM,IBOUND)
        ELSE
C
C3A------READ UNSTRUCTURED HEADER
          CALL UBDSV4Urd(KSTP,KPER,TEXT,NAUX,BNDAUX,IBNDCB,NODES,
     1                 NBNDGW,IOUT,DELT,PERTIM,TOTIM,IBOUND)
        ENDIF
C
C4------READ CLN BUDGET FILE HEADER IF NOT SAME UNIT NUMBER AS GWF
        IF(ICLNCB.EQ.IBNDCB) RETURN
C
        IF(ICLNCB.GT.0)
     1    CALL UBDSV4Urd(KSTP,KPER,TEXT,NAUX,BNDAUX,ICLNCB,NCLNNDS,
     2                 NBNDCLN,IOUT,DELT,PERTIM,TOTIM,IBOUND)
C
C5------RETURN
      RETURN
      END
C------------------------------------------------------------------------------
      SUBROUTINE UBDSVHDRQRTrd(IUNSTR,KSTP,KPER,IOUT,IBNDCB,ICLNCB,
     1  NODES,NCLNNDS,NCOL,NROW,NLAY,NBND,NBNDVL,NAUX,IBOUND,
     2  TEXT,BNDAUX,DELT,PERTIM,TOTIM,IQRTFL,NodQRT,MXRTCELLS,BND)
C     ******************************************************************
C     READ THE COMPACT BUDGET HEADER FOR EITHER STRUCTURED OR
C     UNSTRUCTURED MODELS.
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
        IMPLICIT NONE
        INTEGER,INTENT(IN)::IUNSTR
        INTEGER,INTENT(IN)::KSTP
        INTEGER,INTENT(IN)::KPER
        INTEGER,INTENT(IN)::IOUT
        INTEGER,INTENT(IN)::IBNDCB
        INTEGER,INTENT(IN)::ICLNCB
        INTEGER,INTENT(IN)::NODES
        INTEGER,INTENT(IN)::NCLNNDS
        INTEGER,INTENT(IN)::NCOL
        INTEGER,INTENT(IN)::NROW
        INTEGER,INTENT(IN)::NLAY
        INTEGER,INTENT(IN)::NBND
        INTEGER,INTENT(IN)::NBNDVL
        INTEGER,INTENT(IN)::NAUX
        INTEGER,DIMENSION(NODES+NCLNNDS),INTENT(IN)::IBOUND
        CHARACTER(LEN=16),INTENT(IN)::TEXT
        CHARACTER(LEN=16),DIMENSION(NAUX),INTENT(IN)::BNDAUX
        DOUBLE PRECISION,INTENT(IN)::DELT,PERTIM,TOTIM
        REAL,DIMENSION(NBNDVL,NBND),INTENT(IN)::BND
        INTEGER, DIMENSION (MXRTCELLS),INTENT(IN):: NODQRT
        INTEGER::NBNDGW,NBNDCLN,L,N,IRT,IQRTFL,NUMRT,JJ,INR,
     *           MXRTCELLS
C     ------------------------------------------------------------------
C
C1------SEPARATE GW BOUNDARIES AND CLN BOUNDARIES IF NEEDED
      IF(ICLNCB.EQ.IBNDCB) THEN
        NBNDGW = NBND
      ELSE
C1------COUNT GW BOUNDARIES AND CLN BOUNDARIES
        NBNDGW = 0
        NBNDCLN = 0
        IRT = 0
        DO L=1,NBND
          N=BND(1,L)
          IF(N.GT.NODES) THEN
            NBNDCLN = NBNDCLN + 1
          ELSE
            NBNDGW = NBNDGW + 1
          ENDIF
          IF(IQRTFL.GT.0) THEN
            NumRT = BND(5,L)
            IF (NumRT.EQ.0) CYCLE
            DO JJ = 1,NumRT
              IRT = IRT + 1
              INR = NodQRT(IRT)
              IF(INR.GT.NODES) THEN
                NBNDCLN = NBNDCLN + 1
              ELSE
                NBNDGW = NBNDGW + 1
              ENDIF
            ENDDO
          ENDIF
        ENDDO
      ENDIF
C
C3------IF STRUCTURED, READ STRUCTURED HEADER, OTHERWISE READ
C       UNSTRUTURED HEADER
        IF(IUNSTR.EQ.0)THEN
C
C3A------READ STRUCTURED HEADER
          CALL UBDSV4rd(KSTP,KPER,TEXT,NAUX,BNDAUX,IBNDCB,NCOL,NROW,
     1                NLAY,NBNDGW,IOUT,DELT,PERTIM,TOTIM,IBOUND)
        ELSE
C
C3A------READ UNSTRUCTURED HEADER
          CALL UBDSV4Urd(KSTP,KPER,TEXT,NAUX,BNDAUX,IBNDCB,NODES,
     1                 NBNDGW,IOUT,DELT,PERTIM,TOTIM,IBOUND)
        ENDIF
C
C4------READ CLN BUDGET FILE HEADER IF NOT SAME UNIT NUMBER AS GWF
        IF(ICLNCB.EQ.IBNDCB) RETURN
C
        IF(ICLNCB.GT.0)
     1    CALL UBDSV4Urd(KSTP,KPER,TEXT,NAUX,BNDAUX,ICLNCB,NCLNNDS,
     2                 NBNDCLN,IOUT,DELT,PERTIM,TOTIM,IBOUND)
C
C5------RETURN
      RETURN
      END
C---------------------------------------------------------------------------
      SUBROUTINE UBDSVRECrd(IUNSTR,N,NODES,NCLNNDS,IBNDCB,ICLNCB,NBNDVL,
     1  LAUX,NAUX,Q,VAL,IBOUND,NCOL,NROW,NLAY)
C     ******************************************************************
C     READ A SINGLE RECORD IN THE COMPACT BUDGET FORM.  USE THE NODE
C     NUMBER TO DETERMINE IF IT IS A GW OR CLN NODE AND WRITE THE
C     RECORD TO THE APPROPRIATE GW OR CLN BUDGET FILE.
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
        IMPLICIT NONE
        INTEGER,INTENT(IN)::IUNSTR
        INTEGER,INTENT(IN)::N
        INTEGER,INTENT(IN)::NODES
        INTEGER,INTENT(IN)::NCLNNDS
        INTEGER,INTENT(IN)::IBNDCB
        INTEGER,INTENT(IN)::ICLNCB
        INTEGER,INTENT(IN)::NBNDVL
        INTEGER,INTENT(IN)::LAUX
        INTEGER,INTENT(IN)::NAUX
        REAL,INTENT(IN)::Q
        REAL,DIMENSION(NBNDVL),INTENT(IN)::VAL
        INTEGER,DIMENSION(NODES+NCLNNDS),INTENT(IN)::IBOUND
        INTEGER,INTENT(IN)::NCOL
        INTEGER,INTENT(IN)::NROW
        INTEGER,INTENT(IN)::NLAY
        INTEGER::IL,IJ,IR,IC
C     ------------------------------------------------------------------
C
C1------IF CLN NODE, READ TO CLN BUDGET
        IF(N.GT.NODES) THEN
          IF(ICLNCB.GT.0) THEN
            CALL UBDSVBUrd(ICLNCB,NCLNNDS,N-NODES,Q,VAL,NBNDVL,NAUX,
     1                 LAUX,IBOUND)
          ENDIF
        ELSE
C
C2------READ STRUCTURED GW RECORD
          IF(IUNSTR.EQ.0)THEN
            IL = (N-1) / (NCOL*NROW) + 1
            IJ = N - (IL-1)*NCOL*NROW
            IR = (IJ-1)/NCOL + 1
            IC = IJ - (IR-1)*NCOL
            CALL UBDSVBrd(IBNDCB,NCOL,NROW,IC,IR,IL,Q,
     1                  VAL,NBNDVL,NAUX,LAUX,IBOUND,NLAY)
          ELSE
C
C3------READ UNSTRUCTURED GW RECORD
            CALL UBDSVBUrd(IBNDCB,NODES,N,Q,VAL,NBNDVL,NAUX,LAUX,IBOUND)
          ENDIF
        ENDIF
C
C4------RETURN
      RETURN
      END



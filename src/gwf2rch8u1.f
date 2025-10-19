      MODULE GWFRCHMODULE
        INTEGER, SAVE, POINTER                 ::NRCHOP,IRCHCB,MXNDRCH
        INTEGER,SAVE, POINTER ::NPRCH,IRCHPF,INIRCH,NIRCH,mxznrch,ISELEV
        INTEGER, SAVE,POINTER::IPONDOPT,IRTSOPT,IRTSRD,INRTS,ICONCRCHOPT
        REAL,    SAVE,   DIMENSION(:),  ALLOCATABLE ::RECH,SELEV,RECHSV
        INTEGER, SAVE,   DIMENSION(:),  ALLOCATABLE ::IRCH,iznrch
        REAL,    SAVE,   DIMENSION(:),  ALLOCATABLE ::rtsrch
        REAL,    SAVE,   DIMENSION(:),  ALLOCATABLE      ::RCHF
        REAL,    SAVE,   DIMENSION(:,:),  ALLOCATABLE      ::RCHCONC  !DIMENSION NEQS,NCONCRCH
        INTEGER, SAVE,   DIMENSION(:),  ALLOCATABLE ::IRCHCONC     !DIMENSION MCOMP
        DOUBLE PRECISION, SAVE, POINTER :: tstartrch,tendrch,factrrch,
     1        TIMRCH
       END MODULE GWFRCHMODULE


      SUBROUTINE GWF2RCH8U1AR(IN,INBCT)
C     ******************************************************************
C     ALLOCATE ARRAY STORAGE FOR RECHARGE
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,      ONLY:IOUT,NCOL,NROW,IFREFM,NODLAY,IUNSTR
      USE GWFBASMODULE, ONLY: IATS
      USE GWFRCHMODULE,ONLY:NRCHOP,IRCHCB,NPRCH,IRCHPF,RECH,IRCH,
     1  MXNDRCH,INIRCH,NIRCH,SELEV,iznrch,mxznrch,ISELEV,IPONDOPT,
     1  IRTSOPT,ICONCRCHOPT,INRTS,TIMRCH,IRTSRD,RECHSV,RCHF,
     1  RCHCONC,IRCHCONC,ICONCRCHOPT
      USE GWTBCTMODULE, ONLY: MCOMPT
C
      CHARACTER*400 LINE
      CHARACTER*4 PTYP
C     ------------------------------------------------------------------
C
C1-------ALLOCATE SCALAR VARIABLES.
      ALLOCATE(NRCHOP,IRCHCB,MXNDRCH)
      ALLOCATE(NPRCH,IRCHPF,INIRCH,NIRCH)
      ALLOCATE(INRTS,TIMRCH,IRTSRD)
      IRTSRD=0
C
C2------IDENTIFY PACKAGE.
      IRCHPF=0
      WRITE(IOUT,1)IN
    1 FORMAT(1X,/1X,'RCH -- RECHARGE PACKAGE, VERSION 7, 5/2/2005',
     1' INPUT READ FROM UNIT ',I4)
C
C3------READ NRCHOP AND IRCHCB.
      CALL URDCOM(IN,IOUT,LINE)
      CALL UPARARRAL(IN,IOUT,LINE,NPRCH)
      IF(IFREFM.EQ.0) THEN
         READ(LINE,'(2I10)') NRCHOP,IRCHCB
      ELSE
         LLOC=1
         CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NRCHOP,R,IOUT,IN)
         CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IRCHCB,R,IOUT,IN)
      END IF
C
C3B------READ KEYWORD OPTIONS SEEPELEV, RTS AND RECHARGE CONC.
      ALLOCATE(IPONDOPT,IRTSOPT,ICONCRCHOPT)
      IPONDOPT=0
      IRTSOPT=0
      ICONCRCHOPT=0
      LLOC=1
   10 CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,N,R,IOUT,IN)
C3B1------FOR SEEPAGE-FACE ELEVATION
      IF(LINE(ISTART:ISTOP).EQ.'SEEPELEV') THEN
        WRITE(IOUT,13)
   13   FORMAT(1X,'SEEPAGE-FACE ELEVATIONS WILL BE READ.',
     1      '  VARIABLE INSELEV REQUIRED IN RECORD 5.')
        IPONDOPT = 1
      END IF
C3B2------FOR RTS
      IF(LINE(ISTART:ISTOP).EQ.'RTS') THEN
C3B2A-----CHECK TO SEE IF ATS IS ON. OR ELSE WRITE WARNING AND STOP
        IF(IATS.EQ.0)THEN
          WRITE(IOUT,15)
          STOP
        ENDIF
15      FORMAT(1X,'TRANSIENT RECHARGE NEEDS ADAPTIVE TIME-STEPPING.',
     1     'STOPPING')
C3B2B------SET OPTION, AND READ MAXIMUM NUMBER OF ZONES OF TRANSIENT RCH.
        ALLOCATE(MXZNRCH)
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,MXZNRCH,R,IOUT,INOC)
        WRITE(IOUT,14)MXZNRCH
   14   FORMAT(1X,'TRANSIENT RECHARGE FILE WITH',I8,' ZONES WILL BE',
     1      ' READ. RECHARGE ZONE INDICES WILL BE READ FROM RCH FILE.'
     2         /1X,107('-'))
        IRTSOPT = 1
      END IF
C3BC------FOR CONCENTRATION OF RECHARGE
      IF(LINE(ISTART:ISTOP).EQ.'CONCENTRATION' .OR.
     1   LINE(ISTART:ISTOP).EQ.'CONC' ) THEN
        WRITE(IOUT,16)
   16   FORMAT(1X,'SPECIES CONCENTRATION WILL BE READ.',
     1      '  ARRAY IRCHCONC REQUIRED TO INDICATE SPECIES')
        ICONCRCHOPT = 1
      END IF
      IF(LLOC.LT.200) GO TO 10
C3C------READ NUMBER OF RECHARGE NODES IF UNSTRUCTURED AND NRCHOP=2
      IF(IUNSTR.EQ.1.AND.NRCHOP.EQ.2)THEN
        READ(IN,*) MXNDRCH
      ELSE
        MXNDRCH = NODLAY(1)
      ENDIF
C3D-----ALLOCATE ZONAL ARRAY
        IF(IRTSOPT.EQ.1)THEN
          ALLOCATE(iznrch(mxndrch))
        ELSE
          ALLOCATE(iznrch(1))
        ENDIF
C3D-----ALLOCATE CONCENTRATION ARRAY
        IF(ICONCRCHOPT.GT.0)THEN
          ALLOCATE(IRCHCONC(MCOMPT))
C3D1------READ INDEX ARRAY FOR COMPONENTS WHOSE CONC IS READ
          READ(IN,*)(IRCHCONC(I),I=1,MCOMPT)
          ICONCRCH = 0
          DO II=1,MCOMPT
            ICONCRCH = ICONCRCH + IRCHCONC(II)
          ENDDO
          ALLOCATE(RCHCONC(mxndrch,ICONCRCH))
C3D1------READ ARRAY OF COMPONENT NUMBERS
        ELSE
          ALLOCATE(IRCHCONC(1))
          ALLOCATE(RCHCONC(1,1))
        ENDIF
C
C4------CHECK TO SEE THAT OPTION IS LEGAL.
      IF(NRCHOP.LT.1.OR.NRCHOP.GT.3) THEN
        WRITE(IOUT,8) NRCHOP
    8   FORMAT(1X,'ILLEGAL RECHARGE OPTION CODE (NRCHOP = ',I5,
     &       ') -- SIMULATION ABORTING')
        CALL USTOP(' ')
      END IF
C
C5------OPTION IS LEGAL -- PRINT OPTION CODE.
      IF(NRCHOP.EQ.1) WRITE(IOUT,201)
  201 FORMAT(1X,'OPTION 1 -- RECHARGE TO TOP LAYER')
      IF(NRCHOP.EQ.2) WRITE(IOUT,202)
  202 FORMAT(1X,'OPTION 2 -- RECHARGE TO ONE SPECIFIED NODE IN EACH',
     1     ' VERTICAL COLUMN')
      IF(NRCHOP.EQ.3) WRITE(IOUT,203)
  203 FORMAT(1X,'OPTION 3 -- RECHARGE TO HIGHEST ACTIVE NODE IN',
     1     ' EACH VERTICAL COLUMN')
C
C6------IF CELL-BY-CELL FLOWS ARE TO BE SAVED, THEN PRINT UNIT NUMBER.
      IF(IRCHCB.GT.0) WRITE(IOUT,204) IRCHCB
  204 FORMAT(1X,'CELL-BY-CELL FLOWS WILL BE SAVED ON UNIT ',I4)
C
C7------ALLOCATE SPACE FOR THE RECHARGE (RECH) AND INDICATOR (IRCH)
C7------ARRAYS.
      ALLOCATE (RECH(MXNDRCH))
      ALLOCATE (IRCH(MXNDRCH))
C8------IF TRANSPORT IS ACTIVE THEN ALLOCATE ARRAY TO STORE FLUXES
      IF(INBCT.GT.0)THEN
        ALLOCATE (RCHF(MXNDRCH))
      ENDIF
C--------ALLOCATE SPACE TO SAVE SEEPAGE-FACE INFORMATION
      ALLOCATE (ISELEV)
      ISELEV = 0
      INSELEV = 0
      IF (IPONDOPT.GT.0) THEN
        ALLOCATE(SELEV(mxndrch))
        DO I=1,MXNDRCH
          SELEV(I) = 1.0E20
        ENDDO
      ELSE
          ALLOCATE(SELEV(1))
      ENDIF
C---------ALLOCATE SPACE TO SAVE ORIGINAL RECH ARRAY FROM STRESS PERIODS
      IF(IRTSOPT.GT.0)THEN
        ALLOCATE (RECHSV(MXNDRCH))
      ELSE
        ALLOCATE (RECHSV(1))
      ENDIF
C
C8------READ NAMED PARAMETERS
      WRITE(IOUT,5) NPRCH
    5 FORMAT(1X,//1X,I5,' Recharge parameters')
      IF(NPRCH.GT.0) THEN
         DO 20 K=1,NPRCH
         CALL UPARARRRP(IN,IOUT,N,0,PTYP,1,1,0)
         IF(PTYP.NE.'RCH') THEN
            WRITE(IOUT,7)
    7       FORMAT(1X,'Parameter type must be RCH')
            CALL USTOP(' ')
         END IF
   20    CONTINUE
      END IF
C
C9------RETURN
      RETURN
      END
      SUBROUTINE GWF2RCH8U1RP(IN,IURTS,KPER)
C     ******************************************************************
C     READ RECHARGE DATA FOR STRESS PERIOD
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,      ONLY:IOUT,NCOL,NROW,NLAY,IFREFM,DELR,DELC,
     1  NODLAY,AREA,IUNSTR,NODES
      USE GWFRCHMODULE,ONLY:NRCHOP,NPRCH,IRCHPF,RECH,IRCH,INIRCH,NIRCH,
     *  SELEV,iznrch,mxznrch,ISELEV,IPONDOPT,IRTSOPT,RECHSV,ICONCRCHOPT,
     *  tstartrch,tendrch,factrrch,RTSRCH,INRTS,IRTSRD,TIMRCH,
     *  RCHCONC,IRCHCONC
      USE GWTBCTMODULE, ONLY: MCOMPT
      REAL, DIMENSION(:,:),ALLOCATABLE  ::TEMP
      INTEGER, DIMENSION(:,:),ALLOCATABLE  ::ITEMP
C
      CHARACTER*24 ANAME(5)
      CHARACTER(LEN=200) line
C
      DATA ANAME(1) /'    RECHARGE LAYER INDEX'/
      DATA ANAME(2) /'                RECHARGE'/
      DATA ANAME(3) /'                   SELEV'/
      DATA ANAME(4) /'                  iznrch'/
      DATA ANAME(5) /'                    CONC'/
C     ------------------------------------------------------------------
      ALLOCATE (TEMP(NCOL,NROW))
      ALLOCATE (ITEMP(NCOL,NROW))
C2------IDENTIFY PACKAGE.
      WRITE(IOUT,1)IN
    1 FORMAT(1X,/1X,'RCH -- RECHARGE PACKAGE, VERSION 7, 5/2/2005',
     1' INPUT READ FROM UNIT ',I4)
C
C2------READ FLAGS SHOWING WHETHER DATA IS TO BE REUSED.
      lloc = 1
      iniznrch=0
      INSELEV=0
      INCONC=0
      CALL URDCOM(In, Iout, line)
C3------GET OPTIONS FIRST
   10 CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,N,R,IOUT,IN)
      IF(LINE(ISTART:ISTOP).EQ.'INRCHZONES') THEN
C3B------READ KEYWORD OPTION FOR RTS ZONES TO BE READ.
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,INIZNRCH,R,IOUT,INOC)
        WRITE(IOUT,14) INIZNRCH
14      FORMAT(/1X,'FLAG FOR INPUT OF RTS ZONES (INIZNRCH) = ',
     1        I8)
      ELSEIF(LINE(ISTART:ISTOP).EQ.'INSELEV') THEN
C3C------IS KEWORD OPTION FOR SEEPAGE ELEVATION TO BE READ
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,INSELEV,R,IOUT,INOC)
        WRITE(IOUT,15) INSELEV
15      FORMAT(/1X,'FLAG FOR INPUT OF SEEPAGE ELEVATIONS (INSELEV) = ',
     1        I8)
       ELSEIF(LINE(ISTART:ISTOP).EQ.'INCONC') THEN
C3C------IS KEWORD OPTION FOR CONCENTRATION TO BE READ
        INCONC = 1
        WRITE(IOUT,16) INCONC
16      FORMAT(/1X,'FLAG FOR INPUT OF CONCENTRATIONS (INCONC) = ',
     1        I8)
      END IF
      IF(LLOC.LT.200) GO TO 10
      LLOC = 1
C3D------READ FLAGS
      IF(IFREFM.EQ.0)THEN
        IF(NRCHOP.EQ.2) THEN
          READ(LINE,'(2I10)') INRECH,INIRCH
        ELSE
          READ(LINE,'(I10)') INRECH
          INIRCH = NODLAY(1)
        ENDIF
      ELSE
        IF(NRCHOP.EQ.2) THEN
          CALL URWORD(line, lloc, istart, istop, 2, inrech, r, Iout, In)
          CALL URWORD(line, lloc, istart, istop, 2, inirch, r, Iout, In)
        ELSE
          CALL URWORD(line, lloc, istart, istop, 2, inrech, r, Iout, In)
          INIRCH = NODLAY(1)
        ENDIF
      END IF
      IF(INIRCH.GE.0) NIRCH = INIRCH
      IF(INSELEV.GE.0) ISELEV = INSELEV
C
C3------TEST INRECH TO SEE HOW TO DEFINE RECH.
      IF(INRECH.LT.0) THEN
C
C3A-----INRECH<0, SO REUSE RECHARGE ARRAY FROM LAST STRESS PERIOD.
        WRITE(IOUT,3)
    3   FORMAT(1X,/1X,'REUSING RECH FROM LAST STRESS PERIOD')
      ELSE
        IF(IUNSTR.EQ.0)THEN
C
C3B-----INRECH=>0, SO READ RECHARGE RATE.
          IF(NPRCH.EQ.0) THEN
C
C3B1--------THERE ARE NO PARAMETERS, SO READ RECH USING U2DREL.
            CALL U2DREL(TEMP,ANAME(2),NROW,NCOL,0,IN,IOUT)
          ELSE
C3B2--------DEFINE RECH USING PARAMETERS.  INRECH IS THE NUMBER OF
C3B2--------PARAMETERS TO USE THIS STRESS PERIOD.
            CALL PRESET('RCH')
            WRITE(IOUT,33)
   33       FORMAT(1X,///1X,
     1      'RECH array defined by the following parameters:')
            IF(INRECH.EQ.0) THEN
              WRITE(IOUT,34)
   34         FORMAT(' ERROR: When parameters are defined for the RCH',
     &      ' Package, at least one parameter',/,' must be specified',
     &      ' each stress period -- STOP EXECUTION (GWF2RCH8U1RPLL)')
              CALL USTOP(' ')
            END IF
            CALL UPARARRSUB2(TEMP,NCOL,NROW,0,INRECH,IN,IOUT,'RCH',
     1            ANAME(2),'RCH',IRCHPF)
          END IF
          N=0
          DO I=1,NROW
          DO J=1,NCOL
            N=N+1
            RECH(N)=TEMP(J,I)
          ENDDO
          ENDDO
        ELSE ! READ RECH FOR UNSTRUCTURED GRID
C3B-------INRECH=>0, SO READ RECHARGE RATE.
          IF(NPRCH.EQ.0) THEN
C
C3B1--------THERE ARE NO PARAMETERS, SO READ RECH USING U2DREL.
            CALL U2DREL(RECH,ANAME(2),1,NIRCH,0,IN,IOUT)
          ELSE
C
C3B2--------DEFINE RECH USING PARAMETERS.  INRECH IS THE NUMBER OF
C3B2--------PARAMETERS TO USE THIS STRESS PERIOD.
            CALL PRESET('RCH')
            WRITE(IOUT,33)
            IF(INRECH.EQ.0) THEN
              WRITE(IOUT,34)
              CALL USTOP(' ')
            END IF
            CALL UPARARRSUB2(RECH,NIRCH,1,0,INRECH,IN,IOUT,'RCH',
     1            ANAME(2),'RCH',IRCHPF)
          END IF
        ENDIF
      ENDIF
C
C5------IF NRCHOP=2 THEN A LAYER INDICATOR ARRAY IS NEEDED.  TEST INIRCH
C5------TO SEE HOW TO DEFINE IRCH.
        IF(NRCHOP.EQ.2) THEN
          IF(INIRCH.LT.0) THEN
C
C5A---------INIRCH<0, SO REUSE LAYER INDICATOR ARRAY FROM LAST STRESS PERIOD.
            WRITE(IOUT,2)
    2       FORMAT(1X,/1X,'REUSING IRCH FROM LAST STRESS PERIOD')
          ELSE
C
C5B---------INIRCH=>0, SO CALL U2DINT TO READ LAYER INDICATOR ARRAY(IRCH)
            IF(IUNSTR.EQ.0)THEN
              CALL U2DINT(ITEMP,ANAME(1),NROW,NCOL,0,IN,IOUT)
              N=0
              DO 57 IR=1,NROW
              DO 57 IC=1,NCOL
                N=N+1
                IF(ITEMP(IC,IR).LT.1 .OR. ITEMP(IC,IR).GT.NLAY) THEN
                  WRITE(IOUT,56) IC,IR,ITEMP(IC,IR)
   56             FORMAT(/1X,'INVALID LAYER NUMBER IN IRCH FOR COLUMN',
     1            I4,'  ROW',I4,'  :',I4)
                 CALL USTOP(' ')
                END IF
                IRCH(N) = (ITEMP(IC,IR)-1)*NROW*NCOL + (IR-1)*NCOL + IC
   57         CONTINUE
              NIRCH = NROW*NCOL
            ELSE
              CALL U2DINT(IRCH,ANAME(1),1,NIRCH,0,IN,IOUT)
C----------------------------------------------------
C ------------CHECK FOR IRCH BEING LARGER THAN NODES
              IFLAG = 0
              DO I=1,NIRCH
                IF(IRCH(I).GT.NODES)THEN
                  IFLAG = IRCH(I)
                  GO TO 112
                ENDIF
              ENDDO
112           CONTINUE
C ------------WRITE MESSAGE AND STOP IF IEVT IS LARGER THAN NODES
              IF(IFLAG.GT.0)THEN
                WRITE(IOUT,75)IFLAG,NODES
75              FORMAT('INDEX NODE NO.',I10,
     1          ', LARGER THAN TOTAL GWF NODES (',I10,'), STOPPING')
                STOP
              ENDIF
C----------------------------------------------------
            END IF
          END IF
        ELSE ! NRCHOP IS NOT 2 SO SET TOP LAYER OF NODES IN IRCH
          DO I=1,NIRCH
            IRCH(I) = I
          ENDDO
        END IF
C
C-------IF RECHARGE IS READ THEN MULTIPLY BY AREA TO GIVE FLUX
        IF(INRECH.GE.0) THEN
C
C4--------MULTIPLY RECHARGE RATE BY CELL AREA TO GET VOLUMETRIC RATE.
          DO 50 NN=1,NIRCH
            N = IRCH(NN)
            RECH(NN)=RECH(NN)*AREA(N)
   50     CONTINUE
        END IF
C----------------------------------------------------------------
C----------RECHARGE ZONES
      IF(IRTSOPT.EQ.0) GO TO 101
      IF(INiznrch.LE.0) THEN
C
C3A-----INiznrch=<0, SO REUSE iznrch ARRAY FROM LAST STRESS PERIOD.
        WRITE(IOUT,5)
    5   FORMAT(1X,/1X,'REUSING iznrch FROM LAST STRESS PERIOD')
      ELSEif(INiznrch.gt.0)then
C3B-----READ IZNRCH ARRAY AND FIRST TIME OF RTS FILE AT KPER=1
        mxznrch = iniznrch
        IF(IUNSTR.EQ.0)THEN
          CALL U2DINT(iznrch,ANAME(4),NROW,NCOL,0,IN,IOUT)
        ELSE
          CALL U2DINT(iznrch,ANAME(4),1,NIRCH,0,IN,IOUT)
        ENDIF
C3C-------READ FIRST LINE OF RTS FILE AT KPER=1
        IF(KPER.EQ.1)THEN
          inrts = IURTS
          allocate(tstartrch,tendrch,factrrch,rtsrch(mxznrch))
         read(inrts,*)tstartrch,tendrch,factrrch,(rtsrch(i),i=1,mxznrch)
         write(iout,7)tstartrch,tendrch,factrrch,(rtsrch(i),i=1,mxznrch)
7        format(2x,'*** RTS read - Tstart, Tend, Factor, Rts(mxznrch)'/
     1     5x,200g15.7)
C3D-------SET FLAGS FOR RTS AND ATS
          TIMRCH = TENDRCH
          IRTSRD = 0
        ENDIF
      ENDIF
C-----------------------------------------------------------------
C4--------APPLY RTS TO RECHARGE ARRAY IF RECHARGE OR ZONES CHANGE
        IF(INRECH.GE.0.OR.INIZNRCH.GT.0)THEN
c---------save original stress-period RECH in RECHSV array for later use
          IF(INRECH.GE.0)THEN
            DO NN=1,NIRCH
              RECHSV(NN) = RECH(NN)
            ENDDO
          ENDIF
c---------Add RTS recharge to RECH already on nodes
          DO 52 NN=1,NIRCH
          N = IRCH(NN)
          izr = iznrch(n)
          if(izr.ge.1.and.izr.le.mxznrch)
     *    RECH(NN)=RECHSV(NN) + rtsrch(izr)*AREA(N)*factrrch
   52     CONTINUE
          WRITE(IOUT,6)
6         FORMAT(2X,'*** RECH ARRAY UPDATED FROM RTS FILE ***')
        ENDIF
101   CONTINUE
C----------------------------------------------------------------
C----------UNCONFINED RECHARGE WITHOUT PONDING
      IF(IPONDOPT.EQ.0) GO TO 102
      IF(INSELEV.LE.0) THEN
C
C3A-----INSELEV<0, SO REUSE SELEV ARRAY FROM LAST STRESS PERIOD.
        WRITE(IOUT,4)
    4   FORMAT(1X,/1X,'REUSING SELEV FROM LAST STRESS PERIOD')
      ELSEif(INSELEV.gt.0)then
        IF(IUNSTR.EQ.0)THEN
          CALL U2DREL(SELEV,ANAME(3),NROW,NCOL,0,IN,IOUT)
        ELSE
          CALL U2DREL(SELEV,ANAME(3),1,NIRCH,0,IN,IOUT)
        ENDIF
      ENDIF
102   CONTINUE
C----------------------------------------------------------------
C----------CONCENTRATION OF RECHARGE FOR TRANSPORT
      IF(ICONCRCHOPT.EQ.0) GO TO 103
      IF(INCONC.LE.0) THEN
C
C3A-----INCONC<0, SO REUSE CONCENTRATION ARRAY FROM LAST STRESS PERIOD.
        WRITE(IOUT,8)
    8   FORMAT(1X,/1X,'REUSING CONCENTRATION FROM LAST STRESS PERIOD')
      ELSEif(INCONC.gt.0)then
        ICONCRCH = 0
        DO II=1,MCOMPT
!          WRITE(IOUT,*)(' READING FOR COMPONENT NUMBER',MCOMPT)
          WRITE(IOUT,*) ' READING FOR COMPONENT NUMBER',MCOMPT  !kkz - remove parentheses per JCH (alternative is to use a FORMAT)
          IF(IRCHCONC(II).NE.0)THEN
            ICONCRCH = ICONCRCH + 1
            IF(IUNSTR.EQ.0)THEN
             CALL U2DREL(RCHCONC(1,ICONCRCH),ANAME(5),NROW,NCOL,
     1             0,IN,IOUT)
            ELSE
             CALL U2DREL(RCHCONC(1,ICONCRCH),ANAME(5),1,NIRCH,0,IN,IOUT)
            ENDIF
          ENDIF
        ENDDO
      ENDIF
103   CONTINUE
C---------------------------------------------------------------
      DEALLOCATE(TEMP)
      DEALLOCATE(ITEMP)
C6------RETURN
      RETURN
      END
      SUBROUTINE GWF2RCH8U1FM(KPER)
C     ******************************************************************
C     SUBTRACT RECHARGE FROM RHS
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,      ONLY:NCOL,NROW,NLAY,IBOUND,RHS,IA,JA,JAS,NODLAY,
     1  IVC,AMAT,HNEW,ISSFLG
      USE GWFRCHMODULE,ONLY:NRCHOP,RECH,IRCH,NIRCH,SELEV,ISELEV
      DOUBLE PRECISION RECHFLUX,acoef,eps,pe,hd,rch,rcheps
C     ------------------------------------------------------------------
C
      ISS=ISSFLG(KPER)
      eps = 1.0e-5
C3------FILL RECH ON ACTIVE NODE.
      DO 10 NN=1,NIRCH
        N = IRCH(NN)
        RECHFLUX = RECH(NN)
C---------------------------------------------------------
C-------FIND TOP-MOST ACTIVE NODE IF NOT N
        IF(NRCHOP.EQ.3.AND.IBOUND(N).EQ.0)THEN
          CALL FIRST_ACTIVE_BELOW(N)
        ENDIF
C---------------------------------------------------------
C3A--------IF CELL IS VARIABLE HEAD, apply recharge as newton raphson
        IF(IBOUND(N).GT.0) then
          IF(ISELEV.GT.0)THEN
            pe = SELEV(nn)
            hd = hnew(n)
            call realrech(pe,rechflux,hd,rch)
            hd = hd + eps
            call realrech(pe,rechflux,hd,rcheps)
            acoef = (rcheps - rch) / eps
            amat(ia(n)) = amat(ia(n)) + acoef
            RHS(N)=RHS(N)-RCH + acoef * hnew(n)
          ELSE
            RHS(N)=RHS(N)-RECHFLUX
          ENDIF
        ENDIF
   10 CONTINUE
C
C6------RETURN
      RETURN
      END
C----------------------------------------------
      subroutine FIRST_ACTIVE_BELOW(N)
C     ******************************************************************
C     FIND FIRST ACTIVE NODE BELOW NODE N
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,      ONLY:IBOUND,IA,JA,JAS,NODLAY,IVC,IVSD,NODES
C----------------------------------------------
      IF(IVSD.EQ.-1)THEN
C-------VERTICALLY ADJACENT NODE BELOW IS N + NNDLAY
        NNDLAY = NODLAY(1)
3       CONTINUE
        JJ = N + NNDLAY
        IF(JJ.GT.NODES) RETURN
        IF(IBOUND(JJ).NE.0)THEN ! on const head node is ok (no effect)
          N = JJ
          GO TO 2
        ELSE
          N = JJ
          GO TO 3
        ENDIF
2       CONTINUE
      ELSE
C-------FIND VERTICALLY ADJACENT NODE BELOW
4       CONTINUE
        DO II = IA(N)+1,IA(N+1)-1
          JJ = JA(II)
          IIS = JAS(II)
          IF(IVC(IIS).EQ.1.AND.JJ.GT.N)THEN !VERTICAL DIRECTION DOWN
            IF(JJ.GT.NODES) RETURN
            IF(IBOUND(JJ).NE.0)THEN ! on const head node is ok (no effect)
              N = JJ
              GO TO 5
            ELSE
              N = JJ
              GO TO 4
            ENDIF
          ENDIF
        ENDDO
5       CONTINUE
      ENDIF
C
C------RETURN
      RETURN
      END
C----------------------------------------------
      subroutine realrech(pe,rechflux,hd,rch)
      DOUBLE PRECISION RECHFLUX,pe,hd,rch,depth,
     *  slope,epsilon
c----------------------------------------------
C------MAKE RCH GO TO ZERO OVER A HEAD INCREASE OF 0.01
      slope = -ABS(RECHFLUX/0.01)
      slope = - 1.0e4
      epsilon = 1.0e-3
      depth = hd - pe
      if(depth.lt.0)then
        rch = rechflux
cc      else
cc        rch = slope * depth**2 + rechflux
      elseif(depth.lt.epsilon)then
        rch = slope/(2.0*epsilon)*depth**2+rechflux
      else
        rch = rechflux + slope*(depth-epsilon/2.0)
      endif
C6------RETURN
      RETURN
      END
c----------------------------------------------
      SUBROUTINE GWF2RCH8U1BD(KSTP,KPER,INBCT)
C     ******************************************************************
C     CALCULATE VOLUMETRIC BUDGET FOR RECHARGE
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
C2------CLEAR THE RATE ACCUMULATORS.
      ZERO=0.
      RATIN=ZERO
      RATOUT=ZERO
C
C3------CLEAR THE BUFFER & SET FLAG FOR SAVING CELL-BY-CELL FLOW TERMS.
      DO 2 N=1,NODES
      BUFF(N)=ZERO
    2 CONTINUE
      IF(INBCT.GT.0)THEN
        DO N=1,NIRCH
          RCHF(N) = ZERO
        ENDDO
      ENDIF
      IBD=0
      IF(IRCHCB.GT.0) IBD=ICBCFL
      ALLOCATE(IBUFF(NIRCH))
C
C5------PROCESS EACH RECHARGE CELL LOCATION.
        DO 10 NN=1,NIRCH
        N = IRCH(NN)
        RECHFLUX = RECH(NN)
C---------------------------------------------------------
C-------FIND TOP-MOST ACTIVE NODE IF NOT N
        IF(NRCHOP.EQ.3.AND.IBOUND(N).EQ.0)THEN
          CALL FIRST_ACTIVE_BELOW(N)
        ENDIF
        IBUFF(NN) = N
C---------------------------------------------------------
C5A-----IF CELL IS VARIABLE HEAD, THEN DO BUDGET FOR IT.
        IF(IBOUND(N).GT.0) THEN
          IF(ISELEV.GT.0)THEN
            pe = SELEV(nn)
            hd = hnew(n)
            call realrech(pe,rechflux,hd,rch)
            QQ=rch
          ELSE
            QQ = RECHFLUX
          ENDIF
          Q=QQ
C
C5B-----ADD RECH TO BUFF.
          BUFF(N)=QQ
          IF(INBCT.GT.0) RCHF(NN) = Q
          FMBE(N) = FMBE(N) + QQ
C
C5C-----IF RECH POSITIVE ADD IT TO RATIN, ELSE ADD IT TO RATOUT.
          IF(Q.GE.ZERO) THEN
            RATIN=RATIN+QQ
          ELSE
            RATOUT=RATOUT-QQ
          END IF
        END IF
   10   CONTINUE
C
C8------IF CELL-BY-CELL FLOW TERMS SHOULD BE SAVED, CALL APPROPRIATE
C8------UTILITY MODULE TO WRITE THEM.
      IF(IUNSTR.EQ.0)THEN
        IF(IBD.EQ.1) CALL UBUDSV(KSTP,KPER,TEXT,IRCHCB,BUFF,NCOL,NROW,
     1                          NLAY,IOUT)
        IF(IBD.EQ.2) THEN
          ALLOCATE(ITEMP(NCOL,NROW))
          N=0
          DO I=1,NROW
            DO J=1,NCOL
              N=N+1
              ITEMP(J,I)= (IBUFF(N)-1) / (NCOL*NROW) + 1
            ENDDO
          ENDDO
          CALL UBDSV3(KSTP,KPER,TEXT,IRCHCB,BUFF,ITEMP,NRCHOP,
     1                NCOL,NROW,NLAY,IOUT,DELT,PERTIM,TOTIM,IBOUND)
          DEALLOCATE(ITEMP)
        ENDIF
      ELSE
        IF(IBD.EQ.1) CALL UBUDSVU(KSTP,KPER,TEXT,IRCHCB,BUFF,NODES,
     1                          IOUT,PERTIM,TOTIM)
        IF(IBD.EQ.2) CALL UBDSV3U(KSTP,KPER,TEXT,IRCHCB,BUFF,IBUFF,
     1           NIRCH,NRCHOP,NODES,IOUT,DELT,PERTIM,TOTIM,IBOUND)
      ENDIF
C
C9------MOVE TOTAL RECHARGE RATE INTO VBVL FOR PRINTING BY BAS1OT.
      ROUT=RATOUT
      RIN=RATIN
      VBVL(4,MSUM)=ROUT
      VBVL(3,MSUM)=RIN
C
C10-----ADD RECHARGE FOR TIME STEP TO RECHARGE ACCUMULATOR IN VBVL.
      VBVL(2,MSUM)=VBVL(2,MSUM)+ROUT*DELT
      VBVL(1,MSUM)=VBVL(1,MSUM)+RIN*DELT
C
C11-----MOVE BUDGET TERM LABELS TO VBNM FOR PRINT BY MODULE BAS_OT.
      VBNM(MSUM)=TEXT
C
C12-----INCREMENT BUDGET TERM COUNTER.
      MSUM=MSUM+1
c
      DEALLOCATE(IBUFF)
C
C13-----RETURN
      RETURN
      END
C----------------------------------------------------------------
      SUBROUTINE GWF2RCH8U1DA(INBCT)
C  Deallocate RCH DATA
      USE GWFRCHMODULE
      INTEGER ALLOC_ERR
C
        DEALLOCATE(MXZNRCH, STAT = ALLOC_ERR)
        DEALLOCATE(NRCHOP, STAT = ALLOC_ERR)
        DEALLOCATE(IRCHCB, STAT = ALLOC_ERR)
        DEALLOCATE(NPRCH, STAT = ALLOC_ERR)
        DEALLOCATE(IRCHPF, STAT = ALLOC_ERR)
        DEALLOCATE(RECH, STAT = ALLOC_ERR)
        DEALLOCATE(IRCH, STAT = ALLOC_ERR)
        DEALLOCATE(NRCHOP,IRCHCB,MXNDRCH, STAT = ALLOC_ERR)
        DEALLOCATE(NPRCH,IRCHPF,INIRCH,NIRCH, STAT = ALLOC_ERR)
        DEALLOCATE(INRTS,TIMRCH,IRTSRD, STAT = ALLOC_ERR)
        DEALLOCATE(IPONDOPT,IRTSOPT, STAT = ALLOC_ERR)
        DEALLOCATE(iznrch, STAT = ALLOC_ERR)
c        IF(INBCT.GT.0)THEN
        DEALLOCATE(RCHF, STAT = ALLOC_ERR)
c        ENDIF
        DEALLOCATE(ISELEV,SELEV,RECHSV, STAT = ALLOC_ERR)
        DEALLOCATE(IRCHCONC, STAT = ALLOC_ERR);
        DEALLOCATE(RCHCONC, STAT = ALLOC_ERR);
C
      RETURN
      END

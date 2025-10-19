C
C TIME-VARIANT MATERIALS (TVM) PACKAGE
C
C   Allows hydraulic conductivity and storage values to be changed
C   by linear interpolation throughout a transient simulation.
C
C Original author: D. Merrick, HydroAlgorithmics (May 2014).
C
      MODULE TVMMODULE
          INTEGER,SAVE :: ITVMPRINT,NTVMHK,NTVMVKA,NTVMSS,NTVMSY,
     1      NTVMDDFTR
          INTEGER,SAVE,ALLOCATABLE,DIMENSION(:) :: ITVMHK,ITVMVKA
          INTEGER,SAVE,ALLOCATABLE,DIMENSION(:) :: ITVMSS,ITVMSY
          INTEGER,SAVE,ALLOCATABLE,DIMENSION(:) :: ITVMDDFTR
          REAL,SAVE :: TVMLOGBASEHK,TVMLOGBASEHKINV
          REAL,SAVE :: TVMLOGBASEVKA,TVMLOGBASEVKAINV
          REAL,SAVE :: TVMLOGBASESS,TVMLOGBASESSINV
          REAL,SAVE :: TVMLOGBASESY,TVMLOGBASESYINV
          REAL,SAVE :: TVMDDFTR,TVMDDFTRINV
          REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: HKNEWA,VKANEWA
          REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: HKNEWB,VKANEWB
          REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: SSNEWA,SYNEWA
          REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: SSNEWB,SYNEWB
          REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: DDFTRNEWA,DDFTRNEWB
          DOUBLE PRECISION,SAVE,ALLOCATABLE,DIMENSION(:) :: INITPGF
      END MODULE
C
C
      SUBROUTINE TVMU1ARPRE
C     ******************************************************************
C     SAVE INITIAL GEOMETRIC FACTORS FROM PGF ARRAY PRIOR TO THEIR
C     MODIFICATION IN THE LPF PACKAGE
C     ******************************************************************
C
C     SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,ONLY:PGF,NJAS
      USE TVMMODULE
C-------ALLOCATE SPACE FOR SAVED PGF ARRAY
      ALLOCATE(INITPGF(NJAS))
C-------COPY ARRAY
      DO I = 1, NJAS
          INITPGF(I) = PGF(I)
      END DO
C-------RETURN.
      RETURN
      END
C
C
      SUBROUTINE TVMU1AR(IN,IUNITLPF,IUNITSFR,IUNITDPF)
C     ******************************************************************
C     ALLOCATE ARRAY STORAGE FOR TIME-VARIANT MATERIALS
C     ******************************************************************
C
C     SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,ONLY:IOUT,IFREFM,ITRSS,IUNSTR
      USE GWFBCFMODULE,ONLY:IKCFLAG
      USE GWFSFRMODULE,ONLY:ISFROPT
      USE TVMMODULE
      CHARACTER*200 :: LINE
C     ------------------------------------------------------------------
C-------IDENTIFY PACKAGE
      WRITE(IOUT,1) IN
    1 FORMAT(1X,/1X,'TVM -- TIME-VARIANT MATERIALS VERSION 1,',
     1  ' 6/11/2014',/1X,'INPUT READ FROM UNIT ',I4)
C-------ENSURE LPF PACKAGE IS USED
      IF(IUNITLPF.EQ.0) THEN
          WRITE (IOUT, 12) 
  12      FORMAT(1X,'LPF PACKAGE MUST BE ACTIVE TO CHANGE MATERIAL',
     1    1X,'PROPERTIES OF THE POROUS MEDIUM.')
      END IF
C-------ENSURE NODAL HYDRAULIC CONDUCTIVITIES ARE USED
      IF(IKCFLAG.NE.0) THEN
          WRITE (IOUT, *) 'TVM MAY ONLY BE USED WITH NODAL ',
     1        'CONDUCTIVITIES (IKCFLAG = 0 IN LPF PACKAGE).'
          CALL USTOP(' ')
      END IF
C-------WARN ABOUT USAGE IN STEADY-STATE SIMULATION
      IF(ITRSS.EQ.0) THEN
          WRITE (IOUT, *) 'WARNING: TVM PACKAGE IS ACTIVE BUT',
     1        ' SIMULATION IS STEADY-STATE!'
      END IF
C-------WARN ABOUT USAGE WITH STRUCTURED GRID
      IF(IUNSTR.EQ.0) THEN
          WRITE (IOUT, *) 'WARNING: TVM PACKAGE EXPECTS UNSTRUCTURED',
     1        ' CELL NUMBERS, BUT STRUCTURED GRID INPUT IS BEING USED.'
      END IF
C-------WARN ABOUT UNSATURATED SFR OPTIONS
      IF(IUNITSFR.GT.0.AND.ISFROPT.GT.1) THEN
          WRITE (IOUT, *) 'WARNING: SFR PACKAGE IS ACTIVE WITH ISFROPT',
     1        ' GREATER THAN 1. THIS IS UNSUPPORTED BY THE TVM',
     2        ' PACKAGE. UNSATURATED ZONE HYDRAULIC CONDUCTIVITIES AND',
     3        ' RESIDUAL WATER CONTENT VALUES WILL NOT BE UPDATED WITH',
     4        ' TVM MATERIAL PROPERTY CHANGES!'
      END IF
C-------READ AND PRINT COMMENTS, THEN READ FIRST LINE
      CALL URDCOM(IN,IOUT,LINE)
C-------READ PRINT FLAG
      IF(IFREFM.EQ.0) THEN
        READ(LINE,'(I10,5F10.0)') ITVMPRINT,TVMLOGBASEHK,TVMLOGBASEVKA,
     1                            TVMLOGBASESS,TVMLOGBASESY,TVMDDFTR
      ELSE
        READ(LINE,*) ITVMPRINT,TVMLOGBASEHK,TVMLOGBASEVKA,TVMLOGBASESS,
     1               TVMLOGBASESY,TVMDDFTR
      END IF
C-------REFLECT INPUT AND PRE-CALCULATE 1/LOG(TVMLOGBASE)
      !---TVMLOGBASEHK-------------------------------------------------
      IF(TVMLOGBASEHK.LT.0) THEN
          WRITE (IOUT, *)'STEP FUNCTION FOR HK IS USED, TVMLOGBASEHK = '
     1     ,TVMLOGBASEHK
      ELSE IF(TVMLOGBASEHK.EQ.0) THEN
          WRITE (IOUT, *) 'TVM INTERPOLATES HK VALUES LINEARLY',
     1        ' (TVMLOGBASEHK = 0)'
          TVMLOGBASEHKINV = 1.0
      ELSE
          WRITE (IOUT, *) 'TVM INTERPOLATES HK VALUES ',
     1        'LOGARITHMICALLY, TVMLOGBASEHK = ', TVMLOGBASEHK
          TVMLOGBASEHKINV = 1.0 / LOG(TVMLOGBASEHK)
      END IF
      !---TVMLOGBASEVKA------------------------------------------------
      IF(TVMLOGBASEVKA.LT.0) THEN
          WRITE (IOUT, *)'STEP FUNCTION FOR VKA IS USED,', 
     1     ' TVMLOGBASEVKA = ',TVMLOGBASEVKA
      ELSE IF(TVMLOGBASEVKA.EQ.0) THEN
          WRITE (IOUT, *) 'TVM INTERPOLATES VKA VALUES LINEARLY',
     1        ' (TVMLOGBASEVKA = 0)'
          TVMLOGBASEVKAINV = 1.0
      ELSE
          WRITE (IOUT, *) 'TVM INTERPOLATES VKA VALUES ',
     1        'LOGARITHMICALLY, TVMLOGBASEVKA = ', TVMLOGBASEVKA 
          TVMLOGBASEVKAINV = 1.0 / LOG(TVMLOGBASEVKA)
      END IF
      !---TVMLOGBASESS-------------------------------------------------
      IF(TVMLOGBASESS.LT.0) THEN
          WRITE (IOUT, *) 'STEP FUNCTION FOR SS IS USED,', 
     1     ' TVMLOGBASESS = ',TVMLOGBASESS
      ELSE IF(TVMLOGBASESS.EQ.0) THEN
          WRITE (IOUT, *) 'TVM INTERPOLATES SS VALUES LINEARLY',
     1        ' (TVMLOGBASESS = 0)'
          TVMLOGBASESSINV = 1.0
      ELSE
          WRITE (IOUT, *) 'TVM INTERPOLATES SS VALUES ',
     1        'LOGARITHMICALLY, TVMLOGBASESS = ', TVMLOGBASESS
          TVMLOGBASESSINV = 1.0 / LOG(TVMLOGBASESS)
      END IF
      !---TVMLOGBASESY-------------------------------------------------
      IF(TVMLOGBASESY.LT.0) THEN
          WRITE (IOUT, *) 'STEP FUNCTION FOR SY IS USED,', 
     1     ' TVMLOGBASESY = ',TVMLOGBASESY
      ELSE IF(TVMLOGBASESY.EQ.0) THEN
          WRITE (IOUT, *) 'TVM INTERPOLATES SY VALUES LINEARLY',
     1        ' (TVMLOGBASESY = 0)'
          TVMLOGBASESYINV = 1.0
      ELSE
          WRITE (IOUT, *) 'TVM INTERPOLATES SY VALUES ',
     1        'LOGARITHMICALLY, TVMLOGBASESY = ', TVMLOGBASESY
          TVMLOGBASESYINV = 1.0 / LOG(TVMLOGBASESY)
      END IF
      !---TVMDDFTR-------------------------------------------------
      IF(IUNITDPF.EQ.0) GO TO 100 
      IF(TVMDDFTR.LT.0) THEN
          WRITE (IOUT, *) 'STEP FUNCTION FOR DDFTR IS USED,', 
     1     ' TVMDDFTR = ',TVMDDFTR
      ELSE IF(TVMDDFTR.EQ.0) THEN
          WRITE (IOUT, *) 'TVM INTERPOLATES DDFTR VALUES LINEARLY',
     1        ' (TVMDDFTR = 0)'
          TVMDDFTRINV = 1.0
      ELSE
          WRITE (IOUT, *) 'TVM INTERPOLATES DDFTR VALUES ',
     1        'LOGARITHMICALLY, TVMDDFTR = ', TVMDDFTR
          TVMDDFTRINV = 1.0 / LOG(TVMDDFTR)
      END IF
100   CONTINUE   
C-------INIT VARIABLES
      NTVMHK = 0
      NTVMVKA = 0
      NTVMSS = 0
      NTVMSY = 0
      NTVMDDFTR = 0
C-------READ INITIAL CHANGED VALUES (START OF STRESS PERIOD 1)
      CALL TVMU1READ(IN,KPER)
C-------RETURN.
      RETURN
      END
C
C
      SUBROUTINE TVMU1READ(IN,KPER)
C     ******************************************************************
C     READ ONE SET OF VALUES FOR INTERPOLATION
C     ******************************************************************
C
C     SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,ONLY:IOUT,IFREFM,TOP,BOT,AREA
      USE GWFBCFMODULE,ONLY:SC1,SC2,VKA,HK,ISFAC
      USE TVMMODULE
      USE GWFDPFMODULE, ONLY: DDFTR
C     ------------------------------------------------------------------
C-------DEALLOCATE ANY EXISTING ARRAYS PRIOR TO READING NEW ARRAY VALUES
      CALL TVMU1DASP
C-------IDENTIFY PACKAGE
      WRITE(IOUT,*) 'TVM READING FROM UNIT ',IN
C-------READ NUMBER OF VALUES TO CHANGE FOR EACH MATERIAL PROPERTY
      IF(IFREFM.EQ.0) THEN
        READ(IN,'(5I10)') NTVMHK,NTVMVKA,NTVMSS,NTVMSY,NTVMDDFTR
      ELSE
        READ(IN,*) NTVMHK,NTVMVKA,NTVMSS,NTVMSY,NTVMDDFTR
      END IF
C-------WRITE SUMMARY OF THIS SET OF VALUES TO OUT FILE
      IF(ITVMPRINT.GE.1) THEN
c          WRITE(IOUT,*) 'READING ',NTVMHK,' HK, ',NTVMVKA,' VKA, ',
c     1        NTVMSS,' SS, ',NTVMSY,' SY AND',NTVMDDFTR,' DDFTR',
c     1        ' VALUES FOR STRESS PERIOD ', KPER
           WRITE(IOUT,1) NTVMHK,NTVMVKA,NTVMSS,NTVMSY,NTVMDDFTR,KPER
1          FORMAT(1X,'READING',I8,' HK, ',I8,' VKA, ',I8,' SS, ',I8,
     1     ' SY AND ',I8,' DDFTR VALUES FOR STRESS PERIOD ',I8)      
      END IF
C-------ALLOCATE ARRAYS
      IF(NTVMHK.GT.0) THEN
          ALLOCATE(ITVMHK(NTVMHK))
          ALLOCATE(HKNEWA(NTVMHK))
          ALLOCATE(HKNEWB(NTVMHK))
      END IF
      IF(NTVMVKA.GT.0) THEN
          ALLOCATE(ITVMVKA(NTVMVKA))
          ALLOCATE(VKANEWA(NTVMVKA))
          ALLOCATE(VKANEWB(NTVMVKA))
      END IF
      IF(NTVMSS.GT.0) THEN
          ALLOCATE(ITVMSS(NTVMSS))
          ALLOCATE(SSNEWA(NTVMSS))
          ALLOCATE(SSNEWB(NTVMSS))
      END IF
      IF(NTVMSY.GT.0) THEN
          ALLOCATE(ITVMSY(NTVMSY))
          ALLOCATE(SYNEWA(NTVMSY))
          ALLOCATE(SYNEWB(NTVMSY))
      END IF
      IF(NTVMDDFTR.GT.0) THEN
          ALLOCATE(ITVMDDFTR(NTVMDDFTR))
          ALLOCATE(DDFTRNEWA(NTVMDDFTR))
          ALLOCATE(DDFTRNEWB(NTVMDDFTR))
      END IF
C-------READ NEW HK VALUES
      DO I = 1, NTVMHK
          IF(IFREFM.EQ.0) THEN
            READ(IN,'(I10,F10.0)') ITVMHK(I), HKNEWB(I)
          ELSE
            READ(IN,*) ITVMHK(I), HKNEWB(I)
          END IF
          HKNEWA(I) = HK(ITVMHK(I))
          IF(ITVMPRINT.GE.2) THEN
              IF(KPER.GE.1) THEN
                IF(TVMLOGBASEHK.GE.0)THEN  
C                  WRITE(IOUT,*) '  HK FOR CELL ',ITVMHK(I),' WILL BE',
C     1            ' INTERPOLATED FROM ',HKNEWA(I),' TO ',HKNEWB(I)
                  WRITE(IOUT,2) ITVMHK(I),HKNEWA(I),HKNEWB(I)   
 2                FORMAT('  HK FOR CELL ',I8,' WILL BE INTERPOLATED',
     1            1X,'FROM', G15.6,' TO ',G15.6)       
                ELSE
C                  WRITE(IOUT,*) '  HK FOR CELL ',ITVMHK(I),' WILL BE',
C     1            ' SET TO ',HKNEWB(I) 
                  WRITE(IOUT,3) ITVMHK(I),HKNEWB(I)  
3                 FORMAT('  HK FOR CELL ',I8,' WILL BE SET TO ',G15.6)
                ENDIF    
              ELSE
                WRITE(IOUT,4) ITVMHK(I),HKNEWB(I)
4               FORMAT('  HK FOR CELL ',I8,' WILL BE CHANGED TO ',G15.6)
              END IF
          END IF
      END DO
C-------READ NEW VKA VALUES
      DO I = 1, NTVMVKA
          IF(IFREFM.EQ.0) THEN
            READ(IN,'(I10,F10.0)') ITVMVKA(I), VKANEWB(I)
          ELSE
            READ(IN,*) ITVMVKA(I), VKANEWB(I)
          END IF
          VKANEWA(I) = VKA(ITVMVKA(I))
          IF(ITVMPRINT.GE.2) THEN
              IF(KPER.GE.1) THEN
                IF(TVMLOGBASEVKA.GE.0)THEN  
C                  WRITE(IOUT,*) '  VKA FOR CELL ',ITVMVKA(I),' WILL BE',
C     1            ' INTERPOLATED FROM ',VKANEWA(I),' TO ',VKANEWB(I)
                  WRITE(IOUT,5)ITVMVKA(I),VKANEWA(I),VKANEWB(I)
 5                FORMAT(' VKA FOR CELL ',I8,' WILL BE INTERPOLATED',
     1            1X,'FROM', G15.6,' TO ',G15.6)                    
                ELSE
C                  WRITE(IOUT,*) '  VKA FOR CELL ',ITVMVKA(I),' WILL BE',
C     1            ' SET TO ',VKANEWB(I)  
                  WRITE(IOUT,6)ITVMVKA(I),VKANEWB(I)  
6                 FORMAT(' VKA FOR CELL ',I8,' WILL BE SET TO ',G15.6) 
                ENDIF    
              ELSE
C                WRITE(IOUT,*) '  VKA FOR CELL ',ITVMVKA(I),' WILL BE',
C     1            ' CHANGED TO ',VKANEWB(I)
                WRITE(IOUT,7) ITVMVKA(I),VKANEWB(I)   
7               FORMAT(' VKA FOR CELL ',I8,' WILL BE CHANGED TO ',G15.6)
              END IF
          END IF
      END DO
C-------READ NEW SS VALUES (SC1)
      DO I = 1, NTVMSS
          IF(IFREFM.EQ.0) THEN
            READ(IN,'(I10,F10.0)') ITVMSS(I), SSNEWB(I)
          ELSE
            READ(IN,*) ITVMSS(I), SSNEWB(I)
          END IF
          SSNEWB(I) = SSNEWB(I) * AREA(ITVMSS(I))
          IF(ISFAC.EQ.0) THEN
              THICK = TOP(ITVMSS(I)) - BOT(ITVMSS(I))
              SSNEWB(I) = SSNEWB(I) * THICK
          END IF
          SSNEWA(I) = SC1(ITVMSS(I))
          IF(ITVMPRINT.GE.2) THEN
              IF(KPER.GE.1) THEN
                IF(TVMLOGBASESS.GE.0)THEN
C                  WRITE(IOUT,*) '  SC1 FOR CELL ',ITVMSS(I),' WILL BE',
C     1            ' INTERPOLATED FROM ',SSNEWA(I),' TO ',SSNEWB(I)
                  WRITE(IOUT,8) ITVMSS(I),SSNEWA(I),SSNEWB(I) 
 8                FORMAT(' SC1 FOR CELL ',I8,' WILL BE INTERPOLATED',
     1            1X,'FROM', G15.6,' TO ',G15.6) 
                ELSE
C                  WRITE(IOUT,*) '  SC1 FOR CELL ',ITVMSS(I),' WILL BE',
C     1            ' SET TO ',SSNEWB(I)  
                  WRITE(IOUT,9) ITVMSS(I),SSNEWB(I)               
9                 FORMAT(' SC1 FOR CELL ',I8,' WILL BE SET TO ',G15.6) 
                ENDIF
              ELSE
C                WRITE(IOUT,*) '  SC1 FOR CELL ',ITVMVKA(I),' WILL BE',
C     1            ' CHANGED TO ',SSNEWB(I)
                WRITE(IOUT,10)ITVMVKA(I),SSNEWB(I)
10              FORMAT(' SC1 FOR CELL ',I8,' WILL BE CHANGED TO ',G15.6)
              END IF
          END IF
      END DO
C-------READ NEW SY VALUES (SC2)
      DO I = 1, NTVMSY
          IF(IFREFM.EQ.0) THEN
            READ(IN,'(I10,F10.0)') ITVMSY(I), SYNEWB(I)
          ELSE
            READ(IN,*) ITVMSY(I), SYNEWB(I)
          END IF
          SYNEWB(I) = SYNEWB(I) * AREA(ITVMSY(I))
          SYNEWA(I) = SC2(ITVMSY(I))
          IF(ITVMPRINT.GE.2) THEN
              IF(KPER.GE.1) THEN
                IF(TVMLOGBASESY.GE.0)THEN  
C                  WRITE(IOUT,*) '  SC2 FOR CELL ',ITVMSY(I),' WILL BE',
C     1            ' INTERPOLATED FROM ',SYNEWA(I),' TO ',SYNEWB(I)
                  WRITE(IOUT,11) ITVMSY(I),SYNEWA(I),SYNEWB(I)
 11                FORMAT(' SC2 FOR CELL ',I8,' WILL BE INTERPOLATED',
     1            1X,'FROM', G15.6,' TO ',G15.6) 
                ELSE
C                  WRITE(IOUT,*) '  SC2 FOR CELL ',ITVMSY(I),' WILL BE',
C     1            ' SET TO ',SYNEWB(I)  
                  WRITE(IOUT,12) ITVMSY(I),SYNEWB(I)
12                FORMAT(' SC2 FOR CELL ',I8,' WILL BE SET TO ',G15.6)
                ENDIF
              ELSE
C                WRITE(IOUT,*) '  SC2 FOR CELL ',ITVMSY(I),' WILL BE',
C     1            ' CHANGED TO ',SYNEWB(I)
                WRITE(IOUT,13) ITVMSY(I),SYNEWB(I)
13              FORMAT(' SC2 FOR CELL ',I8,' WILL BE CHANGED TO ',G15.6)
              END IF
          END IF
      END DO
C-------READ NEW DDFTR VALUES 
      DO I = 1, NTVMDDFTR
          IF(IFREFM.EQ.0) THEN
            READ(IN,'(I10,F10.0)') ITVMDDFTR(I), DDFTRNEWB(I)
          ELSE
            READ(IN,*) ITVMDDFTR(I), DDFTRNEWB(I)
          END IF
          THICK = TOP(ITVMDDFTR(I)) - BOT(ITVMDDFTR(I))
          DDFTRNEWB(I) = DDFTRNEWB(I) * AREA(ITVMDDFTR(I)) * THICK
          DDFTRNEWA(I) = DDFTR(ITVMDDFTR(I))
          IF(ITVMPRINT.GE.2) THEN
              IF(KPER.GE.1) THEN
                IF(TVMDDFTR.GE.0)THEN
C                 WRITE(IOUT,*) '  DDFTR FOR CELL ',ITVMDDFTR(I),' WILL',
C     1         ' BE INTERPOLATED FROM ',DDFTRNEWA(I),' TO ',DDFTRNEWB(I)
                 WRITE(IOUT,14) ITVMDDFTR(I),DDFTRNEWA(I),DDFTRNEWB(I)
14                FORMAT(' DDFTR FOR CELL ',I8,' WILL BE INTERPOLATED',
     1            1X,'FROM', G15.6,' TO ',G15.6) 
                ELSE
C                 WRITE(IOUT,*) '  DDFTR FOR CELL ',ITVMDDFTR(I),' WILL',
C     1            ' BE SET TO ',DDFTRNEWB(I)  
                 WRITE(IOUT,15) ITVMDDFTR(I),DDFTRNEWB(I) 
15               FORMAT(' DDFTR FOR CELL ',I8,' WILL BE SET TO ',G15.6)
                ENDIF
              ELSE
C                WRITE(IOUT,*) '  DDFTR FOR CELL ',ITVMDDFTR(I),' WILL',
C     1            ' BE CHANGED TO ',DDFTRNEWB(I)
                WRITE(IOUT,16)ITVMDDFTR(I),DDFTRNEWB(I)  
16            FORMAT(' DDFTR FOR CELL ',I8,' WILL BE CHANGED TO ',G15.6)
              END IF
          END IF
      END DO      
C-------RETURN
      RETURN
      END
C
C
      SUBROUTINE TVMU1RP(IN,IUNITCLN,IUNITSFR,IUNITDPF,KPER)
C     ******************************************************************
C     READ MATERIAL PROPERTY VALUES FOR CURRENT STRESS PERIOD
C     ******************************************************************
      USE GLOBAL,ONLY:IOUT,NLAY,NJAS,PGF
      USE GWFBCFMODULE,ONLY:SC1,SC2,VKA,HK
      USE TVMMODULE
      USE GWFDPFMODULE, ONLY: DDFTR
C-------APPLY FINAL VALUES FROM PREVIOUS STRESS PERIOD, IF ANY
C-------CHANGE HK VALUES
      IF(NTVMHK.GT.0) THEN
          IF(ITVMPRINT.GE.1) THEN
            IF(KPER.GT.1) THEN
C              WRITE(IOUT,*) 'TVM CHANGING ',NTVMHK,' HK VALUES AT ',
C     1          ' END OF STRESS PERIOD ',(KPER - 1)
              WRITE(IOUT,1)NTVMHK,(KPER - 1) 
1             FORMAT(' TVM CHANGING ',I8,' HK VALUES AT END OF STRESS',
     1        1X,'PERIOD ',I8)
            ELSE
C              WRITE(IOUT,*) 'TVM CHANGING ',NTVMHK,' HK VALUES AT ',
C     1          ' START OF STRESS PERIOD 1'
              WRITE(IOUT,2)NTVMHK   
2             FORMAT(' TVM CHANGING ',I8,' HK VALUES AT START OF STRESS'
     1        ,1X,'PERIOD  1')                 
            END IF
          END IF
          DO I = 1, NTVMHK
              HK(ITVMHK(I)) = HKNEWB(I)
              IF(ITVMPRINT.GE.2) THEN
                IF(KPER.GT.1) THEN
C                  WRITE(IOUT,*) '  HK FOR CELL ',ITVMHK(I),' CHANGED',
C     1              ' TO ',HK(ITVMHK(I)),' AT END OF STRESS PERIOD ',
C     2              (KPER - 1)
                  WRITE(IOUT,3) ITVMHK(I),HK(ITVMHK(I)),(KPER - 1)  
3                FORMAT('  HK FOR CELL ',I8,' CHANGED TO ',G15.6,
     1           ' AT END OF STRESS PERIOD ',I8) 
                ELSE
C                  WRITE(IOUT,*) '  HK FOR CELL ',ITVMHK(I),' CHANGED',
C     1             ' TO ',HK(ITVMHK(I)),' AT START OF STRESS PERIOD 1'
                  WRITE(IOUT,4)ITVMHK(I),HK(ITVMHK(I))               
4                 FORMAT('  HK FOR CELL ',I8,' CHANGED TO ',G15.6,
     1            ' AT START OF STRESS PERIOD 1')
                END IF
              END IF
          END DO
      END IF
C-------CHANGE VKA VALUES
      IF(NTVMVKA.GT.0) THEN
          IF(ITVMPRINT.GE.1) THEN
            IF(KPER.GT.1) THEN
C              WRITE(IOUT,*) 'TVM CHANGING ',NTVMVKA,' VKA VALUES AT ',
C     1          ' END OF STRESS PERIOD ',(KPER - 1)
              WRITE(IOUT,11)NTVMVKA,(KPER - 1)              
11            FORMAT(' TVM CHANGING ',I8,' VKA VALUES AT END OF STRESS',
     1        1X,'PERIOD ',I8)              
            ELSE
C              WRITE(IOUT,*) 'TVM CHANGING ',NTVMVKA,' VKA VALUES AT ',
C     1          ' START OF STRESS PERIOD 1'
              WRITE(IOUT,12) NTVMVKA      
12            FORMAT(' TVM CHANGING ',I8,' VKA VALUES AT START OF'
     1        ,1X,'STRESS PERIOD  1')                
            END IF
          END IF
          DO I = 1, NTVMVKA
              VKA(ITVMVKA(I)) = VKANEWB(I)
              IF(ITVMPRINT.GE.2) THEN
                IF(KPER.GT.1) THEN
C                  WRITE(IOUT,*) '  VKA FOR CELL ',ITVMVKA(I),' CHANGED',
C     1              ' TO ',VKA(ITVMVKA(I)),' AT END OF STRESS PERIOD ',
C     2              (KPER - 1)
                  WRITE(IOUT,13) ITVMVKA(I),VKA(ITVMVKA(I)),(KPER - 1)
13                FORMAT('  VKA FOR CELL ',I8,' CHANGED TO ',G15.6,
     1           ' AT END OF STRESS PERIOD ',I8) 
                ELSE
C                  WRITE(IOUT,*) '  VKA FOR CELL ',ITVMVKA(I),' CHANGED',
C     1             ' TO ',VKA(ITVMVKA(I)),' AT START OF STRESS PERIOD 1'
                  WRITE(IOUT,14) ITVMVKA(I),VKA(ITVMVKA(I))            
14                FORMAT('  VKA FOR CELL ',I8,' CHANGED TO ',G15.6,
     1            ' AT START OF STRESS PERIOD 1')                  
                END IF
              END IF
          END DO
      END IF
C-------CHANGE SS VALUES (SC1)
      IF(NTVMSS.GT.0) THEN
          IF(ITVMPRINT.GE.1) THEN
            IF(KPER.GT.1) THEN
C              WRITE(IOUT,*) 'TVM CHANGING ',NTVMSS,' SC1 VALUES AT ',
C     1            ' END OF STRESS PERIOD ',(KPER - 1)
               WRITE(IOUT,21) NTVMSS,(KPER - 1)       
21            FORMAT(' TVM CHANGING ',I8,' SC1 VALUES AT END OF STRESS',
     1        1X,'PERIOD ',I8)                
            ELSE
C              WRITE(IOUT,*) 'TVM CHANGING ',NTVMSS,' SC1 VALUES AT ',
C     1            ' START OF STRESS PERIOD 1'
              WRITE(IOUT,22) NTVMSS    
22            FORMAT(' TVM CHANGING ',I8,' SC1 VALUES AT START OF'
     1        ,1X,'STRESS PERIOD  1')               
            END IF
          END IF
          DO I = 1, NTVMSS
              SC1(ITVMSS(I)) = SSNEWB(I)
              IF(ITVMPRINT.GE.2) THEN
                IF(KPER.GT.1) THEN
C                  WRITE(IOUT,*) '  SC1 FOR CELL ',ITVMSS(I),' CHANGED',
C     1              ' TO ',SC1(ITVMSS(I)),' AT END OF STRESS PERIOD ',
C     2              (KPER - 1)
                  WRITE(IOUT,23) ITVMSS(I),SC1(ITVMSS(I)),(KPER - 1)  
23                FORMAT('  SC1 FOR CELL ',I8,' CHANGED TO ',G15.6,
     1           ' AT END OF STRESS PERIOD ',I8)                   
                ELSE
C                  WRITE(IOUT,*) '  SC1 FOR CELL ',ITVMSS(I),' CHANGED',
C     1              ' TO ',SC1(ITVMSS(I)),' AT START OF STRESS PERIOD 1'
                  WRITE(IOUT,24) ITVMSS(I),SC1(ITVMSS(I))           
24                FORMAT('  SC1 FOR CELL ',I8,' CHANGED TO ',G15.6,
     1            ' AT START OF STRESS PERIOD 1')                  
                END IF
              END IF
          END DO
      END IF
C-------CHANGE SY VALUES (SC2)
      IF(NTVMSY.GT.0) THEN
          IF(ITVMPRINT.GE.1) THEN
            IF(KPER.GT.1) THEN
C              WRITE(IOUT,*) 'TVM CHANGING ',NTVMSY,' SC2 VALUES AT ',
C     1            ' END OF STRESS PERIOD ',(KPER - 1)
              WRITE(IOUT,31)NTVMSY,(KPER - 1)         
31            FORMAT(' TVM CHANGING ',I8,' SC2 VALUES AT END OF STRESS',
     1        1X,'PERIOD ',I8)              
            ELSE
C              WRITE(IOUT,*) 'TVM CHANGING ',NTVMSY,' SC2 VALUES AT ',
C     1            ' START OF STRESS PERIOD 1'
              WRITE(IOUT,32) NTVMSY   
32            FORMAT(' TVM CHANGING ',I8,' SC2 VALUES AT START OF'
     1        ,1X,'STRESS PERIOD  1')              
            END IF
          END IF
          DO I = 1, NTVMSY
              SC2(ITVMSY(I)) = SYNEWB(I)
              IF(ITVMPRINT.GE.2) THEN
                IF(KPER.GT.1) THEN
C                  WRITE(IOUT,*) '  SC2 FOR CELL ',ITVMSY(I),' CHANGED',
C     1              ' TO ',SC2(ITVMSY(I)),' AT END OF STRESS PERIOD ',
C     2              (KPER - 1)
                  WRITE(IOUT,33)ITVMSY(I),SC2(ITVMSY(I)),(KPER - 1)
33                FORMAT('  SC2 FOR CELL ',I8,' CHANGED TO ',G15.6,
     1           ' AT END OF STRESS PERIOD ',I8)                    
                ELSE
C                  WRITE(IOUT,*) '  SC2 FOR CELL ',ITVMSY(I),' CHANGED',
C     1              ' TO ',SC2(ITVMSY(I)),' AT START OF STRESS PERIOD 1'
                  WRITE(IOUT,34)ITVMSY(I),SC2(ITVMSY(I))  
34                FORMAT('  SC2 FOR CELL ',I8,' CHANGED TO ',G15.6,
     1            ' AT START OF STRESS PERIOD 1')                    
                END IF
              END IF
          END DO
      END IF
C-------CHANGE DDFTR VALUES 
      IF(IUNITDPF.GT.0 .AND. NTVMDDFTR.GT.0) THEN
          IF(ITVMPRINT.GE.1) THEN
            IF(KPER.GT.1) THEN
C              WRITE(IOUT,*) 'TVM CHANGING ',NTVMDDFTR,' DDFTR VALUES ',
C     1            ' AT END OF STRESS PERIOD ',(KPER - 1)
              WRITE(IOUT,41) NTVMDDFTR,(KPER - 1)      
41            FORMAT(' TVM CHANGING ',I8,' DDFTR VALUES AT END OF',
     1        1X,'STRESS PERIOD ',I8)               
            ELSE
C              WRITE(IOUT,*) 'TVM CHANGING ',NTVMDDFTR,' DDFTR VALUES ',
C     1            ' AT START OF STRESS PERIOD 1'
              WRITE(IOUT,42)NTVMDDFTR       
42            FORMAT(' TVM CHANGING ',I8,' DDFTR VALUES AT START OF'
     1        ,1X,'STRESS PERIOD  1')                 
            END IF
          END IF
          DO I = 1, NTVMDDFTR
              DDFTR(ITVMSY(I)) = DDFTRNEWB(I)
              IF(ITVMPRINT.GE.2) THEN
                IF(KPER.GT.1) THEN
C                  WRITE(IOUT,*) '  DDFTR FOR CELL ',ITVMDDFTR(I),
C     1              ' CHANGED TO ',DDFTR(ITVMDDFTR(I)),
C     2             ' AT END OF STRESS PERIOD ', (KPER - 1)
                  WRITE(IOUT,43) ITVMDDFTR(I),DDFTR(ITVMDDFTR(I)),
     1              (KPER - 1)      
43                FORMAT('  DDFTR FOR CELL ',I8,' CHANGED TO ',G15.6,
     1           ' AT END OF STRESS PERIOD ',I8)                   
                ELSE
C                  WRITE(IOUT,*) '  DDFTR FOR CELL ',ITVMDDFTR(I),
C     1              ' CHANGED TO ',DDFTR(ITVMDDFTR(I)),
C     2              ' AT START OF STRESS PERIOD 1'
                  WRITE(IOUT,44)ITVMDDFTR(I),DDFTR(ITVMDDFTR(I))
44                FORMAT('  DDFTR FOR CELL ',I8,' CHANGED TO ',G15.6,
     1            ' AT START OF STRESS PERIOD 1')                   
                END IF
              END IF
          END DO
      END IF      
C-------RECOMPUTE PGF ARRAY VALUES IF ANY HK OR VKA VALUES HAVE CHANGED
      IF(NTVMHK.GT.0.OR.NTVMVKA.GT.0) THEN
C----------REPLACE PGF VALUES WITH INITIAL GEOMETRIC FACTORS
          DO I = 1, NJAS
              PGF(I) = INITPGF(I)
          END DO
C----------LPF UPDATE (ASSUMED: LPF PACKAGE IS ACTIVE)
          CALL SGWF2LPFU1N
          IF(NLAY.GT.1) THEN
              CALL SGWF2LPFU1VCOND
          END IF
          CALL FILLPGFH
C----------CLN UPDATE
          IF(IUNITCLN.GT.0) THEN
              CALL SFILLPGF_CLN
          END IF
      END IF
C-------READ NEXT SET OF CHANGED VALUES
      CALL TVMU1READ(IN,KPER)
C-------RETURN
      RETURN
      END
C
C
      SUBROUTINE TVMU1AD(IN,IUNITCLN,IUNITSFR,IUNITDPF,KPER)
C     ******************************************************************
C     SET MATERIAL PROPERTY VALUES AT CURRENT TIME STEP
C     ******************************************************************
      USE GLOBAL,ONLY:IOUT,NLAY,NJAS,PGF
      USE GWFBASMODULE,ONLY:TOTIM
      USE GWFBCFMODULE,ONLY:SC1,SC2,VKA,HK,ISFAC
      USE TVMMODULE
      USE GWFDPFMODULE, ONLY: DDFTR      
C-------NOTHING TO CHANGE IF ALL COUNTS ARE ZERO
      IF(NTVMHK.LE.0.AND.NTVMVKA.LE.0.AND.NTVMSS.LE.0.AND.
     1   NTVMSY.LE.0.AND.NTVMDDFTR.LE.0) THEN
          IF(ITVMPRINT.GE.1) THEN
              WRITE(IOUT,*) 'TVM LEAVING ALL MATERIAL VALUES UNCHANGED',
     1            ' AT TIME ',TOTIM,' AS ALL COUNTS ARE ZERO'
          END IF
          RETURN
      END IF
C-------CHANGE HK VALUES IF NOT STEP FUNCTION
      IF(NTVMHK.GT.0.AND.TVMLOGBASEHK.GT.-0.1) THEN
          IF(ITVMPRINT.GE.1) THEN
C              WRITE(IOUT,*) 'TVM CHANGING ',NTVMHK,' HK VALUES AT ',
C     1            ' TIME ',TOTIM
              WRITE(IOUT,1)NTVMHK,TOTIM
1             FORMAT(' TVM CHANGING ',I8,' HK VALUES AT TIME ',G15.6)  
          END IF
          DO I = 1, NTVMHK
              HK(ITVMHK(I)) = TVMU1INTERP(HKNEWA(I),HKNEWB(I),KPER,
     1                            TVMLOGBASEHK,TVMLOGBASEHKINV)
              IF(ITVMPRINT.GE.2) THEN
C                  WRITE(IOUT,*) '  HK FOR CELL ',ITVMHK(I),' CHANGED',
C     1                ' TO ',HK(ITVMHK(I)),' AT TIME ',TOTIM
                  WRITE(IOUT,2) ITVMHK(I),HK(ITVMHK(I)),TOTIM
2                 FORMAT('  HK FOR CELL ',I8,' CHANGED TO ',G15.6,
     1             ' AT TIME ',G15.6)                   
              END IF
          END DO
      END IF
C-------CHANGE VKA VALUES
      IF(NTVMVKA.GT.0.AND.TVMLOGBASEVKA.GT.-0.1) THEN
          IF(ITVMPRINT.GE.1) THEN
C              WRITE(IOUT,*) 'TVM CHANGING ',NTVMVKA,' VKA VALUES AT ',
C     1            ' TIME ',TOTIM
              WRITE(IOUT,11)NTVMVKA,TOTIM  
11            FORMAT(' TVM CHANGING ',I8,' VKA VALUES AT TIME ',G15.6)
          END IF
          DO I = 1, NTVMVKA
              VKA(ITVMVKA(I)) = TVMU1INTERP(VKANEWA(I),VKANEWB(I),KPER,
     1                              TVMLOGBASEVKA,TVMLOGBASEVKAINV)
              IF(ITVMPRINT.GE.2) THEN
C                  WRITE(IOUT,*) '  VKA FOR CELL ',ITVMVKA(I),' CHANGED',
C     1                ' TO ',VKA(ITVMVKA(I)),' AT TIME ',TOTIM
                  WRITE(IOUT,12) ITVMVKA(I),VKA(ITVMVKA(I)),TOTIM       
12                FORMAT('  VKA FOR CELL ',I8,' CHANGED TO ',G15.6,
     1             ' AT TIME ',G15.6)                   
              END IF
          END DO
      END IF
C-------CHANGE SS VALUES (SC1)
      IF(NTVMSS.GT.0.AND.TVMLOGBASESS.GT.-0.1) THEN
          IF(ITVMPRINT.GE.1) THEN
C              WRITE(IOUT,*) 'TVM CHANGING ',NTVMVKA,' SC1 VALUES AT ',
C     1            ' TIME ',TOTIM
              WRITE(IOUT,21)NTVMVKA,TOTIM
 21            FORMAT(' TVM CHANGING ',I8,' SC1 VALUES AT TIME ',G15.6)
          END IF
          DO I = 1, NTVMSS
              SC1(ITVMSS(I)) = TVMU1INTERP(SSNEWA(I),SSNEWB(I),KPER,
     1                             TVMLOGBASESS,TVMLOGBASESSINV)
              IF(ITVMPRINT.GE.2) THEN
C                  WRITE(IOUT,*) '  SC1 FOR CELL ',ITVMSS(I),' CHANGED',
C     1                ' TO ',SC1(ITVMSS(I)),' AT TIME ',TOTIM
                  WRITE(IOUT,22)ITVMSS(I),SC1(ITVMSS(I)),TOTIM
22                FORMAT('  SC1 FOR CELL ',I8,' CHANGED TO ',G15.6,
     1             ' AT TIME ',G15.6)                  
              END IF
          END DO
      END IF
C-------CHANGE SY VALUES (SC2)
      IF(NTVMSY.GT.0.AND.TVMLOGBASESY.GT.-0.1) THEN
          IF(ITVMPRINT.GE.1) THEN
C              WRITE(IOUT,*) 'TVM CHANGING ',NTVMVKA,' SC2 VALUES AT ',
C     1            ' TIME ',TOTIM
              WRITE(IOUT,31) NTVMVKA,TOTIM    
31            FORMAT(' TVM CHANGING ',I8,' SC2 VALUES AT TIME ',G15.6)
          END IF
          DO I = 1, NTVMSY
              SC2(ITVMSY(I)) = TVMU1INTERP(SYNEWA(I),SYNEWB(I),KPER,
     1                             TVMLOGBASESY,TVMLOGBASESYINV)
              IF(ITVMPRINT.GE.2) THEN
C                  WRITE(IOUT,*) '  SC2 FOR CELL ',ITVMSY(I),' CHANGED',
C     1                ' TO ',SC2(ITVMSY(I)),' AT TIME ',TOTIM
                  WRITE(IOUT,32)ITVMSY(I),SC2(ITVMSY(I)),TOTIM
32                FORMAT('  SC2 FOR CELL ',I8,' CHANGED TO ',G15.6,
     1             ' AT TIME ',G15.6)                   
              END IF
          END DO
      END IF
C-------CHANGE DDFTR VALUES 
      IF(NTVMDDFTR.GT.0.AND.TVMDDFTR.GT.-0.1) THEN
          IF(ITVMPRINT.GE.1) THEN
C              WRITE(IOUT,*) 'TVM CHANGING ',NTVMDDFTR,' DDFTR VALUES ',
C     1            ' AT TIME ',TOTIM
              WRITE(IOUT,41) NTVMDDFTR,TOTIM
41            FORMAT(' TVM CHANGING ',I8,' DDFTR VALUES AT TIME ',G15.6)
          END IF
          DO I = 1, NTVMDDFTR
              DDFTR(ITVMDDFTR(I)) = TVMU1INTERP(DDFTRNEWA(I),
     1             DDFTRNEWB(I),KPER,TVMDDFTR,TVMDDFTRINV)
              IF(ITVMPRINT.GE.2) THEN
C                  WRITE(IOUT,*) '  DDFTR FOR CELL ',ITVMDDFTR(I),
C     1              ' CHANGED TO ',DDFTR(ITVMDDFTR(I)),' AT TIME ',TOTIM
                  WRITE(IOUT,42) ITVMDDFTR(I),DDFTR(ITVMDDFTR(I)),TOTIM
42                FORMAT('  DDFTR FOR CELL ',I8,' CHANGED TO ',G15.6,
     1             ' AT TIME ',G15.6)                   
              END IF
          END DO
      END IF      
C-------RECOMPUTE PGF ARRAY VALUES IF ANY HK OR VKA VALUES HAVE CHANGED
      IF(NTVMHK.GT.0.OR.NTVMVKA.GT.0) THEN
C----------REPLACE PGF VALUES WITH INITIAL GEOMETRIC FACTORS
          DO I = 1, NJAS
              PGF(I) = INITPGF(I)
          END DO
C----------LPF UPDATE (ASSUMED: LPF PACKAGE IS ACTIVE)
          CALL SGWF2LPFU1N
          IF(NLAY.GT.1) THEN
              CALL SGWF2LPFU1VCOND
          END IF
          CALL FILLPGFH
C----------CLN UPDATE
          IF(IUNITCLN.GT.0) THEN
              CALL SFILLPGF_CLN
          END IF
      END IF
C
      RETURN
      END
C
C
      REAL FUNCTION TVMU1INTERP(STARTVALUE, ENDVALUE, KPER, TVMLOGBASE,
     1                          TVMLOGBASEINV)
C     ******************************************************************
C     INTERPOLATE BETWEEN STARTVALUE AND ENDVALUE AT PERTIM IN KPER
C     ******************************************************************
      USE GWFBASMODULE,ONLY:PERTIM
      USE GLOBAL,ONLY:PERLEN
      USE TVMMODULE
C
      REAL :: VALUEA, VALUEB
C
      IF(PERLEN(KPER).LE.0) THEN
          TVMU1INTERP = ENDVALUE
      ELSE
          FRAC = PERTIM / PERLEN(KPER)

          IF(TVMLOGBASE.LE.0) THEN
              VALUEA = STARTVALUE
              VALUEB = ENDVALUE
          ELSE IF(TVMLOGBASE.EQ.10.0) THEN
              VALUEA = LOG10(STARTVALUE)
              VALUEB = LOG10(ENDVALUE)
          ELSE
              VALUEA = LOG(STARTVALUE)*TVMLOGBASEINV
              VALUEB = LOG(ENDVALUE)*TVMLOGBASEINV
          END IF

          TVMU1INTERP = VALUEA + FRAC * (VALUEB - VALUEA)

          IF(TVMLOGBASE.GT.0) THEN
              TVMU1INTERP = TVMLOGBASE**TVMU1INTERP
          END IF
      END IF
C
      RETURN
      END
C
C
      SUBROUTINE TVMU1DASP
C     ******************************************************************
C     DEALLOCATE TVM STRESS PERIOD ARRAYS
C     ******************************************************************
      USE TVMMODULE
C
      IF(NTVMHK.GT.0) THEN
          DEALLOCATE(ITVMHK)
          DEALLOCATE(HKNEWA)
          DEALLOCATE(HKNEWB)
          NTVMHK = 0
      END IF
      IF(NTVMVKA.GT.0) THEN
          DEALLOCATE(ITVMVKA)
          DEALLOCATE(VKANEWA)
          DEALLOCATE(VKANEWB)
          NTVMVKA = 0
      END IF
      IF(NTVMSS.GT.0) THEN
          DEALLOCATE(ITVMSS)
          DEALLOCATE(SSNEWA)
          DEALLOCATE(SSNEWB)
          NTVMSS = 0
      END IF
      IF(NTVMSY.GT.0) THEN
          DEALLOCATE(ITVMSY)
          DEALLOCATE(SYNEWA)
          DEALLOCATE(SYNEWB)
          NTVMSY = 0
      END IF
      IF(NTVMDDFTR.GT.0) THEN
          DEALLOCATE(ITVMDDFTR)
          DEALLOCATE(DDFTRNEWA)
          DEALLOCATE(DDFTRNEWB)
          NTVMDDFTR = 0
      END IF
C
      RETURN
      END
C
C
      SUBROUTINE TVMU1DA
C     ******************************************************************
C     DEALLOCATE TVM ARRAYS
C     ******************************************************************
      USE TVMMODULE
C
      CALL TVMU1DASP
C
      DEALLOCATE(INITPGF)
C
      RETURN
      END

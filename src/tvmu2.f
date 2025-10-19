C
C TIME-VARIANT MATERIALS (TVM) PACKAGE
C
C   Allows hydraulic conductivity and storage values to be changed
C   by linear interpolation throughout a transient simulation.
C
C Original author: D. Merrick, HydroAlgorithmics (May 2014).
C
C Updates:
C
C   January 2025 (Update to version 2): Corrected mass balance reporting for storage changes when only one of Ss or Sy
C       is changed. Added error checking for logarithm bases <= 1 and unsupported ALTSTO flag.
C   September 2017 (Update to version 2): TVM now interoperates correctly with the HFB6 package when both are enabled.
C       Fixed reporting issues when values are changed at the start of the first stress period.
C   January 2016 (Version 2): Storage changes using linear and log interpolation are now incorporated into formulate
C       and budget subroutines to conserve volume. Listing outputs now correctly report number of SC1 and SC2 changes.
C
C
C
      SUBROUTINE TVMU2ARPRE
C     ******************************************************************
C     SAVE INITIAL GEOMETRIC FACTORS FROM PGF ARRAY PRIOR TO THEIR
C     MODIFICATION IN THE LPF PACKAGE
C     ******************************************************************
C
C     SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,ONLY:PGF,NJAS
      USE TVMU2MODULE
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
      SUBROUTINE TVMU2AR(IN,IUNITLPF,IUNITSFR,IUNITDPF,INBCT)
C     ******************************************************************
C     ALLOCATE ARRAY STORAGE FOR TIME-VARIANT MATERIALS
C     ******************************************************************
C
C     SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,ONLY:IOUT,IFREFM,ITRSS,IUNSTR,NODES,NLAY
      USE GWFBCFMODULE,ONLY:IKCFLAG,LAYCON,IALTSTO
      USE GWFSFRMODULE,ONLY:ISFROPT
      USE TVMU2MODULE
      CHARACTER*400 :: LINE
      INTEGER :: KPERDUMMY
C     ------------------------------------------------------------------
C-------IDENTIFY PACKAGE
      WRITE(IOUT,1) IN
1     FORMAT(1X,/1X,'TVM -- TIME-VARIANT MATERIALS VERSION 2,',
     1  ' 9/29/2017',/1X,'INPUT READ FROM UNIT ',I4)
C-------ENSURE LPF PACKAGE IS USED
      IF(IUNITLPF.EQ.0) THEN
          WRITE (IOUT, 12) 
12        FORMAT(1X,'LPF PACKAGE MUST BE ACTIVE TO CHANGE MATERIAL',
     1    1X,'PROPERTIES OF THE POROUS MEDIUM.')
      END IF
C-------ENSURE NODAL HYDRAULIC CONDUCTIVITIES ARE USED
      IF(IKCFLAG.NE.0) THEN
          WRITE (IOUT, *) 'TVM MAY ONLY BE USED WITH NODAL ',
     1        'CONDUCTIVITIES (IKCFLAG = 0 IN LPF PACKAGE).'
          CALL USTOP(' ')
      END IF
C-------WARN ABOUT USAGE IN STEADY-STATE SIMULATION
      IF(ITRSS.NE.1) THEN
          WRITE (IOUT, *) 'WARNING: TVM PACKAGE IS ACTIVE BUT',
     1        ' STEADY-STATE STRESS PERIODS ARE IN USE!'
      END IF
C-------WARN ABOUT USAGE WITH STRUCTURED GRID
      IF(IUNSTR.EQ.0) THEN
          WRITE (IOUT, *) 'WARNING: TVM PACKAGE EXPECTS UNSTRUCTURED',
     1        ' CELL NUMBERS, BUT STRUCTURED GRID INPUT IS BEING USED.'
      END IF
C-------WARN ABOUT UNSATURATED SFR OPTIONS
      IF(IUNITSFR.GT.0) THEN
        IF(ISFROPT.GT.1) THEN
          WRITE (IOUT, *) 'WARNING: SFR PACKAGE IS ACTIVE WITH ISFROPT',
     1        ' GREATER THAN 1. THIS IS UNSUPPORTED BY THE TVM',
     2        ' PACKAGE. UNSATURATED ZONE HYDRAULIC CONDUCTIVITIES AND',
     3        ' RESIDUAL WATER CONTENT VALUES WILL NOT BE UPDATED WITH',
     4        ' TVM MATERIAL PROPERTY CHANGES!'
        END IF
      END IF
C-------ENSURE ALTSTO FLAG IS NOT USED
      IF(IALTSTO.NE.0) THEN
          WRITE (IOUT, *) 'ALTSTO FLAG IS SPECIFIED IN LPF PACKAGE, ',
     1        'WHICH IS NOT SUPPORTED BY TVM.'
          CALL USTOP(' ')
      END IF
C-------DECIDE WHETHER SPECIFIC YIELD IS USED IN THE SIMULATION
      IHAVESC2 = 0
      DO I = 1, NLAY
          L = LAYCON(I)
          IF(L.EQ.2 .OR. L.EQ.3 .OR. L.EQ.4 .OR. L.EQ.5) THEN
              IHAVESC2 = 1
          END IF
      END DO
C-------READ AND PRINT COMMENTS, THEN READ FIRST LINE
      CALL URDCOM(IN,IOUT,LINE)
C-------READ PRINT FLAG
      IF(INBCT.GT.0) THEN !FOR TRANSPORT READ POROSITY
        IF(IFREFM.EQ.0) THEN
          READ(LINE,'(I10,6F10.0)') ITVMPRINT,TVMLOGBASEHK,
     1      TVMLOGBASEVKA,TVMLOGBASESS,TVMLOGBASESY,TVMDDFTR,
     1      TVMLOGBASEPOR
        ELSE
          READ(LINE,*) ITVMPRINT,TVMLOGBASEHK,TVMLOGBASEVKA,
     1      TVMLOGBASESS,TVMLOGBASESY,TVMDDFTR,
     1      TVMLOGBASEPOR
        END IF
      ELSE
        IF(IFREFM.EQ.0) THEN
          READ(LINE,'(I10,5F10.0)') ITVMPRINT,TVMLOGBASEHK,
     1      TVMLOGBASEVKA,TVMLOGBASESS,TVMLOGBASESY,TVMDDFTR
        ELSE
          READ(LINE,*) ITVMPRINT,TVMLOGBASEHK,TVMLOGBASEVKA,
     1      TVMLOGBASESS,TVMLOGBASESY,TVMDDFTR
        END IF
      ENDIF
C-------REFLECT INPUT AND PRE-CALCULATE 1/LOG(TVMLOGBASE)
      !---TVMLOGBASEHK-------------------------------------------------
      IF(TVMLOGBASEHK.LT.0) THEN
          WRITE (IOUT, *)'STEP FUNCTION FOR HK IS USED, TVMLOGBASEHK = '
     1     ,TVMLOGBASEHK
      ELSE IF(TVMLOGBASEHK.EQ.0) THEN
          WRITE (IOUT, *) 'TVM INTERPOLATES HK VALUES LINEARLY',
     1        ' (TVMLOGBASEHK = 0)'
          TVMLOGBASEHKINV = 1.0
      ELSE IF(TVMLOGBASEHK.LE.1.) THEN
          WRITE (IOUT, *) 'TVMLOGBASEHK MUST BE > 1 (LOGARITHM BASE),',
     1        ' < 0 (STEP FUNCTION) OR = 0 (LINEAR).',
     1        ' TVMLOGBASEHK = ', TVMLOGBASEHK
          CALL USTOP(' ')
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
      ELSE IF(TVMLOGBASEVKA.LE.1.) THEN
          WRITE (IOUT, *) 'TVMLOGBASEVKA MUST BE > 1 (LOGARITHM BASE),',
     1        ' < 0 (STEP FUNCTION) OR = 0 (LINEAR).',
     1        ' TVMLOGBASEVKA = ', TVMLOGBASEVKA
          CALL USTOP(' ')
      ELSE
          WRITE (IOUT, *) 'TVM INTERPOLATES VKA VALUES ',
     1        'LOGARITHMICALLY, TVMLOGBASEVKA = ', TVMLOGBASEVKA 
          TVMLOGBASEVKAINV = 1.0 / LOG(TVMLOGBASEVKA)
      END IF
      !---TVMLOGBASESS-------------------------------------------------
      IF(TVMLOGBASESS.LT.0) THEN
          WRITE (IOUT, *) 'TVM DOES NOT SUPPORT STEP FUNCTION FOR SS', 
     1     ' (TVMLOGBASESS = ',TVMLOGBASESS,')'
          CALL USTOP(' ')
      ELSE IF(TVMLOGBASESS.EQ.0) THEN
          WRITE (IOUT, *) 'TVM INTERPOLATES SS VALUES LINEARLY',
     1        ' (TVMLOGBASESS = 0)'
          TVMLOGBASESSINV = 1.0
      ELSE IF(TVMLOGBASESS.LE.1.) THEN
          WRITE (IOUT, *) 'TVMLOGBASESS MUST BE > 1 (LOGARITHM BASE)',
     1        ' OR = 0 (LINEAR).',
     1        ' TVMLOGBASESS = ', TVMLOGBASESS
          CALL USTOP(' ')
      ELSE
          WRITE (IOUT, *) 'TVM INTERPOLATES SS VALUES ',
     1        'LOGARITHMICALLY, TVMLOGBASESS = ', TVMLOGBASESS
          TVMLOGBASESSINV = 1.0 / LOG(TVMLOGBASESS)
      END IF
      !---TVMLOGBASESY-------------------------------------------------
      IF(TVMLOGBASESY.LT.0) THEN
          WRITE (IOUT, *) 'TVM DOES NOT SUPPORT STEP FUNCTION FOR SY', 
     1     ' (TVMLOGBASESY = ',TVMLOGBASESY,')'
          CALL USTOP(' ')
      ELSE IF(TVMLOGBASESY.EQ.0) THEN
          WRITE (IOUT, *) 'TVM INTERPOLATES SY VALUES LINEARLY',
     1        ' (TVMLOGBASESY = 0)'
          TVMLOGBASESYINV = 1.0
      ELSE IF(TVMLOGBASESY.LE.1.) THEN
          WRITE (IOUT, *) 'TVMLOGBASESY MUST BE > 1 (LOGARITHM BASE)',
     1        ' OR = 0 (LINEAR).',
     1        ' TVMLOGBASESY = ', TVMLOGBASESY
          CALL USTOP(' ')
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
      ELSE IF(TVMDDFTR.LE.1.) THEN
          WRITE (IOUT, *) 'TVMDDFTR MUST BE > 1 (LOGARITHM BASE),',
     1        ' < 0 (STEP FUNCTION) OR = 0 (LINEAR).',
     1        ' TVMDDFTR = ', TVMDDFTR
          CALL USTOP(' ')
      ELSE
          WRITE (IOUT, *) 'TVM INTERPOLATES DDFTR VALUES ',
     1        'LOGARITHMICALLY, TVMDDFTR = ', TVMDDFTR
          TVMDDFTRINV = 1.0 / LOG(TVMDDFTR)
      END IF
100   CONTINUE   
      !---TVMLOGBASEPOR-------------------------------------------------
      IF(INBCT.EQ.0) GO TO 200 
      IF(TVMLOGBASEPOR.LT.0) THEN
          WRITE (IOUT, *) 'TVM DOES NOT SUPPORT STEP FUNCTION FOR ', 
     1     'POROSITY (TVMLOGBASEPOR = ',TVMLOGBASEPOR,')'
          CALL USTOP(' ')
      ELSE IF(TVMLOGBASEPOR.EQ.0) THEN
          WRITE (IOUT, *) 'TVM INTERPOLATES POROSITY VALUES LINEARLY',
     1        ' (TVMLOGBASEPOR = 0)'
          TVMLOGBASEPORINV = 1.0
      ELSE IF(TVMLOGBASEPOR.LE.1.) THEN
          WRITE (IOUT, *) 'TVMLOGBASEPOR MUST BE > 1 (LOGARITHM BASE)',
     1        ' OR = 0 (LINEAR).',
     1        ' TVMLOGBASEPOR = ', TVMLOGBASEPOR
          CALL USTOP(' ')
      ELSE
          WRITE (IOUT, *) 'TVM INTERPOLATES POROSITY VALUES ',
     1        'LOGARITHMICALLY, TVMLOGBASEPOR = ', TVMLOGBASEPOR
          TVMLOGBASEPORINV = 1.0 / LOG(TVMLOGBASEPOR)
      END IF
200   CONTINUE   
C-------INIT VARIABLES
      NTVMHK = 0
      NTVMVKA = 0
      NTVMSS = 0
      NTVMSY = 0
      NTVMDDFTR = 0
      NTVMPOR = 0
C-------Allocate memory for old storage values so that storage changes can be formulated mass-conservatively
      ALLOCATE(SSOLD(NODES))
      ALLOCATE(SYOLD(NODES))
      IF(INBCT.GT.0) ALLOCATE(POROLD(NODES))
C-------READ INITIAL CHANGED VALUES (START OF STRESS PERIOD 1)
      KPERDUMMY = 0
      CALL TVMU2READ(IN,KPERDUMMY,INBCT)
C-------RETURN.
      RETURN
      END
C
C
      SUBROUTINE TVMU2READ(IN,KPER,INBCT)
C     ******************************************************************
C     READ ONE SET OF VALUES FOR INTERPOLATION
C     ******************************************************************
C
C     SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,ONLY:IOUT,IFREFM,TOP,BOT,AREA
      USE GWFBCFMODULE,ONLY:SC1,SC2,VKA,HK,ISFAC
      USE TVMU2MODULE
      USE GWFDPFMODULE, ONLY: DDFTR
      USE GWTBCTMODULE, ONLY: PRSITY
C     ------------------------------------------------------------------
C-------DEALLOCATE ANY EXISTING ARRAYS PRIOR TO READING NEW ARRAY VALUES
      CALL TVMU2DASP
C-------IDENTIFY PACKAGE
      WRITE(IOUT,*) 'TVM READING FROM UNIT ',IN
C-------READ NUMBER OF VALUES TO CHANGE FOR EACH MATERIAL PROPERTY
      IF(INBCT.GT.0) THEN !FOR TRANSPORT READ POROSITY
        IF(IFREFM.EQ.0) THEN
          READ(IN,'(6I10)') NTVMHK,NTVMVKA,NTVMSS,NTVMSY,NTVMDDFTR,
     1      NTVMPOR
        ELSE
          READ(IN,*) NTVMHK,NTVMVKA,NTVMSS,NTVMSY,NTVMDDFTR,
     1      NTVMPOR
        END IF
C-------WRITE SUMMARY OF THIS SET OF VALUES TO OUT FILE
        IF(ITVMPRINT.GE.1) THEN
           WRITE(IOUT,21) NTVMHK,NTVMVKA,NTVMSS,NTVMSY,NTVMDDFTR,
     1     NTVMPOR,KPER
21          FORMAT(1X,'READING',I8,' HK, ',I8,' VKA, ',I8,' SS, ',I8,
     1     ' SY ',I8,' DDFTR AND',I8,' POROSITY VALUES FOR STRESS',
     1     ' PERIOD ',I8)      
        END IF
      ELSE
        IF(IFREFM.EQ.0) THEN
          READ(IN,'(5I10)') NTVMHK,NTVMVKA,NTVMSS,NTVMSY,NTVMDDFTR
        ELSE
          READ(IN,*) NTVMHK,NTVMVKA,NTVMSS,NTVMSY,NTVMDDFTR
        END IF
C-------WRITE SUMMARY OF THIS SET OF VALUES TO OUT FILE
        IF(ITVMPRINT.GE.1) THEN
           WRITE(IOUT,1) NTVMHK,NTVMVKA,NTVMSS,NTVMSY,NTVMDDFTR,KPER
1          FORMAT(1X,'READING',I8,' HK, ',I8,' VKA, ',I8,' SS, ',I8,
     1     ' SY AND ',I8,' DDFTR VALUES FOR STRESS PERIOD ',I8)      
        END IF
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
      IF(NTVMSS.GT.0.OR.NTVMSY.GT.0) THEN
          ALLOCATE(ITVMSMAPSS(NTVMSS+NTVMSY))
          ALLOCATE(ITVMSMAPSY(NTVMSS+NTVMSY))
      END IF
      IF(NTVMPOR.GT.0) THEN
          ALLOCATE(ITVMPOR(NTVMPOR))
          ALLOCATE(PORNEWA(NTVMPOR))
          ALLOCATE(PORNEWB(NTVMPOR))
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
                  WRITE(IOUT,2) ITVMHK(I),HKNEWA(I),HKNEWB(I)   
2                 FORMAT('  HK FOR CELL ',I8,' WILL BE INTERPOLATED',
     1            1X,'FROM', G15.6,' TO ',G15.6)       
                ELSE
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
                  WRITE(IOUT,5)ITVMVKA(I),VKANEWA(I),VKANEWB(I)
5                 FORMAT(' VKA FOR CELL ',I8,' WILL BE INTERPOLATED',
     1            1X,'FROM', G15.6,' TO ',G15.6)                    
                ELSE
                  WRITE(IOUT,6)ITVMVKA(I),VKANEWB(I)  
6                 FORMAT(' VKA FOR CELL ',I8,' WILL BE SET TO ',G15.6) 
                ENDIF    
              ELSE
                WRITE(IOUT,7) ITVMVKA(I),VKANEWB(I)   
7               FORMAT(' VKA FOR CELL ',I8,' WILL BE CHANGED TO ',G15.6)
              END IF
          END IF
      END DO
C-------READ NEW SS VALUES (SC1)
      NTVMSCHANGES = 0
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
          NTVMSCHANGES = NTVMSCHANGES + 1
          ITVMSMAPSS(NTVMSCHANGES) = I
          ITVMSMAPSY(NTVMSCHANGES) = 0
          IF(ITVMPRINT.GE.2) THEN
              IF(KPER.GE.1) THEN
                IF(TVMLOGBASESS.GE.0)THEN
                  WRITE(IOUT,8) ITVMSS(I),SSNEWA(I),SSNEWB(I) 
8                 FORMAT(' SC1 FOR CELL ',I8,' WILL BE INTERPOLATED',
     1            1X,'FROM', G15.6,' TO ',G15.6) 
                ELSE
                  WRITE(IOUT,9) ITVMSS(I),SSNEWB(I)               
9                 FORMAT(' SC1 FOR CELL ',I8,' WILL BE SET TO ',G15.6) 
                ENDIF
              ELSE
                WRITE(IOUT,10) ITVMSS(I),SSNEWB(I)
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
C-----------Locate corresponding SS change (if any)
          IMAPPEDS = 0
          IF(NTVMSCHANGES.GE.I) THEN ! Quick evaluation if SS and SY changed cells match in order
            IF(ITVMSS(ITVMSMAPSS(I)).EQ.ITVMSY(I)) THEN
                IMAPPEDS = I
            END IF
          END IF
          IF(IMAPPEDS.EQ.0) THEN     ! Search for corresponding SS change
            DO J = 1, NTVMSS
              IF(J.GT.NTVMSCHANGES) THEN
                EXIT
              END IF
              IF(ITVMSS(ITVMSMAPSS(J)).EQ.ITVMSY(I)) THEN
                IMAPPEDS = J
              END IF
            END DO
          END IF
          IF(IMAPPEDS.EQ.0) THEN     ! No matching SS change, so introduce a new one just for SY
            NTVMSCHANGES = NTVMSCHANGES + 1
            IMAPPEDS = NTVMSCHANGES
            ITVMSMAPSS(IMAPPEDS) = 0
          END IF
          ITVMSMAPSY(IMAPPEDS) = I
C
          IF(ITVMPRINT.GE.2) THEN
              IF(KPER.GE.1) THEN
                IF(TVMLOGBASESY.GE.0)THEN  
                  WRITE(IOUT,11) ITVMSY(I),SYNEWA(I),SYNEWB(I)
11                FORMAT(' SC2 FOR CELL ',I8,' WILL BE INTERPOLATED',
     1            1X,'FROM', G15.6,' TO ',G15.6) 
                ELSE
                  WRITE(IOUT,12) ITVMSY(I),SYNEWB(I)
12                FORMAT(' SC2 FOR CELL ',I8,' WILL BE SET TO ',G15.6)
                ENDIF
              ELSE
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
                 WRITE(IOUT,14) ITVMDDFTR(I),DDFTRNEWA(I),DDFTRNEWB(I)
14                FORMAT(' DDFTR FOR CELL ',I8,' WILL BE INTERPOLATED',
     1            1X,'FROM', G15.6,' TO ',G15.6) 
                ELSE
                 WRITE(IOUT,15) ITVMDDFTR(I),DDFTRNEWB(I) 
15               FORMAT(' DDFTR FOR CELL ',I8,' WILL BE SET TO ',G15.6)
                ENDIF
              ELSE
                WRITE(IOUT,16)ITVMDDFTR(I),DDFTRNEWB(I)  
16            FORMAT(' DDFTR FOR CELL ',I8,' WILL BE CHANGED TO ',G15.6)
              END IF
          END IF
      END DO      
C-------READ NEW POROSITY VALUES 
      DO I = 1, NTVMPOR
          IF(IFREFM.EQ.0) THEN
            READ(IN,'(I10,F10.0)') ITVMPOR(I), PORNEWB(I)
          ELSE
            READ(IN,*) ITVMPOR(I), PORNEWB(I)
          END IF
          PORNEWA(I) = PRSITY(ITVMPOR(I))
          IF(ITVMPRINT.GE.2) THEN
            IF(KPER.GE.1) THEN
              IF(TVMLOGBASEPOR.GE.0)THEN
                WRITE(IOUT,17) ITVMPOR(I),PORNEWA(I),PORNEWB(I)
17              FORMAT(' POROSITY FOR CELL ',I8,' WILL BE INTERPOLATED',
     1                 1X,'FROM', G15.6,' TO ',G15.6) 
              ELSE
                WRITE(IOUT,18) ITVMPOR(I),PORNEWB(I) 
18              FORMAT(' POROSITY FOR CELL ',I8,' WILL BE SET TO ',
     1                 G15.6)
              ENDIF
            ELSE
              WRITE(IOUT,19)ITVMPOR(I),PORNEWB(I)  
19            FORMAT(' POROSITY FOR CELL ',I8,' WILL BE CHANGED TO ',
     1               G15.6)
            END IF
          END IF
      END DO      
C-------RETURN
      RETURN
      END
C
C
      SUBROUTINE TVMU2RP(IN,IUNITCLN,IUNITSFR,IUNITDPF,IUNITHFB,KPER,
     1  INBCT)
C     ******************************************************************
C     READ MATERIAL PROPERTY VALUES FOR CURRENT STRESS PERIOD
C     ******************************************************************
      USE GLOBAL,ONLY:IOUT,NLAY,NJAS,PGF,NODES
      USE GWFBCFMODULE,ONLY:SC1,SC2,VKA,HK
      USE TVMU2MODULE
      USE GWFDPFMODULE, ONLY: DDFTR
      USE GWTBCTMODULE, ONLY: PRSITY
C-------IF A STEP FUNCTION IS USED, STORE THE PREVIOUS SP'S STORAGE VALUES
      IF(TVMLOGBASESS.LT.0) THEN
          DO I = 1, NODES
              SSOLD(I) = SC1(I)
          END DO
      END IF
      IF(TVMLOGBASESY.LT.0.AND.IHAVESC2.NE.0) THEN
          DO I = 1, NODES
              SYOLD(I) = SC2(I)
          END DO
      END IF
C-------APPLY FINAL VALUES FROM PREVIOUS STRESS PERIOD, IF ANY
C-------CHANGE HK VALUES
      IF(NTVMHK.GT.0) THEN
          IF(ITVMPRINT.GE.1) THEN
            IF(KPER.GT.1) THEN
              WRITE(IOUT,1)NTVMHK,(KPER - 1) 
1             FORMAT(' TVM CHANGING ',I8,' HK VALUES AT END OF STRESS',
     1        1X,'PERIOD ',I8)
            ELSE
              WRITE(IOUT,2)NTVMHK   
2             FORMAT(' TVM CHANGING ',I8,' HK VALUES AT START OF STRESS'
     1        ,1X,'PERIOD  1')                 
            END IF
          END IF
          DO I = 1, NTVMHK
              HK(ITVMHK(I)) = HKNEWB(I)
              IF(ITVMPRINT.GE.2) THEN
                IF(KPER.GT.1) THEN
                  WRITE(IOUT,3) ITVMHK(I),HK(ITVMHK(I)),(KPER - 1)  
3                FORMAT('  HK FOR CELL ',I8,' CHANGED TO ',G15.6,
     1           ' AT END OF STRESS PERIOD ',I8) 
                ELSE
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
              WRITE(IOUT,11)NTVMVKA,(KPER - 1)              
11            FORMAT(' TVM CHANGING ',I8,' VKA VALUES AT END OF STRESS',
     1        1X,'PERIOD ',I8)              
            ELSE
              WRITE(IOUT,12) NTVMVKA      
12            FORMAT(' TVM CHANGING ',I8,' VKA VALUES AT START OF'
     1        ,1X,'STRESS PERIOD  1')                
            END IF
          END IF
          DO I = 1, NTVMVKA
              VKA(ITVMVKA(I)) = VKANEWB(I)
              IF(ITVMPRINT.GE.2) THEN
                IF(KPER.GT.1) THEN
                  WRITE(IOUT,13) ITVMVKA(I),VKA(ITVMVKA(I)),(KPER - 1)
13                FORMAT('  VKA FOR CELL ',I8,' CHANGED TO ',G15.6,
     1           ' AT END OF STRESS PERIOD ',I8) 
                ELSE
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
               WRITE(IOUT,21) NTVMSS,(KPER - 1)       
21            FORMAT(' TVM CHANGING ',I8,' SC1 VALUES AT END OF STRESS',
     1        1X,'PERIOD ',I8)                
            ELSE
              WRITE(IOUT,22) NTVMSS    
22            FORMAT(' TVM CHANGING ',I8,' SC1 VALUES AT START OF'
     1        ,1X,'STRESS PERIOD  1')               
            END IF
          END IF
          DO I = 1, NTVMSS
              SC1(ITVMSS(I)) = SSNEWB(I)
              IF(ITVMPRINT.GE.2) THEN
                IF(KPER.GT.1) THEN
                  WRITE(IOUT,23) ITVMSS(I),SC1(ITVMSS(I)),(KPER - 1)  
23                FORMAT('  SC1 FOR CELL ',I8,' CHANGED TO ',G15.6,
     1           ' AT END OF STRESS PERIOD ',I8)                   
                ELSE
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
              WRITE(IOUT,31)NTVMSY,(KPER - 1)         
31            FORMAT(' TVM CHANGING ',I8,' SC2 VALUES AT END OF STRESS',
     1        1X,'PERIOD ',I8)              
            ELSE
              WRITE(IOUT,32) NTVMSY   
32            FORMAT(' TVM CHANGING ',I8,' SC2 VALUES AT START OF'
     1        ,1X,'STRESS PERIOD  1')              
            END IF
          END IF
          DO I = 1, NTVMSY
              SC2(ITVMSY(I)) = SYNEWB(I)
              IF(ITVMPRINT.GE.2) THEN
                IF(KPER.GT.1) THEN
                  WRITE(IOUT,33)ITVMSY(I),SC2(ITVMSY(I)),(KPER - 1)
33                FORMAT('  SC2 FOR CELL ',I8,' CHANGED TO ',G15.6,
     1           ' AT END OF STRESS PERIOD ',I8)                    
                ELSE
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
              WRITE(IOUT,41) NTVMDDFTR,(KPER - 1)      
41            FORMAT(' TVM CHANGING ',I8,' DDFTR VALUES AT END OF',
     1        1X,'STRESS PERIOD ',I8)               
            ELSE
              WRITE(IOUT,42)NTVMDDFTR       
42            FORMAT(' TVM CHANGING ',I8,' DDFTR VALUES AT START OF'
     1        ,1X,'STRESS PERIOD  1')                 
            END IF
          END IF
          DO I = 1, NTVMDDFTR
              DDFTR(ITVMDDFTR(I)) = DDFTRNEWB(I)
              IF(ITVMPRINT.GE.2) THEN
                IF(KPER.GT.1) THEN
                  WRITE(IOUT,43) ITVMDDFTR(I),DDFTR(ITVMDDFTR(I)),
     1              (KPER - 1)      
43                FORMAT('  DDFTR FOR CELL ',I8,' CHANGED TO ',G15.6,
     1           ' AT END OF STRESS PERIOD ',I8)                   
                ELSE
                  WRITE(IOUT,44)ITVMDDFTR(I),DDFTR(ITVMDDFTR(I))
44                FORMAT('  DDFTR FOR CELL ',I8,' CHANGED TO ',G15.6,
     1            ' AT START OF STRESS PERIOD 1')                   
                END IF
              END IF
          END DO
      END IF      
C-------CHANGE POROSITY VALUES (PRSITY)
      IF(NTVMPOR.GT.0) THEN
          IF(ITVMPRINT.GE.1) THEN
            IF(KPER.GT.1) THEN
              WRITE(IOUT,51)NTVMPOR,(KPER - 1)
51            FORMAT(' TVM CHANGING ',I8,' POROSITY VALUES AT END OF ',
     1        1X,'STRESS PERIOD ',I8)
            ELSE
              WRITE(IOUT,52) NTVMPOR
52            FORMAT(' TVM CHANGING ',I8,' POROSITY VALUES AT START OF '
     1        ,1X,'STRESS PERIOD  1')              
            END IF
          END IF
          DO I = 1, NTVMPOR
              PRSITY(ITVMPOR(I)) = PORNEWB(I)
              IF(ITVMPRINT.GE.2) THEN
                IF(KPER.GT.1) THEN
                  WRITE(IOUT,53)ITVMPOR(I),PRSITY(ITVMPOR(I)),(KPER - 1)
53                FORMAT('  POROSITY FOR CELL ',I8,' CHANGED TO ',G15.6,
     1           ' AT END OF STRESS PERIOD ',I8)                    
                ELSE
                  WRITE(IOUT,54)ITVMPOR(I),PRSITY(ITVMPOR(I))  
54                FORMAT('  POROSITY FOR CELL ',I8,' CHANGED TO ',G15.6,
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
C----------HFB UPDATE
          IF(IUNITHFB.GT.0) THEN
              CALL SGWF2HFB7MC(KPER)
          END IF
      END IF
C-------READ NEXT SET OF CHANGED VALUES
      CALL TVMU2READ(IN,KPER,INBCT)
C-------RETURN
      RETURN
      END
C
C
      SUBROUTINE TVMU2AD(IN,IUNITCLN,IUNITSFR,IUNITDPF,IUNITHFB,KPER,
     1  INBCT)
C     ******************************************************************
C     SET MATERIAL PROPERTY VALUES AT CURRENT TIME STEP
C     ******************************************************************
      USE GLOBAL,ONLY:IOUT,NLAY,NJAS,PGF,NODES,AREA,LAYNOD
      USE GWFBASMODULE,ONLY:TOTIM
      USE GWFBCFMODULE,ONLY:SC1,SC2,VKA,HK,ISFAC,LAYCON
      USE TVMU2MODULE
      USE GWFDPFMODULE, ONLY: DDFTR
      USE GWTBCTMODULE, ONLY: PRSITY
C-------IF LINEAR OR LOG INTERPOLATION IS USED (NOT A STEP FUNCTION),
C-------STORE THE PREVIOUS TIME STEP'S STORAGE VALUES
      IF(TVMLOGBASESS.GE.0) THEN
          DO I = 1, NODES
              SSOLD(I) = SC1(I)
          END DO
      END IF
      IF(TVMLOGBASESY.GE.0.AND.IHAVESC2.NE.0) THEN
          DO I = 1, NODES
              SYOLD(I) = SC2(I)
          END DO
      END IF
      IF(INBCT.GT.0) THEN
        IF(TVMLOGBASEPOR.GE.0) THEN
          DO I = 1, NODES
              POROLD(I) = PRSITY(I)
          END DO
        END IF
      END IF
C-------NOTHING TO CHANGE IF ALL COUNTS ARE ZERO
      IF(NTVMHK.LE.0.AND.NTVMVKA.LE.0.AND.NTVMSS.LE.0.AND.
     1   NTVMSY.LE.0.AND.NTVMDDFTR.LE.0.AND.NTVMPOR.LE.0) THEN
          IF(ITVMPRINT.GE.1) THEN
              WRITE(IOUT,*) 'TVM LEAVING ALL MATERIAL VALUES UNCHANGED',
     1            ' AT TIME ',TOTIM,' AS ALL COUNTS ARE ZERO'
          END IF
          RETURN
      END IF
C-------CHANGE HK VALUES IF NOT STEP FUNCTION
      IF(NTVMHK.GT.0.AND.TVMLOGBASEHK.GE.0) THEN
          IF(ITVMPRINT.GE.1) THEN
              WRITE(IOUT,1)NTVMHK,TOTIM
1             FORMAT(' TVM CHANGING ',I8,' HK VALUES AT TIME ',G15.6)  
          END IF
          DO I = 1, NTVMHK
              HK(ITVMHK(I)) = TVMU2INTERP(HKNEWA(I),HKNEWB(I),KPER,
     1                            TVMLOGBASEHK,TVMLOGBASEHKINV)
              IF(ITVMPRINT.GE.2) THEN
                  WRITE(IOUT,2) ITVMHK(I),HK(ITVMHK(I)),TOTIM
2                 FORMAT('  HK FOR CELL ',I8,' CHANGED TO ',G15.6,
     1             ' AT TIME ',G15.6)                   
              END IF
          END DO
      END IF
C-------CHANGE VKA VALUES
      IF(NTVMVKA.GT.0.AND.TVMLOGBASEVKA.GE.0) THEN
          IF(ITVMPRINT.GE.1) THEN
              WRITE(IOUT,11)NTVMVKA,TOTIM  
11            FORMAT(' TVM CHANGING ',I8,' VKA VALUES AT TIME ',G15.6)
          END IF
          DO I = 1, NTVMVKA
              VKA(ITVMVKA(I)) = TVMU2INTERP(VKANEWA(I),VKANEWB(I),KPER,
     1                              TVMLOGBASEVKA,TVMLOGBASEVKAINV)
              IF(ITVMPRINT.GE.2) THEN
                  WRITE(IOUT,12) ITVMVKA(I),VKA(ITVMVKA(I)),TOTIM       
12                FORMAT('  VKA FOR CELL ',I8,' CHANGED TO ',G15.6,
     1             ' AT TIME ',G15.6)                   
              END IF
          END DO
      END IF
C-------CHANGE SS VALUES (SC1)
      IF(NTVMSS.GT.0.AND.TVMLOGBASESS.GE.0) THEN
          IF(ITVMPRINT.GE.1) THEN
              WRITE(IOUT,21)NTVMSS,TOTIM
 21            FORMAT(' TVM CHANGING ',I8,' SC1 VALUES AT TIME ',G15.6)
          END IF
          DO I = 1, NTVMSS
              SC1(ITVMSS(I)) = TVMU2INTERP(SSNEWA(I),SSNEWB(I),KPER,
     1                             TVMLOGBASESS,TVMLOGBASESSINV)
              IF(ITVMPRINT.GE.2) THEN
                  WRITE(IOUT,22)ITVMSS(I),SC1(ITVMSS(I)),TOTIM
22                FORMAT('  SC1 FOR CELL ',I8,' CHANGED TO ',G15.6,
     1             ' AT TIME ',G15.6)                  
              END IF
          END DO
      END IF
C-------CHANGE SY VALUES (SC2)
      IF(NTVMSY.GT.0.AND.TVMLOGBASESY.GE.0) THEN
          IF(ITVMPRINT.GE.1) THEN
              WRITE(IOUT,31) NTVMSY,TOTIM    
31            FORMAT(' TVM CHANGING ',I8,' SC2 VALUES AT TIME ',G15.6)
          END IF
          DO I = 1, NTVMSY
              SC2(ITVMSY(I)) = TVMU2INTERP(SYNEWA(I),SYNEWB(I),KPER,
     1                             TVMLOGBASESY,TVMLOGBASESYINV)
              IF(ITVMPRINT.GE.2) THEN
                  WRITE(IOUT,32)ITVMSY(I),SC2(ITVMSY(I)),TOTIM
32                FORMAT('  SC2 FOR CELL ',I8,' CHANGED TO ',G15.6,
     1             ' AT TIME ',G15.6)                   
              END IF
          END DO
      END IF
C-------CHANGE DDFTR VALUES 
      IF(NTVMDDFTR.GT.0.AND.TVMDDFTR.GE.0) THEN
          IF(ITVMPRINT.GE.1) THEN
              WRITE(IOUT,41) NTVMDDFTR,TOTIM
41            FORMAT(' TVM CHANGING ',I8,' DDFTR VALUES AT TIME ',G15.6)
          END IF
          DO I = 1, NTVMDDFTR
              DDFTR(ITVMDDFTR(I)) = TVMU2INTERP(DDFTRNEWA(I),
     1             DDFTRNEWB(I),KPER,TVMDDFTR,TVMDDFTRINV)
              IF(ITVMPRINT.GE.2) THEN
                  WRITE(IOUT,42) ITVMDDFTR(I),DDFTR(ITVMDDFTR(I)),TOTIM
42                FORMAT('  DDFTR FOR CELL ',I8,' CHANGED TO ',G15.6,
     1             ' AT TIME ',G15.6)                   
              END IF
          END DO
      END IF
C-------CHANGE POROSITY VALUES (PRSITY)
      IF(NTVMPOR.GT.0.AND.TVMLOGBASEPOR.GE.0) THEN
        IF(ITVMPRINT.GE.1) THEN
          WRITE(IOUT,51) NTVMPOR,TOTIM    
51        FORMAT(' TVM CHANGING ',I8,' POROSITY VALUES AT TIME ',G15.6)
        END IF
        DO I = 1, NTVMPOR
          PRSITY(ITVMPOR(I)) = TVMU2INTERP(PORNEWA(I),PORNEWB(I),KPER,
     1                         TVMLOGBASEPOR,TVMLOGBASEPORINV)
          IF(ITVMPRINT.GE.2) THEN
            WRITE(IOUT,52)ITVMPOR(I),PRSITY(ITVMPOR(I)),TOTIM
52          FORMAT('  POROSITY FOR CELL ',I8,' CHANGED TO ',G15.6,
     1             ' AT TIME ',G15.6)                   
          END IF
        END DO
      END IF
C------------------------------------------------------------------
C-----CHECK FOR ERRORS, POROSITY CANNOT BE LESS THAN SPECIFIC YIELD
      IF(NTVMSY.GT.0.AND.TVMLOGBASESY.GE.0) THEN
        IF(INBCT.GT.0) THEN
          DO I = 1, NTVMSY
                N=ITVMSY(I)
                K=LAYNOD(N)
                IF(LAYCON(K).NE.3.AND.LAYCON(K).NE.2.AND.LAYCON(K).NE.4
     1          .AND.LAYCON(K).NE.5)GO TO 90
                IF(PRSITY(N)+1.e-5.LT.SC2(N)/AREA(N)) THEN
                  WRITE(IOUT,56) N
56                FORMAT(5X,'*** Porosity < Sy FOR NODE: ',I10,' ***')
                  STOP
                ENDIF
90              CONTINUE
          END DO
        ENDIF
      END IF
      IF(NTVMPOR.GT.0.AND.TVMLOGBASEPOR.GE.0) THEN
          DO I = 1, NTVMPOR
                N=ITVMPOR(I)
                K=LAYNOD(N)
                IF(LAYCON(K).NE.3.AND.LAYCON(K).NE.2.AND.LAYCON(K).NE.4
     1          .AND.LAYCON(K).NE.5)GO TO 91
                IF(PRSITY(N)+1.e-5.LT.SC2(N)/AREA(N)) THEN
                  WRITE(IOUT,57) N
57                FORMAT(5X,'*** Porosity < Sy FOR NODE: ',I10,' ***')
                  STOP
                ENDIF
91              CONTINUE
          END DO
      END IF
C------------------------------------------------------------------
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
C----------HFB UPDATE
          IF(IUNITHFB.GT.0) THEN
              CALL SGWF2HFB7MC(KPER)
          END IF
      END IF
C
      RETURN
      END
C
C
      SUBROUTINE TVMU2FM(KITER,KSTP,KPER)
C     ******************************************************************
C     ADJUST LPF FORMULATION OF EQUATIONS TO ACCOUNT FOR STORAGE CHANGES
C     ******************************************************************
      USE GLOBAL,ONLY:IOUT,NLAY,NJAS,IBOUND,NODLAY,AMAT,RHS,HOLD,HNEW,
     1    TOP,BOT,IA,So,Sn
      USE GWFBASMODULE,ONLY:DELT
      USE GWFBCFMODULE,ONLY:SC1,SC2,LAYCON
      USE TVMU2MODULE
      DOUBLE PRECISION RHO,RHO1,RHO2,THCK,SOLD,SNEW,SOLDADJ,TOTTHICK,DS
      DOUBLE PRECISION SSPREV,SSCUR,SYPREV,SYCUR,SATO,SATN,BBOT,TTOP,HD,
     * TLED,RHOOLD,TP,BT,EPS
C
      TLED = 1.0D0 / DELT
C-------Adjust storage changes from LPF FM---------------------------------------------------------
      DO 200 I = 1, NTVMSCHANGES
C-------Retrieve mappings to SS and SY changes
      ISSINDEX = ITVMSMAPSS(I)
      ISYINDEX = ITVMSMAPSY(I)
      IF(ISSINDEX.EQ.0) THEN
          N = ITVMSY(ISYINDEX)
      ELSE
          N = ITVMSS(ISSINDEX)
      END IF
C-------Skip external cells
      IF(IBOUND(N).LE.0) GO TO 200
C-------Storage change volume correction is not supported with a step function in this version of TVM
      IF(TVMLOGBASESS.LT.0) THEN
!          IF(KPER.EQ.1.OR.KSTP.NE.1) THEN
              ISSINDEX = 0
!          END IF
      END IF
      IF(TVMLOGBASESY.LT.0) THEN
!          IF(KPER.EQ.1.OR.KSTP.NE.1) THEN
              ISYINDEX = 0
!          END IF
      END IF
      IF(ISSINDEX.LE.0.AND.ISYINDEX.LE.0) THEN
          GO TO 200
      END IF
C-------Retrieve previous and current SS, SY values
      SSCUR = SC1(N)
      IF(IHAVESC2.NE.0) THEN
          SYCUR = SC2(N)
      END IF
      IF(ISSINDEX.GT.0) THEN
          SSPREV = SSOLD(N)
      ELSE
          SSPREV = SSCUR
      END IF
      IF(ISYINDEX.GT.0) THEN
          SYPREV = SYOLD(N)
      ELSE
          SYPREV = SYCUR
      END IF
C-------Find layer of node
      DO K=1, NLAY
          NNDLAY = NODLAY(K)
          NSTRT = NODLAY(K-1)+1
          IF(N.GE.NSTRT.AND.N.LE.NNDLAY)THEN
              EXIT
          END IF
      END DO
C6A------SEE IF THIS LAYER IS CONVERTIBLE OR NON-CONVERTIBLE.
      IF(LAYCON(K).EQ.3 .OR. LAYCON(K).EQ.2) GO TO 150
      IF(LAYCON(K).EQ.4 .OR. LAYCON(K).EQ.5) GO TO 160
C6B------NON-CONVERTIBLE LAYER, SO USE PRIMARY STORAGE
C--------Reverse LPF's RHS addition and reapply using the old Ss value
      IF(ISSINDEX.GT.0) THEN
          RHO=SSCUR*TLED
          RHOOLD=SSPREV*TLED
          RHS(N)=RHS(N)+RHO*HOLD(N)-RHOOLD*HOLD(N)
      END IF
      GO TO 200
C
C6C------A CONVERTIBLE LAYER, SO CHECK OLD AND NEW HEADS TO DETERMINE
C6C------WHEN TO USE PRIMARY AND SECONDARY STORAGE
  150 CONTINUE
C
      TP=TOP(N)
      BT=BOT(N)
      RHO2=SYCUR*TLED
      RHO1=SSCUR*TLED
C
C6C1-----FIND STORAGE FACTOR AT START OF TIME STEP.
      SOLD = RHO2
      SOLDADJ = SYPREV*TLED
      IF(HOLD(N).GT.TP) THEN
          SOLD = RHO1
          SOLDADJ = SSPREV*TLED
      END IF
C
C6C2-----FIND STORAGE FACTOR AT END OF TIME STEP.
      SNEW = RHO2
      IF(HNEW(N).GT.TP) THEN
          SNEW = RHO1
      END IF
C
C--------Reverse HOLD contribution to LPF's RHS and AMAT additions
      AMAT(IA(N)) = AMAT(IA(N)) + SNEW
      RHS(N) = RHS(N) + SOLD * (HOLD(N) - TP) + SNEW * TP
C
C--------First add the volume of water that would be released from storage
C--------under previous storage conditions if head were to drop to cell bottom
      IF(HOLD(N).GT.TP) THEN                                  !-----Cell was fully confined at start of time step
          RHS(N) = RHS(N) - SYPREV * TLED * (TP - BT)             ! Remove old volume controlled by specific yield:   - (Previous SC2) * (Cell thickness)
          RHS(N) = RHS(N) - SSPREV * TLED * (HOLD(N) - TP)        ! Remove old volume controlled by specific storage: - (Previous SC1) * ((Previous head) - (Cell top))
!          RHS(N) = RHS(N) - SSPREV * TLED * (HOLD(N) - BT)        ! OR Remove old volume controlled by specific storage: - (Previous SC1) * ((Previous head) - (Cell bottom))
      ELSE                                                    !-----Cell was unconfined at start of time step
          RHS(N) = RHS(N) - SYPREV * TLED * (HOLD(N) - BT)        ! Remove old volume controlled by specific yield:   - (Previous SC2) * ((Previous head) - (Cell bottom))
      END IF
C
C--------Now subtract the volume of water that would enter into storage under new
C--------storage conditions if head were to rise from cell bottom to its new value
      IF(HNEW(N).GT.TP) THEN                                  !-----Cell is currently fully confined
          RHS(N) = RHS(N) + RHO2 * (TP - BT)                      ! Add new volume controlled by specific yield:      + (Current SC2) * (Cell thickness)
          AMAT(IA(N)) = AMAT(IA(N)) - RHO1                        ! Add new volume controlled by specific storage:    + (Current SC1) * ((Current head) - (Cell top))
          RHS(N) = RHS(N) - (RHO1 * TP)                           ! ...
!          AMAT(IA(N)) = AMAT(IA(N)) - RHO1                        ! OR Add new volume controlled by specific storage:    + (Current SC1) * ((Current head) - (Cell bottom))
!          RHS(N) = RHS(N) - (RHO1 * BT)                           ! ...
      ELSE                                                    !-----Cell is currently unconfined
          AMAT(IA(N)) = AMAT(IA(N)) - RHO2                        ! Add new volume controlled by specific yield:      + (Current SC2) * ((Current head) - (Cell bottom))
          RHS(N) = RHS(N) - RHO2 * BT                             ! ...
      END IF
C
      GO TO 200
C
C6D------A CONVERTIBLE LAYER WITH CONTINUOUS FUNCTION,
C6D------SO USE NEWTON EXPANSION OF STORGE TERM
  160 CONTINUE
C
C--------IF THE CELL IS EXTERNAL THEN SKIP IT.
      IF(IBOUND(N).LE.0) GO TO 200
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
C-------Reverse HOLD and SATO contributions to LPF's AMAT and RHS additions
C-------and reapply using old storage values
      IF(ISSINDEX.GT.0) THEN
          AMAT(IA(N)) = AMAT(IA(N))
     1               - DS*SSCUR*HOLD(N)*TLED
          RHS(N)  = RHS(N)
     1               + SATN*SSCUR*HOLD(N)*TLED
     2               - SATO*SSPREV*HOLD(N)*TLED
     3               - DS*SSCUR*HNEW(N)*HOLD(N)*TLED
      END IF
      IF(ISYINDEX.GT.0) THEN
          RHS(N)  = RHS(N)
     1               + TOTTHICK*SYCUR*SATO*TLED
     2               - TOTTHICK*SYPREV*SATO*TLED
      END IF
    ! IF(ISSINDEX.GT.0) THEN
    !     AMAT(IA(N)) = AMAT(IA(N))
    !1               - DS*SSCUR*HOLD(N)*TLED
    !2               + DS*SSPREV*HOLD(N)*TLED
    !     RHS(N)  = RHS(N)
    !1               - DS*SSCUR*HNEW(N)*HOLD(N)*TLED
    !2               + DS*SSPREV*HNEW(N)*HOLD(N)*TLED
    !3               + SATN*SSCUR*HOLD(N)*TLED
    !4               - SATN*SSPREV*HOLD(N)*TLED
    ! END IF
    ! IF(ISYINDEX.GT.0) THEN
    !     RHS(N)  = RHS(N)
    !1               + TOTTHICK*SYCUR*SATO*TLED
    !2               - TOTTHICK*SYPREV*SATO*TLED
    ! END IF
C
  200 CONTINUE
C----------------------------------------------------------------------------
C
      RETURN
      END
C
C
      SUBROUTINE TVMU2BDS(KSTP,KPER,STOIN,STOUT)
C     ******************************************************************
C     ADJUST LPF BUDGET TERMS TO ACCOUNT FOR STORAGE CHANGES
C     ******************************************************************
      USE GLOBAL,ONLY:IOUT,NLAY,NJAS,IBOUND,NODLAY,HOLD,HNEW,
     1    TOP,BOT,IA,So,Sn,BUFF,FLOWJA
      USE GWFBASMODULE,ONLY:DELT,MSUM,VBVL,VBNM
      USE GWFBCFMODULE,ONLY:SC1,SC2,LAYCON
      USE TVMU2MODULE
C
      INTEGER I,ISSINDEX,ISYINDEX,N,KPER,KSTP,K,NNDLAY,NSTRT
      DOUBLE PRECISION STOIN,STOUT,STRG,STRGADJ,THCK,BBOT,TTOP,TP,BT
      DOUBLE PRECISION RHO,RHO1,RHO2,RHO1ADJ,RHO2ADJ,TOTTHICK,SOLD,SNEW
      DOUBLE PRECISION SSPREV,SSCUR,SYPREV,SYCUR,TLED,STRGADD,STRGREM
C
      TLED = 1.0D0 / DELT
C-------Adjust storage changes from LPF BD---------------------------------------------------------
      DO 300 I = 1, NTVMSCHANGES
C-------Retrieve mappings to SS and SY changes
      ISSINDEX = ITVMSMAPSS(I)
      ISYINDEX = ITVMSMAPSY(I)
      IF(ISSINDEX.EQ.0) THEN
          N = ITVMSY(ISYINDEX)
      ELSE
          N = ITVMSS(ISSINDEX)
      END IF
C-------Skip external cells
      IF(IBOUND(N).LE.0) GO TO 300
C-------Storage change volume correction is not supported with a step function in this version of TVM
      IF(TVMLOGBASESS.LT.0) THEN
!          IF(KPER.EQ.1.OR.KSTP.NE.1) THEN
              ISSINDEX = 0
!          END IF
      END IF
      IF(TVMLOGBASESY.LT.0) THEN
!          IF(KPER.EQ.1.OR.KSTP.NE.1) THEN
              ISYINDEX = 0
!          END IF
      END IF
      IF(ISSINDEX.LE.0.AND.ISYINDEX.LE.0) THEN
          GO TO 300
      END IF
C-------Retrieve previous and current SS, SY values
      SSCUR = SC1(N)
      SYCUR = SC2(N)
      IF(ISSINDEX.GT.0) THEN
          SSPREV = SSOLD(N)
      ELSE
          SSPREV = SSCUR
      END IF
      IF(ISYINDEX.GT.0) THEN
          SYPREV = SYOLD(N)
      ELSE
          SYPREV = SYCUR
      END IF
C-------Find layer of node
      DO K=1, NLAY
          NNDLAY = NODLAY(K)
          NSTRT = NODLAY(K-1)+1
          IF(N.GE.NSTRT.AND.N.LE.NNDLAY)THEN
              EXIT
          END IF
      END DO
C-------COMBINED COMPUTATION FOR LAYCON=4 AND LAYCON=5
      IF(LAYCON(K).EQ.4.OR.LAYCON(K).EQ.5) GO TO 333
C
C7-----CHECK LAYER TYPE TO SEE IF ONE STORAGE CAPACITY OR TWO.
      IF(LAYCON(K).NE.3.AND.LAYCON(K).NE.2) GO TO 285
C
C7A----TWO STORAGE CAPACITIES.
      TP = TOP(N)
      BT=BOT(N)
      RHO2 = SC2(N) * TLED
      RHO1 = SC1(N) * TLED
      SOLD = RHO2
      IF(HOLD(N).GT.TP) THEN
          SOLD = RHO1
      END IF
      SNEW=RHO2
      IF(HNEW(N).GT.TP) THEN
          SNEW = RHO1
      END IF
      STRG = SOLD * (HOLD(N) - TP) + SNEW * TP - SNEW * HNEW(N)
C
C--------First add the volume of water that would be released from storage
C--------under previous storage conditions if head were to drop to cell bottom
      IF(HOLD(N).GT.TP) THEN                                  !-----Cell was fully confined at start of time step
          STRGADJ = SYPREV * TLED * (TP - BT)
     1               + SSPREV * TLED * (HOLD(N) - TP)
!          STRGADJ = SSPREV * TLED * (HOLD(N) - BT)
      ELSE                                                    !-----Cell was unconfined at start of time step
          STRGADJ = SYPREV * TLED * (HOLD(N) - BT)
      END IF
C
C--------Now subtract the volume of water that would enter into storage under new
C--------storage conditions if head were to rise from cell bottom to its new value
      IF(HNEW(N).GT.TP) THEN                                  !-----Cell is currently fully confined
          STRGADJ = STRGADJ - (RHO2 * (TP - BT))
     1               - (RHO1 * HNEW(N))
     2               + (RHO1 * TP)
!          STRGADJ = STRGADJ - (RHO1 * HNEW(N))
!    1                + (RHO1 * BT)
      ELSE                                                    !-----Cell is currently unconfined
          STRGADJ = STRGADJ - (RHO2 * HNEW(N))
     1               + (RHO2 * BT)
      END IF
C
      GO TO 288
C
C7B----ONE STORAGE CAPACITY.
  285 RHO = SC1(N) * TLED
      STRG = RHO * HOLD(N) - RHO * HNEW(N)
      STRGADJ = SSPREV * TLED * HOLD(N) - SSCUR * TLED * HNEW(N)
      GO TO 288
333   CONTINUE
C7C-----COMBINED COMPUTATION FOR LAYCON=4 AND LAYCON=5
      SOLD = So(N)
      BBOT = BOT(N)
      TTOP = TOP(N)
      TOTTHICK = TTOP - BBOT
      CALL SAT_THIK(N,HNEW(N),TOTTHICK,BBOT,THCK,K,TTOP)
      SNEW = THCK
      RHO2 = SC2(N) * TLED * (SOLD - SNEW) * TOTTHICK ! LPF budget calc for SC2 contribution
      RHO1 = SC1(N) * TLED * SNEW * (HOLD(N) - HNEW(N)) ! LPF budget calc for SC1 contribution
      STRGREM = 0.
      STRGADD = 0.
      IF(ISYINDEX.GT.0) THEN
          ! SC2 was changed: TVM removes old storage volume and adds new storage volume
          RHO2ADJ = TOTTHICK * TLED * (SYPREV * SOLD - SYCUR * SNEW)
          STRGREM = STRGREM + RHO2
          STRGADD = STRGADD + RHO2ADJ
      ELSE
          ! SC2 was not changed: use same contribution as LPF package
          RHO2ADJ = RHO2
      END IF
      IF(ISSINDEX.GT.0) THEN
          ! SC1 was changed: TVM removes old storage volume and adds new storage volume
          RHO1ADJ = TLED * SOLD * SSPREV * HOLD(N)
     1               - TLED * SNEW * SSCUR * HNEW(N)
          STRGREM = STRGREM + RHO1
          STRGADD = STRGADD + RHO1ADJ
      ELSE
          ! SC1 was not changed: use same contribution as LPF package
          RHO1ADJ = RHO1
      END IF
      STRG = RHO1 + RHO2
      STRGADJ = RHO1ADJ + RHO2ADJ
C
C------Reverse LPF's BUFF and FLOWJA contributions and re-apply with new values
  288 BUFF(N) = STRGADJ
      FLOWJA(IA(N)) = -STRGADJ
      IF(STRG.LT.0.) THEN
          STOUT = STOUT + STRGREM
      ELSE
          STOIN = STOIN - STRGREM
      END IF
      IF(STRGADJ.LT.0.) THEN
          STOUT = STOUT - STRGADD
      ELSE
          STOIN = STOIN + STRGADD
      END IF
C
  300 CONTINUE
C----------------------------------------------------------------------------
C
      RETURN
      END
C
C
      REAL FUNCTION TVMU2INTERP(STARTVALUE, ENDVALUE, KPER, TVMLOGBASE,
     1                          TVMLOGBASEINV)
C     ******************************************************************
C     INTERPOLATE BETWEEN STARTVALUE AND ENDVALUE AT PERTIM IN KPER
C     ******************************************************************
      USE GWFBASMODULE,ONLY:PERTIM
      USE GLOBAL,ONLY:PERLEN
      USE TVMU2MODULE
C
      REAL :: VALUEA, VALUEB
C
      IF(PERLEN(KPER).LE.0) THEN
          TVMU2INTERP = ENDVALUE
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

          TVMU2INTERP = VALUEA + FRAC * (VALUEB - VALUEA)

          IF(TVMLOGBASE.GT.0) THEN
              TVMU2INTERP = TVMLOGBASE**TVMU2INTERP
          END IF
      END IF
C
      RETURN
      END
C
C
      SUBROUTINE TVMU2DASP
C     ******************************************************************
C     DEALLOCATE TVM STRESS PERIOD ARRAYS
C     ******************************************************************
      USE TVMU2MODULE
C
      IF(NTVMSS.GT.0.OR.NTVMSY.GT.0) THEN
          DEALLOCATE(ITVMSMAPSS)
          DEALLOCATE(ITVMSMAPSY)
          NTVMSCHANGES = 0
      END IF
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
      IF(NTVMPOR.GT.0) THEN
          DEALLOCATE(ITVMPOR)
          DEALLOCATE(PORNEWA)
          DEALLOCATE(PORNEWB)
          NTVMPOR = 0
      END IF
C
      RETURN
      END
C
C
      SUBROUTINE TVMU2DA
C     ******************************************************************
C     DEALLOCATE TVM ARRAYS
C     ******************************************************************
      USE TVMU2MODULE
C-------Deallocate per-SP arrays
      CALL TVMU2DASP
C-------Deallocate saved PGF array
      DEALLOCATE(INITPGF)
C-------Deallocate previous-SP storage value arrays if a step function was used
      DEALLOCATE(SSOLD)
      DEALLOCATE(SYOLD)
      IF(ALLOCATED(POROLD)) DEALLOCATE(POROLD)
C
      RETURN
      END

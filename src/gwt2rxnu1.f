C ---------------------------------------------------------------------
C     GENERALIZED MULTI-SPECIES KINETIC REACTION MODULE (REACTION1)
C     ACTIVATED BY THE 'REACTION1' KEYWORD IN THE BCT PACKAGE.  READS A
C     PER-REACTION LUMPED KINETIC CONSTANT, THE PARTICIPATING SPECIES
C     (REACTANTS TO THE LEFT AND PRODUCTS TO THE RIGHT OF THE REACTION
C     EQUATION), THEIR STOICHIOMETRIC COEFFICIENTS AND PER-SPECIES
C     MOLECULAR WEIGHTS.  THE REACTION IS APPLIED AS AN OPERATOR-SPLIT
C     STEP AFTER TRANSPORT, TO BOTH THE MOBILE (BCT) AND THE IMMOBILE
C     DUAL-POROSITY (DPT) DOMAIN CONCENTRATIONS.
C ---------------------------------------------------------------------
      SUBROUTINE GWT2RXN1AR(IN)
C     ******************************************************************
C     ALLOCATE ARRAYS AND READ DATA FOR THE KINETIC REACTION MODULE
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,      ONLY:IOUT,NEQS,NODES,IDPT,IFREFM, IUNIT
      USE GWTBCTMODULE,ONLY:NTCOMP
      USE GWTRXNMODULE,ONLY:IRXNCB,NRXN,NRXPART,AKRXN,NPARTL,NPARTR,
     1  IPARTL,IPARTR,ASTOL,ASTOR,AMWCOMP,RXNQQ,RXNQQIM,RXNCONCO,
     2  RXNCONCOIM
C
      CHARACTER*24 ANAME(7)
      CHARACTER*300 LINE
      INTEGER I,IR,IC,K,LLOC,ISTART,ISTOP
      REAL R
      DATA ANAME(2)/'AKRXN - KINETIC CONST '/
      DATA ANAME(3)/'REACTANT COMP INDEX   '/
      DATA ANAME(4)/'REACT STOICHIOMETRY L '/
      DATA ANAME(5)/'PRODUCT COMP INDEX    '/
      DATA ANAME(6)/'PROD STOICHIOMETRY R  '/
      DATA ANAME(7)/'AMWCOMP - MOL WEIGHT  '/
C     ------------------------------------------------------------------
C1------ALLOCATE SCALARS.
      ALLOCATE(IRXNCB,NRXN,NRXPART)
      WRITE(IOUT,1)IN
    1 FORMAT(1X,/1X,'REACTION1 -- GENERALIZED KINETIC REACTION MODULE,',
     1  1X,'INPUT READ FROM BCT UNIT',I4)
C
C1A-----READ IRXNCB (C-B-C FLUX FLAG OR UNIT NUMBER) AND NRXN.
      READ(IN,'(A)') LINE
      IF(IFREFM.EQ.0) THEN
        READ(LINE,'(I10,I10)') IRXNCB,NRXN
      ELSE
        LLOC=1
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IRXNCB,R,IOUT,IN)
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NRXN,R,IOUT,IN)
      ENDIF
C
      IF(NRXN.LT.1) THEN
        WRITE(IOUT,*)' NRXN MUST BE AT LEAST 1 FOR REACTION1; STOPPING'
        CALL USTOP(' ')
      ENDIF
      WRITE(IOUT,2) IRXNCB,NRXN
    2 FORMAT(1X,'C-B-C REACTION FLUX FLAG OR UNIT (IRXNCB)     =',I5
     1  /1X,'NUMBER OF REACTION EQUATIONS (NRXN)           =',I5)
C
C2------ALLOCATE REACTION ARRAYS.  A SIDE CANNOT HAVE MORE PARTICIPANTS
C2------THAN THE TOTAL NUMBER OF COMPONENTS, SO NRXPART =< NTCOMP.
      NRXPART = NTCOMP
      ALLOCATE(AKRXN(NRXN))
      ALLOCATE(NPARTL(NRXN),NPARTR(NRXN))
      ALLOCATE(IPARTL(NRXN,NRXPART),IPARTR(NRXN,NRXPART))
      ALLOCATE(ASTOL(NRXN,NRXPART),ASTOR(NRXN,NRXPART))
      ALLOCATE(AMWCOMP(NTCOMP))
      AKRXN = 0.0D0
      NPARTL = 0
      NPARTR = 0
      IPARTL = 0
      IPARTR = 0
      ASTOL = 0.0
      ASTOR = 0.0
      AMWCOMP = 0.0
      K = 0
C
C3------READ EACH REACTION EQUATION (LEFT REACTANTS, RIGHT PRODUCTS).
      DO IR = 1,NRXN
        READ(IN,*) AKRXN(IR)
        IF(AKRXN(IR).LE.0.0D0) THEN
          WRITE(IOUT,*)' AKRXN MUST BE > 0 FOR REACTION',IR
          CALL USTOP(' ')
        ENDIF
C ......LEFT SIDE OF EQUATION       
        READ(IN,*) NPARTL(IR)
        IF(NPARTL(IR).LT.1.OR.NPARTL(IR).GT.NRXPART) THEN
          WRITE(IOUT,*)' INVALID NUMBER OF LEFT SPECIES FOR REACTION',IR
          CALL USTOP(' ')
        ENDIF
        CALL U1DINT(IPARTL(IR,1),ANAME(3),NPARTL(IR),K,IN,IOUT)
        CALL U1DREL(ASTOL(IR,1),ANAME(4),NPARTL(IR),K,IN,IOUT)
C ......RIGHT SIDE OF EQUATION         
        READ(IN,*) NPARTR(IR)
        IF(NPARTR(IR).LT.1.OR.NPARTR(IR).GT.NRXPART) THEN
          WRITE(IOUT,*)' INVALID NUMBER OF RIGHT SPECIES FOR REACTION',
     1      IR
          CALL USTOP(' ')
        ENDIF
        CALL U1DINT(IPARTR(IR,1),ANAME(5),NPARTR(IR),K,IN,IOUT)
        CALL U1DREL(ASTOR(IR,1),ANAME(6),NPARTR(IR),K,IN,IOUT)
C
C3A-------ECHO AND VALIDATE THE REACTION EQUATION.
        WRITE(IOUT,3) IR
    3   FORMAT(1X,/1X,'REACTION EQUATION NUMBER',I4)
        WRITE(IOUT,10) AKRXN(IR)
   10   FORMAT(1X,'  KINETIC CONSTANT (AKRXN)            =',1PG15.6)
        WRITE(IOUT,4) NPARTL(IR)
    4   FORMAT(1X,'  NUMBER OF REACTANT (LEFT) SPECIES  =',I5)
C ......LEFT SIDE OF EQUATION     
        DO I = 1,NPARTL(IR)
          IC = IPARTL(IR,I)
          IF(IC.LT.1.OR.IC.GT.NTCOMP) THEN
            WRITE(IOUT,*)' REACTANT COMPONENT INDEX OUT OF RANGE:',IC
            CALL USTOP(' ')
          ENDIF
          IF(ASTOL(IR,I).LE.0.0) THEN
            WRITE(IOUT,*)' REACTANT STOICHIOMETRY MUST BE > 0 FOR COMP',
     1        IC
            CALL USTOP(' ')
          ENDIF
          WRITE(IOUT,5) IC,ASTOL(IR,I)
    5     FORMAT(1X,'    REACTANT COMPONENT',I5,'   COEFFICIENT',G15.6)
        ENDDO
C ......RIGHT SIDE OF EQUATION 
        WRITE(IOUT,6) NPARTR(IR)
    6   FORMAT(1X,'  NUMBER OF PRODUCT (RIGHT) SPECIES  =',I5)
        DO I = 1,NPARTR(IR)
          IC = IPARTR(IR,I)
          IF(IC.LT.1.OR.IC.GT.NTCOMP) THEN
            WRITE(IOUT,*)' PRODUCT COMPONENT INDEX OUT OF RANGE:',IC
            CALL USTOP(' ')
          ENDIF
          IF(ASTOR(IR,I).LE.0.0) THEN
            WRITE(IOUT,*)' PRODUCT STOICHIOMETRY MUST BE > 0 FOR COMP',
     1        IC
            CALL USTOP(' ')
          ENDIF
          WRITE(IOUT,7) IC,ASTOR(IR,I)
    7     FORMAT(1X,'    PRODUCT COMPONENT ',I5,'   COEFFICIENT',G15.6)
        ENDDO
      ENDDO
C
C4------READ MOLECULAR WEIGHT [KG/MOL] FOR EACH COMPONENT.
      CALL U1DREL(AMWCOMP,ANAME(7),NTCOMP,K,IN,IOUT)
      WRITE(IOUT,8)
    8 FORMAT(1X,/1X,'MOLECULAR WEIGHT (AMWCOMP) BY COMPONENT:')
      DO I = 1,NTCOMP
        WRITE(IOUT,9) I,AMWCOMP(I)
    9   FORMAT(1X,'    COMPONENT',I5,'   AMWCOMP =',G15.6)
      ENDDO
C
C4A-----EVERY PARTICIPATING COMPONENT NEEDS A POSITIVE MOLECULAR WEIGHT
C4A-----SO MASS <-> MOLE CONVERSION IS WELL DEFINED.
      DO IR = 1,NRXN
        DO I = 1,NPARTL(IR)
          IF(AMWCOMP(IPARTL(IR,I)).LE.0.0) THEN
            WRITE(IOUT,*)' AMWCOMP MUST BE > 0 FOR REACTANT COMPONENT',
     1        IPARTL(IR,I)
            CALL USTOP(' ')
          ENDIF
        ENDDO
        DO I = 1,NPARTR(IR)
          IF(AMWCOMP(IPARTR(IR,I)).LE.0.0) THEN
            WRITE(IOUT,*)' AMWCOMP MUST BE > 0 FOR PRODUCT COMPONENT',
     1        IPARTR(IR,I)
            CALL USTOP(' ')
          ENDIF
        ENDDO
      ENDDO
C
C5------ALLOCATE PER-NODE REACTION MASS-CHANGE ACCUMULATORS FOR BUDGET.
      ALLOCATE(RXNQQ(NEQS,NTCOMP))
      RXNQQ = 0.0D0
      ALLOCATE(RXNCONCO(NEQS,NTCOMP))
      RXNCONCO = 0.0D0
C5A-----FOR DUAL POROSITY TRANSPORT      
      IF(IUNIT(30).GT.0) THEN
        ALLOCATE(RXNQQIM(NODES,NTCOMP))
        RXNQQIM = 0.0D0
        ALLOCATE(RXNCONCOIM(NODES,NTCOMP))
        RXNCONCOIM = 0.0D0
      ENDIF
C
C6------RETURN
      RETURN
      END
C ---------------------------------------------------------------------
      SUBROUTINE GWT2RXN1FM(KSTP,KPER)
C     ******************************************************************
C     APPLY THE KINETIC REACTION TO THE MOBILE-DOMAIN CONCENTRATIONS
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,      ONLY:NODES,NEQS,AREA,TOP,BOT,Sn
      USE CLN1MODULE,  ONLY:ACLNNDS
      USE GWFBASMODULE,ONLY:DELT
      USE GWTBCTMODULE,ONLY:CONC,ICBUND,NTCOMP,PRSITY,MCOMPT
      USE GWTRXNMODULE,ONLY:RXNQQ
C
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: CWRK,CWRKO,FVOL
      DOUBLE PRECISION ALENG,VOLU,VFAC,FWAT
C     ------------------------------------------------------------------
      ALLOCATE(CWRK(NTCOMP),CWRKO(NTCOMP),FVOL(NTCOMP))
C
C1------RESET THE PER-NODE MASS-CHANGE ACCUMULATOR.
      RXNQQ = 0.0D0
C
C2------LOOP OVER ALL TRANSPORT NODES.
      DO N = 1,NEQS
        IF(ICBUND(N).EQ.0) CYCLE
C2A-------CELL VOLUME (GWF OR CLN NODE).
        IF(N.LE.NODES) THEN
          ALENG = TOP(N) - BOT(N)
        ELSE
          ALENG = ACLNNDS(N-NODES,4)
        ENDIF
        VOLU = AREA(N) * ALENG
C2B-------VOLUME FRACTION EACH SPECIES IS STORED PER.  MOBILE SPECIES
C2B-------ARE PER UNIT WATER VOLUME, IMCOMP SPECIES (IC > MCOMPT) ARE PER
C2B-------UNIT TOTAL CELL VOLUME.  THE REACTION EXTENT IS IN MOLES PER
C2B-------UNIT CELL VOLUME, SO FVOL CONVERTS IT TO EACH BASIS.
        IF(N.LE.NODES) THEN
          FWAT = PRSITY(N) * Sn(N)
        ELSE
          FWAT = 1.0D0
        ENDIF
        IF(FWAT.LE.0.0D0) CYCLE
        DO IC = 1,NTCOMP
          IF(IC.GT.MCOMPT) THEN
            FVOL(IC) = 1.0D0
          ELSE
            FVOL(IC) = FWAT
          ENDIF
        ENDDO
C2C-------GATHER CONCENTRATIONS, REACT, THEN STORE MASS CHANGE.
        DO IC = 1,NTCOMP
          CWRK(IC)  = CONC(N,IC)
          CWRKO(IC) = CONC(N,IC)
        ENDDO
        CALL SGWT2RXN1REACT(CWRK,NTCOMP,DELT,FVOL)
        DO IC = 1,NTCOMP
          CONC(N,IC) = CWRK(IC)
          VFAC = VOLU * FVOL(IC)
          RXNQQ(N,IC) = (CWRK(IC) - CWRKO(IC)) * VFAC
        ENDDO
      ENDDO
C
      DEALLOCATE(CWRK,CWRKO,FVOL)
C3------RETURN
      RETURN
      END
C ---------------------------------------------------------------------
      SUBROUTINE GWT2RXN1FMIM(KSTP,KPER)
C     ******************************************************************
C     APPLY THE KINETIC REACTION TO THE IMMOBILE (DUAL-POROSITY) DOMAIN
C     CONCENTRATIONS
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,      ONLY:NODES,AREA,TOP,BOT
      USE GWFBASMODULE,ONLY:DELT
      USE GWFDPFMODULE,ONLY:SnIM
      USE GWTBCTMODULE,ONLY:NTCOMP,MCOMPT
      USE GWTDPTMODULE,ONLY:CONCIM,ICBUNDIM,PRSITYIM
      USE GWTRXNMODULE,ONLY:RXNQQIM
C
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: CWRK,CWRKO,FVOL
      DOUBLE PRECISION ALENG,VOLU,VFAC,FWAT
C     ------------------------------------------------------------------
      ALLOCATE(CWRK(NTCOMP),CWRKO(NTCOMP),FVOL(NTCOMP))
C
C1------RESET THE PER-NODE MASS-CHANGE ACCUMULATOR.
      RXNQQIM = 0.0D0
C
C2------LOOP OVER ALL POROUS-MATRIX NODES (NO CLN IN IMMOBILE DOMAIN).
      DO N = 1,NODES
        IF(ICBUNDIM(N).EQ.0) CYCLE
        ALENG = TOP(N) - BOT(N)
        VOLU = AREA(N) * ALENG
C2A-------VOLUME FRACTION EACH SPECIES IS STORED PER (SEE GWT2RXN1FM).
        FWAT = PRSITYIM(N) * SnIM(N)
        IF(FWAT.LE.0.0D0) CYCLE
        DO IC = 1,NTCOMP
          IF(IC.GT.MCOMPT) THEN
            FVOL(IC) = 1.0D0
          ELSE
            FVOL(IC) = FWAT
          ENDIF
        ENDDO
        DO IC = 1,NTCOMP
          CWRK(IC)  = CONCIM(N,IC)
          CWRKO(IC) = CONCIM(N,IC)
        ENDDO
        CALL SGWT2RXN1REACT(CWRK,NTCOMP,DELT,FVOL)
        DO IC = 1,NTCOMP
          CONCIM(N,IC) = CWRK(IC)
          VFAC = VOLU * FVOL(IC)
          RXNQQIM(N,IC) = (CWRK(IC) - CWRKO(IC)) * VFAC
        ENDDO
      ENDDO
C
      DEALLOCATE(CWRK,CWRKO,FVOL)
C3------RETURN
      RETURN
      END
C ---------------------------------------------------------------------
      SUBROUTINE SGWT2RXN1REACT(CMASS,NC,DT,FVOL)
C     ******************************************************************
C     ADVANCE ALL REACTION EQUATIONS OVER ONE TIME STEP AT A SINGLE
C     CELL.  CMASS HOLDS COMPONENT CONCENTRATIONS ON A MASS BASIS AND
C     IS UPDATED IN PLACE.  CONVERSION TO/FROM MOLAR UNITS USES AMWCOMP.
C     FVOL(IC) IS THE FRACTION OF THE CELL VOLUME THAT COMPONENT IC IS
C     STORED PER (WATER-FILLED FRACTION FOR MOBILE SPECIES, 1 FOR
C     SPECIES CARRIED PER UNIT TOTAL VOLUME).  THE REACTION EXTENT IS
C     IN MOLES PER UNIT TOTAL CELL VOLUME SO THAT THE MASS CONSUMED AND
C     GENERATED OBEYS STOICHIOMETRY ACROSS THE DIFFERENT BASES.
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GWTRXNMODULE,ONLY:NRXN,NPARTL,NPARTR,IPARTL,IPARTR,
     1  ASTOL,ASTOR,AMWCOMP
C
      INTEGER NC,IR,I,IC
      DOUBLE PRECISION CMASS(NC),DT,FVOL(NC)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: CMOL
C     ------------------------------------------------------------------
      ALLOCATE(CMOL(NC))
C
C1------CONVERT MASS -> MOLAR (AMWCOMP <= 0 MEANS TREAT AS MOLAR ALREADY).
      DO IC = 1,NC
        IF(AMWCOMP(IC).GT.0.0) THEN
          CMOL(IC) = CMASS(IC) / AMWCOMP(IC)
        ELSE
          CMOL(IC) = CMASS(IC)
        ENDIF
      ENDDO
C
C2------APPLY EACH REACTION SEQUENTIALLY (OPERATOR SPLIT).
      DO IR = 1,NRXN
        IF(NPARTL(IR).EQ.2) THEN
C2A-------TWO REACTANTS: USE ANALYTICAL BILINEAR SOLUTION.
          CALL SGWT2RXN1EXACT(IR,CMOL,NC,DT,FVOL)
        ELSE
C2B-------GENERAL CASE: EXPLICIT SUB-STEPPED INTEGRATION.
          CALL SGWT2RXN1RK(IR,CMOL,NC,DT,FVOL)
        ENDIF
      ENDDO
C
C3------CONVERT MOLAR -> MASS.
      DO IC = 1,NC
        IF(AMWCOMP(IC).GT.0.0) THEN
          CMASS(IC) = CMOL(IC) * AMWCOMP(IC)
        ELSE
          CMASS(IC) = CMOL(IC)
        ENDIF
      ENDDO
C
      DEALLOCATE(CMOL)
C4------RETURN
      RETURN
      END
C ---------------------------------------------------------------------
      SUBROUTINE SGWT2RXN1EXACT(IR,CMOL,NC,DT,FVOL)
C     ******************************************************************
C     ANALYTICAL UPDATE FOR A BILINEAR TWO-REACTANT REACTION
C         R = AKRXN(IR) * X1 * X2 ,   DX1/DT = -S1*R ,  DX2/DT = -S2*R
C     WITH INVARIANT J = S1*X2 - S2*X1 (SEE AU_EXTRACTION_REACTIONS).
C     PRODUCTS ARE INCREMENTED BY THEIR STOICHIOMETRY TIMES THE EXTENT.
C     STOICHIOMETRY IS DIVIDED BY FVOL SO THAT THE EXTENT IS IN MOLES
C     PER UNIT TOTAL CELL VOLUME FOR EVERY PARTICIPANT.
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GWTRXNMODULE,ONLY:AKRXN,NPARTL,NPARTR,IPARTL,IPARTR,
     1  ASTOL,ASTOR
C
      INTEGER IR,NC,I1,I2,IC,I
      DOUBLE PRECISION CMOL(NC),DT,FVOL(NC)
      DOUBLE PRECISION S1,S2,X1,X2,X1N,BIGJ,E,DENOM,XI,EPS
C     ------------------------------------------------------------------
      EPS = 1.0D-20
      I1 = IPARTL(IR,1)
      I2 = IPARTL(IR,2)
      S1 = ASTOL(IR,1) / FVOL(I1)
      S2 = ASTOL(IR,2) / FVOL(I2)
      X1 = CMOL(I1)
      X2 = CMOL(I2)
C
C1------NO REACTION IF EITHER REACTANT IS ABSENT.
      IF(X1.LE.0.0D0.OR.X2.LE.0.0D0) RETURN
C
C2------INVARIANT AND ANALYTICAL SOLUTION FOR THE PRIMARY REACTANT X1.
      BIGJ = S1*X2 - S2*X1
      IF(ABS(BIGJ).GT.EPS) THEN
        E = EXP(-AKRXN(IR)*BIGJ*DT)
        DENOM = S1*X2 - S2*X1*E
        IF(ABS(DENOM).LT.EPS) THEN
          X1N = X1
        ELSE
          X1N = BIGJ*X1*E / DENOM
        ENDIF
      ELSE
C2A-------DEGENERATE CASE J ~ 0: DX1/DT = -AKRXN(IR)*S2*X1**2.
        X1N = X1 / (1.0D0 + AKRXN(IR)*S2*X1*DT)
      ENDIF
      IF(X1N.LT.0.0D0) X1N = 0.0D0
C
C3------REACTION EXTENT (MOLES PER UNIT S1) AND CLAMP TO AVAILABLE X2.
      XI = (X1 - X1N) / S1
      IF(XI.LT.0.0D0) XI = 0.0D0
      IF(S2*XI.GT.X2) XI = X2 / S2
C
C4------UPDATE REACTANTS.
      CMOL(I1) = X1 - S1*XI
      CMOL(I2) = X2 - S2*XI
      IF(CMOL(I1).LT.0.0D0) CMOL(I1) = 0.0D0
      IF(CMOL(I2).LT.0.0D0) CMOL(I2) = 0.0D0
C
C5------UPDATE PRODUCTS.
      DO I = 1,NPARTR(IR)
        IC = IPARTR(IR,I)
        CMOL(IC) = CMOL(IC) + ASTOR(IR,I)/FVOL(IC)*XI
      ENDDO
C
C6------RETURN
      RETURN
      END
C ---------------------------------------------------------------------
      SUBROUTINE SGWT2RXN1RK(IR,CMOL,NC,DT,FVOL)
C     ******************************************************************
C     GENERAL N-REACTANT INTEGRATION (HEUN WITH INTERNAL SUB-STEPS) FOR
C     THE RATE LAW  R = AKRXN(IR) * PRODUCT(X_LEFT),  DX_I/DT = -S_I * R.
C     USED WHEN A REACTION HAS OTHER THAN TWO REACTANTS, OR WHEN THE
C     CALLER FORCES THE NUMERICAL PATH FOR A TWO-REACTANT EQUATION.
C     SLE AND SRE ARE THE STOICHIOMETRIES DIVIDED BY FVOL, SO THE EXTENT
C     IS IN MOLES PER UNIT TOTAL CELL VOLUME FOR EVERY PARTICIPANT.
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GWTRXNMODULE,ONLY:AKRXN,NPARTL,NPARTR,IPARTL,IPARTR,
     1  ASTOL,ASTOR
C
      INTEGER IR,NC,ISUB,NSUB,I,IC
      DOUBLE PRECISION CMOL(NC),DT,FVOL(NC)
      DOUBLE PRECISION HSUB,RRATE,RRATE2,XIMAX,XI,REL,RELMAX,EPS
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: CMOLP,SLE,SRE
C     ------------------------------------------------------------------
      EPS = 1.0D-20
C
C1------NO REACTION IF ANY REACTANT IS ABSENT (MATCHES SGWT2RXN1EXACT).
      DO I = 1,NPARTL(IR)
        IF(CMOL(IPARTL(IR,I)).LE.0.0D0) RETURN
      ENDDO
C
C2------EFFECTIVE STOICHIOMETRY ON THE TOTAL-CELL-VOLUME EXTENT BASIS.
      ALLOCATE(SLE(NPARTL(IR)),SRE(NPARTR(IR)))
      DO I = 1,NPARTL(IR)
        SLE(I) = ASTOL(IR,I) / FVOL(IPARTL(IR,I))
      ENDDO
      DO I = 1,NPARTR(IR)
        SRE(I) = ASTOR(IR,I) / FVOL(IPARTR(IR,I))
      ENDDO
C
C3------SUB-STEPS FROM MAX FRACTIONAL REACTANT CHANGE OVER DT.
C3------(RRATE ALREADY INCLUDES AKRXN(IR); DO NOT MULTIPLY BY AKRXN AGAIN.)
      CALL SGWT2RXN1RATE(IR,CMOL,NC,RRATE)
      RELMAX = 0.0D0
      IF(RRATE.GT.0.0D0) THEN
        DO I = 1,NPARTL(IR)
          IC = IPARTL(IR,I)
          IF(CMOL(IC).GT.EPS .AND. SLE(I).GT.0.0D0) THEN
            REL = SLE(I)*RRATE*DT / CMOL(IC)
            IF(REL.GT.RELMAX) RELMAX = REL
          ENDIF
        ENDDO
      ENDIF
      NSUB = 1
      IF(RELMAX.GT.0.0D0) NSUB = 1 + INT(10.0D0*RELMAX)
      IF(NSUB.LT.1) NSUB = 1
      IF(NSUB.GT.1000) NSUB = 1000
      HSUB = DT / DBLE(NSUB)
      ALLOCATE(CMOLP(NC))
C
C4------HEUN (PREDICTOR-CORRECTOR) SUB-STEPS WITH EXTENT CLAMPING.
      DO ISUB = 1,NSUB
        CALL SGWT2RXN1RATE(IR,CMOL,NC,RRATE)
        IF(RRATE.LE.0.0D0) GO TO 20
C4A-------PREDICTOR EXTENT AND STATE.
        XI = RRATE * HSUB
        XIMAX = 1.0D30
        DO I = 1,NPARTL(IR)
          IC = IPARTL(IR,I)
          IF(SLE(I).GT.0.0D0) THEN
            IF(CMOL(IC)/SLE(I).LT.XIMAX) THEN
              XIMAX = CMOL(IC)/SLE(I)
            ENDIF
          ENDIF
        ENDDO
        IF(XI.GT.XIMAX) XI = XIMAX
        IF(XI.LT.0.0D0) XI = 0.0D0
        DO IC = 1,NC
          CMOLP(IC) = CMOL(IC)
        ENDDO
        DO I = 1,NPARTL(IR)
          IC = IPARTL(IR,I)
          CMOLP(IC) = CMOL(IC) - SLE(I)*XI
          IF(CMOLP(IC).LT.0.0D0) CMOLP(IC) = 0.0D0
        ENDDO
        DO I = 1,NPARTR(IR)
          IC = IPARTR(IR,I)
          CMOLP(IC) = CMOLP(IC) + SRE(I)*XI
        ENDDO
C4B-------CORRECTOR: AVERAGE RATE AT START AND PREDICTED STATE.
        CALL SGWT2RXN1RATE(IR,CMOLP,NC,RRATE2)
        XI = 0.5D0*(RRATE + RRATE2)*HSUB
        XIMAX = 1.0D30
        DO I = 1,NPARTL(IR)
          IC = IPARTL(IR,I)
          IF(SLE(I).GT.0.0D0) THEN
            IF(CMOL(IC)/SLE(I).LT.XIMAX) THEN
              XIMAX = CMOL(IC)/SLE(I)
            ENDIF
          ENDIF
        ENDDO
        IF(XI.GT.XIMAX) XI = XIMAX
        IF(XI.LT.0.0D0) XI = 0.0D0
        DO I = 1,NPARTL(IR)
          IC = IPARTL(IR,I)
          CMOL(IC) = CMOL(IC) - SLE(I)*XI
          IF(CMOL(IC).LT.0.0D0) CMOL(IC) = 0.0D0
        ENDDO
        DO I = 1,NPARTR(IR)
          IC = IPARTR(IR,I)
          CMOL(IC) = CMOL(IC) + SRE(I)*XI
        ENDDO
      ENDDO
20    CONTINUE
C
      DEALLOCATE(CMOLP,SLE,SRE)
C5------RETURN
      RETURN
      END
C ---------------------------------------------------------------------
      SUBROUTINE SGWT2RXN1RATE(IR,CMOL,NC,RRATE)
C     ******************************************************************
C     EVALUATE THE BILINEAR (FIRST-ORDER IN EACH REACTANT) RATE
C         R = AKRXN(IR) * PRODUCT(X_LEFT)
C     ANY NON-POSITIVE REACTANT FORCES R = 0 (NO REACTION).
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GWTRXNMODULE,ONLY:AKRXN,NPARTL,IPARTL
C
      INTEGER IR,NC,I,IC
      DOUBLE PRECISION CMOL(NC),RRATE
C     ------------------------------------------------------------------
      RRATE = AKRXN(IR)
      DO I = 1,NPARTL(IR)
        IC = IPARTL(IR,I)
        IF(CMOL(IC).LE.0.0D0) THEN
          RRATE = 0.0D0
          RETURN
        ENDIF
        RRATE = RRATE * CMOL(IC)
      ENDDO
      IF(RRATE.LT.0.0D0) RRATE = 0.0D0
C
      RETURN
      END
C ---------------------------------------------------------------------
      SUBROUTINE GWT2RXN1BD(KSTP,KPER,ICOMP)
C     ******************************************************************
C     MASS BALANCE OF THE KINETIC REACTION IN THE MOBILE DOMAIN.
C     BUDGET IS ACCUMULATED FOR EVERY COMPONENT 1..NTCOMP (INCLUDING
C     IMCOMP SPECIES WITH ICOMP > MCOMPT).
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,      ONLY:IOUT,NCOL,NROW,NLAY,BUFF,NODES,IUNSTR,NEQS,
     1  INCLN
      USE CLN1MODULE,  ONLY:NCLNNDS,ICLNMB
      USE GWFBASMODULE,ONLY:ICBCFL,DELT,PERTIM,TOTIM
      USE GWTBCTMODULE,ONLY:ICBUND,MSUMT,VBVLT,VBNMT,NTCOMP
      USE GWTRXNMODULE,ONLY:IREACT,IRXNCB,RXNQQ
C
      CHARACTER*16 TEXT(2)
      DOUBLE PRECISION RATIN,RATOUT,QQ,ZERO
      INTEGER IBD
      DATA TEXT(1) /'   MASS REACTION'/
      DATA TEXT(2) /'CLN MASS REACTN '/
C     ------------------------------------------------------------------
C1------RETURN IF REACTION IS OFF OR COMPONENT INDEX IS OUT OF RANGE.
      IF(IREACT.LE.0) RETURN
      IF(ICOMP.LT.1.OR.ICOMP.GT.NTCOMP) RETURN
C
C2------INITIALIZE ACCUMULATORS AND CELL-BY-CELL FLAG.
      ZERO = 0.0D0
      RATIN = ZERO
      RATOUT = ZERO
      IBD = 0
      IF(IRXNCB.LT.0 .AND. ICBCFL.NE.0) IBD = -1
      IF(IRXNCB.GT.0) IBD = ICBCFL
C
C3------CLEAR THE BUFFER.
      DO N = 1,NEQS
        BUFF(N) = ZERO
      ENDDO
C
C4------ACCUMULATE THE REACTION MASS RATE FOR THIS COMPONENT.
      DO N = 1,NEQS
        IF(ICBUND(N).EQ.0) CYCLE
        QQ = RXNQQ(N,ICOMP) / DELT
        BUFF(N) = BUFF(N) + QQ
        IF(QQ.GE.ZERO) THEN
          RATIN = RATIN + QQ
        ELSE
          RATOUT = RATOUT - QQ
        ENDIF
      ENDDO
C
C5------SAVE CELL-BY-CELL FLOWS IF REQUESTED.
      IF(IBD.GE.1) THEN
        IF(IUNSTR.EQ.0) THEN
          CALL UBUDSV(KSTP,KPER,TEXT(1),IRXNCB,BUFF,NCOL,NROW,
     1                          NLAY,IOUT)
        ELSE
          CALL UBUDSVU(KSTP,KPER,TEXT(1),IRXNCB,BUFF,NODES,
     1                          IOUT,PERTIM,TOTIM)
        ENDIF
        IF(INCLN.GT.0.AND.ICLNMB.NE.0) THEN
          CALL UBUDSVU(KSTP,KPER,TEXT(2),ICLNMB,BUFF(NODES+1:NEQS),
     1                 NCLNNDS,IOUT,PERTIM,TOTIM)
        ENDIF
      ENDIF
C
C6------MOVE RATES, VOLUMES AND LABELS INTO THE PRINT ARRAYS.
      VBVLT(3,MSUMT,ICOMP) = RATIN
      VBVLT(4,MSUMT,ICOMP) = RATOUT
      VBVLT(1,MSUMT,ICOMP) = VBVLT(1,MSUMT,ICOMP) + RATIN*DELT
      VBVLT(2,MSUMT,ICOMP) = VBVLT(2,MSUMT,ICOMP) + RATOUT*DELT
      VBNMT(MSUMT,ICOMP) = TEXT(1)
      MSUMT = MSUMT + 1
C
C7------RETURN
      RETURN
      END
C ---------------------------------------------------------------------
      SUBROUTINE GWT2RXN1BDIM(KSTP,KPER,ICOMP)
C     ******************************************************************
C     MASS BALANCE OF THE KINETIC REACTION IN THE IMMOBILE DOMAIN.
C     BUDGET IS ACCUMULATED FOR EVERY COMPONENT 1..NTCOMP (INCLUDING
C     IMCOMP SPECIES WITH ICOMP > MCOMPT).
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,      ONLY:IOUT,NCOL,NROW,NLAY,BUFF,NODES,IUNSTR
      USE GWFBASMODULE,ONLY:ICBCFL,DELT,PERTIM,TOTIM
      USE GWTBCTMODULE,ONLY:MSUMT,VBVLT,VBNMT,NTCOMP
      USE GWTDPTMODULE,ONLY:ICBUNDIM,IDPTCB
      USE GWTRXNMODULE,ONLY:IREACT,RXNQQIM
C
      CHARACTER*16 TEXT
      DOUBLE PRECISION RATIN,RATOUT,QQ,ZERO
      INTEGER IBD
      DATA TEXT /'IM MASS REACTION'/
C     ------------------------------------------------------------------
C1------RETURN IF REACTION IS OFF OR COMPONENT INDEX IS OUT OF RANGE.
      IF(IREACT.LE.0) RETURN
      IF(ICOMP.LT.1.OR.ICOMP.GT.NTCOMP) RETURN
C
C2------INITIALIZE ACCUMULATORS AND CELL-BY-CELL FLAG.
      ZERO = 0.0D0
      RATIN = ZERO
      RATOUT = ZERO
      IBD = 0
      IF(IDPTCB.LT.0 .AND. ICBCFL.NE.0) IBD = -1
      IF(IDPTCB.GT.0) IBD = ICBCFL
C
C3------CLEAR THE BUFFER (POROUS MATRIX NODES ONLY).
      DO N = 1,NODES
        BUFF(N) = ZERO
      ENDDO
C
C4------ACCUMULATE THE REACTION MASS RATE FOR THIS COMPONENT.
      DO N = 1,NODES
        IF(ICBUNDIM(N).EQ.0) CYCLE
        QQ = RXNQQIM(N,ICOMP) / DELT
        BUFF(N) = BUFF(N) + QQ
        IF(QQ.GE.ZERO) THEN
          RATIN = RATIN + QQ
        ELSE
          RATOUT = RATOUT - QQ
        ENDIF
      ENDDO
C
C5------SAVE CELL-BY-CELL FLOWS IF REQUESTED.
      IF(IBD.GE.1) THEN
        IF(IUNSTR.EQ.0) THEN
          CALL UBUDSV(KSTP,KPER,TEXT,IDPTCB,BUFF,NCOL,NROW,NLAY,IOUT)
        ELSE
          CALL UBUDSVU(KSTP,KPER,TEXT,IDPTCB,BUFF,NODES,
     1                          IOUT,PERTIM,TOTIM)
        ENDIF
      ENDIF
C
C6------MOVE RATES, VOLUMES AND LABELS INTO THE PRINT ARRAYS.
      VBVLT(3,MSUMT,ICOMP) = RATIN
      VBVLT(4,MSUMT,ICOMP) = RATOUT
      VBVLT(1,MSUMT,ICOMP) = VBVLT(1,MSUMT,ICOMP) + RATIN*DELT
      VBVLT(2,MSUMT,ICOMP) = VBVLT(2,MSUMT,ICOMP) + RATOUT*DELT
      VBNMT(MSUMT,ICOMP) = TEXT
      MSUMT = MSUMT + 1
C
C7------RETURN
      RETURN
      END
C ---------------------------------------------------------------------
      SUBROUTINE SGWT2RXN1PART(ICOMP,LPART)
C     ******************************************************************
C     RETURN LPART=1 IF COMPONENT ICOMP APPEARS IN ANY REACTION EQUATION
C     (AS EITHER A REACTANT OR A PRODUCT), OTHERWISE LPART=0.
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GWTRXNMODULE,ONLY:NRXN,NPARTL,NPARTR,IPARTL,IPARTR
C
      INTEGER ICOMP,LPART,IR,I
C     ------------------------------------------------------------------
      LPART = 0
      DO IR = 1,NRXN
        DO I = 1,NPARTL(IR)
          IF(IPARTL(IR,I).EQ.ICOMP) LPART = 1
        ENDDO
        DO I = 1,NPARTR(IR)
          IF(IPARTR(IR,I).EQ.ICOMP) LPART = 1
        ENDDO
      ENDDO
C
      RETURN
      END
C ---------------------------------------------------------------------
      SUBROUTINE GWT2RXN1DA
C     ******************************************************************
C     DEALLOCATE THE KINETIC REACTION MODULE ARRAYS
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GWTRXNMODULE
      INTEGER ALLOC_ERR
C     ------------------------------------------------------------------
      DEALLOCATE(NPARTL, STAT = ALLOC_ERR)
      DEALLOCATE(NPARTR, STAT = ALLOC_ERR)
      DEALLOCATE(IPARTL, STAT = ALLOC_ERR)
      DEALLOCATE(IPARTR, STAT = ALLOC_ERR)
      DEALLOCATE(ASTOL, STAT = ALLOC_ERR)
      DEALLOCATE(ASTOR, STAT = ALLOC_ERR)
      DEALLOCATE(AMWCOMP, STAT = ALLOC_ERR)
      DEALLOCATE(RXNQQ, STAT = ALLOC_ERR)
      DEALLOCATE(RXNQQIM, STAT = ALLOC_ERR)
      DEALLOCATE(RXNCONCO, STAT = ALLOC_ERR)
      DEALLOCATE(RXNCONCOIM, STAT = ALLOC_ERR)
      DEALLOCATE(IRXNCB, STAT = ALLOC_ERR)
      DEALLOCATE(NRXN, STAT = ALLOC_ERR)
      DEALLOCATE(NRXPART, STAT = ALLOC_ERR)
      DEALLOCATE(AKRXN, STAT = ALLOC_ERR)
C
      RETURN
      END

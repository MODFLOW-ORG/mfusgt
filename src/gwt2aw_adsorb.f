C -------------------------------------------------------------------------------------            
      SUBROUTINE AW_ADSORB1AL (IN)
C     ******************************************************************
C     ALLOCATE ARRAYS FOR ADSORPTION ON AIR-WATER INTERFACE AND READ CONSTANTS 
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,   ONLY:IOUT,NCOL,NROW,NLAY,IBOUND,BUFF,NODES,IUNSTR 
      USE GWFBASMODULE,ONLY:MSUM,ICBCFL,IAUXSV,DELT,PERTIM,TOTIM,DELT
      USE GWTBCTMODULE, ONLY: MCOMP
      USE AW_ADSORBMODULE  
      CHARACTER*300 LINE      
      CHARACTER*24 ANAME 
      DATA ANAME /'TABULAR AWI FN. ZONE MAP'/      
C     ------------------------------------------------------------------
C      
      ALLOCATE(AREA_AWI(NODES),AREA_AWIO(NODES))
      ALLOCATE(AK_AWI(NODES,MCOMP)) 
C ---------------------------------------------------------------------
C1A ---------ALLOCATE SPACE AND READ ZONE MAP IF EITHER K_AWI OR AW_AREA ARE TABULAR
C   --------USE SAME ZONE MAP IF EITHER (OR BOTH) ARE TABULAR
      IF(IAREA_FN.EQ.5. OR. IKAWI_FN.EQ.4) THEN 
        ITAB_AWI = 1
        ALLOCATE(NAZONES, NATABROWS)
C1B -------READ NUMBER OF ZONES AND NUMBER OF ROWS IN ZONE TABLES        
        READ(IN,*) NAZONES, NATABROWS
C1C -------DIMENSION AND READ ZONE MAP          
        ALLOCATE(IAWIZONMAP(NODES))
        CALL U1DINT(IAWIZONMAP,ANAME,NODES,0,IN,IOUT)  
      ENDIF      
C1 -------FOR THE DIFFERENT AREA COMPUTATION OPTIONS            
      IF(IAREA_FN.EQ.1) THEN 
C1A ------READ AMAX, COMPUTE FROM A = AMAX * (1-Sw)
        ALLOCATE (AWAMAX(NODES))
      ELSEIF(IAREA_FN.EQ.4) THEN 
C1B ------READ X2, X1, X0 COMPUTE FROM A = X2*Sw^2 + X1*Sw + X0
        ALLOCATE (AWAREA_X2(NODES),AWAREA_X1(NODES),AWAREA_X0(NODES))
      ELSEIF(IAREA_FN.EQ.2) THEN 
C1C ------READ GRAIN DIAMETER IN ARRAY AWAMAX, COMPUTE AMAX FROM 3.9 * d^(-1.2)
        ALLOCATE (AWAMAX(NODES))  
      ELSEIF (IAREA_FN.EQ.3) THEN  
C1D ------READ RHOW * GRAVITY / SIGMA COMPUTE AMAX FROM ROG_SIGMA * POROSITY
        ALLOCATE (ROG_SIGMA)
        ALLOCATE (AWAMAX(NODES))
        CALL URDCOM(IN,IOUT,LINE)
        LLOC=1
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,I,ROG_SIGMA,IOUT,INDIS)
        WRITE(IOUT,10) ROG_SIGMA
10      FORMAT(2X,'SP GRAVITY / SURFACE TENSION    (ROG_SIGMA) =',G15.6)
      ELSEIF (IAREA_FN.EQ.5) THEN 
C1E ------ALLOCATE TABULAR FUNCTION ARRAY FOR AREA VERSUS SATURATION        
        ALLOCATE(AWI_AREA_TAB(2,NATABROWS,NAZONES))
C              
      ENDIF    
C ---------------------------------------------------
C2 -------ALLOCATE ARRAYS FOR THE LANGMUIR ADSORPTION ISOTHERM 
      ALLOCATE (ALANGAW(NODES,MCOMP),BLANGAW(NODES,MCOMP))
C ----------------------------------------------------
C3A -------READ SIGMA / RT IF IKAWI_FN IS 3 (BRUSSEAU FORMULATION)        
      IF(IKAWI_FN.EQ.3) THEN 
        ALLOCATE (SIGMA_RT)      
        CALL URDCOM(IN,IOUT,LINE)
        LLOC=1
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,I,SIGMA_RT,IOUT,INDIS)
        WRITE(IOUT,11) SIGMA_RT
11      FORMAT(2X,'SURFACE TENSION / RT             (SIGMA_RT) =',G15.6)
      ELSEIF(IKAWI_FN.EQ.4) THEN 
C3B -------ALLOCATE TABULAR FUNCTIN ARRAYS FOR K_AWI VERSUS CONCENTRATION FOR EACH SPECIES    
        ALLOCATE(AWI_KAWI_TAB(2,NATABROWS,NAZONES,MCOMP))         
C          
      ENDIF
C ---------------------------------------------------------------------
C5 -------READ TABULAR FUNCTION ARRAY FOR AREA VERSUS SATURATION
      IF (IAREA_FN.EQ.5) THEN 
        DO IZON = 1,NAZONES
          WRITE(IOUT,65) IZON  
65        FORMAT(/10X,'AWI AREA-SAT TABLE NUMBER',I5/10X,27('-')/
     *  10X,'SATURATION',7X,'AREA')
          DO ITROWS = 1,NATABROWS  
            READ (IN,*) (AWI_AREA_TAB(I,ITROWS,IZON),I=1,2)  
            WRITE(IOUT,66)(AWI_AREA_TAB(I,ITROWS,IZON),I=1,2) 
66          FORMAT(10X,E14.6,3X,E14.6) 
          ENDDO  
        ENDDO            
      ENDIF
C
C6-------RETURN
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE AW_ADSORB1RP1 (IN)
C     ******************************************************************
C     READ AND PREPARE AREA ARRAYS FOR ADSORPTION ON AIR-WATER INTERFACE 
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,   ONLY:IOUT,NCOL,NROW,NLAY,IBOUND,BUFF,NODES,IUNSTR,
     1                   NODLAY 
      USE GWFBASMODULE,ONLY:MSUM,ICBCFL,IAUXSV,DELT,PERTIM,TOTIM,DELT
      USE GWTBCTMODULE, ONLY: MCOMP,PRSITY
      USE AW_ADSORBMODULE  
      REAL,    DIMENSION(:,:),    ALLOCATABLE ::TEMP      
      CHARACTER*24 ANAME(5)
      DATA ANAME(1) /'                 AW_AMAX'/
      DATA ANAME(2) /'               AWAREA_X2'/
      DATA ANAME(3) /'               AWAREA_X1'/
      DATA ANAME(4) /'               AWAREA_X0'/
      DATA ANAME(5) /'               GRAIN_DIA'/     
C     ------------------------------------------------------------------
C     
C1 ------FOR STRUCTURED GRID SETTINGS
      IF(IUNSTR.EQ.0) THEN 
C          
        ALLOCATE(TEMP(NCOL,NROW))          
        IF(IAREA_FN.EQ.1) THEN   
C1A--------READ AMAX ARRAY
          DO K=1,NLAY
            KK=K
            CALL U2DREL(TEMP(1,1),ANAME(1),NROW,NCOL,KK,IN,IOUT)
            DO I=1,NROW
            DO J=1,NCOL
              N=J+(I-1)*NCOL+(K-1)*NROW*NCOL
              AWAMAX(N) = TEMP(J,I)
            ENDDO
            ENDDO
          ENDDO      
        ELSEIF(IAREA_FN.EQ.4) THEN 
C1B------READ ARRAYS FOR X2, X1 AND X0
          DO K=1,NLAY
            KK=K
            CALL U2DREL(TEMP(1,1),ANAME(2),NROW,NCOL,KK,IN,IOUT)
            DO I=1,NROW
            DO J=1,NCOL
              N=J+(I-1)*NCOL+(K-1)*NROW*NCOL
              AWAREA_X2(N) = TEMP(J,I)
            ENDDO
            ENDDO
          ENDDO  
C          
          DO K=1,NLAY
            KK=K
            CALL U2DREL(TEMP(1,1),ANAME(3),NROW,NCOL,KK,IN,IOUT)
            DO I=1,NROW
            DO J=1,NCOL
              N=J+(I-1)*NCOL+(K-1)*NROW*NCOL
              AWAREA_X1(N) = TEMP(J,I)
            ENDDO
            ENDDO
          ENDDO           
C          
          DO K=1,NLAY
            KK=K
            CALL U2DREL(TEMP(1,1),ANAME(4),NROW,NCOL,KK,IN,IOUT)
            DO I=1,NROW
            DO J=1,NCOL
              N=J+(I-1)*NCOL+(K-1)*NROW*NCOL
              AWAREA_X0(N) = TEMP(J,I)
            ENDDO
            ENDDO
          ENDDO           
        ELSEIF(IAREA_FN.EQ.2) THEN 
C1C------READ ARRAY FOR GRAIN DIAMETER AND CONVERT TO AMAX 
          DO K=1,NLAY
            KK=K
            CALL U2DREL(TEMP(1,1),ANAME(5),NROW,NCOL,KK,IN,IOUT)
            DO I=1,NROW
            DO J=1,NCOL
              N=J+(I-1)*NCOL+(K-1)*NROW*NCOL
              AWAMAX(N) = TEMP(J,I)
              AWAMAX(N) = 3.9 / AWAMAX(N)**1.2
            ENDDO
            ENDDO
          ENDDO
        ELSEIF(IAREA_FN.EQ.3) THEN 
C1C------ AMAX IS ROG_SIGMA TIMES POROSITY
          DO K=1,NLAY
            KK=K
            DO I=1,NROW
            DO J=1,NCOL
              N=J+(I-1)*NCOL+(K-1)*NROW*NCOL
              AWAMAX(N) = ROG_SIGMA * PRSITY(N)
            ENDDO
            ENDDO
          ENDDO                  
        ENDIF  
        DEALLOCATE (TEMP)
C----------------------------------------------------------------------          
C2 ------ELE FOR UNSTRUCTURED GRID SETTINGS
      ELSEIF(IUNSTR.EQ.1) THEN 
C          
        IF(IAREA_FN.EQ.1) THEN   
C2A--------READ AMAX ARRAY
          DO K = 1,NLAY
            NNDLAY = NODLAY(K)
            NSTRT = NODLAY(K-1)+1
            NDSLAY = NNDLAY - NODLAY(K-1)
            CALL U1DREL(AWAMAX(NSTRT),ANAME(1),NDSLAY,K,IN,IOUT)
          ENDDO            
        ELSEIF(IAREA_FN.EQ.4) THEN 
C2B------READ ARRAYS FOR X2, X1 AND X0
          DO K = 1,NLAY
            NNDLAY = NODLAY(K)
            NSTRT = NODLAY(K-1)+1
            NDSLAY = NNDLAY - NODLAY(K-1)
            CALL U1DREL(AWAREA_X2(NSTRT),ANAME(2),NDSLAY,K,IN,IOUT)
          ENDDO
C
          DO K = 1,NLAY
            NNDLAY = NODLAY(K)
            NSTRT = NODLAY(K-1)+1
            NDSLAY = NNDLAY - NODLAY(K-1)
            CALL U1DREL(AWAREA_X1(NSTRT),ANAME(3),NDSLAY,K,IN,IOUT)
          ENDDO          
          DO K = 1,NLAY
            NNDLAY = NODLAY(K)
            NSTRT = NODLAY(K-1)+1
            NDSLAY = NNDLAY - NODLAY(K-1)
            CALL U1DREL(AWAREA_X0(NSTRT),ANAME(4),NDSLAY,K,IN,IOUT)
          ENDDO                   
        ELSEIF(IAREA_FN.EQ.2) THEN 
C2C------READ ARRAY FOR GRAIN DIAMETER AND CONVERT TO AMAX 
          DO K = 1,NLAY
            NNDLAY = NODLAY(K)
            NSTRT = NODLAY(K-1)+1
            NDSLAY = NNDLAY - NODLAY(K-1)
            CALL U1DREL(AWAMAX(NSTRT),ANAME(5),NDSLAY,K,IN,IOUT)
          ENDDO
          DO N = 1,NODES
            AWAMAX(N) = 3.9 / AWAMAX(N)**1.2   
          ENDDO 
      ELSEIF(IAREA_FN.EQ.3) THEN 
C2C------ AMAX IS ROG_SIGMA TIMES POROSITY 
          DO N = 1,NODES
            AWAMAX(N) = ROG_SIGMA * PRSITY(N)   
          ENDDO              
        ENDIF  
      ENDIF
C
C3-------RETURN
      RETURN
      END 
C-----------------------------------------------------------------------
      SUBROUTINE AW_ADSORB1RP2 (IN,ICOMP)
C     ******************************************************************
C     READ AND PREPARE lANGMUIR ARRAYS FOR ADSORPTION ON AIR-WATER INTERFACE 
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,   ONLY:IOUT,NCOL,NROW,NLAY,IBOUND,BUFF,NODES,IUNSTR,
     1                   NODLAY 
      USE GWFBASMODULE,ONLY:MSUM,ICBCFL,IAUXSV,DELT,PERTIM,TOTIM,DELT
      USE GWTBCTMODULE, ONLY: MCOMP
      USE AW_ADSORBMODULE  
      REAL,    DIMENSION(:,:),    ALLOCATABLE ::TEMP     
      REAL*8 AAA,BEE
      CHARACTER*24 ANAME(5)
      DATA ANAME(1) /'                ALANG_AW'/
      DATA ANAME(2) /'                BLANG_AW'/      
C     ------------------------------------------------------------------
C ----SKIP READING ARRAYS ALANG AND BLANG IF TABULAR INPUT FOR ISOTHERM     
      IF(IKAWI_FN.EQ.4) GO TO 10
C1 ------FOR STRUCTURED GRID SETTINGS
      IF(IUNSTR.EQ.0) THEN 
C          
        ALLOCATE(TEMP(NCOL,NROW)) 
C1A--------READ LANGMUIR A AND B ARRAYS
        DO K=1,NLAY
          KK=K
          CALL U2DREL(TEMP(1,1),ANAME(1),NROW,NCOL,KK,IN,IOUT)
          DO I=1,NROW
          DO J=1,NCOL
            N=J+(I-1)*NCOL+(K-1)*NROW*NCOL
            ALANGAW(N,ICOMP) = TEMP(J,I)
          ENDDO
          ENDDO
        ENDDO                  
C
        DO K=1,NLAY
          KK=K
          CALL U2DREL(TEMP(1,1),ANAME(2),NROW,NCOL,KK,IN,IOUT)
          DO I=1,NROW
          DO J=1,NCOL
            N=J+(I-1)*NCOL+(K-1)*NROW*NCOL
            BLANGAW(N,ICOMP) = TEMP(J,I)
          ENDDO
          ENDDO
        ENDDO    
        DEALLOCATE (TEMP) 
C------------------------------------------------------
C2 ------ELSE FOR UNSTRUCTURED GRID SETTINGS
      ELSEIF(IUNSTR.EQ.1) THEN 
C          
C2A--------READ LANGMUIR A AND B ARRAYS
        DO K = 1,NLAY
          NNDLAY = NODLAY(K)
          NSTRT = NODLAY(K-1)+1
          NDSLAY = NNDLAY - NODLAY(K-1)
          CALL U1DREL(ALANGAW(NSTRT,ICOMP),ANAME(1),NDSLAY,K,IN,IOUT)
        ENDDO 
C
        DO K = 1,NLAY
          NNDLAY = NODLAY(K)
          NSTRT = NODLAY(K-1)+1
          NDSLAY = NNDLAY - NODLAY(K-1)
          CALL U1DREL(BLANGAW(NSTRT,ICOMP),ANAME(2),NDSLAY,K,IN,IOUT)
        ENDDO 
      ENDIF  
C ---------------------------------------------------------------
10    CONTINUE
C ---------------------------------------------------------------
C3----------PREPARE LANGMUIR A AND B ARRAYS         
        IF(IKAWI_FN.EQ.1) THEN  
C3A -------A AND B ARE READ NOTHING TO PREPARE         
C          
        ELSEIF(IKAWI_FN.EQ.2) THEN   
C3B--------CMAX AND KL WERE READ. COMPUTE LANGMUIR A AND B 
C
          DO N=1,NODES
            ALANGAW(N,ICOMP) = ALANGAW(N,ICOMP) + BLANGAW(N,ICOMP)   
          ENDDO                           
        ELSEIF(IKAWI_FN.EQ.3) THEN   
C3C--------A_AW AND B_AW WERE READ. COMPUTE LANGMUIR A AND B FROM BRUSSEAU FORMULATION
C
          DO N=1,NODES
            BAW = BLANGAW(N,ICOMP)
            AAW = ALANGAW(N,ICOMP)
            BEE = 1.0 / AAW
            AAA = SIGMA_RT * BAW / AAW
            BLANGAW(N,ICOMP) = BEE 
            ALANGAW(N,ICOMP) = AAA
          ENDDO                           
        ELSEIF(IKAWI_FN.EQ.4) THEN                   
C3D--------TABULAR FUNCTION FOR K_AWI VERSUS CONCENTRATION IN WATER
          DO IZON = 1,NAZONES
          WRITE(IOUT,65) IZON, ICOMP  
65        FORMAT(/10X,'KAW-CONC TABLE NUMBER',I5/10X,27('-')/
     *  10X,'CONC',7X,'K_AW',3X,'FOR SPECIES', I5)
          DO ITROWS = 1,NATABROWS  
            READ (IN,*) (AWI_KAWI_TAB(I,ITROWS,IZON,ICOMP),I=1,2)  
            WRITE(IOUT,66)(AWI_KAWI_TAB(I,ITROWS,IZON,ICOMP),I=1,2) 
66          FORMAT(10X,E14.6,3X,E14.6) 
          ENDDO  
        ENDDO
        ENDIF 
C
C4-------RETURN
      RETURN
      END 
C-----------------------------------------------------------------------
      SUBROUTINE AWI_AREA 
C     ******************************************************************
C     COMPUTE AREA OF A-W INTERFACE FROM SATURATION OF WATER FOR VARIOUS OPTIONS
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,   ONLY:IOUT,NCOL,NROW,NLAY,IBOUND,BUFF,NODES,IUNSTR,
     1                   SN,SO  
      USE GWFBASMODULE,ONLY:MSUM,ICBCFL,IAUXSV,DELT,PERTIM,TOTIM,DELT
      USE GWTBCTMODULE, ONLY: MCOMP
      USE AW_ADSORBMODULE, ONLY: IAREA_FN, AREA_AWI,AREA_AWIO,AWAMAX,
     1  IAWIZONMAP,NATABROWS,NAZONES,AWI_AREA_TAB, AWAREA_X2,AWAREA_X1, 
     1  AWAREA_X0 
      DOUBLE PRECISION SW,SWO,TABAWI
C     ------------------------------------------------------------------
C1 -----FOR ALL NODES     
      DO N = 1,NODES
        SW = SN(N)
        SWO = SO(N)
        IF(IAREA_FN. EQ.1. OR. IAREA_FN. EQ.2. OR. IAREA_FN. EQ. 3) THEN
C2 -------WHEN AMAX IS USED          
          AREA_AWI(N) = AWAMAX(N) * (1.0 - SW)   
          AREA_AWIO(N) = AWAMAX(N) * (1.0 - SWO) 
        ELSEIF(IAREA_FN.EQ.4) THEN   
C3 -------WHEN X2, X1 AND X0 ARE USED 
          AREA_AWI(N)=AWAREA_X2(N)*SW*SW + AWAREA_X1(N)*SW+AWAREA_X0(N) 
         AREA_AWIO(N)=AWAREA_X2(N)*SWO*SWO+AWAREA_X1(N)*SWO+AWAREA_X0(N)
        ELSEIF(IAREA_FN.EQ.5) THEN   
C4- ------PICK FROM TABLE 
          IZON = IAWIZONMAP(N) 
          AREA_AWI(N) = tabAWI(SW,AWI_AREA_TAB,izon,natabrows,nazones) 
          AREA_AWIO(N) = tabAWI(SWO,AWI_AREA_TAB,izon,natabrows,nazones) 
        ENDIF   
      ENDDO     
C
C5-------RETURN
      RETURN
      END           
C-----------------------------------------------------------------------
      SUBROUTINE AWI_KAWI (ICOMP)
C     ******************************************************************
C     COMPUTE PARTITION COEFF OF A-W INTERFACE FROM CONC FOR LANGMUIR AND TABLE,
c     FILL STORAGE TERM ON A-W INTERFACE INTO MATRIX AND RHS       
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,   ONLY:IOUT,NCOL,NROW,NLAY,IBOUND,BUFF,NODES,IUNSTR,
     1                   AMAT,RHS,AREA,TOP,BOT,IA
      USE GWFBASMODULE,ONLY:MSUM,ICBCFL,IAUXSV,DELT,PERTIM,TOTIM,DELT
      USE GWTBCTMODULE, ONLY: MCOMP,ICBUND,CONC,CONCO
      USE AW_ADSORBMODULE, ONLY: IKAWI_FN, AK_AWI,ALANGAW,BLANGAW,
     1  AREA_AWI,AREA_AWIO,IAWIZONMAP,NATABROWS,NAZONES,AWI_KAWI_TAB  
      DOUBLE PRECISION CW,CWO,CWEPS,DC,ALANG,VOLU,VODT,ADSTERM,RT,DT,
     1  DTERM,RTERM,TABAWI,FL
C     ------------------------------------------------------------------        
C1-------INITIALIZE AND COMPUTE FOR ALL ACTIVE TRANSPORT NODES       
      DO N = 1,NODES 
        IF(ICBUND(N). EQ. 0) CYCLE 
        DTERM = 0.0
        RTERM = 0.0
C2-------GET VOL / DELT, AND NEW AND OLD CONCENTRATIONS        
        ALENG = TOP(N) - BOT(N)
        VOLU = AREA(N) * ALENG
        VODT = VOLU / DELT       
        CW = CONC(N,ICOMP)
        IF(CW.LT.0.0) CW = 0.0
        CWO = CONCO(N,ICOMP)
        IF(CWO.LT.0.0) CWO = 0.0
C3-------COMPUTE FOR DIFFERENT CASES AND FILL IN MATRIX        
        IF(IKAWI_FN. EQ.1. OR. IKAWI_FN. EQ.2. OR. IKAWI_FN. EQ. 3) THEN
C4 -------USE LANGMUIR ISOTHERM FILLED AS NEWTON     
          ADSTERM = ALANGAW(N,ICOMP) * VODT  
          FL = BLANGAW(N,ICOMP)
C          
          RT = AREA_AWI(N) * CW/(1.0 + FL*CW) - 
     1         AREA_AWIO(N) * CWO/(1.0 + FL*CWO)
          RT = RT *ADSTERM
          ADSTERM = ADSTERM * AREA_AWI(N) 
          DT = (1+FL*CW)*ADSTERM - ADSTERM*CW*FL ! ANALYTIC DERIV OF U/V = (VdU - UdV)/V^2
          DT = DT / (1.0 + FL*CW)**2         
        ELSEIF(IKAWI_FN.EQ.4) THEN  
C5 --------FILL TABULAR ISOTHERM TERMS ON LHS AND RHS
          IZON = IAWIZONMAP(N) 
          TAB_KN = tabAWI
     1     (CW,AWI_KAWI_TAB(1,1,1,ICOMP),izon,natabrows,nazones)
          TAB_KO = tabAWI
     1     (CWO,AWI_KAWI_TAB(1,1,1,ICOMP),izon,natabrows,nazones)
          ADSTERM = VODT 
          RT = TAB_KN * CW * AREA_AWI(N) -
     1         TAB_KO * CWO *  AREA_AWIO(N)
          RT = RT * ADSTERM
          EPS = 1.0E-4
          CWEPS = CW + EPS
          TAB_KEPS = tabAWI
     1     (CWEPS,AWI_KAWI_TAB(1,1,1,ICOMP),izon,natabrows,nazones)
          DT = (TAB_KEPS * CWEPS - TAB_KN * CW) / EPS
          DT = DT * ADSTERM * AREA_AWI(N)
        ENDIF    
C6 -------FILL DT AND RT TERMS ONTO LHS AND RHS
          DTERM = DTERM - DT
          RTERM = RTERM - DT * CW + RT    
C7 ---------FILL LHS AND RHS TERMS INTO MATRIX DIAGONAL AND RHS VECTOR        
        IPIV = IA(N)
        AMAT(IPIV) = AMAT(IPIV) +  DTERM
        RHS(N) = RHS(N) +  RTERM
       ENDDO         
C
C8-------RETURN
      RETURN
      END           
C-----------------------------------------------------------------------
      SUBROUTINE AWI_GETKAWI (ICOMP)
C     ******************************************************************
C     COMPUTE PARTITION COEFF OF A-W INTERFACE FROM CONC FOR LANGMUIR AND TABLE,
c     FILL STORAGE TERM ON A-W INTERFACE INTO MATRIX AND RHS       
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,   ONLY:IOUT,NCOL,NROW,NLAY,IBOUND,BUFF,NODES,IUNSTR,
     1                   AMAT,RHS,AREA,TOP,BOT,IA
      USE GWFBASMODULE,ONLY:MSUM,ICBCFL,IAUXSV,DELT,PERTIM,TOTIM,DELT
      USE GWTBCTMODULE, ONLY: MCOMP,ICBUND,CONC,CONCO
      USE AW_ADSORBMODULE, ONLY: IKAWI_FN, AK_AWI,ALANGAW,BLANGAW,
     1  IAWIZONMAP,NATABROWS,NAZONES,AWI_KAWI_TAB  
      DOUBLE PRECISION CW,CWO,CWEPS,DC,ALANG,VOLU,VODT,ADSTERM,RT,DT,
     1  DTERM,RTERM,TABAWI,FL
C     ------------------------------------------------------------------        
C1-------INITIALIZE AND COMPUTE FOR ALL ACTIVE TRANSPORT NODES       
      DO N = 1,NODES 
        IF(ICBUND(N). EQ. 0) CYCLE 
C2-------GET POSITIVE CW        
        CW = CONC(N,ICOMP)
        IF(CW.LT.0.0) CW = 0.0
C3-------COMPUTE FOR DIFFERENT CASES         
        IF(IKAWI_FN. EQ.1. OR. IKAWI_FN. EQ.2. OR. IKAWI_FN. EQ. 3) THEN
C4 -------USE LANGMUIR ISOTHERM FILLED AS NEWTON     
          ADSTERM = ALANGAW(N,ICOMP) 
          FL = BLANGAW(N,ICOMP)
C          
          AK_AWI(N,ICOMP) = ADSTERM/(1.0 + FL*CW)        
        ELSEIF(IKAWI_FN.EQ.4) THEN  
C5 --------FILL TABULAR ISOTHERM TERMS ON LHS AND RHS
          IZON = IAWIZONMAP(N) 
          TAB_KN = tabAWI
     1     (CW,AWI_KAWI_TAB(1,1,1,ICOMP),izon,natabrows,nazones)
          
          AK_AWI(N,ICOMP) = TAB_KN
        ENDIF    
       ENDDO         
C
C4-------RETURN
      RETURN
      END           
C-----------------------------------------------------------------------
      SUBROUTINE AW_ADSORB1BD(KSTP,KPER,ICOMP,ISS)
C     ******************************************************************
C     COMPUTE MASS BALANCE TERM FOR AIR-WATER INTERFACE ADSORPTION 
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,   ONLY:IOUT,NCOL,NROW,NLAY,IBOUND,BUFF,NODES,IUNSTR,
     1                   AREA,BOT,TOP,NODLAY 
      USE GWFBASMODULE,ONLY:MSUM,ICBCFL,IAUXSV,DELT,PERTIM,TOTIM,DELT
      USE GWTBCTMODULE, ONLY: MCOMP,VBNMT,VBVLT,CONC,CONCO,ICBUND,
     1  IBCTCB,ICT,MSUMT
      USE AW_ADSORBMODULE, ONLY: IKAWI_FN, ALANGAW, BLANGAW, AREA_AWI,
     1  IAWIZONMAP, AWI_KAWI_TAB, natabrows, nazones, AREA_AWIO 
C
      CHARACTER*16 TEXT
      DOUBLE PRECISION RATIN,RATOUT,QQ,VODT,ADSTERM,FL,CW,CWO,
     1  ALENG,VOLU,TAB_KO,TAB_KN,TABAWI
      DATA TEXT /'AD A-W INTERFACE'/
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
C3------LOOP THROUGH EACH NODE AND CALCULATE AIR-WATER-INTERFACE ADSORPTION
      DO 100 N=1,NODES
C
C4--------IF THE CELL IS INACTIVE, IGNORE IT.
        IF(ICBUND(N).EQ.0)GO TO 99
C
C5--------COMPUTE MASS
        ALENG = TOP(N) - BOT(N)
        VOLU = AREA(N) * ALENG
        VODT = VOLU / DELT
        QQ = 0.0
        Q = 0.0       
C-----------------------------------------------------------------------------
      IF(ICT.EQ.0)THEN  !----------WATER PHASE CONCENTRATION FORMULATION
C----------------------------------------------------------------------------
        CW = CONC(N,ICOMP)
        IF(CW.LT.0.0) CW = 0.0
        CWO = CONCO(N,ICOMP)
        IF(CWO.LT.0.0) CWO = 0.0
        IF(IKAWI_FN. EQ.1. OR. IKAWI_FN. EQ.2. OR. IKAWI_FN. EQ. 3) THEN
C6 -------USE LANGMUIR ISOTHERM      
          ADSTERM = ALANGAW(N,ICOMP) * VODT  
          FL = BLANGAW(N,ICOMP)
          QQ = ADSTERM * (AREA_AWI(N)*CW/(1.0+FL*CW) -
     1                    AREA_AWIO(N)*CWO/(1.0+FL*CWO))  
        ELSEIF(IKAWI_FN.EQ.4) THEN  
C7 --------USE TABULAR ISOTHERM 
          IZON = IAWIZONMAP(N) 
          TAB_KN = tabAWI
     1     (CW,AWI_KAWI_TAB(1,1,1,ICOMP),izon,natabrows,nazones)
          TAB_KO = tabAWI
     1     (CWO,AWI_KAWI_TAB(1,1,1,ICOMP),izon,natabrows,nazones)
          ADSTERM = VODT 
          QQ = TAB_KN * CW * AREA_AWI(N) - TAB_KO * CWO * AREA_AWI(N)
          QQ = QQ * ADSTERM
          ENDIF                
        ENDIF
C
        QQ = - QQ  ! STORAGE TERM NEGATIVE IS INFLOW AS PER MODFLOW CONVENTION
        Q = QQ
C
C8------PRINT FLOW RATE IF REQUESTED.
        IF(IBD.LT.0) THEN
         IF(IBDLBL.EQ.0) WRITE(IOUT,61) TEXT,KPER,KSTP
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
      double precision function tabAWI 
     *       (pvar,retcrvs,izon,natabrows,nazones)
C     ******************************************************************
C     compute value of tabAWI from table for given primaru variable pvar
C     ******************************************************************      
C           
C       specifications
C     ------------------------------------------------------------------
      real*8 pvar, hd, a, b, hir, fn ,Fir 
      dimension retcrvs(2,natabrows,nazones)
C     ------------------------------------------------------------------
c       ! set at first location if pvar is less than at first location
      hd = retcrvs(1,1,izon)
      if(pvar. le. hd) then 
        tabAWI =  retcrvs(2,1,izon) 
        return
      endif    
c      ! set at last location if pvar is greater than at last location 
      hd = retcrvs(1,natabrows,izon)
      if(pvar. ge. hd) then 
        tabAWI =  retcrvs(2,natabrows,izon) 
        return
      endif          
c      ! find location in table and interpolate for anything in between 
cc      ifound = 0
      do irow = 1,natabrows
         hir = retcrvs(1,irow,izon)
         if (pvar. gt. hir) go to 10 
cc         ifound = 1
         ir = irow
         irm1 = irow - 1
         a = retcrvs(2,irow,izon) - retcrvs(2,irm1,izon)
         b = hir - retcrvs(1,irm1,izon)
         fn = a / b * (pvar - hir)
         fir = retcrvs(2,irow,izon)
         tabAWI = Fir + fn
         go to 20
10       continue         
      enddo    
20    continue       
cc      ! pick end value when overshooting
cc      if(ifound.eq.0)then
cc        tabAWI = retcrvs(iv,natabrows,izon) 
cc      endif
c      
      return
      end function tabAWI        
C------------------------------------------------------------------------
      SUBROUTINE GWT2AWI_ADSORB1DA
      USE AW_ADSORBMODULE
      INTEGER ALLOC_ERR            
C
      DEALLOCATE(NAZONES, STAT = ALLOC_ERR)
      DEALLOCATE(NATABROWS, STAT = ALLOC_ERR)
      DEALLOCATE(IAWIZONMAP, STAT = ALLOC_ERR)
      DEALLOCATE(ROG_SIGMA, STAT = ALLOC_ERR)
      DEALLOCATE(SIGMA_RT, STAT = ALLOC_ERR)
      DEALLOCATE(AWAMAX, STAT = ALLOC_ERR)
      DEALLOCATE(AWAREA_X2, STAT = ALLOC_ERR)
      DEALLOCATE(AWAREA_X1, STAT = ALLOC_ERR)
      DEALLOCATE(AWAREA_X0, STAT = ALLOC_ERR)
      DEALLOCATE(AREA_AWI, STAT = ALLOC_ERR)
      DEALLOCATE(AREA_AWIO, STAT = ALLOC_ERR)
      DEALLOCATE(ALANGAW, STAT = ALLOC_ERR)
      DEALLOCATE(BLANGAW, STAT = ALLOC_ERR)
      DEALLOCATE(AK_AWI, STAT = ALLOC_ERR)
      DEALLOCATE(AWI_AREA_TAB, STAT = ALLOC_ERR)
      DEALLOCATE(AWI_KAWI_TAB, STAT = ALLOC_ERR) 
C
C16------RETURN
      RETURN
      END  
C----------------------------------------------------------------------
      SUBROUTINE GWT2DCYAW(N,ICOMP,DTERMD,RTERMD,VOLU)
C     ******************************************************************
C     FORMULATE DECAY TERM FOR EACH NODE AND COMPONENT ON A-W INTERFACE IN MOBILE DOMAIN
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,   ONLY:IOUT,NCOL,NROW,NLAY,IBOUND,TOP,NODLAY,ISYM,
     1 HNEW,AMAT,RHS,HOLD,ISSFLG,IA,JA,NJA,NODES,AREA,BOT,TOP,Sn,So
      USE GWTBCTMODULE, ONLY: ICT,CONC,ZODRAW,FODRAW,IZOD,IFOD
      USE GWFDPFMODULE, ONLY: PHIF,SnIM,SoIM
      USE CLN1MODULE, ONLY: ACLNNDS
      USE GWFBASMODULE,ONLY:DELT
      USE AW_ADSORBMODULE, ONLY: IKAWI_FN, AK_AWI,ALANGAW,BLANGAW,
     1  AREA_AWI,AREA_AWIO,IAWIZONMAP,NATABROWS,NAZONES,AWI_KAWI_TAB
      DOUBLE PRECISION ADSTERM,CW,FL,CWO,DT,RT,ALENG,VOLU,CEPS,
     1  X,Y,EPSS,EPS,CT,CWEPS,QA,QEPS,DQ,OMP,DTERMD,RTERMD,tabAWI
C     ------------------------------------------------------------------
        CW = CONC(N,ICOMP)
        IF(CW.LT.0.0) CW = 0.0

C---------ZERO ORDER DECAY ON A-W INTERFACE - APPLY NEWTON EXPANSION OF SUPPLY-DEMAND FUNCTION
          IF (IZOD.EQ.4.AND.IKAWI_FN.NE.0) THEN
C3-------FOR ZERO ORDER DECAY ON AIR-WATER INTERFACE USE LANGMUIR ISOTHERM FILLED AS NEWTON
            CT = - VOLU * ZODRAW(N,ICOMP) * AREA_AWI(N)
            EPS = 0.01
            CEPS = MAX(0.0,CW)
            X = CEPS /EPS
            CALL SMOOTH(X,Y)
            QA =  CT * Y
C-----------CALCULATE DQ/DH
            EPSS = 0.001 * EPS
            CEPS = MAX(0.0,CW+EPSS)
            X = (CEPS)/EPS
            CALL SMOOTH(X,Y)
            QEPS = CT * Y
            DQ = (QEPS - QA) / EPSS
            DTERMD = DTERMD + DQ
            RTERMD = RTERMD - QA + DQ*CW
          ENDIF
C -----------------------------------------------------------------------------
C---------FIRST ORDER DECAY ON A-W INTERFACE
          IF(IFOD.EQ.4.AND.IKAWI_FN.NE.0)THEN
            CT = -ALANGAW(N,ICOMP)*VOLU*FODRAW(N,ICOMP)*AREA_AWI(N)
            IF(IKAWI_FN. EQ.1.OR.IKAWI_FN.EQ.2.OR.IKAWI_FN.EQ. 3) THEN
C--------------FOR FUNCTIONAL CASES
C--------------FOR NON-LINEAR ADSORPTION FILL AS NEWTON
              ETA = BLANGAW(N,ICOMP)
              QA =  CT * CW / (1.0+ETA*CW)
C--------------CALCULATE DQ/DH
              EPSS = 0.00001
              CEPS = CW + EPSS
              QEPS =  CT * CEPS / (1.0 + ETA*CEPS)
              DQ = (QEPS - QA) / EPSS
              DTERMD = DTERMD + DQ
              RTERMD = RTERMD - QA + DQ*CW
            ELSEIF(IKAWI_FN.EQ.4) THEN
C5 --------FILL TABULAR ISOTHERM TERMS ON LHS AND RHS
              IZON = IAWIZONMAP(N)
              TAB_KN = tabAWI
     1          (CW,AWI_KAWI_TAB(1,1,1,ICOMP),izon,natabrows,nazones)
C
              EPS = 1.0E-4
              CEPS = CW + EPS
              TAB_KEPS = tabAWI
     1          (CEPS,AWI_KAWI_TAB(1,1,1,ICOMP),izon,natabrows,nazones)
              DT = (TAB_KEPS * CEPS - TAB_KN * CW) / EPS
              DT = DT * VOLU * AREA_AWI(N)
              DTERMD = DTERMD + DT
              RT = - TAB_KN * CW * AREA_AWI(N)
              RT = RT * VOLU
              RTERMD = RTERMD - RT + DT*CW
            ENDIF
          ENDIF
C
      RETURN
      END
C----------------------------------------------------------------------
      SUBROUTINE GWT2GENAW(N,ICOMP,DTERMD,RTERMD,VOLU,NPAREN)
C     ******************************************************************
C     FORMULATE DECAY TERM FOR EACH NODE AND COMPONENT ON A-W INTERFACE IN MOBILE DOMAIN
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,   ONLY:IOUT,NCOL,NROW,NLAY,IBOUND,TOP,NODLAY,ISYM,
     1 HNEW,AMAT,RHS,HOLD,ISSFLG,IA,JA,NJA,NODES,AREA,BOT,TOP,Sn,So
      USE GWTBCTMODULE, ONLY: ICT,STOTIO,JPARENT,NPARENT,SPTLRCT,CONC,
     1 ZODRAW,FODRAW,IZOD,IFOD,ISPRCT      
      USE GWFDPFMODULE, ONLY: PHIF,SnIM,SoIM
      USE CLN1MODULE, ONLY: ACLNNDS
      USE GWFBASMODULE,ONLY:DELT
      USE AW_ADSORBMODULE, ONLY: IKAWI_FN, AK_AWI,ALANGAW,BLANGAW,
     1  AREA_AWI,AREA_AWIO,IAWIZONMAP,NATABROWS,NAZONES,AWI_KAWI_TAB
      DOUBLE PRECISION ADSTERM,CW,FL,CWO,DT,RT,ALENG,VOLU,CEPS,
     1  X,Y,EPSS,EPS,CT,CWEPS,QA,QEPS,DQ,OMP,DTERMD,RTERMD,tabAWI,STOIT
C     ------------------------------------------------------------------
        IPAREN = JPARENT(ICOMP,NPAREN)
        IF(ISPRCT.EQ.0)THEN
          STOIT  = STOTIO(ICOMP,NPAREN)
        ELSE
          STOIT = SPTLRCT(N,ICOMP,NPAREN)
        ENDIF
        CW = CONC(N,IPAREN)
        IF(CW.LT.0.0) CW = 0.0
C----------------------------------------------------------------------
C---------ZERO ORDER DECAY ON A-W INTERFACE - APPLY NEWTON EXPANSION OF SUPPLY-DEMAND FUNCTION
          IF (IZOD.EQ.4.AND.IKAWI_FN.NE.0) THEN
C3-------FOR ZERO ORDER DECAY ON AIR-WATER INTERFACE USE LANGMUIR ISOTHERM FILLED AS NEWTON
            CT =  VOLU * ZODRAW(N,IPAREN) * AREA_AWI(N) * STOIT
            EPS = 0.01
            CEPS = MAX(0.0,CW)
            X = CEPS /EPS
            CALL SMOOTH(X,Y)
            QA =  CT * Y
C-----------CALCULATE DQ/DH
            EPSS = 0.001 * EPS
            CEPS = MAX(0.0,CW+EPSS)
            X = (CEPS)/EPS
            CALL SMOOTH(X,Y)
            QEPS = CT * Y
            DQ = (QEPS - QA) / EPSS
            DTERMD = DTERMD + DQ
            RTERMD = RTERMD - QA + DQ*CW
          ENDIF
C -----------------------------------------------------------------------------
C---------FIRST ORDER DECAY ON A-W INTERFACE
          IF(IFOD.EQ.4.AND.IKAWI_FN.NE.0)THEN
            CT = ALANGAW(N,ICOMP) * VOLU * FODRAW(N,IPAREN)
     *        * AREA_AWI(N) * STOIT
            IF(IKAWI_FN. EQ.1.OR.IKAWI_FN.EQ.2.OR.IKAWI_FN.EQ. 3) THEN
C--------------FOR FUNCTIONAL CASES
C--------------FOR NON-LINEAR ADSORPTION FILL AS NEWTON
              ETA = BLANGAW(N,ICOMP)
              QA =  CT * CW / (1.0+ETA*CW)
C--------------CALCULATE DQ/DH
              EPSS = 0.00001
              CEPS = CW + EPSS
              QEPS =  CT * CEPS / (1.0 + ETA*CEPS)
              DQ = (QEPS - QA) / EPSS
              DTERMD = DTERMD + DQ
              RTERMD = RTERMD - QA + DQ*CW
            ELSEIF(IKAWI_FN.EQ.4) THEN
C5 --------FILL TABULAR ISOTHERM TERMS ON LHS AND RHS
              IZON = IAWIZONMAP(N)
              TAB_KN = tabAWI
     1          (CW,AWI_KAWI_TAB(1,1,1,ICOMP),izon,natabrows,nazones)
C
              EPS = 1.0E-4
              CEPS = CW + EPS
              TAB_KEPS = tabAWI
     1          (CEPS,AWI_KAWI_TAB(1,1,1,ICOMP),izon,natabrows,nazones)
              DT = (TAB_KEPS * CEPS - TAB_KN * CW) / EPS
              DT = -DT * VOLU * AREA_AWI(N) * STOIT
              DTERMD = DTERMD + DT
              RT = - TAB_KN * CW * AREA_AWI(N)
              RT = - RT * VOLU  * STOIT
              RTERMD = RTERMD - RT + DT*CW
            ENDIF
          ENDIF
C
      RETURN
      END
C----------------------------------------------------------------------
      SUBROUTINE GWT2DCYAW1BD(N,ICOMP,QQ,VOLU)
C     ******************************************************************
C     FORMULATE DECAY TERM FOR EACH NODE AND COMPONENT ON A-W INTERFACE IN MOBILE DOMAIN
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,   ONLY:IOUT,NCOL,NROW,NLAY,IBOUND,TOP,NODLAY,ISYM,
     1 HNEW,AMAT,RHS,HOLD,ISSFLG,IA,JA,NJA,NODES,AREA,BOT,TOP,Sn,So
      USE GWTBCTMODULE, ONLY: ICT,CONC,ZODRAW,FODRAW,IZOD,IFOD
      USE GWFDPFMODULE, ONLY: PHIF,SnIM,SoIM
      USE CLN1MODULE, ONLY: ACLNNDS
      USE GWFBASMODULE,ONLY:DELT
      USE AW_ADSORBMODULE, ONLY: IKAWI_FN, AK_AWI,ALANGAW,BLANGAW,
     1  AREA_AWI,AREA_AWIO,IAWIZONMAP,NATABROWS,NAZONES,AWI_KAWI_TAB
      DOUBLE PRECISION ADSTERM,CW,FL,CWO,DT,RT,ALENG,VOLU,CEPS,
     1  X,Y,EPSS,EPS,CT,CWEPS,QQ,QEPS,DQ,OMP,tabAWI
C     ------------------------------------------------------------------
        CW = CONC(N,ICOMP)
        IF(CW.LT.0.0) CW = 0.0

C---------ZERO ORDER DECAY ON A-W INTERFACE - APPLY NEWTON EXPANSION OF SUPPLY-DEMAND FUNCTION
          IF (IZOD.EQ.4.AND.IKAWI_FN.NE.0) THEN
C3-------FOR ZERO ORDER DECAY ON AIR-WATER INTERFACE USE LANGMUIR ISOTHERM FILLED AS NEWTON
            CT = - VOLU * ZODRAW(N,ICOMP) * AREA_AWI(N)
            EPS = 0.01
            CEPS = MAX(0.0,CW)
            X = CEPS /EPS
            CALL SMOOTH(X,Y)
            QQ =  CT * Y
          ENDIF
C -----------------------------------------------------------------------------
C---------FIRST ORDER DECAY ON A-W INTERFACE
          IF(IFOD.EQ.4.AND.IKAWI_FN.NE.0)THEN
            CT = -ALANGAW(N,ICOMP)*VOLU*FODRAW(N,ICOMP)*AREA_AWI(N)
            IF(IKAWI_FN. EQ.1.OR.IKAWI_FN.EQ.2.OR.IKAWI_FN.EQ. 3) THEN
C--------------FOR FUNCTIONAL CASES
C--------------FOR NON-LINEAR ADSORPTION FILL AS NEWTON
              ETA = BLANGAW(N,ICOMP)
              QQ =  CT * CW / (1.0+ETA*CW)
            ELSEIF(IKAWI_FN.EQ.4) THEN
C5 --------FILL TABULAR ISOTHERM TERMS ON LHS AND RHS
              IZON = IAWIZONMAP(N)
              TAB_KN = tabAWI
     1          (CW,AWI_KAWI_TAB(1,1,1,ICOMP),izon,natabrows,nazones)
C
              RT = - TAB_KN * CW * AREA_AWI(N)
              RT = RT * VOLU
              QQ = QQ - RT
            ENDIF
          ENDIF
C
      RETURN
      END
C----------------------------------------------------------------------
      SUBROUTINE GWT2GENAW1BD(N,ICOMP,QQ,VOLU,NPAREN)
C     ******************************************************************
C     FORMULATE DECAY TERM FOR EACH NODE AND COMPONENT ON A-W INTERFACE IN MOBILE DOMAIN
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,   ONLY:IOUT,NCOL,NROW,NLAY,IBOUND,TOP,NODLAY,ISYM,
     1 HNEW,AMAT,RHS,HOLD,ISSFLG,IA,JA,NJA,NODES,AREA,BOT,TOP,Sn,So
      USE GWTBCTMODULE, ONLY: ICT,STOTIO,JPARENT,NPARENT,SPTLRCT,CONC,
     1  ZODRAW,FODRAW,IZOD,IFOD,ISPRCT 
      USE GWFDPFMODULE, ONLY: PHIF,SnIM,SoIM
      USE CLN1MODULE, ONLY: ACLNNDS
      USE GWFBASMODULE,ONLY:DELT
      USE AW_ADSORBMODULE, ONLY: IKAWI_FN, AK_AWI,ALANGAW,BLANGAW,
     1  AREA_AWI,AREA_AWIO,IAWIZONMAP,NATABROWS,NAZONES,AWI_KAWI_TAB
      DOUBLE PRECISION ADSTERM,CW,FL,CWO,DT,RT,ALENG,VOLU,CEPS,
     1  X,Y,EPSS,EPS,CT,CWEPS,QQ,OMP,tabAWI,STOIT
C     ------------------------------------------------------------------
        IPAREN = JPARENT(ICOMP,NPAREN)
        IF(ISPRCT.EQ.0)THEN
          STOIT  = STOTIO(ICOMP,NPAREN)
        ELSE
          STOIT = SPTLRCT(N,ICOMP,NPAREN)
        ENDIF
        CW = CONC(N,IPAREN)
        IF(CW.LT.0.0) CW = 0.0
C----------------------------------------------------------------------
C---------ZERO ORDER DECAY ON A-W INTERFACE - APPLY NEWTON EXPANSION OF SUPPLY-DEMAND FUNCTION
          IF (IZOD.EQ.4.AND.IKAWI_FN.NE.0) THEN
C3-------FOR ZERO ORDER DECAY ON AIR-WATER INTERFACE USE LANGMUIR ISOTHERM FILLED AS NEWTON
            CT =  VOLU * ZODRAW(N,IPAREN) * AREA_AWI(N) * STOIT
            EPS = 0.01
            CEPS = MAX(0.0,CW)
            X = CEPS /EPS
            CALL SMOOTH(X,Y)
            QQ =  CT * Y
          ENDIF
C -----------------------------------------------------------------------------
C---------FIRST ORDER DECAY ON A-W INTERFACE
          IF(IFOD.EQ.4.AND.IKAWI_FN.NE.0)THEN
            CT = ALANGAW(N,ICOMP) * VOLU * FODRAW(N,IPAREN)
     *        * AREA_AWI(N) * STOIT
            IF(IKAWI_FN. EQ.1.OR.IKAWI_FN.EQ.2.OR.IKAWI_FN.EQ. 3) THEN
C--------------FOR FUNCTIONAL CASES
C--------------FOR NON-LINEAR ADSORPTION FILL AS NEWTON
              ETA = BLANGAW(N,ICOMP)
              QQ =  CT * CW / (1.0+ETA*CW)
            ELSEIF(IKAWI_FN.EQ.4) THEN
C5 --------FILL TABULAR ISOTHERM TERMS ON LHS AND RHS
              IZON = IAWIZONMAP(N)
              TAB_KN = tabAWI
     1          (CW,AWI_KAWI_TAB(1,1,1,ICOMP),izon,natabrows,nazones)
C
              RT = - TAB_KN * CW * AREA_AWI(N)
              RT = -RT * VOLU  * STOIT
              QQ = QQ - RT
            ENDIF
          ENDIF
C
      RETURN
      END
      
 
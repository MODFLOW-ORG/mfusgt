C -------------------------------------------------------------------------------------      
C          
      SUBROUTINE AW_ADSORBIM1AL (IN) 
C     ******************************************************************
C     ALLOCATE ARRAYS FOR ADSORPTION ON AIR-WATER INTERFACE AND READ CONSTANTS 
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,   ONLY:IOUT,NCOL,NROW,NLAY,IBOUND,BUFF,NODES,IUNSTR 
      USE GWFBASMODULE,ONLY:MSUM,ICBCFL,IAUXSV,DELT,PERTIM,TOTIM,DELT
      USE GWTBCTMODULE, ONLY: MCOMP
      USE AW_ADSORBIMMODULE  
      CHARACTER*300 LINE      
      CHARACTER*24 ANAME 
      DATA ANAME /'TABULAR AWI FN. ZONE MAP'/      
C     ------------------------------------------------------------------
C      
      ALLOCATE(AREA_AWIIM(NODES),AREA_AWIIMO(NODES))
      ALLOCATE(AK_AWIIM(NODES,MCOMP)) 
C ---------------------------------------------------------------------
C1A ---------ALLOCATE SPACE AND READ ZONE MAP IF EITHER K_AWI OR AW_AREA ARE TABULAR
C   --------USE SAME ZONE MAP IF EITHER (OR BOTH) ARE TABULAR
      IF(IAREA_FNIM.EQ.5. OR. IKAWI_FNIM.EQ.4) THEN 
        ITAB_AWIIM = 1
        ALLOCATE(NAZONESIM, NATABROWSIM)
C1B -------READ NUMBER OF ZONES AND NUMBER OF ROWS IN ZONE TABLES        
        READ(IN,*) NAZONESIM, NATABROWSIM
C1C -------DIMENSION AND READ ZONE MAP          
        ALLOCATE(IAWIZONMAPIM(NODES))
        CALL U1DINT(IAWIZONMAPIM,ANAME,NODES,0,IN,IOUT)  
      ENDIF      
C1 -------FOR THE DIFFERENT AREA COMPUTATION OPTIONS            
      IF(IAREA_FNIM.EQ.1) THEN 
C1A ------READ AMAX, COMPUTE FROM A = AMAX * (1-Sw)
        ALLOCATE (AWAMAXIM(NODES))
      ELSEIF(IAREA_FNIM.EQ.4) THEN 
C1B ------READ X2, X1, X0 COMPUTE FROM A = X2*Sw^2 + X1*Sw + X0
        ALLOCATE (AWAREA_X2IM(NODES),AWAREA_X1IM(NODES),
     1    AWAREA_X0IM(NODES))
      ELSEIF(IAREA_FNIM.EQ.3) THEN 
C1C ------READ GRAIN DIAMETER IN ARRAY AWAMAXIM, COMPUTE AMAX FROM 3.9 * d^(-1.2)
        ALLOCATE (AWAMAXIM(NODES))  
      ELSEIF (IAREA_FN.EQ.2) THEN  
C1D ------READ RHOW * GRAVITY / SIGMA COMPUTE AMAX FROM ROG_SIGMA * POROSITY
        ALLOCATE (ROG_SIGMAIM)
        CALL URDCOM(IN,IOUT,LINE)
        LLOC=1
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,I,ROG_SIGMAIM,IOUT,INDIS)
        WRITE(IOUT,10) ROG_SIGMAIM
10      FORMAT(2X,'SP GRAVITY * GRAVITY / SURFACE TENSION ',1X,
     1    '  (ROG_SIGMA) =',G15.6)
      ELSEIF (IAREA_FNIM.EQ.5) THEN 
C1D ------ALLOCATE TABULAR FUNCTION ARRAY FOR AREA VERSUS SATURATION        
        ALLOCATE(AWI_AREA_TABIM(2,NATABROWSIM,NAZONESIM))
C              
      ENDIF    
C ---------------------------------------------------
C2 -------ALLOCATE ARRAYS FOR THE LANGMUIR ADSORPTION ISOTHERM 
      ALLOCATE (ALANGAWIM(NODES,MCOMP),BLANGAWIM(NODES,MCOMP))
C ----------------------------------------------------
      IF(IKAWI_FNIM.EQ.3) THEN 
        ALLOCATE (SIGMA_RTIM)      
        CALL URDCOM(IN,IOUT,LINE)
        LLOC=1
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,I,SIGMA_RTIM,IOUT,INDIS)
        WRITE(IOUT,11) SIGMA_RTIM
11      FORMAT(2X,'SURFACE TENSION / RT           (SIGMA_RTIM) =',G15.6)
      ELSEIF(IKAWI_FNIM.EQ.4) THEN 
C3B -------ALLOCATE TABULAR FUNCTIN ARRAYS FOR K_AWI VERSUS CONCENTRATION FOR EACH SPECIES    
        ALLOCATE(AWI_KAWI_TABIM(2,NATABROWSIM,NAZONESIM,MCOMP))         
C          
      ENDIF
C ---------------------------------------------------------------------
C5 -------READ TABULAR FUNCTION ARRAY FOR AREA VERSUS SATURATION
      IF (IAREA_FNIM.EQ.5) THEN 
        DO IZON = 1,NAZONESIM
          WRITE(IOUT,65) IZON  
65        FORMAT(/10X,'AWI AREA-SAT TABLE NUMBER',I5/10X,27('-')/
     *  10X,'SATURATION',7X,'AREA')
          DO ITROWS = 1,NATABROWSIM  
            READ (IN,*) (AWI_AREA_TABIM(I,ITROWS,IZON),I=1,2)  
            WRITE(IOUT,66)(AWI_AREA_TABIM(I,ITROWS,IZON),I=1,2) 
66          FORMAT(10X,E14.6,3X,E14.6) 
          ENDDO  
        ENDDO            
      ENDIF
C
C6-------RETURN
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE AW_ADSORBIM1RP1 (IN)
C     ******************************************************************
C     READ AND PREPARE AREA ARRAYS FOR ADSORPTION ON AIR-WATER INTERFACE 
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,   ONLY:IOUT,NCOL,NROW,NLAY,IBOUND,BUFF,NODES,IUNSTR,
     1                   NODLAY 
      USE GWFBASMODULE,ONLY:MSUM,ICBCFL,IAUXSV,DELT,PERTIM,TOTIM,DELT
      USE GWTBCTMODULE, ONLY: MCOMP
      USE GWTDPTMODULE, ONLY: PRSITYIM
      USE AW_ADSORBIMMODULE  
      REAL,    DIMENSION(:,:),    ALLOCATABLE ::TEMP      
      CHARACTER*24 ANAME(5)
      DATA ANAME(1) /'               AW_AMAXIM'/
      DATA ANAME(2) /'             AWAREA_X2IM'/
      DATA ANAME(3) /'             AWAREA_X1IM'/
      DATA ANAME(4) /'             AWAREA_X0IM'/
      DATA ANAME(5) /'             GRAIN_DIA'/     
C     ------------------------------------------------------------------
C     
C1 ------FOR STRUCTURED GRID SETTINGS
      IF(IUNSTR.EQ.0) THEN 
C          
        ALLOCATE(TEMP(NCOL,NROW))          
        IF(IAREA_FNIM.EQ.1) THEN   
C1A--------READ AMAX ARRAY
          DO K=1,NLAY
            KK=K
            CALL U2DREL(TEMP(1,1),ANAME(1),NROW,NCOL,KK,IN,IOUT)
            DO I=1,NROW
            DO J=1,NCOL
              N=J+(I-1)*NCOL+(K-1)*NROW*NCOL
              AWAMAXIM(N) = TEMP(J,I)
            ENDDO
            ENDDO
          ENDDO      
        ELSEIF(IAREA_FNIM.EQ.4) THEN 
C1B------READ ARRAYS FOR X2, X1 AND X0
          DO K=1,NLAY
            KK=K
            CALL U2DREL(TEMP(1,1),ANAME(2),NROW,NCOL,KK,IN,IOUT)
            DO I=1,NROW
            DO J=1,NCOL
              N=J+(I-1)*NCOL+(K-1)*NROW*NCOL
              AWAREA_X2IM(N) = TEMP(J,I)
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
              AWAREA_X1IM(N) = TEMP(J,I)
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
              AWAREA_X0IM(N) = TEMP(J,I)
            ENDDO
            ENDDO
          ENDDO           
        ELSEIF(IAREA_FNIM.EQ.2) THEN 
C1C------READ ARRAY FOR GRAIN DIAMETER AND CONVERT TO AMAX 
          DO K=1,NLAY
            KK=K
            CALL U2DREL(TEMP(1,1),ANAME(5),NROW,NCOL,KK,IN,IOUT)
            DO I=1,NROW
            DO J=1,NCOL
              N=J+(I-1)*NCOL+(K-1)*NROW*NCOL
              AWAMAXIM(N) = TEMP(J,I)
              AWAMAXIM(N) = 3.9 / AWAMAXIM(N)**1.2
            ENDDO
            ENDDO
          ENDDO        
        ELSEIF(IAREA_FNIM.EQ.3) THEN 
C1C------ AMAX IS ROG_SIGMA TIMES POROSITY
          DO K=1,NLAY
            KK=K
            DO I=1,NROW
            DO J=1,NCOL
              N=J+(I-1)*NCOL+(K-1)*NROW*NCOL
              AWAMAXIM(N) = ROG_SIGMAIM * PRSITYIM(N)
            ENDDO
            ENDDO
          ENDDO                         
        ENDIF  
        DEALLOCATE (TEMP)
C----------------------------------------------------------------------          
C2 ------ELE FOR UNSTRUCTURED GRID SETTINGS
      ELSEIF(IUNSTR.EQ.1) THEN 
C          
        IF(IAREA_FNIM.EQ.1) THEN   
C2A--------READ AMAX ARRAY
          DO K = 1,NLAY
            NNDLAY = NODLAY(K)
            NSTRT = NODLAY(K-1)+1
            NDSLAY = NNDLAY - NODLAY(K-1)
            CALL U1DREL(AWAMAXIM(NSTRT),ANAME(1),NDSLAY,K,IN,IOUT)
          ENDDO            
        ELSEIF(IAREA_FNIM.EQ.4) THEN 
C2B------READ ARRAYS FOR X2, X1 AND X0
          DO K = 1,NLAY
            NNDLAY = NODLAY(K)
            NSTRT = NODLAY(K-1)+1
            NDSLAY = NNDLAY - NODLAY(K-1)
            CALL U1DREL(AWAREA_X2IM(NSTRT),ANAME(2),NDSLAY,K,IN,IOUT)
          ENDDO
C
          DO K = 1,NLAY
            NNDLAY = NODLAY(K)
            NSTRT = NODLAY(K-1)+1
            NDSLAY = NNDLAY - NODLAY(K-1)
            CALL U1DREL(AWAREA_X1IM(NSTRT),ANAME(3),NDSLAY,K,IN,IOUT)
          ENDDO          
          DO K = 1,NLAY
            NNDLAY = NODLAY(K)
            NSTRT = NODLAY(K-1)+1
            NDSLAY = NNDLAY - NODLAY(K-1)
            CALL U1DREL(AWAREA_X0IM(NSTRT),ANAME(4),NDSLAY,K,IN,IOUT)
          ENDDO                   
        ELSEIF(IAREA_FNIM.EQ.2) THEN 
C2C------READ ARRAY FOR GRAIN DIAMETER AND CONVERT TO AMAX 
          DO K = 1,NLAY
            NNDLAY = NODLAY(K)
            NSTRT = NODLAY(K-1)+1
            NDSLAY = NNDLAY - NODLAY(K-1)
            CALL U1DREL(AWAMAXIM(NSTRT),ANAME(5),NDSLAY,K,IN,IOUT)
          ENDDO
          DO N = 1,NODES
            AWAMAXIM(N) = 3.9 / AWAMAXIM(N)**1.2   
          ENDDO   
        ELSEIF(IAREA_FNIM.EQ.3) THEN 
C2C------ AMAX IS ROG_SIGMA TIMES POROSITY 
          DO N = 1,NODES
            AWAMAXIM(N) = ROG_SIGMAIM * PRSITYIM(N)   
          ENDDO                        
        ENDIF  
      ENDIF
C
C3-------RETURN
      RETURN
      END 
C-----------------------------------------------------------------------
      SUBROUTINE AW_ADSORBIM1RP2 (IN,ICOMP)
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
      USE AW_ADSORBIMMODULE  
      REAL,    DIMENSION(:,:),    ALLOCATABLE ::TEMP     
      REAL*8 AAA,BEE
      CHARACTER*24 ANAME(5)
      DATA ANAME(1) /'              ALANG_AWIM'/
      DATA ANAME(2) /'              BLANG_AWIM'/      
C     ------------------------------------------------------------------
C ----SKIP READING ARRAYS ALANG AND BLANG IF TABULAR INPUT FOR ISOTHERM     
      IF(IKAWI_FNIM.EQ.4) GO TO 10
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
            ALANGAWIM(N,ICOMP) = TEMP(J,I)
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
            BLANGAWIM(N,ICOMP) = TEMP(J,I)
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
          CALL U1DREL(ALANGAWIM(NSTRT,ICOMP),ANAME(1),NDSLAY,K,IN,IOUT)
        ENDDO 
C
        DO K = 1,NLAY
          NNDLAY = NODLAY(K)
          NSTRT = NODLAY(K-1)+1
          NDSLAY = NNDLAY - NODLAY(K-1)
          CALL U1DREL(BLANGAWIM(NSTRT,ICOMP),ANAME(2),NDSLAY,K,IN,IOUT)
        ENDDO 
      ENDIF  
C ---------------------------------------------------------------
10    CONTINUE
C ---------------------------------------------------------------
C3----------PREPARE LANGMUIR A AND B ARRAYS         
        IF(IKAWI_FNIM.EQ.1) THEN  
C3A -------A AND B ARE READ NOTHING TO PREPARE         
C          
        ELSEIF(IKAWI_FNIM.EQ.2) THEN   
C3B--------CMAX AND KL WERE READ. COMPUTE LANGMUIR A AND B 
C
          DO N=1,NODES
            ALANGAWIM(N,ICOMP) = ALANGAWIM(N,ICOMP) + BLANGAWIM(N,ICOMP)
          ENDDO                           
        ELSEIF(IKAWI_FNIM.EQ.3) THEN   
C3C--------A_AW AND B_AW WERE READ. COMPUTE LANGMUIR A AND B FROM BRUSSEAU FORMULATION
C
          DO N=1,NODES
            BAW = BLANGAWIM(N,ICOMP)
            AAW = ALANGAWIM(N,ICOMP)
            BEE = 1.0 / AAW
            AAA = SIGMA_RTIM * BAW / AAW
            BLANGAWIM(N,ICOMP) = BEE 
            ALANGAWIM(N,ICOMP) = AAA
          ENDDO                           
        ELSEIF(IKAWI_FNIM.EQ.4) THEN                   
C3D--------TABULAR FUNCTION FOR K_AWI VERSUS CONCENTRATION IN WATER
          DO IZON = 1,NAZONESIM
          WRITE(IOUT,65) IZON, ICOMP  
65        FORMAT(/10X,'KAW-CONCIM TABLE NUMBER',I5/10X,27('-')/
     *  10X,'CONCIM',7X,'K_AW',3X,'FOR SPECIES', I5)
          DO ITROWS = 1,NATABROWSIM  
            READ (IN,*) (AWI_KAWI_TABIM(I,ITROWS,IZON,ICOMP),I=1,2)  
            WRITE(IOUT,66)(AWI_KAWI_TABIM(I,ITROWS,IZON,ICOMP),I=1,2) 
66          FORMAT(10X,E14.6,3X,E14.6) 
          ENDDO  
        ENDDO
        ENDIF 
C
C4-------RETURN
      RETURN
      END 
C-----------------------------------------------------------------------
      SUBROUTINE AWI_AREAIM 
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
      USE GWFDPFMODULE, ONLY: PHIF,SnIM,SoIM
      USE AW_ADSORBIMMODULE, ONLY: IAREA_FNIM, AREA_AWIIM,AWAMAXIM,
     1  IAWIZONMAPIM,NATABROWSIM,NAZONESIM,AWI_AREA_TABIM, 
     1  AWAREA_X2IM,AWAREA_X1IM, AWAREA_X0IM,AREA_AWIIMO 
      USE AW_ADSORBMODULE, ONLY: IAREA_FN, AREA_AWI,AWAMAX,
     1  IAWIZONMAP,NATABROWS,NAZONES,AWI_AREA_TAB, 
     1  AWAREA_X2,AWAREA_X1, AWAREA_X0,AREA_AWIO 
      DOUBLE PRECISION SW,SWO,TABAWI, SWIM,SWIMO
C     ------------------------------------------------------------------
C1 -----FOR ALL NODES     
      DO N = 1,NODES
        SW = SN(N)
        SWIM = SnIM(N)
        SWO = So(N)
        SWIOM = SoIM(N)        
        IF(IAREA_FNIM.EQ.1.OR.IAREA_FNIM.EQ.2.OR.IAREA_FNIM. EQ. 3) THEN
C2 -------WHEN AMAX IS USED          
          AREA_AWIIM(N) = AWAMAXIM(N) * (1.0 - SWIM)  * (1.0 - PHIF(N))
         AREA_AWIIMO(N) = AWAMAXIM(N) * (1.0 - SWIMO)  * (1.0 - PHIF(N))
        ELSEIF(IAREA_FNIM.EQ.4) THEN   
C3 -------WHEN X2, X1 AND X0 ARE USED 
          AREA_AWIIM(N) = (AWAREA_X2IM(N)*SWIM*SWIM +  
     1     AWAREA_X1IM(N) *SWIM + AWAREA_X0IM(N)) * (1.0 - PHIF(N))
          AREA_AWIIMO(N) = (AWAREA_X2IM(N)*SWIMO*SWIMO +  
     1     AWAREA_X1IM(N) *SWIMO + AWAREA_X0IM(N)) * (1.0 - PHIF(N))          
        ELSEIF(IAREA_FNIM.EQ.5) THEN   
C4- ------PICK FROM TABLE 
          IZON = IAWIZONMAPIM(N) 
          AREA_AWIIM(N) = tabAWI(SWIM,AWI_AREA_TABIM,izon,NATABROWSIM,
     1      NAZONESIM)* (1.0 - PHIF(N)) 
          AREA_AWIIMO(N) = tabAWI(SWIMO,AWI_AREA_TABIM,izon,NATABROWSIM,
     1      NAZONESIM)* (1.0 - PHIF(N))           
        ENDIF   
C ----------------------------------------------------------------------
C5 -------ADJUST AREAS OF MOBILE DOMAIN BY PHIF
        IF(IAREA_FN.EQ.1.OR.IAREA_FN.EQ.2.OR.IAREA_FN. EQ. 3) THEN
C6 -------WHEN AMAX IS USED          
          AREA_AWI(N) = AWAMAX(N) * SW  * PHIF(N)
          AREA_AWIO(N) = AWAMAX(N) * SWO  * PHIF(N)
        ELSEIF(IAREA_FN.EQ.4) THEN   
C7 -------WHEN X2, X1 AND X0 ARE USED 
          AREA_AWI(N) = (AWAREA_X2(N)*SW*SW + AWAREA_X1(N) *SW + 
     1      AWAREA_X0(N)) * PHIF(N)
          AREA_AWIO(N) = (AWAREA_X2(N)*SWO*SWO + AWAREA_X1(N) *SWO + 
     1      AWAREA_X0(N)) * PHIF(N)          
        ELSEIF(IAREA_FN.EQ.5) THEN   
C8 -------PICK FROM TABLE 
          IZON = IAWIZONMAP(N) 
          AREA_AWI(N) = tabAWI(SW,AWI_AREA_TAB,izon,NATABROWS,
     1      NAZONES)* PHIF(N) 
          AREA_AWIO(N) = tabAWI(SWO,AWI_AREA_TAB,izon,NATABROWS,
     1      NAZONES)* PHIF(N)           
        ENDIF   
      ENDDO     
C
C9-------RETURN
      RETURN
      END           
C-----------------------------------------------------------------------
      SUBROUTINE AWI_KAWIIM (ICOMP)
C     ******************************************************************
C     COMPUTE PARTITION COEFF OF A-W INTERFACE FROM CONCIM FOR LANGMUIR AND TABLE,
c     FILL STORAGE TERM ON A-W INTERFACE INTO MATRIX AND RHS       
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,   ONLY:IOUT,NCOL,NROW,NLAY,IBOUND,BUFF,NODES,IUNSTR,
     1                   AMAT,RHS,AREA,TOP,BOT,IA
      USE GWFBASMODULE,ONLY:MSUM,ICBCFL,IAUXSV,DELT,PERTIM,TOTIM,DELT
      USE GWTDPTMODULE, ONLY: ICBUNDIM,CONCIM,CONCOIM,DIADDT,RDDT     
      USE AW_ADSORBIMMODULE, ONLY: IKAWI_FNIM, AK_AWIIM,ALANGAWIM,
     1  BLANGAWIM,AREA_AWIIM,IAWIZONMAPIM,NATABROWSIM,NAZONESIM,
     1  AWI_KAWI_TABIM,AREA_AWIIMO  
      DOUBLE PRECISION CW,CWO,CWEPS,DC,ALANG,VOLU,VODT,ADSTERM,RT,DT,
     1  DTERM,RTERM,TABAWI
C     ------------------------------------------------------------------        
C1-------INITIALIZE AND COMPUTE FOR ALL ACTIVE TRANSPORT NODES       
      DO N = 1,NODES 
        IF(ICBUNDIM(N). LE. 0) CYCLE 
        DTERM = 0.0
        RTERM = 0.0
C2-------GET VOL / DELT, AND NEW AND OLD CONCIMENTRATIONS        
        ALENG = TOP(N) - BOT(N)
        VOLU = AREA(N) * ALENG
        VODT = VOLU / DELT       
        CW = CONCIM(N,ICOMP)
        IF(CW.LT.0.0) CW = 0.0
        CWO = CONCOIM(N,ICOMP)
        IF(CWO.LT.0.0) CWO = 0.0
C3-------COMPUTE FOR DIFFERENT CASES AND FILL IN MATRIX        
        IF(IKAWI_FNIM.EQ.1.OR.IKAWI_FNIM.EQ.2.OR.IKAWI_FNIM. EQ. 3) THEN
C4 -------USE LANGMUIR ISOTHERM FILLED AS NEWTON     
          ADSTERM = ALANGAWIM(N,ICOMP) * VODT  
          FL = BLANGAWIM(N,ICOMP)
C          
          RT = AREA_AWIIM(N)*CW/(1.0 + FL*CW) - 
     1         AREA_AWIIMO(N)*CWO/(1.0 + FL*CWO)
          RT = RT *ADSTERM
          ADSTERM = ADSTERM * AREA_AWIIM(N)
          DT = (1+FL*CW)*ADSTERM - ADSTERM*CW*FL ! ANALYTIC DERIV OF U/V = (VdU - UdV)/V^2
          DT = DT / (1.0 + FL*CW)**2         
        ELSEIF(IKAWI_FNIM.EQ.4) THEN  
C5 --------FILL TABULAR ISOTHERM TERMS ON LHS AND RHS
          IZON = IAWIZONMAPIM(N) 
          TAB_KN = tabAWI
     1     (CW,AWI_KAWI_TABIM(1,1,1,ICOMP),izon,NATABROWSIM,NAZONESIM)
          TAB_KO = tabAWI
     1     (CWO,AWI_KAWI_TABIM(1,1,1,ICOMP),izon,NATABROWSIM,NAZONESIM)
          ADSTERM = VODT 
          RT = AREA_AWIIM(N)*TAB_KN * CW - AREA_AWIIMO(N)*TAB_KO * CWO
          RT = RT * ADSTERM
          EPS = 1.0E-4
          CWEPS = CW + EPS
          TAB_KEPS = tabAWI
     1    (CWEPS,AWI_KAWI_TABIM(1,1,1,ICOMP),izon,NATABROWSIM,NAZONESIM)
          DT = (TAB_KEPS * CWEPS - TAB_KN * CW) / EPS
          DT = DT * ADSTERM* AREA_AWIIM(N)
        ENDIF    
C6 -------FILL DT AND RT TERMS ONTO LHS AND RHS
        DTERM = DTERM - DT
        RTERM = RTERM - DT * CW + RT    
C7 ---------FILL LHS AND RHS TERMS INTO MATRIX DIAGONAL AND RHS VECTOR        
        DIADDT(N)=DIADDT(N) + DTERM
        RDDT(N)=RDDT(N) + RTERM 
       ENDDO         
C
C8-------RETURN
      RETURN
      END           
C-----------------------------------------------------------------------
      SUBROUTINE AWI_GETKAWIIM (ICOMP)
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
      USE GWTDPTMODULE, ONLY: ICBUNDIM,CONCIM,CONCOIM,DIADDT,RDDT     
      USE AW_ADSORBIMMODULE, ONLY: IKAWI_FNIM, AK_AWIIM,ALANGAWIM,
     1  BLANGAWIM,IAWIZONMAPIM,NATABROWSIM,NAZONESIM,
     1  AWI_KAWI_TABIM  
      DOUBLE PRECISION CW,CWO,CWEPS,DC,ALANG,VOLU,VODT,ADSTERM,RT,DT,
     1  DTERM,RTERM,TABAWI
C     ------------------------------------------------------------------        
C1-------INITIALIZE AND COMPUTE FOR ALL ACTIVE TRANSPORT NODES       
      DO N = 1,NODES 
        IF(ICBUNDIM(N). EQ. 0) CYCLE 
C2-------GET VOL        
        ALENG = TOP(N) - BOT(N)
        VOLU = AREA(N) * ALENG
        CW = CONCIM(N,ICOMP)
        IF(CW.LT.0.0) CW = 0.0
C3-------COMPUTE FOR DIFFERENT CASES         
        IF(IKAWI_FNIM.EQ.1.OR.IKAWI_FNIM.EQ.2.OR.IKAWI_FNIM.EQ.3) THEN
C4 -------USE LANGMUIR ISOTHERM FILLED AS NEWTON     
          ADSTERM = ALANGAWIM(N,ICOMP) * VOLU
          FL = BLANGAWIM(N,ICOMP)
C          
          AK_AWIIM(N,ICOMP) = ADSTERM/(1.0 + FL*CW)        
        ELSEIF(IKAWI_FNIM.EQ.4) THEN  
C5 --------FILL TABULAR ISOTHERM TERMS ON LHS AND RHS
          IZON = IAWIZONMAPIM(N) 
          TAB_KN = tabAWI
     1     (CW,AWI_KAWI_TABIM(1,1,1,ICOMP),izon,natabrowsIM,nazonesIM)
          
          AK_AWIIM(N,ICOMP) = TAB_KN * VOLU
        ENDIF    
       ENDDO         
C
C4-------RETURN
      RETURN
      END           
C-----------------------------------------------------------------------
      SUBROUTINE AW_ADSORBIM1BD(KSTP,KPER,ICOMP,ISS)
C     ******************************************************************
C     COMPUTE MASS BALANCE TERM FOR AIR-WATER INTERFACE ADSORPTION 
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,   ONLY:IOUT,NCOL,NROW,NLAY,IBOUND,BUFF,NODES,IUNSTR,
     1                   AREA,BOT,TOP,NODLAY 
      USE GWFBASMODULE,ONLY:MSUM,ICBCFL,IAUXSV,DELT,PERTIM,TOTIM,DELT
      USE GWTBCTMODULE, ONLY: MCOMP,VBNMT,VBVLT,IBCTCB,ICT,MSUMT
      USE GWTDPTMODULE, ONLY: ICBUNDIM,CONCIM,CONCOIM
      USE AW_ADSORBIMMODULE, ONLY: IKAWI_FNIM, ALANGAWIM, BLANGAWIM, 
     1  AREA_AWIIM,IAWIZONMAPIM, AWI_KAWI_TABIM, NATABROWSIM, NAZONESIM,
     1  AREA_AWIIMO 
C
      CHARACTER*16 TEXT
      DOUBLE PRECISION RATIN,RATOUT,QQ,VODT,ADSTERM,FL,CW,CWO,
     1  ALENG,VOLU,TAB_KO,TAB_KN,TABAWI
      DATA TEXT /'AD A-W INTRFC.IM'/
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
        IF(ICBUNDIM(N).EQ.0)GO TO 99
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
        CW = CONCIM(N,ICOMP)
        IF(CW.LT.0.0) CW = 0.0
        CWO = CONCOIM(N,ICOMP)
        IF(CWO.LT.0.0) CWO = 0.0
        IF(IKAWI_FNIM.EQ.1.OR.IKAWI_FNIM.EQ.2.OR.IKAWI_FNIM. EQ. 3) THEN
C6 -------USE LANGMUIR ISOTHERM      
          ADSTERM = ALANGAWIM(N,ICOMP) * VODT  
          FL = BLANGAWIM(N,ICOMP)
          QQ = ADSTERM * (AREA_AWIIM(N)*CW/(1.0+FL*CW) - 
     1                    AREA_AWIIMO(N)*CWO/(1.0+FL*CWO))  
        ELSEIF(IKAWI_FNIM.EQ.4) THEN  
C7 --------USE TABULAR ISOTHERM 
          IZON = IAWIZONMAPIM(N) 
          TAB_KN = tabAWI
     1     (CW,AWI_KAWI_TABIM(1,1,1,ICOMP),izon,NATABROWSIM,NAZONESIM)
          TAB_KO = tabAWI
     1     (CWO,AWI_KAWI_TABIM(1,1,1,ICOMP),izon,NATABROWSIM,NAZONESIM)
          ADSTERM = VODT 
          QQ = AREA_AWIIM(N)*TAB_KN * CW - AREA_AWIIMO(N)*TAB_KO * CWO
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
      SUBROUTINE DPT2AWI_ADSORB1DA
      USE AW_ADSORBIMMODULE
      INTEGER ALLOC_ERR            
C
      DEALLOCATE(NAZONESIM, STAT = ALLOC_ERR)
      DEALLOCATE(NATABROWSIM, STAT = ALLOC_ERR)
      DEALLOCATE(IAWIZONMAPIM, STAT = ALLOC_ERR)
      DEALLOCATE(ROG_SIGMAIM, STAT = ALLOC_ERR)
      DEALLOCATE(SIGMA_RTIM, STAT = ALLOC_ERR)
      DEALLOCATE(AWAMAXIM, STAT = ALLOC_ERR)
      DEALLOCATE(AWAREA_X2IM, STAT = ALLOC_ERR)
      DEALLOCATE(AWAREA_X1IM, STAT = ALLOC_ERR)
      DEALLOCATE(AWAREA_X0IM, STAT = ALLOC_ERR)
      DEALLOCATE(AREA_AWIIM, STAT = ALLOC_ERR)
      DEALLOCATE(AREA_AWIIMO, STAT = ALLOC_ERR)
      DEALLOCATE(ALANGAWIM, STAT = ALLOC_ERR)
      DEALLOCATE(BLANGAWIM, STAT = ALLOC_ERR)
      DEALLOCATE(AK_AWIIM, STAT = ALLOC_ERR)
      DEALLOCATE(AWI_AREA_TABIM, STAT = ALLOC_ERR)
      DEALLOCATE(AWI_KAWI_TABIM, STAT = ALLOC_ERR) 
C
C16------RETURN
      RETURN
      END  
C----------------------------------------------------------------------
      SUBROUTINE GWT2DPT1DCYAWIM(N,ICOMP,DTERMD,RTERMD,VOLU)
C     ******************************************************************
C     FORMULATE DECAY TERM FOR EACH NODE AND COMPONENT ON A-W INTERFACE IN IMMOBILE DOMAIN
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,   ONLY:IOUT,NCOL,NROW,NLAY,IBOUND,TOP,NODLAY,ISYM,
     1 HNEW,AMAT,RHS,HOLD,ISSFLG,IA,JA,NJA,NODES,AREA,BOT,TOP,Sn,So
      USE GWTBCTMODULE, ONLY: ICT
      USE GWFDPFMODULE, ONLY: PHIF,SnIM,SoIM
      USE CLN1MODULE, ONLY: ACLNNDS
      USE GWFBASMODULE,ONLY:DELT
      USE GWTDPTMODULE, ONLY: CONCIM,CONCOIM,IADSORBIM,PRSITYIM,
     * ADSORBIM,FLICHIM,IZODIM,IFODIM,ZODRAWIM,FODRAWIM,
     1  ZODRSIM,FODRSIM
      USE AW_ADSORBMODULE, ONLY: IKAWI_FN, AK_AWI,ALANGAW,BLANGAW,
     1  AREA_AWI,AREA_AWIO,IAWIZONMAP,NATABROWS,NAZONES,AWI_KAWI_TAB  
      USE AW_ADSORBIMMODULE, ONLY: ALANGAWIM,BLANGAWIM,AK_AWIIM
      DOUBLE PRECISION ADSTERM,CW,FL,CWO,DT,RT,ALENG,VOLU,CEPS,
     1  X,Y,EPSS,EPS,CT,CWEPS,QA,QEPS,DQ,OMP,DTERMD,RTERMD,tabAWI
C     ------------------------------------------------------------------      
        CW = CONCIM(N,ICOMP)
        IF(CW.LT.0.0) CW = 0.0
      
C---------ZERO ORDER DECAY ON A-W INTERFACE - APPLY NEWTON EXPANSION OF SUPPLY-DEMAND FUNCTION
          IF (IZODIM.GE.2.AND.IKAWI_FN.NE.0) THEN 
C3-------FOR ZERO ORDER DECAY ON AIR-WATER INTERFACE USE LANGMUIR ISOTHERM FILLED AS NEWTON
            CT = - VOLU * ZODRAWIM(N,ICOMP) * AREA_AWI(N)
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
          IF(IFODIM.GE.2.AND.IKAWI_FN.NE.0)THEN
            CT = -ALANGAWIM(N,ICOMP)*VOLU*FODRAWIM(N,ICOMP)*AREA_AWI(N)
            IF(IKAWI_FN. EQ.1.OR.IKAWI_FN.EQ.2.OR.IKAWI_FN.EQ. 3) THEN
C--------------FOR FUNCTIONAL CASES
C--------------FOR NON-LINEAR ADSORPTION FILL AS NEWTON
              ETA = BLANGAWIM(N,ICOMP)
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
      SUBROUTINE GWT2DPT1GENAWIM(N,ICOMP,DTERMD,RTERMD,VOLU,NPAREN)
C     ******************************************************************
C     FORMULATE DECAY TERM FOR EACH NODE AND COMPONENT ON A-W INTERFACE IN IMMOBILE DOMAIN
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,   ONLY:IOUT,NCOL,NROW,NLAY,IBOUND,TOP,NODLAY,ISYM,
     1 HNEW,AMAT,RHS,HOLD,ISSFLG,IA,JA,NJA,NODES,AREA,BOT,TOP,Sn,So
      USE GWTBCTMODULE, ONLY: ICT,STOTIO,JPARENT,NPARENT,ISPRCT
      USE GWFDPFMODULE, ONLY: PHIF,SnIM,SoIM
      USE CLN1MODULE, ONLY: ACLNNDS
      USE GWFBASMODULE,ONLY:DELT
      USE GWTDPTMODULE, ONLY: CONCIM,CONCOIM,IADSORBIM,PRSITYIM,
     * ADSORBIM,FLICHIM,IZODIM,IFODIM,ZODRAWIM,FODRAWIM,
     1  ZODRSIM,FODRSIM, SPTLRCTIM 
      USE AW_ADSORBMODULE, ONLY: IKAWI_FN, AK_AWI,ALANGAW,BLANGAW,
     1  AREA_AWI,AREA_AWIO,IAWIZONMAP,NATABROWS,NAZONES,AWI_KAWI_TAB  
      USE AW_ADSORBIMMODULE, ONLY: ALANGAWIM,BLANGAWIM,AK_AWIIM
      DOUBLE PRECISION ADSTERM,CW,FL,CWO,DT,RT,ALENG,VOLU,CEPS,
     1  X,Y,EPSS,EPS,CT,CWEPS,QA,QEPS,DQ,OMP,DTERMD,RTERMD,tabAWI,STOIT
C     ------------------------------------------------------------------      
        IPAREN = JPARENT(ICOMP,NPAREN)
        IF(ISPRCT.EQ.0)THEN
          STOIT  = STOTIO(ICOMP,NPAREN)
        ELSE
          STOIT = SPTLRCTIM(N,ICOMP,NPAREN)
        ENDIF
        CW = CONCIM(N,IPAREN)
        IF(CW.LT.0.0) CW = 0.0        
C----------------------------------------------------------------------        
C---------ZERO ORDER DECAY ON A-W INTERFACE - APPLY NEWTON EXPANSION OF SUPPLY-DEMAND FUNCTION 
          IF (IZODIM.EQ.4.AND.IKAWI_FN.NE.0) THEN 
C3-------FOR ZERO ORDER DECAY ON AIR-WATER INTERFACE USE LANGMUIR ISOTHERM FILLED AS NEWTON
            CT =  VOLU * ZODRAWIM(N,IPAREN) * AREA_AWI(N) * STOIT
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
          IF(IFODIM.EQ.4.AND.IKAWI_FN.NE.0)THEN
            CT = ALANGAWIM(N,ICOMP) * VOLU * FODRAWIM(N,IPAREN)
     *        * AREA_AWI(N) * STOIT
            IF(IKAWI_FN. EQ.1.OR.IKAWI_FN.EQ.2.OR.IKAWI_FN.EQ. 3) THEN
C--------------FOR FUNCTIONAL CASES
C--------------FOR NON-LINEAR ADSORPTION FILL AS NEWTON
              ETA = BLANGAWIM(N,ICOMP)
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
              DT = DT * VOLU * AREA_AWI(N) * STOIT
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
      SUBROUTINE GWT2DPT1DCYAWIM1BD(N,ICOMP,QQ,VOLU)
C     ******************************************************************
C     FORMULATE DECAY TERM FOR EACH NODE AND COMPONENT ON A-W INTERFACE IN IMMOBILE DOMAIN
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,   ONLY:IOUT,NCOL,NROW,NLAY,IBOUND,TOP,NODLAY,ISYM,
     1 HNEW,AMAT,RHS,HOLD,ISSFLG,IA,JA,NJA,NODES,AREA,BOT,TOP,Sn,So
      USE GWTBCTMODULE, ONLY: ICT
      USE GWFDPFMODULE, ONLY: PHIF,SnIM,SoIM
      USE CLN1MODULE, ONLY: ACLNNDS
      USE GWFBASMODULE,ONLY:DELT
      USE GWTDPTMODULE, ONLY: CONCIM,CONCOIM,IADSORBIM,PRSITYIM,
     * ADSORBIM,FLICHIM,IZODIM,IFODIM,ZODRAWIM,FODRAWIM,
     1  ZODRSIM,FODRSIM
      USE AW_ADSORBMODULE, ONLY: IKAWI_FN, AK_AWI,ALANGAW,BLANGAW,
     1  AREA_AWI,AREA_AWIO,IAWIZONMAP,NATABROWS,NAZONES,AWI_KAWI_TAB  
      USE AW_ADSORBIMMODULE, ONLY: ALANGAWIM,BLANGAWIM,AK_AWIIM
      DOUBLE PRECISION ADSTERM,CW,FL,CWO,DT,RT,ALENG,VOLU,CEPS,
     1  X,Y,EPSS,EPS,CT,CWEPS,QQ,QEPS,DQ,OMP,tabAWI
C     ------------------------------------------------------------------      
        CW = CONCIM(N,ICOMP)
        IF(CW.LT.0.0) CW = 0.0
      
C---------ZERO ORDER DECAY ON A-W INTERFACE - APPLY NEWTON EXPANSION OF SUPPLY-DEMAND FUNCTION
          IF (IZODIM.EQ.4.AND.IKAWI_FN.NE.0) THEN 
C3-------FOR ZERO ORDER DECAY ON AIR-WATER INTERFACE USE LANGMUIR ISOTHERM FILLED AS NEWTON
            CT = - VOLU * ZODRAWIM(N,ICOMP) * AREA_AWI(N)
            EPS = 0.01
            CEPS = MAX(0.0,CW)
            X = CEPS /EPS
            CALL SMOOTH(X,Y)
            QQ =  CT * Y
          ENDIF
C -----------------------------------------------------------------------------
C---------FIRST ORDER DECAY ON A-W INTERFACE
          IF(IFODIM.EQ.4.AND.IKAWI_FN.NE.0)THEN
            CT = -ALANGAWIM(N,ICOMP)*VOLU*FODRAWIM(N,ICOMP)*AREA_AWI(N)
            IF(IKAWI_FN. EQ.1.OR.IKAWI_FN.EQ.2.OR.IKAWI_FN.EQ. 3) THEN
C--------------FOR FUNCTIONAL CASES
C--------------FOR NON-LINEAR ADSORPTION FILL AS NEWTON
              ETA = BLANGAWIM(N,ICOMP)
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
      SUBROUTINE GWT2DPT1GENAWIM1BD(N,ICOMP,QQ,VOLU,NPAREN)
C     ******************************************************************
C     FORMULATE DECAY TERM FOR EACH NODE AND COMPONENT ON A-W INTERFACE IN IMMOBILE DOMAIN
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,   ONLY:IOUT,NCOL,NROW,NLAY,IBOUND,TOP,NODLAY,ISYM,
     1 HNEW,AMAT,RHS,HOLD,ISSFLG,IA,JA,NJA,NODES,AREA,BOT,TOP,Sn,So
      USE GWTBCTMODULE, ONLY: ICT,STOTIO,JPARENT,NPARENT,ISPRCT
      USE GWFDPFMODULE, ONLY: PHIF,SnIM,SoIM
      USE CLN1MODULE, ONLY: ACLNNDS
      USE GWFBASMODULE,ONLY:DELT
      USE GWTDPTMODULE, ONLY: CONCIM,CONCOIM,IADSORBIM,PRSITYIM,
     * ADSORBIM,FLICHIM,IZODIM,IFODIM,ZODRAWIM,FODRAWIM,
     1  ZODRSIM,FODRSIM, SPTLRCTIM 
      USE AW_ADSORBMODULE, ONLY: IKAWI_FN, AK_AWI,ALANGAW,BLANGAW,
     1  AREA_AWI,AREA_AWIO,IAWIZONMAP,NATABROWS,NAZONES,AWI_KAWI_TAB  
      USE AW_ADSORBIMMODULE, ONLY: ALANGAWIM,BLANGAWIM,AK_AWIIM
      DOUBLE PRECISION ADSTERM,CW,FL,CWO,DT,RT,ALENG,VOLU,CEPS,
     1  X,Y,EPSS,EPS,CT,CWEPS,QQ,OMP,tabAWI,STOIT
C     ------------------------------------------------------------------      
        IPAREN = JPARENT(ICOMP,NPAREN)
        IF(ISPRCT.EQ.0)THEN
          STOIT  = STOTIO(ICOMP,NPAREN)
        ELSE
          STOIT = SPTLRCTIM(N,ICOMP,NPAREN)
        ENDIF
        CW = CONCIM(N,IPAREN)
        IF(CW.LT.0.0) CW = 0.0        
C----------------------------------------------------------------------        
C---------ZERO ORDER DECAY ON A-W INTERFACE - APPLY NEWTON EXPANSION OF SUPPLY-DEMAND FUNCTION 
          IF (IZODIM.EQ.4.AND.IKAWI_FN.NE.0) THEN 
C3-------FOR ZERO ORDER DECAY ON AIR-WATER INTERFACE USE LANGMUIR ISOTHERM FILLED AS NEWTON
            CT =  VOLU * ZODRAWIM(N,IPAREN) * AREA_AWI(N) * STOIT
            EPS = 0.01
            CEPS = MAX(0.0,CW)
            X = CEPS /EPS
            CALL SMOOTH(X,Y)
            QQ =  CT * Y
          ENDIF
C -----------------------------------------------------------------------------
C---------FIRST ORDER DECAY ON A-W INTERFACE
          IF(IFODIM.EQ.4.AND.IKAWI_FN.NE.0)THEN
            CT = ALANGAWIM(N,ICOMP) * VOLU * FODRAWIM(N,IPAREN)
     *        * AREA_AWI(N) * STOIT
            IF(IKAWI_FN. EQ.1.OR.IKAWI_FN.EQ.2.OR.IKAWI_FN.EQ. 3) THEN
C--------------FOR FUNCTIONAL CASES
C--------------FOR NON-LINEAR ADSORPTION FILL AS NEWTON
              ETA = BLANGAWIM(N,ICOMP)
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
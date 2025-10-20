C----------------------------------------------------------------------------------------
      SUBROUTINE SCLN2COND1RP
C     ******************************************************************
C      ALLOCATE SPACE AND READ PROPERTIES FOR CONDUIT TYPE CLNs
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE CLN1MODULE
      USE GLOBAL, ONLY: IUNIT,IOUT,NEQS,NODES,NROW,NCOL,IFREFM,IUNSTR,
     *                  INCLN
      DOUBLE PRECISION PERIF,AREAF
      CHARACTER*400 LINE
C----------------------------------------------------------------------------------------
C12------ALLOCATE SPACE FOR CONDUIT TYPE CLNs AND PREPARE TO REFLECT INPUT TO LISTING FILE
      ALLOCATE (ACLNCOND(NCONDUITYP,5))
      ALLOCATE (BHEPROP(NCONDUITYP,4))
      BHEPROP = 0.0
      IF(IBHETYP.EQ.0) THEN 
        WRITE(IOUT,23)
23      FORMAT(/20X,' CONDUIT NODE INFORMATION'/
     1  20X,40('-')/5X,'CONDUIT NODE',8X,'RADIUS',3X,'CONDUIT SAT K',
     1  /5X,12('-'),8X,6('-'),3X,13('-'))
      ELSE
        WRITE(IOUT,24)
24      FORMAT(/20X,' CONDUIT NODE INFORMATION'/
     1  20X,40('-')/5X,'CONDUIT NODE',8X,'RADIUS',3X,'CONDUIT SAT K',
     1  5X,'COND. PIPE', 6X,'THICK PIPE',5X,'COND. FLUID',
     1  5X,'CONV. COEFF',
     1  /5X,12('-'),4X,10('-'),3X,13('-'))          
      ENDIF  
C13------READ CONDUIT PROPERTIES FOR EACH CONDUIT TYPE
      DO I=1,NCONDUITYP
        CALL URDCOM(INCLN,IOUT,LINE)
        IF(IFREFM.EQ.0) THEN
          IF(IBHETYP.EQ.0)THEN
            READ(LINE,'(I10,2F10.3)') IFNO,FRAD,CONDUITK
            LLOC=71
          ELSE 
            READ(LINE,'(I10,6F10.3)')IFNO,FRAD,CONDUITK,TCOND,TTHK,
     *         TCFLUID,TCOEF
            LLOC=111              
          ENDIF    
        ELSE
          LLOC=1
          CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IFNO,R,IOUT,INCLN)
          CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,I,FSRAD,IOUT,INCLN)
          CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,I,CONDUITK,IOUT,INCLN)
          FRAD = FSRAD
          IF(IBHETYP.EQ.1)THEN
            CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,I,TCOND,IOUT,INCLN)  
            CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,I,TTHK,IOUT,INCLN)
            CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,I,TCFLUID,IOUT,INCLN)  
            CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,I,TCOEF,IOUT,INCLN)
          ENDIF  
        END IF
C
C14--------FILL PROPERTY ARRAYS WITH READ AND PREPARE INFORMATION
        ACLNCOND(I,1) = IFNO
        ACLNCOND(I,2) = FRAD
        ACLNCOND(I,3) = CONDUITK
        CALL CLNA(IFNO,AREAF)
        ACLNCOND(I,4) = AREAF
        CALL CLNP(I,PERIF)
        ACLNCOND(I,5) = PERIF
        IF(IBHETYP.EQ.1)THEN
          BHEPROP(I,1) = TCOND
          BHEPROP(I,2) = TTHK
          BHEPROP(I,3) = TCFLUID
          BHEPROP(I,4) = TCOEF
          WRITE(IOUT,26) IFNO,FRAD,CONDUITK,TCOND,TTHK,TCFLUID,TCOEF
26        FORMAT(5X,I10,6(1X,E15.6))          
        ELSE
        WRITE(IOUT,25)IFNO,FRAD,CONDUITK
25      FORMAT(5X,I10,2(1X,E15.6))
        ENDIF          
      ENDDO
C--------RETURN
      RETURN
      END
C----------------------------------------------------------------------------      
      SUBROUTINE SCLN2REC1RP
C     ******************************************************************
C      ALLOCATE SPACE AND READ PROPERTIES FOR RECTAQNGULAR TYPE CLNs
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE CLN1MODULE
      USE GLOBAL, ONLY: IUNIT,IOUT,NEQS,NODES,NROW,NCOL,IFREFM,IUNSTR,
     *                  INCLN
      DOUBLE PRECISION PERIF,AREAF
      CHARACTER*400 LINE
C----------------------------------------------------------------------------------------
C12------ALLOCATE SPACE FOR CONDUIT TYPE CLNs AND PREPARE TO REFLECT INPUT TO LISTING FILE
      ALLOCATE (ACLNREC(NRECTYP,6))
      WRITE(IOUT,23)
23    FORMAT(/20X,' RECTANGULAR CLN SECTION INFORMATION'/
     1  20X,40('-')/5X,'RECTANGULAR NODE',8X,'LENGTH',8X,'HEIGHT',3X,
     1 'CONDUIT SAT K' /5X,12('-'),8X,6('-'),8X,6('-'),3X,13('-'))
C13------READ RECTANGULAR GEOMETRY PROPERTIES FOR EACH RECTANGULAR TYPE
      DO I=1,NRECTYP
        CALL URDCOM(INCLN,IOUT,LINE)
        IF(IFREFM.EQ.0) THEN
          READ(LINE,'(I10,3F10.3)') IFNO,FLENGTH,FHEIGHT,CONDUITK
          LLOC=71
        ELSE
          LLOC=1
          CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IFNO,R,IOUT,INCLN)
          CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,I,FSW,IOUT,INCLN)
          CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,I,FSH,IOUT,INCLN)
          CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,I,CONDUITK,IOUT,INCLN)
          FWIDTH = FSW
          FHEIGHT = FSH
        END IF
C
C14--------FILL PROPERTY ARRAYS WITH READ AND PREPARE INFORMATION
        IFTOTNO = IFNO + NCONDUITYP
        ACLNREC(I,1) = IFTOTNO 
        ACLNREC(I,2) = FWIDTH
        ACLNREC(I,3) = FHEIGHT
        ACLNREC(I,4) = CONDUITK
        CALL CLNA(IFTOTNO,AREAF)
        ACLNREC(I,5) = AREAF
        CALL CLNP(IFTOTNO,PERIF)
        ACLNREC(I,6) = PERIF
        WRITE(IOUT,24)IFNO,FWIDTH,FHEIGHT,CONDUITK
24      FORMAT(5X,I10,3(1X,E15.6))
      ENDDO
C--------RETURN
      RETURN
      END
C----------------------------------------------------------------------------      
CADD      ALLOCATE SPACE AND READ PROPERTIES FOR OTHER TYPES OF CLNs HERE
C----------------------------------------------------------------------------
      SUBROUTINE CLNA(IC,AREAF)
C--------COMPUTE X-SECTIONAL FLOW AREA FOR NODE
      USE CLN1MODULE, ONLY:  ACLNCOND,NCONDUITYP,ACLNREC,NRECTYP
      DOUBLE PRECISION AREAF,RADFSQ
C--------------------------------------------------------------------------------------
      PI = 3.1415926
      IF(IC.LE.NCONDUITYP)THEN
C1-------CLN NODE IS A CONDUIT
        RADFSQ = ACLNCOND(IC,2)**2
        AREAF = PI * RADFSQ
      ELSEIF(IC.GT.NCONDUITYP.AND.IC.LE.NCONDUITYP+NRECTYP)THEN
        ICL = IC - NCONDUITYP
        AREAF = ACLNREC(ICL,2) * ACLNREC(ICL,3)
      ELSEIF(IC.GT.NCONDUITYP+NRECTYP)then   
C2------ADD COMPUTATION FOR AREA FOR OTHER CLN TYPES HERE
CADD      ADD COMPUTATION FOR AREA FOR OTHER TYPES OF CLNs HERE
      ENDIF
C7------RETURN
      RETURN
      END
C--------------------------------------------------------------------------------------
      SUBROUTINE CLNAGET(IC,AREAF)
C--------GET X-SECTIONAL FLOW AREA FOR NODE FROM RESPECTIVRE ARRAY
      USE CLN1MODULE, ONLY:  ACLNCOND,NCONDUITYP,ACLNREC,NRECTYP
      DOUBLE PRECISION AREAF,RADFSQ
C--------------------------------------------------------------------------------------
      IF(IC.LE.NCONDUITYP)THEN
C1-------CLN NODE IS A CONDUIT
        AREAF = ACLNCOND(IC,4)
      ELSEIF(IC.GT.NCONDUITYP.AND.IC.LE.NCONDUITYP+NRECTYP)THEN 
        ICL = IC - NCONDUITYP
        AREAF = ACLNREC(ICL,5)
      ELSEIF(IC.GT.NCONDUITYP+NRECTYP)then   
C2------ADD COMPUTATION FOR AREA FOR OTHER CLN TYPES HERE
CADD      ADD COMPUTATION FOR AREA FOR OTHER TYPES OF CLNs HERE
      ENDIF
C7------RETURN
      RETURN
      END
C--------------------------------------------------------------------------------------
      SUBROUTINE CLNK(IC1,NC1,IC2,NC2,FK)
C--------COMPUTE EFFECTIVE LEAKANCE (CONSTANT TERM) FOR NODE
      USE CLN1MODULE, ONLY:  ACLNCOND,NCONDUITYP,ACLNNDS,GRAV,VISK,
     1  IFLINCLN,ACLNREC,NRECTYP
      DOUBLE PRECISION FK,RADFSQ1,RADFSQ2,AR1,AR2,PERI1,PERI2,
     1  EL,CK1,CK2,CK
C--------------------------------------------------------------------------------------
C1 ---GET LOCAL VARIABLES FOR AREA, PERIMETER AND CONDUCTANCE TERM FOR EACH NODE      
      EL = (ACLNNDS(NC1,4) + ACLNNDS(NC2,4)) * 0.5  ! CONNECTION LENGTH IS HALF OF CLN CELL LENGTH
      CALL GETAPK12(IC1,IC2,AR1,AR2,PERI1,PERI2,CK1,CK2)
      IFLIN = IABS(IFLINCLN(NC1))
C ------------------------------------------------------------------------------------------
      IF(IFLIN.EQ.1)THEN 
C1-------CLN NODE IS LAMINAR FLOW USING HAGEN-POISEUILLE EQUATION          
        IF(VISK.LT.1.0E-10.AND.GRAV.LT.1.0E-10)THEN
C --------THIS PART IS FOR BACKWARD COMPATIBILITY WHERE CONDUITK HAD THE GRAV AND VISCOSITY TERM
          CK1 = CK1 * 4.0 * (AR1/PERI1)**2
          CK2 = CK2 * 4.0 * (AR2/PERI2)**2
        ELSE
          CK1 = 0.5*GRAV/VISK*(AR1/PERI1)**2
          CK2 = 0.5*GRAV/VISK*(AR2/PERI2)**2
        ENDIF  
        CK = 2.0 *CK1*CK2 / (CK1 + CK2)         
        FK = CK / EL  
      ELSEIF(IFLIN.EQ.2) THEN  
C2-------CLN NODE IS TURBULENT FLOW USING DARCY-WEISBACH EQUATION 
        FK = SQRT(32.0*GRAV/EL)       
      ELSEIF(IFLIN.EQ.3) THEN             
C3-------CLN NODE IS TURBULENT FLOW USING HAZEN-WILLIAMS EQUATION
        CK = 2.0 *CK1*CK2 / (CK1 + CK2)
        FK = 0.849 * CK / EL**0.54    
      ELSEIF(IFLIN.EQ.4) THEN  
C4-------CLN NODE IS TURBULENT FLOW USING MANNINGS EQUATION
        CK = 2.0 *CK1*CK2 / (CK1 + CK2)
        FK = 1.0/CK/SQRT(EL)   
      ENDIF  
C7------RETURN
      RETURN
      END
C ------------------------------------------------------------------------------      
      SUBROUTINE GETAPK12(IC1,IC2,AR1,AR2,PERI1,PERI2,CK1,CK2)
C     ******************************************************************
C--------RETRIEVE AREA, PERIMETER AND CONDUCTANCE TERM FOR VARIOUS GEOMETRIES
C     ******************************************************************      
      USE CLN1MODULE, ONLY:  ACLNCOND,NCONDUITYP,ACLNNDS,GRAV,VISK,
     1  IFLINCLN,ACLNREC,NRECTYP
      DOUBLE PRECISION AR1,AR2,PERI1,PERI2,CK1,CK2      
      
C1A ---FOR NODE 1      
      IF(IC1.LE.NCONDUITYP)THEN
C1A1---------CLN1 NODE IS A CONDUIT          
        AR1 = ACLNCOND(IC1,4) 
        PERI1 = ACLNCOND(IC1,5)   
        CK1 = ACLNCOND(IC1,3)
      ELSEIF(IC1.GT.NCONDUITYP.AND.IC1.LE.NCONDUITYP+NRECTYP)THEN
C1A2---------CLN1 NODE IS A RECTANGULAR SECTION
        ILC1 = IC1 - NCONDUITYP
        AR1 = ACLNREC(ILC1,5) 
        PERI1 = ACLNREC(ILC1,6)  
        CK1 = ACLNREC(ILC1,4)
      ELSEIF(IC1.GT.NCONDUITYP+NRECTYP)THEN
C1A3--------ADD COMPUTATION FOR LAMINAR K FOR OTHER CLN TYPES HERE
CADD--     ADD COMPUTATION FOR LAMINAR K FOR OTHER CLN TYPES HERE        
      ENDIF
C1B ---FOR NODE2      
      IF(IC2.LE.NCONDUITYP)THEN
C1B1---------CLN1 NODE IS A CONDUIT          
        AR2 = ACLNCOND(IC2,4) 
        PERI2 = ACLNCOND(IC2,5)         
        CK2 = ACLNCOND(IC2,3)
      ELSEIF(IC1.GT.NCONDUITYP.AND.IC1.LE.NCONDUITYP+NRECTYP)THEN
C1B2---------CLN1 NODE IS A RECTANGULAR SECTION
        ILC2 = IC2 - NCONDUITYP
        AR2 = ACLNREC(ILC2,5) 
        PERI2 = ACLNREC(ILC2,6)  
        CK2 = ACLNREC(ILC2,4)
      ELSEIF(IC1.GT.NCONDUITYP+NRECTYP)THEN
C1B3--------ADD COMPUTATION FOR LAMINAR K FOR OTHER CLN TYPES HERE
CADD--     ADD COMPUTATION FOR LAMINAR K FOR OTHER CLN TYPES HERE          
      ENDIF  
C7------RETURN
      RETURN
      END      
C------------------------------------------------------------------------
      SUBROUTINE TURBFUNC(NC1,NC2,IFLIN,HD,TURB)
C     ******************************************************************
C     COMPUTE THE TURBULENT GRADIENT FUNCTION TERM FOR THE VARIOUS FORMULATIONS.
C     ******************************************************************
      USE CLN1MODULE, ONLY:  ACLNCOND,NCONDUITYP,ACLNNDS,GRAV,IFLINCLN,
     1   VISK
      USE GLOBAL, ONLY: HNEW
      DOUBLE PRECISION TURB,FK,RADFSQ,AR,PERI,EL,CONDUITK,HD,ALPHA,BETA,
     1  SQTERM,EPCLIP,PETAM1,PETAM3,PETA,AREPB,APBORE,C1,C2,
     1  AR1,AR2,PERI1,PERI2,CK1,CK2 
      DATA EPCLIP/ 0.01 /      
ccsp      DATA EPCLIP/ 0.1 /      
C-----------------------------------------------------------------------------      
C
C1------RETURN FOR LINEAR FLOW      
      IF(IFLIN.EQ.1) RETURN
C-----------------------------------------------------------------------------  
C2------FOR DARCY-WEISBACH EQUATION
      IF(IFLIN.EQ.2)THEN
C2A-----COMPUTE TERMS ALPHA AND BETA WITHIN THE LOG
        EL = (ACLNNDS(NC1,4)+ACLNNDS(NC2,4)) * 0.5  ! CONNECTION LENGTH IS HALF OF CLN CELL LENGTH SUMS    
        IC1 = ACLNNDS(NC1,2)
        IC2 = ACLNNDS(NC2,2)
        CALL GETAPK12(IC1,IC2,AR1,AR2,PERI1,PERI2,CK1,CK2)
        AR = (AR1 + AR2) * 0.5   ! AVERAGE AREA OF CONNECTED CELLS
        PERI = (PERI1 + PERI2) * 0.5 ! AVERAGE PERIMETER OF CONNECTED CELLS
        CONDUITK = (CK1 + CK2) * 0.5 ! AVERAGE ROUGHNESS OF CONNECTED CELLS
        ALPHA = CONDUITK * PERI / (14.84 * AR)
        SQTERM = (AR/PERI)**3 * GRAV / EL
        BETA = 0.222 * VISK / SQRT(SQTERM)          
        IF(HD.GT.EPCLIP)THEN
          TURB = HD**(-0.5) * LOG10 (ALPHA + BETA/SQRT(HD))
        ELSE
C2B-------FLATTEN AND SMOOTH THE FUNCTION BEFORE IT GROWS VERY LARGE  
          AREPB = ALPHA * SQRT(EPCLIP) + BETA
          APBORE = ALPHA + BETA/SQRT(EPCLIP)
          C1 = -AREPB * LOG(APBORE) - BETA
          C1 = C1 / (9.21034 * EPCLIP**2.5 * AREPB)
          C2 = LOG10(APBORE)/SQRT(EPCLIP) - C1*EPCLIP**2
          TURB = C1*HD**2 + C2
        ENDIF
C-------------------------------------------------------------------------------        
      ELSEIF(IFLIN.EQ.3.OR.IFLIN.EQ.4)THEN    
C3------FOR POWER-OF-HD STYLE EQUATIONS
C-------------------------------------------------------------------------------          
        IF(IFLIN.EQ.3)THEN 
C3A------FOR HAZEN-WILLIAMS EQUATION            
          PETA = 0.54 
        ELSEIF(IFLIN.EQ.4)THEN    
C3B------FOR MANNINGS EQUATION
          PETA = 0.5
        ENDIF
C--------------------------------------------------------------------
        PETAM1 = PETA-1.0
        IF(HD.GT.EPCLIP)THEN
          TURB = HD**PETAM1
        ELSE
C3C-------FLATTEN AND SMOOTH THE FUNCTION BEFORE IT GROWS VERY LARGE
         PETAM3 = PETA - 3.0   
         TURB= 0.5*(PETAM1*EPCLIP**PETAM3*HD**2 - EPCLIP**PETAM1*PETAM3)
C3C-------straighten FUNCTION BEFORE IT GROWS VERY LARGE
csp         PETAM2 = PETA - 2.0   
csp         TURB= PETAM1 * EPCLIP**PETAM2 * HD - PETAM2 * EPCLIP**PETAM1 
        ENDIF
      ENDIF
      TURB = ABS(TURB)
C
C4-----RETURN.
      RETURN
      END
C--------------------------------------------------------------------------------------
      SUBROUTINE CLNR(IC,FRAD)
C--------COMPUTE RADIUS FOR CONNECTION OF CLN SEGMENT TO 3-D GRID WITH THEIM EQUATION
      USE CLN1MODULE, ONLY: ACLNCOND,NCONDUITYP,ACLNREC,NRECTYP
      DOUBLE PRECISION FRAD
C--------------------------------------------------------------------------------------
      IF(IC.LE.NCONDUITYP)THEN
C1-------CLN NODE IS A CONDUIT
        FRAD = ACLNCOND(IC,2)
      ELSEIF(IC.GT.NCONDUITYP.AND.IC.LE.NCONDUITYP+NRECTYP)THEN
C1B---------CLN1 NODE IS A RECTANGULAR SECTION: FRAD IS HALF OF HYDRAULIC RADIUS
        ICL = IC - NCONDUITYP
        FRAD = 0.5 * ACLNREC(ICL,4)/ACLNREC(ICL,5)
      ELSEIF(IC.GT.NCONDUITYP+NRECTYP)THEN
C2------ADD COMPUTATION FOR EFFECTIVE RADIUS FOR OTHER CLN TYPES HERE
CADD     ADD COMPUTATION FOR RADIUS FOR OTHER CLN TYPES HERE
      ENDIF
C7------RETURN
      RETURN
      END
C--------------------------------------------------------------------------------------
      SUBROUTINE CLNP(IC,FPER)
C--------COMPUTE EFFECTIVE PERIMETER FOR CONNECTION OF CLN SEGMENT TO 3-D GRID
      USE CLN1MODULE, ONLY: ACLNCOND,NCONDUITYP,ACLNREC,NRECTYP
      DOUBLE PRECISION FPER
C--------------------------------------------------------------------------------------
      PI = 3.1415926
      IF(IC.LE.NCONDUITYP)THEN
C1-------CLN NODE IS A CONDUIT
        FPER = 2 * PI * ACLNCOND(IC,2)
      ELSEIF(IC.GT.NCONDUITYP.AND.IC.LE.NCONDUITYP+NRECTYP)THEN
C1B---------CLN1 NODE IS A RECTANGULAR SECTION
        ICL = IC - NCONDUITYP
        FPER = 2 * (ACLNREC(ICL,2) + ACLNREC(ICL,3))
      ELSEIF(IC.GT.NCONDUITYP+NRECTYP)THEN
C2------ADD COMPUTATION FOR PERIMETER FOR OTHER CLN TYPES HERE
CADD      ADD COMPUTATION FOR PERIMETER FOR OTHER CLN TYPES HERE      
      ENDIF
C7------RETURN
      RETURN
      END
C--------------------------------------------------------------------------------
      SUBROUTINE CLNPGET(IC,FPER)
C--------GET EFFECTIVE TOTAL PERIMETER FOR CONNECTION OF CLN SEGMENT TO 3-D GRID FROM ARRAYS
      USE CLN1MODULE, ONLY: ACLNCOND,NCONDUITYP,ACLNREC,NRECTYP
      DOUBLE PRECISION FPER
C--------------------------------------------------------------------------------------
      IF(IC.LE.NCONDUITYP)THEN
C1-------CLN NODE IS A CONDUIT
        FPER = ACLNCOND(IC,5)
      ELSEIF(IC.GT.NCONDUITYP.AND.IC.LE.NCONDUITYP+NRECTYP)THEN
C1B---------CLN1 NODE IS A RECTANGULAR SECTION
        ICL = IC - NCONDUITYP
        FPER = ACLNREC(ICL,6) 
      ELSEIF(IC.GT.NCONDUITYP+NRECTYP)THEN
C2------ADD COMPUTATION FOR PERIMETER FOR OTHER CLN TYPES HERE
CADD      ADD COMPUTATION FOR PERIMETER FOR OTHER CLN TYPES HERE      
      ENDIF
C7------RETURN
      RETURN
      END 
C-------------------------------------------------------------------------------
      SUBROUTINE CLNPW(ICLN,HD,PERIW,IGWCLN)
C--------COMPUTE WETTED X-SECTIONAL PERIMETER FOR HORIZONTAL CLN CELL
      USE CLN1MODULE, ONLY:  ACLNCOND,NCONDUITYP,ACLNNDS,ACLNREC,NRECTYP
      USE GLOBAL, ONLY: NODES,HNEW
      DOUBLE PRECISION PERIW,RADFSQ,HD,FRAD,BBOT,DEPTH,
     1  AFAC,BFAC,CFAC,DTOP,EFAC
      INTEGER IGWCLN ! INDEX OF CONNECTION = 0 FOR CLN-CLN; =1 FOR GW-CLN
C--------------------------------------------------------------------------------------
      N = ACLNNDS(ICLN,1)
      IC = ACLNNDS(ICLN,2)
      BBOT = ACLNNDS(ICLN,5)
      DEPTH = HD - BBOT      
      IF(IC.LE.NCONDUITYP)THEN
C1-------CLN NODE IS A CONDUIT
        PI = 3.1415926
        FRAD = ACLNCOND(IC,2)
        IF(DEPTH.LE.0)THEN
          PERIW = 0.0
        ELSEIF(DEPTH.LE.FRAD)THEN
          PERIW = 2.0*FRAD*ACOS((FRAD-DEPTH)/FRAD)
        ELSEIF(DEPTH.LE.2.0*FRAD)THEN
          PERIW = 2.0*FRAD*(PI - ACOS((DEPTH-FRAD)/FRAD))
        ELSE
          PERIW = 2* PI *FRAD
        ENDIF
      ELSEIF(IC.GT.NCONDUITYP.AND.IC.LE.NCONDUITYP+NRECTYP)THEN
C1 -------CLN NODE IS A RECTANGULAR SECTION
        ICL = IC - NCONDUITYP
        HEIGHT = ACLNREC(ICL,3)
        WIDTH = ACLNREC(ICL,2)
        EPSILON = 0.01 
        IF(DEPTH.LE.0.0) DEPTH = 0.0
        IF(DEPTH.LT.HEIGHT)THEN
          PERIW = ACLNREC(ICL,2) + 2.0* DEPTH
cc          PERIW = ACLNREC(ICL,2)   ! Wide rectangle assumption *************************************************************
        ELSEIF(DEPTH.LT.(HEIGHT + EPSILON))THEN 
          DTOP = DEPTH - HEIGHT
          EFAC = WIDTH + 2*HEIGHT
          CFAC = 2.0
          BFAC = (3*WIDTH - 4 * EPSILON) / EPSILON**2
          AFAC = -(2.0*BFAC*EPSILON + 2.0) / (3 * EPSILON**2)
          PERIW = AFAC*DTOP**3 + BFAC*DTOP**2 + CFAC*DTOP + EFAC
        ELSE    
          PERIW = ACLNREC(ICL,2) + 2.0 * HEIGHT   + WIDTH  
        ENDIF   
C ------FOR GW-CLN CONNECTION, NEED TO SMOOTHEN RECTANGLE TO ZERO WHEN DEPTH COMES TO ZERO
        IF(IGWCLN.EQ.1)THEN
C5----------PROVIDE SMOOTH PATCHING FUNCTION WHEN DEPTH DROPS TO ZERO
          EPSILON = 0.01
          FPATCH = 1.0D0
          HR = DEPTH
          IF(HR.LT.EPSILON)THEN
            FPATCH = 2.0D0*(HR/EPSILON)**3 + 3.0D0*(HR/EPSILON)**2 
            PERIW = PERIW * FPATCH 
          ENDIF            
        ENDIF
      ELSEIF(IC.GT.NCONDUITYP)THEN
C2------ADD COMPUTATION FOR WETTED X-SECTIONAL PERIMETER FOR OTHER CLN TYPES HERE
CADD      ADD COMPUTATION FOR WETTED X-SECTIONAL PERIMETER FOR OTHER CLN TYPES HERE    
      ENDIF
C
C5------RETURN.
      RETURN
      END
C-------------------------------------------------------------------------------
      SUBROUTINE CLNAW(ICLN,HD,AREAW)
C--------COMPUTE WETTED X-SECTIONAL FLOW AREA FOR NODE
      USE CLN1MODULE, ONLY:  ACLNCOND,NCONDUITYP,ACLNNDS,ACLNREC,NRECTYP
      USE GLOBAL, ONLY: NODES,HNEW
      DOUBLE PRECISION AREAW,RADFSQ,HD,FRAD,BBOT,DEPTH
C--------------------------------------------------------------------------------------
      N = ACLNNDS(ICLN,1)
      IC = ACLNNDS(ICLN,2)
      BBOT = ACLNNDS(ICLN,5)
      DEPTH = HD - BBOT
      IF(IC.LE.NCONDUITYP)THEN
C1-------CLN NODE IS A CONDUIT
        PI = 3.1415926
        FRAD = ACLNCOND(IC,2)
        IF(DEPTH.LE.0)THEN
          AREAW = 0.0
        ELSEIF(DEPTH.LE.FRAD)THEN
          AREAW = FRAD*FRAD*ACOS((FRAD-DEPTH)/FRAD) - (FRAD-DEPTH)*
     1         SQRT(FRAD*FRAD - (FRAD-DEPTH)**2)
        ELSEIF(DEPTH.LE.2.0*FRAD)THEN
          AREAW = FRAD*FRAD*(PI - ACOS((DEPTH-FRAD)/FRAD))
     1    - (FRAD-DEPTH) * SQRT(FRAD*FRAD - (FRAD-DEPTH)**2)
        ELSE
          AREAW = PI *FRAD*FRAD
        ENDIF
      ELSEIF(IC.GT.NCONDUITYP.AND.IC.LE.NCONDUITYP+NRECTYP)THEN
C1 -------CLN NODE IS A RECTANGULAR SECTION
        ICL = IC - NCONDUITYP
        HEIGHT = ACLNREC(ICL,3)
        IF(DEPTH.LE.0)THEN
          AREAW = 0.0  
        ELSEIF(DEPTH.LT.HEIGHT)THEN
          AREAW = ACLNREC(ICL,2) * DEPTH
        ELSE
          AREAW = ACLNREC(ICL,2) * HEIGHT   
        ENDIF 
      ELSEIF(IC.GT.NCONDUITYP+NRECTYP)THEN
C2------ADD COMPUTATION FOR WETTED X-SECTIONAL FLOW AREA FOR OTHER CLN TYPES HERE
CADD      ADD COMPUTATION FOR WETTED X-SECTIONAL FLOW AREA FOR OTHER CLN TYPES HERE
      ENDIF
C
C5------RETURN.
      RETURN
      END
C-----------------------------------------------------------------------


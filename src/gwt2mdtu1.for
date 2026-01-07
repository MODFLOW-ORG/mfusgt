      module GWTmdTMODULE
c
c      new input variables for USG-Transport 6/7/20
c----------these are variable by element in USG-Transport, but are constant
c          here except for decay, which varies in the space/time zones
c          keep this simple for now
        integer, save, pointer :: IMDTCB,IMDTCF,imdtcf2,IMULTIMD
        integer, save, allocatable, dimension (:) ::
     1    mdflag, mdtop,mdbot,mdembed
        real, save, allocatable, dimension (:) :: VOLFRACMD, DIFFLENMD,
     &     PORMD, TORTMD, RHOBMD 
        real*8, save, allocatable, dimension (:) :: amd1,amd2
         real, save, allocatable, dimension (:,:) :: KDMD,DECAYMD,
     1    YIELDMD,DIFFMD
         real*8, save, allocatable, dimension (:,:) :: alpha,retardm,
     1    alambdar,aiold1,aiold2,ainew1,ainew2,d,delt1,delt2,gam1,gam2,
     1    beta1,beta2,aa1,aa2,bb1,bb2,p1,p2,q1,q2,  cn,co         
         DOUBLE PRECISION, SAVE, POINTER :: TSHIFTMD
         CHARACTER (LEN=200), SAVE :: CROOTNAME_MD
c
      end module GWTmdTMODULE
C-------------------------------------------------------------------------------
      subroutine gwt2mdtu1ar(IN,NPER)
C     ******************************************************************
C     INITIALIZE VARIABLES AND READ DATA FOR matrix duffusion
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL, ONLY:IFREFM,IOUT,NODES,NLAY,NODLAY,ImdT,IDPF,
     1    IBOUND,TOP,BOT,AREA,ISSFLG,PGF,IA,JAS,JA,Sn,ISSFLG
      use GWTmdTMODULE
      USE GWTBCTMODULE, ONLY: MCOMP,ICBUND,ZODRW,ZODRS,FODRW,FODRS,
     1  PRSITY,ADSORB,IZOD,IFOD,IADSORB,CINACT,MCOMPT,NTCOMP,IHEAT,
     1  ICHAIN,ISPRCT,NPARENT
      USE GWFBCFMODULE, ONLY: SC1,SC2,HK,LAYCON
      USE GWFDPFMODULE, ONLY: PHIF,SnIM,SoIM,SC1IM,SC2IM,HNEWIM
      CHARACTER*200 LINE
      CHARACTER*16 TEXT, ISTATUS
      CHARACTER*400 CFILENAME
      DOUBLE PRECISION MASLOCW,MASLOCS,MASLOCC,REFHD
      INTEGER IFRAHK
C
      REAL, DIMENSION(:),ALLOCATABLE  ::TEMPC
      CHARACTER*24 ANAME(12),CNAME
      DATA ANAME(1) /'    mdT MATRIX TYPE FLAG'/
      DATA ANAME(2) /'mdT MATRIX VOLUME FRACTN'/
      DATA ANAME(3) /'     mdT matrix POROSITY'/
      DATA ANAME(4) /'        mdT BULK DENSITY'/
      DATA ANAME(5) /'    mdT DIFFUSION LENGTH'/
      DATA ANAME(6) /'   mdT MATRIX TORTUOSITY'/
      DATA ANAME(7) /'     mdT ADSORPTION COEF'/
      DATA ANAME(8) /'        mdT MATRIX DECAY'/
      DATA ANAME(9) /'  mdT MATRIX YIELD COEFF'/
      DATA ANAME(10) /' mdT MATRIX DIFFSN COEFF'/
      DATA ANAME(11) /'  mdT MATRIX AI1MD COEFF'/
      DATA ANAME(12) /'  mdT MATRIX AI2MD COEFF'/
      INCLUDE 'openspec.inc'
C     ------------------------------------------------------------------
C
C1------IDENTIFY PACKAGE
      ImdT = 1 ! FLAG FOR DUAL POROSITY TRANSPORT IS ON
      WRITE(IOUT,1) IN
    1 FORMAT(1X,/1X,'mdT -- MATRIX DIFFUSION TRANSPORT PACKAGE, ',
     1'VERSION 1, 06/26/2020',/,9X,'INPUT READ FROM UNIT',I3)
C
C-----------------------------------------------------------------------
C2------ALLOCATE VARIABLES
      ALLOCATE(IMDTCB,IMDTCF,IMDTCF2,IMULTIMD)
      ALLOCATE(TSHIFTMD)
      TSHIFTMD = 0.0
      IMULTIMD = 0
C
C-----------------------------------------------------------------------
C3------READ GENERAL MDT INFORMATION AND FLAGS
      IF(IFREFM.EQ.0) THEN
        READ(IN,2)imdtcb,IMDTCF
        LLOC=11
      ELSE
        CALL URDCOM(IN,IOUT,LINE)
        LLOC=1
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IMDTCB,R,IOUT,IN)
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IMDTCF,R,IOUT,IN)
      ENDIF
2     FORMAT(I10)
      IMDTCF2 = IMDTCF
C
C3B-----PRINT VALUES
C
      IF(ImdTCB.GT.0) WRITE(IOUT,9) ImdTCB
    9 FORMAT(1X,'MATRIX DIFFUSION CELL-BY-CELL MASS FLUX WILL BE SAVED',
     1 1X,'ON UNIT ImdTCB  =',I3)
C      
      IF(ImdTCF.GT.0)WRITE(IOUT,11) imdtcf 
   11 FORMAT(1X,'MD COEFFICIENTS WILL BE SAVED ON UNIT ImdTCF = ',I3)  
C-------------------------------------------------------------------------------
C3C--------GET OPTIONS FOR WHEN FLOW IS NOT DUAL POROSITY
      IF(IDPF.EQ.0)THEN
        IFRAHK = 0
100     CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,I,R,IOUT,IN)
        IF(LINE(ISTART:ISTOP).EQ.'FRAHK') THEN
C3D-------SET FLAG FOR HK TO BE FRACTURE-DOMAIN VALUE.
          IFRAHK = 1
        ELSEIF(LINE(ISTART:ISTOP).EQ.'FRADARCY') THEN
C3D-------SET FLAG FOR HK TO BE FRACTURE-DOMAIN VALUE.
          IFRAHK = 2  
        ELSEIF(LINE(ISTART:ISTOP).EQ.'TSHIFTMD') THEN
C3D-----READ KEYWORD OPTION FOR NON-ZERO STARTING TIME OF MD
          CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,I,TS,IOUT,IN)
          TSHIFTMD = TS
          WRITE(IOUT,37) TS
37        FORMAT(1X,'TIME SHIFT AMOUNT (TSHIFTMD) = ', F10.3)
        ELSEIF(LINE(ISTART:ISTOP).EQ.'SEPARATE_AI2') THEN
C3D-----READ KEYWORD OPTION FOR NON-ZERO STARTING TIME OF MD
          CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IMDTCF2,R,IOUT,IN)
          WRITE(IOUT,38) IMDTCF2
38        FORMAT(1X,'MD COEFFICIENT AI2 WILL BE SAVED SEPARATELY ON',
     1       1X,'UNIT ImdTCF2 (SEPARATE_AI2) = ',I3)           
        ELSEIF(LINE(ISTART:ISTOP).EQ.'MULTIFILE_MD') THEN
C3B------READ KEYWORD OPTION TO OUTPUT BINARY CONC IN SEPARATE FILES.
          IMULTIMD = 1
          CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,N,R,IOUT,IN)
          CROOTNAME_MD = LINE(ISTART:ISTOP)
          WRITE(IOUT,39) CROOTNAME_MD
39        FORMAT(1X,'ROOT FILE NAME FOR AI1 & AI2 OUTPUT = ', A20) 
        ENDIF
        IF(LLOC.LT.200) GO TO 100
C3H----------PRINT OPTION FLAG VALUES
        IF(IFRAHK.EQ.1) WRITE(IOUT,12) IFRAHK
12      FORMAT(1X,'CONDUCTANCE ARE FOR FRACTURE VOLUME; (IFRAHK)',12X,
     1   '=',I4)
        IF(IFRAHK.EQ.0) WRITE(IOUT,13) IFRAHK
13      FORMAT(1X,'CONDUCTANCE ARE FOR TOTAL (FRACTURE + MATRIX)',
     1            1X,'DOMAIN; (IFRAHK)  =',I4)
        ENDIF
C---------------------------------------------------------------------------------------
C4i----------FOR MULTIFILE OPEN FILES
        IF(IMULTIMD.NE.0)THEN 
          IF(IMDTCF. EQ.0) THEN
            WRITE(IOUT,*)'NEED POSITIVE VALUE FOR IMDTCF IN MDT INPUT'
            STOP
          ENDIF  
          DO ICOMP = 1,NTCOMP
            CALL GET_TEXTAI1(ICOMP,TEXT)
            CFILENAME = trim(CROOTNAME_MD) // TRIM(TEXT) //'.AI1'
!            OPEN(ISPCUN, FILE = CFILENAME, FORM = 'BINARY',  
            OPEN(IMDTCF, FILE = CFILENAME, FORM = FORM,     ! - use FORM from openspec.inc rather than hardwiring BINARY or UNFORMATTED
     1       STATUS = 'REPLACE')
            CLOSE (IMDTCF)
          ENDDO
          IF(IMDTCF2.NE.0)THEN
C4ii--------OPEN SEPARATE FILE FOR AI2              
            DO ICOMP = 1,NTCOMP
              CALL GET_TEXTAI2(ICOMP,TEXT)
              CFILENAME = trim(CROOTNAME_MD) // TRIM(TEXT) //'.AI2'
!            OPEN(ISPCUN, FILE = CFILENAME, FORM = 'BINARY',  
              OPEN(IMDTCF, FILE = CFILENAME, FORM = FORM,     ! - use FORM from openspec.inc rather than hardwiring BINARY or UNFORMATTED
     1         STATUS = 'REPLACE')
              CLOSE (IMDTCF)
            ENDDO
          ENDIF  
        ENDIF
C-----------------------------------------------------------------------------------
C4-----ALLOCATE ARRAYS AND INITIALIZE
c
      allocate (mdflag(nodes),mdtop(nodes),mdbot(nodes),mdembed(nodes))
      allocate (VOLFRACMD(nodes), DIFFLENMD(nodes), PORMD(nodes),
     &  TORTMD(nodes), RHOBMD(nodes), amd1(nodes),amd2(nodes),
     &  KDMD(nodes,mcompt),DECAYMD(nodes,mcompt),YIELDMD(nodes,mcompt),
     &  DIFFMD(nodes,mcompt),alpha(nodes,mcompt),retardm(nodes,mcompt),
     &  alambdar(nodes,mcompt),
     &  aiold1(nodes,mcompt),aiold2(nodes,mcompt),ainew1(nodes,mcompt),
     &  ainew2(nodes,mcompt),d(nodes,mcompt),delt1(nodes,mcompt),
     &  delt2(nodes,mcompt),gam1(nodes,mcompt),gam2(nodes,mcompt),
     &  beta1(nodes,mcompt),beta2(nodes,mcompt),aa1(nodes,mcompt),
     &  aa2(nodes,mcompt),bb1(nodes,mcompt),bb2(nodes,mcompt),
     &  p1(nodes,mcompt),p2(nodes,mcompt),q1(nodes,mcompt),
     &  q2(nodes,mcompt) ) 
      ainew1 = 0.
      ainew2 = 0. 
C-----------------------------------------------------------------------
C5A-----READ mdflag ARRAY: 
c     7 combinations of matrix diffusion: no md, md on bottom, embedded md
c     md on top, md on top and bottom, embedded md and on bottom, embedded md and on top, and
c     embedded md and md on top and bottom.  These correspond to mdflags 1-7, respectively
c      
        DO K = 1,NLAY
          KK = K
          NNDLAY = NODLAY(K)
          NSTRT = NODLAY(K-1)+1
          NDSLAY = NNDLAY - NODLAY(K-1)
          CALL U1DINT(mdflag(NSTRT),ANAME(1),NDSLAY,K,IN,IOUT)
        ENDDO
C-----------------------------------------------------------------------
C6------READ MATERIAL PROPERTIES FOR MATRIX DOMAIN.
C-----------------------------------------------------------------------
C6A-----READ HIGH k MATERIAL VOLUME FRACTION INTO ARRAY VOLFRACMDMD IF FLOW WAS SINGLE DOMAIN
      IF(IDPF.EQ.0)THEN
        DO K = 1,NLAY
          KK = K
          NNDLAY = NODLAY(K)
          NSTRT = NODLAY(K-1)+1
          NDSLAY = NNDLAY - NODLAY(K-1)
          CALL U1DREL(VOLFRACMD(NSTRT),ANAME(2),NDSLAY,K,IN,IOUT)
        ENDDO
      ELSE
C6B-----COMPUTE VOLFRACMDMD FROM PHIF READ IN DURING DPF
        DO N = 1,NODES
          VOLFRACMD(N) = PHIF(M)
        ENDDO
      ENDIF
c
C6C-----READ MATRIX DOMAIN POROSITY INTO ARRAY PORMD
      DO K = 1,NLAY
        NNDLAY = NODLAY(K)
        NSTRT = NODLAY(K-1)+1
        NDSLAY = NNDLAY - NODLAY(K-1)
        CALL U1DREL(PORMD(NSTRT),ANAME(3),NDSLAY,K,IN,IOUT)
      ENDDO
c
C6C-----READ MATRIX DOMAIN BULK DENSITY INTO ARRAY RHOBMD
      DO K = 1,NLAY
        NNDLAY = NODLAY(K)
        NSTRT = NODLAY(K-1)+1
        NDSLAY = NNDLAY - NODLAY(K-1)
        CALL U1DREL(RHOBMD(NSTRT),ANAME(4),NDSLAY,K,IN,IOUT)
      ENDDO
c
C6C-----READ MATRIX DOMAIN DIFFUSION LENGTH INTO ARRAY DIFFLENMD
      DO K = 1,NLAY
        NNDLAY = NODLAY(K)
        NSTRT = NODLAY(K-1)+1
        NDSLAY = NNDLAY - NODLAY(K-1)
        CALL U1DREL(DIFFLENMD(NSTRT),ANAME(5),NDSLAY,K,IN,IOUT)
      ENDDO
      do n=1,nodes
        if(difflenmd(n).lt.1.0e-15) difflenmd(n) = 1.0e-15 
      enddo 
c
C6C-----READ MATRIX DOMAIN TORTUOSITY INTO ARRAY TORTMD
      DO K = 1,NLAY
        NNDLAY = NODLAY(K)
        NSTRT = NODLAY(K-1)+1
        NDSLAY = NNDLAY - NODLAY(K-1)
        CALL U1DREL(TORTMD(NSTRT),ANAME(6),NDSLAY,K,IN,IOUT)
      ENDDO
C ---------------------------------------------------------------------
C7------READ SPECIES DEPENDENT ARRAYS FOR ALL LAYERS.
C-----------------------------------------------------------------------
      ALLOCATE(TEMPC(Nodes))
      DO ICOMP=1,NTCOMP
        WRITE(IOUT,20) ICOMP
20      FORMAT(80('-')/1X,'THE FOLLOWING ARRAYS ARE READ FOR SPECIES ',
     *    1X,'NUMBER ',I3/80('-'))
C
C7A-----READ MATRIX DOMAIN ADSORPTION COEFFICIENT INTO ARRAY KDMD
        DO K=1,NLAY
          NNDLAY = NODLAY(K)
          NSTRT = NODLAY(K-1)+1
          NDSLAY = NNDLAY - NODLAY(K-1)
          CALL U1DREL(TEMPC(NSTRT),ANAME(7),NDSLAY,K,IN,IOUT)
        ENDDO
C
        DO N = 1, NODES
          KDMD(N,ICOMP) = TEMPC(N)
        ENDDO
C
C7B-----READ MATRIX DOMAIN DECAY INTO ARRAY DECAYMD
        DO K=1,NLAY
          NNDLAY = NODLAY(K)
          NSTRT = NODLAY(K-1)+1
          NDSLAY = NNDLAY - NODLAY(K-1)
          CALL U1DREL(TEMPC(NSTRT),ANAME(8),NDSLAY,K,IN,IOUT)
        ENDDO
C
        DO N = 1, NODES
          DECAYMD(N,ICOMP) = TEMPC(N)
          if(decaymd(n,icomp).lt.1.0e-15) decaymd(n,icomp) = 1.0e-15
        ENDDO
C
C7C-----READ MATRIX DOMAIN YIELD COEFFICIENT INTO ARRAY YIELDMD
        DO K=1,NLAY
          NNDLAY = NODLAY(K)
          NSTRT = NODLAY(K-1)+1
          NDSLAY = NNDLAY - NODLAY(K-1)
          CALL U1DREL(TEMPC(NSTRT),ANAME(9),NDSLAY,K,IN,IOUT)
        ENDDO
C
        DO N = 1, NODES
          YIELDMD(N,ICOMP) = TEMPC(N)
        ENDDO
C
C7D-----READ MATRIX DOMAIN DIFFUSION COEFFICIENT INTO ARRAY DIFFMD
        DO K=1,NLAY
          NNDLAY = NODLAY(K)
          NSTRT = NODLAY(K-1)+1
          NDSLAY = NNDLAY - NODLAY(K-1)
          CALL U1DREL(TEMPC(NSTRT),ANAME(10),NDSLAY,K,IN,IOUT)
        ENDDO
C
        DO N = 1, NODES
          DIFFMD(N,ICOMP) = TEMPC(N)
        ENDDO
C
C7E-----READ AI1MD AND AI2MD IF SHIFT TIME IS NON-ZERO
        IF(TSHIFTMD. GT. 1.0E-10) THEN
C7F-----READ MATRIX DOMAIN COEFFICIENT AIOLD1 INTO ARRAY AIOLD1
          DO K=1,NLAY
            NNDLAY = NODLAY(K)
            NSTRT = NODLAY(K-1)+1
            NDSLAY = NNDLAY - NODLAY(K-1)
            CALL U1DREL(TEMPC(NSTRT),ANAME(11),NDSLAY,K,IN,IOUT)
          ENDDO
C
          DO N = 1, NODES
            AIOLD1(N,ICOMP) = TEMPC(N)
          ENDDO
C
C7G-----READ MATRIX DOMAIN COEFFICIENT AIOLD2 INTO ARRAY AIOLD2
          DO K=1,NLAY
            NNDLAY = NODLAY(K)
            NSTRT = NODLAY(K-1)+1
            NDSLAY = NNDLAY - NODLAY(K-1)
            CALL U1DREL(TEMPC(NSTRT),ANAME(12),NDSLAY,K,IN,IOUT)
          ENDDO
C
          DO N = 1, NODES
            AIOLD2(N,ICOMP) = TEMPC(N)
          ENDDO
        ELSE
C7H ------SET AI1MD AND AI2MD TO ZERO IF SHIFT TIME IS ZERO          
          DO N = 1, NODES
            AIOLD1(N,ICOMP) = 0.0
            AIOLD2(N,ICOMP) = 0.0
          ENDDO
        ENDIF  
      ENDDO
      DEALLOCATE(TEMPC)
C-----------------------------------------------------------------------
C8------PREPARE ARRAYS OF PARAMETERS
C-----------------------------------------------------------------------
C
C9------IF FLOW IS DUAL POROSITY THEN CONDUCTANCE AND STORAGE TERMS ARE ADJUSTED IN DPF. OTHERWISE,
C ------ADJUST CONDUCTANCE AND STORAGE HERE, IF ALSO IFRAHK=1 (PROPERTIES WERE INPUT FOR FRACTURE DOMAIN)
      IF(IDPF.EQ.0)THEN
C9A------SCALE STORAGE COEFFICIENTS BY MOBILE FRACTION IF TRANSIENT FLOW
        ITISTR = 0
        DO NPE=1,NPER
          IF(ISSFLG(NPE).EQ.0) THEN
            ITISTR = 1
            GO TO 123
          ENDIF
        ENDDO
 123    CONTINUE
C
        IF(ITISTR.EQ.0) GO TO 211
C9B-------SCALE SC1 AND SC2 ONLY IF IFRAHK = 1 AND GIVEN FOR FRACTURE WHERE FLOW OCCURS
        IF(IFRAHK.GE.1)THEN
        DO N=1,NODES
          SC1(N) = SC1(N) * VOLFRACMD(N)
        ENDDO
        NCNVRT=0
        DO K = 1,NLAY
          IF(LAYCON(K).NE.0)THEN
            NCNVRT=1
          ENDIF
        ENDDO
        IF(NCNVRT.GT.0)THEN
          DO N=1,NODES
            SC2(N) = SC2(N) * VOLFRACMD(N)
          ENDDO
          IF(IDPF.NE.0)THEN
            DO N=1,NODES
              SC2IM(N) = SC2IM(N) * (1.0 - VOLFRACMD(N))
            ENDDO
C          ELSE
C            ALLOCATE(SC2IM(NODES))
C            DO N=1,NODES
C              SC2IM(N) = 0.0
C            ENDDO
          ENDIF
        ENDIF
        ELSE
C9C ----SPECIFIC YIELD IS FOR TOTAL DOMAIN AND EFFECTIVE POROSITY IS FOR FRACTURE
C9C-----CHECK FOR ERRORS, POROSITY CANNOT BE LESS THAN SPECIFIC YIELD IN FRACTURE
        IEFFPOR = 0
        DO K = 1,NLAY
        NNDLAY = NODLAY(K)
        NSTRT = NODLAY(K-1)+1
        NDSLAY = NNDLAY - NODLAY(K-1)
        DO N = NSTRT,NNDLAY
          IF(LAYCON(K).NE.3.AND.LAYCON(K).NE.2.AND.LAYCON(K).NE.4
     1       .AND.LAYCON(K).NE.5)GO TO 90
          IF(ICBUND(N).EQ.0.OR.IBOUND(N).EQ.0)GO TO 90
          SYIF = SC2(N)/AREA(N) / VOLFRACMD(N)
          IF(PRSITY(N)+1.e-5.LT.SYIF) THEN
              PRSITY(N) = SYIF+1.e-5
              IEFFPOR = 1
c            WRITE(IOUT,56) N, PRSITY(N),SC2(N)/AREA(N),VOLFRACMD(N)
c56          FORMAT(5X,'*** Porosity * phif < Sy FOR NODE: ',I10,' ***'/
c     1      9X,'POROSITY IS ',G10.4,' SC2 IS',G10.4,' PHIF IS',G10.4)
c            STOP
          ENDIF
90        CONTINUE
        ENDDO
        ENDDO
        IF(IEFFPOR.EQ.1) THEN
          WRITE(IOUT,57)
57        FORMAT(5X,'*** EFFECTIVE POROSITY INPUT IS ASSUMED TO BE FOR',
     1    1X,'TOTAL DOMAIN AND IS SCALED FOR FRACTURE DOMAIN ***' /
     2    9X,'POROSITY OF FRACURE HAS TO BE LARGER THAN ITS SPECIFIC',
     2    1X,'YIELD WHICH WAS INPUT FOR TOTAL DOMAIN (IFRACHK = 0)')
        ENDIF
        ENDIF
  211   CONTINUE
C9D-----SCALE CONDUCTIVITY IF PROVIDED FOR FRACTURES ONLY
        IF(IFRAHK.GE.1) THEN
        DO K=1,NLAY
          NNDLAY = NODLAY(K)
          NSTRT = NODLAY(K-1)+1
C9D1---------loop over all nodes within each layer
          DO N=NSTRT,NNDLAY
C9D2------------loop over all connections of node N and fill upper triangle with PGF term
            DO II = IA(N)+1,IA(N+1)-1
              JJ = JA(II)
C9D3--------------only for upper triangle of porous medium nodes
              IF(JJ.LE.N.OR.JJ.GT.NODES) CYCLE
              IIS = JAS(II)
              IF(IFRAHK.EQ.1) THEN
              PGF(IIS) = PGF(IIS) * 0.5*(VOLFRACMD(N)+VOLFRACMD(JJ))
              ELSE
              PGF(IIS) = PGF(IIS) / (0.5*(VOLFRACMD(N)+VOLFRACMD(JJ)))
              ENDIF    
            ENDDO
          ENDDO
          ENDDO         
        DO N = 1,NODES
          IF(IFRAHK.EQ.1)THEN  
            HK(N) = HK(N) * VOLFRACMD(N)
          ELSE
            HK(N) = HK(N) / VOLFRACMD(N)
          ENDIF
        ENDDO
        ENDIF
      ENDIF
C-----------------------------------------------------------------------
C10--------INITIALIZE MOBILE DOMAIN VARIABLES      
C-----------------------------------------------------------------------
C11----FOR ZERO ORDER DECAY IN WATER
      IF(IZOD.EQ.1.OR.IZOD.EQ.3)THEN
        DO ICOMP = 1,MCOMP
        DO N = 1,NODES
          ZODRW(N,ICOMP) = ZODRW(N,ICOMP) * VOLFRACMD(N)
        ENDDO
        ENDDO
      ENDIF
C12----FOR ZERO ORDER DECAY IN SOIL
      IF((IZOD.EQ.2.OR.IZOD.EQ.3).AND.IADSORB.NE.0)THEN
        DO ICOMP = 1,MCOMP
        DO N = 1,NODES
          ZODRS(N,ICOMP) = ZODRS(N,ICOMP) * VOLFRACMD(N)
        ENDDO
        ENDDO
      ENDIF
C13----FOR FIRST ORDER DECAY IN WATER
      IF(IFOD.EQ.1.OR.IFOD.EQ.3)THEN
        DO ICOMP = 1,MCOMP
        DO N = 1,NODES
          FODRW(N,ICOMP) = FODRW(N,ICOMP) * VOLFRACMD(N)
        ENDDO
        ENDDO
      ENDIF
C14----FOR FIRST ORDER DECAY ON SOIL DOES NOT NEED VOLFRACMD SINCE ADSORB TERM HAS IT
C      IF((IFOD.EQ.2.OR.IFOD.EQ.3).AND.IADSORB.NE.0)THEN
C        DO ICOMP = 1,MCOMP
C        DO N = 1,NODES
C          FODRS(N,ICOMP) = FODRS(N,ICOMP) * VOLFRACMD(N)
C        ENDDO
C        ENDDO
C      ENDIF
C15----FOR POROSITY FRACTIONED INTO MOBILE AND IMMOBILE DOMAINS
      DO N = 1,NODES
        PRSITY(N) = PRSITY(N) * VOLFRACMD(N)
      ENDDO
C16----FOR KD VALUES ON SOIL
      IF(IADSORB.NE.0)THEN
        DO ICOMP = 1,MCOMP
        DO N = 1,NODES
          ADSORB(N,ICOMP) = ADSORB(N,ICOMP) * VOLFRACMD(N)
        ENDDO
        ENDDO
      ENDIF
C---------------------------------------------------------------------            
C17--------INITIALIZE MD VARIABLES      
      CALL matrixdsr1 
C18------RETURN
      RETURN
      END
c---------------------------------------------------------------------------------
      subroutine matrixdsr1
c      for MODFLOW-USG-T  falta 6/8/20
c
c      this subroutine is called at the beginning of a simulation before time-stepping
c
c     returns mdtop,mdbot,mdembed,amd1,amd2, and modifies VOLFRACMD if no embedded martix DIFFMDusion
c       also returns some constants, retardm, alpha, and alambdar that are used in other matrix
c       DIFFMDusion routines
c
c     compute matrix DIFFMDusion flags and DIFFMDusion areas from mdflag
c     mdtop is for matrix DIFFMDusion in top aquitard
c     mdbot is for matrix DIFFMDusion in bottom aquitard
c     mdembed is for embedded matrix DIFFMDusion
c     amd1 is area for semi-infinite DIFFMDusion into aquitard(s)
c     amd2 is embedded interfacial area in a gridblock, computed using VOLFRACMD, DIFFLENMD
c        and the gridblock volume
c
c
c     volume fraction is one if there is no embedded md; amd=area for top and bottom
c--------------------------------------------------------------------------------------
c        specifications:
c
      use GWTmdTMODULE
      USE GLOBAL, ONLY: AREA,TOP,BOT,nodes
      USE GWTBCTMODULE, ONLY: MCOMPt,icbund
      IMPLICIT REAL*8 (A-H,O-Z)      
c-------------------------------------------------------------------------------------
c
c      element loop
      do 200 n=1, nodes
c
       if(icbund(n).le.0) cycle 
       mdtop(n)=0
       mdbot(n)=0
       mdembed(n)=0
       amd1(n)=area(n)
       aleng=top(n)-bot(n)
       volu=area(n)*aleng
c     use volume balance to compute interfacial area for embedded
c     this is the area in the gridblock, length**2
       amd2(n)=((1.-VOLFRACMD(n))*volu)/DIFFLENMD(n)
c
       if (mdflag(n).eq.0) then
          VOLFRACMD(n)=1.
       endif
       if (mdflag(n).eq.1) then
          mdbot(n)=1
          VOLFRACMD(n)=1.
       endif
       if (mdflag(n).eq.2) then
          mdembed(n)=1
       endif
       if (mdflag(n).eq.3) then
          mdtop(n)=1
          VOLFRACMD(n)=1.
       endif
       if (mdflag(n).eq.4) then
          mdtop(n)=1
          mdbot(n)=1
          VOLFRACMD(n)=1.
       endif
       if (mdflag(n).eq.5) then
          mdbot(n)=1
          mdembed(n)=1
       endif
       if (mdflag(n).eq.6) then
          mdtop(n)=1
          mdembed(n)=1
       endif
       if (mdflag(n).eq.7) then
          mdtop(n)=1
          mdembed(n)=1
          mdbot(n)=1
       endif
c        component loop
c       calculate some constants used in V-W formulation
       do 100 m=1,mcompt
        retardm(n,m)=1+RHOBMD(n)*KDMD(n,m)/PORMD(n)
        alpha(n,m)=TORTMD(n)*DIFFMD(n,m)/retardm(n,m)
        alambdar(n,m)=DECAYMD(n,m)/retardm(n,m)
 100   continue
c
 200  continue
c
      return
      end
c
c
c
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine matrixdsr2(atime)
c
c         called at the beginning of every time-step before interation
c         atime is the elapsed time during the simulation (at the end of the current time-step)
c         calculate time-dependent matrix DIFFMDusion geometry parameters for 4 species
c         these are computed for both the embedded and semi-infinite geometries
c         falta 6/8/20
c
c         returns d,delt1,delt2,gam1,gam2,beta1,beta2
c------------------------------------------------------------------------------------
      use GWTmdTMODULE
      USE GWTBCTMODULE, ONLY: MCOMPt,icbund
      use GLOBAL, only: nodes      
      IMPLICIT REAL*8 (A-H,O-Z)      
      DOUBLE PRECISION ATIME
c-------------------------------------------------------------------------------------
c
      do 50, n=1,nodes
c  for no matrix DIFFMDusion skip calculations
       if(icbund(n).le.0) cycle           
       If (mdflag(n).eq.0) GoTo 50
       do 70, m=1,mcompt
c
c     DIFFMDusion penetration depth, add new limitation for cases with decay
        d(n,m)=(alpha(n,m)*atime)**0.5/2.
        ssdecay=(TORTMD(n)*DIFFMD(n,m)/DECAYMD(n,m))**.5
        if (d(n,m).gt.ssdecay) then
            d(n,m)=ssdecay
        endif
        d2=d(n,m)*d(n,m)
        d3=d(n,m)*d(n,m)*d(n,m)
c       calculate aquitard/low K zone mass integral parameters delt, gam, beta
c       these depend on the limits of integration specified in mdflag
c      first calculate these for  the infinite aquitard case
        delt1(n,m)=d(n,m)
        gam1(n,m)=d2
        beta1(n,m)=2.*d3
c        additional terms for finite embedded heterogeneity
        expterm=Exp(-DIFFLENMD(n)/d(n,m))
        delt2(n,m)=delt1(n,m)-d(n,m)*expterm
        gam2(n,m)=gam1(n,m)-d(n,m)*DIFFLENMD(n)*expterm-d2*expterm
        beta2(n,m)=beta1(n,m)-DIFFLENMD(n)*DIFFLENMD(n)*d(n,m)*
     &  expterm-2.*d2*DIFFLENMD(n)*
     &  expterm-2.*d3*expterm
c
  70   continue
  50  continue
c
      return
      end
c
c
c
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine matrixdsr3(n,m,dt,cnewt,cnewtm_1,coldt)
c
c      Called in innermost loop at every iteration of every time-step
c      this is inside the element and species loop (n,m)
c      compute key parameters for the semi-analytical matrix DIFFMDusion method
c      used to calculate RHS constant and diagonal coefficient prior to calling
c      linear equation solver
c
c      dt is the current time-step size
c
c      n is the element index, m is the species index (parent=m, daughter =m+1, etc)
c
c      uses the current and previous concentrations, cnewt=cnew(n,m); cnewtm_1=cnew(n,m-1)
c        coldt=cold(n,m)
c
c      this subroutine returns the variables aa1(n,m), aa2(n,m), bb1(n,m) and bb2(n,m)
c         that are used to compute the matrix DIFFMDusion mass flux into the high K zone
c
c      falta 6/8/20
c
c----------------------------------------------------------------------------
      use GWTmdTMODULE
      IMPLICIT REAL*8 (A-H,O-Z)      
      DOUBLE PRECISION DT,CNEWT,CNEWTM_1,COLDT      
c-------------------------------------------------------------------------------------
c
c**********    parameters for semi-infinite geometry  ****************
c
      a1=beta1(n,m)+alambdar(n,m)*dt*beta1(n,m)
      b1=gam1(n,m)+alpha(n,m)*dt+alambdar(n,m)*dt*gam1(n,m)
      e1=delt1(n,m)-alpha(n,m)*dt/d(n,m)+
     &alambdar(n,m)*dt*delt1(n,m)
c      compute parameters aa and bb used for V-W flux
c***************************************************************************
       aatop=-e1-a1/(2.*alpha(n,m)*dt)+a1/
     &(2.*d(n,m)*d(n,m))-a1*alambdar(n,m)/(2.*alpha(n,m))
       aabot=a1/d(n,m)+b1
       aa1(n,m)=aatop/aabot
       bbtop=aiold1(n,m)+a1*coldt/(2.*alpha(n,m)*dt)
c term for daughter production in matrix
       If (m.gt.1) Then
         bbtop=bbtop+YIELDMD(n,m)*DECAYMD(n,m-1)*
     &dt*ainew1(n,m-1)/retardm(n,m)
         bbtop=bbtop+a1*YIELDMD(n,m)*DECAYMD(n,m-1)*
     &cnewtm_1/(2.*alpha(n,m)*retardm(n,m))
       End If
c
       bbbot=a1/d(n,m)+b1
       bb1(n,m)=bbtop/bbbot
c       print *, 'aa(1)=',aa(1),'bb(1,1,1,1)=',bb(1,1,1,1)
c
c*************** parameters for finite case for embedded heterogeneity **********
c
      a2=beta2(n,m)+alambdar(n,m)*dt*beta2(n,m)
      b2=gam2(n,m)+alpha(n,m)*dt+alambdar(n,m)*dt*gam2(n,m)
      e2=delt2(n,m)-alpha(n,m)*dt/d(n,m)
     &+alambdar(n,m)*dt*delt2(n,m)
c       compute parameters aa and bb needed for V-W flux
c****************************************************************************
       aatop=-e2-a2/(2.*alpha(n,m)*dt)+a2/
     &(2.*d(n,m)*d(n,m))-a2*alambdar(n,m)/(2.*alpha(n,m))
       aabot=a2/d(n,m)+b2
       aa2(n,m)=aatop/aabot
       bbtop=aiold2(n,m)+a2*coldt/(2.*alpha(n,m)*dt)
c term for daughter production in matrix
       If (m.gt.1) Then
         bbtop=bbtop+YIELDMD(n,m)*DECAYMD(n,m-1)*
     &dt*ainew2(n,m-1)/retardm(n,m)
         bbtop=bbtop+a2*YIELDMD(n,m)*DECAYMD(n,m-1)*
     &cnewtm_1/(2.*alpha(n,m)*retardm(n,m))
       End If
c
       bbbot=a2/d(n,m)+b2
       bb2(n,m)=bbtop/bbbot
c
      return
      end
cc
c
c
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine matrixdsr4rhs(n,m,RTERMD)
c
c       compute contribution to RHS of matrix equations from matrix DIFFMDusion mass flux
c       falta 6/8/20
c       n is element index, m is component index
c       returns RTERMD which is added to RHS in USG-Transport.  This is both element and
c         component specific
c----------------------------------------------------------------------------
      use GWTmdTMODULE
      IMPLICIT REAL*8 (A-H,O-Z)      
      DOUBLE PRECISION RTERMD
c-------------------------------------------------------------------------------------
c
      anum2=0.
c         embedded heterogenity, mdembed=1
      If (mdembed(n).eq.1) Then
          anum2=amd2(n)*PORMD(n)*TORTMD(n)*DIFFMD(n,m)*bb2(n,m)
      End If
c        semi-infinite in bottom layer, mdbot=1
      If(mdbot(n).eq.1)  Then
          anum2=anum2+amd1(n)*PORMD(n)*TORTMD(n)*DIFFMD(n,m)*bb1(n,m)
      End If
c        semi-infinite in top layer, mdtop=1
      If(mdtop(n).eq.1)  Then
          anum2=anum2+amd1(n)*PORMD(n)*TORTMD(n)*DIFFMD(n,m)*bb1(n,m)
      End If
      RTERMD=anum2
c
      return
      end
c
c
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine matrixdsr5diag(n,m,DTERMD)
c       This is the contribution to the diagonal coefficient in the matrix
c       equation.  In MODFLOW-USG-Transport this is called DTERMD, and it is added
c       to AMAT(IPIV).  This will need to be both element and component specific
c       n is the element index and m is the component indes
c       subroutine returns DTERMD  falta 6/8/20
c
c----------------------------------------------------------------------------
      use GWTmdTMODULE
      IMPLICIT REAL*8 (A-H,O-Z)      
      DOUBLE PRECISION DTERMD
c-------------------------------------------------------------------------------------
c
      denom2 = 0.
c    embedded heterogenity, mdembed=1
      If (mdembed(n).eq.1) Then
          denom2=amd2(n)*PORMD(n)*TORTMD(n)*DIFFMD(n,m)/d(n,m)-
     &    amd2(n)*PORMD(n)*TORTMD(n)*DIFFMD(n,m)*aa2(n,m)
      End If
c    semi-infinite on bottom mdbot=1
      If (mdbot(n).eq.1) Then
          denom2=denom2+amd1(n)*PORMD(n)*TORTMD(n)*DIFFMD(n,m)/d(n,m)-
     &amd1(n)*PORMD(n)*TORTMD(n)*DIFFMD(n,m)*aa1(n,m)
      End If
c    semi-infinite on top mdtop=1
      If (mdtop(n).eq.1) Then
          denom2=denom2+amd1(n)*PORMD(n)*TORTMD(n)*DIFFMD(n,m)/d(n,m)-
     &amd1(n)*PORMD(n)*TORTMD(n)*DIFFMD(n,m)*aa1(n,m)
      End If
      DTERMD=denom2
c
      return
      end
c
c
c
c
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c
      subroutine matrixdsr6int(dt)
c
c        This subroutine is called after the linear equations have been solved at the current
c        iteration, but before convergence of the iterative solution.
c        It uses the new concentrations, CONC(n,m) and CONC(n,m-1) and the old concentration CONCO(n,m)
c
c        it returns the updated concentration integrals ainew1(n,m) and ainew2(n,m)
c
c        it also returns the fitting parameters p1(n,m), p2(n,m), q1(n,m) and q2(n,m)
c           these are not needed outside of this subroutine, but can be output if desired to compute
c           concentration profiles in the low K material at the current time.
c
c         compute V-W p and q to compute mass integral at ainew at current timestep
c        falta 6/8/20
c
c----------------------------------------------------------------------------
      use GWTmdTMODULE
      USE GWTBCTMODULE, ONLY: MCOMPt,CONC,CONCO,icbund
      use GLOBAL, only: nodes
      IMPLICIT REAL*8 (A-H,O-Z)      
      DOUBLE PRECISION DT
c-------------------------------------------------------------------------------------
c
      do 100 n=1, nodes
       if(icbund(n).le.0) cycle           
       do 200 m=1,mcompt
c   embedded heterogenity, mdembed=1
        If (mdembed(n).eq.1) Then
          p2(n,m)=aa2(n,m)*CONC(n,m)+bb2(n,m)
          qtop1=(CONC(n,m)-CONCO(n,m))*d(n,m)*d(n,m)/(alpha(n,m)*dt)
     &    -CONC(n,m)
          qtop2=2.*d(n,m)*p2(n,m)
     &    +alambdar(n,m)*CONC(n,m)*d(n,m)*d(n,m)/alpha(n,m)
c        term for daughter production in matrix
          If (m.gt.1) Then
            qtop2=qtop2-YIELDMD(n,m)*DECAYMD(n,m-1)*
     &      CONC(n,m-1)*d(n,m)*d(n,m)/(alpha(n,m)*retardm(n,m))
          End If
          q2(n,m)=(qtop1+qtop2)/(2.*d(n,m)*d(n,m))
          ainew2(n,m)=delt2(n,m)*CONC(n,m)
     &    +gam2(n,m)*p2(n,m)+beta2(n,m)*q2(n,m)
        End If
c
c    semi-infinite in bottom, mdbot=1
        If (mdbot(n).eq.1) Then
          p1(n,m)=aa1(n,m)*CONC(n,m)+bb1(n,m)
          qtop1=(CONC(n,m)-CONCO(n,m))*d(n,m)*d(n,m)/
     &    (alpha(n,m)*dt)-CONC(n,m)
          qtop2=2.*d(n,m)*p1(n,m)
     &    +alambdar(n,m)*CONC(n,m)*d(n,m)*d(n,m)/alpha(n,m)
c            term for daughter production in matrix
           If (m.gt.1) Then
            qtop2=qtop2-YIELDMD(n,m)*DECAYMD(n,m-1)
     &      *CONC(n,m-1)*d(n,m)*d(n,m)/(alpha(n,m)*retardm(n,m))
          End If
          q1(n,m)=(qtop1+qtop2)/(2.*d(n,m)*d(n,m))
          ainew1(n,m)=delt1(n,m)*CONC(n,m)+gam1(n,m)*p1(n,m)
     &    +beta1(n,m)*q1(n,m)
        End If
c
c    semi-infinite in top, mdtop=1
        If (mdtop(n).eq.1) Then
          p1(n,m)=aa1(n,m)*CONC(n,m)+bb1(n,m)
          qtop1=(CONC(n,m)-CONCO(n,m))*d(n,m)*d(n,m)/
     &    (alpha(n,m)*dt)-CONC(n,m)
          qtop2=2.*d(n,m)*p1(n,m)
     &    +alambdar(n,m)*CONC(n,m)*d(n,m)*d(n,m)/alpha(n,m)
c            term for daughter production in matrix
           If (m.gt.1) Then
            qtop2=qtop2-YIELDMD(n,m)*DECAYMD(n,m-1)
     &      *CONC(n,m-1)*d(n,m)*d(n,m)/(alpha(n,m)*retardm(n,m))
          End If
          q1(n,m)=(qtop1+qtop2)/(2.*d(n,m)*d(n,m))
          ainew1(n,m)=delt1(n,m)*CONC(n,m)+gam1(n,m)*p1(n,m)
     &    +beta1(n,m)*q1(n,m)
        End If
c
 200   continue
 100  continue
c
      return
      end
c
c
c
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c
      subroutine matrixdsr7updt
c
      use GWTmdTMODULE
      USE GWTBCTMODULE, ONLY: MCOMPt,icbund
      use GLOBAL, only: nodes      
      IMPLICIT REAL*8 (A-H,O-Z)      
c
c               update aiold1 and aiold2 upon convergence at each time-step
c               falta 6/10/20
c
c                  update VW concentration integrals
      do 100 n=1,nodes
       if(icbund(n).le.0) cycle           
        do 200 m=1,mcompt
          If (mdembed (n).eq.1) Then
           aiold2(n,m)=ainew2(n,m)
          End If
          If (mdbot(n).eq.1) Then
           aiold1(n,m)=ainew1(n,m)
          End If
          If (mdtop(n).eq.1) Then
           aiold1(n,m)=ainew1(n,m)
          End If
c
  200   continue
  100 continue
c
c--------this is where the concentration integrals ainew1 ainew2 could be written to output at the end of a simulation
c
      return
      end
c
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c
      subroutine matrixdsr8 (kstp,kper,icomp)
      use GWTmdTMODULE
      USE GLOBAL,   ONLY:IOUT,NCOL,NROW,NLAY,IBOUND,BUFF,NODES,IUNSTR,
     1  ISSFLG
      USE GWFBASMODULE,ONLY:MSUM,ISPCFL,IAUXSV,DELT,PERTIM,TOTIM
      USE GWTBCTMODULE, ONLY: MSUMT,VBVLT,VBNMT,ICT,MCOMPt,CONC,CONCO,
     1  ICBUND
      IMPLICIT REAL*8 (A-H,O-Z)            
C
      CHARACTER*16 TEXT
      DOUBLE PRECISION RATIN,RATOUT,QQ,VODT,ADSTERM,FL,CW,CWO,ALENG,
     *  DTERMS,RTERMS,VOLU
      DATA TEXT /' MD MASS STORAGE'/      
C     ------------------------------------------------------------------
C
C1------CLEAR RATIN AND RATOUT ACCUMULATORS, AND SET CELL-BY-CELL
C1------BUDGET FLAG.
      ZERO=0.
      RATIN=ZERO
      RATOUT=ZERO
      QQ = ZERO
      IBD=0
      ISS = ISSFLG(KPER)
      IF(ImdTCB.LT.0 .AND. ISPCFL.NE.0) IBD=-1
      IF(ImdTCB.GT.0) IBD=ISPCFL
      IBDLBL=0
C
C3------CLEAR THE BUFFER.
      DO 50 N=1,NODES
      BUFF(N)=ZERO
50    CONTINUE
C
C5------LOOP THROUGH EACH NODE AND CALCULATE STORAGE      
c
      DO N=1,NODES
        IF(ICBUND(N).LE.0) CYCLE      
c
c==========================================================================================
c-------------compute mass flow of component icomp into or out of matrix domain falta 7/26/20
c-------------If mass flow (qq) is positive, mass is flowing out of matrix into the high K zone
c-------------if mass flow is negative, mass is flowing into matrix out of high K zone
c
        ICOMPM1 = ICOMP-1
        IF(ICOMP.EQ.1) ICOMPM1 = ICOMP 
        CALL matrixdsr3 (N,ICOMP,DELT,CONC(N,ICOMP),
     1    CONC(N,ICOMPM1),CONCO(N,ICOMP))
        CALL matrixdsr4rhs(n,ICOMP,RTERMD)
        CALL matrixdsr5diag(n,ICOMP,DTERMD)
        qq = -DTERMD *conc(n,icomp) + RTERMD          
c=============================================================================  
c
c       if(mdembed(n).eq.1) then
c        qq=amd2(n)*pormd(n)*tortmd(n)*diffmd(n,icomp)*((aa2(n,icomp)
c     &  -1./d(n,icomp))*conc(n,icomp)+bb2(n,icomp))
c        endif
c       if(mdbot(n).eq.1) then
c        qq=qq+amd1(n)*pormd(n)*tortmd(n)*diffmd(n,icomp)*((aa1(n,icomp)
c     1   -1./d(n,icomp))*conc(n,icomp)+bb1(n,icomp))
c       endif
c       if(mdtop(n).eq.1) then
c        qq=qq+amd1(n)*pormd(n)*tortmd(n)*diffmd(n,icomp)*((aa1(n,icomp)
c     1   -1./d(n,icomp))*conc(n,icomp)+bb1(n,icomp))
c       endif
c  
c=============================================================================        
C-------------compute mass storage rate for md in matrix for various mdflag options
c        if(mdflag(n).eq.2) then ! embedded
c          qq = aiold2(n,icomp) * amd2(n) * pormd(n) * retardm(n,icomp) ! (New - Old / delt)
c        elseif(mdflag(n).eq.3) then !whatever
c          qq =  aiold1(n,icomp) * amd1(n) * pormd(n) * retardm(n,icomp)
c        elseif(mdflag(n).eq.4) then !whatever
c
c        endif
c=============================================================================
        Q = QQ
C
C5D-----PRINT FLOW RATE IF REQUESTED.
        IF(IBD.LT.0) THEN
           IF(IBDLBL.EQ.0) WRITE(IOUT,61) TEXT,KPER,KSTP
   61      FORMAT(1X,/1X,A,'   PERIOD ',I4,'   STEP ',I3)
          IF(IUNSTR.EQ.0.AND.N.LE.NODES)THEN
            IL = (N-1) / (NCOL*NROW) + 1
            IJ = N - (IL-1)*NCOL*NROW
            IR = (IJ-1)/NCOL + 1
            IC = IJ - (IR-1)*NCOL
            WRITE(IOUT,62) L,IL,IR,IC,Q
   62       FORMAT(1X,'CBC  ',I6,'   LAYER ',I3,'   ROW ',I5,'   COL ',
     1         I5, '   FLUX ',1PG15.6)
          ELSE
            WRITE(IOUT,63) L,N,Q
   63      FORMAT(1X,'CBC  ',I6,'    NODE ',I8,'   FLUX ',1PG15.6)
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
C5G-----POSITIVE FLOW RATE. ADD IT TO RATIN
          RATIN=RATIN+QQ
        ELSE
C
C5H-----NEGATIVE FLOW RATE. ADD IT TO RATOUT
          RATOUT=RATOUT-QQ
        END IF
C
      ENDDO
C
C6------IF CELL-BY-CELL FLOWS WILL BE SAVED AS A 3-D ARRAY,
C6------CALL UBUDSV TO SAVE THEM.
      IF(IBD.EQ.1)THEN
        IF(IUNSTR.EQ.0)THEN
          CALL UBUDSV(KSTP,KPER,TEXT,ImdTCB,BUFF,NCOL,NROW,
     1                          NLAY,IOUT)
        ELSE
          CALL UBUDSVU(KSTP,KPER,TEXT,ImdTCB,BUFF,NODES,
     1                          IOUT,PERTIM,TOTIM)
        ENDIF
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
      return
      end    
C---------------------------------------------------------------------------------
      SUBROUTINE GWT2mdt1OT(KSTP,KPER)
C     ******************************************************************
C     OUTPUT coefficients a1md and a2md FOR EACH SPECIES 
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,      ONLY:ITMUNI,IOUT,IUNSTR,ISSFLG,PERLEN,INCLN
      USE GWFBASMODULE,ONLY:DELT,PERTIM,TOTIM,ISPCFL,IBUDFL,IATS,TMINAT,
     1  NPTIMES,NPSTPS,ITIMOT,TIMOT,ISPCFLAT
      USE GWTBCTMODULE, ONLY: MCOMP,VBVLT,VBNMT,MSUMT,MCOMPT,NTCOMP
      USE GWTmdtMODULE, ONLY: IMDTCF      
C     ------------------------------------------------------------------
C1A---RETURN IF CONCENTRATION PRINT FLAG IS OFF OR IF MDT COEFFICIENTS ARE NOT WRITTEN 
      IF(ISPCFL.EQ.0) RETURN   
      IF(IMDTCF.EQ.0) RETURN      
C
C1------CLEAR PRINTOUT FLAG (IPFLG)      
      IPFLG=0
      WRITE(IOUT,1)
1     FORMAT(/1X,'OUTPUT OF MATRIX DIFFUSION COEFFICIENTS'/1X,40('-')) 
      DO ICOMP = 1,NTCOMP 
C2------WRITE MD COEFFICIENTS IN ACCORDANCE WITH FLAGS IN IOFLG.
        WRITE(IOUT,2)ICOMP
2       FORMAT(/5X,'MATRIX DIFFUSION OUTPUT FOR COMPONENT SPECIES',
     *    1X,'NUMBER',I5/5X,60('-'))
C
        IF(IUNSTR.EQ.0)THEN ! WRITE M2K5 STYLE FOR STRUCTURED GRID
          CALL SGWT2mdt1C(KSTP,KPER,IPFLG,ICOMP)
        ELSE
          CALL SGWT2mdt1CU(KSTP,KPER,IPFLG,ICOMP)
        ENDIF
        IPFLG = 1
      ENDDO        
C-------------------------------------------------------------
C
        CALL SGWF2BAS7T(KSTP,KPER,DELT,PERTIM,TOTIM,ITMUNI,IOUT)
        WRITE(IOUT,101)
  101   FORMAT('1') 
C
C5------RETURN
      RETURN
      END
C-----------------------------------------------------------------------  
C----------------------------------------------------------------------
      SUBROUTINE SGWT2mdt1C(KSTP,KPER,IPFLG,ICOMP)
C     ******************************************************************
C     PRINT AND RECORD AI1 AND AI2 COEFFICIENTS FOR STRUCTURED GWF GRID
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,      ONLY:NCOL,NROW,NLAY,IXSEC,HNEW,NODLAY,
     1                      IBOUND,IOUT
      USE GWFBASMODULE,ONLY:PERTIM,TOTIM,IOFLG,ISPCFM,LBHDSV,
     2                      CSPCFM,IOFLG
      USE GWTmdtMODULE, ONLY: AINEW1, AINEW2, MDFLAG, IMDTCF,IMDTCF2,
     1                      IMULTIMD,CROOTNAME_MD 
C
      REAL,   SAVE,    DIMENSION(:,:,:),    ALLOCATABLE ::BUFF1, BUFF2
      CHARACTER*16 TEXT1, TEXT2,FORMC,ISTATUS
      CHARACTER*400 CFILENAME
      INCLUDE 'openspec.inc'   
C     ------------------------------------------------------------------
      ALLOCATE(BUFF1(NCOL,NROW,NLAY),BUFF2(NCOL,NROW,NLAY))
      CALL GET_TEXTAI1(ICOMP,TEXT1)
      CALL GET_TEXTAI2(ICOMP,TEXT2)
C
C1------FOR EACH LAYER MOVE AI1 TO BUFF IF PRINT OR SAVE IS REQUESTED.
      DO 59 K=1,NLAY
C
C2------IS CONC NEEDED FOR THIS LAYER?
      KL=K
      IF(IXSEC.NE.0) KL=1
      IF(IOFLG(KL,5).EQ.0 .AND. IOFLG(KL,6).EQ.0) GO TO 59
C
C3------MOVE CONCIM TO BUFF FOR THE LAYER.
      DO 58 I=1,NROW
      DO 58 J=1,NCOL
      N = (K-1)*NROW*NCOL + (I-1)*NCOL + J
      BUFF1(J,I,K)=AINEW1(N,ICOMP)
      BUFF2(J,I,K)=AINEW2(N,ICOMP)
   58 CONTINUE
   59 CONTINUE
C
C4------FOR EACH LAYER: DETERMINE IF AI1 AND AI2 SHOULD BE PRINTED.
C4------IF SO THEN CALL ULAPRS OR ULAPRW TO PRINT AI1 AND AI2.
      IF(IMDTCF.LT.0) THEN
         IF(IXSEC.EQ.0) THEN
C4B------PRING AI1               
           DO 69 K=1,NLAY
             KK=K
             IF(IOFLG(K,5).EQ.0) GO TO 69
             IF(ISPCFM.LT.0) CALL ULAPRS(BUFF1(1,1,K),TEXT1,KSTP,KPER,
     1               NCOL,NROW,KK,-ISPCFM,IOUT)
             IF(ISPCFM.GE.0) CALL ULAPRW(BUFF1(1,1,K),TEXT1,KSTP,KPER,
     1               NCOL,NROW,KK,ISPCFM,IOUT)
   69      CONTINUE
C4B------PRING AI2           
          DO 68 K=1,NLAY
           KK=K
           IF(IOFLG(K,5).EQ.0) GO TO 68
           IF(ISPCFM.LT.0) CALL ULAPRS(BUFF2(1,1,K),TEXT2,KSTP,KPER,
     1               NCOL,NROW,KK,-ISPCFM,IOUT)
           IF(ISPCFM.GE.0) CALL ULAPRW(BUFF2(1,1,K),TEXT2,KSTP,KPER,
     1               NCOL,NROW,KK,ISPCFM,IOUT)
           IPFLG=1
   68      CONTINUE
C
C4C-----PRINT AI1 AND AI2 FOR CROSS SECTION.
         ELSE
           IF(IOFLG(1,5).NE.0) THEN
C4D------PRINT AI1               
             IF(ISPCFM.LT.0) CALL ULAPRS(BUFF1,TEXT1,KSTP,KPER,
     1                 NCOL,NLAY,-1,-ISPCFM,IOUT)
             IF(ISPCFM.GE.0) CALL ULAPRW(BUFF1,TEXT1,KSTP,KPER,
     1                 NCOL,NLAY,-1,ISPCFM,IOUT)
C4E------PRING AI2             
             IF(ISPCFM.LT.0) CALL ULAPRS(BUFF2,TEXT2,KSTP,KPER,
     1                 NCOL,NLAY,-1,-ISPCFM,IOUT)
             IF(ISPCFM.GE.0) CALL ULAPRW(BUFF2,TEXT2,KSTP,KPER,
     1                 NCOL,NLAY,-1,ISPCFM,IOUT)
             IPFLG=1             
           END IF
         END IF
      END IF
C
C5------FOR EACH LAYER: DETERMINE IF AI1 AND AI2 SHOULD BE SAVED ON DISK.
C5------IF SO THEN CALL ULASAV OR ULASV2 TO SAVE AI1 AND AI2.
      IFIRST=1
      IF(IMDTCF.LE.0) GO TO 80
C--------------------------------------------------------------------------------------------------------------------
C5A ----FOR MULTIFILE OUTPUT OPTION, ADD COMPONENT SPECIES TO ROOT NAME AND OPEN FILE ON ISPCFM AND GO TO END OF FILE
      IF(IMULTIMD.NE.0)THEN
        CFILENAME = trim(CROOTNAME_MD) // TRIM(TEXT1) //'.AI1'
        ISTATUS = 'OLD'
C
        OPEN(ImdtCF, FILE = CFILENAME, FORM = FORM,
     1  STATUS = ISTATUS, POSITION = 'APPEND')
C -------IF AI2 IS WRITTEN TO A SEPARATE FILE THEN OPEN THAT FILE ALSO        
        IF(IMDTCF2. NE. 0) THEN 
          CFILENAME = trim(CROOTNAME_MD) // TRIM(TEXT2) //'.AI2'  
          OPEN(ImdtCF2, FILE = CFILENAME, FORM = FORM,
     1    STATUS = ISTATUS, POSITION = 'APPEND')          
        ENDIF    
      ENDIF
C--------------------------------------------------------------------------------------------------------------------        
      IF(IXSEC.EQ.0) THEN
C5A------SAVE AI1           
        DO 79 K=1,NLAY
        NSTRT = NODLAY(K-1)+1
        KK=K
        IF(IOFLG(K,6).EQ.0) GO TO 79
        IF(IFIRST.EQ.1) WRITE(IOUT,74) IMDTCF,KSTP,KPER
   74   FORMAT(1X,/1X,'AI1 COEFFICIENTS WILL BE SAVED ON UNIT ',
     1      I4,' AT END OF TIME STEP ',I3,', STRESS PERIOD ',I4)
        IFIRST=0       
        IF(CSPCFM.EQ.' ') THEN
           CALL ULASAV(BUFF1(1,1,K),TEXT1,KSTP,KPER,PERTIM,TOTIM,NCOL,
     1                NROW,KK,IMDTCF)
        ELSE
           CALL ULASV2(BUFF1(1,1,K),TEXT1,KSTP,KPER,PERTIM,TOTIM,NCOL,
     1                NROW,KK,IMDTCF,CSPCFM,LBHDSV,MDFLAG(NSTRT))
        END IF        
   79 CONTINUE
C5B------SAVE AI2      
        IFIRST = 1
        DO 78 K=1,NLAY
        NSTRT = NODLAY(K-1)+1
        KK=K
        IF(IOFLG(K,6).EQ.0) GO TO 78
        IF(IFIRST.EQ.1) WRITE(IOUT,73) IMDTCF2,KSTP,KPER
   73   FORMAT(1X,/1X,'AI2 COEFFICIENTS WILL BE SAVED ON UNIT ',
     1      I4,' AT END OF TIME STEP ',I3,', STRESS PERIOD ',I4)
        IFIRST=0
C5B------SAVE AI2
        IF(CSPCFM.EQ.' ') THEN
           CALL ULASAV(BUFF2(1,1,K),TEXT2,KSTP,KPER,PERTIM,TOTIM,NCOL,
     1                NROW,KK,IMDTCF2)
        ELSE
           CALL ULASV2(BUFF2(1,1,K),TEXT2,KSTP,KPER,PERTIM,TOTIM,NCOL,
     1                NROW,KK,IMDTCF2,CSPCFM,LBHDSV,MDFLAG(NSTRT))
        END IF        
   78   CONTINUE      
C
C5C-----SAVE CONC FOR CROSS SECTION.
      ELSE
        IF(IOFLG(1,6).NE.0) THEN
C5D------SAVE AI1           
          WRITE(IOUT,74) IMDTCF,KSTP,KPER            
          IF(CSPCFM.EQ.' ') THEN
             CALL ULASAV(BUFF1,TEXT1,KSTP,KPER,PERTIM,TOTIM,NCOL,
     1                NLAY,-1,IMDTCF)
          ELSE
             CALL ULASV2(BUFF1,TEXT1,KSTP,KPER,PERTIM,TOTIM,NCOL,
     1                  NLAY,-1,IMDTCF,CSPCFM,LBHDSV,MDFLAG)
          END IF
C5E------SAVE AI2 
          WRITE(IOUT,73) IMDTCF,KSTP,KPER
          IF(CSPCFM.EQ.' ') THEN
             CALL ULASAV(BUFF2,TEXT2,KSTP,KPER,PERTIM,TOTIM,NCOL,
     1                NLAY,-1,IMDTCF2)
          ELSE
             CALL ULASV2(BUFF2,TEXT2,KSTP,KPER,PERTIM,TOTIM,NCOL,
     1                  NLAY,-1,IMDTCF2,CSPCFM,LBHDSV,MDFLAG)
          END IF          
        END IF
      END IF
80    CONTINUE
      DEALLOCATE(BUFF1, BUFF2)
C
C6------RETURN.
      RETURN
      END
C----------------------------------------------------------------------
      SUBROUTINE SGWT2mdt1CU(KSTP,KPER,IPFLG,ICOMP)
C     ******************************************************************
C     PRINT AND RECORD AI1 AND AI2 COEFFICIENTS FOR UNSTRUCTURED GWF GRID
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,      ONLY:NCOL,NROW,NLAY,IXSEC,HNEW,NODLAY,
     1                      IBOUND,IOUT,NODES,BUFF
      USE GWFBASMODULE,ONLY:PERTIM,TOTIM,IOFLG,ISPCFM,LBHDSV,
     2                      CSPCFM,IOFLG
      USE GWTmdtMODULE, ONLY: AINEW1, AINEW2, MDFLAG, IMDTCF,IMDTCF2,
     1                        IMULTIMD,CROOTNAME_MD
      REAL,   SAVE,    DIMENSION(:),    ALLOCATABLE ::BUFF2
C
      CHARACTER*16 TEXT1, TEXT2,ISTATUS 
      CHARACTER*400 CFILENAME
      INCLUDE 'openspec.inc'   
C     ------------------------------------------------------------------
C
      ALLOCATE (BUFF2(NODES))
      CALL GET_TEXTAI1(ICOMP,TEXT1)
      CALL GET_TEXTAI2(ICOMP,TEXT2)
C1------FOR EACH LAYER MOVE CONC TO BUFF IF PRINT OR SAVE IS REQUESTED.
      DO 59 K=1,NLAY
C
C2------IS CONC NEEDED FOR THIS LAYER?
      KL=K
      IF(IXSEC.NE.0) KL=1
      IF(IOFLG(KL,5).EQ.0 .AND. IOFLG(KL,6).EQ.0) GO TO 59
C
C3------MOVE CONC TO BUFF FOR THE LAYER.
      NNDLAY = NODLAY(K)
      NSTRT = NODLAY(K-1)+1
      DO 58 N=NSTRT,NNDLAY
      BUFF(N)=AINEW1(N,ICOMP)
      BUFF2(N)=AINEW2(N,ICOMP)
   58 CONTINUE
   59 CONTINUE
C
C4------FOR EACH LAYER: DETERMINE IF COEFFICIENTS AI1 AND AI2 SHOULD BE PRINTED.
C4------IF SO THEN CALL ULAPRU TO PRINT CONC.
      IF(IMDTCF.LT.0) THEN
         IF(IXSEC.EQ.0) THEN
C4A-----PRINT AINEW1             
           DO 69 K=1,NLAY
           KK=K
           IF(IOFLG(K,5).EQ.0) GO TO 69
           NNDLAY = NODLAY(K)
           NSTRT = NODLAY(K-1)+1
           CALL ULAPRU(BUFF,TEXT1,KSTP,KPER,
     1           NSTRT,NNDLAY,KK,IABS(ISPCFM),IOUT,PERTIM,TOTIM,NODES)
           IPFLG=1
   69      CONTINUE
C4B-----PRINT AINEW2            
           DO 68 K=1,NLAY
           KK=K
           IF(IOFLG(K,5).EQ.0) GO TO 68
           NNDLAY = NODLAY(K)
           NSTRT = NODLAY(K-1)+1
           CALL ULAPRU(BUFF2,TEXT2,KSTP,KPER,
     1           NSTRT,NNDLAY,KK,IABS(ISPCFM),IOUT,PERTIM,TOTIM,NODES)
           IPFLG=1
   68      CONTINUE           
C
C4C-----PRINT AI1 AND AI2 FOR CROSS SECTION.
         ELSE
           IF(IOFLG(1,5).NE.0) THEN
           IF(ISPCFM.NE.0) CALL ULAPRU(BUFF,TEXT1,KSTP,KPER,
     1           NSTRT,NNDLAY,-1,IABS(ISPCFM),IOUT,PERTIM,TOTIM,NODES)
           IF(ISPCFM.NE.0) CALL ULAPRU(BUFF2,TEXT2,KSTP,KPER,
     1           NSTRT,NNDLAY,-1,IABS(ISPCFM),IOUT,PERTIM,TOTIM,NODES)
             IPFLG=1
C
           END IF
         END IF
      END IF
C
C5------FOR EACH LAYER: DETERMINE IF CONC SHOULD BE SAVED ON DISK.
C5------IF SO THEN CALL ULASAV OR ULASV2 TO SAVE AI1 AND AI2 COEFFICIENTS.
      IFIRST=1
      IF(ImdtCF.LE.0) GO TO 80
C--------------------------------------------------------------------------------------------------------------------
C5A ----FOR MULTIFILE OUTPUT OPTION, ADD COMPONENT SPECIES TO ROOT NAME AND OPEN FILE ON ISPCFM AND GO TO END OF FILE
      IF(IMULTIMD.NE.0)THEN
        CFILENAME = trim(CROOTNAME_MD) // TRIM(TEXT1) //'.AI1'
        ISTATUS = 'OLD'
C
        OPEN(ImdtCF, FILE = CFILENAME, FORM = FORM,
     1  STATUS = ISTATUS, POSITION = 'APPEND')
C -------IF AI2 IS WRITTEN TO A SEPARATE FILE THEN OPEN THAT FILE ALSO        
        IF(IMDTCF2. NE. 0) THEN 
          CFILENAME = trim(CROOTNAME_MD) // TRIM(TEXT2) //'.AI2'  
          OPEN(ImdtCF2, FILE = CFILENAME, FORM = FORM,
     1    STATUS = ISTATUS, POSITION = 'APPEND')          
        ENDIF    
      ENDIF
C--------------------------------------------------------------------------------------------------------------------              
      IF(IXSEC.EQ.0) THEN
C5A-----SAVE AI1          
        DO 79 K=1,NLAY
        KK=K
        IF(IOFLG(K,6).EQ.0) GO TO 79
        NNDLAY = NODLAY(K)
        NSTRT = NODLAY(K-1)+1
        IF(IFIRST.EQ.1) WRITE(IOUT,74) ImdtCF,KSTP,KPER
   74   FORMAT(1X,/1X,'AI1 COEFFICIENT WILL BE SAVED ON UNIT ',I8,
     1      ' AT END OF TIME STEP ',I8,', STRESS PERIOD ',I8)
        IFIRST=0
        IF(CSPCFM.EQ.' ') THEN
           CALL ULASAVU(BUFF,TEXT1,KSTP,KPER,PERTIM,TOTIM,NSTRT,
     1                NNDLAY,KK,ImdtCF,NODES)
        ELSE
           CALL ULASV2U(BUFF,TEXT1,KSTP,KPER,PERTIM,TOTIM,NSTRT,
     1            NNDLAY,KK,ImdtCF,CSPCFM,LBHDSV,MDFLAG(NSTRT),NODES)
        END IF
        IPFLG=1
   79   CONTINUE
C5B-----SAVE AI2
        IFIRST=1
        DO 78 K=1,NLAY
        KK=K
        IF(IOFLG(K,6).EQ.0) GO TO 78
        NNDLAY = NODLAY(K)
        NSTRT = NODLAY(K-1)+1
        IF(IFIRST.EQ.1) WRITE(IOUT,73) ImdtCF2,KSTP,KPER
   73   FORMAT(1X,/1X,'AI2 COEFFICIENT WILL BE SAVED ON UNIT ',I8,
     1      ' AT END OF TIME STEP ',I8,', STRESS PERIOD ',I8)
        IFIRST=0
        IF(CSPCFM.EQ.' ') THEN
           CALL ULASAVU(BUFF2,TEXT2,KSTP,KPER,PERTIM,TOTIM,NSTRT,
     1                NNDLAY,KK,ImdtCF2,NODES)
        ELSE
           CALL ULASV2U(BUFF2,TEXT2,KSTP,KPER,PERTIM,TOTIM,NSTRT,
     1            NNDLAY,KK,ImdtCF2,CSPCFM,LBHDSV,MDFLAG(NSTRT),NODES)
        END IF
        IPFLG=1
   78   CONTINUE        
C
C5C-----SAVE COEFFICIENTS FOR CROSS SECTION.
      ELSE
        IF(IOFLG(1,6).NE.0) THEN
C5D-----SAVE AI1            
          WRITE(IOUT,74) ImdtCF,KSTP,KPER
          IF(CSPCFM.EQ.' ') THEN
             CALL ULASAVU(BUFF,TEXT1,KSTP,KPER,PERTIM,TOTIM,NSTRT,
     1                NNDLAY,-1,ImdtCF,NODES)
          ELSE
             CALL ULASV2U(BUFF,TEXT1,KSTP,KPER,PERTIM,TOTIM,NSTRT,
     1                  NNDLAY,-1,ImdtCF,CSPCFM,LBHDSV,MDFLAG,NODES)
          END IF
C5E-----SAVE AI2          
          WRITE(IOUT,73) ImdtCF2,KSTP,KPER
          IF(CSPCFM.EQ.' ') THEN
             CALL ULASAVU(BUFF2,TEXT2,KSTP,KPER,PERTIM,TOTIM,NSTRT,
     1                NNDLAY,-1,ImdtCF2,NODES)
          ELSE
             CALL ULASV2U(BUFF2,TEXT2,KSTP,KPER,PERTIM,TOTIM,NSTRT,
     1                  NNDLAY,-1,ImdtCF2,CSPCFM,LBHDSV,MDFLAG,NODES)
          END IF          
          IPFLG=1
        END IF
      END IF
C
C6------RETURN.
   80 CONTINUE
      DEALLOCATE (BUFF2)
      RETURN
C
      END
C --------------------------------------------------------------------------
      SUBROUTINE GET_TEXTAI1(ICOMP,TEXT)
C     ******************************************************************
C     SELECT TEXT FOR OUTPUT FILE DEPENDING ON SPECIES NUMBER
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      
      USE GWTBCTMODULE, ONLY: MCOMP
      CHARACTER*16 TEXT
      CHARACTER*10 FORMATSTR  !kkz for a formatting string for multiple component numbers below
C     ------------------------------------------------------------------
C
      IF(MCOMP.EQ.1) THEN
        TEXT = '     COEFF AI1MD'  
      ELSE
        !kkz - write species number with leading zero if necessary
        IF(ICOMP.LT.10) THEN
          FORMATSTR='(A15,I1)'
          WRITE(TEXT,FORMATSTR)' COEFF AI1MD 0',ICOMP
        ELSE
          FORMATSTR='(A14,I2)'
          WRITE(TEXT,FORMATSTR)' COEFF AI1MD ',ICOMP
        ENDIF
        ! IF(ICOMP.EQ.1) THEN  
          ! TEXT = ' COEFF AI1MD 01'
        ! ELSEIF(ICOMP.EQ.2) THEN 
          ! TEXT = ' COEFF AI1MD 02'
        ! ELSEIF(ICOMP.EQ.3) THEN             
          ! TEXT = ' COEFF AI1MD 03'
        ! ELSEIF(ICOMP.EQ.4) THEN             
          ! TEXT = ' COEFF AI1MD 04'
        ! ELSEIF(ICOMP.EQ.5) THEN             
          ! TEXT = ' COEFF AI1MD 05'
        ! ELSEIF(ICOMP.EQ.6) THEN             
          ! TEXT = ' COEFF AI1MD 06'
        ! ELSEIF(ICOMP.EQ.7) THEN             
          ! TEXT = ' COEFF AI1MD 07'
        ! ELSEIF(ICOMP.EQ.8) THEN             
          ! TEXT = ' COEFF AI1MD 08'
        ! ELSEIF(ICOMP.EQ.9) THEN             
          ! TEXT = ' COEFF AI1MD 09'
        ! ELSEIF(ICOMP.EQ.10) THEN             
         ! TEXT = ' COEFF AI1MD 10'
        ! ELSEIF(ICOMP.EQ.11) THEN             
          ! TEXT = ' COEFF AI1MD 11'
        ! ELSEIF(ICOMP.EQ.12) THEN               
          ! TEXT = ' COEFF AI1MD 12'
        ! ELSEIF(ICOMP.EQ.13) THEN   
          ! TEXT = ' COEFF AI1MD 13'
        ! ELSEIF(ICOMP.EQ.14) THEN             
          ! TEXT = ' COEFF AI1MD 14'
        ! ELSEIF(ICOMP.EQ.15) THEN               
          ! TEXT = ' COEFF AI1MD 15'
        ! ELSEIF(ICOMP.EQ.16) THEN             
          ! TEXT = ' COEFF AI1MD 16'
        ! ELSEIF(ICOMP.EQ.17) THEN             
          ! TEXT = ' COEFF AI1MD 17'
        ! ELSEIF(ICOMP.EQ.18) THEN               
          ! TEXT = ' COEFF AI1MD 18'
        ! ELSEIF(ICOMP.EQ.19) THEN             
          ! TEXT = ' COEFF AI1MD 19'
        ! ELSEIF(ICOMP.EQ.20) THEN             
          ! TEXT = ' COEFF AI1MD 20'
        ! ENDIF             
      ENDIF            
C
C6------RETURN.
      RETURN
      END
C --------------------------------------------------------------------------
      SUBROUTINE GET_TEXTAI2(ICOMP,TEXT)
C     ******************************************************************
C     SELECT TEXT FOR OUTPUT FILE DEPENDING ON SPECIES NUMBER
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      
      USE GWTBCTMODULE, ONLY: MCOMP
      CHARACTER*16 TEXT
      CHARACTER*10 FORMATSTR  !kkz for a formatting string for multiple component numbers below
C     ------------------------------------------------------------------
C
      IF(MCOMP.EQ.1) THEN
        TEXT = '     COEFF AI2MD'  
      ELSE
        !kkz - write species number with leading zero if necessary
        IF(ICOMP.LT.10) THEN
          FORMATSTR='(A15,I1)'
          WRITE(TEXT,FORMATSTR)' COEFF AI2MD 0',ICOMP
        ELSE
          FORMATSTR='(A14,I2)'
          WRITE(TEXT,FORMATSTR)' COEFF AI2MD ',ICOMP
        ENDIF
        ! IF(ICOMP.EQ.1) THEN  
          ! TEXT = ' COEFF AI2MD 01'
        ! ELSEIF(ICOMP.EQ.2) THEN 
          ! TEXT = ' COEFF AI2MD 02'
        ! ELSEIF(ICOMP.EQ.3) THEN             
          ! TEXT = ' COEFF AI2MD 03'
        ! ELSEIF(ICOMP.EQ.4) THEN             
          ! TEXT = ' COEFF AI2MD 04'
        ! ELSEIF(ICOMP.EQ.5) THEN             
          ! TEXT = ' COEFF AI2MD 05'
        ! ELSEIF(ICOMP.EQ.6) THEN             
          ! TEXT = ' COEFF AI2MD 06'
        ! ELSEIF(ICOMP.EQ.7) THEN             
          ! TEXT = ' COEFF AI2MD 07'
        ! ELSEIF(ICOMP.EQ.8) THEN             
          ! TEXT = ' COEFF AI2MD 08'
        ! ELSEIF(ICOMP.EQ.9) THEN             
          ! TEXT = ' COEFF AI2MD 09'
        ! ELSEIF(ICOMP.EQ.10) THEN             
          ! TEXT = ' COEFF AI2MD 10'
        ! ELSEIF(ICOMP.EQ.11) THEN             
          ! TEXT = ' COEFF AI2MD 11'
        ! ELSEIF(ICOMP.EQ.12) THEN               
          ! TEXT = ' COEFF AI2MD 12'
        ! ELSEIF(ICOMP.EQ.13) THEN   
          ! TEXT = ' COEFF AI2MD 13'
        ! ELSEIF(ICOMP.EQ.14) THEN             
          ! TEXT = ' COEFF AI2MD 14'
        ! ELSEIF(ICOMP.EQ.15) THEN               
          ! TEXT = ' COEFF AI2MD 15'
        ! ELSEIF(ICOMP.EQ.16) THEN             
          ! TEXT = ' COEFF AI2MD 16'
        ! ELSEIF(ICOMP.EQ.17) THEN             
          ! TEXT = ' COEFF AI2MD 17'
        ! ELSEIF(ICOMP.EQ.18) THEN               
          ! TEXT = ' COEFF AI2MD 18'
        ! ELSEIF(ICOMP.EQ.19) THEN             
          ! TEXT = ' COEFF AI2MD 19'
        ! ELSEIF(ICOMP.EQ.20) THEN             
          ! TEXT = ' COEFF AI2MD 20'
        ! ENDIF             
      ENDIF            
C
C6------RETURN.
      RETURN
      END




C     Last change:  SP  Aug 25,2022   
C     *****************************************************************************************
C     MAIN CODE FOR UpdateState for IHM-USG
C     Created from UpdateState for IHM-96 - same functionality but on 
C     unstructured grids      
c**********************************************************************************************
      subroutine UpdateState(StressPeriod,IOUT,IUIHM)
c
c This subroutine is executed at the start of each stress period in MODFLOW-USG. The purpose is to 
c pause MODFLOW's execution in order to update some input files by the interprocessors in the integrated
c model.  
c      
c It starts by closing the input and output files that are not archived. Then a 'READ' statement is used to hold
c the model's execution.  After something is read, the I/O files are re-opened.
c
c StessPeriod       Stress period number. Only used for identification.
c IOUT              Output unit number for standarad modflow output
c IUIHM             Output unit number for debug file of IHM-USG      
c
c**********************************************************************************************
      ! implicit none
!      integer      ITLKCB,IWELCB,IDRNCB,IHEDUN,IDDNUN,
!     &             IRIVCB,IEVTCB,IGHBCB,IRCHCB,IRESCB,ISTCB1,ISTCB2,
!     &             IIBSCB,IFHBCB,IOUT,StressPeriod,nCol,nRow,nLay,fp
      
      USE GLOBAL, ONLY: NODES, AREA, NIUNIT, IUNIT
      USE GWFBCFMODULE, ONLY: SC2
      USE GWFRIVMODULE, ONLY:IRIVCB      
      USE NAMEFILEMODULE 
C
      INTEGER      lUnit(11),cnt,unitno
      LOGICAL      ftest(100),fouttest(100),FILETEST
      character*64 FNAME(40),FOUTNAME(40),FFORM(40),FOUTFORM(40),
     &             FILEACCESS(40),FOUTACCESS(40),SY_fname
      INTEGER      fp,StressPeriod
C-----------------------------------------------------------------------
1087  format(A8,A12,A7,A30,A11,I4,A16,A12)

c1      ! initialize all file tests to be false
      do cnt=1,100
         ftest(cnt) = .FALSE.
      end do
c2      ! write stress period information to screen and output listing file
      if(StressPeriod .eq. 1) then  ! first stress period - just finished launching
            write (*,*) "(MODFLOW Successfully Launched)"
            write (IOUT,*) "(MODFLOW Successfully Launched)"
c     
                IF(IUIHM.GT.0)
     *    WRITE(IUIHM,*)'(MODFLOW Successfully Launched)'
      else
            write (*,*) "(MODFLOW StressPeriod ",StressPeriod, ")"
            write (IOUT,*) "(MODFLOW StressPeriod ",StressPeriod, ")"
                      IF(IUIHM.GT.0)
     *    WRITE(IUIHM,*)"(MODFLOW StressPeriod ",StressPeriod, ")"
      endif
c
      write (*,*) "(MODFLOW Ready to Run for Stress Period:",
     &        StressPeriod, ")"
      write (IOUT,*) "(MODFLOW Ready to Run for Stress Period:",
     &        StressPeriod, ")"
                      IF(IUIHM.GT.0)
     *    WRITE(IUIHM,*) "(MODFLOW Ready to Run for Stress Period:",
     &        StressPeriod, ")"
c-------------------------------------------------------------------------------
c3 Close all open packages to be updated and not archived
      do cnt=1,nfiles
         IF(IARCVs(cnt) .EQ. 0) then    ! if NOT archiving, close this file before the stress period calculations start
             write(IOUT,*)' From 1 - Closing Unit ',IUS(cnt)
             call Close_if_opened(IUS(cnt), CNT,ftest(cnt))
             write(IOUT,1087)' Closed ',fmtargs(cnt),
     &        ' File: ',fnames(cnt),
     &        ' on Unit: ',ius(cnt),' with Access: ',accargs(cnt)
         endif
      end do
c------------------------------------------------------------------------------
c4 Send Message - wait for files to be updated

      write (*,*)"(MODFLOW--Need to Update Files)"
c
c4a Wait for Message - files updated (type and name)
       READ(*,*)
c------------------------------------------------------------------------------
c5 Open Files which have been updated in NON-archive mode; archive files are kept open so MF-USG can read from correct location

      do cnt=1,nfiles
         IF(IARCVs(cnt) .EQ. 0) then    ! if NOT archived, the file has been over-written so open it back and reposition for read files
             call Open_if_opened(IUS(cnt), CNT, ftest(cnt))
             write(IOUT,1087)' Opened ',fmtargs(cnt),
     &        ' File: ',fnames(cnt),
     &        ' on Unit: ',ius(cnt),' with Access: ',accargs(cnt)
         endif
      end do
c------------------------------------------------------------------------------
c6 update Specific Yield - if needed
      IF(SYIU .GT. 0)  call Update_Specific_Yield(
     &             SC2(1),NODES,AREA,StressPeriod,IOUT)      
c7 return
      Return 
      End

c**********************************************************************************************
      Subroutine Close_if_opened(unitno, CNT,filetest)
c     This subroutine closes a file connected to a unit number.  
c**********************************************************************************************
      USE NAMEFILEMODULE
      implicit none
      LOGICAL      filetest
      integer cnt,unitno,iflen
      CHARACTER*300 FNAME
c------------------------------------------------------------------------------
      fname = fnames (cnt)
      iflen = iflens(cnt)
      inquire(unit=unitno,NAME=fname(1:iflen),opened=filetest,
     &        FORM=fmtargs(cnt), ACCESS=accargs(cnt))
      IF(filetest) then
          close (unit=ius(cnt))
      endif
      RETURN
c
      End
c**********************************************************************************************
      Subroutine Open_if_opened(unitno, CNT, filetest)
c     This subroutine opens a file that used to be connected to a unit number. 
c     Also, for files with READ access, the file is positioned after 1 line to read SP information
c**********************************************************************************************
      USE NAMEFILEMODULE
      use GLOBAL, only: iout
      implicit none
      character*64 fileform,fileaccess
      integer      unitno,fileposition,int1,int2,cnt,iflen,m1,i1 
      LOGICAL      filetest, NotOpened
      CHARACTER*300 FNAME
      character*200 line
c------------------------------------------------------------------------------
C1 -------OPEN FILE IF IT WAS CLOSED FOR NON-ARCHIVE MODE
      IF(filetest) THEN        
        fname = fnames (cnt)
        iflen = iflens(cnt)
        open(unit=unitno,file=fname(1:iflen),status='old',
     &         form=fmtargs(cnt),err = 10,share='DENYNONE',
     &         ACCESS=accargs(cnt),ACTION=FILACTS(cnt))
C---------------------------------------------------------------
C2 -------POSITION READ FILES AFTER 1 DATA LINE FOR STRESS PERIOD INFORMATION        
        IF(FILACTS(cnt). EQ . 'READ ') THEN
C2A -------IF READ FILE THEN POSITION CORRECTLY TO READ NEXT STRESS PERIOD VALUES          
          IF(FMTARGS(cnt). NE. 'BINARY ') THEN 
C
C3 -------READ ASCII BOUNDARY FILE HEADER LINE SO AS TO POSITION FILE AT STRESS PERIOD
C---------EXCEPT FOR THE SYF FILE               
            IF(UNITNO.NE.SYIU) CALL URDCOM(unitno,IOUT,LINE)   !skip first line of data (urdcom removes comments)
          ENDIF 
        ENDIF
C ---------------------------------------------------------------
      ENDIF
C4 ----RETURN IF DONE      
      RETURN
C5 ----REPORT ERROR CONDITION AND STOP IF FILE DID NOT OPEN      
   10 continue
      WRITE(IOUT,55)UNITNO
55    FORMAT(1X,'COULD NOT OPEN FILE ON UNIT',I5,'STOPPING')
      STOP
      End

c**********************************************************************************************
      subroutine Update_Specific_Yield(SC2,NODES,AREA,StressPeriod,IOUT)
c This subroutine updates the specific yield values for the water table layer
c by reading values from a file and replacing the contenets of the array.
c The file should be formatted to be read by the U2DREL MODFLOW subroutine.
c
c SY_fname          File name containing the new values for the SY array
c SY_unit           Unit number on which the file will be opened
c StressPeriod      Stress period number - only used for identification purposes
c IOUT              Output unit number for standarad modflow output
c
c**********************************************************************************************
c
      use GLOBAL, only: nodlay,NROW,NCOL,IUNSTR,NODLAY,NLAY,ISYALL
      USE GWFBCFMODULE, ONLY: LAYCON 
      USE NAMEFILEMODULE 
      REAL,  DIMENSION(:,:),    ALLOCATABLE ::TEMP
      integer      SY_unit,NODES,StressPeriod,ios
      real         SC2(NODES)
      double precision AREA(NODES) 
      CHARACTER*24 ANAME
      DATA ANAME /'  SECONDARY STORAGE COEF'/
c------------------------------------------------------------------------------
c
c1 write indentifcation message to MODFLOW standard output
      WRITE(IOUT,*)' Updating Specific Yield Array for Stress Period ',
     &             StressPeriod
c2 -----set layer looping depending on ISYALL option 
      IF(ISYALL.EQ.0) THEN 
        KLAY = 1
      ELSE
        KLAY = NLAY
      ENDIF
C3 ------fill sc2 array     
      IF(IUNSTR.EQ.0)THEN
        ALLOCATE(TEMP(NCOL,NROW))  
        DO K = 1,KLAY
          IF(LAYCON(K).EQ.2. OR. LAYCON(K).EQ.3. OR. 
     1       LAYCON(K).EQ.4) THEN           
            CALL U2DREL(TEMP(1,1),ANAME,NROW,NCOL,K,SYIU,IOUT)
            DO I=1,NROW
              DO J=1,NCOL
                N=J+(I-1)*NCOL+(K-1)*NROW*NCOL
                SC2(N) = TEMP(J,I)
              ENDDO
            ENDDO
          ENDIF  
        ENDDO
        DEALLOCATE(TEMP)
      ELSE
        DO K = 1,KLAY  
          IF(LAYCON(K).EQ.2. OR. LAYCON(K).EQ.3. OR. 
     1       LAYCON(K).EQ.4) THEN  
            NNDLAY = NODLAY(K)
            NSTRT = NODLAY(K-1)+1
            NDSLAY = NNDLAY - NODLAY(K-1)
            CALL U1DREL(SC2(NSTRT),ANAME,NDSLAY,K,SYIU,IOUT)
          ENDIF
        ENDDO  
      ENDIF  
c4 -----multiply SC2 by area to preprocess the array for FM routines      
      DO K = 1,KLAY  
          IF(LAYCON(K).EQ.2. OR. LAYCON(K).EQ.3. OR. 
     1       LAYCON(K).EQ.4) THEN    
          NNDLAY = NODLAY(K)
          NSTRT = NODLAY(K-1)+1
          NDSLAY = NNDLAY - NODLAY(K-1)      
          do INODE=NSTRT,NNDLAY
             SC2(INODE)=SC2(INODE)*AREA(INODE) 
          end do
        ENDIF  
      ENDDO  
      WRITE(IOUT,*)' Done with Updating Specific Yield Array'
c4 done - 
C
      end subroutine
c--------------------------------------------------------------------------

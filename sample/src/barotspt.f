      PROGRAM BAROTSPT
      IMPLICIT NONE
C
      INTEGER    MKDM,MXT
      PARAMETER (MKDM=99, MXT=999)  ! big enough for all reasonable cases
C
      REAL         YP(MKDM)
      REAL*8       SUMYP
C
C**********
C*
C  1) CALCULATE BAROTROPOC TRANSPORT FROM HYCOM LAYER TRANSPORTS.
C
C  2) FROM THE TRANSPORT FILE HEADER:
C
C       TSPTLT - 25 CHARACTER TITLE FOR EACH SECTION OR REGION
C                  ='@0': SKIP THIS SECTION
C                  ='@+': ADD TO NEXT SECTION
C                  ='@-': SUBTRACT FROM NEXT SECTION
C
C     MULTIPLE SECTIONS CAN BE ADDED TOGETHER BY SETTING TSPTLT TO
C      '@+' OR '@-' FOR SEVERAL SECTIONS IN SEQUENCE .  IN SUCH CASES 
C      THE SUM PROPAGATES TO THE FIRST FOLLOWING SECTION THAT DOES NOT 
C      HAVE TSPTLT STARTING WITH '@'.
C
C  4) THE LAYER TRANSPORT FILE IS INPUT ON UNIT 10, AND THE RESULTS
C      OUTPUT ON UNIT 20.
C
C  5) ALAN J. WALLCRAFT, SEPT 2019.
C*
C**********
C
      LOGICAL       LPOSITIVE
      INTEGER       I,ICOUNT,ISUM,J,JJ,K,KB,KDM,IEXPT,MSECT,NPAGE
      REAL          DAYF,DAYL,DAYFQ,YQ
      CHARACTER*240 FLNM_T,FLNM_B,CLINE
C
      INTEGER       LP
      COMMON/LINEPR/LP
C
      LP = 6
C
C --- 'flnm_t' = name of transport section sample input file
C --- 'flnm_b' = name of transport section  baro output file
C
      read (*,'(a)') flnm_t
      write (lp,'(2a)') 'transport  input file: ',
     &                  flnm_t(1:len_trim(flnm_t))
      call flush(lp)
      read (*,'(a)') flnm_b
      write (lp,'(2a)') 'transport output file: ',
     &                  flnm_b(1:len_trim(flnm_b))
      write(lp,*)
C
C     READ AND WRITE TRANSPORT FILE HEADER RECORDS.
C
      OPEN(UNIT=20,FILE=FLNM_B,FORM='formatted',
     &      STATUS='new',ACTION='write')
C
      OPEN(UNIT=10,FILE=FLNM_T,FORM='formatted',
     &      STATUS='old',ACTION='read')
C
       READ(10,'(a)')      CLINE   !line 1
      WRITE(20,'(a)') TRIM(CLINE)
       READ(10,'(a)')      CLINE   !line 2
       READ(CLINE,6401) KDM,IEXPT
      WRITE(20,   3001)   1,IEXPT
       READ(10,'(a)')      CLINE   !line 3
       READ(CLINE,6402) DAYF,DAYL,DAYFQ
      WRITE(20,'(a)') TRIM(CLINE)
       READ(10,'(a)')      CLINE   !line 4
       READ(CLINE,6403) MSECT
      WRITE(20,'(a)') TRIM(CLINE)
C
      IF    (KDM .GT. MKDM) THEN
        WRITE(LP,9100)
        STOP
      ENDIF
C
C     READ TRANSPORT LOCATIONS AND NAMES
C
      DO 10 J= 1,MSECT
         READ(10,'(a)')      CLINE
        WRITE(20,'(a)') TRIM(CLINE)
   10 CONTINUE
C
C     SUM TRANSPORTS.
C
      DO 20 I=1,99999
        DO J= 1,MSECT
          DO K= 1,KDM,10
            READ(10,6420,END=220) YP(K:MIN(K+9,KDM))
          ENDDO
          SUMYP = 0.D0
          DO K= 1,KDM
            SUMYP = SUMYP + YP(K)
          ENDDO
          WRITE(20,3002) J,SUMYP
        ENDDO
   20 CONTINUE
  220 CONTINUE
C
      WRITE(LP,'(/A/)') 'NORMAL EXIT FROM BAROTSPT'
      CALL FLUSH(LP)
      STOP
C
 1000 FORMAT(/ ' MODEL EXPERIMENT:',I5)
 1050 FORMAT(1X,I5,1X,I2,'-LAYER ',A,' STATISTICS')
 1075 FORMAT(/ T33,'LONGITUDE',
     +         T54,'LATITUDE',
     +         T71,'LAYER',
     +         T82,'MEAN (SV)',
     +         T95,'STD DEV' /)
 1100 FORMAT(  T2, A25,': ',2(F7.2,A1,1X),2X,2(F7.2,A1,1X),
     +         T72,I2.2,
     +         T81,F7.2,
     +         T94,F7.2)
 1105 FORMAT(  T2, A25,': ',2(F7.2,A1,1X),2X,2(F7.2,A1,1X),
     +         T71,I2.2,'-',I2.2,
     +         T81,F7.2,
     +         T94,F7.2)
 1107 FORMAT(  T2, A25,': ',2(F7.2,A1,1X),2X,2(F7.2,A1,1X),
     +         T71,' +VE ',
     +         T81,F7.2,
     +         T94,F7.2)
 1110 FORMAT(  T2, A25,
     +         T72,I2.2,
     +         T81,F7.2,
     +         T94,F7.2)
 1115 FORMAT(  T2, A25,
     +         T71,I2.2,'-',I2.2,
     +         T81,F7.2,
     +         T94,F7.2)
 1117 FORMAT(  T2, A25,
     +         T71,' +VE ',
     +         T81,F7.2,
     +         T94,F7.2)
 1120 FORMAT(  T71,'----------','----------','----------')
 1200 FORMAT(1X,F7.1,' DAY AVERAGE OF',I5,' SAMPLES BETWEEN',
     +       ' MODEL DAYS',F8.1,' AND',F8.1 )
 1400 FORMAT(///)
 1500 FORMAT(A1 ///)
 3001 format(1x,'An',i3,' layer HYCOM experiment with label ',i4,'.')
 3002 format(i4,2x,f10.4)
 6401 FORMAT(3X,I3,35X,I4)
 6402 FORMAT(15X,F10.2,3X,F10.2,19X,F6.2)
 6403 FORMAT(I4)
 6415 FORMAT(6X,4F9.2,2X,A25)
 6420 FORMAT(6X,10F10.4)
 6600 FORMAT(1X)
 9100 FORMAT(10X,'**ERROR:  KDM .GT. MKDM**')
      END

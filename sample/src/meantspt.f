      PROGRAM MNTSPT
      IMPLICIT NONE
C
      INTEGER    MKDM,MXT
      PARAMETER (MKDM=99, MXT=999)  ! big enough for all reasonable cases
C
      CHARACTER*1  CLONG(2,MXT),CLAT(2,MXT)
      CHARACTER*25 TSPTLT(MXT)
      REAL         XYSECT(4,MXT), YP(MXT,MKDM)
      REAL         AVGTPT(MXT,MKDM), SUMTPT(MXT,MKDM)
      REAL         VARTPT(MXT,MKDM), STDTPT(MXT,MKDM)
      REAL         AVGTOT(MXT,2), TMPTOT(MXT), POSTOT(MXT)
      REAL         VARTOT(MXT,2), STDTOT(MXT,2)
C
C**********
C*
C  1) CALCULATE MEAN AND VARIABILITY VALUES FOR HYCOM TRANSPORTS.
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
C  4) THE TRANSPORT FILE IS INPUT ON UNIT 10, AND THE RESULTS
C      OUTPUT ON UNIT 20.
C
C  5) JOSEPH E. METZGER AND ALAN J. WALLCRAFT, JUNE 1991 AND MAY 1993.
C     MODIFIED FOR HYCOM BY ALAN J. WALLCRAFT, JULY 2001.
C*
C**********
C
      INTEGER    MKDMMXT
      PARAMETER (MKDMMXT=MKDM*MXT)
C
      LOGICAL      LPOSITIVE
      INTEGER      I,ICOUNT,ISUM,J,JJ,K,KB,KDM,IEXPT,MSECT,NPAGE
      INTEGER      KOUT,LAYBOT(0:99)
      REAL         DAYF,DAYL,DAYFQ,YQ
      CHARACTER*240 FLNM_T,FLNM_M
C
      INTEGER       LP
      COMMON/LINEPR/LP
C
      DATA SUMTPT / MKDMMXT*0.0 /
      DATA VARTPT / MKDMMXT*0.0 /
      DATA VARTOT /    MXT*0.0,    MXT*0.0 /
      DATA AVGTOT /    MXT*0.0,    MXT*0.0 /
C
      DATA TSPTLT / MXT*' ' /
C
      LP = 6
C
C --- 'flnm_t' = name of transport section sample input file
C --- 'flnm_m' = name of transport section  mean output file
C --- 'kout  ' = number of layer combinations to output (<=kdm)
C ---             <0 to output positive net instead of sum of layers 1:kdm-1
C
      read (*,'(a)') flnm_t
      write (lp,'(2a)') 'transport  input file: ',
     &                  flnm_t(1:len_trim(flnm_t))
      call flush(lp)
      read (*,'(a)') flnm_m
      write (lp,'(2a)') 'transport output file: ',
     &                  flnm_m(1:len_trim(flnm_m))
      call flush(lp)
      call blkini(kout,  'kout  ')
      write(lp,*)
C
      LPOSITIVE = KOUT .LT. 0
      KOUT = ABS(KOUT)
C
      IF     (LPOSITIVE) THEN
        WRITE(LP,'(a)')
     &    ' +VE row is vertical sum of all net positive mean transports'
        WRITE(LP,*)
      ENDIF
C
C     READ TRANSPORT FILE HEADER RECORDS.
C
      OPEN(UNIT=10,FILE=FLNM_T,FORM='formatted',
     &      STATUS='old',ACTION='read')
C
      READ(10,6400) KDM,IEXPT,DAYF,DAYL,DAYFQ,MSECT
C
      IF    (KDM .GT. MKDM) THEN
        WRITE(LP,9100)
        STOP
      ENDIF
C
      NPAGE = 56/(KOUT+4)
C
      OPEN(UNIT=20,FILE=FLNM_M,FORM='formatted',
     &      STATUS='new',ACTION='write')
C
      WRITE(20,1000) IEXPT
      WRITE(20,1050) MSECT, KDM, 'TRANSPORT SECTION'
      CALL FLUSH(20)
      WRITE(LP,1000) IEXPT
      WRITE(LP,1050) MSECT, KDM, 'TRANSPORT SECTION'
      CALL FLUSH(LP)
C
      IF     (KOUT.EQ.KDM) THEN
C
C ---   output all layers.
C
        DO K= 0,KOUT
          LAYBOT(K) = K
        ENDDO
      ELSEIF (KOUT.GE.3 .AND. KOUT.LT.KDM) THEN
        LAYBOT(0) = 0
        DO K= 1,KOUT
C
C ---     'laybot' = last layer in next combination of layers
C
          call blkini(laybot(k),  'laybot')
C
          IF     (LAYBOT(K).LE.LAYBOT(K-1)) THEN
            WRITE(LP,'(A)') 'ERROR - LAYBOT MUST BE ASCENDING'
            STOP
          ELSEIF (K.EQ.KOUT .AND. LAYBOT(K).NE.KDM) THEN
            WRITE(LP,'(A)') 'ERROR - LAST LAYBOT MUST BE KDM'
            STOP
          ELSEIF (K.NE.KOUT .AND. LAYBOT(K).GE.KDM) THEN
            WRITE(LP,'(A)') 'ERROR - LAYBOT REACHED KDM TOO SOON'
            STOP
          ENDIF
        ENDDO
      ELSE
        WRITE(LP,'(A)') 'ERROR - BAD KOUT VALUE'
        CALL FLUSH(LP)
        STOP
      ENDIF
C
C     READ TRANSPORT LOCATIONS AND NAMES
C
      DO 10 J= 1,MSECT
        READ(10,6415) (XYSECT(I,J),I=1,4),TSPTLT(J)
   10 CONTINUE
C
C     PREPARE LOCATION INFO FOR PRINTING.
C
      ISUM = 0
      DO 15 J= 1,MSECT
        CLONG(1,J) = 'E'
        CLONG(2,J) = 'E'
        CLAT (1,J) = 'N'
        CLAT (2,J) = 'N'
C
        IF (XYSECT(1,J) .GE. 360.0) THEN
          XYSECT(1,J) = XYSECT(1,J) - 360.0
        ELSE IF (XYSECT(1,J) .GT. 180.0) THEN
          CLONG(1,J) = 'W'
          XYSECT(1,J) = XYSECT(1,J) - 180.0
          XYSECT(1,J) = 180.0 - XYSECT(1,J)
        ENDIF
        IF (XYSECT(2,J) .GE. 360.0) THEN
          XYSECT(2,J) = XYSECT(2,J) - 360.0
        ELSE IF (XYSECT(2,J) .GT. 180.0) THEN
          CLONG(2,J) = 'W'
          XYSECT(2,J) = XYSECT(2,J) - 180.0
          XYSECT(2,J) = 180.0 - XYSECT(2,J)
        ENDIF
        IF (XYSECT(3,J) .LT. 0.0) THEN
          CLAT(1,J) = 'S'
          XYSECT(3,J) = ABS(XYSECT(3,J))
        ENDIF
        IF (XYSECT(4,J) .LT. 0.0) THEN
          CLAT(2,J) = 'S'
          XYSECT(4,J) = ABS(XYSECT(4,J))
        ENDIF
C
C       WEED OUT SUMMED SECTIONS.
C
        IF     (TSPTLT(J)(1:1).EQ.'@') THEN
          IF     (TSPTLT(J)(2:2).NE.'0') THEN
            ISUM = ISUM + 1
          ENDIF
        ELSE
          IF     (ISUM.NE.0) THEN
            ISUM = 0
            CLONG(1,J) = 'X'
          ENDIF
        ENDIF
   15 CONTINUE
C
C     SUM TRANSPORTS.
C
      ICOUNT =  0
      DO 20 I=1,99999
        ICOUNT = ICOUNT + 1
        DO J= 1,MSECT
          DO K= 1,KDM,10
            READ(10,6420,END=220) YP(J,K:MIN(K+9,KDM))
          ENDDO
        ENDDO
        DO K= 1,KOUT
          IF     (LAYBOT(K-1)+1.NE.LAYBOT(K)) THEN
C
C           COMBINE LAYERS.
C
            DO J= 1,MSECT
              YP(J,LAYBOT(K)) = SUM(YP(J,LAYBOT(K-1)+1:LAYBOT(K)))
              YP(J,LAYBOT(K-1)+1:LAYBOT(K)-1) = 0.0
            ENDDO
          ENDIF
        ENDDO
        DO 30 K=1,KDM
C
C         DO ANY REQUIRED SUMS.
C
          ISUM = 0
          YQ   = 0.0
          DO 31 J= 1,MSECT
            IF     (TSPTLT(J)(1:2).EQ.'@+') THEN
              ISUM = ISUM + 1
              YQ   = YQ + YP(J,K)
            ELSEIF (TSPTLT(J)(1:2).EQ.'@-') THEN
              ISUM = ISUM + 1
              YQ   = YQ - YP(J,K)
            ELSEIF (TSPTLT(J)(1:1).NE.'@') THEN
              IF     (ISUM.NE.0) THEN
                YP(J,K) = YP(J,K) + YQ
                ISUM    = 0
                YQ      = 0.0
              ENDIF
            ENDIF
   31     CONTINUE
C
          IF     (ICOUNT.GT.0) THEN
            DO 35 J=1,MSECT
              SUMTPT(J,K) = SUMTPT(J,K) + YP(J,K)
   35       CONTINUE
          ENDIF
   30   CONTINUE
   20 CONTINUE
  220 CONTINUE
      ICOUNT = ICOUNT - 1
C
C     AVERAGE TRANSPORT FOR EACH LAYER
C
      DO 50 J=1,MSECT
        DO 55 K=1,KDM
          AVGTPT(J,K) = SUMTPT(J,K) / FLOAT(ICOUNT)
   55   CONTINUE
   50 CONTINUE
C
C     SUM ALL LAYERS AND ALL BUT LAST LAYER COMBINATION.
C
      DO J=1,MSECT
        IF     (LPOSITIVE) THEN
          DO K=1,KDM
            AVGTOT(J,1) = AVGTOT(J,1) + MAX( 0.0, AVGTPT(J,K) )
            AVGTOT(J,2) = AVGTOT(J,2) +           AVGTPT(J,K)
          ENDDO !i
        ELSE
          DO K=1,KDM-1
            AVGTOT(J,1) = AVGTOT(J,1) + AVGTPT(J,K)
          ENDDO !i
          AVGTOT(J,2) = AVGTOT(J,1) + AVGTPT(J,KDM)
        ENDIF
      ENDDO !j

      DAYL = DAYF + ((ICOUNT-1) * DAYFQ)
      WRITE(20,1200) ICOUNT*DAYFQ,ICOUNT,DAYF,DAYL
      CALL FLUSH(20)
      WRITE(LP,1200) ICOUNT*DAYFQ,ICOUNT,DAYF,DAYL
      CALL FLUSH(LP)
C
C     COMPUTE STANDARD DEVIATION
C
      REWIND 10
      DO 56 I=1,4
        READ(10,6600)
   56 CONTINUE
      DO 57 I=1,MSECT
        READ(10,6600) 
   57 CONTINUE

      DO 60 I=1,ICOUNT
        DO J= 1,MSECT
          DO K= 1,KDM,10
            READ(10,6420) YP(J,K:MIN(K+9,KDM))
          ENDDO
        ENDDO
        DO K= 1,KOUT
          IF     (LAYBOT(K-1)+1.NE.LAYBOT(K)) THEN
            DO J= 1,MSECT
              YP(J,LAYBOT(K)) = SUM(YP(J,LAYBOT(K-1)+1:LAYBOT(K)))
              YP(J,LAYBOT(K-1)+1:LAYBOT(K)-1) = 0.0
            ENDDO
          ENDIF
        ENDDO
        DO 65 K=1,KDM
          ISUM = 0
          YQ   = 0.0
          DO 66 J= 1,MSECT
            IF     (TSPTLT(J)(1:2).EQ.'@+') THEN
              ISUM = ISUM + 1
              YQ   = YQ + YP(J,K)
            ELSEIF (TSPTLT(J)(1:2).EQ.'@-') THEN
              ISUM = ISUM + 1
              YQ   = YQ - YP(J,K)
            ELSEIF (TSPTLT(J)(1:1).NE.'@') THEN
              IF     (ISUM.NE.0) THEN
                YP(J,K) = YP(J,K) + YQ
                ISUM    = 0
                YQ      = 0.0
              ENDIF
            ENDIF
   66     CONTINUE
C
          DO 67 J=1,MSECT
            IF     (K.EQ.1) THEN
              TMPTOT(J) = 0.0
              POSTOT(J) = 0.0
            ENDIF
            IF     (AVGTPT(J,K).GE.0.0) THEN
              POSTOT(J) = POSTOT(J) + YP(J,K)
            ENDIF
            TMPTOT(J) = TMPTOT(J) + YP(J,K)
   67     CONTINUE
          DO 70 J=1,MSECT
            VARTPT(J,K) = VARTPT(J,K) + ((YP(J,K)-AVGTPT(J,K))**2)
            IF     (.NOT. LPOSITIVE .AND. K.EQ.KDM-1) THEN
              VARTOT(J,1) = VARTOT(J,1) + 
     +                        ((TMPTOT(J) - AVGTOT(J,1))**2)
            ELSEIF (      LPOSITIVE .AND. K.EQ.KDM)   THEN
              VARTOT(J,1) = VARTOT(J,1) + 
     +                        ((POSTOT(J) - AVGTOT(J,1))**2)
            ENDIF
            IF     (K.EQ.KDM) THEN
              VARTOT(J,2) = VARTOT(J,2) + 
     +                        ((TMPTOT(J) - AVGTOT(J,2))**2)
            ENDIF
   70     CONTINUE
   65   CONTINUE
   60 CONTINUE
C
      WRITE(20,1075)
C
      JJ = 1
      DO 90 J=1,MSECT
        IF     (JJ.GT.1 .AND. MOD(JJ,NPAGE).EQ.1) THEN
          WRITE(20,1500) CHAR(12)
          WRITE(20,1075)
        ENDIF
        IF     (TSPTLT(J)(1:1).NE.'@') THEN
          DO 95 KB=1,KOUT
            K = LAYBOT(KB)
            STDTPT(J,K) = SQRT(VARTPT(J,K)/FLOAT(ICOUNT))
            IF     (LAYBOT(KB-1)+1.EQ.K) THEN
              IF     (CLONG(1,J) .EQ.'X') THEN
                WRITE(20,1110) TSPTLT(J),K,
     +                         AVGTPT(J,K),STDTPT(J,K)
              ELSE
                WRITE(20,1100) TSPTLT(J),
     +                         XYSECT(1,J),CLONG(1,J),
     +                         XYSECT(2,J),CLONG(2,J),
     +                         XYSECT(3,J), CLAT(1,J),
     +                         XYSECT(4,J), CLAT(2,J),K,
     +                         AVGTPT(J,K),STDTPT(J,K)
              ENDIF
            ELSE
              IF     (CLONG(1,J) .EQ.'X') THEN
                WRITE(20,1115) TSPTLT(J),
     +                         LAYBOT(KB-1)+1,K,
     +                         AVGTPT(J,K),STDTPT(J,K)
              ELSE
                WRITE(20,1105) TSPTLT(J),
     +                         XYSECT(1,J),CLONG(1,J),
     +                         XYSECT(2,J),CLONG(2,J),
     +                         XYSECT(3,J), CLAT(1,J),
     +                         XYSECT(4,J), CLAT(2,J),
     +                         LAYBOT(KB-1)+1,K,
     +                         AVGTPT(J,K),STDTPT(J,K)
              ENDIF
            ENDIF
   95     CONTINUE
          STDTOT(J,1) = SQRT(VARTOT(J,1)/FLOAT(ICOUNT))
          STDTOT(J,2) = SQRT(VARTOT(J,2)/FLOAT(ICOUNT))
          WRITE(20,1120)
          IF     (LPOSITIVE) THEN
            IF     (CLONG(1,J) .EQ.'X') THEN
              WRITE(20,1117) TSPTLT(J),
     +                       AVGTOT(J,1),STDTOT(J,1)
            ELSE
              WRITE(20,1107) TSPTLT(J),
     +                       XYSECT(1,J),CLONG(1,J),
     +                       XYSECT(2,J),CLONG(2,J),
     +                       XYSECT(3,J), CLAT(1,J),
     +                       XYSECT(4,J), CLAT(2,J),
     +                       AVGTOT(J,1),STDTOT(J,1)
            ENDIF
          ELSE
            IF     (CLONG(1,J) .EQ.'X') THEN
              WRITE(20,1115) TSPTLT(J),
     +                       1,LAYBOT(KOUT-1),
     +                       AVGTOT(J,1),STDTOT(J,1)
            ELSE
              WRITE(20,1105) TSPTLT(J),
     +                       XYSECT(1,J),CLONG(1,J),
     +                       XYSECT(2,J),CLONG(2,J),
     +                       XYSECT(3,J), CLAT(1,J),
     +                       XYSECT(4,J), CLAT(2,J),
     +                       1,LAYBOT(KOUT-1),
     +                       AVGTOT(J,1),STDTOT(J,1)
            ENDIF
          ENDIF !lpositive:else
          IF     (CLONG(1,J) .EQ.'X') THEN
            WRITE(20,1115) TSPTLT(J),
     +                     1,KDM,
     +                     AVGTOT(J,2),STDTOT(J,2)
          ELSE
            WRITE(20,1105) TSPTLT(J),
     +                     XYSECT(1,J),CLONG(1,J),
     +                     XYSECT(2,J),CLONG(2,J),
     +                     XYSECT(3,J), CLAT(1,J),
     +                     XYSECT(4,J), CLAT(2,J),
     +                     1,KDM,
     +                     AVGTOT(J,2),STDTOT(J,2)
          ENDIF            
          WRITE(20,*)
          CALL FLUSH(20)
C
          JJ = JJ + 1
        ENDIF
   90 CONTINUE
C
      WRITE(20,1400)
      CALL FLUSH(20)
C
      WRITE(LP,'(/A/)') 'NORMAL EXIT FROM MEANTSPT'
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
 6400 FORMAT( / 3X,I3,35X,I4 / 15X,F10.2,3X,F10.2,19X,F6.2 / I4)
 6415 FORMAT(6X,4F9.2,2X,A25)
 6420 FORMAT(6X,10F10.4)
 6600 FORMAT(1X)
 9100 FORMAT(10X,'**ERROR:  KDM .GT. MKDM**')
      END

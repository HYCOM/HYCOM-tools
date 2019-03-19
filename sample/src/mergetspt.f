      PROGRAM MRGTSPT
      IMPLICIT NONE
C
C**********
C*
C  1) MERGE TWO HYCOM TRANSPORT FILES.
C*
C**********
C
      INTEGER       I
      INTEGER       KDM1,IEXPT1,MSECT1,
     &              KDM2,IEXPT2,MSECT2
      REAL          DAYF1,DAYL1,DAYFQ1,
     &              DAYF2,DAYL2,DAYFQ2
      CHARACTER*240  FLNM_1,FLNM_2,FLNM_O
      CHARACTER*120 CLINE1,CLINE2
C
      INTEGER       LP
      COMMON/LINEPR/LP
C
      LP = 6
C
C --- 'flnm_1' = name of 1st    transport section sample  input file
C --- 'flnm_2' = name of 2nd    transport section sample  input file
C --- 'flnm_o' = name of merged transport section sample output file
C
      read (*,'(a)') flnm_1
      write (lp,'(2a)') 'transport 1st input file: ',
     &                  flnm_1(1:len_trim(flnm_1))
      call flush(lp)
      read (*,'(a)') flnm_2
      write (lp,'(2a)') 'transport 2nd input file: ',
     &                  flnm_2(1:len_trim(flnm_2))
      call flush(lp)
      read (*,'(a)') flnm_o
      write (lp,'(2a)') 'transport    output file: ',
     &                  flnm_o(1:len_trim(flnm_o))
C
C     READ TRANSPORT FILE HEADER RECORDS.
C
      OPEN(UNIT=11,FILE=FLNM_1,FORM='formatted',
     &      STATUS='old',ACTION='read')
      OPEN(UNIT=12,FILE=FLNM_2,FORM='formatted',
     &      STATUS='old',ACTION='read')
C
      READ(11,6400) KDM1,IEXPT1,DAYF1,DAYL1,DAYFQ1,MSECT1
      READ(12,6400) KDM2,IEXPT2,DAYF2,DAYL2,DAYFQ2,MSECT2
C
      IF     (KDM1  .NE.KDM2   .OR.
     &        IEXPT1.NE.IEXPT2 .OR.
     &        DAYFQ1.NE.DAYFQ2 .OR.
     &        MSECT1.NE.MSECT2     ) THEN
        WRITE(LP,*)
        WRITE(LP,*) '***** ERROR - MISMATCHED HEADERS:'
        WRITE(LP,*) 'KDM   = ',KDM1,  KDM2
        WRITE(LP,*) 'IEXPT = ',IEXPT1,IEXPT2
        WRITE(LP,*) 'DAYFQ = ',DAYFQ1,DAYFQ2
        WRITE(LP,*) 'MSECT = ',MSECT1,MSECT2
        WRITE(LP,*)
        CALL FLUSH(LP)
        STOP
      ELSEIF (ABS( DAYF2 - (DAYL1+DAYFQ1) ).GT.DAYFQ1*0.001) THEN
        WRITE(LP,*)
        WRITE(LP,*) '***** ERROR - NON-SEQUENTIAL SAMPLE FILES'
        WRITE(LP,*)
        WRITE(LP,*) 'DAYF2 IS ',DAYF2,' BUT SHOULD BE ',DAYL1+DAYFQ1
        WRITE(LP,*)
        CALL FLUSH(LP)
        STOP
      ENDIF
C
      OPEN(UNIT=20,FILE=FLNM_O,FORM='formatted',
     &      STATUS='new',ACTION='write')
      WRITE(20,3000) KDM1,IEXPT1,DAYF1,DAYL2,DAYFQ1,MSECT1
C
      DO I= 1,MSECT1
        READ(11,'(A)') CLINE1
        READ(12,'(A)') CLINE2
        IF     (CLINE1.NE.CLINE2) THEN
          WRITE(LP,*)
          WRITE(LP,*) '***** ERROR - DIFFERING SECTIONS'
          WRITE(LP,*)
          WRITE(LP,*) CLINE1
          WRITE(LP,*) CLINE2
          WRITE(LP,*)
          CALL FLUSH(LP)
          STOP
        ENDIF
        WRITE(20,'(A)') CLINE1(1:LEN_TRIM(CLINE1))
      ENDDO
      CALL FLUSH(20)
C
C     OUTPUT ALL SAMPLES FROM THE 1ST FILE.
C
      DO
        READ(11,'(A)',IOSTAT=I) CLINE1
        IF     (I.NE.0) THEN
          EXIT
        ENDIF
        WRITE(20,'(A)') CLINE1(1:LEN_TRIM(CLINE1))
      ENDDO
      CALL FLUSH(20)
C
C     OUTPUT ALL SAMPLES FROM THE 2ND FILE.
C
      DO
        READ(12,'(A)',IOSTAT=I) CLINE2
        IF     (I.NE.0) THEN
          EXIT
        ENDIF
        WRITE(20,'(A)') CLINE2(1:LEN_TRIM(CLINE2))
      ENDDO
      CALL FLUSH(20)
C
      WRITE(LP,'(/A/)') 'NORMAL EXIT FROM MERGETSPT'
      CALL FLUSH(LP)
      STOP
C
 3000 format(1x,'Transport Sections' /
     +       1x,'An',i3,' layer HYCOM experiment with label ',i4,'.' /
     +       1x,'From model day',f10.2,' to',f10.2,
     +          ', with values every',f6.2,' days.' /
     +       i4,' lines of transports at locations:')
 6400 FORMAT( / 3X,I3,35X,I4 / 15X,F10.2,3X,F10.2,19X,F6.2 / I4)
      END

      PROGRAM HYCOM_SKILL
      IMPLICIT NONE
C
      CHARACTER*256    CHDFI1,CHDFI2,CHDFO,CHDFB, CLABEL
      INTEGER          IDM,JDM,ITLSDS1,LSTSDS1,INCSDS1,ITLSDS2,INCSDS2
      REAL             RMSE,NOISE,POWER
C
C**********
C*
C  1) ERROR STATISTICS BETWEEN TWO SEQUENCES OF HYCOM .a FIELDS.
C     BASED ON SDSCORR, WITH SOME COMMENTS STILL USING HDF NAMES.
C
C     CONDITIONAL AND UNCONDITIONAL BIAS EACH LIMITED TO 1.0,
C     SO MINIMUM SKILL SCORE IS -2.0.
C
C  2) ARGUMENT INPUT (FROM COMMAND LINE, SEE ARGS AND USAGE)
C
C      CHDFI1  - 1ST  INPUT HYCOM FILENAME
C      CHDFI2  - 2ND  INPUT HYCOM FILENAME
C      CHDFO   -     OUTPUT HYCOM FILENAME (MUST END IN .a)
C      IDM,JDM - ARRAY DIMENSIONS
C      CLABEL  - PARTIAL OUTPUT SDS LABEL (SEE 4)
C      RMSE    - TARGET RMS ERROR VALUE
C                  > 0.0 SETS FLOOR ON SS (NO NOISE)
C      NOISE   - RMS OF INDEPENDENT NOISE
C                  > 0.0 MODIFIED SS (NO RMSE)
C                 USE EXACTLY ONE OF RMSE AND NOISE
C      POWER   - RAISE INPUT FIELDS TO THIS POWER (0.0...1.0)
C                  =-1.0 TAKE -1/     OF THE INPUT
C                  =-0.5 TAKE -1/SQRT OF THE INPUT
C                  = 0.0 TAKE LOG  OF THE INPUT
C                  = 0.5 TAKE SQRT OF THE INPUT
C                  = 1.0 USE THE INPUT AS IS (DEFAULT)
C      ITLSDS1 - RELATIVE POSITION OF 1ST FIRST SDS SELECTED
C      LSTSDS1 - RELATIVE POSITION OF 1ST LAST  SDS SELECTED
C      INCSDS1 - 1ST INCREMENT BETWEEN SELECTED SDS'S
C                  IF INCSDS1.GT.0, MUST HAVE ITLSDS1.LE.LSTSDS1
C                  IF INCSDS1.LT.0, MUST HAVE ITLSDS1.GE.LSTSDS1
C      ITLSDS2 - RELATIVE POSITION OF 2ND FIRST SDS SELECTED
C      INCSDS2 - 2ND INCREMENT BETWEEN SELECTED SDS'S
C                  =0 TO PERPETUALLY REPEAT ITLSDS2
C
C  3) EIGHT SDS'S ARE PRODUCED.
C      A) MEAN ERROR (TITLE: MEAN ERROR for ...)
C      B) RMS ERROR  (TITLE: RMS ERROR for ...)
C      C) NORMALIZED MEAN ERROR (TITLE: 1/RMSE * MEAN ERROR for ...)
C      D) NORMALIZED RMS ERROR  (TITLE: 1/RMSE * RMS ERROR for ...)
C      E) COND.   BIAS (TITLE: Bcond for ...)
C      F) UNCOND. BIAS (TITLE: Buncond for ...)
C      G) SKILL SCORE  (TITLE: Skill Score for ...)
C      H) CORRELATION COEFFICIENT (TITLE: Correlation Coeff. for ...)
C
C     CONDITIONAL AND UNCONDITIONAL BIAS ATE EACH LIMITED TO 1.0,
C      SO MINIMUM SKILL SCORE IS -2.0.
C
C  4) THE FIRST OCCURENCE OF '@1' OR '@2' WITHIN 'CLABEL' WILL BE
C      REPLACED BY THE NON-BLANK PORTION OF THE 1ST OR 2ND INPUT SDS
C      LABEL RESPECTIVELY.
C
C  5) ALAN J. WALLCRAFT, NRL, JUNE 1998.
C*
C**********
C
      INTEGER  OTKMREAL
      EXTERNAL OTKMREAL
C
      REAL       SPVAL
      PARAMETER (SPVAL=2.0**100)
C
      INTEGER       ISDS,JSDS
      INTEGER       I,IOS,N,NDSDS,NPAD,NRECL
      REAL          QRMSE,RN,SDS1A,SDS2A,
     +              SDSAAX,SDSBBX,SDSDDX,SMIN,SMAX
C
      REAL              :: PAD(4096)
      REAL, ALLOCATABLE :: SDS1(:),SDSAA(:),SDSDD(:),
     +                     SDS2(:),SDSBB(:),SDSAB(:),
     +                    SDS1M(:),
     +                    SDS2M(:)
C
C     COMMAND LINE PROCESSING.
C
      CALL ARGS(CHDFI1,CHDFI2,CHDFO,CHDFB,
     +          IDM,JDM, CLABEL, RMSE, NOISE, POWER,
     +          ITLSDS1,LSTSDS1,INCSDS1,ITLSDS2,INCSDS2)
C
C     ALLOCATE SPACE FOR SEVERAL COMPLETE ARRAYS.
C
      NDSDS=IDM*JDM
*     write(6,*) 'ndsds = ',ndsds
*     call flush(6)
      ALLOCATE(  SDS1(NDSDS),SDSAA(NDSDS),SDSDD(NDSDS),
     +           SDS2(NDSDS),SDSBB(NDSDS),SDSAB(NDSDS),
     +          SDS1M(NDSDS),
     +          SDS2M(NDSDS) )
C
      DO I= 1,NDSDS
        SDSDD(I) = 0.0
        SDS1M(I) = 0.0
        SDS2M(I) = 0.0
      ENDDO
C
C     OPEN FILES.
C
      NPAD = 4096 - MOD(IDM*JDM,4096)
      IF     (NPAD.EQ.4096) THEN
        NPAD = 0
      ENDIF
      IF     (NPAD.EQ.0) THEN
        INQUIRE( IOLENGTH=NRECL) SDS1
      ELSE
        INQUIRE( IOLENGTH=NRECL) SDS1,PAD(1:NPAD) 
      ENDIF
C
      OPEN(UNIT=11, FILE=CHDFI1, FORM='UNFORMATTED', STATUS='OLD',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CHDFI1)
        CALL EXIT(3)
      ENDIF
      OPEN(UNIT=12, FILE=CHDFI2, FORM='UNFORMATTED', STATUS='OLD',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CHDFI2)
        CALL EXIT(3)
      ENDIF
      OPEN(UNIT=21, FILE=CHDFO,  FORM='UNFORMATTED', STATUS='NEW',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CHDFO)
        CALL EXIT(3)
      ENDIF
      OPEN(UNIT=22, FILE=CHDFB,  FORM='FORMATTED', STATUS='NEW',
     +         IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CHDFB)
        CALL EXIT(3)
      ENDIF
C
C     HEADER.
C
      WRITE(22,'(a)')      'hycom_skill:'
      WRITE(22,'(5x,a,a)')       'chdfi1 = ',trim(chdfi1)
      WRITE(22,'(5x,a,a)')       'chdfi2 = ',trim(chdfi2)
      WRITE(22,'(5x,a,a)')       'chdfo  = ',trim(chdfo)
      WRITE(22,'(5x,a,a)')       'clabel = ',trim(clabel)
      WRITE(22,'(5x,a,2i5)')     'size   = ',idm,jdm
      IF     (RMSE.NE.0.0) THEN
        WRITE(22,'(5x,a,1pg16.8)') 'rmse   = ',rmse
      ELSE
        WRITE(22,'(5x,a,1pg16.8)') 'noise  = ',noise
      ENDIF
      WRITE(22,'(5x,a,f5.1)')    'power  = ',power
      WRITE(22,'(5x,a,3i5)')     'sds1   = ',itlsds1,lstsds1,incsds1
      WRITE(22,'(5x,a,2i5)')     'sds2   = ',itlsds2,incsds2
C
C     MEAN LOOP.
C
      N    = 0
      JSDS = ITLSDS2 - INCSDS2
      DO 110 ISDS= ITLSDS1,LSTSDS1,INCSDS1
        N = N + 1
        WRITE(6,6200) ISDS
*       CALL FLUSH(6)
C
        IF     (ISDS.EQ.ITLSDS1 .OR. INCSDS2.NE.0) THEN
C
C         READ IN SDS FROM 2ND INPUT FILE, IF REQUIRED.
C
          JSDS = JSDS + INCSDS2
          READ(12,REC=JSDS,IOSTAT=IOS) SDS2(1:NDSDS)
#ifdef ENDIAN_IO
          CALL ENDIAN_SWAP(SDS2,NDSDS)
#endif
          IF     (IOS.NE.0) THEN
            write(6,*) 'Error: can''t read ',TRIM(CHDFI2)
            write(6,*) '       record number: ',JSDS
            CALL EXIT(1)
          ENDIF
        ENDIF
C
        CALL SDSPOW(SDS2,NDSDS, POWER)
        SMIN =  SPVAL
        SMAX = -SPVAL
        DO I= 1,NDSDS
          IF     (SDS2(I).NE.SPVAL) THEN
            SMIN = MIN( SDS2(I), SMIN )
            SMAX = MAX( SDS2(I), SMAX )
          ENDIF
        ENDDO
        WRITE(6,'(a,1p2g16.8)') 'SDS2: min,max =',SMIN,SMAX
*       CALL FLUSH(6)
C
C       READ IN SDS FROM IST INPUT FILE.
C
        READ(11,REC=ISDS,IOSTAT=IOS) SDS1(1:NDSDS)
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(SDS1,NDSDS)
#endif
        IF     (IOS.NE.0) THEN
          write(6,*) 'Error: can''t read ',TRIM(CHDFI1)
          write(6,*) '       record number: ',ISDS
          CALL EXIT(1)
        ENDIF
C
        CALL SDSPOW(SDS1,NDSDS, POWER)
        SMIN =  SPVAL
        SMAX = -SPVAL
        DO I= 1,NDSDS
          IF     (SDS1(I).NE.SPVAL) THEN
            SMIN = MIN( SDS1(I), SMIN )
            SMAX = MAX( SDS1(I), SMAX )
          ENDIF
        ENDDO
        WRITE(6,'(a,1p2g16.8)') 'SDS1: min,max =',SMIN,SMAX
*       CALL FLUSH(6)
C
C       PERFORM SUMS.
C
        DO I= 1,NDSDS
          IF     (SDS1M(I).EQ.SPVAL) THEN
C           do nothing
          ELSEIF (SDS1( I).NE.SPVAL .AND.
     +            SDS2( I).NE.SPVAL      ) THEN
            SDS1M(I) = SDS1M(I) +  SDS1(I)
            SDS2M(I) = SDS2M(I) +  SDS2(I)
            SDSDD(I) = SDSDD(I) + (SDS1(I) - SDS2(I))**2
          ELSE
            SDS1M(I) = SPVAL
            SDS2M(I) = SPVAL
            SDSDD(I) = SPVAL
          ENDIF
        ENDDO
  110 CONTINUE
C
C     FORM MEAN AND RMS ERROR.
C
      RN = 1.0/N
      DO I= 1,NDSDS
        IF     (SDS1M(I).NE.SPVAL) THEN
          SDS1M(I) =      SDS1M(I)*RN
          SDS2M(I) =      SDS2M(I)*RN
          SDSBB(I) = SDS1M(I) - SDS2M(I)
          SDSDD(I) = SQRT(SDSDD(I)*RN)
        ELSE
          SDSBB(I) = SPVAL
        ENDIF
      ENDDO
C
C     OUTPUT MEAN ERROR.
C
      SMIN =  SPVAL
      SMAX = -SPVAL
      DO I= 1,NDSDS
        IF     (SDSBB(I).NE.SPVAL) THEN
          SMIN = MIN( SDSBB(I), SMIN )
          SMAX = MAX( SDSBB(I), SMAX )
        ENDIF
      ENDDO
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(SDSBB,NDSDS)
#endif
      WRITE(21,REC=1,IOSTAT=IOS) SDSBB(1:NDSDS)
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(SDSBB,NDSDS)
#endif
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 
        WRITE(6,*) 'ERROR - WRITE FAILED (RECORD 1).'
        WRITE(6,*) 
        CALL EXIT(1)
      ENDIF
      WRITE(22,'(a,1p2g16.8)')
     +   '           Mean Error: min,max =',SMIN,SMAX
      WRITE( 6,'(a,1p2g16.8)')
     +   '           Mean Error: min,max =',SMIN,SMAX
*     CALL FLUSH(6)
C
C     OUTPUT RMS ERROR.
C
      SMIN =  SPVAL
      SMAX = -SPVAL
      DO I= 1,NDSDS
        IF     (SDSDD(I).NE.SPVAL) THEN
          SMIN = MIN( SDSDD(I), SMIN )
          SMAX = MAX( SDSDD(I), SMAX )
        ENDIF
      ENDDO
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(SDSDD,NDSDS)
#endif
      WRITE(21,REC=2,IOSTAT=IOS) SDSDD(1:NDSDS)
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(SDSDD,NDSDS)
#endif
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 
        WRITE(6,*) 'ERROR - WRITE FAILED (RECORD 1).'
        WRITE(6,*) 
        CALL EXIT(1)
      ENDIF
      WRITE(22,'(a,1p2g16.8)')
     +   '            RMS Error: min,max =',SMIN,SMAX
      WRITE( 6,'(a,1p2g16.8)')
     +   '            RMS Error: min,max =',SMIN,SMAX
*     CALL FLUSH(6)
C
C     FORM NORMALIZED MEAN AND RMS ERROR.
C
      IF     (RMSE.GT.0.0) THEN
        QRMSE = 1.0/RMSE
      ELSE
        QRMSE = 1.0
      ENDIF
      DO I= 1,NDSDS
        IF     (SDSBB(I).NE.SPVAL) THEN
          SDSBB(I) = ABS(SDSBB(I))*QRMSE
          SDSDD(I) =     SDSDD(I) *QRMSE
        ENDIF
      ENDDO
C
C     OUTPUT NORMALIZED MEAN ERROR.
C
      SMIN =  SPVAL
      SMAX = -SPVAL
      DO I= 1,NDSDS
        IF     (SDSBB(I).NE.SPVAL) THEN
          SMIN = MIN( SDSBB(I), SMIN )
          SMAX = MAX( SDSBB(I), SMAX )
        ENDIF
      ENDDO
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(SDSBB,NDSDS)
#endif
      WRITE(21,REC=3,IOSTAT=IOS) SDSBB(1:NDSDS)
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 
        WRITE(6,*) 'ERROR - WRITE FAILED (RECORD 1).'
        WRITE(6,*) 
        CALL EXIT(1)
      ENDIF
      WRITE(22,'(a,1p2g16.8)')
     +   'Normalized Mean Error: min,max =',SMIN,SMAX
      WRITE( 6,'(a,1p2g16.8)') 
     +   'Normalized Mean Error: min,max =',SMIN,SMAX
*     CALL FLUSH(6)
C
C     OUTPUT RMS ERROR.
C
      SMIN =  SPVAL
      SMAX = -SPVAL
      DO I= 1,NDSDS
        IF     (SDSDD(I).NE.SPVAL) THEN
          SMIN = MIN( SDSDD(I), SMIN )
          SMAX = MAX( SDSDD(I), SMAX )
        ENDIF
      ENDDO
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(SDSDD,NDSDS)
#endif
      WRITE(21,REC=4,IOSTAT=IOS) SDSDD(1:NDSDS)
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 
        WRITE(6,*) 'ERROR - WRITE FAILED (RECORD 1).'
        WRITE(6,*) 
        CALL EXIT(1)
      ENDIF
      WRITE(22,'(a,1p2g16.8)') 
     +   'Normalized RMS  Error: min,max =',SMIN,SMAX
      WRITE( 6,'(a,1p2g16.8)') 
     +   'Normalized RMS  Error: min,max =',SMIN,SMAX
*     CALL FLUSH(6)
C
C     STANDARD DEVIATION LOOP.
C
      DO I= 1,NDSDS
        IF     (SDS1M(I).NE.SPVAL) THEN
          SDSAA(I) = 0.0
          SDSBB(I) = 0.0
          SDSAB(I) = 0.0
        ELSE
          SDSAA(I) = SPVAL
          SDSBB(I) = SPVAL
          SDSAB(I) = SPVAL
        ENDIF
      ENDDO
C
      N    = 0
      JSDS = ITLSDS2 - INCSDS2
      DO 210 ISDS= ITLSDS1,LSTSDS1,INCSDS1
        N = N + 1
        IF     (ISDS.EQ.ITLSDS1) THEN
          WRITE(6,6000) ISDS
        ELSE
          WRITE(6,6100) ISDS
        ENDIF
*       CALL FLUSH(6)
C
        IF     (ISDS.EQ.ITLSDS1 .OR. INCSDS2.NE.0) THEN
C
C         READ IN SDS FROM 2ND INPUT FILE, IF REQUIRED.
C
          JSDS = JSDS + INCSDS2
          READ(12,REC=JSDS,IOSTAT=IOS) SDS2(1:NDSDS)
#ifdef ENDIAN_IO
          CALL ENDIAN_SWAP(SDS2,NDSDS)
#endif
          IF     (IOS.NE.0) THEN
            write(6,*) 'Error: can''t read ',TRIM(CHDFI2)
            write(6,*) '       record number: ',JSDS
            CALL EXIT(1)
          ENDIF
        ENDIF
C
        CALL SDSPOW(SDS2,NDSDS, POWER)
C
C       READ IN SDS FROM IST INPUT FILE.
C
        READ(11,REC=ISDS,IOSTAT=IOS) SDS1(1:NDSDS)
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(SDS1,NDSDS)
#endif
        IF     (IOS.NE.0) THEN
          write(6,*) 'Error: can''t read ',TRIM(CHDFI1)
          write(6,*) '       record number: ',ISDS
          CALL EXIT(1)
        ENDIF
C
        CALL SDSPOW(SDS1,NDSDS, POWER)
C
C       PERFORM SUMS.
C
        DO I= 1,NDSDS
          IF     (SDS1M(I).NE.SPVAL) THEN
            SDS1A = SDS1(I) - SDS1M(I)
            SDS2A = SDS2(I) - SDS2M(I)
C
            SDSAA(I) = SDSAA(I) + SDS1A**2
            SDSBB(I) = SDSBB(I) + SDS2A**2
            SDSAB(I) = SDSAB(I) + SDS1A*SDS2A
          ENDIF
        ENDDO
  210 CONTINUE
C
C     FORM ALL REQUIRED FIELDS
C
      RN = 1.0/N
      DO I= 1,NDSDS
        IF     (SDS1M(I).NE.SPVAL) THEN
          SDSAA(I) = SQRT(SDSAA(I)*RN)
          SDSBB(I) = SQRT(SDSBB(I)*RN)
          SDSAB(I) =      SDSAB(I)*RN
          IF     ((SDSAA(I)*SDSBB(I)).NE.0.0) THEN
            SDSDD(I) =  SDSAB(I)/(SDSAA(I)*SDSBB(I))
          ELSE
C
C           CORRELATION NOT DEFINABLE, RETURN 0
C
            SDSDD(I) = 0.0
          ENDIF
C
C         NOISE SETS LEVEL OF INSIGNIFICANT DIFFERENCES, AND 
C         ALLOWS A NON-ZERO RESULT EVEN WHEN SDSAA*SDSBB IS ZERO.
C
          IF     (NOISE.GT.0.0) THEN
            SDSAAX   = SQRT(SDSAA(I)**2 + NOISE**2)
            SDSBBX   = SQRT(SDSBB(I)**2 + NOISE**2)
            SDSDDX   = (SDSDD(I)*SDSAB(I) + NOISE**2) / (SDSAAX*SDSBBX)
          ELSE
            SDSAAX   =      SDSAA(I)
            SDSBBX   = MAX( SDSBB(I), RMSE )
            SDSDDX   =      SDSDD(I)
          ENDIF
          SDS1M(I) = (    (SDS1M(I)  - SDS2M(I))/SDSBBX)**2
          SDS2M(I) = ( ABS(SDSDDX)   -    SDSAAX/SDSBBX)**2
          SDS1M(I) = MIN( SDS1M(I), 1.0 )  !clip at 1
          SDS2M(I) = MIN( SDS2M(I), 1.0 )  !clip at 1
          SDSAA(I) = SDSDDX**2 - SDS1M(I) - SDS2M(I)  !minimum -2
        ENDIF
      ENDDO
      WRITE(6,*)
C
C     OUTPUT CONDITIONAL BIAS.
C
      SMIN =  SPVAL
      SMAX = -SPVAL
      DO I= 1,NDSDS
        IF     (SDS2M(I).NE.SPVAL) THEN
          SMIN = MIN( SDS2M(I), SMIN )
          SMAX = MAX( SDS2M(I), SMAX )
        ENDIF
      ENDDO
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(SDS2M,NDSDS)
#endif
      WRITE(21,REC=5,IOSTAT=IOS) SDS2M(1:NDSDS)
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 
        WRITE(6,*) 'ERROR - WRITE FAILED (RECORD 1).'
        WRITE(6,*) 
        CALL EXIT(1)
      ENDIF
      WRITE(22,'(a,1p2g16.8)')
     +   '     Conditional Bias: min,max =',SMIN,SMAX
      WRITE( 6,'(a,1p2g16.8)')
     +   '     Conditional Bias: min,max =',SMIN,SMAX
*     CALL FLUSH(6)
C
C     OUTPUT UNCONDITIONAL BIAS.
C
      SMIN =  SPVAL
      SMAX = -SPVAL
      DO I= 1,NDSDS
        IF     (SDS1M(I).NE.SPVAL) THEN
          SMIN = MIN( SDS1M(I), SMIN )
          SMAX = MAX( SDS1M(I), SMAX )
        ENDIF
      ENDDO
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(SDS1M,NDSDS)
#endif
      WRITE(21,REC=6,IOSTAT=IOS) SDS1M(1:NDSDS)
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 
        WRITE(6,*) 'ERROR - WRITE FAILED (RECORD 1).'
        WRITE(6,*) 
        CALL EXIT(1)
      ENDIF
      WRITE(22,'(a,1p2g16.8)')
     +   '   Unconditional Bias: min,max =',SMIN,SMAX
      WRITE( 6,'(a,1p2g16.8)')
     +   '   Unconditional Bias: min,max =',SMIN,SMAX
*     CALL FLUSH(6)
C
C     OUTPUT SKILL SCORE BIAS.
C
      SMIN =  SPVAL
      SMAX = -SPVAL
      DO I= 1,NDSDS
        IF     (SDSAA(I).NE.SPVAL) THEN
          SMIN = MIN( SDSAA(I), SMIN )
          SMAX = MAX( SDSAA(I), SMAX )
        ENDIF
      ENDDO
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(SDSAA,NDSDS)
#endif
      WRITE(21,REC=7,IOSTAT=IOS) SDSAA(1:NDSDS)
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 
        WRITE(6,*) 'ERROR - WRITE FAILED (RECORD 1).'
        WRITE(6,*) 
        CALL EXIT(1)
      ENDIF
      WRITE(22,'(a,1p2g16.8)')
     +   '     Skill Score Bias: min,max =',SMIN,SMAX
      WRITE( 6,'(a,1p2g16.8)')
     +   '     Skill Score Bias: min,max =',SMIN,SMAX
*     CALL FLUSH(6)
C
C     OUTPUT R**2 BIAS.
C
      SMIN =  SPVAL
      SMAX = -SPVAL
      DO I= 1,NDSDS
        IF     (SDSDD(I).NE.SPVAL) THEN
          SMIN = MIN( SDSDD(I), SMIN )
          SMAX = MAX( SDSDD(I), SMAX )
        ENDIF
      ENDDO
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(SDSDD,NDSDS)
#endif
      WRITE(21,REC=8,IOSTAT=IOS) SDSDD(1:NDSDS)
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 
        WRITE(6,*) 'ERROR - WRITE FAILED (RECORD 1).'
        WRITE(6,*) 
        CALL EXIT(1)
      ENDIF
      WRITE(22,'(a,1p2g16.8)')
     +   '   Correlation Coeff.: min,max =',SMIN,SMAX
      WRITE( 6,'(a,1p2g16.8)')
     +   '   Correlation Coeff.: min,max =',SMIN,SMAX
*     CALL FLUSH(6)
C
C     OUTPUT STANDARD DEVIATION OF 'OBSERVATIONS'.
C
      SMIN =  SPVAL
      SMAX = -SPVAL
      DO I= 1,NDSDS
        IF     (SDSBB(I).NE.SPVAL) THEN
          SMIN = MIN( SDSBB(I), SMIN )
          SMAX = MAX( SDSBB(I), SMAX )
        ENDIF
      ENDDO
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(SDSBB,NDSDS)
#endif
      WRITE(21,REC=9,IOSTAT=IOS) SDSBB(1:NDSDS)
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 
        WRITE(6,*) 'ERROR - WRITE FAILED (RECORD 1).'
        WRITE(6,*) 
        CALL EXIT(1)
      ENDIF
      WRITE(22,'(a,1p2g16.8)')
     +   '   Standard Deviation: min,max =',SMIN,SMAX
      WRITE( 6,'(a,1p2g16.8)')
     +   '   Standard Deviation: min,max =',SMIN,SMAX
*     CALL FLUSH(6)
C
      CLOSE(11)
      CLOSE(12)
      CLOSE(21)
      CLOSE(22)
C
      CALL EXIT(0)
C
 6000 FORMAT(' Started ISDS = ',I5,$)
 6100 FORMAT(                   I5,$)
 6200 FORMAT(' Started ISDS = ',I5)
      END
      SUBROUTINE SDSPOW(SDS,NDSDS, POWER)
      IMPLICIT NONE
C
      INTEGER NDSDS
      REAL    SDS(NDSDS), POWER
C
C     RAISE THE INPUT TO A POWER.
C
C     POWER - RAISE INPUT FIELDS TO THIS POWER
C               =-1.0 TAKE -1/     OF THE INPUT
C               =-0.5 TAKE -1/SQRT OF THE INPUT
C               = 0.0 TAKE LOG  OF THE INPUT
C               = 0.5 TAKE SQRT OF THE INPUT
C               = 1.0 USE THE INPUT AS IS (DEFAULT)
C
      REAL       SPVAL
      PARAMETER (SPVAL=2.0**100)
C
      INTEGER I
      REAL    SMAX,SMIN
C
      IF     (POWER.NE.1.0) THEN
        SMIN =  SPVAL
        SMAX = -SPVAL
        DO I= 1,NDSDS
          IF     (SDS(I).NE.SPVAL) THEN
            SMIN = MIN( SDS(I), SMIN )
            SMAX = MAX( SDS(I), SMAX )
          ENDIF
        ENDDO
        IF     (SMAX.LT.0.0) THEN
          DO I= 1,NDSDS
            SDS(I) = - SDS(I)
          ENDDO
          SMIN = -SMAX
          SMAX = -SMIN
        ENDIF
        IF     (SMIN.LE.0.0) THEN
          WRITE(6,*) 
          WRITE(6,*) 'ERROR IN SDSPOW - SDS NOT ONE SIGN'
          WRITE(6,*) 'SMIN,SMAX = ',SMIN,SMAX
          WRITE(6,*)
          STOP
        ENDIF
        IF     (POWER.EQ.0.0) THEN
          DO I= 1,NDSDS
            IF     (SDS(I).NE.SPVAL) THEN
              SDS(I) = LOG(SDS(I))
            ENDIF
          ENDDO
        ELSEIF (POWER.EQ.0.5) THEN
          DO I= 1,NDSDS
            IF     (SDS(I).NE.SPVAL) THEN
              SDS(I) = SQRT(SDS(I))
            ENDIF
          ENDDO
        ELSEIF (POWER.GT.0.0) THEN
          DO I= 1,NDSDS
            IF     (SDS(I).NE.SPVAL) THEN
              SDS(I) = SDS(I)**POWER
            ENDIF
          ENDDO
        ELSEIF (POWER.EQ.-1.0) THEN
          DO I= 1,NDSDS
            IF     (SDS(I).NE.SPVAL) THEN
              SDS(I) = -1.0/SDS(I)
            ENDIF
          ENDDO
        ELSEIF (POWER.EQ.-0.5) THEN
          DO I= 1,NDSDS
            IF     (SDS(I).NE.SPVAL) THEN
              SDS(I) = -SQRT(SDS(I))
            ENDIF
          ENDDO
        ELSE
          DO I= 1,NDSDS
            IF     (SDS(I).NE.SPVAL) THEN
              SDS(I) = -(SDS(I)**POWER)
            ENDIF
          ENDDO
        ENDIF
      ENDIF
      RETURN
C     END OF SDSPOW
      END
      SUBROUTINE ARGS(CHDFI1,CHDFI2,CHDFO,CHDFB,
     +                IDM,JDM, CLABEL, RMSE, NOISE, POWER,
     +                ITLSDS1,LSTSDS1,INCSDS1,
     +                ITLSDS2,        INCSDS2)
C
      IMPLICIT NONE
      CHARACTER*256 CHDFI1,CHDFI2,CHDFO,CHDFB, CLABEL
      INTEGER       IDM,JDM,ITLSDS1,LSTSDS1,INCSDS1,ITLSDS2,INCSDS2
      REAL          RMSE,NOISE,POWER
C
C     COMMAND LINE ARGUMENT PROCESSING.
C
      INTEGER          MAXF,  MAXC
      PARAMETER       (MAXF=4,MAXC=512)
C
      CHARACTER*512    CARG
      INTEGER          NARGC, I, II, LSTSDS2, NB, NVALS
      REAL             RVAL(MAXF)
      INTEGER          IVAL(MAXF)
      CHARACTER*(MAXC) CVAL(MAXF)
C
      INTEGER  IARGC
      EXTERNAL IARGC
C
      NARGC = IARGC()
      IF     (NARGC.LT.1) THEN
        CALL USAGE
      ENDIF
C
      CHDFI1  =  ' '
      CHDFI2  =  ' '
      CHDFO   =  ' '
      CLABEL  =  ' '
      RMSE    =  0.0
      NOISE   =  0.0
      POWER   =  1.0
      IDM     =  0
      JDM     =  0
      ITLSDS1 =  0
      LSTSDS1 =  0
      INCSDS1 =  0
      ITLSDS2 =  0
      LSTSDS2 =  0
      INCSDS2 =  0
C
      I = 0
      DO 110 II= 1,99
        I = I + 1
        IF     (I.GT.NARGC) THEN
          GOTO 1110
        ENDIF
C
        CALL OTKARG(I,CARG)
        NB = LEN_TRIM(CARG)
C
        IF     (CARG.EQ.'-i1') THEN
C
C         1ST INPUT FILENAME.
C
          I = I + 1
          IF     (I.GT.NARGC) THEN
            WRITE(6,*) 
            WRITE(6,*) 'Option -i1 should be followed by the name of'
            WRITE(6,*) ' existing HDF file.'
            WRITE(6,*) 
            CALL USAGE
          ELSEIF (CHDFI1.NE.' ') THEN
            WRITE(6,*) 
            WRITE(6,*) 'Please specify option -i1 exactly once'
            WRITE(6,*) 
            CALL USAGE
          ENDIF
          CALL OTKARG(I,CARG)
          NB    = LEN_TRIM(CARG)
          CHDFI1 = CARG(1:NB)
        ELSEIF (CARG.EQ.'-i2') THEN
C
C         2ND INPUT FILENAME.
C
          I = I + 1
          IF     (I.GT.NARGC) THEN
            WRITE(6,*) 
            WRITE(6,*) 'Option -i2 should be followed by the name of'
            WRITE(6,*) ' existing HDF file.'
            WRITE(6,*) 
            CALL USAGE
          ELSEIF (CHDFI2.NE.' ') THEN
            WRITE(6,*) 
            WRITE(6,*) 'Please specify option -i2 exactly once'
            WRITE(6,*) 
            CALL USAGE
          ENDIF
          CALL OTKARG(I,CARG)
          NB    = LEN_TRIM(CARG)
          CHDFI2 = CARG(1:NB)
        ELSEIF (CARG.EQ.'-o') THEN
C
C         OUTPUT FILENAME.
C
          I = I + 1
          IF     (I.GT.NARGC) THEN
            WRITE(6,*) 
            WRITE(6,*) 'Option -o should be followed by filename'
            WRITE(6,*) 
            CALL USAGE
          ELSEIF (CHDFO.NE.' ') THEN
            WRITE(6,*) 
            WRITE(6,*) 'Please specify option -o exactly once'
            WRITE(6,*) 
            CALL USAGE
          ENDIF
          CALL OTKARG(I,CARG)
          NB    = LEN_TRIM(CARG)
          CHDFO = CARG(1:NB)
        ELSEIF (CARG.EQ.'-l') THEN
C
C         TITLE.
C
          I = I + 1
          IF     (I.GT.NARGC) THEN
            WRITE(6,*) 
            WRITE(6,*) 'Option -l should be followed by label string'
            WRITE(6,*) 
            CALL USAGE
          ELSEIF (CLABEL.NE.' ') THEN
            WRITE(6,*) 
            WRITE(6,*) 'Please specify option -l exactly once'
            WRITE(6,*) 
            CALL USAGE
          ENDIF
          CALL OTKARG(I,CARG)
          NB     = LEN_TRIM(CARG)
          CLABEL = CARG(1:NB)
        ELSEIF (CARG.EQ.'-size') THEN
C
C         ARRAY SIZE.
C
          I = I + 1
          CALL OTKPARG(I,MAXF,NVALS,RVAL,IVAL,CVAL,MAXC)
          IF     (NVALS.NE.2) THEN
            WRITE(6,*) 
            WRITE(6,*) 'Usage error for option -size'
            WRITE(6,*) NVALS,' values were specified, '
            WRITE(6,*) 'where 2 were expected'
            WRITE(6,*) 
            CALL USAGE 
          ELSEIF (IDM.NE.0) THEN
            WRITE(6,*) 
            WRITE(6,*) 'Please specify option -size exactly once'
            WRITE(6,*) 
            CALL USAGE
          ENDIF
          IDM = IVAL(1) 
          JDM = IVAL(2) 
        ELSEIF (CARG.EQ.'-sds1') THEN
C
C         WHICH SDS'S.
C
          I = I + 1
          CALL OTKPARG(I,MAXF,NVALS,RVAL,IVAL,CVAL,MAXC)
          IF     (NVALS.NE.3) THEN
            WRITE(6,*) 
            WRITE(6,*) 'Usage error for option -sds1'
            WRITE(6,*) NVALS,' values were specified, '
            WRITE(6,*) 'where 3 were expected'
            WRITE(6,*) 
            CALL USAGE 
          ELSEIF (ITLSDS1.NE.0) THEN
            WRITE(6,*) 
            WRITE(6,*) 'Please specify option -sds1 exactly once'
            WRITE(6,*) 
            CALL USAGE
          ENDIF
          ITLSDS1 = IVAL(1) 
          LSTSDS1 = IVAL(2) 
          INCSDS1 = IVAL(3) 
          IF     ( INCSDS1.EQ.      0        .OR.
     +             ITLSDS1.LT.      1        .OR.
     +             LSTSDS1.LT.      1        .OR.
     +            (INCSDS1.GT.      0 .AND.
     +             ITLSDS1.GT.LSTSDS1      ) .OR.
     +            (INCSDS1.LT.      0 .AND.
     +             ITLSDS1.LT.LSTSDS1      )     ) THEN
            WRITE(6,*) 
            WRITE(6,*) 'error - bad itlsds1,lstsds1,incsds1'
            WRITE(6,*) 
            CALL USAGE
          ENDIF
        ELSEIF (CARG.EQ.'-sds2') THEN
C
C         WHICH SDS'S.
C
          I = I + 1
          CALL OTKPARG(I,MAXF,NVALS,RVAL,IVAL,CVAL,MAXC)
          IF     (NVALS.NE.2) THEN
            WRITE(6,*) 
            WRITE(6,*) 'Usage error for option -sds2'
            WRITE(6,*) NVALS,' values were specified, '
            WRITE(6,*) 'where 2 were expected'
            WRITE(6,*) 
            CALL USAGE 
          ELSEIF (ITLSDS2.NE.0) THEN
            WRITE(6,*) 
            WRITE(6,*) 'Please specify option -sds2 exactly once'
            WRITE(6,*) 
            CALL USAGE
          ELSEIF (ITLSDS1.EQ.0) THEN
            WRITE(6,*) 
            WRITE(6,*) 'Must have -sds1 before -sds2'
            WRITE(6,*) 
            CALL USAGE
          ENDIF
          ITLSDS2 = IVAL(1) 
          INCSDS2 = IVAL(2) 
          LSTSDS2 = ITLSDS2 + INCSDS2*((LSTSDS1 - ITLSDS1)/INCSDS1)
          IF     ( ITLSDS2.LT.      1        .OR.
     +             LSTSDS2.LT.      1        .OR.
     +            (INCSDS2.GT.      0 .AND.
     +             ITLSDS2.GT.LSTSDS2      ) .OR.
     +            (INCSDS2.LT.      0 .AND.
     +             ITLSDS2.LT.LSTSDS2      )     ) THEN
            WRITE(6,*) 
            WRITE(6,*) 'error - bad itlsds2,incsds2'
            WRITE(6,*) 
            CALL USAGE
          ENDIF
        ELSEIF (CARG.EQ.'-rmse') THEN
C
C         STANDARD RMSE VALUE.
C
          I = I + 1
          CALL OTKPARG(I,MAXF,NVALS,RVAL,IVAL,CVAL,MAXC)
          IF     (NVALS.NE.1) THEN
            WRITE(6,*) 
            WRITE(6,*) 'Usage error for option -rmse'
            WRITE(6,*) NVALS,' values were specified, '
            WRITE(6,*) 'where 1 was expected'
            WRITE(6,*) 
            CALL USAGE 
          ELSEIF (RMSE.NE.0.0) THEN
            WRITE(6,*) 
            WRITE(6,*) 'Please specify option -rmse exactly once'
            WRITE(6,*) 
            CALL USAGE
          ELSEIF (NOISE.NE.0.0) THEN
            WRITE(6,*) 
            WRITE(6,*) 'Please specify either -noise or -rmse'
            WRITE(6,*) 
            CALL USAGE
          ENDIF
          RMSE = RVAL(1) 
          IF     (RMSE.LE.0.0) THEN
            WRITE(6,*) 
            WRITE(6,*) 'error - rmse must be > 0'
            WRITE(6,*) 
            CALL USAGE
          ENDIF
        ELSEIF (CARG.EQ.'-noise') THEN
C
C         STANDARD NOISE VALUE.
C
          I = I + 1
          CALL OTKPARG(I,MAXF,NVALS,RVAL,IVAL,CVAL,MAXC)
          IF     (NVALS.NE.1) THEN
            WRITE(6,*) 
            WRITE(6,*) 'Usage error for option -noise'
            WRITE(6,*) NVALS,' values were specified, '
            WRITE(6,*) 'where 1 was expected'
            WRITE(6,*) 
            CALL USAGE 
          ELSEIF (NOISE.NE.0.0) THEN
            WRITE(6,*) 
            WRITE(6,*) 'Please specify option -noise at most once'
            WRITE(6,*) 
            CALL USAGE
          ELSEIF (RMSE .NE.0.0) THEN
            WRITE(6,*) 
            WRITE(6,*) 'Please specify either -noise or -rmse'
            WRITE(6,*) 
            CALL USAGE
          ENDIF
          NOISE = RVAL(1) 
          IF     (NOISE.LE.0.0) THEN
            WRITE(6,*) 
            WRITE(6,*) 'error - noise must be > 0'
            WRITE(6,*) 
            CALL USAGE
          ENDIF
        ELSEIF (CARG.EQ.'-power') THEN
C
C         POWER VALUE.
C
          I = I + 1
          CALL OTKPARG(I,MAXF,NVALS,RVAL,IVAL,CVAL,MAXC)
          IF     (NVALS.NE.1) THEN
            WRITE(6,*) 
            WRITE(6,*) 'Usage error for option -power'
            WRITE(6,*) NVALS,' values were specified, '
            WRITE(6,*) 'where 1 was expected'
            WRITE(6,*) 
            CALL USAGE 
          ELSEIF (POWER.NE.1.0) THEN
            WRITE(6,*) 
            WRITE(6,*) 'Please specify option -power at most once'
            WRITE(6,*) 
            CALL USAGE
          ENDIF
          POWER = RVAL(1) 
        ELSE
          WRITE(6,*) 
          WRITE(6,*) 'Bad command line option'
          WRITE(6,*) 
          CALL USAGE
        ENDIF
  110 CONTINUE
 1110 CONTINUE
C
C     PRINTOUT ARGUMENTS.
C
      WRITE(6,*)
      WRITE(6,*) 'chdfi1 = ',trim(chdfi1)
      WRITE(6,*) 'chdfi2 = ',trim(chdfi2)
      WRITE(6,*) 'chdfo  = ',trim(chdfo)
      WRITE(6,*) 'size   = ',idm,jdm
      WRITE(6,*) 'clabel = ',trim(clabel)
      WRITE(6,*) 'rmse   = ',rmse
      WRITE(6,*) 'noise  = ',noise
      WRITE(6,*) 'power  = ',power
      WRITE(6,*) 'sds1   = ',itlsds1,lstsds1,incsds1
      WRITE(6,*) 'sds2   = ',itlsds2,incsds2
      WRITE(6,*)
C
C     TEST FOR NON-OPTIONAL COMMAND LINE ARGUMENTS.
C
      IF     (CHDFI1   .EQ. ' ') THEN
        WRITE(6,*) 
        WRITE(6,*) 'error - incomplete command line (chdfi1).'
        WRITE(6,*) 
        CALL USAGE 
      ENDIF
      IF     (CHDFI2   .EQ. ' ') THEN
        WRITE(6,*) 
        WRITE(6,*) 'error - incomplete command line (chdfi2).'
        WRITE(6,*) 
        CALL USAGE 
      ENDIF
      IF     (CHDFO    .EQ. ' ') THEN
        WRITE(6,*) 
        WRITE(6,*) 'error - incomplete command line (chdfo).'
        WRITE(6,*) 
        CALL USAGE 
      ENDIF
      IF     (IDM      .EQ. 0) THEN
        WRITE(6,*) 
        WRITE(6,*) 'error - incomplete command line (size).'
        WRITE(6,*) 
        CALL USAGE 
      ENDIF
      IF     (RMSE     .EQ.0.0 .AND.
     +        NOISE    .EQ.0.0      ) THEN
        WRITE(6,*) 
        WRITE(6,*) 'error - incomplete command line (rmse or noise).'
        WRITE(6,*) 
        CALL USAGE 
      ENDIF
      IF     (ITLSDS1  .EQ. 0) THEN
        WRITE(6,*) 
        WRITE(6,*) 'error - incomplete command line (sds1).'
        WRITE(6,*) 
        CALL USAGE 
      ENDIF
      IF     (ITLSDS2  .EQ. 0) THEN
        WRITE(6,*) 
        WRITE(6,*) 'error - incomplete command line (sds2).'
        WRITE(6,*) 
        CALL USAGE 
      ENDIF
C
C     FORM CHDFB FROM CHDFO.
C
      II = LEN_TRIM(CHDFO)
      IF     (CHDFO(II-1:II).EQ.'.a') THEN
        CHDFB = CHDFO(1:II-2)//'.b'
      ELSE
        WRITE(6,*) 
        WRITE(6,*) 'error - chdfo must end with .a'
        WRITE(6,*) 
        CALL USAGE 
      ENDIF
      RETURN
C     END OF ARGS.
      END
      SUBROUTINE USAGE
      WRITE(6,*) ' '
      WRITE(6,*) 'SYNOPSIS'
      WRITE(6,*) ' '
      WRITE(6,10) 
     + 'hycom_skill -i1 <chdfi1> -i2 <chdfi2> -o <chdfo> -l <clabel>'
      WRITE(6,10) 
     + '            -size "<idm>,<jdm>"'
      WRITE(6,10) 
     + '            -rmse <rmse> -noise <noise> -power <p>'
      WRITE(6,10) 
     + '            -sds1 "<itlsds1>,<lstsds1>,<incsds1>"'
      WRITE(6,10) 
     + '            -sds2 "<itlsds2>,<incsds2>"'
      WRITE(6,*) ' '
      WRITE(6,*) 'OPTIONS'
      WRITE(6,10) ' '
      WRITE(6,10) '  <chdfi1>  is the 1st input hycom .a file'
      WRITE(6,10) '  <chdfi2>  is the 2nd input hycom .a file'
      WRITE(6,10) '  <chdfo>   is the    output hycom .a file'
      WRITE(6,10) '  <clabel>  is the title'
      WRITE(6,10) '  <idm>     is the 1st array dimension'
      WRITE(6,10) '  <jdm>     is the 2nd array dimension'
      WRITE(6,10) '  <rmse>    is the target rmse value'
      WRITE(6,10) '             > 0.0 -> noise=0.0; sets floor on SS'
      WRITE(6,10) '  <noise>   is the rms of independent noise'
      WRITE(6,10) '             > 0.0 -> rmse=0.0;  modified R^2 and SS'
      WRITE(6,10) '             use exactly one of -rmse and -noise'
      WRITE(6,10) '  <p>       is the power to raise fields to'
      WRITE(6,10) '             =-1.0 -> take -1/      of the fields'
      WRITE(6,10) '             =-0.5 -> take -1/sqrt  of the fields'
      WRITE(6,10) '             = 0.0 -> take the log  of the fields'
      WRITE(6,10) '             = 0.5 -> take the sqrt of the fields'
      WRITE(6,10) '             = 1.0 -> use fields as is (default)'
      WRITE(6,10) '  <itlsds1> is the first field to use from chdfi1'
      WRITE(6,10) '  <lstsds1> is the last  field to use from chdfi1'
      WRITE(6,10) '  <incsds1> is the increment between field ''s'
      WRITE(6,10) '             >0 -> must have itlsds1<=lstsds1'
      WRITE(6,10) '             <0 -> must have itlsds1>=lstsds1'
      WRITE(6,10) '  <itlsds2> is the first field to use from chdfi2'
      WRITE(6,10) '  <incsds2> is the increment between field''s'
      WRITE(6,10) '             =0 -> always use itlsds2'
      WRITE(6,10) ' '
      WRITE(6,10) ' '    
      WRITE(6,10) 
     + '@(#)AUTHOR A. J. Wallcraft, NRL'
      WRITE(6,10) ' '
C
   10 FORMAT(5X,A)
C
*     CALL FLUSH(6)
      CALL EXIT(1)
      END
      SUBROUTINE OTKPARG(IARG,M,NVALS,RVALS,IVALS,CVALS,MC)
      IMPLICIT NONE
C
C ------ This subroutine provides a menu-independent interface
C        for parsing command line argument IARG as a sequence
C        of values.  The calling sequence is similar to that
C        of PARSE, except an integer parameter takes the place
C        of a character string and it's length (as in PARSE).
C
      INTEGER IARG,NVALS,M,IVALS(M),MC
      REAL RVALS(M)
      CHARACTER*(*) CVALS(M)
C
      INTEGER MAXCH, I, NB
      PARAMETER (MAXCH=512)
      CHARACTER*(MAXCH) STRING
C
      DOUBLE PRECISION DVALS(20)
C
      CALL OTKARG(IARG,STRING)
C
C ------ Pass only the "non-blank" length of the character
C        string to be parsed.
C
      NB = LEN_TRIM(STRING)
      CALL PARSE(STRING,NB,M,NVALS,DVALS,IVALS,CVALS,MC)
C
C ------ Move the double precision (local) values into
C        the single precision subroutine arguments. On the
C        Cray this is redundant.
C
      DO I=1,NVALS
        RVALS(I) = DVALS(I)
      ENDDO
C
      RETURN
      END
      SUBROUTINE OTKARG(I,C)
      implicit none
C
C ------ Return command line argument 'I' in the character
C        variable 'C'
C
      character*(*) C
      integer i
C
      call getarg(i,c)
C
      return
      end

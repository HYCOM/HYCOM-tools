      PROGRAM FZONAL
      IMPLICIT NONE
C
C  hycom_zonal - Usage:  hycom_zonal file.a idm jdm k [if il] [jc js jw]
C
C                 prints the min, max, mean and rms of zonal j-extents
C                 from the k-th (1:idm,1:jdm) array in file.a
C                 if - first meridional point in sum (default 1)
C                 il - last  meridional point in sum (default idm)
C                 jc - center of one zonal region    (default 1)
C                 js - spacing between zonal centers (default 1)
C                 jw - width of each zonal region    (default 1)
C
C  Use hycom_zonal_lat for true latitude extents that allow for
C  variable grid spacing.
C
C  file.a is assumed to contain idm*jdm 32-bit IEEE real values for
C   each array, in standard f77 element order, followed by padding
C   to a multiple of 4096 32-bit words, but otherwise with no control
C   bytes/words, and input values of 2.0**100 indicating a data void.
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  January 2001.
C
      REAL*4, ALLOCATABLE :: A(:,:)
      REAL*4              :: PAD(4096)
      INTEGER IOS
      INTEGER      IARGC
      INTEGER      NARG
      CHARACTER*240 CARG
C
      INTEGER      IDM,JDM,K,NPAD,IF,IL,JC,JS,JW
      CHARACTER*240 CFILE
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.4) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,5000) IDM
        CALL GETARG(3,CARG)
        READ(CARG,5000) JDM
        CALL GETARG(4,CARG)
        READ(CARG,5000) K
        IF = 1
        IL = IDM
        JC = 1
        JS = 1
        JW = 1
      ELSEIF (NARG.EQ.6) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,5000) IDM
        CALL GETARG(3,CARG)
        READ(CARG,5000) JDM
        CALL GETARG(4,CARG)
        READ(CARG,5000) K
        CALL GETARG(5,CARG)
        READ(CARG,5000) IF
        CALL GETARG(6,CARG)
        READ(CARG,5000) IL
        JC = 1
        JS = 1
        JW = 1
      ELSEIF (NARG.EQ.7) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,5000) IDM
        CALL GETARG(3,CARG)
        READ(CARG,5000) JDM
        CALL GETARG(4,CARG)
        READ(CARG,5000) K
        IF = 1
        IL = IDM
        CALL GETARG(5,CARG)
        READ(CARG,5000) JC
        CALL GETARG(6,CARG)
        READ(CARG,5000) JS
        CALL GETARG(7,CARG)
        READ(CARG,5000) JW
      ELSE
        WRITE(6,*) 
     +   'Usage:  hycom_zonal file.a idm jdm k [if il] [jc js jw]'
        CALL EXIT(1)
      ENDIF
C
      NPAD = 4096 - MOD(IDM*JDM,4096)
      IF     (NPAD.EQ.4096) THEN
        NPAD = 0
      ENDIF
C
      ALLOCATE( A(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_zonal: could not allocate ',
     +             IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
C
      CALL ZONAL(A,IDM,JDM,PAD,NPAD, K,IF,IL,JC,JS,JW, CFILE)
      CALL EXIT(0)
 5000 FORMAT(I4)
      END
      SUBROUTINE ZONAL(A,IDM,JDM, PAD,NPAD, K,IF,IL,JC,JS,JW, CFILE)
      IMPLICIT NONE
C
      REAL*4     SPVAL
      PARAMETER (SPVAL=2.0**100)
C
      CHARACTER*240 CFILE
      INTEGER      IDM,JDM,NPAD,K,IF,IL,JC,JS,JW
      REAL*4       A(IDM,JDM),PAD(NPAD)
C
C     MOST OF WORK IS DONE HERE.
C
#ifdef sun
      INTEGER      IR_ISNAN
C
#endif
      CHARACTER*18 CASN
      INTEGER      LEN_TRIM
      INTEGER      I,J,JF,JL,JZ,JZF,JZL,IOS,NRECL
      REAL*8       SUMA,SUMB,SUMR
      REAL*4       MINA,MAXA
#ifdef CRAY
      INTEGER*8    IU8,IOS8
#endif
C
      INQUIRE( IOLENGTH=NRECL) A,PAD
#ifdef CRAY
#ifdef t3e
      IF     (MOD(NRECL,4096).EQ.0) THEN
        WRITE(CASN,8000) NRECL/4096
 8000   FORMAT('-F cachea:',I4.4,':1:0')
        IU8 = 11
        CALL ASNUNIT(IU8,CASN,IOS8)
        IF     (IOS8.NE.0) THEN
          write(6,*) 'Error: can''t asnunit 11'
          write(6,*) 'ios  = ',ios8
          write(6,*) 'casn = ',casn
          CALL EXIT(5)
        ENDIF
      ENDIF
#else
      CALL ASNUNIT(11,'-F syscall -N ieee',IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t asnunit 11'
        write(6,*) 'ios = ',ios
        CALL EXIT(5)
      ENDIF
#endif
#endif
      OPEN(UNIT=11, FILE=CFILE, FORM='UNFORMATTED', STATUS='OLD',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',CFILE(1:LEN_TRIM(CFILE))
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
C
        READ(11,REC=K,IOSTAT=IOS) A
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
        IF     (IOS.NE.0) THEN
          IF     (K.EQ.1) THEN
            WRITE(6,*) 'can''t read ',CFILE(1:LEN_TRIM(CFILE))
            CALL EXIT(4)
          ENDIF
        ENDIF
C
        WRITE(6,6000) 
     +    CFILE(1:LEN_TRIM(CFILE)),
     +    IDM,JDM,K,
     +    IF,IL,
     +    JC,JS,JW
C
        JF = JC - JS*((JC-  1)/JS)
        JL = JC + JS*((JDM-JC)/JS)
        DO JZ= JF,JL,JS
          JZF = MAX( JZ-JW/2,   1)
          JZL = MIN( JZ+JW/2, JDM)
          SUMA = 0.0
          SUMB = 0.0
          MINA =  HUGE(MINA)
          MAXA = -HUGE(MAXA)
          DO J= JZF,JZL
            DO I= IF,MIN(IL,IDM)
              IF     (A(I,J).NE.SPVAL) THEN
                MINA = MIN( MINA, A(I,J) )
                MAXA = MAX( MAXA, A(I,J) )
                SUMA = SUMA + A(I,J)
                SUMB = SUMB + 1.D0
              ENDIF
            ENDDO
            DO I= IDM+1,IL
              IF     (A(I-IDM,J).NE.SPVAL) THEN
                MINA = MIN( MINA, A(I-IDM,J) )
                MAXA = MAX( MAXA, A(I-IDM,J) )
                SUMA = SUMA + A(I-IDM,J)
                SUMB = SUMB + 1.D0
              ENDIF
            ENDDO
          ENDDO
          IF     (SUMB.NE.0.0) THEN
            SUMA = SUMA/SUMB
          ENDIF
C
          SUMR = 0.0
          DO J= JZF,JZL
            DO I= IF,MIN(IL,IDM)
              IF     (A(I,J).NE.SPVAL) THEN
                SUMR = SUMR + (A(I,J)-SUMA)**2
              ENDIF
            ENDDO
            DO I= IDM+1,IL
              IF     (A(I-IDM,J).NE.SPVAL) THEN
                SUMR = SUMR + (A(I-IDM,J)-SUMA)**2
              ENDIF
            ENDDO
          ENDDO
          IF     (SUMB.NE.0.0) THEN
            SUMR = SQRT(SUMR/SUMB)
          ENDIF
          WRITE(6,'(3I5,1p4g16.8)') JZF,JZL,JZ,SUMA,SUMR,MINA,MAXA
        ENDDO
      RETURN
 6000 FORMAT(
     +   '# hycom_zonal file.a idm jdm k if il jc js jw' /
     +   '# hycom_zonal ',A   /
     +   '#             ',8I5 /
     +   '#  JF   JL   JZ     MEAN            RMS',
     +           '             MIN            MAX')
 
      END

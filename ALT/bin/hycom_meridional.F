      PROGRAM FMERID
      IMPLICIT NONE
C
C  hycom_meridional - Usage:  hycom_meridional file.a idm jdm k [jf jl] [ic is iw]
C
C                 prints the mean and rms of meridional extents from
C                 the k-th (1:idm,1:jdm) array in file.a
C                 jf - first meridional point in sum (default 1)
C                 jl - last  meridional point in sum (default idm)
C                 ic - center of one meridional region    (default 1)
C                 is - spacing between meridional centers (default 1)
C                 iw - width of each meridional region    (default 1)
C
C  Use hycom_meridional_lon for true longitude extents that allow for
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
      INTEGER      IDM,JDM,K,NPAD,JF,JL,IC,IS,IW
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
        JF = 1
        JL = JDM
        IC = 1
        IS = 1
        IW = 1
      ELSEIF (NARG.EQ.6) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,5000) IDM
        CALL GETARG(3,CARG)
        READ(CARG,5000) JDM
        CALL GETARG(4,CARG)
        READ(CARG,5000) K
        CALL GETARG(5,CARG)
        READ(CARG,5000) JF
        CALL GETARG(6,CARG)
        READ(CARG,5000) JL
        IC = 1
        IS = 1
        IW = 1
      ELSEIF (NARG.EQ.7) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,5000) IDM
        CALL GETARG(3,CARG)
        READ(CARG,5000) JDM
        CALL GETARG(4,CARG)
        READ(CARG,5000) K
        JF = 1
        JL = JDM
        CALL GETARG(5,CARG)
        READ(CARG,5000) IC
        CALL GETARG(6,CARG)
        READ(CARG,5000) IS
        CALL GETARG(7,CARG)
        READ(CARG,5000) IW
      ELSE
        WRITE(6,*) 
     +   'Usage:  hycom_meridional file.a idm jdm k [if il] [jc js jw]'
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
        WRITE(6,*) 'Error in hycom_meridional: could not allocate ',
     +             IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
C
      CALL MERID(A,IDM,JDM,PAD,NPAD, K,JF,JL,IC,IS,IW, CFILE)
      CALL EXIT(0)
 5000 FORMAT(I4)
      END
      SUBROUTINE MERID(A,IDM,JDM, PAD,NPAD, K,JF,JL,IC,IS,IW, CFILE)
      IMPLICIT NONE
C
      REAL*4     SPVAL
      PARAMETER (SPVAL=2.0**100)
C
      CHARACTER*240 CFILE
      INTEGER      IDM,JDM,NPAD,K,JF,JL,IC,IS,IW
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
      INTEGER      I,J,IF,IL,IZ,IZF,IZL,IOS,NRECL
      REAL*8       SUMA,SUMB,SUMR
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
     +    JF,JL,
     +    IC,IS,IW
C
        IF = IC - IS*((IC-  1)/IS)
        IL = IC + IS*((IDM-IC)/IS)
        DO IZ= IF,IL,IS
          IZF = MAX( IZ-IW/2,   1)
          IZL = MIN( IZ+IW/2, IDM)
          SUMA = 0.0
          SUMB = 0.0
          DO I= IZF,IZL
            DO J= JF,JL
              IF     (A(I,J).NE.SPVAL) THEN
                SUMA = SUMA + A(I,J)
                SUMB = SUMB + 1.D0
              ENDIF
            ENDDO
          ENDDO
          IF     (SUMB.NE.0.0) THEN
            SUMA = SUMA/SUMB
          ENDIF
C
          SUMR = 0.0
          DO I= IZF,IZL
            DO J= JF,JL
              IF     (A(I,J).NE.SPVAL) THEN
                SUMR = SUMR + (A(I,J)-SUMA)**2
              ENDIF
            ENDDO
          ENDDO
          IF     (SUMB.NE.0.0) THEN
            SUMR = SQRT(SUMR/SUMB)
          ENDIF
          WRITE(6,'(3I5,2x,1p2g16.8,a)') IZF,IZL,IZ,SUMA,SUMR
        ENDDO
      RETURN
 6000 FORMAT(
     +   '# hycom_meridional file.a idm jdm k if il jc js jw' /
     +   '# hycom_meridional ',A   /
     +   '#             ',8I5 /
     +   '#  IF   IL   IZ       MEAN            RMS')
      END

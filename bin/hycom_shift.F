      PROGRAM FSHIFT
      IMPLICIT NONE
C
C  hycom_shift - Usage:  hycom_shift f1.a f2.a idm jdm [ishift jshift [ip jp]]
C
C                 compares each (1:idm,1:jdm) array in f1.a and f2.a
C                 assuming the latter is periodically shifted by
C                 ishift and jshift in the 1st and 2nd dimensions.
C                 a1(ip,jp) is also printed if ip and jp are present
C
C                 primary use is for testing periodic cases, but
C                 can be used to compare bit-for-bit identity of 
C                 any two ".a" files by setting [ij]shift=0.
C
C                 ishift and jshift default to 0 if not present.
C
C  f[12].a is assumed to contain idm*jdm 32-bit IEEE real values for
C   each array, in standard f77 element order, followed by padding
C   to a multiple of 4096 32-bit words, but otherwise with no control
C   bytes/words, and input values of 2.0**100 indicating a data void.
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  August 2001.
C
      REAL*4, ALLOCATABLE :: A1(:,:),A2(:,:)
      REAL*4              :: PAD(4096)
C
      INTEGER      IOS
      INTEGER      IARGC
      INTEGER      NARG
      CHARACTER*240 CARG
C
      INTEGER      IDM,JDM,NPAD,ISHIFT,JSHIFT,IP,JP
      CHARACTER*240 CFILE1,CFILE2
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.4) THEN
        CALL GETARG(1,CFILE1)
        CALL GETARG(2,CFILE2)
        CALL GETARG(3,CARG)
        READ(CARG,5000) IDM
        CALL GETARG(4,CARG)
        READ(CARG,5000) JDM
        ISHIFT = 0
        JSHIFT = 0
        IP     = 0
        JP     = 0
      ELSEIF (NARG.EQ.6) THEN
        CALL GETARG(1,CFILE1)
        CALL GETARG(2,CFILE2)
        CALL GETARG(3,CARG)
        READ(CARG,5000) IDM
        CALL GETARG(4,CARG)
        READ(CARG,5000) JDM
        CALL GETARG(5,CARG)
        READ(CARG,5000) ISHIFT
        CALL GETARG(6,CARG)
        READ(CARG,5000) JSHIFT
        IP = 0
        JP = 0
      ELSEIF (NARG.EQ.8) THEN
        CALL GETARG(1,CFILE1)
        CALL GETARG(2,CFILE2)
        CALL GETARG(3,CARG)
        READ(CARG,5000) IDM
        CALL GETARG(4,CARG)
        READ(CARG,5000) JDM
        CALL GETARG(5,CARG)
        READ(CARG,5000) ISHIFT
        CALL GETARG(6,CARG)
        READ(CARG,5000) JSHIFT
        CALL GETARG(7,CARG)
        READ(CARG,5000) IP
        CALL GETARG(8,CARG)
        READ(CARG,5000) JP
      ELSE
        WRITE(6,*) 'Usage: ' //
     +   'hycom_shift f1.a f2.a idm jdm [ishift jshift [ip jp]]'
        CALL EXIT(1)
      ENDIF
C
      NPAD = 4096 - MOD(IDM*JDM,4096)
      IF     (NPAD.EQ.4096) THEN
        NPAD = 0
      ENDIF
C
      ALLOCATE( A1(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_shift: could not allocate ',
     +             IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( A2(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_shift: could not allocate ',
     +             IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
C
      CALL SHIFT(A1,A2,IDM,JDM,PAD,NPAD,
     +           ISHIFT,JSHIFT,IP,JP, CFILE1,CFILE2)
      CALL EXIT(0)
 5000 FORMAT(I4)
      END
      SUBROUTINE SHIFT(A1,A2,IDM,JDM, PAD,NPAD,
     +                 ISHIFT,JSHIFT,IP,JP, CFILE1,CFILE2)
      IMPLICIT NONE
C
      REAL*4     SPVAL
      PARAMETER (SPVAL=2.0**100)
C
      CHARACTER*240 CFILE1,CFILE2
      INTEGER      IDM,JDM,NPAD,ISHIFT,JSHIFT,IP,JP
      REAL*4       A1(IDM,JDM),A2(IDM,JDM),PAD(NPAD)
C
C     MOST OF WORK IS DONE HERE.
C
#ifdef sun
      INTEGER      IR_ISNAN
C
#endif
      CHARACTER*18 CIPJP,CASN
      INTEGER      LEN_TRIM
      INTEGER      I1,I2,J1,J2,K,IOS,KBAD,NBAD,NRECL
      REAL*4       AMN,AMX
#ifdef CRAY
      INTEGER*8    IU8,IOS8
#endif
C
      IF     (MIN(IP,JP).GT.0) THEN
        WRITE(CIPJP,6000) IP,JP
 6000   FORMAT('    (',I4,',',I4,') = ')
      ENDIF
C
      INQUIRE( IOLENGTH=NRECL) A1,PAD
#ifdef CRAY
#ifdef t3e
      IF     (MOD(NRECL,4096).EQ.0) THEN
        WRITE(CASN,8000) NRECL/4096
 8000   FORMAT('-F cachea:',I4.4,':1:0')
        IU8 = 11
        CALL ASNUNIT(IU8,CASN,IOS8)
        IF     (IOS8.NE.0) THEN
          write(6,*) 'Error: can''t asnunit ',IU8
          write(6,*) 'ios  = ',ios8
          write(6,*) 'casn = ',casn
          CALL EXIT(5)
        ENDIF
        IU8 = 12
        CALL ASNUNIT(IU8,CASN,IOS8)
        IF     (IOS8.NE.0) THEN
          write(6,*) 'Error: can''t asnunit ',IU8
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
      CALL ASNUNIT(12,'-F syscall -N ieee',IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t asnunit 12'
        write(6,*) 'ios = ',ios
        CALL EXIT(5)
      ENDIF
#endif
#endif
      OPEN(UNIT=11, FILE=CFILE1, FORM='UNFORMATTED', STATUS='OLD',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',CFILE1(1:LEN_TRIM(CFILE1))
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
      OPEN(UNIT=12, FILE=CFILE2, FORM='UNFORMATTED', STATUS='OLD',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',CFILE2(1:LEN_TRIM(CFILE2))
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
C
      KBAD = 0
      DO 110 K= 1,9999
        READ(11,REC=K,IOSTAT=IOS) A1
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(A1,IDM*JDM)
#endif
        IF     (IOS.NE.0) THEN
          IF     (K.EQ.1) THEN
            WRITE(6,*) 'can''t read ',CFILE1(1:LEN_TRIM(CFILE1))
            CALL EXIT(4)
          ELSE
            GOTO 1110
          ENDIF
        ENDIF
        READ(12,REC=K,IOSTAT=IOS) A2
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(A2,IDM*JDM)
#endif
        IF     (IOS.NE.0) THEN
          IF     (K.EQ.1) THEN
            WRITE(6,*) 'can''t read ',CFILE2(1:LEN_TRIM(CFILE2))
            CALL EXIT(4)
          ELSE
            WRITE(6,*) 'files are not the same length'
            CALL EXIT(4)
          ENDIF
        ENDIF
        NBAD = 0
        AMN  =  SPVAL
        AMX  = -SPVAL
        DO 210 J1= 1,JDM
          J2 = MOD( J1-1+JSHIFT+JDM, JDM) + 1
          DO 212 I1= 1,IDM
            I2 = MOD( I1-1+ISHIFT+IDM, IDM) + 1
C
            IF     (A1(I1,J1).NE.A2(I2,J2)) THEN
*             IF     (NBAD.EQ.0) THEN
*               WRITE(6,'(a,2i5,1pg16.8,2i5,1p2g16.8)')
*    &            'A1,A2 = ',I1,J1,A1(I1,J1),
*    &                       I2,J2,A2(I2,J2),A1(I1,J1)-A2(I2,J2)
*             ENDIF
              NBAD = NBAD + 1
            ENDIF
C
#ifdef sun
            IF     (IR_ISNAN(A1(I1,J1)).NE.1) THEN
              IF     (A1(I1,J1).NE.SPVAL) THEN
                AMN = MIN( AMN, A1(I1,J1) )
                AMX = MAX( AMX, A1(I1,J1) )
              ENDIF
            ENDIF
#else
            IF     (A1(I1,J1).NE.SPVAL) THEN
              AMN = MIN( AMN, A1(I1,J1) )
              AMX = MAX( AMX, A1(I1,J1) )
            ENDIF
#endif
  212     CONTINUE
  210   CONTINUE
        IF     (NBAD.EQ.0) THEN
          IF     (MIN(IP,JP).LE.0) THEN
            WRITE(6,'(a,1p2g16.8)')
     &         'min, max = ',AMN,AMX
          ELSE
            WRITE(6,'(a,1p2g16.8,a,g16.8)')
     &        'min, max = ',AMN,AMX,CIPJP,A1(IP,JP)
          ENDIF
        ELSE
          KBAD = KBAD + 1
          IF     (MIN(IP,JP).LE.0) THEN
            WRITE(6,'(a,1p2g16.8,a,I9)')
     &         'min, max = ',AMN,AMX,
     &         '  no. diff =',NBAD
          ELSE
            WRITE(6,'(a,1p2g16.8,a,g16.8,a,I9)')
     &        'min, max = ',AMN,AMX,CIPJP,A1(IP,JP),
     &         '  no. diff =',NBAD
          ENDIF
        ENDIF
  110 CONTINUE
 1110 CONTINUE
      WRITE(6,*) 
      WRITE(6,*) K-1,' FIELDS PROCESSED, ',KBAD,' WITH DIFFS'
      WRITE(6,*) 
      RETURN
      END

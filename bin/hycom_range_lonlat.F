      PROGRAM RANGELL
      IMPLICIT NONE
C
C  hycom_range_lonlat - Usage:  hycom_range_lonlat file.a [grid.a]
C
C                 prints lonlat location and min/max value for each
C                 (1:idm,1:jdm) array in file.a
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
      REAL*4, ALLOCATABLE :: PLAT(:,:),PLON(:,:), A(:,:)
      REAL*4              :: PAD(4096)
      INTEGER IOS
      INTEGER       IARGC
      INTEGER       NARG
      CHARACTER*240 CARG
C
      INTEGER       IDM,JDM,NPAD
      CHARACTER*240 CFILE,CFILEG,CFILEB
      CHARACTER*6   CVARIN
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.1) THEN
        CALL GETARG(1,CFILE)
        CFILEG = 'regional.grid.a'
      ELSEIF (NARG.EQ.2) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CFILEG)
      ELSE
        WRITE(6,*) 'Usage: hycom_range_lonlat file.a [grid.a]'
        CALL EXIT(1)
      ENDIF
C
C     GET IDM,JDM FROM grid.b.
C
      CFILEB = CFILEG(1:LEN_TRIM(CFILEG)-1) // 'b'
      OPEN(UNIT=22,FILE=CFILEB,FORM='FORMATTED',
     &     STATUS='OLD',ACTION='READ')
      READ( 22,*) IDM,CVARIN
      IF (CVARIN.NE.'idm   ') THEN
        WRITE(6,*) 'hycom_range_lonlat: bad header file ',
     &             CFILEB(1:LEN_TRIM(CFILEB))
        CALL EXIT(2)
      ENDIF
      READ( 22,*) JDM,CVARIN
      IF (CVARIN.NE.'jdm   ') THEN
        WRITE(6,*) 'hycom_range_lonlat: bad header file ',
     &             CFILEB(1:LEN_TRIM(CFILEB))
        CALL EXIT(2)
      ENDIF
      CLOSE(UNIT=22)

C
      NPAD = 4096 - MOD(IDM*JDM,4096)
      IF     (NPAD.EQ.4096) THEN
        NPAD = 0
      ENDIF
C
      ALLOCATE( PLON(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_range_lonlat: could not allocate ',
     +             IDM*JDM,' words for PLON'
        CALL EXIT(3)
      ENDIF
      ALLOCATE( PLAT(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_range_lonlat: could not allocate ',
     +             IDM*JDM,' words for PLAT'
        CALL EXIT(3)
      ENDIF
      ALLOCATE( A(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_range_lonlat: could not allocate ',
     +             IDM*JDM,' words for A'
        CALL EXIT(2)
      ENDIF
C
      CALL RANGE(A,PLON,PLAT,IDM,JDM,PAD,NPAD, CFILE,CFILEG)
      CALL EXIT(0)
      END
      SUBROUTINE RANGE(A,PLON,PLAT,IDM,JDM,PAD,NPAD, CFILE,CFILEG)
      IMPLICIT NONE
C
      REAL*4     SPVAL
      PARAMETER (SPVAL=2.0**100)
C
      CHARACTER*240 CFILE,CFILEG
      INTEGER      IDM,JDM,NPAD
      REAL*4       A(IDM,JDM),PAD(NPAD)
      REAL*4       PLON(IDM,JDM),PLAT(IDM,JDM)
C
C     MOST OF WORK IS DONE HERE.
C
#ifdef sun
      INTEGER      IR_ISNAN
C
#endif
      CHARACTER*99 CLINE
      CHARACTER*30 CIJMAX,CIJMIN,CASN
      INTEGER      I,IMN,IMX,J,JMN,JMX,K,IOS,NRECL
      REAL*4       AMN,AMX,X,Y
#ifdef CRAY
      INTEGER*8    IU8,IOS8
#endif
      INQUIRE( IOLENGTH=NRECL) A,PAD
C
      OPEN(UNIT=21, FILE=CFILEG, FORM='UNFORMATTED', STATUS='OLD',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILEG)
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
      READ(21,REC=1,IOSTAT=IOS) PLON
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(PLON,IDM*JDM)
#endif
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'can''t read ',TRIM(CFILEG)
        CALL EXIT(4)
      ENDIF
      READ(21,REC=2,IOSTAT=IOS) PLAT
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(PLAT,IDM*JDM)
#endif
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'can''t read ',TRIM(CFILEG)
        CALL EXIT(4)
      ENDIF
      CLOSE(UNIT=21)
C
      OPEN(UNIT=11, FILE=CFILE, FORM='UNFORMATTED', STATUS='OLD',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILE)
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
C
      DO 110 K= 1,9999
        READ(11,REC=K,IOSTAT=IOS) A
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
        IF     (IOS.NE.0) THEN
          IF     (K.EQ.1) THEN
            WRITE(6,*) 'can''t read ',TRIM(CFILE)
            CALL EXIT(4)
          ELSE
            GOTO 1110
          ENDIF
        ENDIF
        AMN =  SPVAL
        IMN =  0
        JMN =  0
        AMX = -SPVAL
        IMX =  0
        JMX =  0
        DO 210 J= 1,JDM
          DO 212 I= 1,IDM
#ifdef sun
            IF     (IR_ISNAN(A(I,J)).NE.1) THEN
              IF     (A(I,J).NE.SPVAL) THEN
                IF     (A(I,J).LT.AMN) THEN
                  AMN = A(I,J)
                  IMN = I
                  JMN = J
                ENDIF
                IF     (A(I,J).GT.AMX) THEN
                  AMX = A(I,J)
                  IMX = I
                  JMX = J
                ENDIF
              ENDIF
            ENDIF
#else
            IF     (A(I,J).NE.SPVAL) THEN
              IF     (A(I,J).LT.AMN) THEN
                AMN = A(I,J)
                IMN = I
                JMN = J
              ENDIF
              IF     (A(I,J).GT.AMX) THEN
                AMX = A(I,J)
                IMX = I
                JMX = J
              ENDIF
            ENDIF
#endif
  212     CONTINUE
  210   CONTINUE
C
        Y = PLAT(IMN,JMN)
        X = PLON(IMN,JMN)
        X = MOD(X+180.0,360.0)-180.0
        IF     (X.LT.0.0) THEN
          X = X + 360.0
        ENDIF
        WRITE(CIJMIN,6100) 'min',X,Y
C
        Y = PLAT(IMX,JMX)
        X = PLON(IMX,JMX)
        X = MOD(X+180.0,360.0)-180.0
        IF     (X.LT.0.0) THEN
          X = X + 360.0
        ENDIF
        WRITE(CIJMAX,6100) 'max',X,Y
 6100   FORMAT(A,' at (',F6.2,',',F6.2,') = ')
        WRITE(CLINE,'(a,1pg16.8,2x,a,g16.8)') CIJMIN,AMN,CIJMAX,AMX
        WRITE(6,'(a)') TRIM(CLINE)
  110 CONTINUE
 1110 CONTINUE
      CLOSE(UNIT=11)
      WRITE(6,*) 
      WRITE(6,*) K-1,' FIELDS PROCESSED'
      WRITE(6,*) 
      RETURN
      END

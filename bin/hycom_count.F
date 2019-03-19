      PROGRAM FMASS
      IMPLICIT NONE
C
C  hycom_count     - Usage:  hycom_count     file.a [[grid.a] mask.a]
C  hycom_count_dim - Usage:  hycom_count_dim file.a [mask.a] idm jdm
C
C                 prints the count of sea points of all fields in file.a
C
C                 grid.a is a hycom grid file, default regional.grid.a.
C                 Note that the corresponding grid.b must also exist.
C
C                 mask.a contains an optional mask array, that turns
C                 off sampling anywhere it is 2.0**100 (data void).
C                 Note that grid.a must be specified if mask.a is used.
C
C                 for hycom_count idm,jdm are taken from grid.a.
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
      REAL*4, ALLOCATABLE :: A(:,:),AMSK(:,:)
      REAL*4              :: PAD(4096)
      INTEGER       IOS,L
      INTEGER       IARGC
      INTEGER       NARG
      CHARACTER*240 CARG
C
      LOGICAL       LDIM
      INTEGER       IDM,JDM,NPAD
      REAL*4        CELL
      CHARACTER*6   CVARIN
      CHARACTER*240 CFILE,CFILEG,CFILEM,CFILEB
C
C     READ ARGUMENTS.
C
      CALL GETARG(0,CARG)
      L = LEN_TRIM(CARG)
      LDIM = CARG(L-3:L).EQ.'_dim'

      NARG = IARGC()
C
      IF     (LDIM) THEN
        IF     (NARG.EQ.3) THEN
          CALL GETARG(1,CFILE)
          CALL GETARG(2,CARG)
          READ(CARG,*) IDM
          CALL GETARG(3,CARG)
          READ(CARG,*) JDM
          CFILEM = 'NONE'
        ELSEIF (NARG.EQ.4) THEN
          CALL GETARG(1,CFILE)
          CALL GETARG(2,CFILEG)
          CALL GETARG(3,CARG)
          READ(CARG,*) IDM
          CALL GETARG(4,CARG)
          READ(CARG,*) JDM
        ELSE
          WRITE(6,*) 
     +     'Usage:  hycom_count_dim file.a [grid.a] idm jdm'
          CALL EXIT(1)
        ENDIF
      ELSE
        IF     (NARG.EQ.1) THEN
          CALL GETARG(1,CFILE)
          CFILEG = 'regional.grid.a'
          CFILEM = 'NONE'
        ELSEIF (NARG.EQ.2) THEN
          CALL GETARG(1,CFILE)
          CALL GETARG(2,CFILEG)
          CFILEM = 'NONE'
        ELSEIF (NARG.EQ.3) THEN
          CALL GETARG(1,CFILE)
          CALL GETARG(2,CFILEG)
          CALL GETARG(3,CFILEM)
        ELSE
          WRITE(6,*) 
     +     'Usage:  hycom_count file.a [grid.a] [mask.a]'
          CALL EXIT(1)
        ENDIF
      ENDIF
C
C     GET IDM,JDM FROM grid.b.
C
      IF     (.NOT.LDIM) THEN
        CFILEB = CFILEG(1:LEN_TRIM(CFILEG)-1) // 'b'
C
        OPEN(UNIT=11,FILE=CFILEB,FORM='FORMATTED',
     &       STATUS='OLD',ACTION='READ')
C
        READ( 11,*) IDM,CVARIN
        IF (CVARIN.NE.'idm   ') THEN
          WRITE(6,*) 'hycom_count: bad header file ',
     &               CFILEB(1:LEN_TRIM(CFILEB))
          CALL EXIT(2)
        ENDIF
        READ( 11,*) JDM,CVARIN
        IF (CVARIN.NE.'jdm   ') THEN
          WRITE(6,*) 'hycom_count: bad header file ',
     &                 CFILEB(1:LEN_TRIM(CFILEB))
          CALL EXIT(2)
        ENDIF
C
        CLOSE(UNIT=11)
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
     +             IDM*JDM,' words for A'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( AMSK(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_zonal: could not allocate ',
     +             IDM*JDM,' words for AMSK'
        CALL EXIT(2)
      ENDIF
C
      CALL MCNT(A,AMSK,IDM,JDM,PAD,NPAD, 
     +          CFILE,CFILEG,CFILEM)
      CALL EXIT(0)
 5000 FORMAT(I4)
      END
      SUBROUTINE MCNT(A,AMSK,IDM,JDM, PAD,NPAD,
     +                CFILE,CFILEG,CFILEM)
      IMPLICIT NONE
C
      REAL*4     SPVAL
      PARAMETER (SPVAL=2.0**100)
C
      CHARACTER*240 CFILE,CFILEG,CFILEM
      INTEGER      IDM,JDM,NPAD
      REAL*4       A(IDM,JDM),AMSK(IDM,JDM),PAD(NPAD)
C
C     MOST OF WORK IS DONE HERE.
C
#ifdef sun
      INTEGER      IR_ISNAN
C
#endif
      CHARACTER*18 CASN
      INTEGER      LEN_TRIM
      INTEGER      COUNT,I,J,K,IOS,NRECL
#ifdef CRAY
      INTEGER*8    IU8,IOS8
#endif
C
C     INPUT MASK ARRAY.
C
      INQUIRE( IOLENGTH=NRECL) A,PAD
C
      IF     (CFILEM.EQ.'NONE') THEN
        AMSK(:,:) = 0.0
      ELSE
        OPEN(UNIT=11, FILE=CFILEM, FORM='UNFORMATTED', STATUS='OLD',
     +           ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
        IF     (IOS.NE.0) THEN
          write(6,*) 'Error: can''t open ',TRIM(CFILEM)
          write(6,*) 'ios   = ',ios
          write(6,*) 'nrecl = ',nrecl
          CALL EXIT(3)
        ENDIF
C
        READ(11,REC=1,IOSTAT=IOS) AMSK
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(AMSK,IDM*JDM)
#endif
        IF     (IOS.NE.0) THEN
          WRITE(6,*) 'can''t read ',TRIM(CFILEM)
          CALL EXIT(4)
        ENDIF
      ENDIF
C
      CLOSE(UNIT=11)
C
C     OPEN INPUT FILE.
C
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
        write(6,*) 'Error: can''t open ',TRIM(CFILE)
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
C
C     LOOP OVER ALL FIELDS.
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
C
        COUNT = 0
        DO J= 1,JDM
          DO I= 1,IDM
            IF     (   A(I,J).NE.SPVAL .AND.
     +              AMSK(I,J).NE.SPVAL      ) THEN
              COUNT = COUNT + 1
            ENDIF
          ENDDO
        ENDDO
        WRITE(6,'(a,i8,f8.3)')
     &    'count,percent = ',COUNT,100.0*REAL(COUNT)/REAL(IDM*JDM)
  110 CONTINUE
 1110 CONTINUE
      WRITE(6,*)
      WRITE(6,*) K-1,' FIELDS PROCESSED'
      WRITE(6,*)
C
      CLOSE(UNIT=11)
C
      RETURN
      END

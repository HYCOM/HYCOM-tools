      PROGRAM HYCOM_SSTICE
      IMPLICIT NONE
C
C  hycom_sstice - Usage:  hycom_sstice sst.a ice.a idm jdm ssti.a [maxice]
C
C                  Outputs SST with -1.8 under ice
C
C  ice is zero for no ice and maxice (default 1.0) for all ice
C
C  *.a is assumed to contain idm*jdm 32-bit IEEE real values for
C   each array, in standard f77 element order, followed by padding
C   to a multiple of 4096 32-bit words, but otherwise with no control
C   bytes/words, and input values of 2.0**100 indicating a data void.
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  July 2005.
C
      REAL*4, ALLOCATABLE :: SST(:,:),ICE(:,:)
      REAL*4              :: PAD(4096)
      INTEGER       IOS
      INTEGER       IARGC
      INTEGER       NARG
      CHARACTER*240 CARG
C
      REAL*4        MAXICE
      INTEGER       IDM,JDM,NPAD
      CHARACTER*240 CFILE1,CFILE2,CFILEO
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.6) THEN
        CALL GETARG(1,CFILE1)
        CALL GETARG(2,CFILE2)
        CALL GETARG(3,CARG)
        READ(CARG,*) IDM
        CALL GETARG(4,CARG)
        READ(CARG,*) JDM
        CALL GETARG(5,CFILEO)
        CALL GETARG(6,CARG)
        READ(CARG,*) MAXICE
      ELSEIF (NARG.EQ.5) THEN
        CALL GETARG(1,CFILE1)
        CALL GETARG(2,CFILE2)
        CALL GETARG(3,CARG)
        READ(CARG,*) IDM
        CALL GETARG(4,CARG)
        READ(CARG,*) JDM
        CALL GETARG(5,CFILEO)
        MAXICE = 1.0
      ELSE
        WRITE(6,*)
     &    'Usage: hycom_sstice sst.a ice.a idm jdm ssti.a [maxice]'
        CALL EXIT(1)
      ENDIF
C
      NPAD = 4096 - MOD(IDM*JDM,4096)
      IF     (NPAD.EQ.4096) THEN
        NPAD = 0
      ENDIF
C
      ALLOCATE( SST(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_sstice: could not allocate 1st ',
     +             IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( ICE(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_sstice: could not allocate 2nd ',
     +             IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
C
      CALL SSTI(SST,ICE,IDM,JDM,PAD,NPAD, 
     &          MAXICE, CFILE1,CFILE2,CFILEO)
      CALL EXIT(0)
      END
      SUBROUTINE SSTI(SST,ICE,IDM,JDM,PAD,NPAD,
     &                MAXICE, CFILE1,CFILE2,CFILEO)
      IMPLICIT NONE
C
      REAL*4     SPVAL
      PARAMETER (SPVAL=2.0**100)
C
      CHARACTER*240 CFILE1,CFILE2,CFILEO
      INTEGER       IDM,JDM,NPAD
      REAL*4        SST(IDM,JDM),ICE(IDM,JDM),PAD(NPAD)
      REAL*4        MAXICE
C
C     MOST OF WORK IS DONE HERE.
C
#ifdef sun
      INTEGER      IR_ISNAN
C
#endif
      CHARACTER*18 CASN
      LOGICAL      LMULT
      INTEGER      LEN_TRIM
      INTEGER      I,J,K,IOS,NRECL
      REAL*4       AMN,AMX,RICE,S
#ifdef CRAY
      INTEGER*8    IU8,IOS8
#endif
C
      INQUIRE( IOLENGTH=NRECL) SST,PAD
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
        IU8 = 12
        CALL ASNUNIT(IU8,CASN,IOS8)
        IF     (IOS8.NE.0) THEN
          write(6,*) 'Error: can''t asnunit 12'
          write(6,*) 'ios  = ',ios8
          write(6,*) 'casn = ',casn
          CALL EXIT(5)
        ENDIF
        IU8 = 21
        CALL ASNUNIT(IU8,CASN,IOS8)
        IF     (IOS8.NE.0) THEN
          write(6,*) 'Error: can''t asnunit 21'
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
      CALL ASNUNIT(21,'-F syscall -N ieee',IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t asnunit 21'
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
      OPEN(UNIT=21, FILE=CFILEO, FORM='UNFORMATTED', STATUS='NEW',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',CFILEO(1:LEN_TRIM(CFILEO))
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
C
      RICE = 1.0/MAXICE
      DO 110 K= 1,9999
        READ(11,REC=K,IOSTAT=IOS) SST
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(SST,IDM*JDM)
#endif
        IF     (IOS.NE.0) THEN
          IF     (K.EQ.1) THEN
            WRITE(6,*) 'can''t read ',CFILE1(1:LEN_TRIM(CFILE1))
            CALL EXIT(4)
          ELSE
            GOTO 1110
          ENDIF
        ENDIF
        READ(12,REC=K,IOSTAT=IOS) ICE
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(ICE,IDM*JDM)
#endif
        IF     (IOS.NE.0) THEN
          IF     (K.EQ.1) THEN
            WRITE(6,*) 'can''t read ',CFILE2(1:LEN_TRIM(CFILE2))
            CALL EXIT(4)
          ELSE
            GOTO 1110
          ENDIF
        ENDIF
C
        AMN =  SPVAL
        AMX = -SPVAL

        DO 210 J= 1,JDM
          DO 212 I= 1,IDM
#ifdef sun
            IF     (IR_ISNAN(SST(I,J)).NE.1) THEN
              IF     (SST(I,J).NE.SPVAL .AND.
     &                ICE(I,J).NE.SPVAL      ) THEN
                S        = ICE(I,J)*RICE  !0.0 to 1.0
                SST(I,J) = (1.0-S)*MAX(SST(I,J),-1.8) - 1.8*S
                AMN = MIN( AMN, SST(I,J) )
                AMX = MAX( AMX, SST(I,J) )
              ELSE
                SST(I,J) = SPVAL
              ENDIF
            ENDIF
#else
            IF     (SST(I,J).NE.SPVAL .AND.
     &              ICE(I,J).NE.SPVAL      ) THEN
              S        = MAX( 0.0, MIN( 1.0, ICE(I,J)*RICE ) )
              SST(I,J) = (1.0-S)*MAX(SST(I,J),-1.8) - 1.8*S
              AMN = MIN( AMN, SST(I,J) )
              AMX = MAX( AMX, SST(I,J) )
            ELSE
              SST(I,J) = SPVAL
            ENDIF
#endif
  212     CONTINUE
  210   CONTINUE
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(SST,IDM*JDM)
#endif
        WRITE(21,REC=K,IOSTAT=IOS) SST
        WRITE(6,'(a,1p2g16.8)')
     &     'ssti: min, max = ',AMN,AMX
  110 CONTINUE
 1110 CONTINUE
      WRITE(6,*) 
      WRITE(6,*) K-1,' FIELDS PROCESSED'
      WRITE(6,*) 
C
      CLOSE(11)
      CLOSE(12)
      CLOSE(21)
C
      RETURN
      END

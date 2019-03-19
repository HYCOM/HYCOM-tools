      PROGRAM HYCOM_CFL
      IMPLICIT NONE
C
C  hycom_cfl      - Usage:  hycom_cfl     speed fout.a [grid.a]
C  hycom_cfl_ts   - Usage:  hycom_cfl_ts  speed fout.a [grid.a]
C  hycom_cfl_spd  - Usage:  hycom_cfl_spd ts fout.a [grid.a]
C
C                 Returns, in fout.a, by default and for "_fs,
C                   the maximum time step at each grid point for the 
C                   input advective speed to be stable
C                 Returns, in fout.a, for "_spd,
C                   the maximum speed at each grid point for the 
C                   input time step to be stable
C
C  fout.a will contain idm*jdm 32-bit IEEE real values for
C   each array, in standard f77 element order, followed by padding
C   to a multiple of 4096 32-bit words, but otherwise with no control
C   bytes/words, and input values of 2.0**100 indicating a data void.
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  COAPS (FSU), October 2018
C
      REAL*4, ALLOCATABLE :: STEP(:,:),PSCX(:,:),PSCY(:,:)
      REAL*4              :: PAD(4096)
      INTEGER       IOS,L
      INTEGER       IARGC
      INTEGER       NARG
      CHARACTER*240 CARG
C
      LOGICAL       LSPD
      INTEGER       IDM,JDM,NPAD
      REAL*4        SPEED
      CHARACTER*6   CVARIN
      CHARACTER*240 CFILEG,CFILEB,CFILEO
C
C     READ ARGUMENTS.
C
      CALL GETARG(0,CARG)
      L = LEN_TRIM(CARG)
      LSPD = CARG(L-3:L).EQ.'_spd'
C
      NARG = IARGC()
C
      IF     (NARG.EQ.3) THEN
        CALL GETARG(1,CARG)
        READ(CARG,*)   SPEED   !speed or t
        CALL GETARG(2,CFILEO)
        CALL GETARG(3,CFILEG)
      ELSEIF (NARG.EQ.2) THEN
        CALL GETARG(1,CARG)
        READ(CARG,*)   SPEED   !speed or t
        CALL GETARG(2,CFILEO)
        CFILEG = 'regional.grid.a'
      ELSEIF (LSPD) THEN
        WRITE(6,*) 
     +   'Usage:  hycom_cfl_spd ts fout.a [grid.a]'
        CALL EXIT(1)
      ELSE
        WRITE(6,*) 
     +   'Usage:  hycom_cfl speed fout.a [grid.a]'
        CALL EXIT(1)
      ENDIF
C
C     GET IDM,JDM FROM grid.b.
C
      CFILEB = CFILEG(1:LEN_TRIM(CFILEG)-1) // 'b'
C
      OPEN(UNIT=11,FILE=CFILEB,FORM='FORMATTED',
     &     STATUS='OLD',ACTION='READ')
C
      READ( 11,*) IDM,CVARIN
      IF (CVARIN.NE.'idm   ') THEN
        WRITE(6,*) 'hycom_cfl: bad header file ',
     &             CFILEB(1:LEN_TRIM(CFILEB))
        CALL EXIT(2)
      ENDIF
      READ( 11,*) JDM,CVARIN
      IF (CVARIN.NE.'jdm   ') THEN
        WRITE(6,*) 'hycom_cfl: bad header file ',
     &             CFILEB(1:LEN_TRIM(CFILEB))
        CALL EXIT(2)
      ENDIF
C
      CLOSE(UNIT=11)
C
      NPAD = 4096 - MOD(IDM*JDM,4096)
      IF     (NPAD.EQ.4096) THEN
        NPAD = 0
      ENDIF
C
      ALLOCATE( PSCX(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_cfl: could not allocate ',
     +             IDM*JDM,' words for PSCX'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( PSCY(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_cfl: could not allocate ',
     +             IDM*JDM,' words for PSCY'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( STEP(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_cfl: could not allocate ',
     +             IDM*JDM,' words for STEP'
        CALL EXIT(2)
      ENDIF
C
      IF     (LSPD) THEN
        CALL CFL_SPD(STEP,PSCX,PSCY,IDM,JDM,PAD,NPAD,  !speed,pscx,pscy....
     +               SPEED, CFILEG,CFILEO)             !step, cfileg,cfileo
      ELSE
        CALL CFL_TS( STEP,PSCX,PSCY,IDM,JDM,PAD,NPAD, 
     +               SPEED, CFILEG,CFILEO)
      ENDIF
      CALL EXIT(0)
 5000 FORMAT(I4)
      END
      SUBROUTINE CFL_SPD(SPEED,PSCX,PSCY,IDM,JDM,PAD,NPAD,        
     +                   STEP, CFILEG,CFILEO)
      IMPLICIT NONE
C
      REAL*4     SPVAL
      PARAMETER (SPVAL=2.0**100)
C
      CHARACTER*240 CFILEG,CFILEO
      INTEGER       IDM,JDM,NPAD
      REAL*4        PSCX( IDM,JDM),PSCY(IDM,JDM),
     +              SPEED(IDM,JDM),PAD(NPAD)
      REAL*4        STEP
C
C     MOST OF WORK IS DONE HERE.
C
#ifdef sun
      INTEGER      IR_ISNAN
C
#endif
      CHARACTER*18 CASN
      INTEGER      I,J,IOS,NRECL
      REAL*4       AMN,AMX
#ifdef CRAY
      INTEGER*8    IU8
#endif
C
C     INPUT GRID ARRAYS.
C
      IF     (NPAD.EQ.0) THEN
        INQUIRE( IOLENGTH=NRECL) PSCX
      ELSE
        INQUIRE( IOLENGTH=NRECL) PSCX,PAD
        PAD(:) = SPVAL
      ENDIF

      OPEN(UNIT=11, FILE=CFILEG, FORM='UNFORMATTED', STATUS='OLD',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILEG)
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
C
      READ(11,REC=10,IOSTAT=IOS) PSCX  ! pscx
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(PSCX,IDM*JDM)
#endif
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'can''t read pscx from ',TRIM(CFILEG)
        CALL EXIT(4)
      ENDIF
C
      READ(11,REC=11,IOSTAT=IOS) PSCY  ! pscy
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(PSCY,IDM*JDM)
#endif
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'can''t read pscy from ',TRIM(CFILEG)
        CALL EXIT(4)
      ENDIF
C
      CLOSE(UNIT=11)
C
C     OPEN OUTPUT FILE.
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
      OPEN(UNIT=21, FILE=CFILEO, FORM='UNFORMATTED', STATUS='NEW',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILEO)
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
C
C --- CALCULATE THE SPEED AT THE CFL LIMIT FOR THE TIMESTEP AT EACH GRIDPOINT.
C
      AMN  =  SPVAL
      AMX  = -SPVAL
      DO J= 1,JDM
        DO I= 1,IDM
          SPEED(I,J) = 0.5*MIN(PSCX(I,J),PSCY(I,J))/STEP
          AMN  = MIN( AMN, SPEED(I,J) )
          AMX  = MAX( AMX, SPEED(I,J) )
        ENDDO
      ENDDO
C
      IF     (NPAD.EQ.0) THEN
        WRITE(21,REC=1,IOSTAT=IOS) SPEED
      ELSE
        WRITE(21,REC=1,IOSTAT=IOS) SPEED,PAD
      ENDIF
      WRITE(6,'(a,1p2g16.8)')
     &     'CFL SPEED: min, max = ',AMN,AMX
C
      CLOSE(UNIT=21)
C
      RETURN
      END
      SUBROUTINE CFL_TS(STEP,PSCX,PSCY,IDM,JDM,PAD,NPAD,        
     +                  SPEED, CFILEG,CFILEO)
      IMPLICIT NONE
C
      REAL*4     SPVAL
      PARAMETER (SPVAL=2.0**100)
C
      CHARACTER*240 CFILEG,CFILEO
      INTEGER       IDM,JDM,NPAD
      REAL*4        PSCX(IDM,JDM),PSCY(IDM,JDM),
     +              STEP(IDM,JDM),PAD(NPAD)
      REAL*4        SPEED
C
C     MOST OF WORK IS DONE HERE.
C
#ifdef sun
      INTEGER      IR_ISNAN
C
#endif
      CHARACTER*18 CASN
      INTEGER      I,J,IOS,NRECL
      REAL*4       AMN,AMX
#ifdef CRAY
      INTEGER*8    IU8
#endif
C
C     INPUT GRID ARRAYS.
C
      IF     (NPAD.EQ.0) THEN
        INQUIRE( IOLENGTH=NRECL) PSCX
      ELSE
        INQUIRE( IOLENGTH=NRECL) PSCX,PAD
        PAD(:) = SPVAL
      ENDIF

      OPEN(UNIT=11, FILE=CFILEG, FORM='UNFORMATTED', STATUS='OLD',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILEG)
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
C
      READ(11,REC=10,IOSTAT=IOS) PSCX  ! pscx
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(PSCX,IDM*JDM)
#endif
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'can''t read pscx from ',TRIM(CFILEG)
        CALL EXIT(4)
      ENDIF
C
      READ(11,REC=11,IOSTAT=IOS) PSCY  ! pscy
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(PSCY,IDM*JDM)
#endif
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'can''t read pscy from ',TRIM(CFILEG)
        CALL EXIT(4)
      ENDIF
C
      CLOSE(UNIT=11)
C
C     OPEN OUTPUT FILE.
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
      OPEN(UNIT=21, FILE=CFILEO, FORM='UNFORMATTED', STATUS='NEW',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILEO)
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
C
C --- CALCULATE THE CFL TIMESTEP AT EACH GRIDPOINT.
C
      AMN  =  SPVAL
      AMX  = -SPVAL
      DO J= 1,JDM
        DO I= 1,IDM
          STEP(I,J) = 0.5*MIN(PSCX(I,J),PSCY(I,J))/SPEED
          AMN  = MIN( AMN, STEP(I,J) )
          AMX  = MAX( AMX, STEP(I,J) )
        ENDDO
      ENDDO
C
      IF     (NPAD.EQ.0) THEN
        WRITE(21,REC=1,IOSTAT=IOS) STEP
      ELSE
        WRITE(21,REC=1,IOSTAT=IOS) STEP,PAD
      ENDIF
      WRITE(6,'(a,1p2g16.8)')
     &     'CFL TIMESTEP: min, max = ',AMN,AMX
C
      CLOSE(UNIT=21)
C
      RETURN
      END

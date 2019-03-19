      PROGRAM HYCOM_ICEFREEDAY
      IMPLICIT NONE
C
C  hycom_icefreeday - Usage:  hycom_icefreeday ext.a idm jdm day1 dayi free.a
C
C           Outputs 1 (1:idm,1:jdm) field, representing the first
C           ice free (value 0.0) day in the multi-record ext.a
C
C  ext.a is assumed to contain idm*jdm 32-bit IEEE real values for
C   each array, in standard f77 element order, followed by padding
C   to a multiple of 4096 32-bit words, but otherwise with no control
C   bytes/words, and input values of 2.0**100 indicating a data void.
C   it has 1.0 for sea ice and 0.0 for no sea ice.
C
C  if record number n is the first ice free record at a location, 
C   free.a is set to day1+(n-1)*dayi there.
C  the result will be data void where all records have sea ice.
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  June 1015.
C
      REAL*4, ALLOCATABLE :: EXT(:,:),FREE(:,:)
      REAL*4              :: PAD(4096)
      INTEGER       IOS
      INTEGER       IARGC
      INTEGER       NARG
      CHARACTER*240 CARG
C
      INTEGER       IDM,JDM,NPAD
      REAL          DAY1,DAYI
      CHARACTER*240 CFILEI,CFILEO
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.6) THEN
        CALL GETARG(1,CFILEI)
        CALL GETARG(2,CARG)
        READ(CARG,*) IDM
        CALL GETARG(3,CARG)
        READ(CARG,*) JDM
        CALL GETARG(4,CARG)
        READ(CARG,*) DAY1
        CALL GETARG(5,CARG)
        READ(CARG,*) DAYI
        CALL GETARG(6,CFILEO)
      ELSE
        WRITE(6,'(2a)')
     &    'Usage:  ',
     &    'hycom_icefreeday ext.a idm jdm day1 dayi free.a'
        CALL EXIT(1)
      ENDIF
C
      NPAD = 4096 - MOD(IDM*JDM,4096)
      IF     (NPAD.EQ.4096) THEN
        NPAD = 0
      ENDIF
C
      ALLOCATE( EXT(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_icefreeday: could not allocate ',
     +             IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( FREE(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_icefreeday: could not allocate ',
     +             IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
C
      CALL FREEDAY(EXT,FREE,IDM,JDM,PAD,NPAD,
     &             DAY1,DAYI, CFILEI,CFILEO)
      CALL EXIT(0)
      END
      SUBROUTINE FREEDAY(EXT,FREE,IDM,JDM,PAD,NPAD,
     &                   DAY1,DAYI, CFILEI,CFILEO)
      IMPLICIT NONE
C
      REAL*4     SPVAL
      PARAMETER (SPVAL=2.0**100)
C
      CHARACTER*240 CFILEI,CFILEO
      INTEGER       IDM,JDM,NPAD
      REAL          DAY1,DAYI
      REAL*4        EXT(IDM,JDM),FREE(IDM,JDM),PAD(NPAD)
C
C     MOST OF WORK IS DONE HERE.
C
#ifdef sun
      INTEGER      IR_ISNAN
C
#endif
      CHARACTER*18 CASN
      INTEGER      I,J,IOS,IR,NRECL
      REAL*4       AMN,AMX
#ifdef CRAY
      INTEGER*8    IU8,IOS8
#endif
C
      IF     (NPAD.EQ.0) THEN
        INQUIRE(IOLENGTH=NRECL) EXT
      ELSE
        INQUIRE(IOLENGTH=NRECL) EXT,PAD
        PAD(:) = SPVAL
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(PAD,NPAD)
#endif
      ENDIF
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
      CALL ASNUNIT(21,'-F syscall -N ieee',IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t asnunit 21'
        write(6,*) 'ios = ',ios
        CALL EXIT(5)
      ENDIF
#endif
#endif
      OPEN(UNIT=11, FILE=CFILEI, FORM='UNFORMATTED', STATUS='OLD',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILEI)
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
      OPEN(UNIT=21, FILE=CFILEO, FORM='UNFORMATTED', STATUS='NEW',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILEO)
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
C
      DO J= 1,JDM
        DO I= 1,IDM
          FREE(I,J) = SPVAL
        ENDDO
      ENDDO
C
      DO IR= 1,999999
        READ(11,REC=IR,IOSTAT=IOS) EXT
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(EXT,IDM*JDM)
#endif
        IF     (IOS.NE.0) THEN
          IF     (IR.EQ.1) THEN
            WRITE(6,*) 'can''t read ',TRIM(CFILEI)
            CALL EXIT(4)
          ELSE
            EXIT
          ENDIF
        ENDIF
        DO J= 1,JDM
          DO I= 1,IDM
            IF     (EXT(I,J).EQ.0.0 .AND. FREE(I,J).EQ.SPVAL) THEN
              FREE(I,J) = DAY1 + (IR-1)*DAYI
            ENDIF
          ENDDO !i
        ENDDO !j
      ENDDO !ir
C
      AMN =  SPVAL
      AMX = -SPVAL
      DO J= 1,JDM
        DO I= 1,IDM
          IF     (FREE(I,J).NE.SPVAL) THEN
            AMN = MIN( AMN, FREE(I,J) )
            AMX = MAX( AMX, FREE(I,J) )
          ENDIF
        ENDDO
      ENDDO
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(FREE,IDM*JDM)
#endif
      IF     (NPAD.EQ.0) THEN
        WRITE(21,REC=1,IOSTAT=IOS) FREE
      ELSE
        WRITE(21,REC=1,IOSTAT=IOS) FREE,PAD
      ENDIF
      WRITE(6,'(a,2f10.4)') 'ice free day min, max = ',AMN,AMX
      WRITE(6,*)
      WRITE(6,*) IR-1,' FIELDS PROCESSED'
      WRITE(6,*)
C
      CLOSE(UNIT=11)
      CLOSE(UNIT=21)
C
      RETURN
      END

      PROGRAM HYCOM_PERTURBATION_SCALE
      IMPLICIT NONE
C
C  hycom_perturbation_scale - Usage:  hycom_perturbation_scale hscl fout.a
C
C                 creates a HYCOM .[ab] file containing scale
C                 factor for hycom_perturbation
C
C                 a positive hscl is the gaussian length scale in m
C                 a negative hscl is minus the gaussian length scale in deg
C                 when converting from degrees to m, the minimum length
C                 is clipped half the maximum length scale
C
C                 idm,jdm are taken from grid.a
C
C  the grid is assumed to be p-grid global with an arctic bi-polar patch.
C  always use mask.a for closed domains.
C
C  fout.a will contain idm*jdm 32-bit IEEE real values for the array,
C   in standard f77 element order, followed by padding to a multiple
C   of 4096 32-bit words, but otherwise with no control bytes/words,
C   and values of 2.0**100 indicating a data void.
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  March 2013.
C  Based on a program by Matthew J. Carrier.
C
      REAL*4,  ALLOCATABLE :: HSCL(:,:)
      REAL*4               :: PAD(4096)
      INTEGER       IOS
      INTEGER       IARGC
      INTEGER       NARG
      CHARACTER*240 CARG
C
      LOGICAL       LDEBUG
      INTEGER       IDM,JDM,NPAD,NSAMPLE,NFLD
      REAL*4        HHSCL
      CHARACTER*6   CVARIN
      CHARACTER*240 CFILEG,CFILEM,CFILEB,CFILEO
C
C     READ ARGUMENTS.
C
      CALL GETARG(0,CARG)
      LDEBUG = .FALSE.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.2) THEN
        CALL GETARG(1,CARG)
        READ(CARG,*)    HHSCL
        CALL GETARG(2,CFILEO)
      ELSE
        WRITE(6,*) 'Usage:  hycom_perturbation_scale '//
     +    'hscl fout.a'
        CALL EXIT(1)
      ENDIF
C
      CFILEG = 'regional.grid.a'
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
        WRITE(6,*) 'hycom_perturbation_scale: bad header file ',
     &             TRIM(CFILEB)
        CALL EXIT(2)
      ENDIF
      READ( 11,*) JDM,CVARIN
      IF (CVARIN.NE.'jdm   ') THEN
        WRITE(6,*) 'hycom_perturbation_scale: bad header file ',
     &             TRIM(CFILEB)
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
      ALLOCATE( HSCL(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*)
     +   'Error in hycom_perturbation_scale: could not allocate ',
     +   IDM*JDM,' words for HSCL'
        CALL EXIT(2)
      ENDIF
C
      CALL PERTURBATION(HSCL,IDM,JDM,PAD,NPAD, HHSCL, CFILEG,CFILEO,
     &                  LDEBUG)
      CALL EXIT(0)
 5000 FORMAT(I4)
      END
      SUBROUTINE PERTURBATION(HSCL,IDM,JDM, PAD,NPAD,
     +                        HHSCL, CFILEG,CFILEO,
     +                        LDEBUG)
      IMPLICIT NONE
C
      REAL*4     SPVAL
      PARAMETER (SPVAL=2.0**100)
C
      CHARACTER*240 CFILE,CFILEG,CFILEO
      LOGICAL       LDEBUG
      INTEGER       IDM,JDM,NPAD
      REAL*4        HSCL(IDM,JDM),PAD(NPAD),
     +              HHSCL,DEG2RAD
C
C     MOST OF WORK IS DONE HERE.
C
#ifdef sun
      INTEGER      IR_ISNAN
C
#endif
      CHARACTER*18 CASN
      INTEGER      I,II,IOS,J,K,KS,NRECL
      REAL         SCL,V,X,Y
      REAL*4       AMX,AMN
#ifdef CRAY
      INTEGER*8    IU8,IOS8
#endif
C
C     INPUT GRID ARRAYS.
C
      INQUIRE( IOLENGTH=NRECL) HSCL,PAD
      OPEN(UNIT=11, FILE=CFILEG, FORM='UNFORMATTED', STATUS='OLD',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILEG)
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
C
      IF     (HHSCL.GE.0.0) THEN
        HSCL(:,:) = HHSCL
        AMN = HHSCL
        AMX = HHSCL
      ELSE
        READ(11,REC=2,IOSTAT=IOS) HSCL  ! plat
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(HSCL,IDM*JDM)
#endif
        IF     (IOS.NE.0) THEN
          WRITE(6,*) 'can''t read record 2 of ',TRIM(CFILEG)
          CALL EXIT(4)
        ENDIF
        DEG2RAD = 4.D0*ATAN(1.D0)/180.D0  !PI/180
        AMN =  SPVAL
        AMX = -SPVAL
        DO J= 1,JDM
          DO I= 1,IDM
            HSCL(I,J) = -HHSCL * 111.2E3 * 
     &                    SQRT(ABS(COS(HSCL(I,J)*DEG2RAD)))
            AMX = MAX( AMX, HSCL(I,J) )
            AMN = MIN( AMN, HSCL(I,J) )
          ENDDO !i
        ENDDO !j
        AMN = MAX( AMN, 0.5*AMX )  !clip any small values
        DO J= 1,JDM
          DO I= 1,IDM
            HSCL(I,J) = MAX( HSCL(I,J), AMN )
          ENDDO !i
        ENDDO !j
      ENDIF !hhscl>0:else
        write(6,*) 'hscl = ',AMN,AMX
        call flush(6)
C
      CLOSE(UNIT=11)
C
C     OUTPUT FILE.
C
      OPEN(UNIT=21, FILE=CFILEO, FORM='UNFORMATTED', STATUS='NEW',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILEO)
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
C
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(R2D,IDM*JDM)
#endif
        WRITE(21,REC=1,IOSTAT=IOS) HSCL
        WRITE(6,'(a,1p2g16.8)')
     &     'min, max = ',AMN,AMX
        CALL FLUSH(6)
      WRITE(6,*)
      WRITE(6,*) 1,' FIELDS PROCESSED'
      WRITE(6,*)
C
      CLOSE(21)
C
      RETURN
      END

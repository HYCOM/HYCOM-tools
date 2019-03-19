      PROGRAM HYCOM_BOTFRIC
      IMPLICIT NONE
C
C hycom_botfric - Usage: hycom_botfric depth.a idm jdm z0 cbmin cb.a
C
C                 Calculates roughness-based bottom friction coefficient
C
C  depth.a is assumed to contain idm*jdm 32-bit IEEE real values
C   for each array, in standard f77 element order, followed by padding
C   to a multiple of 4096 32-bit words, but otherwise with no control
C   bytes/words, and input values of 2.0**100 indicating a data void.
C
C  depth.a contains the bottom depth in m
C  z0      is the bottom roughness in m
C  cbmin   is the minimum quadratic bottom friction coeficient
C  cb.a    will contain the roughness-based bottom friction coefficient
C
C  Uses NCOM's formula for 2-D cases:
C    cb(i,j) = max( cbmin, (vonk/log(0.5*depth(i,j)/z0))**2 )
C
C  The same scheme, in 3-D, is discussed in Appendix 2 of 
C   Oey, L.-Y., 2006:  An OGCM with movable land-sea boundaries.
C   Ocean Modeling, 2006, 13, 176-195. 
C  It is using the "Law of the Wall", see:
C   Schlichting, H., 1968. Boundary Layer Theory. McGraw-Hill, 748 pp. 
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  December 2007.
C
      REAL*4, ALLOCATABLE :: DEPTH(:,:),CBP(:,:)
      REAL*4              :: PAD(4096)
      INTEGER IOS
      INTEGER      IARGC
      INTEGER      NARG
      CHARACTER*240 CARG
C
      INTEGER       IDM,JDM,NPAD
      REAL*4        Z0,CBMIN
      CHARACTER*240 CFILE1,CFILEO
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.6) THEN
        CALL GETARG(1,CFILE1)
        CALL GETARG(2,CARG)
        READ(CARG,*) IDM
        CALL GETARG(3,CARG)
        READ(CARG,*) JDM
        CALL GETARG(4,CARG)
        READ(CARG,*) Z0
        CALL GETARG(5,CARG)
        READ(CARG,*) CBMIN
        CALL GETARG(6,CFILEO)
      ELSE
        WRITE(6,*)
     &  'Usage:hycom_botfric depth.a idm jdm z0 cbmin cb.a'
        CALL EXIT(1)
      ENDIF
C
      NPAD = 4096 - MOD(IDM*JDM,4096)
      IF     (NPAD.EQ.4096) THEN
        NPAD = 0
      ENDIF
*
*     WRITE(6,'(a,1pg16.8)')
*    &   'THBASE = ',THBASE
*     WRITE(6,'(a,i6)')
*    &   'NPAD   = ',NPAD   
C
      ALLOCATE(DEPTH(IDM,JDM), 
     &           CBP(IDM,JDM), 
     &         STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_expr: could not allocate ',
     +             2*IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
C
      CALL BF(DEPTH,CBP,IDM,JDM,PAD,NPAD, 
     &        Z0,CBMIN, CFILE1,CFILEO)
      CALL EXIT(0)
      END
      SUBROUTINE BF(DEPTH,CBP,IDM,JDM,PAD,NPAD, 
     &              Z0,CBMIN, CFILE1,CFILEO)
      IMPLICIT NONE
C
      REAL*4     SPVAL
      PARAMETER (SPVAL=2.0**100)
C
      CHARACTER*240 CFILE1,CFILEO
      INTEGER      IDM,JDM,NPAD
      REAL*4       Z0,CBMIN
      REAL*4       DEPTH(IDM,JDM),
     &               CBP(IDM,JDM), PAD(NPAD)
C
C     MOST OF WORK IS DONE HERE.
C
      CHARACTER*18 CASN
      INTEGER      I,J,IOS,NRECL
      REAL*4       VONK,AMN,AMX
#ifdef CRAY
      INTEGER*8    IU8,IOS8
#endif
C
      PAD(:) = SPVAL
      INQUIRE( IOLENGTH=NRECL) CBP,PAD
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
      OPEN(UNIT=11, FILE=CFILE1, FORM='UNFORMATTED', STATUS='OLD',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILE1)
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
c
      READ(11,REC=1,IOSTAT=IOS) DEPTH
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(DEPTH,IDM*JDM)
#endif
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'can''t read DEPTH from ',TRIM(CFILE1)
        CALL EXIT(4)
      ENDIF
C
      CLOSE(11)
*
*     CALL MINMAX(DEPTH,IDM,JDM, AMN,AMX)
*     WRITE(6,'(a,1p2g16.8)')
*    &   'DEPTH - min, max = ',AMN,AMX
C
      VONK = 0.4
      AMN  =  SPVAL
      AMX  = -SPVAL
      DO J= 1,JDM
        DO I= 1,IDM
          IF     (DEPTH(I,J).NE.SPVAL) THEN
            CBP(I,J) = MAX( CBMIN,
     &                      VONK**2/LOG(0.5*DEPTH(I,J)/Z0)**2 )
            AMN = MIN( AMN, CBP(I,J) )
            AMX = MAX( AMX, CBP(I,J) )
          ELSE
            CBP(I,J) = SPVAL
          ENDIF
        ENDDO
      ENDDO
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(CBP,IDM*JDM)
#endif
      WRITE(21,REC=1,IOSTAT=IOS) CBP,PAD
      WRITE(6,'(a,1p2g16.8)')
     &     'CB - min, max = ',AMN,AMX
      CLOSE(21)
C
      RETURN
      END
      SUBROUTINE MINMAX(ATH,IDM,JDM, AMN,AMX)
      IMPLICIT NONE
C
      REAL*4     SPVAL
      PARAMETER (SPVAL=2.0**100)
C
      INTEGER      IDM,JDM
      REAL*4       ATH(IDM,JDM),AMN,AMX
C
      INTEGER I,J
C
      AMN =  SPVAL
      AMX = -SPVAL
      DO J= 1,JDM
        DO I= 1,IDM
          IF     (ATH(I,J).NE.SPVAL) THEN
            AMN = MIN( AMN, ATH(I,J) )
            AMX = MAX( AMX, ATH(I,J) )
          ENDIF
        ENDDO
      ENDDO
      RETURN
      END

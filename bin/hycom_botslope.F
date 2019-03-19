      PROGRAM HYCOM_BOTSLOPE
      IMPLICIT NONE
C
C hycom_botslope - Usage: hycom_botslope depth.a idm jdm slope.a
C
C                 Calculates the bottom slope
C
C  depth.a is assumed to contain idm*jdm 32-bit IEEE real values
C   for each array, in standard f77 element order, followed by padding
C   to a multiple of 4096 32-bit words, but otherwise with no control
C   bytes/words, and input values of 2.0**100 indicating a data void.
C
C  depth.a contains the bottom depth in m
C  slope.a will contain the bottom slope
C
C  Returns the maximum value of N,S,E,W slopes (all positive).
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  April 2003.
C
      REAL*4, ALLOCATABLE :: DEPTH(:,:),SLOPE(:,:)
      REAL*4              :: PAD(4096)
      INTEGER       IOS
      INTEGER       IARGC
      INTEGER       NARG
      CHARACTER*240 CARG
C
      INTEGER       IDM,JDM,NPAD
      CHARACTER*240 CFILE1,CFILEO
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.4) THEN
        CALL GETARG(1,CFILE1)
        CALL GETARG(2,CARG)
        READ(CARG,*) IDM
        CALL GETARG(3,CARG)
        READ(CARG,*) JDM
        CALL GETARG(4,CFILEO)
      ELSE
        WRITE(6,*)
     &  'Usage:hycom_botslope depth.a idm jdm slope.a'
        CALL EXIT(1)
      ENDIF
C
      NPAD = 4096 - MOD(IDM*JDM,4096)
      IF     (NPAD.EQ.4096) THEN
        NPAD = 0
      ENDIF
*
*     WRITE(6,'(a,i6)')
*    &   'NPAD   = ',NPAD   
C
      ALLOCATE(DEPTH(IDM,JDM), 
     &         SLOPE(IDM,JDM), 
     &         STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_expr: could not allocate ',
     +             2*IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
C
      CALL BS(DEPTH,SLOPE,IDM,JDM,PAD,NPAD, CFILE1,CFILEO)
      CALL EXIT(0)
      END
      SUBROUTINE BS(DEPTH,SLOPE,IDM,JDM,PAD,NPAD, CFILE1,CFILEO)
      IMPLICIT NONE
C
      REAL*4     SPVAL
      PARAMETER (SPVAL=2.0**100)
C
      CHARACTER*240 CFILE1,CFILEO
      INTEGER      IDM,JDM,NPAD
      REAL*4       DEPTH(IDM,JDM),
     &             SLOPE(IDM,JDM), PAD(NPAD)
C
C     MOST OF WORK IS DONE HERE.
C
      CHARACTER*18 CASN
      INTEGER      I,II,J,JJ,IOS,NRECL
      REAL*4       SN,SS,SE,SW,AMN,AMX
#ifdef CRAY
      INTEGER*8    IU8,IOS8
#endif
C
      PAD(:) = SPVAL
      INQUIRE( IOLENGTH=NRECL) SLOPE,PAD
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
      AMN  =  SPVAL
      AMX  = -SPVAL
      DO J= 1,JDM
        DO I= 1,IDM
          IF     (DEPTH(I,J).NE.SPVAL) THEN
            II = MIN(I+1,IDM)
            JJ = J
            IF     (DEPTH(II,JJ).NE.SPVAL) THEN
              SN = ABS( DEPTH(I,J) - DEPTH(II,JJ) ) /
     &                ( DEPTH(I,J) + DEPTH(II,JJ) )
            ELSE
              SN = 0.0
            ENDIF
            II = MAX(I-1,1)
            JJ = J
            IF     (DEPTH(II,JJ).NE.SPVAL) THEN
              SS = ABS( DEPTH(I,J) - DEPTH(II,JJ) ) /
     &                ( DEPTH(I,J) + DEPTH(II,JJ) )
            ELSE
              SS = 0.0
            ENDIF
            II = I
            JJ = MIN(J+1,JDM)
            IF     (DEPTH(II,JJ).NE.SPVAL) THEN
              SE = ABS( DEPTH(I,J) - DEPTH(II,JJ) ) /
     &                ( DEPTH(I,J) + DEPTH(II,JJ) )
            ELSE
              SE = 0.0
            ENDIF
            II = I
            JJ = MAX(J-1,1)
            IF     (DEPTH(II,JJ).NE.SPVAL) THEN
              SW = ABS( DEPTH(I,J) - DEPTH(II,JJ) ) /
     &                ( DEPTH(I,J) + DEPTH(II,JJ) )
            ELSE
              SW = 0.0
            ENDIF
C
            SLOPE(I,J) = MAX( SN, SS, SE, SW )
            AMN = MIN( AMN, SLOPE(I,J) )
            AMX = MAX( AMX, SLOPE(I,J) )
          ELSE
            SLOPE(I,J) = SPVAL
          ENDIF
        ENDDO
      ENDDO
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(SLOPE,IDM*JDM)
#endif
      WRITE(21,REC=1,IOSTAT=IOS) SLOPE,PAD
      WRITE(6,'(a,1p2g16.8)')
     &     'slope - min, max = ',AMN,AMX
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

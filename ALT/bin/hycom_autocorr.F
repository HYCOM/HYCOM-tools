      PROGRAM HYCOM_AUTOCORR
      IMPLICIT NONE
C
C  hycom_autocorr - Usage:  hycom_autocorr fin.a idm jdm itlrec increc numrec lagmax fauto.a [lperiod]
C
C                 Outputs the autocorrelations at time lags from 1 to lagmax
C
C  see also hycom_autocorr_lag
C
C  fin*.a is assumed to contain idm*jdm 32-bit IEEE real values for
C   each array, in standard f77 element order, followed by padding
C   to a multiple of 4096 32-bit words, but otherwise with no control
C   bytes/words, and input values of 2.0**100 indicating a data void.
C
C  With lperiod set, record itlrec+(numrec)*increc is record itlrec
C      and the highest record read is itlrec+(numrec-1)*increc.
C  Otherwise, the last record read is itlrec+(numrec-1+lagmax)*increc.
C
C  This version only needs five idm*jdm arrays, but this is at the cost
C   of being very I/O intensive.
C
C  This version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  February 2001.
C
      REAL*4, ALLOCATABLE :: A(:,:),AM(:,:),AV(:,:),AL(:,:),
     &                       ALAG(:,:)
      REAL*4              :: PAD(4096)
      INTEGER       IOS,L
      INTEGER       IARGC
      INTEGER       NARG
      CHARACTER*240 CARG
C
      LOGICAL       LPERIOD
      INTEGER       IDM,JDM,ITLREC,INCREC,LAGMAX,NUMREC,NPAD
      CHARACTER*240 CFILE1,CFILEO
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.8 .OR. NARG.EQ.9) THEN
        CALL GETARG(1,CFILE1)
        CALL GETARG(2,CARG)
        READ(CARG,*) IDM
        CALL GETARG(3,CARG)
        READ(CARG,*) JDM
        CALL GETARG(4,CARG)
        READ(CARG,*) ITLREC
        CALL GETARG(5,CARG)
        READ(CARG,*) INCREC
        CALL GETARG(6,CARG)
        READ(CARG,*) NUMREC
        CALL GETARG(7,CARG)
        READ(CARG,*) LAGMAX
        CALL GETARG(8,CFILEO)
        LPERIOD = NARG.EQ.9
      ELSE
        WRITE(6,'(3a)')
     &    'Usage:  ',
     &    'hycom_autocorr ',
     &    'fin.a idm jdm itlrec increc numrec lagmax fauto.a [lperiod]'
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
        WRITE(6,*) 'Error in hycom_autocorr: could not allocate ',
     +             IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( AM(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_autocorr: could not allocate ',
     +             IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( AV(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_autocorr: could not allocate ',
     +             IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( AL(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_autocorr: could not allocate ',
     +             IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( ALAG(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_autocorr: could not allocate ',
     +             IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
C
      CALL AUTO(A,AM,AV,AL,ALAG,IDM,JDM,LAGMAX,PAD,NPAD,
     &          ITLREC,INCREC,NUMREC, CFILE1,CFILEO, LPERIOD)
      CALL EXIT(0)
      END
      SUBROUTINE AUTO(A,AM,AV,AL,ALAG,IDM,JDM,LAGMAX,PAD,NPAD,
     &                ITLREC,INCREC,NUMREC, CFILE1,CFILEO, LPERIOD)
      IMPLICIT NONE
C
      REAL*4     SPVAL
      PARAMETER (SPVAL=2.0**100)
C
      CHARACTER*240 CFILE1,CFILEO
      LOGICAL       LPERIOD
      INTEGER       IDM,JDM,LAGMAX,NPAD,ITLREC,INCREC,NUMREC
      REAL*4        A(IDM,JDM),AM(IDM,JDM),AV(IDM,JDM),AL(IDM,JDM),
     &              ALAG(IDM,JDM),PAD(NPAD)
C
C     MOST OF WORK IS DONE HERE.
C
#ifdef sun
      INTEGER      IR_ISNAN
C
#endif
      CHARACTER*18 CASN
      INTEGER      I,J,IOS,IOUT,IR,IRL,L,NR,NRL,NRECL,NUMR
      REAL*4       AMN,AMX,RNUMR
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
C
      DO J= 1,JDM
        DO I= 1,IDM
          AM(I,J) = 0.0
          AV(I,J) = 0.0
        ENDDO !i
      ENDDO !j
C
      DO NR= 1,NUMREC
        IR = ITLREC + INCREC*(NR-1)
        READ(11,REC=IR,IOSTAT=IOS) A
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
        IF     (IOS.NE.0) THEN
          WRITE(6,*) 'can''t read record ',IR,' of ',
     &               TRIM(CFILE1)
          CALL EXIT(4)
        ENDIF
        DO J= 1,JDM
          DO I= 1,IDM
#ifdef sun
            IF     (IR_ISNAN(A(I,J)).NE.1) THEN
              IF     (A(I,J).NE.SPVAL) THEN
                AM(I,J) = AM(I,J) + A(I,J)
              ENDIF
            ENDIF
#else
            IF     (A(I,J).NE.SPVAL) THEN
              AM(I,J) = AM(I,J) + A(I,J)
            ENDIF
#endif
          ENDDO !i
        ENDDO !j
      ENDDO !nr
C
      RNUMR = 1.0/NUMREC
      AMN =  SPVAL
      AMX = -SPVAL
      DO J= 1,JDM
        DO I= 1,IDM
          IF     (A(I,J).NE.SPVAL) THEN
            AM(I,J) = AM(I,J)*RNUMR
            AMN = MIN( AMN, AM(I,J) )
            AMX = MAX( AMX, AM(I,J) )
          ELSE
            AM(I,J) = SPVAL
          ENDIF
        ENDDO !i
      ENDDO !j
C
      DO L= 1,LAGMAX
        DO J= 1,JDM
          DO I= 1,IDM
            ALAG(I,J) = 0.0
          ENDDO !i
        ENDDO !j
        DO NR= 1,NUMREC
          IR = ITLREC + INCREC*(NR-1)
          READ(11,REC=IR,IOSTAT=IOS) A
#ifdef ENDIAN_IO
          CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
          IF     (IOS.NE.0) THEN
            WRITE(6,*) 'can''t read record ',IR,' of ',
     &                 TRIM(CFILE1)
            CALL EXIT(4)
          ENDIF
          IF      (L.EQ.1) THEN
            DO J= 1,JDM
              DO I= 1,IDM
#ifdef sun
                IF     (IR_ISNAN(A(I,J)).NE.1) THEN
                  IF     (A(I,J).NE.SPVAL) THEN
                    AV(I,J) = AV(I,J) + (A(I,J)-AM(I,J))**2
                  ENDIF
                ENDIF
#else
                IF     (A(I,J).NE.SPVAL) THEN
                  AV(I,J) = AV(I,J) + (A(I,J)-AM(I,J))**2
                ENDIF
#endif
              ENDDO !i
            ENDDO !j
          ENDIF !l==1
          IF     (LPERIOD) THEN
            NRL = MOD(NR+L-1,NUMREC)+1
          ELSE
            NRL =     NR+L
          ENDIF
          IRL = ITLREC + INCREC*(NRL-1)
          READ(11,REC=IRL,IOSTAT=IOS) AL
#ifdef ENDIAN_IO
          CALL ENDIAN_SWAP(AL,IDM*JDM)
#endif
          IF     (IOS.NE.0) THEN
            WRITE(6,*) 'can''t read record ',IRL,' of ',
     &                 TRIM(CFILE1)
            CALL EXIT(4)
          ENDIF
          DO J= 1,JDM
            DO I= 1,IDM
#ifdef sun
              IF     (IR_ISNAN(A(I,J)).NE.1) THEN
                IF     (A(I,J).NE.SPVAL) THEN
                  ALAG(I,J) = ALAG(I,J) + (A( I,J)-AM(I,J))*
     &                                    (AL(I,J)-AM(I,J))
                ENDIF
              ENDIF
#else
              IF     (A(I,J).NE.SPVAL) THEN
                ALAG(I,J) = ALAG(I,J) + (A( I,J)-AM(I,J))*
     &                                  (AL(I,J)-AM(I,J))
              ENDIF
#endif
            ENDDO !i
          ENDDO !j
        ENDDO !nr
C
        AMN =  SPVAL
        AMX = -SPVAL
        DO J= 1,JDM
          DO I= 1,IDM
            IF     (A(I,J).NE.SPVAL) THEN
              IF     (AV(I,J).NE.0.0) THEN
                ALAG(I,J) = ALAG(I,J)/AV(I,J)
              ELSE
                ALAG(I,J) = 0.0
              ENDIF
              AMN = MIN( AMN, ALAG(I,J) )
              AMX = MAX( AMX, ALAG(I,J) )
            ENDIF
          ENDDO !i
        ENDDO !j
        WRITE(21,REC=L,IOSTAT=IOS) ALAG
        WRITE(6,'(a,1p2g16.8)') 'ALAG: min, max = ',AMN,AMX
      ENDDO !l
C
      CLOSE(UNIT=11)
      CLOSE(UNIT=21)
C
      RETURN
      END

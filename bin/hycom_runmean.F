      PROGRAM FMEAN
      IMPLICIT NONE
C
C  hycom_runmean      - Usage:  hycom_runmean      fin.a idm jdm numsum itlrec increc numout fout.a
C  hycom_runmean_half - Usage:  hycom_runmean_half fin.a idm jdm numsum itlrec increc numout fout.a
C
C                 Outputs means of numsum consecutive fields.
C                 The n-th mean starts with input field itlrec+(n-1)*increc.
C                 There are numout means produced, so the input file must
C                 contain at least itlrec+(numout)*increc+numsum-1 fields.
C                 If numout=0 the maximum number of means is produced.
C                 For hycom_runmean_half the 1st and last fields in the
C                  mean are multiplied by 0.5.
C
C  One typical usage is to have increc=numsum, e.g. to produce a daily
C  mean from 6-hrly fields (increc=numsum=4).
C
C  A special case that is treated differently is: increc=1 and numout
C  set to the number of records in fin.a.  Then the data set is assumed
C  to be periodic, i.e. record n maps to record mod(n-1+numout,numout)+1
C  which is between 1 and numout.  In this case, itlrec must be 
C  1-(numsum+1)/2 (<1) to produce the running mean of the entire data set.
C
C  fin*.a is assumed to contain idm*jdm 32-bit IEEE real values for
C   each array, in standard f77 element order, followed by padding
C   to a multiple of 4096 32-bit words, but otherwise with no control
C   bytes/words, and input values of 2.0**100 indicating a data void.
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  February 2001.
C
      REAL*4, ALLOCATABLE :: A(:,:),AM(:,:)
      REAL*4              :: PAD(4096)
      INTEGER       IOS,L
      INTEGER       IARGC
      INTEGER       NARG
      CHARACTER*240 CARG
C
      LOGICAL       LHALF
      INTEGER       IDM,JDM,NUMSUM,ITLREC,INCREC,NUMOUT,NPAD
      CHARACTER*240 CFILE1,CFILEO
C
C     READ ARGUMENTS.
C
      CALL GETARG(0,CARG)
      L = LEN_TRIM(CARG)
      LHALF = CARG(L-4:L).EQ.'_half'
C
      NARG = IARGC()
C
      IF     (NARG.EQ.8) THEN
        CALL GETARG(1,CFILE1)
        CALL GETARG(2,CARG)
        READ(CARG,*) IDM
        CALL GETARG(3,CARG)
        READ(CARG,*) JDM
        CALL GETARG(4,CARG)
        READ(CARG,*) NUMSUM
        CALL GETARG(5,CARG)
        READ(CARG,*) ITLREC
        CALL GETARG(6,CARG)
        READ(CARG,*) INCREC
        CALL GETARG(7,CARG)
        READ(CARG,*) NUMOUT
        CALL GETARG(8,CFILEO)
      ELSEIF (LHALF) THEN
        WRITE(6,'(3a)')
     &    'Usage:  ',
     &    'hycom_runmean_half ',
     &    'fin.a idm jdm numsum itlrec increc numout fout.a'
        CALL EXIT(1)
      ELSE
        WRITE(6,'(3a)')
     &    'Usage:  ',
     &    'hycom_runmean ',
     &    'fin.a idm jdm numsum itlrec increc numout fout.a'
        CALL EXIT(1)
      ENDIF
C
      IF     (ITLREC.LT.1) THEN  !special case, periodic input
        IF     (ITLREC.NE.1-(numsum+1)/2 .OR.
     &          INCREC.NE.1 .OR. NUMOUT.LE.0) THEN
          IF (LHALF) THEN
            WRITE(6,'(3a)')
     &        'Usage (periodic):  ',
     &        'hycom_runmean_half ',
     &        'fin.a idm jdm numsum 1-(numsum+1)/2 1 numrec fout.a'
            CALL EXIT(0)
          ELSE
            WRITE(6,'(3a)')
     &        'Usage (periodic):  ',
     &        'hycom_runmean ',
     &        'fin.a idm jdm numsum 1-(numsum+1)/2 1 numrec fout.a'
            CALL EXIT(0)
          ENDIF
        ENDIF
      ENDIF
C
      NPAD = 4096 - MOD(IDM*JDM,4096)
      IF     (NPAD.EQ.4096) THEN
        NPAD = 0
      ENDIF
C
      ALLOCATE( A(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_runmean: could not allocate ',
     +             IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( AM(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_runmean: could not allocate ',
     +             IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
C
      CALL MEAN(A,AM,IDM,JDM,PAD,NPAD,
     &          NUMSUM,ITLREC,INCREC,NUMOUT, LHALF, CFILE1,CFILEO)
      CALL EXIT(0)
      END
      SUBROUTINE MEAN(A,AM,IDM,JDM,PAD,NPAD,
     &                NUMSUM,ITLREC,INCREC,NUMOUT, LHALF, CFILE1,CFILEO)
      IMPLICIT NONE
C
      REAL*4     SPVAL
      PARAMETER (SPVAL=2.0**100)
C
      CHARACTER*240 CFILE1,CFILEO
      LOGICAL      LHALF
      INTEGER      IDM,JDM,NPAD,NUMSUM,ITLREC,INCREC,NUMOUT
      REAL*4       A(IDM,JDM),AM(IDM,JDM),PAD(NPAD)
C
C     MOST OF WORK IS DONE HERE.
C
#ifdef sun
      INTEGER      IR_ISNAN
C
#endif
      CHARACTER*18 CASN
      INTEGER      LEN_TRIM
      INTEGER      I,J,IOS,IOUT,IR,NR,NRECL,NUMR
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
        write(6,*) 'Error: can''t open ',CFILE1(1:LEN_TRIM(CFILE1))
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
      IF     (NUMOUT.EQ.0) THEN
        NUMOUT = 99999
      ENDIF
C
      DO 110 IOUT= 1,NUMOUT
        DO J= 1,JDM
          DO I= 1,IDM
            AM(I,J) =  0.0
          ENDDO
        ENDDO
C
        DO NR= 1,NUMSUM
          IR = ITLREC + INCREC*(IOUT-1) + NR-1
          IF     (ITLREC.LT.1) THEN
            IR = MOD(IR-1+NUMOUT,NUMOUT)+1
          ENDIF
          READ(11,REC=IR,IOSTAT=IOS) A
#ifdef ENDIAN_IO
          CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
          IF     (IOS.NE.0) THEN
            IF     (NUMOUT.EQ.99999 .AND. IOUT.NE.1) THEN
              NUMOUT = IOUT-1
              GOTO 1110
            ELSE
              WRITE(6,*) 'can''t read record ',IR,' of ',
     &                   CFILE1(1:LEN_TRIM(CFILE1))
              CALL EXIT(4)
            ENDIF
          ENDIF
          IF     (LHALF .AND. (NR.EQ.1 .OR. NR.EQ.NUMSUM)) THEN
            DO J= 1,JDM
              DO I= 1,IDM
#ifdef sun
                IF     (IR_ISNAN(A(I,J)).NE.1) THEN
                  IF     (A(I,J).NE.SPVAL) THEN
                    AM(I,J) = AM(I,J) + 0.5*A(I,J)
                  ENDIF
                ENDIF
#else
                IF     (A(I,J).NE.SPVAL) THEN
                  AM(I,J) = AM(I,J) + 0.5*A(I,J)
                ENDIF
#endif
              ENDDO !i
            ENDDO !j
          ELSE
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
          ENDIF !lhalf
        ENDDO !nr
C
        IF     (LHALF) THEN
          RNUMR = 1.0/(NUMSUM-1)
        ELSE
          RNUMR = 1.0/NUMSUM
        ENDIF
        AMN =  SPVAL
        AMX = -SPVAL
        DO J= 1,JDM
          DO I= 1,IDM
            IF     (A(I,J).NE.SPVAL) THEN
              A(I,J) = AM(I,J)*RNUMR
              AMN = MIN( AMN, A(I,J) )
              AMX = MAX( AMX, A(I,J) )
            ELSE
              A(I,J) = SPVAL
            ENDIF
          ENDDO
        ENDDO
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
        WRITE(21,REC=IOUT,IOSTAT=IOS) A
        WRITE(6,'(a,1p2g16.8)') 'min, max = ',AMN,AMX
  110 CONTINUE
 1110 CONTINUE
      IF     (LHALF) THEN
        WRITE(6,*)
        WRITE(6,*) NUMOUT,' END-WEIGHTED MEANS PRODUCED'
        WRITE(6,*)
      ELSE
        WRITE(6,*)
        WRITE(6,*) NUMOUT,' MEANS PRODUCED'
        WRITE(6,*)
      ENDIF
C
      CLOSE(UNIT=11)
      CLOSE(UNIT=21)
C
      RETURN
      END

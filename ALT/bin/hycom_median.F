      PROGRAM HYCOM_MEDIAN
      IMPLICIT NONE
C
C  hycom_median - Usage:  hycom_median fin.a idm jdm itlrec increc numrec fout.a
C
C                 Outputs a single (1:idm,1:jdm) field, representing the
C                 median of input fields itl+(n-1)*inc for n=1:numrec.
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
      REAL*4, ALLOCATABLE :: A(:,:,:),AM(:,:)
      REAL*4              :: PAD(4096)
      INTEGER      IOS,L
      INTEGER      IARGC
      INTEGER      NARG
      CHARACTER*240 CARG
C
      INTEGER       IDM,JDM,ITLREC,INCREC,NUMREC,ITEST,JTEST,NPAD
      CHARACTER*240 CFILE1,CFILEO
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.7) THEN
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
        CALL GETARG(7,CFILEO)
        ITEST = 0
        JTEST = 0
      ELSEIF (NARG.EQ.9) THEN  !undocumented debug option
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
        CALL GETARG(7,CFILEO)
        CALL GETARG(8,CARG)
        READ(CARG,*) ITEST
        CALL GETARG(9,CARG)
        READ(CARG,*) JTEST
      ELSE
        WRITE(6,'(2a)')
     &    'Usage:  ',
     &    'hycom_median fin.a idm jdm itlrec increc numrec fout.a'
        CALL EXIT(1)
      ENDIF
C
      IF     (NUMREC.LE.0) THEN
        WRITE(6,'(a)')
     &    'Error in hycom_median: numrec must be > 0'
        CALL EXIT(3)
      ENDIF
C
      NPAD = 4096 - MOD(IDM*JDM,4096)
      IF     (NPAD.EQ.4096) THEN
        NPAD = 0
      ENDIF
C
      ALLOCATE( A(IDM,JDM,NUMREC), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_median: could not allocate ',
     +             IDM*JDM*NUMREC,' words'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( AM(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_median: could not allocate ',
     +             IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
C
      CALL MEDIAN(A,AM,IDM,JDM,PAD,NPAD,
     &            ITLREC,INCREC,NUMREC, ITEST,JTEST, 
     &            CFILE1,CFILEO)
      CALL EXIT(0)
      END
      SUBROUTINE MEDIAN(A,AM,IDM,JDM,PAD,NPAD,
     &                  ITLREC,INCREC,NUMREC, ITEST,JTEST,
     &                  CFILE1,CFILEO)
      IMPLICIT NONE
C
      REAL*4     SPVAL
      PARAMETER (SPVAL=2.0**100)
C
      CHARACTER*240 CFILE1,CFILEO
      INTEGER      IDM,JDM,NPAD,ITLREC,INCREC,NUMREC,ITEST,JTEST
      REAL*4       A(IDM,JDM,NUMREC),AM(IDM,JDM),PAD(NPAD)
C
C     MOST OF WORK IS DONE HERE.
C
#ifdef sun
      INTEGER      IR_ISNAN
C
#endif
      CHARACTER*18 CASN
      INTEGER      LEN_TRIM
      INTEGER      I,J,K,KK,KS,IOS,IR,NR,NRECL
      REAL*4       AMN,AMX,SMALL,TEMP
#ifdef CRAY
      INTEGER*8    IU8,IOS8
#endif
C
      INQUIRE( IOLENGTH=NRECL) AM,PAD
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
C     READ IN ALL RECORDS.
C
      DO NR= 1,NUMREC
        IR = ITLREC + INCREC*(NR-1)
        READ(11,REC=IR,IOSTAT=IOS) A(:,:,NR)
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(A(1,1,NR),IDM*JDM)
#endif
        IF     (IOS.NE.0) THEN
          IF     (NR.EQ.1) THEN
            WRITE(6,*) 'can''t read ',TRIM(CFILE1)
            CALL EXIT(4)
          ELSE
            WRITE(6,*) TRIM(CFILE1),' is too short'
            CALL EXIT(4)
          ENDIF
        ENDIF
      ENDDO !nr
C
      CLOSE(UNIT=11)
C
      IF     (ITEST.GT.0 .AND. JTEST.GT.0) THEN
         WRITE(6,*) 'a  = ',a(itest,jtest,:)
      ENDIF
C
      AMN =  SPVAL
      AMX = -SPVAL
      DO J= 1,JDM
        DO I= 1,IDM
          IF     (A(I,J,1).EQ.SPVAL) THEN
            AM(I,J) = SPVAL
          ELSE
C
C           FIND MEDIAN VALUE USING SELECTION SORT.
C
            DO KS=1,NUMREC/2+1
              K=KS
              SMALL=A(I,J,KS)
              DO KK=KS+1,NUMREC
                IF     (A(I,J,KK).LT.SMALL) THEN
                  K=KK
                  SMALL=A(I,J,KK)
                ENDIF
              ENDDO !kk
              IF     (I.EQ.ITEST .AND. J.EQ.JTEST) THEN
                WRITE(6,*) '  ks,k,a.ks,a.k = ',
     &                        ks,k,a(i,j,ks),a(i,j,k)
              ENDIF
              TEMP=A(I,J,KS)
              A(I,J,KS)=A(I,J,K)
              A(I,J,K) =TEMP
            ENDDO !ks
C
            IF     (MOD(NUMREC,2).EQ.1) THEN !odd
              AM(I,J) =      A(I,J,NUMREC/2+1)
            ELSE !even
              AM(I,J) = 0.5*(A(I,J,NUMREC/2+1) + A(I,J,NUMREC/2))
            ENDIF !odd:even
            AMN = MIN( AMN, AM(I,J) )
            AMX = MAX( AMX, AM(I,J) )
          ENDIF !spval:else
        ENDDO !1
      ENDDO !j
C
      IF     (ITEST.GT.0 .AND. JTEST.GT.0) THEN
         WRITE(6,*) 'a  = ',a(itest,jtest,:)
         WRITE(6,*) 'am = ',am(itest,jtest)
      ENDIF
C
C     OUTPUT MEDIAN ARRAY.
C
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(AM,IDM*JDM)
#endif
      WRITE(21,REC=1,IOSTAT=IOS) AM
      WRITE(6,'(a,1p2g16.8)') 'min, max = ',AMN,AMX
C
      CLOSE(UNIT=21)
C
      RETURN
      END

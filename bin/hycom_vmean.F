      PROGRAM FVMEAN
      IMPLICIT NONE
C
C  hycom_vmean - Usage:  hycom_vmean fin.a idm jdm kdm itlrec increc numrec fout.a
C  hycom_vcnt  - Usage:  hycom_vcnt  fin.a idm jdm kdm itlrec increc numrec fout.a
C  hycom_vsdev - Usage:  hycom_vsdev fin.a idm jdm kdm itlrec increc numrec fout.a
C  hycom_vmin  - Usage:  hycom_vmin  fin.a idm jdm kdm itlrec increc numrec fout.a
C  hycom_vmax  - Usage:  hycom_vmax  fin.a idm jdm kdm itlrec increc numrec fout.a
C  hycom_vamax - Usage:  hycom_vamax fin.a idm jdm kdm itlrec increc numrec fout.a
C
C           Outputs kdm (1:idm,1:jdm) fields, representing
C           the mean/sdev/min/max/absmax of fields itl+(k-1)*inc+(n-1)*kdm*inc
C           for k=1:kdm and n=1:numrec (or n=1:e-o-f if numrec=0).
C           corresponding input fields can have different data voids,
C           hycom_vcnt outputs the number of non-voids
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
      REAL*4, ALLOCATABLE :: A(:,:),AM(:,:,:),AC(:,:,:),AQ(:,:,:)
      REAL*4              :: PAD(4096)
      INTEGER      IOS,L
      INTEGER      IARGC
      INTEGER      NARG
      CHARACTER*240 CARG
C
      INTEGER      IDM,JDM,KDM,ITLREC,INCREC,NUMREC,NPAD,ITYPE
      CHARACTER*240 CFILE1,CFILEO
C
C     READ ARGUMENTS.
C
      CALL GETARG(0,CARG)
      L = LEN_TRIM(CARG)
*     WRITE(6,"(4a)") TRIM(CARG),'"',CARG(L-4:L),'"'
      IF     (CARG(L-5:L).EQ.'_vmean') THEN
        ITYPE=1
      ELSEIF (CARG(L-4:L).EQ.'_vmax') THEN
        ITYPE=2
      ELSEIF (CARG(L-4:L).EQ.'_amax') THEN
        ITYPE=-2
      ELSEIF (CARG(L-4:L).EQ.'_vmin') THEN
        ITYPE=3
      ELSEIF (CARG(L-5:L).EQ.'_vsdev') THEN
        ITYPE=4
      ELSEIF (CARG(L-4:L).EQ.'_vcnt') THEN
        ITYPE=5
      ELSE
        WRITE(6,'(2a)')
     &    'Usage:  ',
     &  'hycom_vmean or hycom_vsdev or hycom_v[a]max or hycom_vmin ...'
        CALL EXIT(1)
      ENDIF
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
        READ(CARG,*) KDM
        CALL GETARG(5,CARG)
        READ(CARG,*) ITLREC
        CALL GETARG(6,CARG)
        READ(CARG,*) INCREC
        CALL GETARG(7,CARG)
        READ(CARG,*) NUMREC
        CALL GETARG(8,CFILEO)
      ELSE
        WRITE(6,'(2a)')
     &    'Usage:  ',
     &    'hycom_vmean fin.a idm jdm kdm itlrec increc numrec fout.a'
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
        WRITE(6,*) 'Error in hycom_vmean: could not allocate ',
     +             IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( AM(IDM,JDM,KDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_vmean: could not allocate ',
     +             IDM*JDM*KDM,' words'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( AC(IDM,JDM,KDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_vmean: could not allocate ',
     +             IDM*JDM*KDM,' words'
        CALL EXIT(2)
      ENDIF
      IF     (ITYPE.EQ.4) THEN
        ALLOCATE( AQ(IDM,JDM,KDM), STAT=IOS )
        IF     (IOS.NE.0) THEN
          WRITE(6,*) 'Error in hycom_vmean: could not allocate ',
     +               IDM*JDM*KDM,' words'
          CALL EXIT(2)
        ENDIF
      ENDIF
C
      CALL MEAN(A,AM,AC,AQ,IDM,JDM,KDM,PAD,NPAD,
     &          ITLREC,INCREC,NUMREC, ITYPE, CFILE1,CFILEO)
      CALL EXIT(0)
      END
      SUBROUTINE MEAN(A,AC,AM,AQ,IDM,JDM,KDM,PAD,NPAD,
     &                ITLREC,INCREC,NUMREC, ITYPE, CFILE1,CFILEO)
      IMPLICIT NONE
C
      REAL*4     SPVAL
      PARAMETER (SPVAL=2.0**100)
C
      CHARACTER*240 CFILE1,CFILEO
      INTEGER      IDM,JDM,KDM,NPAD,ITLREC,INCREC,NUMREC,ITYPE
      REAL*4       A( IDM,JDM),AM(IDM,JDM,KDM),
     &                         AC(IDM,JDM,KDM),
     &                         AQ(IDM,JDM,KDM),PAD(NPAD)
C
C     MOST OF WORK IS DONE HERE.
C
#ifdef sun
      INTEGER      IR_ISNAN
C
#endif
      CHARACTER*18 CASN
      INTEGER      LEN_TRIM
      INTEGER      I,J,K,IOS,IR,NR,NRECL,NUMR
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
      DO K= 1,KDM
        DO J= 1,JDM
          DO I= 1,IDM
            AC(I,J,K) = 0.0
            IF     (ITYPE.EQ.1 .OR. ITYPE.EQ.-2) THEN
              AM(I,J,K) =  0.0
            ELSEIF (ITYPE.EQ.4) THEN
              AM(I,J,K) =  0.0
              AQ(I,J,K) =  0.0
            ELSEIF (ITYPE.EQ.2) THEN
              AM(I,J,K) = -HUGE(AM(1,1,1))
            ELSEIF (ITYPE.EQ.3) THEN
              AM(I,J,K) =  HUGE(AM(1,1,1))
            ENDIF
          ENDDO
        ENDDO
      ENDDO
C
      IF     (NUMREC.EQ.0) THEN
        NUMR = 99999
      ELSE
        NUMR = NUMREC
      ENDIF
      DO 110 NR= 1,NUMR
        DO K= 1,KDM
          IR = ITLREC + INCREC*(KDM*(NR-1) + K-1)
          READ(11,REC=IR,IOSTAT=IOS) A
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
          IF     (IOS.NE.0) THEN
            IF     (NR.EQ.1 .AND. K.EQ.1) THEN
              WRITE(6,*) 'can''t read ',CFILE1(1:LEN_TRIM(CFILE1))
              CALL EXIT(4)
            ELSEIF (NUMREC.EQ.0) THEN
              IF     (K.EQ.1) THEN
                NUMREC = NR -1
                GOTO 1110
              ELSE
                WRITE(6,*) CFILE1(1:LEN_TRIM(CFILE1)),
     +                     ' not a multiple of ',KDM,' records long'
                CALL EXIT(4)
              ENDIF
            ELSE
              WRITE(6,*) CFILE1(1:LEN_TRIM(CFILE1)),' is too short'
              CALL EXIT(4)
            ENDIF
          ENDIF
          DO J= 1,JDM
            DO I= 1,IDM
#ifdef sun
              IF     (IR_ISNAN(A(I,J)).NE.1) THEN
                IF     (A(I,J).NE.SPVAL) THEN
                  AC(I,J,K) = AC(I,J,K) + 1.0
                  IF     (ITYPE.EQ.1) THEN
                    AM(I,J,K) = AM(I,J,K) + A(I,J)
                  ELSEIF (ITYPE.EQ.4) THEN
                    AM(I,J,K) = AM(I,J,K) + A(I,J)
                    AQ(I,J,K) = AQ(I,J,K) + A(I,J)**2
                  ELSEIF (ITYPE.EQ.2) THEN
                    AM(I,J,K) = MAX(AM(I,J,K), A(I,J))
                  ELSEIF (ITYPE.EQ.-2) THEN
                    IF     (ABS(A(I,J)).GT.ABS(AM(I,J,K))) THEN
                      AM(I,J,K) = A(I,J)
                    ENDIF
                  ELSEIF (ITYPE.EQ.3) THEN
                    AM(I,J,K) = MIN(AM(I,J,K), A(I,J))
                  ENDIF
                ENDIF
              ENDIF
#else
              IF     (A(I,J).NE.SPVAL) THEN
                AC(I,J,K) = AC(I,J,K) + 1.0
                IF     (ITYPE.EQ.1) THEN
                  AM(I,J,K) = AM(I,J,K) + A(I,J)
                ELSEIF (ITYPE.EQ.4) THEN
                  AM(I,J,K) = AM(I,J,K) + A(I,J)
                  AQ(I,J,K) = AQ(I,J,K) + A(I,J)**2
                ELSEIF (ITYPE.EQ.2) THEN
                  AM(I,J,K) = MAX(AM(I,J,K), A(I,J))
                ELSEIF (ITYPE.EQ.-2) THEN
                  IF     (ABS(A(I,J)).GT.ABS(AM(I,J,K))) THEN
                    AM(I,J,K) = A(I,J)
                  ENDIF
                ELSEIF (ITYPE.EQ.3) THEN
                  AM(I,J,K) = MIN(AM(I,J,K), A(I,J))
                ENDIF
              ENDIF
#endif
            ENDDO
          ENDDO
        ENDDO
  110 CONTINUE
 1110 CONTINUE
C
      DO K= 1,KDM
        AMN =  SPVAL
        AMX = -SPVAL
        DO J= 1,JDM
          DO I= 1,IDM
            IF     (AC(I,J,K).NE.0.0) THEN
              IF     (ITYPE.EQ.1) THEN
                RNUMR  = 1.0/AC(I,J,K)
                A(I,J) = AM(I,J,K)*RNUMR
              ELSEIF (ITYPE.EQ.4) THEN
                RNUMR  = 1.0/AC(I,J,K)
                A(I,J) = SQRT( MAX( 0.0,
     &                               AQ(I,J,K)*RNUMR -
     &                              (AM(I,J,K)*RNUMR)**2 ) )
              ELSEIF (ITYPE.EQ.5) THEN
                A(I,J) = AC(I,J,K)
              ELSE !2,3
                A(I,J) = AM(I,J,K)
              ENDIF
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
        WRITE(21,REC=K,IOSTAT=IOS) A
        WRITE(6,'(a,1p2g16.8)') 'min, max = ',AMN,AMX
      ENDDO
      WRITE(6,*)
      WRITE(6,*) NUMREC,' BY ',KDM,' FIELDS PROCESSED'
      WRITE(6,*)
C
      CLOSE(UNIT=11)
      CLOSE(UNIT=21)
C
      RETURN
      END

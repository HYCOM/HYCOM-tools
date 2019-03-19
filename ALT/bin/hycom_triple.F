      PROGRAM HYCOM_TRIPLE
      IMPLICIT NONE
C
C  hycom_triple - Usage:  hycom_triple fin.a idm jdm [i1 j1 idm2 jdm2] fout.a
C
C                 Outputs a 3x finer grid using bilinear interpolation
C
C                 i1,  j1   = fin(1,1) maps to fout(i1,j1); 1,1 by default
C                 idm2,jdm2 = idm*3,jdm*3 by default
C
C                 Often used to "unwind" hycom_thirdbox
C
C  fin*.a is assumed to contain idm*jdm 32-bit IEEE real values for
C   each array, in standard f77 element order, followed by padding
C   to a multiple of 4096 32-bit words, but otherwise with no control
C   bytes/words, and input values of 2.0**100 indicating a data void.
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  January 2001.
C
      REAL*4, ALLOCATABLE :: A(:,:),B(:,:),A2(:,:),B2(:,:)
      REAL*4              :: PAD(4096),PAD2(4096)
      INTEGER      IOS,IOS2,L
      INTEGER      IARGC
      INTEGER      NARG
      CHARACTER*240 CARG
C
      INTEGER      IDM,JDM,I1,J1,IDM2,JDM2,NPAD,NPAD2
      CHARACTER*240 CFILE,CFILEO
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.4) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,*) IDM
        CALL GETARG(3,CARG)
        READ(CARG,*) JDM
        I1   = 1
        J1   = 1
        IDM2 = IDM*3
        JDM2 = JDM*3
        CALL GETARG(4,CFILEO)
      ELSEIF (NARG.EQ.8) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,*) IDM
        CALL GETARG(3,CARG)
        READ(CARG,*) JDM
        CALL GETARG(4,CARG)
        READ(CARG,*) I1
        CALL GETARG(5,CARG)
        READ(CARG,*) J1
        CALL GETARG(6,CARG)
        READ(CARG,*) IDM2
        CALL GETARG(7,CARG)
        READ(CARG,*) JDM2
        CALL GETARG(8,CFILEO)
      ELSE
        WRITE(6,*)
     &'Usage: hycom_triple fin.a idm jdm [i1 j1 idm2 jdm2] fout.a'
        CALL EXIT(1)
      ENDIF
C
      IF     (I1+3*IDM-3.GT.IDM2+2 .OR.
     &        J1+3*JDM-3.GT.JDM2+2     ) THEN
        WRITE(6,*) 'Error in hycom_triple: idm2 or jdm2 too large'
        CALL EXIT(1)
      ENDIF
C
      NPAD = 4096 - MOD(IDM*JDM,4096)
      IF     (NPAD.EQ.4096) THEN
        NPAD = 0
      ENDIF
      NPAD2 = 4096 - MOD(IDM2*JDM2,4096)
      IF     (NPAD2.EQ.4096) THEN
        NPAD2 = 0
      ENDIF
C
      ALLOCATE( A(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_triple: could not allocate 1st ',
     +             IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( B(0:IDM+1,0:JDM+1), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_triple: could not allocate 2nd ',
     +             IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( A2(IDM2,JDM2), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_triple: could not allocate 3rd ',
     +             IDM2*JDM2,' words'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( B2(IDM2+I1,2), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_triple: could not allocate last ',
     +             IDM2*2,' words'
        CALL EXIT(2)
      ENDIF
C
      CALL TRIPLE(A,B,IDM,JDM,A2,B2,IDM2,JDM2,I1,J1,
     &             PAD,NPAD,PAD2,NPAD2, CFILE,CFILEO)
      CALL EXIT(0)
      END
      SUBROUTINE TRIPLE(A,B,IDM,JDM,A2,B2,IDM2,JDM2,I1,J1,
     &                  PAD,NPAD,PAD2,NPAD2, CFILE,CFILEO)
      IMPLICIT NONE
C
      REAL*4     SPVAL
      PARAMETER (SPVAL=2.0**100)
C
      CHARACTER*240 CFILE,CFILEO
      INTEGER      IDM,JDM,NPAD,IDM2,JDM2,NPAD2,I1,J1
      REAL*4       A(IDM,JDM),B(0:IDM+1,0:JDM+1),PAD(NPAD)
      REAL*4       A2(IDM2,JDM2),B2(IDM2+I1,2),PAD2(NPAD2)
C
C     MOST OF WORK IS DONE HERE.
C
#ifdef CRAY
      INTEGER*8    IU8,IOS8
#endif
      CHARACTER*18 CASN
      LOGICAL      LPERIOD,LARCTIC
      INTEGER      I,II,IQ,J,JJ,JQ,K,IOS,NRECL,NRECL2
      REAL*4       AMN,AMX,THIRD,THIRD2
      REAL*8       RC,RS
C
      THIRD  = 1.0/3.0
      THIRD2 = 1.0 - THIRD
C
      INQUIRE( IOLENGTH=NRECL)  A, PAD
      INQUIRE( IOLENGTH=NRECL2) A2,PAD2
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
      IF     (MOD(NRECL2,4096).EQ.0) THEN
        WRITE(CASN,8000) NRECL2/4096
 8000   FORMAT('-F cachea:',I4.4,':1:0')
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
      OPEN(UNIT=11, FILE=CFILE, FORM='UNFORMATTED', STATUS='OLD',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILE)
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
      OPEN(UNIT=21, FILE=CFILEO, FORM='UNFORMATTED', STATUS='NEW',
     +         ACCESS='DIRECT', RECL=NRECL2, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILEO)
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
C
      DO 110 K= 1,9999999
        READ(11,REC=K,IOSTAT=IOS) A
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
        IF     (IOS.NE.0) THEN
          IF     (K.EQ.1) THEN
            WRITE(6,*) 'can''t read ',TRIM(CFILE)
            CALL EXIT(4)
          ELSE
            GOTO 1110
          ENDIF
        ENDIF
C
        LPERIOD = MINVAL(A(IDM,:)).NE.SPVAL
        LARCTIC = MINVAL(A(:,JDM)).NE.SPVAL
*       write(6,*) 'larctic = ',larctic
*       write(6,*) 'lperiod = ',lperiod
C
        DO J= 1,JDM
          DO I= 1,IDM
            B(I,J) = A(I,J)
          ENDDO
        ENDDO
        IF     (LARCTIC) THEN
          DO I= 1,IDM
            B(I,0) = SPVAL
          ENDDO
          DO J= JDM+1,JDM+1
            JJ = JDM-1-(J-JDM)
            DO I= 1,IDM
              II = IDM-MOD(I-1,IDM)
              B(I,J) = B(II,JJ)
            ENDDO !i
          ENDDO !j
        ELSE
          DO I= 1,IDM
            B(I,    0) = SPVAL
            B(I,JDM+1) = SPVAL
          ENDDO
        ENDIF
        IF     (LPERIOD) THEN
          DO J= 0,JDM+1
            B(    0,J) = B(IDM,J)
            B(IDM+1,J) = B(  1,J)
          ENDDO
        ELSE
          DO J= 0,JDM+1
            B(    0,J) = SPVAL
            B(IDM+1,J) = SPVAL
          ENDDO
        ENDIF
C
        AMN =  SPVAL
        AMX = -SPVAL
        DO J= 1,JDM2
          DO I= 1,IDM2
            A2(I,j) = SPVAL
          ENDDO !i
        ENDDO !j
        J = 1
          DO I= 1,IDM
            II =  I1+3*I-3
*           write(6,*) 'I,II = ',I,II
            IF     (B(I,J).NE.SPVAL .AND. B(I+1,J).NE.SPVAL) THEN
              B2(II,  2) =        B(I,J)
              B2(II+1,2) = THIRD2*B(I,J) + THIRD *B(I+1,J)
              B2(II+2,2) = THIRD *B(I,J) + THIRD2*B(I+1,J)
            ELSEIF (B(I,  J).NE.SPVAL) THEN
              B2(II,  2) = B(I,J)
              B2(II+1,2) = B(I,J)
              B2(II+2,2) = B(I,J)
            ELSEIF (B(I+1,J).NE.SPVAL) THEN
              B2(II,  2) = B(I+1,J)
              B2(II+1,2) = B(I+1,J)
              B2(II+2,2) = B(I+1,J)
            ELSE
              B2(II,  2) = SPVAL
              B2(II+1,2) = SPVAL
              B2(II+2,2) = SPVAL
            ENDIF
          ENDDO !i
          IF     (LPERIOD) THEN
            DO II= 1,I1-1
*             write(6,*) 'I,II = ',0,II
              B2(II,2) = B2(II+IDM2,2)
            ENDDO !ii
          ENDIF !lperiod
          DO II= 1,IDM2
            IF     (B2(II,2).NE.SPVAL) THEN
              DO JJ= 1,J1
                A2(II,JJ) = B2(II,2)
              ENDDO !jj
              AMX = MAX( AMX, A2(II,J1) )
              AMN = MIN( AMN, A2(II,J1) )
            ELSE
              DO JJ= 1,J1
                A2(II,JJ) = SPVAL
              ENDDO !jj
            ENDIF
          ENDDO !ii
*       write(6,*) 'J,JJ = ',1,J1,AMX
        DO J= 2,JDM
          JJ = J1+3*J-3
          IF     (JJ-2.GT.JDM2) THEN
            EXIT
          ENDIF
*         write(6,*) 'J,JJ = ',J,JJ,AMX
          B2(:,1) = B2(:,2)
          DO I= 1,IDM
            II =  I1+3*I-3
            IF     (B(I,J).NE.SPVAL .AND. B(I+1,J).NE.SPVAL) THEN
              B2(II,  2) =        B(I,J)
              B2(II+1,2) = THIRD2*B(I,J) + THIRD *B(I+1,J)
              B2(II+2,2) = THIRD *B(I,J) + THIRD2*B(I+1,J)
            ELSEIF (B(I,  J).NE.SPVAL) THEN
              B2(II,  2) = B(I,J)
              B2(II+1,2) = B(I,J)
              B2(II+2,2) = B(I,J)
            ELSEIF (B(I+1,J).NE.SPVAL) THEN
              B2(II,  2) = B(I+1,J)
              B2(II+1,2) = B(I+1,J)
              B2(II+2,2) = B(I+1,J)
            ELSE
              B2(II,  2) = SPVAL
              B2(II+1,2) = SPVAL
              B2(II+2,2) = SPVAL
            ENDIF
          ENDDO !i
          IF     (LPERIOD) THEN
            DO II= 1,I1-1
              B2(II,2) = B2(II+IDM2,2)
            ENDDO !ii
          ENDIF !lperiod
          DO II= 1,IDM2
            IF     (B2(II,1).NE.SPVAL .AND. B2(II,2).NE.SPVAL) THEN
              A2(II,MIN(JDM2,JJ))   =        B2(II,2)
              A2(II,MIN(JDM2,JJ-1)) = THIRD2*B2(II,2) + THIRD *B2(II,1)
              A2(II,MIN(JDM2,JJ-2)) = THIRD *B2(II,2) + THIRD2*B2(II,1)
              AMX = MAX( AMX,
     &                   A2(II,MIN(JDM2,JJ)  ),
     &                   A2(II,MIN(JDM2,JJ-1)),
     &                   A2(II,MIN(JDM2,JJ-2)) )
              AMN = MIN( AMN,
     &                   A2(II,MIN(JDM2,JJ)  ),
     &                   A2(II,MIN(JDM2,JJ-1)),
     &                   A2(II,MIN(JDM2,JJ-2)) )
            ELSEIF (B2(II,2).NE.SPVAL) THEN
              A2(II,MIN(JDM2,JJ))   = B2(II,2)
              A2(II,MIN(JDM2,JJ-1)) = B2(II,2)
              A2(II,MIN(JDM2,JJ-2)) = B2(II,2)
              AMX = MAX( AMX, A2(II,JJ-2) )
              AMN = MIN( AMN, A2(II,JJ-2) )
            ELSEIF (B2(II,1).NE.SPVAL) THEN
              A2(II,MIN(JDM2,JJ))   = B2(II,1)
              A2(II,MIN(JDM2,JJ-1)) = B2(II,1)
              A2(II,MIN(JDM2,JJ-2)) = B2(II,1)
              AMX = MAX( AMX, A2(II,JJ-2) )
              AMN = MIN( AMN, A2(II,JJ-2) )
            ELSEIF (B2(II,2).NE.SPVAL) THEN
              A2(II,MIN(JDM2,JJ))   = SPVAL
              A2(II,MIN(JDM2,JJ-1)) = SPVAL
              A2(II,MIN(JDM2,JJ-2)) = SPVAL
            ENDIF
          ENDDO !ii
        ENDDO !j
        IF     (LARCTIC) THEN
          J  = JDM2
          JJ = JDM2-1
*         write(6,*) 'J,JJ = ',J,JJ,AMX
          DO I= 1,IDM2
            II = IDM2-MOD(I-1,IDM2)
            A2(I,J) = A2(II,JJ)
          ENDDO !i
       ENDIF !larctic
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(A2,IDM2*JDM2)
#endif
        WRITE(21,REC=K,IOSTAT=IOS) A2
        WRITE(6,'(a,1p2g16.8)')
     &     'min, max = ',AMN,AMX
  110 CONTINUE
 1110 CONTINUE
      WRITE(6,*) 
      WRITE(6,*) K-1,' FIELDS PROCESSED (IDM2,JDM2 = ',IDM2,JDM2,')'
      WRITE(6,*) 
C
      CLOSE(11)
      CLOSE(21)
C
      RETURN
      END

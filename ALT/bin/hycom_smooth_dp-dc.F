      PROGRAM FEXPR
      IMPLICIT NONE
C
C  hycom_smooth_dp-dc - Usage:  hycom_smooth_dp-dc fin.a idm jdm [s] fout.a
C
C                 Outputs a smoothed version of each input field.
C
C  The smoother is "delsq_plus - delsq_cross", from:
C  Killworth, P.D., D. Stainforth. D.J. Webb, S.M. Paterson, 1991,
C  The development of a free-surface Bryan-Cox-Semtner ocean model,
C  JPO, 21, 1333-1348.
C
C  fin*.a is assumed to contain idm*jdm 32-bit IEEE real values for
C   each array, in standard f77 element order, followed by padding
C   to a multiple of 4096 32-bit words, but otherwise with no control
C   bytes/words, and input values of 2.0**100 indicating a data void.
C
C  s is a scale factor betweem 0.0 and 1.0.
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  January 2001.
C
      REAL*4, ALLOCATABLE :: A(:,:),B(:,:),C(:,:)
      REAL*4              :: PAD(4096)
      INTEGER IOS
      INTEGER      IARGC
      INTEGER      NARG
      CHARACTER*240 CARG
C
      INTEGER       IDM,JDM,NPAD
      REAL*4        SS
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
        SS = 1.0
        CALL GETARG(4,CFILEO)
      ELSEIF (NARG.EQ.5) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,*) IDM
        CALL GETARG(3,CARG)
        READ(CARG,*) JDM
        CALL GETARG(4,CARG)
        READ(CARG,*) SS
        CALL GETARG(5,CFILEO)
      ELSE
        WRITE(6,*)
     &    'Usage: hycom_smooth_dp-dc fin.a fin2.a idm jdm [s] fout.a'
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
        WRITE(6,*)
     +    'Error in hycom_smooth_dp-dc: could not allocate 1st ',
     +    IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( B(0:IDM+1,0:JDM+1), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 
     +    'Error in hycom_smooth_dp-dc: could not allocate 2nd ',
     +    IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( C(-1:IDM+2,-1:JDM+2), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 
     +    'Error in hycom_smooth_dp-dc: could not allocate 3rd ',
     +    IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
C
      CALL SMTH(A,B,C,IDM,JDM,PAD,NPAD, SS, CFILE,CFILEO)
      CALL EXIT(0)
      END
      SUBROUTINE SMTH(A,B,C,IDM,JDM,PAD,NPAD, SS, CFILE,CFILEO)
      IMPLICIT NONE
C
      REAL*4     SPVAL
      PARAMETER (SPVAL=2.0**100)
C
      CHARACTER*240 CFILE,CFILEO
      INTEGER      IDM,JDM,NPAD
      REAL*4       A(IDM,JDM),B( 0:IDM+1, 0:JDM+1),
     &                        C(-1:IDM+2,-1:JDM+2),PAD(NPAD),SS
C
C     MOST OF WORK IS DONE HERE.
C
#ifdef CRAY
      INTEGER*8    IU8,IOS8
#endif
      CHARACTER*18 CASN
      INTEGER      LEN_TRIM
      INTEGER      I,IQ,ISM,J,JQ,K,IOS,NRECL
      REAL*4       AMN,AMX,QC,RS,DC1,DC2,DC3,DC4,DP1,DP2,DP3,DP4
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
      OPEN(UNIT=11, FILE=CFILE, FORM='UNFORMATTED', STATUS='OLD',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',CFILE(1:LEN_TRIM(CFILE))
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
      DO 110 K= 1,9999
        READ(11,REC=K,IOSTAT=IOS) A
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
        IF     (IOS.NE.0) THEN
          IF     (K.EQ.1) THEN
            WRITE(6,*) 'can''t read ',CFILE(1:LEN_TRIM(CFILE))
            CALL EXIT(4)
          ELSE
            GOTO 1110
          ENDIF
        ENDIF
C
C       SMOOTHING PASS.
C
        DO J= 1,JDM
          DO I= 1,IDM
            C(I,J) = A(I,J)
          ENDDO
        ENDDO
        DO J= 1,JDM
          C(   -1,J) = SPVAL
          C(    0,J) = SPVAL
          C(IDM+1,J) = SPVAL
          C(IDM+2,J) = SPVAL
        ENDDO
        DO I=-1,IDM+2
          C(I,   -1) = SPVAL
          C(I,    0) = SPVAL
          C(I,JDM+1) = SPVAL
          C(I,JDM+2) = SPVAL
        ENDDO
C
        DO J= 0,JDM+1
          DO I= 0,IDM+1
            IF     (C(I,J).NE.SPVAL) THEN
              B(I,J) = C(I,J)
            ELSE
              RS = 0.0
              QC = 0.0
              IF     (C(I+1,J).NE.SPVAL) THEN
                RS = RS + C(I+1,J)
                QC = QC + 1.0
              ENDIF
              IF     (C(I-1,J).NE.SPVAL) THEN
                RS = RS + C(I-1,J)
                QC = QC + 1.0
              ENDIF
              IF     (C(I,J+1).NE.SPVAL) THEN
                RS = RS + C(I,J+1)
                QC = QC + 1.0
              ENDIF
              IF     (C(I,J-1).NE.SPVAL) THEN
                RS = RS + C(I,J-1)
                QC = QC + 1.0
              ENDIF
              if     (rs.gt.0.5*spval) then
                write(6,*) 'error - i,j,b = ',
     &            i,j,b(i,j),b(i+1,j),b(i-1,j),b(i,j+1),b(i,j-1)
              endif
              IF     (QC.NE.0.0) THEN
                B(I,J) = RS/QC
              ELSE
                B(I,J) = SPVAL
              ENDIF
            ENDIF
          ENDDO
        ENDDO
C
        AMN =  SPVAL
        AMX = -SPVAL
        DO J= 1,JDM
          DO I= 1,IDM
            IF     ( C(I,J).NE.SPVAL) THEN
              DP1 = (B(I+1,J) - B(I,J))
              DP2 = (B(I-1,J) - B(I,J))
              DP3 = (B(I,J+1) - B(I,J))
              DP4 = (B(I,J-1) - B(I,J))
              if     (b(i+1,j).eq.spval) then
                write(6,*) 'error - i,j,dp1,b = ',
     &            i,j,dp1,b(i+1,j),b(i,j)
              endif
              if     (b(i-1,j).eq.spval) then
                write(6,*) 'error - i,j,dp2,b = ',
     &            i,j,dp2,b(i-1,j),b(i,j)
              endif
              if     (b(i,j+1).eq.spval) then
                write(6,*) 'error - i,j,dp3,b = ',
     &            i,j,dp3,b(i,j+1),b(i,j)
              endif
              if     (b(i,j-1).eq.spval) then
                write(6,*) 'error - i,j,dp4,b = ',
     &            i,j,dp4,b(i,j-1),b(i,j)
              endif
              IF     ( C(I+1,J+1).NE.SPVAL) THEN
                DC1 = (B(I+1,J+1) - B(I,J))
              ELSE
                DC1 = 0.0
              ENDIF
              IF     ( C(I-1,J+1).NE.SPVAL) THEN
                DC2 = (B(I-1,J+1) - B(I,J))
              ELSE
                DC2 = 0.0
              ENDIF
              IF     ( C(I+1,J-1).NE.SPVAL) THEN
                DC3 = (B(I+1,J-1) - B(I,J))
              ELSE
                DC3 = 0.0
              ENDIF
              IF     ( C(I-1,J-1).NE.SPVAL) THEN
                DC4 = (B(I-1,J-1) - B(I,J))
              ELSE
                DC4 = 0.0
              ENDIF
              A(I,J) = A(I,J) + SS*(      (DP1+DP2+DP3+DP4)
     &                               -0.5*(DC1+DC2+DC3+DC4) )
              AMX = MAX( AMX, A(I,J) )
              AMN = MIN( AMN, A(I,J) )
            ENDIF
          ENDDO
        ENDDO
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
        WRITE(21,REC=K,IOSTAT=IOS) A
        WRITE(6,'(a,1p2g16.8)')
     &     'min, max = ',AMN,AMX
  110 CONTINUE
 1110 CONTINUE
      WRITE(6,*) 
      WRITE(6,*) K-1,' FIELDS PROCESSED'
      WRITE(6,*) 
C
      CLOSE(11)
      CLOSE(21)
C
      RETURN
      END

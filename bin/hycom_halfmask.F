      PROGRAM HYCOM_HALFMASK
      IMPLICIT NONE
C
C  hycom_halfmask     - Usage:  hycom_halfmask     fin.a idm jdm bl bh fout.a [inc]
C  hycom_halfmask_inv - Usage:  hycom_halfmask_inv fin.a idm jdm bl bh fout.a [inc]
C
C                 Partially masks fields in fin where they are between bl and bh
C                 For hycom_halfmask_inv, where they are not between bl and bh
C
C  fin.a and are assumed to contain idm*jdm 32-bit IEEE real values
C   for each array, in standard f77 element order, followed by padding
C   to a multiple of 4096 32-bit words, but otherwise with no control
C   bytes/words, and input values of 2.0**100 indicating a data void.
C
C  the fields in fout.a will be copies of those in fin.a except that
C   they will additionally have data voids in 3 out of 4 locations where
C   the original values were between bl and bh inclusive (hycom_halfmask)
C   or where the original values were not between bl and bh inclusive 
C   (hycom_halfmask_inv).  The values at points (1:idm:2,1:jdm:2) are
C   unchanged.
C
C  if inc is present inc^2-1 out of inc^2 locations will have data voids
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  August 2002.
C
      REAL*4, ALLOCATABLE :: A(:,:)
      REAL*4              :: PAD(4096)
      INTEGER       IOS,L
      INTEGER       IARGC
      INTEGER       NARG
      CHARACTER*240 CARG
C
      INTEGER       IDM,JDM,INC,NPAD,ITYPE
      REAL*4        BL,BH,BN
      CHARACTER*240 CFILE1,CFILEO
C
C     READ ARGUMENTS.
C
      CALL GETARG(0,CARG)
      L = LEN_TRIM(CARG)
*     WRITE(6,"(4a)") TRIM(CARG),'"',CARG(L-4:L),'"'
      IF     (CARG(L-13:L).EQ.'hycom_halfmask') THEN
        ITYPE=1
      ELSEIF (CARG(L-17:L).EQ.'hycom_halfmask_inv') THEN
        ITYPE=2
      ELSE
        WRITE(6,'(2a)')
     &    'Usage:  ',
     &    'hycom_halfmask or hycom_halfmask_inv ...'
        CALL EXIT(1)
      ENDIF

      NARG = IARGC()
C
      IF     (NARG.EQ.6) THEN
        CALL GETARG(1,CFILE1)
        CALL GETARG(2,CARG)
        READ(CARG,*) IDM
        CALL GETARG(3,CARG)
        READ(CARG,*) JDM
        CALL GETARG(4,CARG)
        READ(CARG,*) BL
        CALL GETARG(5,CARG)
        READ(CARG,*) BH
        CALL GETARG(6,CFILEO)
        INC = 2
        BN  = 2.0**100  !data void marker
      ELSEIF (NARG.EQ.7) THEN
        CALL GETARG(1,CFILE1)
        CALL GETARG(2,CARG)
        READ(CARG,*) IDM
        CALL GETARG(3,CARG)
        READ(CARG,*) JDM
        CALL GETARG(4,CARG)
        READ(CARG,*) BL
        CALL GETARG(5,CARG)
        READ(CARG,*) BH
        CALL GETARG(6,CFILEO)
        CALL GETARG(7,CARG)
        READ(CARG,*) INC
        BN = 2.0**100  !data void marker
      ELSEIF (NARG.EQ.8) THEN  !undocumented option
        CALL GETARG(1,CFILE1)
        CALL GETARG(2,CARG)
        READ(CARG,*) IDM
        CALL GETARG(3,CARG)
        READ(CARG,*) JDM
        CALL GETARG(4,CARG)
        READ(CARG,*) BL
        CALL GETARG(5,CARG)
        READ(CARG,*) BH
        CALL GETARG(6,CFILEO)
        CALL GETARG(7,CARG)
        READ(CARG,*) INC
        CALL GETARG(8,CARG)
        READ(CARG,*) BN
      ELSE
        WRITE(6,*)
     &    'Usage: hycom_halfmask fin.a idm jdm bl bh fout.a [inc]'
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
        WRITE(6,*) 'Error in hycom_halfmask: could not allocate 1st ',
     +             IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
C
      CALL HALF(A,IDM,JDM,PAD,NPAD, INC,BL,BH,BN, ITYPE, CFILE1,CFILEO)
      CALL EXIT(0)
      END
      SUBROUTINE HALF(A,IDM,JDM,PAD,NPAD,
     +                INC,BL,BH,BN, ITYPE, CFILE1,CFILEO)
      IMPLICIT NONE
C
      REAL*4     SPVAL
      PARAMETER (SPVAL=2.0**100)
C
      CHARACTER*240 CFILE1,CFILEO
      INTEGER      IDM,JDM,NPAD,INC,ITYPE
      REAL*4       A(IDM,JDM),PAD(NPAD), BL,BH,BN
C
C     MOST OF WORK IS DONE HERE.
C
#ifdef sun
      INTEGER      IR_ISNAN
C
#endif
      CHARACTER*18 CASN
      INTEGER      I,J,K,IOS,N,NRECL
      REAL*4       AMN,AMX
#ifdef CRAY
      INTEGER*8    IU8,IOS8
#endif
C
      IF     (NPAD.EQ.0) THEN
        INQUIRE( IOLENGTH=NRECL) A
      ELSE
        INQUIRE( IOLENGTH=NRECL) A,PAD
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
      DO 110 K= 1,9999
        READ(11,REC=K,IOSTAT=IOS) A
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
        IF     (IOS.NE.0) THEN
          IF     (K.EQ.1) THEN
            WRITE(6,*) 'can''t read ',TRIM(CFILE1)
            CALL EXIT(4)
          ELSE
            GOTO 1110
          ENDIF
        ENDIF
        AMN =  SPVAL
        AMX = -SPVAL
        N   =  0
        DO J= 1,JDM
          DO I= 1,IDM
            IF     (MOD(J,INC).EQ.1 .AND. MOD(I,INC).EQ.1) THEN
C ---         leave alone every (INC**2)-th point, usually INC=2
#ifdef sun
              IF     (IR_ISNAN(A(I,J)).NE.1) THEN
                IF     (A(I,J).NE.SPVAL) THEN
                  AMN = MIN( AMN, A(I,J) )
                  AMX = MAX( AMX, A(I,J) )
                ENDIF
              ENDIF
#else
              IF     (A(I,J).NE.SPVAL) THEN
                AMN = MIN( AMN, A(I,J) )
                AMX = MAX( AMX, A(I,J) )
              ENDIF
#endif
              CYCLE
            ENDIF
#ifdef sun
            IF     (IR_ISNAN(A(I,J)).NE.1) THEN
              IF     (A(I,J).NE.SPVAL) THEN
                IF     (ITYPE.EQ.1) THEN !hycom_halfmask
                  IF     (A(I,J).GE.BL .AND.
     +                    A(I,J).LE.BH      ) THEN
                    A(I,J) = BN
                    N = N+1
                  ELSE
                    AMN = MIN( AMN, A(I,J) )
                    AMX = MAX( AMX, A(I,J) )
                  ENDIF
                ELSE !hycom_halfmask_inv
                  IF     (A(I,J).GE.BL .AND.
     +                    A(I,J).LE.BH      ) THEN
                    AMN = MIN( AMN, A(I,J) )
                    AMX = MAX( AMX, A(I,J) )
                  ELSE
                    A(I,J) = BN
                    N = N+1
                  ENDIF
                ENDIF !itype
              ENDIF
            ENDIF
#else
            IF     (A(I,J).NE.SPVAL) THEN
              IF     (ITYPE.EQ.1) THEN !hycom_halfmask
                IF     (A(I,J).GE.BL .AND.
     +                  A(I,J).LE.BH      ) THEN
                  A(I,J) = BN
                  N = N+1
                ELSE
                  AMN = MIN( AMN, A(I,J) )
                  AMX = MAX( AMX, A(I,J) )
                ENDIF
              ELSE !hycom_halfmask_inv
                IF     (A(I,J).GE.BL .AND.
     +                  A(I,J).LE.BH      ) THEN
                  AMN = MIN( AMN, A(I,J) )
                  AMX = MAX( AMX, A(I,J) )
                ELSE
                  A(I,J) = BN
                  N = N+1
                ENDIF
              ENDIF !itype
            ENDIF
#endif
          ENDDO !i
        ENDDO !j
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
        IF     (NPAD.EQ.0) THEN
          WRITE(21,REC=K,IOSTAT=IOS) A
        ELSE
          WRITE(21,REC=K,IOSTAT=IOS) A,PAD
        ENDIF
        WRITE(6,'(a,1p2g16.8,i12)')
     &     'min, max, masked = ',AMN,AMX,N
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

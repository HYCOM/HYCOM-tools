      PROGRAM FMASS
      IMPLICIT NONE
C
C  hycom_xward     - Usage:  hycom_xward uin.a vin.a xout.a yout.a [grid.a]
C                            hycom_xward uvin.a xyout.a [grid.a]
C                            hycom_eward xin.a yin.a eout.a nout.a [grid.a]
C                            hycom_eward xyin.a uvout.a [grid.a]
C
C                 rotates p-grid vectors between [en]-ward and [xy]-ward
C
C                 grid.a is a hycom grid file, default regional.grid.a.
C                 Note that the corresponding grid.b must also exist.
C
C                 idm,jdm are taken from grid.a and the input array
C                 is assumed to be on the p-grid (i.e. plat from grid.a
C                 represents the field's latitude).
C
C  *in.a is assumed to contain idm*jdm 32-bit IEEE real values for
C   each array, in standard f77 element order, followed by padding
C   to a multiple of 4096 32-bit words, but otherwise with no control
C   bytes/words, and input values of 2.0**100 indicating a data void.
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  November 2008.
C
      REAL*4, ALLOCATABLE :: U(:,:),V(:,:),PANG(:,:)
      REAL*4              :: PAD(4096)
      INTEGER       IOS,L
      INTEGER       IARGC
      INTEGER       NARG
      CHARACTER*240 CARG
C
      LOGICAL       XWARD,YFILE
      INTEGER       IDM,JDM,NPAD
      REAL*4        CELL
      CHARACTER*6   CVARIN
      CHARACTER*240 CFILE1,CFILE2,CFILE3,CFILE4,CFILEG,CFILEB
C
C     READ ARGUMENTS.
C
      CALL GETARG(0,CARG)
      XWARD = CARG.EQ.'hycom_xward'
C
      NARG = IARGC()
C
      IF     (NARG.EQ.2) THEN
        CALL GETARG(1,CFILE1)
        CALL GETARG(2,CFILE3)
        CFILEG = 'regional.grid.a'
        YFILE  = .false.  !no separate y-ward file
      ELSEIF (NARG.EQ.3) THEN
        CALL GETARG(1,CFILE1)
        CALL GETARG(2,CFILE3)
        CALL GETARG(3,CFILEG)
        YFILE  = .false.  !no separate y-ward file
      ELSEIF (NARG.EQ.4) THEN
        CALL GETARG(1,CFILE1)
        CALL GETARG(2,CFILE2)
        CALL GETARG(3,CFILE3)
        CALL GETARG(4,CFILE4)
        CFILEG = 'regional.grid.a'
        YFILE  = .true.   !   separate y-ward file
      ELSEIF (NARG.EQ.5) THEN
        CALL GETARG(1,CFILE1)
        CALL GETARG(2,CFILE2)
        CALL GETARG(3,CFILE3)
        CALL GETARG(4,CFILE4)
        CALL GETARG(5,CFILEG)
        YFILE  = .true.   !   separate y-ward file
      ELSEIF (XWARD) THEN
        WRITE(6,*) 
     +   'Usage:  hycom_xward uin.a vin.a xout.a yout.a [grid.a]'
        CALL EXIT(1)
      ELSE
        WRITE(6,*) 
     +   'Usage:  hycom_eward xin.a yin.a eout.a nout.a [grid.a]'
        CALL EXIT(1)
      ENDIF
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
        WRITE(6,*) 'hycom_xward: bad header file ',
     &             TRIM(CFILEB)
        CALL EXIT(2)
      ENDIF
      READ( 11,*) JDM,CVARIN
      IF (CVARIN.NE.'jdm   ') THEN
        WRITE(6,*) 'hycom_xward: bad header file ',
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
      ALLOCATE( U(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_xward: could not allocate ',
     +             IDM*JDM,' words for U'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( V(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_xward: could not allocate ',
     +             IDM*JDM,' words for V'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( PANG(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_xward: could not allocate ',
     +             IDM*JDM,' words for PANG'
        CALL EXIT(2)
      ENDIF
C
      CALL WARD(U,V,PANG,IDM,JDM,PAD,NPAD, 
     +          XWARD,YFILE, CFILE1,CFILE2,CFILE3,CFILE4,CFILEG)
      CALL EXIT(0)
 5000 FORMAT(I4)
      END
      SUBROUTINE WARD(U,V,PANG,IDM,JDM,PAD,NPAD,
     +                XWARD,YFILE, CFILE1,CFILE2,CFILE3,CFILE4,CFILEG)
      IMPLICIT NONE
C
      REAL*4     SPVAL
      PARAMETER (SPVAL=2.0**100)
C
      CHARACTER*240 CFILE1,CFILE2,CFILE3,CFILE4,CFILEG
      LOGICAL       XWARD,YFILE
      INTEGER       IDM,JDM,NPAD
      REAL*4        U(IDM,JDM),V(IDM,JDM),PANG(IDM,JDM),PAD(NPAD)
C
C     MOST OF WORK IS DONE HERE.
C
#ifdef sun
      INTEGER      IR_ISNAN
C
#endif
      CHARACTER*18 CASN
      INTEGER      I,J,K,IOS,NRECL
      REAL*4       UMN,UMX,VMN,VMX,X,Y
#ifdef CRAY
      INTEGER*8    IU8,IOS8
#endif
C
C     INPUT GRID ARRAY.
C
      OPEN(UNIT=11, FILE=CFILEG, FORM='UNFORMATTED', STATUS='OLD',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILEG)
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
C
      READ(11,REC=9,IOSTAT=IOS) PANG
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(PANG, IDM*JDM)
#endif
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'can''t read ',TRIM(CFILEG)
        CALL EXIT(4)
      ENDIF
C
      CLOSE(UNIT=11)
C
C     OPEN INPUT FILES.
C
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
        IU8 = 12
        CALL ASNUNIT(IU8,CASN,IOS8)
        IF     (IOS8.NE.0) THEN
          write(6,*) 'Error: can''t asnunit 12'
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
        IU8 = 22
        CALL ASNUNIT(IU8,CASN,IOS8)
        IF     (IOS8.NE.0) THEN
          write(6,*) 'Error: can''t asnunit 22'
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
      CALL ASNUNIT(12,'-F syscall -N ieee',IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t asnunit 12'
        write(6,*) 'ios = ',ios
        CALL EXIT(5)
      ENDIF
      CALL ASNUNIT(21,'-F syscall -N ieee',IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t asnunit 21'
        write(6,*) 'ios = ',ios
        CALL EXIT(5)
      ENDIF
      CALL ASNUNIT(22,'-F syscall -N ieee',IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t asnunit 22'
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
      OPEN(UNIT=21, FILE=CFILE3, FORM='UNFORMATTED', STATUS='NEW',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILE3)
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
      IF     (YFILE) THEN
        OPEN(UNIT=12, FILE=CFILE2, FORM='UNFORMATTED', STATUS='OLD',
     +           ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
        IF     (IOS.NE.0) THEN
          write(6,*) 'Error: can''t open ',TRIM(CFILE2)
          write(6,*) 'ios   = ',ios
          write(6,*) 'nrecl = ',nrecl
          CALL EXIT(3)
        ENDIF
        OPEN(UNIT=22, FILE=CFILE4, FORM='UNFORMATTED', STATUS='NEW',
     +           ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
        IF     (IOS.NE.0) THEN
          write(6,*) 'Error: can''t open ',TRIM(CFILE4)
          write(6,*) 'ios   = ',ios
          write(6,*) 'nrecl = ',nrecl
          CALL EXIT(3)
        ENDIF
      ENDIF
C
C     LOOP OVER ALL FIELDS.
C
      DO K= 1,9999
        IF     (YFILE) THEN
          READ(11,REC=K,IOSTAT=IOS) U
#ifdef ENDIAN_IO
          CALL ENDIAN_SWAP(U,IDM*JDM)
#endif
          IF     (IOS.NE.0) THEN
            IF     (K.EQ.1) THEN
              WRITE(6,*) 'can''t read ',TRIM(CFILE1)
              CALL EXIT(4)
            ELSE
              EXIT !k
            ENDIF
          ENDIF !ios
          READ(12,REC=K,IOSTAT=IOS) V
#ifdef ENDIAN_IO
          CALL ENDIAN_SWAP(V,IDM*JDM)
#endif
          IF     (IOS.NE.0) THEN
            IF     (K.EQ.1) THEN
              WRITE(6,*) 'can''t read ',TRIM(CFILE2)
              CALL EXIT(4)
            ELSE
              EXIT !k
            ENDIF
          ENDIF !ios
        ELSE
          READ(11,REC=2*K-1,IOSTAT=IOS) U
#ifdef ENDIAN_IO
          CALL ENDIAN_SWAP(U,IDM*JDM)
#endif
          IF     (IOS.NE.0) THEN
            IF     (K.EQ.1) THEN
              WRITE(6,*) 'can''t read ',TRIM(CFILE1)
              CALL EXIT(4)
            ELSE
              WRITE(6,*) 'odd number of records in ',TRIM(CFILE1)
              CALL EXIT(4)
            ENDIF
          ENDIF !ios
          READ(11,REC=2*K,IOSTAT=IOS) V
#ifdef ENDIAN_IO
          CALL ENDIAN_SWAP(V,IDM*JDM)
#endif
          IF     (IOS.NE.0) THEN
            EXIT !k
          ENDIF !ios
        ENDIF
C
        UMN  =  SPVAL
        UMX  = -SPVAL
        VMN  =  SPVAL
        VMX  = -SPVAL
        DO J= 1,JDM
          DO I= 1,IDM
            IF     (U(I,J).NE.SPVAL .AND.
     +              V(I,J).NE.SPVAL      ) THEN
              X = U(I,J)
              Y = V(I,J)
              IF     (XWARD) THEN 
C               EAST-NORTH TO X-Y
                U(I,J) = COS(-PANG(I,J))*X +
     &                   SIN( PANG(I,J))*Y
                V(I,J) = COS(-PANG(I,J))*Y -
     &                   SIN( PANG(I,J))*X
              ELSE
C               X-Y TO EAST-NORTH
                U(I,J) = COS( PANG(I,J))*X +
     &                   SIN(-PANG(I,J))*Y
                V(I,J) = COS( PANG(I,J))*Y -
     &                   SIN(-PANG(I,J))*X
              ENDIF
              UMN  = MIN( UMN, U(I,J) )
              UMX  = MAX( UMX, U(I,J) )
              VMN  = MIN( VMN, V(I,J) )
              VMX  = MAX( VMX, V(I,J) )
            ENDIF
          ENDDO
        ENDDO
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(U,IDM*JDM)
        CALL ENDIAN_SWAP(V,IDM*JDM)
#endif
        IF     (YFILE) THEN
          WRITE(21,REC=  K,  IOSTAT=IOS) U
          WRITE(22,REC=  K,  IOSTAT=IOS) V
        ELSE
          WRITE(21,REC=2*K-1,IOSTAT=IOS) U
          WRITE(21,REC=2*K,  IOSTAT=IOS) V
        ENDIF
        IF     (XWARD) THEN
          WRITE(6,'(a,1p2g16.8)')
     &       'x-ward: min,max = ',UMN,UMX
          WRITE(6,'(a,1p2g16.8)')
     &       'y-ward: min,max = ',VMN,VMX
        ELSE
          WRITE(6,'(a,1p2g16.8)')
     &       'e-ward: min,max = ',UMN,UMX
          WRITE(6,'(a,1p2g16.8)')
     &       'n-ward: min,max = ',VMN,VMX
        ENDIF
      ENDDO !k
C
      CLOSE(UNIT=21)
      CLOSE(UNIT=11)
      IF     (YFILE) THEN
        CLOSE(UNIT=22)
        CLOSE(UNIT=12)
      ENDIF
C
      RETURN
      END

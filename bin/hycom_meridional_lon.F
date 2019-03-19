      PROGRAM FMERID
      IMPLICIT NONE
C
C  hycom_meridional_lon - Usage:  hycom_meridional_lon file.a k cell [grid.a] [mask.a]
C
C                 prints the mean and rms of meridional longitudinal extents
C                 from the k-th (1:idm,1:jdm) array in file.a
C
C                 longitude bands are centered on 0E and are cell degrees
C                 wide, e.g. (N-0.5)*cell to (N+0.5)*cell.
C
C                 grid.a is a hycom grid file, default regional.grid.a.
C                 Note that the corresponding grid.b must also exist.
C
C                 mask.a contains an optional mask array, that turns
C                 off sampling anywhere it is 2.0**100 (data void).
C                 Note that grid.a must be specified if mask.a is used.
C
C                 idm,jdm are taken from grid.a and the input array
C                 is assumed to be on the p-grid (i.e. plon from grid.a
C                 represents the field's longitude).
C
C  file.a is assumed to contain idm*jdm 32-bit IEEE real values for
C   each array, in standard f77 element order, followed by padding
C   to a multiple of 4096 32-bit words, but otherwise with no control
C   bytes/words, and input values of 2.0**100 indicating a data void.
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  January 2001.
C
      REAL*4, ALLOCATABLE :: A(:,:),AMSK(:,:),PLON(:,:),PSC2(:,:)
      REAL*4              :: PAD(4096)
      INTEGER IOS
      INTEGER      IARGC
      INTEGER      NARG
      CHARACTER*240 CARG
C
      INTEGER       IDM,JDM,K,NPAD
      REAL*4        CELL
      CHARACTER*6   CVARIN
      CHARACTER*240 CFILE,CFILEG,CFILEM,CFILEB
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.3) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,5000) K
        CALL GETARG(3,CARG)
        READ(CARG,*)    CELL
        CFILEG = 'regional.grid.a'
        CFILEM = 'NONE'
      ELSEIF (NARG.EQ.4) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,5000) K
        CALL GETARG(3,CARG)
        READ(CARG,*)    CELL
        CALL GETARG(4,CFILEG)
        CFILEM = 'NONE'
      ELSEIF (NARG.EQ.5) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,5000) K
        CALL GETARG(3,CARG)
        READ(CARG,*)    CELL
        CALL GETARG(4,CFILEG)
        CALL GETARG(5,CFILEM)
      ELSE
        WRITE(6,*) 
     +   'Usage:  hycom_meridional_lon file.a k cell [grid.a] [mask.a]'
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
        WRITE(6,*) 'hycom_meridional_lon: bad header file ',
     &             CFILEB(1:LEN_TRIM(CFILEB))
        CALL EXIT(2)
      ENDIF
      READ( 11,*) JDM,CVARIN
      IF (CVARIN.NE.'jdm   ') THEN
        WRITE(6,*) 'hycom_meridional_lon: bad header file ',
     &             CFILEB(1:LEN_TRIM(CFILEB))
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
      ALLOCATE( A(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_meridional: could not allocate ',
     +             IDM*JDM,' words for A'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( AMSK(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_meridional: could not allocate ',
     +             IDM*JDM,' words for AMSK'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( PLON(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_meridional: could not allocate ',
     +             IDM*JDM,' words for PLON'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( PSC2(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_meridional: could not allocate ',
     +             IDM*JDM,' words for PSC2'
        CALL EXIT(2)
      ENDIF
C
      CALL MERID(A,AMSK,PLON,PSC2,IDM,JDM,PAD,NPAD, K,CELL,
     +           CFILE,CFILEG,CFILEM)
      CALL EXIT(0)
 5000 FORMAT(I4)
      END
      SUBROUTINE MERID(A,AMSK,PLON,PSC2,IDM,JDM, PAD,NPAD, K,CELL,
     +                 CFILE,CFILEG,CFILEM)
      IMPLICIT NONE
C
      REAL*4     SPVAL
      PARAMETER (SPVAL=2.0**100)
C
      CHARACTER*240 CFILE,CFILEG,CFILEM
      INTEGER      IDM,JDM,NPAD,K
      REAL*4       A(IDM,JDM),AMSK(IDM,JDM),
     +             PLON(IDM,JDM),PSC2(IDM,JDM),PAD(NPAD),CELL
C
C     MOST OF WORK IS DONE HERE.
C
      REAL*8, ALLOCATABLE :: SUMA(:),SUMB(:),SUMR(:)
C
#ifdef sun
      INTEGER      IR_ISNAN
C
#endif
      CHARACTER*18 CASN
      INTEGER      LEN_TRIM
      INTEGER      I,IF,IL,IZ,IOS,J,NRECL
      REAL*8       SUMT
#ifdef CRAY
      INTEGER*8    IU8,IOS8
#endif
C
C     INPUT ARRAY A.
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
      ENDIF
#else
      CALL ASNUNIT(11,'-F syscall -N ieee',IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t asnunit 11'
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
C
      READ(11,REC=K,IOSTAT=IOS) A
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'can''t read ',CFILE(1:LEN_TRIM(CFILE))
        CALL EXIT(4)
      ENDIF
C
      CLOSE(UNIT=11)
C
C     INPUT MASK ARRAY.
C
      IF     (CFILEM.EQ.'NONE') THEN
        AMSK(:,:) = 0.0
      ELSE
        OPEN(UNIT=11, FILE=CFILEM, FORM='UNFORMATTED', STATUS='OLD',
     +           ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
        IF     (IOS.NE.0) THEN
          write(6,*) 'Error: can''t open ',CFILEM(1:LEN_TRIM(CFILEM))
          write(6,*) 'ios   = ',ios
          write(6,*) 'nrecl = ',nrecl
          CALL EXIT(3)
        ENDIF
C
        READ(11,REC=1,IOSTAT=IOS) AMSK
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(AMSK,IDM*JDM)
#endif
        IF     (IOS.NE.0) THEN
          WRITE(6,*) 'can''t read ',CFILEM(1:LEN_TRIM(CFILEM))
          CALL EXIT(4)
        ENDIF
      ENDIF
C
      CLOSE(UNIT=11)
C
C     INPUT GRID ARRAYS.
C
      OPEN(UNIT=11, FILE=CFILEG, FORM='UNFORMATTED', STATUS='OLD',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',CFILEG(1:LEN_TRIM(CFILEG))
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
C
      READ(11,REC=10,IOSTAT=IOS) PLON  ! pscx
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(PLON,IDM*JDM)
#endif
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'can''t read ',CFILEG(1:LEN_TRIM(CFILEG))
        CALL EXIT(4)
      ENDIF
C
      READ(11,REC=11,IOSTAT=IOS) PSC2  ! pscy
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(PSC2,IDM*JDM)
#endif
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'can''t read ',CFILEG(1:LEN_TRIM(CFILEG))
        CALL EXIT(4)
      ENDIF
C
      PSC2(:,:) = PSC2(:,:)*PLON(:,:)  ! psc2=pscx*pscy
C
      READ(11,REC= 1,IOSTAT=IOS) PLON  ! plon
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(PLON,IDM*JDM)
#endif
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'can''t read ',CFILEG(1:LEN_TRIM(CFILEG))
        CALL EXIT(4)
      ENDIF
C
      CLOSE(UNIT=11)
C
C     CONVERT LONGITUDE TO CELL SPACE.
C
      DO J= 1,JDM
        DO I= 1,IDM
          PLON(I,J) = MOD( PLON(I,J) + 1080.0, 360.0 )
          IF     (PLON(I,J).GT.(NINT(180.0/CELL)+0.5)*CELL) THEN
            PLON(I,J) = PLON(I,J) - 360.0
          ENDIF
          PLON(I,J) = NINT(PLON(I,J)/CELL)
        ENDDO
      ENDDO
      IF = MINVAL(PLON(:,:))
      IL = MAXVAL(PLON(:,:))
C
C     ALLOCATE SUM ARRAYS IN CELL SPACE.
C
      ALLOCATE( SUMA(IF:IL),
     +          SUMB(IF:IL),
     +          SUMR(IF:IL) )
C
C     PRINT A HEADER.
C
      WRITE(6,6000) 
     +  TRIM(CFILE),
     +  K,CELL,
     +  TRIM(CFILEG),
     +  TRIM(CFILEM)
C
C     DO ALL LATITUDE BAND SUMS.
C
      SUMT    = 0.0
      SUMA(:) = 0.0
      SUMB(:) = 0.0
      DO J= 1,JDM
        DO I= 1,IDM
          IF     (   A(I,J).NE.SPVAL .AND.
     +            AMSK(I,J).NE.SPVAL      ) THEN
            IZ = PLON(I,J)
            SUMA(IZ) = SUMA(IZ) + A(I,J)*PSC2(I,J)
            SUMB(IZ) = SUMB(IZ) +        PSC2(I,J)
            SUMT     = SUMT     +        PSC2(I,J)
          ENDIF
        ENDDO
      ENDDO
      DO IZ= IF,IL
        IF     (SUMB(IZ).NE.0.0) THEN
          SUMA(IZ) = SUMA(IZ)/SUMB(IZ)
        ENDIF
      ENDDO
C
C     DO ALL LATITUDE BAND SUM OF SQUARES.
C
      SUMR(:) = 0.0
      DO J= 1,JDM
        DO I= 1,IDM
          IF     (   A(I,J).NE.SPVAL .AND.
     +            AMSK(I,J).NE.SPVAL      ) THEN
            IZ = PLON(I,J)
            SUMR(IZ) = SUMR(IZ) + (A(I,J)-SUMA(IZ))**2 * PSC2(I,J)
          ENDIF
        ENDDO
      ENDDO
      DO IZ= IF,IL
        IF     (SUMB(IZ).NE.0.0) THEN
          SUMR(IZ) = SQRT(SUMR(IZ)/SUMB(IZ))
        ENDIF
        WRITE(6,'(F7.2,F10.5,1p2g16.8)') IZ*CELL,
     +                                   SUMB(IZ)/SUMT,
     +                                   SUMA(IZ),SUMR(IZ)
      ENDDO
      RETURN
 6000 FORMAT(
     +   '# hycom_meridional_lon file.a k cell grid.a mask.a' /
     +   '# hycom_meridional_lon ',A       /
     +   '#                 ',I5,F7.2 /
     +   '#                 ',A       /
     +   '#                 ',A       /
     +   '#   LON      FRAC        MEAN             RMS')
      END

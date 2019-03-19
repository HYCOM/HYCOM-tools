      PROGRAM OM_ISLANDS
      IMPLICIT NONE
C
C  hycom_islands - Usage:  hycom_islands depth.a [len [grid.a]]
C
C                 prints the location of islands in depth.a
C
C                 len is the maximum size of the islands (len by len),
C                 default is 3.
C
C                 grid.a is a hycom grid file, default regional.grid.a.
C                 Note that the corresponding grid.b must also exist.
C
C                 idm,jdm are taken from grid.a
C
C    This program detects an island by looking for NxN boxes (N<=len)
C    surrounded by sea and containing at least one land point.
C    Such points are definately on islands (no false positives),
C    but this algorithm is not guarenteed to find all islands.
C    It will find all single point islands and most "small" islands.
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  January 2001.
C
      LOGICAL, ALLOCATABLE :: L(:,:),M(:,:)
      REAL*4,  ALLOCATABLE :: D(:,:),PLON(:,:),PLAT(:,:)
      REAL*4               :: PAD(4096)
      INTEGER IOS
      INTEGER      IARGC
      INTEGER      NARG
      CHARACTER*240 CARG
C
      INTEGER       IDM,JDM,MXLEN,NPAD
      CHARACTER*6   CVARIN
      CHARACTER*240 CFILE,CFILEG,CFILEB
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.1) THEN
        CALL GETARG(1,CFILE)
        MXLEN  = 3
        CFILEG = 'regional.grid.a'
      ELSEIF (NARG.EQ.2) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,*) MXLEN
        CFILEG = 'regional.grid.a'
      ELSEIF (NARG.EQ.3) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,*) MXLEN
        CALL GETARG(3,CFILEG)
      ELSE
        WRITE(6,*) 
     +   'Usage:  hycom_islands depth.a [len [grid.a]]'
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
        WRITE(6,*) 'hycom_islands: bad header file ',
     &             CFILEB(1:LEN_TRIM(CFILEB))
        CALL EXIT(2)
      ENDIF
      READ( 11,*) JDM,CVARIN
      IF (CVARIN.NE.'jdm   ') THEN
        WRITE(6,*) 'hycom_islands: bad header file ',
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
      ALLOCATE( L(IDM+MXLEN+1,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_zonal: could not allocate ',
     +             (IDM+MXLEN+1)*JDM,' words for L'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( M(IDM+MXLEN+1,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_zonal: could not allocate ',
     +             (IDM+MXLEN+1)*JDM,' words for M'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( D(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_zonal: could not allocate ',
     +             IDM*JDM,' words for D'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( PLON(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_zonal: could not allocate ',
     +             IDM*JDM,' words for PLON'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( PLAT(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_zonal: could not allocate ',
     +             IDM*JDM,' words for PLAT'
        CALL EXIT(2)
      ENDIF
C
      CALL ISLAND(D,L,M,PLON,PLAT,IDM,JDM,PAD,NPAD,
     +            MXLEN, CFILE,CFILEG)
      CALL EXIT(0)
 5000 FORMAT(I4)
      END
      SUBROUTINE ISLAND(D,L,M,PLON,PLAT,IDM,JDM,PAD,NPAD,
     +                  MXLEN, CFILE,CFILEG)
      IMPLICIT NONE
C
      REAL*4     SPVAL
      PARAMETER (SPVAL=2.0**100)
C
      CHARACTER*240 CFILE,CFILEG
      INTEGER      IDM,JDM,NPAD,MXLEN
      LOGICAL      L(IDM+MXLEN+1,JDM),M(IDM+MXLEN+1,JDM)
      REAL*4       D(IDM,JDM),PLON(IDM,JDM),PLAT(IDM,JDM),PAD(NPAD)
C
C     MOST OF WORK IS DONE HERE.
C
#ifdef sun
      INTEGER      IR_ISNAN
C
#endif
      CHARACTER*18 CASN
      INTEGER      LEN_TRIM
      INTEGER      I,II,J,JJ,IOS,NRECL
      INTEGER      COUNTI,COUNTM,COUNTX,LEN,NISLE
      REAL*4       XP,YP
#ifdef CRAY
      INTEGER*8    IU8,IOS8
#endif
C
C     INPUT DEPTH ARRAY.
C
      INQUIRE( IOLENGTH=NRECL) D,PAD
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
      READ(11,REC=1,IOSTAT=IOS) D
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(D,IDM*JDM)
#endif
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'can''t read ',CFILE(1:LEN_TRIM(CFILE))
        CALL EXIT(4)
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
      READ(11,REC=1,IOSTAT=IOS) PLON
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(PLON,IDM*JDM)
#endif
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'can''t read ',CFILEG(1:LEN_TRIM(CFILEG))
        CALL EXIT(4)
      ENDIF
C
      READ(11,REC=2,IOSTAT=IOS) PLAT
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(PLAT,IDM*JDM)
#endif
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'can''t read ',CFILEG(1:LEN_TRIM(CFILEG))
        CALL EXIT(4)
      ENDIF
C
      CLOSE(UNIT=11)
C
C     FORM LAND MASK
C
      DO J= 1,JDM
        DO I= 1,IDM
          L(I,J) = D(I,J).EQ.SPVAL .OR. D(I,J).LE.0.0
          M(I,J) = .FALSE.
        ENDDO !i
        DO I= 1,MXLEN+1
          L(I+IDM,J) = L(I,J)
          M(I+IDM,J) = .FALSE.
        ENDDO !i
      ENDDO !j
C
C     FIND ISLANDS.
C
      WRITE(6,'(a,I4)')
     &  '# hycom_islands with len =',MXLEN
      WRITE(6,'(a)')
     &  '# isle     I     J       LON       LAT  SIZE'
C
      NISLE = 0
      DO LEN= 1,MXLEN
C
      DO J= LEN+2,JDM
        DO I= 1,IDM
          IF     (.NOT.L(I,J)) THEN
            IF     (COUNT(L(I+1:I+LEN+1,J)).EQ.0) THEN
C
C             NORTH IS ALL SEA
C
*             write(6,'(a,4i6)') 'NB',I,J
              IF     (COUNT(L(I:I+LEN+1,J-LEN-1)).EQ.0) THEN
C
C               SOUTH IS ALL SEA
C
                COUNTI = COUNT(L(I+1:I+LEN,J-LEN:J-1))
                COUNTX = COUNT(L(I:I+LEN+1,J-LEN:J-1))
*               write(6,'(a,4i6)') 'NB',I,J,COUNTI,COUNTX
                IF     (COUNTI.NE.0 .AND. COUNTI.EQ.COUNTX) THEN
C
C                 ISLAND, SINCE EAST AND WEST ARE ALL SEA
C
                  COUNTM = COUNT(M(I+1:I+LEN,J-LEN:J-1))
*                 write(6,'(a,4i6)') 'AB',I,J,COUNTM
                  IF     (COUNTI.NE.COUNTM) THEN
C
C                   NEW ISLAND.
C
                    NISLE = NISLE + 1
                    DO JJ= J-LEN,J-1
                      DO II= I+1,I+LEN
                        IF     (L(II,JJ) .AND. .NOT.M(II,JJ)) THEN
                          M(II,JJ) = .TRUE.
                          IF     (II.LE.IDM) THEN
                            YP =     PLAT(II,JJ)
                            XP = MOD(PLON(II,JJ)    +1440.0,360.0)
                          ELSE
                            YP =     PLAT(II-IDM,JJ)
                            XP = MOD(PLON(II-IDM,JJ)+1440.0,360.0)
                          ENDIF
                          IF     (XP.GT.180.0) THEN
                            XP = XP - 360.0
                          ENDIF
                          WRITE(6,'(3I6,2F10.3,I6)')
     &                      NISLE,II,JJ,XP,YP,COUNTI-COUNTM
                        ENDIF !land
                      ENDDO !ii
                    ENDDO !jj
                  ENDIF !new island
                ENDIF !island
              ENDIF !sea south
            ENDIF !sea north
          ENDIF !sea i,j
        ENDDO !i
      ENDDO !j
C
      ENDDO !len
C
      WRITE(6,'(a)')
     &  '# isle     I     J       LON       LAT  SIZE'
      WRITE(6,'(a,I4)')
     &  '# hycom_islands with len =',MXLEN
C
      RETURN
      END

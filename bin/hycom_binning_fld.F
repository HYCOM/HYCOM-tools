      PROGRAM HYCOM_BINNING_FLD
      IMPLICIT NONE
C
C  hycom_binning_fld - Usage:  hycom_binning_fld fin.a cell [grid.a] [mask.a] fout.a
C
C                 creates a HYCOM .[ab] file containing the array of means of
C                 lon/lat bins from all the (1:idm,1:jdm) arrays in fin.a
C                 linearly interpolated back to the original grid.
C
C                 use hycom_binning for plain text output (single array).
C                 use hycom_binning_nc  for netCDF output.
C
C                 longitude bands are centered on 0E and latitude bands
C                 are centered on the equator, both are nominally cell
C                 degrees wide, e.g. (N-0.5)*cell to (N+0.5)*cell.
C                 However, the cell sample width is (approximately) 
C                 scaled by 1/cos(lat) in longitude, i.e. the cell
C                 width in meters is aproximately constant.
C
C                 grid.a is a hycom grid file, default regional.grid.a.
C                 Note that the corresponding grid.b must also exist.
C
C                 mask.a contains an optional mask array, that turns
C                 off sampling anywhere it is 2.0**100 (data void).
C                 Note that grid.a must be specified if mask.a is used.
C
C                 idm,jdm are taken from grid.a and the input array
C                 is assumed to be on the p-grid (i.e. plat from grid.a
C                 represents the field's latitude).
C
C  fin.a is assumed to contain idm*jdm 32-bit IEEE real values for
C   each array, in standard f77 element order, followed by padding
C   to a multiple of 4096 32-bit words, but otherwise with no control
C   bytes/words, and input values of 2.0**100 indicating a data void.
C
C  this version for "serial" Unix systems.
C  Based on hycom_binning.F.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  August 2004.
C
      REAL*4, ALLOCATABLE :: A(:,:),AMSK(:,:),
     +                       PLAT(:,:),PLON(:,:),PSC2(:,:)
      REAL*4              :: PAD(4096)
      INTEGER IOS
      INTEGER      IARGC
      INTEGER      NARG
      CHARACTER*240 CARG
C
      INTEGER       IDM,JDM,NPAD
      REAL*4        CELL
      CHARACTER*6   CVARIN
      CHARACTER*240 CFILE,CFILEG,CFILEM,CFILEB,CFILEO
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.3) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,*)    CELL
        CFILEG = 'regional.grid.a'
        CFILEM = 'NONE'
        CALL GETARG(3,CFILEO)
      ELSEIF (NARG.EQ.4) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,*)    CELL
        CALL GETARG(3,CFILEG)
        CFILEM = 'NONE'
        CALL GETARG(4,CFILEO)
      ELSEIF (NARG.EQ.5) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,*)    CELL
        CALL GETARG(3,CFILEG)
        CALL GETARG(4,CFILEM)
        CALL GETARG(5,CFILEO)
      ELSE
        WRITE(6,*) 
     +   'Usage:  hycom_binning_fld fin.a cell [grid.a] [mask.a] fout.a'
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
        WRITE(6,*) 'hycom_binning_fld: bad header file ',
     &             CFILEB(1:LEN_TRIM(CFILEB))
        CALL EXIT(2)
      ENDIF
      READ( 11,*) JDM,CVARIN
      IF (CVARIN.NE.'jdm   ') THEN
        WRITE(6,*) 'hycom_binning_fld: bad header file ',
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
        WRITE(6,*) 'Error in hycom_binning_fld: could not allocate ',
     +             IDM*JDM,' words for A'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( AMSK(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_binning_fld: could not allocate ',
     +             IDM*JDM,' words for AMSK'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( PLAT(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_binning_fld: could not allocate ',
     +             IDM*JDM,' words for PLAT'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( PLON(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_binning_fld: could not allocate ',
     +             IDM*JDM,' words for PLON'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( PSC2(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_binning_fld: could not allocate ',
     +             IDM*JDM,' words for PSC2'
        CALL EXIT(2)
      ENDIF
C
      CALL BINNING(A,AMSK,PLAT,PLON,PSC2,IDM,JDM,PAD,NPAD,
     +             CELL, CFILE,CFILEG,CFILEM,CFILEO)
      CALL EXIT(0)
 5000 FORMAT(I4)
      END
      SUBROUTINE BINNING(A,AMSK,PLAT,PLON,PSC2,IDM,JDM, PAD,NPAD,
     +                   CELL, CFILE,CFILEG,CFILEM,CFILEO)
      IMPLICIT NONE
C
      REAL*4     SPVAL
      PARAMETER (SPVAL=2.0**100)
C
      CHARACTER*240 CFILE,CFILEG,CFILEM,CFILEO
      INTEGER      IDM,JDM,NPAD
      REAL*4       A(IDM,JDM),AMSK(IDM,JDM),
     +             PLAT(IDM,JDM),PLON(IDM,JDM),PSC2(IDM,JDM),PAD(NPAD),
     +             CELL
C
C     MOST OF WORK IS DONE HERE.
C
      REAL*8,  ALLOCATABLE :: SUMA(:,:),SUMB(:,:)
      REAL,    ALLOCATABLE :: AMEAN(:,:),ALON(:),ALAT(:)
      INTEGER, ALLOCATABLE :: IW(:)
C
#ifdef sun
      INTEGER      IR_ISNAN
C
#endif
      CHARACTER*18 CASN
      INTEGER      LEN_TRIM
      LOGICAL      LGLOBAL
      INTEGER      I,IC,IF,IL,IZ,IP,IQ,J,JC,JF,JL,JZ,IOS,K,NRECL
      REAL*8       SUMAT,SUMBT
      REAL*4       XLON,QD,AMX,AMN,DEG2RAD,MINCOS
      REAL*4       SX,SXM,SY,SYM
      REAL*4       LONMIN,LONMAX
#ifdef CRAY
      INTEGER*8    IU8,IOS8
#endif
C
C     INPUT ARRAY A.
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
      READ(11,REC=10,IOSTAT=IOS) PLAT  ! pscx
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(PLAT,IDM*JDM)
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
      PSC2(:,:) = PSC2(:,:)*PLAT(:,:)  ! psc2=pscx*pscy
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
      READ(11,REC= 2,IOSTAT=IOS) PLAT  ! plat
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
C    ASSUME THAT LONGITUDES ARE EITHER GLOBAL OR COMPACT
C
      LGLOBAL = PLON(IDM,1)-PLON(1,1).GT.350.0
      IF     (LGLOBAL) THEN
        LONMIN  = -0.5*CELL
        LONMAX  = (NINT(360.0/CELL)-1.5)*CELL 
      ELSE
        LONMIN  = MINVAL(PLON(:,:))
        LONMAX  = MAXVAL(PLON(:,:))
      ENDIF
C
C     MOVE PLON WITHIN +360 OF LONMIN, AND
C     CONVERT LONGITUDE TO CELL SPACE.
C
      DO J= 1,JDM
        DO I = 1,IDM
          XLON = PLON(I,J)
          IF     (XLON.LT.LONMIN) THEN
            DO
              XLON = XLON + 360.0
              IF     (XLON.GE.LONMIN) THEN
                EXIT
              ENDIF
            ENDDO
          ELSEIF (XLON.GE.LONMIN+360.0) THEN
            DO
              XLON = XLON - 360.0
              IF     (XLON.LT.LONMIN+360.0) THEN
                EXIT
              ENDIF
            ENDDO
          ENDIF
          PLON(I,J) = (XLON-LONMIN)/CELL + 0.5
        ENDDO !i
      ENDDO !j
      IF = NINT(MINVAL(PLON(:,:)))
      IL = NINT(MAXVAL(PLON(:,:)))
*     write(6,*) "IF,MINLON = ",IF,MINVAL(PLON(:,:))
*     write(6,*) "IL,MAXLON = ",IL,MAXVAL(PLON(:,:))
C
C     CONVERT LATITUDE TO CELL SPACE.
C
      PLAT(:,:) = PLAT(:,:)/CELL
      JF = NINT(MINVAL(PLAT(:,:)))
      JL = NINT(MAXVAL(PLAT(:,:)))
C
      IC = (IF+IL)/2
      JC = (JF+JL)/2
C
C     ALLOCATE SUM ARRAYS IN CELL SPACE.
C
      ALLOCATE(  SUMA(IF:IL,JF:JL),
     +           SUMB(IF:IL,JF:JL) )
      ALLOCATE( AMEAN(IF-1:IL+1,JF-1:JL+1) )
      ALLOCATE( ALON(IF:IL) )
      ALLOCATE( ALAT(JF:JL) )
      ALLOCATE(   IW(JF:JL) )
C
      AMEAN(:,:) = SPVAL
C
C     OUTPUT FILE.
C
      OPEN(UNIT=21, FILE=CFILEO, FORM='UNFORMATTED', STATUS='NEW',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',CFILEO(1:LEN_TRIM(CFILEO))
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
C
C     INPUT FILE
C
      OPEN(UNIT=11, FILE=CFILE, FORM='UNFORMATTED', STATUS='OLD',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',CFILE(1:LEN_TRIM(CFILE))
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
C
C     LONGITUDE CELL WIDTH (1+2*IW)
C
      DEG2RAD = 4.D0*ATAN(1.D0)/180.D0  !PI/180
      MINCOS  = 1.0/REAL(NINT(360.0/CELL)-1)
      DO JZ= JF,JL
        IW(JZ) = INT( 0.5/MAX(  MINCOS,
     &                         ABS(COS(JZ*CELL*DEG2RAD)) ) )
      ENDDO !jz
C
C     ALL RECORDS IN THE FILE
C
      DO K= 1,HUGE(K)
        READ(11,REC=K,IOSTAT=IOS) A
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
        IF     (K.EQ.1 .AND. IOS.NE.0) THEN
          WRITE(6,*) 'can''t read ',CFILE(1:LEN_TRIM(CFILE))
          CALL EXIT(4)
        ELSEIF (IOS.NE.0) THEN
          CLOSE(UNIT=11)
          EXIT
        ENDIF
C
C       DO ALL BINS.
C
        SUMAT     = 0.0
        SUMBT     = 0.0
        SUMA(:,:) = 0.0
        SUMB(:,:) = 0.0
        DO J= 1,JDM
          DO I= 1,IDM
            IF     (   A(I,J).NE.SPVAL .AND.
     +              AMSK(I,J).NE.SPVAL      ) THEN
              IZ = NINT(PLON(I,J))
              JZ = NINT(PLAT(I,J))
              IF     (IZ.LT.IF .OR. IZ.GT.IL .OR.
     +                JZ.LT.JF .OR. JZ.GT.JL     ) THEN
                CYCLE !i
              ENDIF
C
              SUMAT = SUMAT + A(I,J)*PSC2(I,J)
              SUMBT = SUMBT +        PSC2(I,J)
C
C             ALLOW FOR CELL WIDTH
C
              DO IP= IZ-IW(JZ),IZ+IW(JZ)
                IF     (LGLOBAL) THEN
                  IQ = MOD(IP,IL+1)
                ELSEIF (IP.LT.IF .OR. IP.GT.IL) THEN
                  CYCLE !ip
                ELSE
                  IQ = IP
                ENDIF
                SUMA(IQ,JZ) = SUMA(IQ,JZ) + A(I,J)*PSC2(I,J)
                SUMB(IQ,JZ) = SUMB(IQ,JZ) +        PSC2(I,J)
*
*               IF     (IQ.EQ.IC .AND. JZ.EQ.JC) THEN
*                 WRITE(6,'(a,2i5,1pg16.8)') 'debug: ',I,J,A(I,J)
*               ENDIF
              ENDDO !ip
            ENDIF
          ENDDO !j
        ENDDO !i
        DO JZ= JF,JL
          DO IZ= IF,IL
            IF     (SUMB(IZ,JZ).NE.0.0) THEN
              AMEAN(IZ,JZ) = SUMA(IZ,JZ)/SUMB(IZ,JZ)
            ELSE
              AMEAN(IZ,JZ) = SPVAL
            ENDIF
          ENDDO !iz
          IF     (LGLOBAL) THEN
            AMEAN(IF-1,JZ) = AMEAN(IL,JZ)  !periodic
            AMEAN(IL+1,JZ) = AMEAN(IF,JZ)  !periodic
          ENDIF
        ENDDO !jz
        DO IZ= IF,IL
          ALON(IZ) = IZ*CELL
        ENDDO !iz
        DO JZ= JF,JL
          ALAT(JZ) = JZ*CELL
        ENDDO !jz
C
        SUMAT = SUMAT/SUMBT
C
C       INTERPOLATE BACK TO THE ORIGINAL GRID
C
        AMN =  SPVAL
        AMX = -SPVAL
        DO J= 1,JDM
          DO I= 1,IDM
            IF     (   A(I,J).NE.SPVAL .AND.
     +              AMSK(I,J).NE.SPVAL      ) THEN
              IZ = PLON(I,J)
              IF     (PLON(I,J).LT.0.0) THEN
                IZ = IZ - 1
              ENDIF
              JZ = PLAT(I,J)
              IF     (PLAT(I,J).LT.0.0) THEN
                JZ = JZ - 1
              ENDIF
              SX  = PLON(I,J) - IZ
              SY  = PLAT(I,J) - JZ
              SXM = 1.0 - SX
              SYM = 1.0 - SY
              IF     (MAX(AMEAN(IZ,  JZ),  
     +                    AMEAN(IZ+1,JZ),
     +                    AMEAN(IZ,  JZ+1),
     +                    AMEAN(IZ+1,JZ+1) ).NE.SPVAL) THEN  !bilinear
                A(I,J) = SXM*SYM*AMEAN(IZ,  JZ)   +
     +                   SX *SYM*AMEAN(IZ+1,JZ)   +
     +                   SXM*SY *AMEAN(IZ,  JZ+1) +
     +                   SX *SY *AMEAN(IZ+1,JZ+1)
                AMX = MAX( AMX, A(I,J) )
                AMN = MIN( AMN, A(I,J) )
                if     (k.eq.0 .and. mod(i+j,100).eq.0) then
                  write(6,'(a,4i5,4f10.4)')
     +              'i,j,iz,jz,lon,lat,sx,sy    = ',
     +               i,j,iz,jz,plon(i,j),plat(i,j),sx,sy
                  write(6,'(a,4i5,5f10.4)')
     +              'i,j,iz,jz,am,a = ',
     +               i,j,iz,jz,
     +               amean(iz,jz),amean(iz+1,jz),
     +               amean(iz,jz+1),amean(iz+1,jz+1),
     +               a(i,j)
                endif
              ELSE  !distance weighted
                SUMAT = 0.0
                SUMBT = 0.0
                IF     (AMEAN(IZ,  JZ)  .NE.SPVAL) THEN
                  QD = 1.0 / MAX( 0.01, SX**2+SY**2 )
                  SUMAT = SUMAT + SQRT(QD)*AMEAN(IZ,  JZ)
                  SUMBT = SUMBT + SQRT(QD)
                  if     (k.eq.0 .and. mod(i+j,100).eq.0) then
                    write(6,'(a,4i5,5f10.4)')
     +                'i,j,ip,jp,d,qd,qd*am = ',
     +                 i,j,iz,  jz,
     +                 SX,SY,SX**2+SY**2,
     +                 qd,sqrt(qd)*amean(iz,  jz)
                  endif
                ENDIF
                IF     (AMEAN(IZ+1,JZ)  .NE.SPVAL) THEN
                  QD = 1.0 / MAX( 0.01, SXM**2+SY**2 )
                  SUMAT = SUMAT + SQRT(QD)*AMEAN(IZ+1,JZ)
                  SUMBT = SUMBT + SQRT(QD)
                  if     (k.eq.0 .and. mod(i+j,100).eq.0) then
                    write(6,'(a,4i5,5f10.4)')
     +                'i,j,ip,jp,d,qd,qd*am = ',
     +                 i,j,iz+1,jz,
     +                 SXM,SY,SXM**2+SY**2,
     +                 qd,sqrt(qd)*amean(iz+1,jz)
                  endif
                ENDIF
                IF     (AMEAN(IZ,  JZ+1).NE.SPVAL) THEN
                  QD = 1.0 / MAX( 0.01, SX**2+SYM**2 )
                  SUMAT = SUMAT + SQRT(QD)*AMEAN(IZ,  JZ+1)
                  SUMBT = SUMBT + SQRT(QD)
                  if     (k.eq.0 .and. mod(i+j,100).eq.0) then
                    write(6,'(a,4i5,5f10.4)')
     +                'i,j,ip,jp,d,qd,qd*am = ',
     +                 i,j,iz,  jz+1,
     +                 SX,SYM,SX**2+SYM**2,
     +                 qd,sqrt(qd)*amean(iz,  jz+1)
                  endif
                ENDIF
                IF     (AMEAN(IZ+1,JZ+1).NE.SPVAL) THEN
                  QD = 1.0 / MAX( 0.01, SXM**2+SYM**2 )
                  SUMAT = SUMAT + SQRT(QD)*AMEAN(IZ+1,JZ+1)
                  SUMBT = SUMBT + SQRT(QD)
                  if     (k.eq.0 .and. mod(i+j,100).eq.0) then
                    write(6,'(a,4i5,5f10.4)')
     +                'i,j,ip,jp,d,qd,qd*am = ',
     +                 i,j,iz+1,jz+1,
     +                 SXM,SYM,SXM**2+SYM**2,
     +                 qd,sqrt(qd)*amean(iz+1,jz+1)
                  endif
                ENDIF
                IF     (SUMBT.NE.0.0) THEN
                  A(I,J) = SUMAT/SUMBT
                  AMX = MAX( AMX, A(I,J) )
                  AMN = MIN( AMN, A(I,J) )
                  if     (k.eq.0 .and. mod(i+j,100).eq.0) then
                    write(6,'(a,4i5,5f10.4)')
     +                'i,j,iz,jz,lon,lat,sx,sy,bt = ',
     +                 i,j,iz,jz,plon(i,j),plat(i,j),sx,sy,
     +                 sumbt
                    write(6,'(a,4i5,5f10.4)')
     +                'i,j,iz,jz,am,a = ',
     +                 i,j,iz,jz,
     +                 amean(iz,jz),amean(iz+1,jz),
     +                 amean(iz,jz+1),amean(iz+1,jz+1),
     +                 a(i,j)
                endif
                ELSE
                  A(I,J) = SPVAL
                ENDIF
              ENDIF
            ELSE
              A(I,J) = SPVAL
            ENDIF
          ENDDO !j
        ENDDO !i
C
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
        IF     (NPAD.EQ.0) THEN
          WRITE(21,REC=K,IOSTAT=IOS) A
        ELSE
          WRITE(21,REC=K,IOSTAT=IOS) A,PAD
        ENDIF
        WRITE(6,'(a,1p2g16.8)')
     &     'min, max = ',AMN,AMX
      ENDDO !k
      WRITE(6,*)
      WRITE(6,*) K-1,' FIELDS PROCESSED'
      WRITE(6,*)
C
      CLOSE(11)
      CLOSE(21)
C
      RETURN
 6000 FORMAT(
     +   '# hycom_binning_fld fin.a cell grid.a mask.a fout.a' /
     +   '# hycom_binning_fld ',A       /
     +   '#               ',F7.2 /
     +   '#               ',A       /
     +   '#               ',A       /
     +   '#    LON    LAT        MEAN')
      END

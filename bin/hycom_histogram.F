      PROGRAM HYCOM_HISTOGRAM
      IMPLICIT NONE
C
C  hycom_histogram - Usage:  hycom_histogram file.a k hmin hmax hinc [lonmin lonmax latmin latmax] [grid.a [mask.a]]
C
C                 prints the percentage of the domain where the k-th
C                 (1:idm,1:jdm) array in file.a lies in each histogram bin
C
C                 also provides approximate 50, 95 and 99 percentils
C
C                 hmin and hmax define the histogram extent and 
C                 hinc defines the histogram bin size.
C
C                 lonmin lonmax latmin latmax define the extent of
C                 the sampled subregion, default the entire region.
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
      REAL*4, ALLOCATABLE :: A(:,:),AMSK(:,:),
     +                       PLAT(:,:),PLON(:,:),PSC2(:,:)
      REAL*4              :: PAD(4096)
      INTEGER IOS
      INTEGER      IARGC
      INTEGER      NARG
      CHARACTER*240 CARG
C
      INTEGER       IDM,JDM,K,NPAD
      REAL*4        HMIN,HMAX,HINC,LONMIN,LONMAX,LATMIN,LATMAX
      CHARACTER*6   CVARIN
      CHARACTER*240 CFILE,CFILEG,CFILEM,CFILEB
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.5) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,*) K
        CALL GETARG(3,CARG)
        READ(CARG,*)    HMIN
        CALL GETARG(4,CARG)
        READ(CARG,*)    HMAX
        CALL GETARG(5,CARG)
        READ(CARG,*)    HINC
        CFILEG = 'regional.grid.a'
        CFILEM = 'NONE'
        LATMAX = 999.0  ! indicates entire region
      ELSEIF (NARG.EQ.6) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,*) K
        CALL GETARG(3,CARG)
        READ(CARG,*)    HMIN
        CALL GETARG(4,CARG)
        READ(CARG,*)    HMAX
        CALL GETARG(5,CARG)
        READ(CARG,*)    HINC
        CALL GETARG(6,CFILEG)
        CFILEM = 'NONE'
        LATMAX = 999.0  ! indicates entire region
      ELSEIF (NARG.EQ.7) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,*) K
        CALL GETARG(3,CARG)
        READ(CARG,*)    HMIN
        CALL GETARG(4,CARG)
        READ(CARG,*)    HMAX
        CALL GETARG(5,CARG)
        READ(CARG,*)    HINC
        CALL GETARG(6,CFILEG)
        CALL GETARG(7,CFILEM)
        LATMAX = 999.0  ! indicates entire region
      ELSEIF (NARG.EQ.9) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,*) K
        CALL GETARG(3,CARG)
        READ(CARG,*)    HMIN
        CALL GETARG(4,CARG)
        READ(CARG,*)    HMAX
        CALL GETARG(5,CARG)
        READ(CARG,*)    HINC
        CALL GETARG(6,CARG)
        READ(CARG,*)    LONMIN
        CALL GETARG(7,CARG)
        READ(CARG,*)    LONMAX
        CALL GETARG(8,CARG)
        READ(CARG,*)    LATMIN
        CALL GETARG(9,CARG)
        READ(CARG,*)    LATMAX
        CFILEG = 'regional.grid.a'
        CFILEM = 'NONE'
      ELSEIF (NARG.EQ.10) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,*) K
        CALL GETARG(3,CARG)
        READ(CARG,*)    HMIN
        CALL GETARG(4,CARG)
        READ(CARG,*)    HMAX
        CALL GETARG(5,CARG)
        READ(CARG,*)    HINC
        CALL GETARG(6,CARG)
        READ(CARG,*)    LONMIN
        CALL GETARG(7,CARG)
        READ(CARG,*)    LONMAX
        CALL GETARG(8,CARG)
        READ(CARG,*)    LATMIN
        CALL GETARG(9,CARG)
        READ(CARG,*)    LATMAX
        CALL GETARG(10,CFILEG)
        CFILEM = 'NONE'
      ELSEIF (NARG.EQ.11) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,*) K
        CALL GETARG(3,CARG)
        READ(CARG,*)    HMIN
        CALL GETARG(4,CARG)
        READ(CARG,*)    HMAX
        CALL GETARG(5,CARG)
        READ(CARG,*)    HINC
        CALL GETARG(6,CARG)
        READ(CARG,*)    LONMIN
        CALL GETARG(7,CARG)
        READ(CARG,*)    LONMAX
        CALL GETARG(8,CARG)
        READ(CARG,*)    LATMIN
        CALL GETARG(9,CARG)
        READ(CARG,*)    LATMAX
        CALL GETARG(10,CFILEG)
        CALL GETARG(11,CFILEM)
      ELSE
        WRITE(6,*) 
     +   'Usage:  hycom_histogram file.a k hmin hmax hinc ' //
     +   '[lonmin lonmax latmin latmax] [grid.a [mask.a]]'
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
        WRITE(6,*) 'hycom_histogram: bad header file ',
     &             CFILEB(1:LEN_TRIM(CFILEB))
        CALL EXIT(2)
      ENDIF
      READ( 11,*) JDM,CVARIN
      IF (CVARIN.NE.'jdm   ') THEN
        WRITE(6,*) 'hycom_histogram: bad header file ',
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
        WRITE(6,*) 'Error in hycom_zonal: could not allocate ',
     +             IDM*JDM,' words for A'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( AMSK(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_zonal: could not allocate ',
     +             IDM*JDM,' words for AMSK'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( PLAT(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_zonal: could not allocate ',
     +             IDM*JDM,' words for PLAT'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( PLON(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_zonal: could not allocate ',
     +             IDM*JDM,' words for PLON'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( PSC2(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_zonal: could not allocate ',
     +             IDM*JDM,' words for PSC2'
        CALL EXIT(2)
      ENDIF
C
      CALL HISTOGRAM(A,AMSK,PLAT,PLON,PSC2,IDM,JDM,PAD,NPAD,
     +               K,HMIN,HMAX,HINC,
     +               LONMIN,LONMAX,LATMIN,LATMAX,
     +               CFILE,CFILEG,CFILEM)
      CALL EXIT(0)
      END
      SUBROUTINE HISTOGRAM(A,AMSK,PLAT,PLON,PSC2,IDM,JDM, PAD,NPAD,
     +                     K,HMIN,HMAX,HINC,
     +                     LONMIN,LONMAX,LATMIN,LATMAX,
     +                     CFILE,CFILEG,CFILEM)
      IMPLICIT NONE
C
      REAL*4     SPVAL
      PARAMETER (SPVAL=2.0**100)
C
      CHARACTER*240 CFILE,CFILEG,CFILEM
      INTEGER      IDM,JDM,NPAD,K
      REAL*4       A(IDM,JDM),AMSK(IDM,JDM),
     +             PLAT(IDM,JDM),PLON(IDM,JDM),PSC2(IDM,JDM),PAD(NPAD),
     +             HMIN,HMAX,HINC,
     +             LONMIN,LONMAX,LATMIN,LATMAX
C
C     MOST OF WORK IS DONE HERE.
C
      REAL*8, ALLOCATABLE :: SUMH(:)
C
#ifdef sun
      INTEGER      IR_ISNAN
C
#endif
      CHARACTER*18 CASN
      INTEGER      LEN_TRIM
      LOGICAL      OUT_OF_LON_BAND
      LOGICAL      LSUBREG
      INTEGER      I,J,JF,JL,JZ,IOS,NRECL
      REAL*4       Q,HLARGE
      REAL*8       SUMHT,SUMPC,SUMPCM,PCTILE(3)
#ifdef CRAY
      INTEGER*8    IU8,IOS8
#endif
C
C     SUBREGION?
C
      LSUBREG = LATMAX.LT.91.0  ! sub-region specified
C
      IF     (LSUBREG) THEN
C
C       MOVE LONMAX WITHIN +360 OF LONMIN
C
        IF     (LONMAX.LT.LONMIN) THEN
          DO
            LONMAX = LONMAX + 360.0
            IF     (LONMAX.GT.LONMIN) THEN
              EXIT
            ENDIF
          ENDDO
        ELSEIF (LONMAX.GE.LONMIN+360.0) THEN
          DO
            LONMAX = LONMAX - 360.0
            IF     (LONMAX.LT.LONMIN+360.0) THEN
              EXIT
            ENDIF
          ENDDO
        ENDIF
      ENDIF !lsubreg
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
C     ALLOCATE SUM ARRAYS IN HISTOGRAM SPACE.
C
      JF = -1
      JL = INT((HMAX-HMIN+HINC*0.5)/HINC)
      ALLOCATE( SUMH(JF:JL) )
C
C     PRINT A HEADER.
C
      IF     (.NOT. LSUBREG) THEN ! full region
        WRITE(6,6000) 
     +    TRIM(CFILE),
     +    K,HMIN,HMAX,HINC,
     +    TRIM(CFILEG),
     +    TRIM(CFILEM)
      ELSE
        WRITE(6,6100) 
     +    TRIM(CFILE),
     +    K,HMIN,HMAX,HINC,
     +    LONMIN,LONMAX,LATMIN,LATMAX,
     +    TRIM(CFILEG),
     +    TRIM(CFILEM)
      ENDIF
C
C     DO ALL LATITUDE BAND SUMS.
C
      HLARGE  = HMAX + 2.0*HINC  !guard against integer overflow
      SUMHT   = 0.0
      SUMH(:) = 0.0
      DO J= 1,JDM
        DO I= 1,IDM
          IF     (   A(I,J).NE.SPVAL .AND.
     +            AMSK(I,J).NE.SPVAL      ) THEN
            IF     (LSUBREG) THEN
              IF     (PLAT(I,J).LT.LATMIN .OR.
     +                PLAT(I,J).GT.LATMAX     ) THEN
                CYCLE
              ENDIF
              IF     (OUT_OF_LON_BAND(PLON(I,J), LONMIN,LONMAX)) THEN
                CYCLE
              ENDIF
            ENDIF
            JZ = MAX(JF,
     &               MIN( JL,
     &                    INT((MIN(A(I,J),HLARGE)-HMIN)/HINC) ) )
*           if     (jz.eq.jf) then
*             WRITE(6,*) 'I,J,A,JZ = ',I,J,A(I,J),JZ
*           endif
C
            SUMH(JZ) = SUMH(JZ) + PSC2(I,J)
            SUMHT    = SUMHT    + PSC2(I,J)
          ENDIF
        ENDDO
      ENDDO
      JZ = JF
      SUMH(JZ) = 100.0*SUMH(JZ)/SUMHT
      SUMPC    = SUMH(JZ)
      IF     (SUMH(JZ).NE.0.0) THEN
        WRITE(6,'(A11,2F11.4,2F12.5)')
     +    "-Inf",HMIN+(JZ+1)*HINC,HMIN+(JZ+0.5)*HINC,
     +    SUMH(JZ),SUMPC
      ENDIF
      DO JZ= JF+1,JL-1
        SUMH(JZ) = 100.0*SUMH(JZ)/SUMHT
        SUMPC = SUMPC + SUMH(JZ)
        WRITE(6,'(3F11.4,2F12.5)')
     +    HMIN+JZ*HINC,HMIN+(JZ+1)*HINC,HMIN+(JZ+0.5)*HINC,
     +    SUMH(JZ),SUMPC
      ENDDO
      JZ = JL
      IF     (SUMH(JZ).NE.0.0) THEN
        SUMH(JZ) = 100.0*SUMH(JZ)/SUMHT
        SUMPC    = 100.0
        WRITE(6,'(F11.4,A11,F11.4,2F12.5)')
     +    HMIN+JZ*HINC,"+Inf",HMIN+(JZ+0.5)*HINC,
     +    SUMH(JZ),SUMPC
      ENDIF
C
C     PERCENTILS
C
      PCTILE(1) = 50.0
      PCTILE(2) = 95.0
      PCTILE(3) = 99.0
      DO I= 1,3
        SUMPC  =  0.0
        DO JZ = JF,JL-1
          SUMPCM = SUMPC
          SUMPC  = SUMPC + SUMH(JZ)
          IF     (SUMPC.GE.PCTILE(I)) THEN
            Q = (SUMPC-PCTILE(I))/(SUMPC-SUMPCM)
            WRITE(6,'(A,F6.1,A,F11.4)')
     &        '## ',PCTILE(I),' percentil =',HMIN+JZ*HINC+(1.0-Q)*HINC
            EXIT
          ENDIF
        ENDDO !jz
      ENDDO !i
      RETURN
 6000 FORMAT(
     +   '# hycom_histogram file.a k hmin hmax hinc grid.a mask.a' /
     +   '# hycom_histogram ',A       /
     +   '#                 ',I5,3F11.4 /
     +   '#                 ',A       /
     +   '#                 ',A       /
     +   '#       MIN        MAX    TYPICAL     PERCENT       TOTAL')
 6100 FORMAT(
     +   '# hycom_histogram file.a k hmin hmax hinc ',
     +   'lonmin lonmax latmin latmax grid.a mask.a' /
     +   '# hycom_histogram ',A       /
     +   '#                 ',I5,3F11.4 /
     +   '#                 ',2F9.2,2F8.2 /
     +   '#                 ',A       /
     +   '#                 ',A       /
     +   '#       MIN        MAX    TYPICAL     PERCENT       TOTAL')
      END
      LOGICAL FUNCTION OUT_OF_LON_BAND(LON, LONMIN,LONMAX)
      IMPLICIT NONE
      REAL*4 LON, LONMIN,LONMAX
C
C     IS LON OUTSIDE THE RANGE LONMIN TO LONMAX?
C
      REAL*4 XLON
C
      XLON = LON
      IF     (LON.LT.LONMIN) THEN
        DO
          XLON = XLON + 360.0
          IF     (XLON.GT.LONMIN) THEN
            EXIT
          ENDIF
        ENDDO
      ELSEIF (LON.GE.LONMIN+360.0) THEN
        DO
          XLON = XLON - 360.0
          IF     (XLON.LT.LONMIN+360.0) THEN
            EXIT
          ENDIF
        ENDDO
      ENDIF
C
      OUT_OF_LON_BAND = XLON .GT. LONMAX
      END

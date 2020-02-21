      PROGRAM TRACK_HISTOGRAM
      IMPLICIT NONE
C
C  track_histogram - Usage:  track_histogram file.txt hmin hmax hinc [lonmin lonmax latmin latmax] [latlon]
C
C                 prints the percentage of the domain where the
C                 locations in file.txt lies in each histogram bin
C
C                 also provides approximate 50, 95 and 99 percentils
C
C                 hmin and hmax define the histogram extent and 
C                 hinc defines the histogram bin size.
C
C                 lonmin lonmax latmin latmax define the extent of
C                 the sampled subregion, default the entire region.
C
C                 latlon is a marker
C
C
C  file.txt is assumed to lines that are blank, or start with "#",
C  or are of the form: longitude latitude value
C
C  if latlon is present, the lines are: latitude longitude value
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  November 2015.
C
      REAL*4        A,PLAT,PLON
      INTEGER       IOS
      INTEGER       IARGC
      INTEGER       NARG
      CHARACTER*240 CARG
C
      REAL*4        HMIN,HMAX,HINC,LONMIN,LONMAX,LATMIN,LATMAX
      CHARACTER*6   CVARIN
      CHARACTER*240 CFILE
      LOGICAL       LATLON
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.4) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,*)    HMIN
        CALL GETARG(3,CARG)
        READ(CARG,*)    HMAX
        CALL GETARG(4,CARG)
        READ(CARG,*)    HINC
        LATMAX = 999.0  ! indicates entire region
        LATLON = .FALSE.
      ELSEIF (NARG.EQ.5) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,*)    HMIN
        CALL GETARG(3,CARG)
        READ(CARG,*)    HMAX
        CALL GETARG(4,CARG)
        READ(CARG,*)    HINC
        LATMAX = 999.0  ! indicates entire region
        LATLON = .TRUE.
      ELSEIF (NARG.EQ.9) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,*)    HMIN
        CALL GETARG(3,CARG)
        READ(CARG,*)    HMAX
        CALL GETARG(4,CARG)
        READ(CARG,*)    HINC
        CALL GETARG(5,CARG)
        READ(CARG,*)    LONMIN
        CALL GETARG(6,CARG)
        READ(CARG,*)    LONMAX
        CALL GETARG(7,CARG)
        READ(CARG,*)    LATMIN
        CALL GETARG(8,CARG)
        READ(CARG,*)    LATMAX
        LATLON = .FALSE.
      ELSEIF (NARG.EQ.10) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,*)    HMIN
        CALL GETARG(3,CARG)
        READ(CARG,*)    HMAX
        CALL GETARG(4,CARG)
        READ(CARG,*)    HINC
        CALL GETARG(5,CARG)
        READ(CARG,*)    LONMIN
        CALL GETARG(6,CARG)
        READ(CARG,*)    LONMAX
        CALL GETARG(7,CARG)
        READ(CARG,*)    LATMIN
        CALL GETARG(8,CARG)
        READ(CARG,*)    LATMAX
        LATLON = .TRUE.
      ELSE
        WRITE(6,*) 
     +   'Usage:  track_histogram file.txt hmin hmax hinc ' //
     +   '[lonmin lonmax latmin latmax] [latlon]'
        CALL EXIT(1)
      ENDIF
C
      CALL HISTOGRAM(A,PLAT,PLON,
     +               HMIN,HMAX,HINC,
     +               LONMIN,LONMAX,LATMIN,LATMAX,LATLON,
     +               CFILE)
      CALL EXIT(0)
      END
      SUBROUTINE HISTOGRAM(A,PLAT,PLON,
     +                     HMIN,HMAX,HINC,
     +                     LONMIN,LONMAX,LATMIN,LATMAX,LATLON,
     +                     CFILE)
      IMPLICIT NONE
C
      CHARACTER*240 CFILE
      LOGICAL       LATLON
      REAL*4        A,PLAT,PLON,
     +              HMIN,HMAX,HINC,
     +              LONMIN,LONMAX,LATMIN,LATMAX
C
C     MOST OF WORK IS DONE HERE.
C
      REAL*8, ALLOCATABLE :: SUMH(:)
C
#ifdef sun
      INTEGER      IR_ISNAN
C
#endif
      CHARACTER*240 CLINE
      LOGICAL       OUT_OF_LON_BAND
      LOGICAL       LSUBREG
      INTEGER       I,J,JF,JL,JZ,IOS,NRECL
      REAL*4        Q,HLARGE
      REAL*8        SUMHT,SUMPC,SUMPCM,PCTILE(3)
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
C     ALLOCATE SUM ARRAYS IN HISTOGRAM SPACE.
C
      JF = -1
      JL = INT((HMAX-HMIN+HINC*0.5)/HINC)
      ALLOCATE( SUMH(JF:JL) )
C
C     PRINT A HEADER.
C
      IF     (LATLON) THEN
        IF     (.NOT. LSUBREG) THEN ! full globe
          WRITE(6,6010) 
     +      TRIM(CFILE),
     +      HMIN,HMAX,HINC
        ELSE
          WRITE(6,6110) 
     +      TRIM(CFILE),
     +      HMIN,HMAX,HINC,
     +      LONMIN,LONMAX,LATMIN,LATMAX
        ENDIF
      ELSE
        IF     (.NOT. LSUBREG) THEN ! full globe
          WRITE(6,6000) 
     +      TRIM(CFILE),
     +      HMIN,HMAX,HINC
        ELSE
          WRITE(6,6100) 
     +      TRIM(CFILE),
     +      HMIN,HMAX,HINC,
     +      LONMIN,LONMAX,LATMIN,LATMAX
        ENDIF
      ENDIF !latlon
C
C     OPEN INPUT TRACKS
C
      OPEN(UNIT=11, FILE=CFILE, FORM='FORMATTED', STATUS='OLD',
     +         IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILE)
        write(6,*) 'ios   = ',ios
        CALL EXIT(3)
      ENDIF
C
C     DO ALL LATITUDE BAND SUMS.
C
      HLARGE  = HMAX + 2.0*HINC  !guard against integer overflow
      SUMHT   = 0.0
      SUMH(:) = 0.0
      DO I= 1,HUGE(I)
        READ(11,'(A)',IOSTAT=IOS) CLINE
        IF     (IOS.NE.0) THEN
          EXIT
        ENDIF
        IF     (CLINE     .EQ.' ' .OR.
     +          CLINE(1:1).EQ.'#'     ) THEN
          CYCLE
        ENDIF
        IF     (LATLON) THEN
          READ(CLINE,*) PLAT,PLON,A
        ELSE
          READ(CLINE,*) PLON,PLAT,A
        ENDIF
C
        IF     (LSUBREG) THEN
          IF     (PLAT.LT.LATMIN .OR.
     +            PLAT.GT.LATMAX     ) THEN
            CYCLE
          ENDIF
          IF     (OUT_OF_LON_BAND(PLON, LONMIN,LONMAX)) THEN
            CYCLE
          ENDIF
        ENDIF
        JZ = MAX(JF,
     &           MIN( JL,
     &                INT((MIN(A,HLARGE)-HMIN)/HINC) ) )
*           if     (jz.eq.jf) then
*             WRITE(6,*) 'I,A,JZ = ',I,A,JZ
*           endif
C
        SUMH(JZ) = SUMH(JZ) + 1.0
        SUMHT    = SUMHT    + 1.0
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
     +   '# track_histogram file.txt hmin hmax hinc' /
     +   '# track_histogram ',A       /
     +   '#                 ',3F11.4 /
     +   '#       MIN        MAX    TYPICAL     PERCENT       TOTAL')
 6010 FORMAT(
     +   '# track_histogram file.txt hmin hmax hinc latlon' /
     +   '# track_histogram ',A       /
     +   '#                 ',3F11.4 /
     +   '#       MIN        MAX    TYPICAL     PERCENT       TOTAL')
 6100 FORMAT(
     +   '# track_histogram file.txt hmin hmax hinc ',
     +   'lonmin lonmax latmin latmax' /
     +   '# track_histogram ',A       /
     +   '#                 ',3F11.4 /
     +   '#                 ',2F9.2,2F8.2 /
     +   '#       MIN        MAX    TYPICAL     PERCENT       TOTAL')
 6110 FORMAT(
     +   '# track_histogram file.txt hmin hmax hinc ',
     +   'lonmin lonmax latmin latmax latlon' /
     +   '# track_histogram ',A       /
     +   '#                 ',3F11.4 /
     +   '#                 ',2F9.2,2F8.2 /
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

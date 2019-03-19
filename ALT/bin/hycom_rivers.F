      PROGRAM LOC_RIVERS
      IMPLICIT NONE
C
C  hycom_rivers - Usage:  hycom_rivers f_lon l_lon f_lat l_lat [rivers.d]
C
C                 prints the location and monthly flow of all rivers
C                 in rivers.d that are within the specifed lon/lat box.
C 
C                 see README.ALL.bin.rivers for the format of rivers.d.
C                 rivers.d defaults to {command_name}.d
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  October 2001.
C
      INTEGER      IARGC
      INTEGER      NARG
      CHARACTER*240 CARG
C
      INTEGER      IOS
      REAL         ALAT(2),ALON(2),RFLOW,RLAT,RLON,RLONIN
      CHARACTER*240 CFILE,CRIVR,CFLOW,CTEMP,CEXEC
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.5) THEN
        CALL GETARG(1,CARG)
        READ(CARG,*) ALON(1)
        CALL GETARG(2,CARG)
        READ(CARG,*) ALON(2)
        CALL GETARG(3,CARG)
        READ(CARG,*) ALAT(1)
        CALL GETARG(4,CARG)
        READ(CARG,*) ALAT(2)
        CALL GETARG(5,CFILE)
      ELSEIF (NARG.EQ.4) THEN
        CALL GETARG(1,CARG)
        READ(CARG,*) ALON(1)
        CALL GETARG(2,CARG)
        READ(CARG,*) ALON(2)
        CALL GETARG(3,CARG)
        READ(CARG,*) ALAT(1)
        CALL GETARG(4,CARG)
        READ(CARG,*) ALAT(2)
        CALL GETARG(0,CEXEC)
        CFILE = CEXEC(1:LEN_TRIM(CEXEC)) // '.d'
      ELSE
        WRITE(6,*)
     &    'Usage: hycom_rivers f_lon l_lon f_lat l_lat [rivers.d]'
        CALL EXIT(1)
      ENDIF
C
C     CONVERT TO DEGREES BETWEEN 0E AND 720.
C
      IF     (ALON(1).LT.0.0) THEN
        ALON(1) = ALON(1) + 360.0
      ENDIF
      IF     (ALON(2).LT.ALON(1)) THEN
        ALON(2) = ALON(2) + 360.0
      ENDIF
      IF     (ALON(2).LT.ALON(1)) THEN
        ALON(2) = ALON(2) + 360.0
      ENDIF
C
      OPEN(UNIT=11, FILE=CFILE, FORM='FORMATTED', STATUS='OLD',
     +     IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error: can''t open ',CFILE(1:LEN_TRIM(CFILE))
        WRITE(6,*) 'ios   = ',ios
        CALL EXIT(3)
      ENDIF
      DO
        READ(11,'(a)',IOSTAT=IOS) CRIVR
        IF     (IOS.NE.0) THEN
          EXIT
        ENDIF
        READ(11,'(a)')            CFLOW
        READ(11,'(a)')            CTEMP
        READ(CRIVR,*) RFLOW,RLON,RLAT
        RLONIN = RLON
C
C       IS THIS RIVER IN THE BOX?
C
        IF     (ALAT(1).LE.RLAT .AND. RLAT.LE.ALAT(2)) THEN
          IF     (RLON.LT.ALON(1)) THEN
            RLON = RLON + 360.0
          ENDIF
          IF     (RLON.LT.ALON(1)) THEN
            RLON = RLON + 360.0
          ENDIF
          IF     (ALON(1).LE.RLON .AND. RLON.LE.ALON(2)) THEN
            IF     (RLONIN.GT.180.0) THEN
              WRITE(CRIVR(12:20),'(F9.4)') RLONIN-360.0 !180:360 to -180:0
            ENDIF
            WRITE(6,'(a)')     CRIVR( 1:LEN_TRIM(CRIVR))
            WRITE(6,'(40x,a)') CFLOW( 1: 60)
            WRITE(6,'(40x,a)') CFLOW(61:120)
          ENDIF
        ENDIF
      ENDDO
      CALL EXIT(0)
      END

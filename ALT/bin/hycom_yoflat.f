      PROGRAM FNDLAT
      IMPLICIT NONE
C
C**********
C*
C 1)  MAPS HYCOM GRID TO/FROM LATITUDE
C
C 2)  USES SUN FORTRAN SPECIFIC CARRIAGE CONTROL EDIT DESCIPTOR ('$').
C
C 3)  ALAN J. WALLCRAFT,  NRL,  MAY 2000.
C*
C**********
C
      INTEGER JDM,LOOP,MAPTYPE
      REAL    GRIDSZ,PI,PLAT,PY,QLAT,RADIAN,Y
C
      REAL ALAT,ALAT1,ALATU,DIST,DIST1,DISTU,GRID
      ALAT( DIST1,GRID)=(2.*ATAN(EXP(DIST1*GRID/RADIAN))-PI/2.)
      DIST( ALAT1,GRID)=ALOG(TAN((2.*ALAT1+PI)/4.))*RADIAN/GRID
      ALATU(DIST1,GRID)=DIST1*(GRID/RADIAN)
      DISTU(ALAT1,GRID)=ALAT1*(RADIAN/GRID)
C
      RADIAN=57.2957795
      PI    = 3.1415926536
C
C     SELECT GRID SIZE AND MAPPING DIRECTION.
C
      WRITE(6,'(a)',ADVANCE='NO')
     &   ' Grid size (deg) : '
      READ( 5,*)  GRIDSZ
      WRITE(6,'(a)',ADVANCE='NO')
     &   ' Eq. is at P(?)  : '
      READ(5,*) PY
      WRITE(6,'(a)',ADVANCE='NO')
     &   ' Input latitude (0) or grid (1) : '
      READ( 5,*)  MAPTYPE
C
      DO 11 LOOP=1,99999
C
        IF     (MAPTYPE.EQ.0 .OR. MAPTYPE.EQ.2) THEN
C
C         INPUT LATITUDE.
C
          WRITE(6,'(a)',ADVANCE='NO')
     &       ' Latitude (degN, -99.0 to exit) : '
          READ(5,*) PLAT
          IF     (ABS(PLAT).GT.90.0) THEN
            GOTO 111
          ENDIF
C
          IF     (MAPTYPE.EQ.0) THEN
            Y = PY + DIST( PLAT/RADIAN,GRIDSZ)
          ELSE
            Y = PY + DISTU(PLAT/RADIAN,GRIDSZ)
          ENDIF
          WRITE(6,7100) PLAT,Y
 7100     FORMAT(' Lat =',F8.3,
     +           '  is at P(',F9.3,')' /)
        ELSEIF (MAPTYPE.EQ.1 .OR. MAPTYPE.EQ.3) THEN
C
C         INPUT GRID DISTANCE.
C
          WRITE(6,'(a)',ADVANCE='NO')
     &       ' Index, -1.0 to exit : '
          READ(5,*) Y
          IF     (Y.LE.-1.0) THEN
            GOTO 111
          ENDIF
C
          IF     (MAPTYPE.EQ.1) THEN
            PLAT = ALAT( Y-PY,GRIDSZ)*RADIAN
          ELSE
            PLAT = ALATU(Y-PY,GRIDSZ)*RADIAN
          ENDIF
          WRITE(6,8100) Y,PLAT
 8100     FORMAT(' P(',F9.3,
     +           '), is at Lat =',F8.3 /)
        ENDIF
C
   11 CONTINUE
  111 CONTINUE
      STOP
C     END OF PROGRAM FNDLAT.
      END
      SUBROUTINE IEEE_RETROSPECTIVE()
C
C     DUMMY ROUTINE TO TURN OFF IEEE WARNING MESSAGES ON A SUN.
C
      END

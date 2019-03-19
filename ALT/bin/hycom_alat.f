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
     &   ' Grid size at equator (degrees) : '
      READ( 5,*)  GRIDSZ
      WRITE(6,'(a)',ADVANCE='NO')
     &   ' Latitudinal array extent (jdm) : '
      READ( 5,*)  JDM
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
     &       ' Latitude of P(1) (degN, -99.0 to exit) : '
          READ(5,*) PLAT
          IF     (ABS(PLAT).GT.90.0) THEN
            GOTO 111
          ENDIF
C
          IF     (MAPTYPE.EQ.0) THEN
            Y = DIST( PLAT/RADIAN,GRIDSZ)
          ELSE
            Y = DISTU(PLAT/RADIAN,GRIDSZ)
          ENDIF
          PY = 1.0-Y
          IF     (JDM.LE.1) THEN
            WRITE(6,7100) PLAT,PY
 7100       FORMAT(' P(1) Lat =',F8.3,
     +             '  Eq. is at P(',F9.3,')' /)
          ELSE
            IF     (MAPTYPE.EQ.0) THEN
              QLAT = ALAT( JDM-1.0-PY,GRIDSZ)*RADIAN
            ELSE
              QLAT = ALATU(JDM-1.0-PY,GRIDSZ)*RADIAN
            ENDIF
            WRITE(6,7150) PLAT,PY,QLAT
 7150       FORMAT(' P(1) Lat =',F8.3,
     +             '  Eq. is at P(',F9.3,
     +             '), P(JDM-1) Lat =',F8.3 /)
          ENDIF
        ELSEIF (MAPTYPE.EQ.1 .OR. MAPTYPE.EQ.3) THEN
C
C         INPUT GRID DISTANCE.
C
          WRITE(6,'(a)',ADVANCE='NO')
     &       ' Eq. is at P(?), 1.0 to exit : '
          READ(5,*) PY
          IF     (PY.EQ.1.0) THEN
            GOTO 111
          ENDIF
C
          IF     (MAPTYPE.EQ.1) THEN
            PLAT = ALAT( 1.0-PY,GRIDSZ)*RADIAN
          ELSE
            PLAT = ALATU(1.0-PY,GRIDSZ)*RADIAN
          ENDIF
          IF     (JDM.LE.1) THEN
            WRITE(6,8100) PY,PLAT
 8100       FORMAT(' Eq. is at P(',F9.3,
     +             '), P(1) Lat =',F8.3 /)
          ELSE
            IF     (MAPTYPE.EQ.1) THEN
              QLAT = ALAT( JDM-1.0-PY,GRIDSZ)*RADIAN
            ELSE
              QLAT = ALATU(JDM-1.0-PY,GRIDSZ)*RADIAN
            ENDIF
            WRITE(6,8150) PY,PLAT,QLAT
 8150       FORMAT(' Eq. is at P(',F9.3,
     +             '), P(1) Lat =',F8.3,
     +             ', P(JDM-1) Lat =',F8.3 /)
          ENDIF
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

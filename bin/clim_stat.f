      PROGRAM CLMSTAT
      IMPLICIT NONE
C
      CHARACTER*40     CTITLE
      INTEGER          IWI,JWI,KWI
      REAL             XFIN,YFIN,DXIN,DYIN,ZLEV(99)
C
      INTEGER          IOS,K
      CHARACTER*240     CFILE
C
C**********
C*
C 1)  PRINT NATIVE CLIMATOLOGY FILE STATISTICS.
C
C 2)  CLIM FILE ON UNIT 55, OR USE THE ENVIRONEMENT VARIABLE FOR055.
C
C 3)  ALAN J. WALLCRAFT,  FEBRUARY 1993.
C*
C**********
C
C     OPEN THE FILE.
C
      CFILE = ' '
      CALL GETENV('FOR055',CFILE)
      IF     (CFILE.EQ.' ') THEN
        CFILE = 'fort.55'
      ENDIF
      OPEN(UNIT=55, FILE=CFILE, FORM='UNFORMATTED', STATUS='OLD',
     +     IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        WRITE(0,*) 'clim_stat: cannot open ',TRIM(CFILE)
        CALL EXIT(1)
        STOP
      ENDIF
C
C     READ THE FIRST TWO RECORDS.
C
      READ(UNIT=55,IOSTAT=IOS) CTITLE
      IF     (IOS.NE.0) THEN
        WRITE(0,*) 'clim_stat: cannot read ',TRIM(CFILE)
        CALL EXIT(2)
        STOP
      ENDIF
      READ( 55)     IWI,JWI,KWI,XFIN,YFIN,DXIN,DYIN,
     +              (ZLEV(K),K=1,KWI)
      CLOSE(UNIT=55)
C
C     STATISTICS.
C
      WRITE(6,6000) CTITLE
      WRITE(6,6100) IWI,JWI,KWI,XFIN,YFIN,DXIN,DYIN,
     +              (ZLEV(K),K=1,KWI)
*     CALL FLUSH(6)
      CALL EXIT(0)
      STOP
C
 6000 FORMAT(A40)
 6100 FORMAT(
     +      'IWI,JWI,KWI =',I4,',',I4,',',I2,
     +   3X,'XFIN,YFIN =',F8.3,',',F7.3,
     +   3X,'DXIN,DYIN =',F5.2,',',F4.2 /
     +      'ZLEV = ',10F7.1 / (7X,10F7.1) )
C     END OF CLMSTAT.
      END

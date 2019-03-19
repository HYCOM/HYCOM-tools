      PROGRAM HYCOM_RIVER_ANOM
      IMPLICIT NONE
C
C  hycom_river_anom - Usage: hycom_river_anom climo.txt daily.txt anom.txt
C
C                 calculates daily river anomaly from climatology
C
C   climo.txt contains 12 monthy river values, 1 per line
C   daily.txt contains daily river values, each line is YYYY-MM-DD value
C   anom.txt appends two values to daily.txt, daily climo and anomaly
C
C   A negative daily value (e.g. -999.00) indicates a missing observation,
C    which is replaced by persisting the last good anomaly w.r.t. 
C    climatology or zero anomaly if there is no previous good value.
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  April 2016.
C
      INTEGER      IARGC
      INTEGER      NARG
      CHARACTER*240 CARG
C
      CHARACTER*240 CFILEC,CFILED,CFILEA
      CHARACTER*240 CLINE
C
      REAL          RIVMON(-2:15)
      REAL          RIVA,RIVC,RIVD,RIVP,W0,W1,W2,W3,X,X1
      INTEGER       IDY,IYR,MON, JDY,JYR,JMN, IOS,K
C
      INTEGER MONDAY(0:13)
      DATA    MONDAY / -31,   0,  31,  59,  90, 120, 151,
     &                      181, 212, 243, 273, 304, 334, 365 /
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.3) THEN
        CALL GETARG(1,CFILEC)
        CALL GETARG(2,CFILED)
        CALL GETARG(3,CFILEA)
      ELSE
        WRITE(6,"(a,a)")
     +    'Usage: hycom_river_anom climo.txt daily.txt anom.txt'
        CALL EXIT(1)
      ENDIF
C
C     OPEN ALL FILES.
C
      OPEN(UNIT=11, FILE=CFILEC, FORM='FORMATTED', STATUS='OLD',
     +     IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error: can''t open ',TRIM(CFILEC)
        WRITE(6,*) 'ios   = ',ios
        CALL EXIT(3)
      ENDIF
      OPEN(UNIT=12, FILE=CFILED, FORM='FORMATTED', STATUS='OLD',
     +     IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error: can''t open ',TRIM(CFILED)
        WRITE(6,*) 'ios   = ',ios
        CALL EXIT(4)
      ENDIF
      OPEN(UNIT=21, FILE=CFILEA, FORM='FORMATTED', STATUS='NEW',
     +     IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error: can''t open ',TRIM(CFILEA)
        WRITE(6,*) 'ios   = ',ios
        CALL EXIT(5)
      ENDIF
C
C --- READ MONTHLY VALUES
C
      READ(11,*) RIVMON(1:12)
      CLOSE(11)
C
      RIVMON(-2) = RIVMON(10)
      RIVMON(-1) = RIVMON(11)
      RIVMON( 0) = RIVMON(12)
      RIVMON(13) = RIVMON( 1)
      RIVMON(14) = RIVMON( 2)
      RIVMON(15) = RIVMON( 3)
C
C     READ ALL DAILY RECORDS
C
      RIVP = 0.0
      JYR  = -1
C
      DO K= 1,HUGE(K)
        READ(12,'(a)',IOSTAT=IOS) CLINE
        IF     (IOS.NE.0) THEN
          EXIT
        ENDIF
C
        IF     (CLINE(1:1).EQ."#") THEN
          WRITE(21,'(a)') TRIM(CLINE)  !a comment
          CYCLE
        ENDIF
C
        READ(CLINE(1: 4),*) IYR
        READ(CLINE(6: 7),*) MON
        READ(CLINE(9:10),*) IDY
        READ(CLINE(11:), *) RIVD
C
        IF     (JYR.NE.-1) THEN
C
C ---     THESE TESTS ARE NOT EXHASTIVE, BUT SHOULD CATCH MOST CASES
C
          IF     (IYR.EQ.JYR .AND. MON.EQ.JMN .AND. IDY.NE.JDY+1) THEN
            WRITE(6,*) 'Error: input is not every day'
            WRITE(6,*) 'old iyr,mon,idy =',jyr,jmn,jdy
            WRITE(6,*) 'new iyr,mon,idy =',iyr,mon,idy
            CALL EXIT(7)
          ELSEIF (IYR.EQ.JYR .AND. MON.NE.JMN .AND. IDY.NE.1) THEN
            WRITE(6,*) 'Error: input is not every day'
            WRITE(6,*) 'old iyr,mon,idy =',jyr,jmn,jdy
            WRITE(6,*) 'new iyr,mon,idy =',iyr,mon,idy
            CALL EXIT(7)
          ELSEIF (IYR.NE.JYR .AND. (MON.NE.1 .OR. IDY.NE.1)) THEN
            WRITE(6,*) 'Error: input is not every day'
            WRITE(6,*) 'old iyr,mon,idy =',jyr,jmn,jdy
            WRITE(6,*) 'new iyr,mon,idy =',iyr,mon,idy
            CALL EXIT(7)
          ENDIF
        ENDIF !not first date
        JYR = IYR
        JMN = MON
        JDY = IDY
C
        IF     (IDY.LE.15) THEN
C
C ---      BETWEEN MON-1 AND MON MONTH CENTERS
C
           X  = ((MONDAY(MON-1) + 15.0) -
     &           (MONDAY(MON)   + IDY )   )/
     &          ((MONDAY(MON-1) + 15.0) -
     &           (MONDAY(MON)   + 15.0)   )
           X1 = 1.0-X
C
           W1 = X1*(1.0+X *(1.0-1.5*X ))
           W2 = X *(1.0+X1*(1.0-1.5*X1))
           W0 = -0.5*X *X1*X1
           W3 = -0.5*X1*X *X
C
           RIVC = W0*RIVMON(MON-2) + 
     &            W1*RIVMON(MON-1) + 
     &            W2*RIVMON(MON  ) + 
     &            W3*RIVMON(MON+1)
        ELSE
C
C ---      BETWEEN MON AND MON+1 MONTH CENTERS
C
           X  = ((MONDAY(MON)   + IDY ) -
     &           (MONDAY(MON)   + 15.0)   )/
     &          ((MONDAY(MON+1) + 15.0) -
     &           (MONDAY(MON)   + 15.0)   )
           X1 = 1.0-X
C
           W1 = X1*(1.0+X *(1.0-1.5*X ))
           W2 = X *(1.0+X1*(1.0-1.5*X1))
           W0 = -0.5*X *X1*X1
           W3 = -0.5*X1*X *X
C
           RIVC = W0*RIVMON(MON-1) + 
     &            W1*RIVMON(MON  ) + 
     &            W2*RIVMON(MON+1) + 
     &            W3*RIVMON(MON+2)
        ENDIF
C
        IF     (RIVD.LT.0.0) THEN  !missing observation
          RIVA = RIVP              !persist the last good anomaly
          RIVD = RIVA + RIVC
        ELSE
          RIVA = RIVD - RIVC
          RIVP = RIVA
        ENDIF
C
        WRITE(21,'(a,3f15.2)') CLINE(1:10),RIVD,RIVC,RIVA
      ENDDO !k
      CLOSE(21)
      CLOSE(12)
      END

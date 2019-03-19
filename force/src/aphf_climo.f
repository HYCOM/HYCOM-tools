      PROGRAM APHF_CLIMO
      USE MOD_ZA  ! HYCOM array I/O interface
      IMPLICIT NONE
C
C     FLUX ARRAYS.
C
      INTEGER, ALLOCATABLE :: MSK(:,:)
      REAL*4,  ALLOCATABLE :: FLX(:,:)
C
      CHARACTER PREAMBL(5)*79,CLINE*80
      INTEGER   INDX
C
C     NAMELIST.
C
      REAL*8           TSTART,TMAX
      NAMELIST/AFTIME/ TSTART,TMAX
C
C**********
C*
C 1)  FROM A CLIMATOLOGICAL (SINGLE YEAR REPEATING) FORCING FILE,
C      CREATE AN ACTUAL-YEAR FILE.
C
C 2)  NAMELIST INPUT:
C
C     /AFTIME/
C        TMAX   - TIME OF END   OF CURRENT INTEGRATION     (DAYS)
C        TSTART - TIME OF START OF CURRENT INTEGRATION     (DAYS)
C
C 3)  INPUT:
C        ON UNIT  5:    NAMELIST /AFTIME/
C        ON UNIT 20:    UNFORMATTED MODEL CLIMO. FLUX FILE, SEE (4).
C     OUTPUT:
C        ON UNIT 10:    UNFORMATTED MODEL ACTUAL FLUX FILE, SEE (4).
C
C 4)  THE FLUXES ARE AT EVERY GRID POINT OF THE MODEL'S 'P' GRID.
C     ARRAY SIZE IS 'IDM' BY 'JDM'.
C
C 5)  ALAN J. WALLCRAFT, NRL, SEPTEMBER 2009.
C*
C**********
C
      INTEGER I,IYEAR,IDAY,IHOUR,IOS,KREC,KOUT,KREWIND
      REAL*8  FDAYIN,FDAYOUT,FDAYOFF,FINC
      REAL*4  HMINA,HMINB,HMAXA,HMAXB
C
C --- MODEL ARRAYS.
C
      CALL XCSPMD  !define idm,jdm
      ALLOCATE(  MSK(IDM,JDM),
     +           FLX(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,'(/ a,i9,a /)')
     +    'error - can''t allocate',IDM*JDM*3,' words'
        CALL ZHFLSH(6)
        STOP
      ENDIF
C
      CALL ZAIOST
C
C     NAMELIST INPUT.
C
      CALL ZHOPEN(6, 'FORMATTED', 'UNKNOWN', 0)
C
      TSTART = 0.0
      TMAX   = 0.0
      WRITE(6,*) 'READING /AFTIME/'
      CALL ZHFLSH(6)
      READ( 5,AFTIME)
      WRITE(6,AFTIME)
      WRITE(6,*) 
      CALL ZHFLSH(6)
C
C --- CONVERSION FACTOR FROM CLIMO TO ACTUAL YEAR
C
      CALL FORDAY(TSTART,3, IYEAR,IDAY,IHOUR)
      FDAYOFF = (INT(TSTART) - IDAY) - 1095  
C
C     INITIALIZE ARRAY INPUT AND OUTPUT.
C
      CALL ZAIOPN('NEW', 10)
      CALL ZAIOPN('OLD', 20)
      CALL ZHOPEN(10, 'FORMATTED', 'NEW', 0)
      CALL ZHOPEN(20, 'FORMATTED', 'OLD', 0)
C
      READ( 20,4101) PREAMBL
      WRITE(10,4101) PREAMBL
C
      KREWIND = 0
      KOUT    = 0
      DO KREC= 1,HUGE(KREC)
C
C       READ IN FLUXES.
C
        READ(20,'(A)',IOSTAT=IOS) CLINE
        IF     (IOS.NE.0) THEN
          IF     (KREWIND.LT.10) THEN
            KREWIND = KREWIND + 1
            FDAYOFF = FDAYOFF + 366.D0
            CALL ZAIORW(20)
            REWIND(20)
            READ( 20,4101) PREAMBL
            READ(20,'(A)') CLINE
          ELSE  !probably in an infinate loop
            WRITE(6,'(A)') 
     &        'end of input file'
            EXIT
          ENDIF
        ENDIF
C
        INDX = INDEX(CLINE,'=')
        READ(CLINE(INDX+1:),*) FDAYIN,FINC,HMINB,HMAXB
        FDAYOUT = FDAYIN + FDAYOFF
C
        IF     (FDAYOUT+FINC.LE.TSTART) THEN
          CALL ZAIOSK(20)
          CYCLE
        ENDIF
        IF     (FDAYOUT-FINC.GE.TMAX)   THEN
          WRITE(6,'(A,i5,2f10.3)') 
     &      'rec,dayin,dayout =',krec,fdayin,fdayout
          CALL ZHFLSH(6)
          EXIT
        ENDIF
C
        KOUT = KOUT + 1
        IF     (KOUT.EQ.1) THEN
          WRITE(6,'(A,i5,2f10.3)') 
     &      'rec,dayin,dayout =',krec,fdayin,fdayout
          CALL ZHFLSH(6)
        ENDIF
C
        IF     (HMINB.EQ.HMAXB) THEN  !constant field
          HMINA    = HMINB
          HMAXA    = HMINB
          FLX(:,:) = HMINB
          CALL ZAIOSK(20)
        ELSE
          CALL ZAIORD(FLX,MSK,.FALSE., HMINA,HMAXA, 20)
          IF     (ABS(HMINA-HMINB).GT.ABS(HMINB)*1.E-4 .OR.
     &            ABS(HMAXA-HMAXB).GT.ABS(HMAXB)*1.E-4     ) THEN
            WRITE(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &        'error - .a and .b grid files not consistent (frm):',
     &        '.a,.b min = ',HMINA,HMINB,HMINA-HMINB,
     &        '.a,.b max = ',HMAXA,HMAXB,HMAXA-HMAXB
            CALL ZHFLSH(6)
            STOP
          ENDIF
        ENDIF
C
C       WRITE OUT HYCOM FLUXS.
C
        I = MAX( 1, INDX-25 )  !4112 line at most 80 characters
        CALL ZAIOWR(FLX,MSK,.FALSE., HMINA,HMAXA, 10, .FALSE.)
        WRITE(10,4112) CLINE(I:INDX),FDAYOUT,FINC,HMINA,HMAXA
        CALL ZHFLSH(10)
        WRITE( 6,4112) CLINE(I:INDX),FDAYOUT,FINC,HMINA,HMAXA
        CALL ZHFLSH(6)
      ENDDO !KREC
C
      CALL ZAIOCL(10)
      CLOSE( UNIT=10)
      CALL ZAIOCL(20)
      CLOSE( UNIT=20)
      STOP
C
 4101 FORMAT(A79)
 4112 FORMAT(A,F12.5,F10.6,1P2E16.7)
C     END OF PROGRAM APHF_CLIMO.
      END

      subroutine forday(dtime,yrflag, iyear,iday,ihour)
      implicit none
c
      real*8  dtime
      integer yrflag, iyear,iday,ihour
c
c --- converts model day to "calendar" date (year,julian-day,hour).
c
      real*8  dtim1,day
      integer iyr,nleap
c
      if     (yrflag.eq.0) then
c ---   360 days per model year, starting Jan 16
        iyear =  int((dtime+15.001d0)/360.d0) + 1
        iday  =  mod( dtime+15.001d0 ,360.d0) + 1
        ihour = (mod( dtime+15.001d0 ,360.d0) + 1.d0 - iday)*24.d0
c
      elseif (yrflag.eq.1) then
c ---   366 days per model year, starting Jan 16
        iyear =  int((dtime+15.001d0)/366.d0) + 1
        iday  =  mod( dtime+15.001d0 ,366.d0) + 1
        ihour = (mod( dtime+15.001d0 ,366.d0) + 1.d0 - iday)*24.d0
c
      elseif (yrflag.eq.2) then
c ---   366 days per model year, starting Jan 01
        iyear =  int((dtime+ 0.001d0)/366.d0) + 1
        iday  =  mod( dtime+ 0.001d0 ,366.d0) + 1
        ihour = (mod( dtime+ 0.001d0 ,366.d0) + 1.d0 - iday)*24.d0
c
      elseif (yrflag.eq.3) then
c ---   model day is calendar days since 01/01/1901
        iyr   = (dtime-1.d0)/365.25d0
        nleap = iyr/4
        dtim1 = 365.d0*iyr + nleap + 1.d0
        day   = dtime - dtim1 + 1.d0
        if     (dtim1.gt.dtime) then
          iyr = iyr - 1
        elseif (day.ge.367.d0) then
          iyr = iyr + 1
        elseif (day.ge.366.d0 .and. mod(iyr,4).ne.3) then
          iyr = iyr + 1
        endif
        nleap = iyr/4
        dtim1 = 365.d0*iyr + nleap + 1.d0
c
        iyear =  1901 + iyr
        iday  =  dtime - dtim1 + 1
        ihour = (dtime - dtim1 + 1.d0 - iday)*24.d0
c
      else
        write( 6,*)
        write( 6,*) 'error in forday - unsupported yrflag value'
        write( 6,*)
        stop '(forday)'
      endif
      return
      end

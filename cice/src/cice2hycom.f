      PROGRAM CICE2HYCOM
      USE MOD_ZA  ! HYCOM array I/O interface
      IMPLICIT NONE
C
C     CICE ARRAY.
C
      REAL*4,  ALLOCATABLE :: CICE(:,:),WORK(:,:)
      INTEGER, ALLOCATABLE ::  MSK(:,:)
C
C**********
C*
C 1)  CONVERT A CICE NETCDF FILE TO A HYCOM arche-like .a/.b FILE.
C     INTERPOLATE VELOCITY FIELDS TO THE CELL CENTER FIRST.
C
C 2)  INPUT IS:
C       cice       filename  (existing file, .nc)
C       bathymetry filename  (existing file, .a or .b)
C       hycom      filename  (new file, .a or .b)
C       4-line header
C       experiment number x10
C       field 1 name
C       field 2 name
C          ...
C       field n name
C
C 3)  ALAN J. WALLCRAFT,  NRL, MARCH 2011.
C*
C**********
C
      INTEGER       I,IEXPT,IOS,J
      REAL          WDAY,WYR,FDY, XMIN,XMAX
      REAL*8        WDAY8
      CHARACTER*8   CNAMEH
      CHARACTER*10  CNAME
      CHARACTER*80  CLINE
      CHARACTER*256 CFILE,HFILE,AFILE,BFILE
C
C     GET ARRAY SIZE.
C
      READ(*,'(a)') CFILE
      WRITE(6,6000) ' INPUT:',TRIM(CFILE)
      CALL CREADI(CFILE, WDAY,IDM,JDM)
      WDAY8 = WDAY
      CALL WNDAY(WDAY8, WYR,FDY)
      WRITE(6,6100) '  DATE:',FDY,NINT(WYR)
C
      CALL ZAIOST
C
      ALLOCATE( CICE(IDM,JDM),
     &          WORK(IDM,JDM),
     &           MSK(IDM,JDM) )
C
C     INITIALIZE MASK FROM BATHYMETRY.
C
      READ(*,'(a)') HFILE
      WRITE(6,6000) ' BATHY:',TRIM(HFILE)
      CALL ZHFLSH(6)
C
      AFILE=HFILE(1:LEN_TRIM(HFILE)-2) // ".a"
      CALL ZAIOPF(AFILE,'OLD', 20)
      CALL ZAIORD(CICE,MSK,.FALSE., XMIN,XMAX, 20)
      CALL ZAIOCL(20)
C
      DO J= 1,JDM
        DO I= 1,IDM
          IF     (CICE(I,J).GT.0.5*2.0**100 .OR.
     &            CICE(I,J).LE.0.0              ) THEN
            MSK(I,J) = 0
          ELSE
            MSK(I,J) = 1
          ENDIF
        ENDDO
      ENDDO
C
C     INITIALIZE OUTPUT.
C
      READ(*,'(a)') HFILE
      WRITE(6,6000) 'OUTPUT:',TRIM(HFILE)
      CALL ZHFLSH(6)
C
      AFILE=HFILE(1:LEN_TRIM(HFILE)-2) // ".a"
      BFILE=HFILE(1:LEN_TRIM(HFILE)-2) // ".b"
      CALL ZAIOPF(AFILE,'NEW', 10)
      CALL ZHOPNC(10, BFILE, 'FORMATTED', 'NEW', 0)
C
      DO I= 1,4
        READ(  *,'(a)')   CLINE
        WRITE( 6,'(a80)') CLINE
        WRITE(10,'(a80)') CLINE
        CALL ZHFLSH(10)
      ENDDO
C
      READ(*,*) IEXPT
      WRITE( 6,4116) 22,IEXPT,3,IDM,JDM
      WRITE(10,4116) 22,IEXPT,3,IDM,JDM
      CALL ZHFLSH(10)
C
C     LOOP THROUGH FIELDS
C
      DO
        READ(*,'(a)',IOSTAT=IOS) CNAME
        IF     (IOS.NE.0) THEN
          EXIT
        ENDIF
C
        CALL CREAD1(CFILE,CNAME, CICE,WORK,IDM,JDM)
C
        J = LEN_TRIM(CNAME)
        IF     (J.LE.8) THEN
          CNAMEH = CNAME
        ELSE
          I = INDEX(CNAME,'_')
          IF     (I.EQ.0) THEN
            CNAMEH = CNAME(1:8)
          ELSE
            CNAMEH = CNAME(1:8-(J-I+1)) // CNAME(I:J)
          ENDIF
        ENDIF
        CALL ZAIOWR(CICE,MSK,.TRUE., XMIN,XMAX, 10, .FALSE.)
        WRITE(10,4117) CNAMEH,0,WDAY,0,0.0,XMIN,XMAX
        CALL ZHFLSH(10)
        WRITE( 6,4117) CNAMEH,0,WDAY,0,0.0,XMIN,XMAX
        CALL ZHFLSH(6)
      ENDDO
C
      CALL ZAIOCL(10)
      CLOSE( UNIT=10)
      STOP
C
 4116 format (
     & i5,4x,'''iversn'' = hycom version number x10'/
     & i5,4x,'''iexpt '' = experiment number x10'/
     & i5,4x,'''yrflag'' = days in year flag'/
     & i5,4x,'''idm   '' = longitudinal array size'/
     & i5,4x,'''jdm   '' = latitudinal  array size'/
     & 'field       time step  model day',
     & '  k  dens        min              max')
 4117 format (a8,' =',i11,f11.3,i3,f7.3,1p2e16.7)
 6000 FORMAT(1X,A,2X,A)
 6100 FORMAT(1X,A,2X,F7.2,'/',I4)
C     END OF PROGRAM CICE2HYCOM.
      END
      SUBROUTINE WNDAY(WDAY, YEAR,DAY)
      IMPLICIT NONE
      REAL*8 WDAY
      REAL*4 YEAR,DAY
C
C**********
C*
C  1) CONVERT 'FLUX DAY' INTO JULIAN DAY AND YEAR.
C
C  2) THE 'FLUX DAY' IS THE NUMBER OF DAYS SINCE 001/1901 (WHICH IS 
C      FLUX DAY 1.0).
C     FOR EXAMPLE:
C      A) YEAR=1901.0 AND DAY=1.0, REPRESENTS 0000Z HRS ON 001/1901
C         SO WDAY WOULD BE 1.0.
C      B) YEAR=1901.0 AND DAY=2.5, REPRESENTS 1200Z HRS ON 002/1901
C         SO WDAY WOULD BE 2.5.
C     YEAR MUST BE NO LESS THAN 1901.0, AND NO GREATER THAN 2099.0.
C     NOTE THAT YEAR 2000 IS A LEAP YEAR (BUT 1900 AND 2100 ARE NOT).
C
C  3) ALAN J. WALLCRAFT, PLANNING SYSTEMS INC., FEBRUARY 1993.
C*
C**********
C
      INTEGER IYR,NLEAP
      REAL*8  WDAY1
C
C     FIND THE RIGHT YEAR.
C
      IYR   = (WDAY-1.0)/365.25
      NLEAP = IYR/4
      WDAY1 = 365.0*IYR + NLEAP + 1.0
      DAY   = WDAY - WDAY1 + 1.0
      IF     (WDAY1.GT.WDAY) THEN
        IYR   = IYR - 1
      ELSEIF (DAY.GE.367.0) THEN
        IYR   = IYR + 1
      ELSEIF (DAY.GE.366.0 .AND. MOD(IYR,4).NE.3) THEN
        IYR   = IYR + 1
      ENDIF
      NLEAP = IYR/4
      WDAY1 = 365.0*IYR + NLEAP + 1.0
C
C     RETURN YEAR AND JULIAN DAY.
C
      YEAR = 1901 + IYR
      DAY  = WDAY - WDAY1 + 1.0
      RETURN
C     END OF WNDAY.
      END
      SUBROUTINE CREADI(CFILE, WDAY,IDM,JDM)
      USE netcdf   ! NetCDF fortran 90 interface
      IMPLICIT NONE
C
      CHARACTER*256 CFILE
      REAL          WDAY
      INTEGER       IDM,JDM
C
C**********
C*
C  1)  INITIALIZE ARRAY SIZES FOR READING CICE FIELDS.
C
C      SEE 'FREAD1' FOR READING ACTUAL CICE RECORDS.
C
C  2) ON EXIT:
C      WDAY = SAMPLE DAY (SINCE 1900-12-31 00:00:00)
C      IDM  = LONGIUDINAL ARRAY DIMENSION
C      IDM  = LATITUDINAL ARRAY DIMENSION
C
C  3) ALAN J. WALLCRAFT,  MARCH 2011.
C*
C*********
C
      REAL    TIME(1)
      INTEGER ncFID,ncDID,ncVID
C
      call nchek('nf90_open',
     &            nf90_open(trim(CFILE), nf90_nowrite, ncFID))
C
      call nchek('nf90_inq_dimid-nj',
     &            nf90_inq_dimid(        ncFID, 'nj',ncDID))
      call nchek('nf90_inquire_dimension-nj',
     &            nf90_inquire_dimension(ncFID,      ncDID,
     &                                               len=JDM))
C
      call nchek('nf90_inq_dimid-ni',
     &            nf90_inq_dimid(        ncFID, 'ni',ncDID))
      call nchek('nf90_inquire_dimension-ni',
     &            nf90_inquire_dimension(ncFID,      ncDID,
     &                                               len=IDM))
C
      call nchek('nf90_inq_varid-time',
     &            nf90_inq_varid(ncFID,'time',ncVID))
      call nchek('nf90_get_var-time',
     &            nf90_get_var(  ncFID,        ncVID,TIME(1:1)))
      WDAY = TIME(1)
C
      call nchek("nf90_close",
     &            nf90_close(ncFID))
      RETURN
C     END OF FREADI.
      END
      SUBROUTINE CREAD1(CFILE,CNAME, CICE,WORK,IDM,JDM)
      USE netcdf   ! NetCDF fortran 90 interface
      IMPLICIT NONE
C
      CHARACTER*256 CFILE
      CHARACTER*10  CNAME
      INTEGER       IDM,JDM
      REAL*4        CICE(IDM,JDM),WORK(IDM,JDM)
C
C**********
C*
C  1) READ A SINGLE CICE FIELD.
C
C  2) INTERPOLATE TO THE CELL CENTER IN NECESSARY.
C
C  3) ALAN J. WALLCRAFT,  NRL,  MARCH 2011.
C*
C*********
C
      REAL*4     SPVALH
      PARAMETER (SPVALH=2.0**100)  !HYCOM DATA VOID
C
      CHARACTER*14 COORDS
      INTEGER      ncFID,ncDID,ncVID
      INTEGER      I,IP,J,JP
C
C     READ THE FLUX RECORD.
C
      call nchek('nf90_open',
     &            nf90_open(trim(CFILE), nf90_nowrite, ncFID))
C
      write(6,*) "CNAME  = ",TRIM(CNAME)
      call nchek('nf90_inq_varid-CNAME',
     &            nf90_inq_varid(ncFID,CNAME,ncVID))
      call nchek("nf90_get_att-coordinates",
     &            nf90_get_att(  ncFID,      ncVID,
     &                                       "coordinates",
     &                                       COORDS))
      write(6,*) "COORDS = ",trim(COORDS)
      IF     (COORDS(1:4).EQ.'TLON') THEN
        call nchek('nf90_get_var-CICE',
     &              nf90_get_var(  ncFID,      ncVID,
     &                                         CICE(:,:),
     &                                         (/ 1,1,1 /) ))
        DO J= 1,JDM
          DO I= 1,IDM
            IF     (CICE(I,J).GT.0.9E+30) THEN
              CICE(I,J) = SPVALH
            ENDIF
          ENDDO !i
        ENDDO !j
      ELSEIF (COORDS(1:4).EQ.'ULON') THEN
        call nchek('nf90_get_var-CICE',
     &              nf90_get_var(  ncFID,      ncVID,
     &                                         WORK(:,:),
     &                                         (/ 1,1,1 /) ))
C
C       SET VELOCITY TO ZERO OVER LAND AND WHERE ICE FREE
C
        DO J= 1,JDM
          DO I= 1,IDM
            IF     (WORK(I,J).GT.0.9E+30) THEN
              WORK(I,J) = 0.0
            ENDIF
          ENDDO !i
        ENDDO !j
        DO J= 1,JDM
          JP = MIN(J+1,JDM)
          DO I= 1,IDM
            IP = MOD(I,IDM) + 1
            CICE(I,J) = 0.25*(WORK(I, J ) +
     &                        WORK(IP,J ) +
     &                        WORK(I, JP) +
     &                        WORK(IP,JP)  )
          ENDDO !i
        ENDDO !j
      ELSE
        write(6,*) 'Unkown coordinates'
        write(6,*) trim(COORDS)
        STOP
      ENDIF
      call nchek("nf90_close",
     &            nf90_close(ncFID))
      RETURN
C     END OF FREAD
      END
      subroutine nchek(cnf90,status)
      use netcdf   ! NetCDF fortran 90 interface
      implicit none
c
      character*(*), intent(in) :: cnf90
      integer,       intent(in) :: status
c
c     subroutine to handle NetCDF errors
c
*     if     (.TRUE. ) then !debug
      if     (.FALSE.) then !nodebug
        write(6,'(a)') trim(cnf90)
      endif

      if (status /= nf90_noerr) then
        write(6,'(/a)')   'error from NetCDF library'
        write(6,'(a/)')   trim(cnf90)
        write(6,'(a/)')   trim(nf90_strerror(status))
        stop
      end if
      end subroutine nchek

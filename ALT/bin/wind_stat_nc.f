      PROGRAM WIND_STAT_NC
      USE netcdf   ! NetCDF fortran 90 interface
      IMPLICIT NONE
C
      CHARACTER*40     CTITLE
      INTEGER          IWI,JWI,NREC
      REAL             XFIN,YFIN,DXIN,DYIN
C
      CHARACTER*240    CFILE
      INTEGER          ncFID,ncDID,ncVID
      INTEGER          KREC
      REAL             JDAY,YEAR
      DOUBLE PRECISION TIME_NEXT
C
      REAL,             ALLOCATABLE :: WDAY(:)
      DOUBLE PRECISION, ALLOCATABLE :: TIME(:),WLON(:),WLAT(:)
C
C**********
C*
C 1)  PRINT MODEL WIND FILE STATISTICS.
C
C 2)  NRL-style NetCDF WIND FILE FROM ENVIRONEMENT VARIABLE CDF055.
C
C 3)  ALAN J. WALLCRAFT,  JULY 2010.
C*
C**********
c
      integer          i
      double precision time_i
C
C     OPEN THE FILE.
C
      CFILE = ' '
      CALL GETENV('CDF055',CFILE)
      IF     (CFILE.EQ.' ') THEN
        WRITE(0,*) 'wind_stat_nc: no CDF055 environment variable'
        CALL EXIT(1)
        STOP
      ENDIF
      call nchek('nf90_open',
     &            nf90_open(trim(CFILE), nf90_nowrite, ncFID))
C
C     READ THE TITLE.
C
      call nchek("nf90_get_att-title",
     &            nf90_get_att(ncFID,nf90_global,
     &                               "title",
     &                               CTITLE))
C
C     READ THE AXES.
C
      call nchek('nf90_inq_dimid-MT',
     &            nf90_inq_dimid(        ncFID, 'MT',ncDID))
      call nchek('nf90_inquire_dimension-MT',
     &            nf90_inquire_dimension(ncFID,      ncDID,
     &                                               len=NREC))
      ALLOCATE( TIME(NREC) )
      call nchek('nf90_inq_varid-MT',
     &            nf90_inq_varid(ncFID,'MT',ncVID))
      call nchek('nf90_get_var-MT',
     &            nf90_get_var(  ncFID,     ncVID,TIME(:)))
      call nchek("nf90_get_att-next_MT",
     &            nf90_get_att(ncFID,       ncVID,
     &                                      "next_MT",
     &                                      TIME_NEXT))
C
      do i= 2,nrec
        time_i = time(1) + 
     &             (time_next - time(1))/nrec * (i-1.d0)
        time_i = NINT(time_i*24.d0)/24.d0
        if     (abs(time(i)-time_i).gt.1.d-5) then
          write(6,'(i5,f20.6,f20.16)') i,time(i),time(i)-time_i
        endif
      enddo !i
C
      call nchek('nf90_inq_dimid-Latitude',
     &            nf90_inq_dimid(        ncFID, 'Latitude',ncDID))
      call nchek('nf90_inquire_dimension-Latitude',
     &            nf90_inquire_dimension(ncFID,            ncDID,
     &                                                 len=JWI))
      ALLOCATE( WLAT(JWI) )
      call nchek('nf90_inq_varid-Latitude',
     &            nf90_inq_varid(ncFID,'Latitude',ncVID))
      call nchek('nf90_get_var-Latitude',
     &            nf90_get_var(  ncFID,           ncVID,WLAT(:)))
C
      call nchek('nf90_inq_dimid-Longitude',
     &            nf90_inq_dimid(        ncFID, 'Longitude',ncDID))
      call nchek('nf90_inquire_dimension-Longitude',
     &            nf90_inquire_dimension(ncFID,             ncDID,
     &                                                  len=IWI))
      ALLOCATE( WLON(IWI) )
      call nchek('nf90_inq_varid-Longitude',
     &            nf90_inq_varid(ncFID,'Longitude',ncVID))
      call nchek('nf90_get_var-Longitude',
     &            nf90_get_var(  ncFID,            ncVID,WLON(:)))
      call nchek("nf90_close",
     &            nf90_close(ncFID))
C
      XFIN = WLON(1)
      DXIN = (WLON(IWI) - WLON(1))/(IWI-1.D0)  !assume uniform spaced
      YFIN = WLAT(1)
      IF     (ABS((WLAT(JWI)-WLAT(1))-
     &            (WLAT(  2)-WLAT(1))*(JWI-1.D0)).LT.1.D-2) THEN
*       write(6,*)     (WLAT(JWI)-WLAT(1))
*       write(6,*)     (WLAT(  2)-WLAT(1))*(JWI-1.D0)
*       write(6,*) ABS((WLAT(JWI)-WLAT(1))-
*    &                 (WLAT(  2)-WLAT(1))*(JWI-1.D0))
        DYIN = (WLAT(JWI) - WLAT(1))/(JWI-1.D0)  !assume uniform spaced
      ELSE
        DYIN = 0.0  ! assume a gaussian grid
      ENDIF
      ALLOCATE( WDAY(NREC+1) )
      DO KREC= 1,NREC
        WDAY(KREC) = TIME(KREC)
      ENDDO
      WDAY(NREC+1) = TIME_NEXT
C
C     STATISTICS.
C
      WRITE(6,6000) CTITLE
      WRITE(6,6100) IWI,JWI,XFIN,YFIN,DXIN,DYIN,
     +              NREC,(WDAY(KREC), KREC=1,NREC+1)
*     CALL FLUSH(6)
C
C     SUMMARY.
C
      CALL WNDAY(WDAY(1), YEAR,JDAY)
      IF     (YEAR.LT.1904.5) THEN
        WRITE(6,6200) NREC,JDAY,NINT(YEAR),WDAY(NREC+1)-WDAY(1)
      ELSE
        WRITE(6,6250) NREC,JDAY,NINT(YEAR),WDAY(NREC+1)-WDAY(1)
      ENDIF
      CALL EXIT(0)
      STOP
C
 6000 FORMAT(A40)
 6100 FORMAT(
     +      'IWI,JWI =',I4,',',I4,
     +   3X,'XFIN,YFIN =',F9.3,',',F9.3,
     +   3X,'DXIN,DYIN =',F6.3,',',F6.3 /
     +      'NREC =',I5,5X,'WDAY =' / (8F10.3) )
 6200 FORMAT(I5,' RECORD CLIMATOLOGY STARTING ON',F7.2,'/',I4,
     +   ' COVERING',F9.2,' DAYS')
 6250 FORMAT(I5,' WIND RECORDS STARTING ON',F7.2,'/',I4,
     +   ' COVERING',F9.2,' DAYS')
C     END OF WDSTAT.
      END
      SUBROUTINE WNDAY(WDAY, YEAR,DAY)
      IMPLICIT NONE
      REAL WDAY,YEAR,DAY
C
C**********
C*
C  1) CONVERT 'WIND DAY' INTO JULIAN DAY AND YEAR.
C
C  2) THE 'WIND DAY' IS THE NUMBER OF DAYS SINCE 001/1901 (WHICH IS 
C      WIND DAY 1.0).
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
      REAL    WDAY1
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
        write(6,'(/a)')   'error in profout - from NetCDF library'
        write(6,'(a/)')   trim(cnf90)
        write(6,'(a/)')   trim(nf90_strerror(status))
        stop
      end if
      end subroutine nchek

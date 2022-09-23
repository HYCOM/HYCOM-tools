      PROGRAM STOCH
      USE MOD_ZA  ! HYCOM array I/O interface
      IMPLICIT NONE
C
C     STOCHASTIC ARRAYS.
C
      INTEGER, ALLOCATABLE :: MSK(:,:)
      REAL*4,  ALLOCATABLE :: PLON(:,:),PLAT(:,:)
      REAL*4,  ALLOCATABLE :: VPM(:,:,:),TPM(:,:,:)
C
      LOGICAL   LARCTIC
      CHARACTER PREAMBL(5)*79
C
C     NAMELIST.
C
      CHARACTER*40     CTITLE
      NAMELIST/AFTITL/ CTITLE
C
      REAL*8           FSTART,DTI,SINC,SSTART,WSTART,TSTART,TMAX
      REAL*4           STFR_SPC(3),STFR_TIM(3),STFR_RMS(3)
      NAMELIST/AFTIME/ FSTART,DTI,SINC,SSTART,WSTART,TSTART,TMAX,
     +                 STFR_SPC,   STFR_TIM,   STFR_RMS
C
C**********
C*
C 1)  CREATE T,S,U,V MODEL GRID STOCHASTIC SURFACE FORCING
C      FILES SUITABLE FOR INPUT TO HYCOM OVER THE GIVEN REGION.
C
C 2)  NAMELIST INPUT:
C
C     /AFTITL/
C        CTITLE - ONE (40-CHARACTER) LINE TITLE.
C
C     /AFTIME/
C        FSTART - TIME OF FLUX START, IGNORED HERE         (DAYS)
C        WSTART - TIME OF WIND START, IGNORED HERE         (DAYS)
C        DTI    - NOMINAL TIME STEP FOR STOCHASTIC FORCING (SECS)
C        SINC   - TIME INCREMENT BETWEEN OUTPUT FIELDS     (DAYS)
C        SSTART - TIME OF STOCHASTIC FORCING               (DAYS)
C        TMAX   - TIME OF END   OF CURRENT INTEGRATION     (DAYS)
C        TSTART - TIME OF START OF CURRENT INTEGRATION     (DAYS)
C        
C        STFR_SPC - HOR SPATIAL SCALE FOR STOCHASTIC FORCING (M)
C                    (1:T;2:S;3:UV)
C        STFR_TIM - TEMPORAL    SCALE FOR STOCHASTIC FORCING (HR)
C                    (1:T;2:S;3:UV)
C        STFR_RMS - RMS AMPLITUDE     FOR STOCHASTIC FORCING (UNITS/HR)
C                    (1:T;2:S;3:UV)
C
C     NAMELIST /AFTIME/ IS PATTERNED AFTER /XXTIME/ SO THAT THE
C      MODEL''S STANDARD AWK-BASED RUN SCRIPT CUSTOMIZER CAN ALSO
C      WORK FOR THE FLUX GENERATION SCRIPT.  IN PARTICULAR, 'WSTART'
C      AND 'FSTART' ARE READ IN, BUT NOT USED.
C
C 4)  INPUT:
C        ON UNIT  5:    NAMELIST /AFTITL/, /AFTIME/
C     OUTPUT:
C        ON UNIT 10:    UNFORMATTED MODEL T STOCASTIC FILE, SEE (6).
C        ON UNIT 11:    UNFORMATTED MODEL S STOCASTIC FILE, SEE (6).
C        ON UNIT 12:    UNFORMATTED MODEL U STOCASTIC FILE, SEE (6).
C        ON UNIT 13:    UNFORMATTED MODEL V STOCASTIC FILE, SEE (6).
C
C 6)  THE OUTPUT FIELDS ARE AT EVERY GRID POINT OF THE MODEL'S
C     'P' GRID.  ARRAY SIZE IS 'IDM' BY 'JDM'.
C
C 7)  SEVERAL STATISTICS ARE WRITTEN OUT IN ORDER TO CHECK THE 
C      INTERPOLATION BETWEEN VARIOUS MACHINES.  MIN, MAX, MEAN AND RMS 
C      OF THE ENTIRE BASIN ARE OUTPUT FOR kPAR FOR EACH 
C      RECORD.  NOTE HOWEVER THAT THESE VALUES MAY NOT REPRESENT THE
C      STATISTICS OF THE FLUXS AS SEEN BY THE MODEL, IF THE INPUT FLUX 
C      DATA HAS NON-REALISTIC VALUES OVER LAND.  IT IS UP TO THE USER 
C      TO CHECK THE LOG FILES FOR CONSISTENCY BETWEEN MACHINES.
C
C 8)  ALAN J. WALLCRAFT,  NRL, SEPTEMBER 2017.
C      BASED ON NCOM VERSIONS BY PAUL MARTIN AND CLARK ROWLEY.
C*
C**********
C
      EXTERNAL AVERMS,MINMAX
C
      REAL*4     ZERO,RADIAN
      PARAMETER (ZERO=0.0, RADIAN=57.2957795)
C
      INTEGER NREC,MREC,NSTEP
C
      CHARACTER*80 CLINE
      REAL*8  WDAY8
      REAL*4  HMINA,HMINB,HMAXA,HMAXB
C
      INTEGER I,INIT,J,KREC,KSTEP
      REAL*4  FDY,WYR,DXREF,DTI4,
     +        XMIN,XMAX,XAVE,XRMS
C
C --- MODEL ARRAYS.
C
      CALL XCSPMD  !define idm,jdm
      ALLOCATE(  MSK(IDM,JDM),
     +          PLON(IDM,JDM),
     +          PLAT(IDM,JDM),
     +           TPM(IDM,JDM,2),
     +           VPM(IDM,JDM,2) )
C
C     NAMELIST INPUT.
C
      CALL ZHOPEN(6, 'FORMATTED', 'UNKNOWN', 0)
C
      CTITLE = ' '
      WRITE(6,*) 'READING /AFTITL/'
      CALL ZHFLSH(6)
      READ( 5,AFTITL)
      WRITE(6,AFTITL)
C
      FSTART = 0.0
      DTI    = 0.0
      SINC   = 0.0
      SSTART = 0.0
      WSTART = 0.0
      TSTART = 0.0
      TMAX   = 0.0
      WRITE(6,*) 'READING /AFTIME/'
      CALL ZHFLSH(6)
      READ( 5,AFTIME)
      IF     (SINC.GT.0.03D0) THEN
        SINC = NINT(SINC*24.D0)/24.D0  !assume multiple of 1 hour
      ENDIF
      WRITE(6,AFTIME)
      WRITE(6,*) 
      CALL ZHFLSH(6)
C
      IF     (SSTART.LT.TSTART) THEN
        WRITE(6,*) 'ERROR - SSTART MUST BE AT LEAST TSTART'
        WRITE(6,*)
        CALL ZHFLSH(6)
        STOP
      ENDIF
      IF     (SINC.LE.0.0) THEN
        WRITE(6,*) 'ERROR - SINC MUST BE POSITIVE'
        WRITE(6,*)
        CALL ZHFLSH(6)
        STOP
      ENDIF
      IF     (DTI.LE.0.0) THEN
        WRITE(6,*) 'ERROR - DTI MUST BE POSITIVE'
        WRITE(6,*)
        CALL ZHFLSH(6)
        STOP
      ENDIF
C
C     GRID INPUT.
C
      CALL ZAIOST
C
      CALL ZHOPNC(21, 'regional.grid.b', 'FORMATTED', 'OLD', 0)
      CALL ZAIOPF(    'regional.grid.a', 'OLD', 21)
C
      READ(21,*) ! skip idm
      READ(21,*) ! skip jdm
      READ(21,*) ! skip mapflg
      READ(21,'(A)') CLINE
      I = INDEX(CLINE,'=')
      READ (CLINE(I+1:),*)   HMINB,HMAXB
      CALL ZAIORD(PLON,MSK,.FALSE., HMINA,HMAXA, 21)
      IF     (ABS(HMINA-HMINB).GT.ABS(HMINB)*1.E-4 .OR.
     &        ABS(HMAXA-HMAXB).GT.ABS(HMAXB)*1.E-4     ) THEN
        WRITE(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - .a and .b grid files not consistent (plon):',
     &    '.a,.b min = ',HMINA,HMINB,HMINA-HMINB,
     &    '.a,.b max = ',HMAXA,HMAXB,HMAXA-HMAXB
        CALL ZHFLSH(6)
        STOP
      ENDIF
C
      READ(21,'(A)') CLINE
      I = INDEX(CLINE,'=')
      READ (CLINE(I+1:),*)   HMINB,HMAXB
      CALL ZAIORD(PLAT,MSK,.FALSE., HMINA,HMAXA, 21)
      IF     (ABS(HMINA-HMINB).GT.ABS(HMINB)*1.E-4 .OR.
     &        ABS(HMAXA-HMAXB).GT.ABS(HMAXB)*1.E-4     ) THEN
        WRITE(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - .a and .b grid files not consistent (plat):',
     &    '.a,.b min = ',HMINA,HMINB,HMINA-HMINB,
     &    '.a,.b max = ',HMAXA,HMAXB,HMAXA-HMAXB
        CALL ZHFLSH(6)
        STOP
      ENDIF
C
      CLOSE(UNIT=21)
      CALL ZAIOCL(21)
C
C     CHECK FOR AN ARCTIC BIPOLAR PATCH
C
      LARCTIC = PLON(3*IDM/4+1,JDM-1).EQ.PLON(IDM/4,JDM) .AND.
     &          PLAT(3*IDM/4+1,JDM-1).EQ.PLAT(IDM/4,JDM)
C
      WRITE(6,*)
      WRITE(6,*) 'LARCTIC = ',LARCTIC
      WRITE(6,*)
      CALL ZHFLSH(6)
C
C     INITIALIZE OUTPUT.
C
      WRITE(6,6000) 'OUTPUT:',CTITLE
      CALL ZHFLSH(6)
C
      CALL ZAIOPN('NEW', 10)
      CALL ZAIOPN('NEW', 11)
      CALL ZAIOPN('NEW', 12)
      CALL ZAIOPN('NEW', 13)
C
      CALL ZHOPEN(10, 'FORMATTED', 'NEW', 0)
      CALL ZHOPEN(11, 'FORMATTED', 'NEW', 0)
      CALL ZHOPEN(12, 'FORMATTED', 'NEW', 0)
      CALL ZHOPEN(13, 'FORMATTED', 'NEW', 0)
C
      PREAMBL(1) = CTITLE
      PREAMBL(3) = ' '
      PREAMBL(4) = ' '
      WRITE(PREAMBL(5),'(A,2I5)')
     +        'i/jdm =',
     +       IDM,JDM
      WRITE(PREAMBL(2),'(A,F10.2,A,A,F10.2,A,A,F10.6,A)')
     +        'stfr_spc =',STFR_SPC(1),' m   ',
     +        'stfr_tim =',STFR_TIM(1),' hr   ',
     +        'stfr_rms =',STFR_RMS(1),' degC/hr'
      WRITE(10,4101) PREAMBL
      WRITE(6,*)
      WRITE(6, 4101) PREAMBL
      WRITE(PREAMBL(2),'(A,F10.2,A,A,F10.2,A,A,F10.6,A)')
     +        'stfr_spc =',STFR_SPC(2),' m   ',
     +        'stfr_tim =',STFR_TIM(2),' hr   ',
     +        'stfr_rms =',STFR_RMS(2),' psu/hr'
      WRITE(11,4101) PREAMBL
      WRITE(6,*)
      WRITE(6, 4101) PREAMBL
      WRITE(PREAMBL(2),'(A,F10.2,A,A,F10.2,A,A,F10.6,A)')
     +        'stfr_spc =',STFR_SPC(3),' m   ',
     +        'stfr_tim =',STFR_TIM(3),' hr   ',
     +        'stfr_rms =',STFR_RMS(3),' m/s/hr'
      PREAMBL(3) = 'on u-grid'
      WRITE(12,4101) PREAMBL
      WRITE(6,*)
      WRITE(6, 4101) PREAMBL
      PREAMBL(3) = 'on v-grid'
      WRITE(13,4101) PREAMBL
      WRITE(6,*)
      WRITE(6, 4101) PREAMBL
      WRITE(6,*)
C
C     PROCESS ALL THE STOCHASTIC RECORDS.
C
      NREC  = NINT( (TMAX   - TSTART) / SINC ) + 1
      MREC  = NINT( (SSTART - TSTART) / SINC ) + 1
      NSTEP = NINT((SINC*86400.0d0)/DTI)
      DXREF = (PLAT(IDM/2,JDM/2+1)-PLAT(IDM/2+1,JDM/2))*111200.0
      DTI4  = DTI
C
      WRITE(6,'(a,3i4,f12.3)')
     &  'nrec,mrec,nstep,dxref = ',nrec,mrec,nstep,dxref
      WRITE(6,*)
C
      TPM(:,:,1) = 0.0
      TPM(:,:,2) = 0.0
      VPM(:,:,1) = 0.0
      VPM(:,:,2) = 0.0
C
      DO 810 KREC= 1,NREC
C
        WDAY8 = TSTART + (KREC-1)*SINC
C
        IF     (KREC.EQ.MREC) THEN
C
C         FORM INITIAL STOCHASTIC FIELDS.
C
          init = 0
          call stoch_r(init,idm,jdm,2,dxref,dti4,
     &                 stfr_spc(1),stfr_tim(1),stfr_rms(1),tpm)
          init = 0
          call stoch_v(init,idm,jdm,  dxref,dti4,
     &                 stfr_spc(3),stfr_tim(3),stfr_rms(3),vpm(1,1,1))
          init = 1
          call stoch_v(init,idm,jdm,  dxref,dti4,
     &                 stfr_spc(3),stfr_tim(3),stfr_rms(3),vpm(1,1,2))
        ELSEIF (KREC.GT.MREC) THEN
C
C         TIME ADVANCE STOCHASTIC FIELDS.
C
          init = 1
          do kstep= 1,nstep
            call stoch_r(init,idm,jdm,2,dxref,dti4,
     &                   stfr_spc(1),stfr_tim(1),stfr_rms(1),tpm)
            call stoch_v(init,idm,jdm,  dxref,dti4,
     &                   stfr_spc(3),stfr_tim(3),stfr_rms(3),vpm(1,1,1))
            call stoch_v(init,idm,jdm,  dxref,dti4,
     &                   stfr_spc(3),stfr_tim(3),stfr_rms(3),vpm(1,1,2))
          enddo !istep
        ENDIF !krec
C
C       WRITE OUT STATISTICS.
C
        CALL ARCUPD(TPM(1,1,1),IDM,JDM, LARCTIC,.TRUE.)
        CALL MINMAX(TPM(1,1,1),IDM,JDM, XMIN,XMAX)
        CALL AVERMS(TPM(1,1,1),IDM,JDM, XAVE,XRMS)
        WRITE(6,8100) 'stocht', XMIN,XMAX,XAVE,XRMS
C
        CALL ARCUPD(TPM(1,1,2),IDM,JDM, LARCTIC,.TRUE.)
        CALL MINMAX(TPM(1,1,2),IDM,JDM, XMIN,XMAX)
        CALL AVERMS(TPM(1,1,2),IDM,JDM, XAVE,XRMS)
        WRITE(6,8100) 'stochs', XMIN,XMAX,XAVE,XRMS
C
        CALL ARCUPD(VPM(1,1,1),IDM,JDM, LARCTIC,.TRUE.)
        CALL MINMAX(VPM(1,1,1),IDM,JDM, XMIN,XMAX)
        CALL AVERMS(VPM(1,1,1),IDM,JDM, XAVE,XRMS)
        WRITE(6,8100) 'stochu', XMIN,XMAX,XAVE,XRMS
C
        CALL ARCUPD(VPM(1,1,2),IDM,JDM, LARCTIC,.TRUE.)
        CALL MINMAX(VPM(1,1,2),IDM,JDM, XMIN,XMAX)
        CALL AVERMS(VPM(1,1,2),IDM,JDM, XAVE,XRMS)
        WRITE(6,8100) 'stochv', XMIN,XMAX,XAVE,XRMS
C
C       WRITE OUT HYCOM FIELDS.
C
        CALL ZAIOWR(TPM(1,1,1),MSK,.FALSE., XMIN,XMAX, 10, .FALSE.)
        WRITE(10,4122) 'stocht',WDAY8,SINC,XMIN,XMAX
        CALL ZHFLSH(10)
C
        CALL ZAIOWR(TPM(1,1,2),MSK,.FALSE., XMIN,XMAX, 11, .FALSE.)
        WRITE(11,4122) 'stochs',WDAY8,SINC,XMIN,XMAX
        CALL ZHFLSH(11)
C
        CALL ZAIOWR(VPM(1,1,1),MSK,.FALSE., XMIN,XMAX, 12, .FALSE.)
        WRITE(12,4122) 'stochu',WDAY8,SINC,XMIN,XMAX
        CALL ZHFLSH(12)
C
        CALL ZAIOWR(VPM(1,1,2),MSK,.FALSE., XMIN,XMAX, 13, .FALSE.)
        WRITE(13,4122) 'stochv',WDAY8,SINC,XMIN,XMAX
        CALL ZHFLSH(13)
C
        CALL WNDAY(WDAY8, WYR,FDY)
        WRITE(6,6350) KREC,WDAY8,FDY,NINT(WYR)
        CALL ZHFLSH(6)
  810 CONTINUE
C
      CALL ZAIOCL(10)
      CALL ZAIOCL(11)
      CALL ZAIOCL(12)
      CALL ZAIOCL(13)
      CLOSE( UNIT=10)
      CLOSE( UNIT=11)
      CLOSE( UNIT=12)
      CLOSE( UNIT=13)
C
C     SUMMARY.
C
      CALL WNDAY(TSTART, WYR,FDY)
      WRITE(6,6450) NREC,FDY,NINT(WYR),WDAY8-TSTART
      CALL ZHFLSH(6)
      STOP
C
 4101 FORMAT(A79)
 4112 FORMAT(2X,A,': range = ',1P2E16.7)
 4122 FORMAT(2X,A,': day,span,range =',F12.5,F10.6,1P2E16.7)
 6000 FORMAT(1X,A,2X,A40 //)
 6200 FORMAT(/ 1X,'MIN,MAX I COORDS = ',F8.2,',',F8.2 
     +       / 1X,'MIN,MAX J COORDS = ',F8.2,',',F8.2 /)
 6350 FORMAT(10X,'WRITING STOC RECORD',I5,
     +           '    SDAY =',F10.3,
     +            '  SDATE =',F8.3,'/',I4 /)
 6450 FORMAT(I5,' STOC RECORDS STARTING ON',F8.3,'/',I4,
     +   ' COVERING',F12.5,' DAYS')
 8100 FORMAT(1X,A,': MIN=',F12.8,' MAX=',F12.8,
     +             ' AVE=',F12.8,' RMS=',F12.8)
C     END OF PROGRAM STOCH
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

      subroutine stoch_v(init,n,m,delx,dti,
     &                   stfv_ss,stfv_ts,stfv_rms, stfv)
c
c  created 2014-04-07, Paul J Martin, NRL.
      implicit none
c
c  declare passed variables.
      integer init,n,m
      real*4  delx,dti
      real*4  stfv_ss(1),stfv_ts(1),stfv_rms(1)
      real*4  stfv(n,m)
c
c  compute/update stochastic forcing fields on 1st call for new time step.
*     write(6,*) 'v: init,n,m    =',init,n,m
*     write(6,*) 'v: delx,dti    =',delx,dti
*     call zhflsh(6)
*
      call stoch_field(init,delx,dti,n,m,1,
     &  stfv_ss,stfv_ts,stfv_rms,stfv)
c
      return
      end

      subroutine stoch_r(init,n,m,nr,delx,dti,
     &                   stfr_ss,stfr_ts,stfr_rms, stfr)
c
c  created 2014-04-07, Paul J Martin, NRL.
      implicit none
c
c  declare passed variables.
      integer init,n,m,nr
      real*4  delx,dti
      real*4  stfr_ss(nr),stfr_ts(nr),stfr_rms(nr)
      real*4  stfr(n,m,nr)
c
c  compute/update stochastic forcing fields on 1st call for new time step.
*     write(6,*) 'r: init,n,m,nr =',init,n,m,nr
*     write(6,*) 'r: delx,dti    =',delx,dti
*     call zhflsh(6)
*
      call stoch_field(init,delx,dti,n,m,nr,
     &  stfr_ss,stfr_ts,stfr_rms,stfr)
c
      return
      end

      subroutine stoch_field(init,delx,dti,n,m,nstf,
     &  stf_ss,stf_ts,stf_rms,stf)
c  subroutine to compute stochastic forcing field.
c      init   = integer flag to denote initial calculation of stochastic
c               forcing field (when =0).
c      delx   = approximate grid spacing in m.
c      dti    = time step
c      n,m    = horizontal dimensions of local tile.
c      l      = number of vertical layers +1.
c      nstf   = number of stochastic forcing fields being computed.
c      stf_ss = horizontal spatial scale for stochastic forcing (+m).
c      stf_ts = temporal scale for stochastic forcing (hr).
c      stf_rms= rms amplitude for stochastic forcing (units/hr).
c      stf    = returned stochastic forcing field.
c  created 2014-04-04, Paul J Martin, NRL.
      implicit none
c
c  declare passed variables.
      integer init,n,m,nstf
      real*4  delx,dti
      real*4  stf_ss(nstf),stf_ts(nstf),stf_rms(nstf)
      real*4  stf(n,m,nstf)
c
c  declare local temporary variables.
      integer i,j,ir,iflt,nflt
      real*4  a,b,c,d
      real*4, allocatable :: r1(:,:),r2(:,:)
      real*8  a8,b8,c8,d8
      integer indgrd
      real*4  timed,dum,amult,cint,vscale
c
c  diagnostic printout.
      if (init .eq. 0) then
        call init_random_seed
      endif
c
c  allocate needed local arrays.
      allocate( r1(n,m) )
      allocate( r2(n,m) )
c
c  loop through all the stochastic forcing fields being computed.
      do ir=1,nstf
*       write(6,*) 'ir,delx,dti =',ir,delx,dti
*       write(6,*) 'ir,stf      =',ir,stf_ss(ir),stf_ts(ir),stf_rms(ir)
*       call zhflsh(6)
c
c  calc random 2D field (note:  random_number is a system routine).
c  this routine returns uniformly-distributed random numbers from 0 to 1.
      call random_number(r1)
      r2 = 0.
c
c  average over 3x3 boxes to give an approximate normal distribution.
c  this has small tails because of small number of values (9) averaged.
c  the value "a" is subtracted to provide an ~ mean value of zero.
c  note:  value of "a" used here may need slight tweek to give mean =0.!!
      a=4.5
      do j=1,m
        do i=2,n-1
          r2(i,j)=r1(i-1,j)+r1(i,j)+r1(i+1,j)
        enddo
        r2(1,j)=r2(2,j)
        r2(n,j)=r2(n-1,j)
      enddo
      do i=1,n
        do j=2,m-1
          r1(i,j)=(r2(i,j-1)+r2(i,j)+r2(i,j+1)) - a
        enddo
        r1(i,1)=r1(i,2)
        r1(i,m)=r1(i,m-1)
      enddo
c
c  use spatial filtering to achieve desired spatial correlation scale.
c
c  estimate number of passes of spatial filter needed.
c  note:  this formula may need some refinement.!!
      nflt=int(0.300*(stf_ss(ir)/delx)**2)
c
c  apply spatial filter.
      do iflt=1,nflt
c  apply filter in x-direction.
        do j=1,m
          do i=2,n-1
            r2(i,j)=(r1(i-1,j)+r1(i+1,j)) + 2.*r1(i,j)
          enddo
          r2(1,j)=r2(2,j)
          r2(n,j)=r2(n-1,j)
        enddo
c  apply filter in y-direction.
        a=1./16.0
        do i=1,n
          do j=2,m-1
            r1(i,j)=( (r2(i,j-1)+r2(i,j+1)) + 2.*r2(i,j) )*a
          enddo
          r1(i,1)=r1(i,2)
          r1(i,m)=r1(i,m-1)
        enddo
      enddo
c
c  scale the amplitude of the stochastic forcing field.
c  note:  this formula may need some tuning.!!
      a=stf_rms(ir)*(nflt+1)**0.465/0.86366
      do j=1,m
        do i=1,n
          r1(i,j)=r1(i,j)*a
        enddo
      enddo
c
c  apply temporal correlation.
c  note that the particular temporal weightings (a,b) used here may seem
c  a bit strange, but these will maintain a fairly constant amplitude
c  for the computed stochastic forcing field (stf) with time (which is
c  desireable).
      a=dti/(3600.0*stf_ts(ir))
      c=exp(-a)
      d=sqrt(1. - exp(-2.*a))
      if (init .eq. 0) then
        do j=1,m
          do i=1,n
            stf(i,j,ir)=r1(i,j)
          enddo
        enddo
      else
*       write(6,*) 'ir,a,c,d    =',ir,a,c,d
*       call zhflsh(6)
        do j=1,m
          do i=1,n
            stf(i,j,ir)=c*stf(i,j,ir)+d*r1(i,j)
          enddo
        enddo
      endif
c
      enddo
c
c  deallocate local temporary arrays.
      deallocate(r1,r2)
c
      return
      end

      subroutine init_random_seed
c  subroutine to get random seed values for random number generator 
c  from file = /dev/urandom .
c  note that the intel workstations return n=34 and 10-digit integer 
c  seed values with both + and - signs.
c  created 2014-06-27, Paul J Martin, NRL.
      implicit none
c
      integer n,istat,i,datetime(8)
      integer, allocatable :: iseed(:)
c
      call random_seed(size = n)
*     write(6,*) 'rs: n =',n
*     call zhflsh(6)
      allocate( iseed(n) )
c
      call date_and_time(values=datetime)
      do i= 1,min(8,n)
        iseed(i) = datetime(9-i)
      enddo
      do i= 9,n
        iseed(i) = i
      enddo
*     write(6,*) 'dtime =',datetime
      write(6,*) 'iseed =',iseed(1:n)
      call zhflsh(6)
      call random_seed(put = iseed)
c
      deallocate( iseed )
      return
      end

      PROGRAM CICE_RESTART_MASK
      IMPLICIT NONE
C
C  cice_restart_mask - Usage:  cice_restart_mask rold grid idm jdm [itest jtest ] nc nl rnew
C
C  Changes land/sea mask for a CICE restart
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  April 2015.
C
      LOGICAL, ALLOCATABLE :: TMASK(:,:),UMASK(:,:),ZMASK(:,:)
      REAL*8,  ALLOCATABLE :: A8(:,:),H8(:,:)
C
      INTEGER       IARGC
      INTEGER       NARG
      CHARACTER*240 CARG
C
      INTEGER       IDM,JDM,ITEST,JTEST,NC,NL
      CHARACTER*240 CFILE1,CFILEO,CFILEG
      INTEGER       I,II,J,K,N,ISTEP,IOS,NRECL
      REAL*8        RDAY,FDAY
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.7) THEN
        CALL GETARG(1,CFILE1)
        CALL GETARG(2,CFILEG)
        CALL GETARG(3,CARG)
        READ(CARG,*) IDM
        CALL GETARG(4,CARG)
        READ(CARG,*) JDM
        ITEST = 0
        JTEST = 0
        CALL GETARG(5,CARG)
        READ(CARG,*) NC
        CALL GETARG(6,CARG)
        READ(CARG,*) NL
        CALL GETARG(7,CFILEO)
      ELSEIF (NARG.EQ.9) THEN
        CALL GETARG(1,CFILE1)
        CALL GETARG(2,CFILEG)
        CALL GETARG(3,CARG)
        READ(CARG,*) IDM
        CALL GETARG(4,CARG)
        READ(CARG,*) JDM
        CALL GETARG(5,CARG)
        READ(CARG,*) ITEST
        CALL GETARG(6,CARG)
        READ(CARG,*) JTEST
        CALL GETARG(7,CARG)
        READ(CARG,*) NC
        CALL GETARG(8,CARG)
        READ(CARG,*) NL
        CALL GETARG(9,CFILEO)
      ELSE
        WRITE(6,*)
     &    'Usage:  '//
     &    'cice_restart_mask rold grid idm jdm [itest jtest] nc nl rnew'
        CALL EXIT(1)
      ENDIF
C
      ALLOCATE( A8(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in cice_restart: could not allocate ',
     +             IDM*JDM,' 8-byte words'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( H8(0:IDM+1,0:JDM+1), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in cice_restart: could not allocate ',
     +             (IDM+2)*(JDM+2),' 8-byte words'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( TMASK(IDM,JDM),
     &          UMASK(IDM,JDM),
     &          ZMASK(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in cice_restart: could not allocate ',
     +             3*IDM*JDM,' 4-byte words'
        CALL EXIT(2)
      ENDIF
C
C     TMASK AND UMASK
C
      INQUIRE( IOLENGTH=NRECL) A8(:,:)
      OPEN(UNIT=22, FILE=CFILEG, FORM='UNFORMATTED', STATUS='OLD',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILEG)
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
      WRITE(6,*) 'IDM,JDM = ',IDM,JDM
      WRITE(6,*) 'NRECL   = ',IDM,JDM
      WRITE(6,*) 'reading from ',TRIM(CFILEG)
      A8(:,:) = -1.0d0
      READ( 22,REC=1) A8(:,:)
      CLOSE(22)
      H8(1:IDM,1:JDM) = A8(:,:)
      DO I= 1,IDM
        II = IDM-MOD(I-1,IDM)
        H8(I,0)     = 0.0
        H8(I,JDM+1) = H8(II,JDM)
      ENDDO !i
      DO J= 0,JDM+1
        H8(    0,J) = H8(IDM,J)
        H8(IDM+1,J) = H8(  1,J)
      ENDDO !j
C
      DO J= 1,JDM
        DO I= 1,IDM
          TMASK(I,J) =      H8(I,  J)     .GT. 0.5D0
          UMASK(I,J) = MIN( H8(I,  J), 
     &                      H8(I+1,J), 
     &                      H8(I,  J+1), 
     &                      H8(I+1,J+1) ) .GT. 0.5D0
          ZMASK(I,J) = .TRUE.  !no masking
        ENDDO !j
      ENDDO !i
C
      IF     (MAX(ITEST,JTEST).GT.0) THEN
        WRITE(6,*) 'A8   ',ITEST,JTEST,   A8(ITEST,JTEST)
        WRITE(6,*) 'H8   ',ITEST,JTEST,   H8(ITEST,JTEST)
        WRITE(6,*) 'TMASK',ITEST,JTEST,TMASK(ITEST,JTEST)
        WRITE(6,*) 'UMASK',ITEST,JTEST,UMASK(ITEST,JTEST)
      ENDIF
C
C --- RESTART
C
      OPEN(UNIT=11, FILE=CFILEO, FORM='UNFORMATTED', STATUS='NEW',
     +         IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILEO)
        write(6,*) 'ios   = ',ios
        CALL EXIT(3)
      ENDIF
      OPEN(UNIT=21, FILE=CFILE1, FORM='UNFORMATTED', STATUS='OLD',
     +         IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILE1)
        write(6,*) 'ios   = ',ios
        CALL EXIT(3)
      ENDIF
C
C     HEADER RECORD
C
      READ( 21) ISTEP,RDAY,FDAY
      WRITE(11) ISTEP,RDAY,FDAY
      WRITE(6,'(A,I9,5X,2F20.6)')
     &  'step,time,time_forc = ',ISTEP,RDAY/86400.d0,FDAY/86400.d0
C
C     ARRAY RECORDS
C
      DO N= 1,NC
        write(6,*) 'cat ',n,' min/max area, vol ice, vol snow, Tsfc'
        DO K= 1,4
          CALL R21_W11(A8,TMASK,IDM,JDM,ITEST,JTEST)
        ENDDO !k
      ENDDO !n
      write(6,*) 'min/max eicen for each layer and category'
      DO K= 1,NL*NC
        CALL R21_W11(A8,TMASK,IDM,JDM,ITEST,JTEST)
      ENDDO !k
      write(6,*) 'min/max esnon for each layer and category'
      DO K= 1,NC
        CALL R21_W11(A8,TMASK,IDM,JDM,ITEST,JTEST)
      ENDDO !k
      write(6,*) 'min/max velocity components'
      CALL R21_W11(A8,UMASK,IDM,JDM,ITEST,JTEST)
      CALL R21_W11(A8,UMASK,IDM,JDM,ITEST,JTEST)
      write(6,*) 'radiation fields (no mask)'
      CALL R21_W11(A8,ZMASK,IDM,JDM,ITEST,JTEST)
      CALL R21_W11(A8,ZMASK,IDM,JDM,ITEST,JTEST)
      CALL R21_W11(A8,ZMASK,IDM,JDM,ITEST,JTEST)
      CALL R21_W11(A8,ZMASK,IDM,JDM,ITEST,JTEST)
      CALL R21_W11(A8,ZMASK,IDM,JDM,ITEST,JTEST)
      write(6,*) 'min/max ocean stress components'
      CALL R21_W11(A8,TMASK,IDM,JDM,ITEST,JTEST)
      CALL R21_W11(A8,TMASK,IDM,JDM,ITEST,JTEST)
      write(6,*) 'internal stress components'
      CALL R21_W11(A8,TMASK,IDM,JDM,ITEST,JTEST)
      CALL R21_W11(A8,TMASK,IDM,JDM,ITEST,JTEST)
      CALL R21_W11(A8,TMASK,IDM,JDM,ITEST,JTEST)
      CALL R21_W11(A8,TMASK,IDM,JDM,ITEST,JTEST)
      CALL R21_W11(A8,TMASK,IDM,JDM,ITEST,JTEST)
      CALL R21_W11(A8,TMASK,IDM,JDM,ITEST,JTEST)
      CALL R21_W11(A8,TMASK,IDM,JDM,ITEST,JTEST)
      CALL R21_W11(A8,TMASK,IDM,JDM,ITEST,JTEST)
      CALL R21_W11(A8,TMASK,IDM,JDM,ITEST,JTEST)
      CALL R21_W11(A8,TMASK,IDM,JDM,ITEST,JTEST)
      CALL R21_W11(A8,TMASK,IDM,JDM,ITEST,JTEST)
      CALL R21_W11(A8,TMASK,IDM,JDM,ITEST,JTEST)
      write(6,*) 'ice mask for dynamics'
      CALL R21_W11(A8,UMASK,IDM,JDM,ITEST,JTEST)
      write(6,*) 'min/max sst, frzmlt'
      CALL R21_W11(A8,TMASK,IDM,JDM,ITEST,JTEST)
      CALL R21_W11(A8,TMASK,IDM,JDM,ITEST,JTEST)
C
      CLOSE(11)
      CLOSE(21)
C
      END
      SUBROUTINE R21_W11(A8,AM,IDM,JDM,ITEST,JTEST)
      IMPLICIT NONE
C
      INTEGER IDM,JDM,ITEST,JTEST
      REAL*8  A8(IDM,JDM)
      LOGICAL AM(IDM,JDM)
C
C     READ ONE RECORD ON UNIT 21,
C     APPLY THE MASK
C     WRITE IT OUT ON UNIT 11
C
      INTEGER I,J,IOS
      REAL*8  A8MAX,A8MIN
C
      READ(21,IOSTAT=IOS) A8
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'can''t read unit 21'
        CALL EXIT(4)
      ENDIF
      A8MIN = 0.0
      A8MAX = 0.0
      DO J= 1,JDM
        DO I= 1,IDM
          IF     (.NOT.AM(I,J)) THEN
            A8(I,J) = 0.D0
          ELSE
            A8MIN = MIN( A8MIN, A8(I,J) )
            A8MAX = MAX( A8MAX, A8(I,J) )
          ENDIF
        ENDDO !j
      ENDDO !i
      WRITE(11) A8
      WRITE(6,*) A8MIN,A8MAX
      IF     (MAX(ITEST,JTEST).GT.0) THEN
        WRITE(6,*) ITEST,JTEST,A8(ITEST,JTEST)
      ENDIF
      END

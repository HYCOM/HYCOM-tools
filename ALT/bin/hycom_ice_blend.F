      PROGRAM HYCOM_ICE_BLEND
      IMPLICIT NONE
C
C  hycom_ice_blend - Usage:  hycom_ice_blend ice_a.a ice_m.a idm jdm ice_b.a type [maxice_a maxice_m]
C
C                  Blends sea ice analysis and model fields
C
C  input  ice is zero for no ice and maxice (default 1.0) for all ice.
C  output ice is zero for no ice and 1.0 for all ice.
C  type is the blending algorithm
C    0: ARCc0.08 approach
C    1: GLBb0.08 approach
C    2: ARCc0.08 approach, edge_m = 0.45
C
C  *.a is assumed to contain idm*jdm 32-bit IEEE real values for
C   each array, in standard f77 element order, followed by padding
C   to a multiple of 4096 32-bit words, but otherwise with no control
C   bytes/words, and input values of 2.0**100 indicating a data void.
C
C  if the "ice_a.a" and "ice_m.a" filename do not end with ".a" they
C   are assumed to be raw files, i.e. no padding, with -9999.0 as their
C   data void marker.
C  note that ice_b.a is always a hycom ".a" file.
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  February 2014.
C
      REAL*4, ALLOCATABLE :: ICE_A(:,:),ICE_M(:,:)
      REAL*4              :: PAD(4096)
      INTEGER       IOS
      INTEGER       IARGC
      INTEGER       NARG
      CHARACTER*240 CARG
C
      REAL*4        MAXICE_A,MAXICE_M
      INTEGER       ITYPE,IDM,JDM,NPAD
      CHARACTER*240 CFILEA,CFILEM,CFILEO
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.8) THEN
        CALL GETARG(1,CFILEA)
        CALL GETARG(2,CFILEM)
        CALL GETARG(3,CARG)
        READ(CARG,*) IDM
        CALL GETARG(4,CARG)
        READ(CARG,*) JDM
        CALL GETARG(5,CFILEO)
        CALL GETARG(6,CARG)
        READ(CARG,*) ITYPE
        CALL GETARG(7,CARG)
        READ(CARG,*) MAXICE_A
        CALL GETARG(8,CARG)
        READ(CARG,*) MAXICE_M
      ELSEIF (NARG.EQ.6) THEN
        CALL GETARG(1,CFILEA)
        CALL GETARG(2,CFILEM)
        CALL GETARG(3,CARG)
        READ(CARG,*) IDM
        CALL GETARG(4,CARG)
        READ(CARG,*) JDM
        CALL GETARG(5,CFILEO)
        CALL GETARG(6,CARG)
        READ(CARG,*) ITYPE
        MAXICE_A = 1.0
        MAXICE_M = 1.0
      ELSE
        WRITE(6,'(2a)')
     &   'Usage: hycom_ice_blend ice_a.a ice_m.a idm jdm ice_b.a type',
     &                           ' [maxice_a maxice_m]'
        CALL EXIT(1)
      ENDIF
C
      NPAD = 4096 - MOD(IDM*JDM,4096)
      IF     (NPAD.EQ.4096) THEN
        NPAD = 0
      ENDIF
C
      ALLOCATE( ICE_M(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_ice_blend: could not allocate 1st ',
     +             IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( ICE_A(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_ice_blend: could not allocate 2nd ',
     +             IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
C
      CALL ICE_BLEND(ICE_M,ICE_A,IDM,JDM,PAD,NPAD, 
     &               ITYPE, MAXICE_A,MAXICE_M, CFILEA,CFILEM,CFILEO)
      CALL EXIT(0)
      END
      SUBROUTINE ICE_BLEND(ICE_M,ICE_A,IDM,JDM,PAD,NPAD,
     &                     ITYPE, MAXICE_A,MAXICE_M,
     &                     CFILEA,CFILEM,CFILEO)
      IMPLICIT NONE
C
      REAL*4     SPVAL
      PARAMETER (SPVAL =2.0**100)
C
      CHARACTER*240 CFILEA,CFILEM,CFILEO
      INTEGER       IDM,JDM,NPAD,ITYPE
      REAL*4        ICE_M(IDM,JDM),ICE_A(IDM,JDM),PAD(NPAD)
      REAL*4        MAXICE_A,MAXICE_M
C
C     MOST OF WORK IS DONE HERE.
C
      CHARACTER*18 CASN
      INTEGER      I,J,K,IOS,NRECL,MRECL_A,MRECL_M
      REAL*4       AMN,AMX,RICE_A,RICE_M,S_A,S_M,SPVALA,SPVALM,
     &             EDGE_M,EDGE_A,Q,
     &             edge_om,diff_om,aice_o,aice_m
#ifdef CRAY
      INTEGER*8    IU8,IOS8
#endif
C
      IF     (NPAD.EQ.0) THEN
        INQUIRE(IOLENGTH=NRECL) ICE_M
      ELSE
        INQUIRE(IOLENGTH=NRECL) ICE_M,PAD
        PAD(:) = SPVAL
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(PAD,NPAD)
#endif
      ENDIF
      I = LEN_TRIM(CFILEA)
      IF     (CFILEA(I-1:I).NE.'.a') THEN !raw file
        INQUIRE(IOLENGTH=MRECL_A) ICE_A
        SPVALA = -9999.0
      ELSE !.a file
        MRECL_A = NRECL
        SPVALA  = SPVAL
      ENDIF
      I = LEN_TRIM(CFILEM)
      IF     (CFILEM(I-1:I).NE.'.a') THEN !raw file
        INQUIRE(IOLENGTH=MRECL_M) ICE_M
        SPVALM = -9999.0
      ELSE !.a file
        MRECL_M = NRECL
        SPVALM  = SPVAL
      ENDIF
#ifdef CRAY
#ifdef t3e
      IF     (MOD(MRECL_A,4096).EQ.0) THEN
        WRITE(CASN,8000) MRECL_A/4096
 8000   FORMAT('-F cachea:',I4.4,':1:0')
        IU8 = 11
        CALL ASNUNIT(IU8,CASN,IOS8)
        IF     (IOS8.NE.0) THEN
          write(6,*) 'Error: can''t asnunit 11'
          write(6,*) 'ios  = ',ios8
          write(6,*) 'casn = ',casn
          CALL EXIT(5)
        ENDIF
      ENDIF
      IF     (MOD(MRECL_M,4096).EQ.0) THEN
        WRITE(CASN,8000) MRECL_M/4096
        IU8 = 12
        CALL ASNUNIT(IU8,CASN,IOS8)
        IF     (IOS8.NE.0) THEN
          write(6,*) 'Error: can''t asnunit 12'
          write(6,*) 'ios  = ',ios8
          write(6,*) 'casn = ',casn
          CALL EXIT(5)
        ENDIF
      ENDIF
      IF     (MOD(NRECL,4096).EQ.0) THEN
        WRITE(CASN,8000) NRECL/4096
        IU8 = 21
        CALL ASNUNIT(IU8,CASN,IOS8)
        IF     (IOS8.NE.0) THEN
          write(6,*) 'Error: can''t asnunit 21'
          write(6,*) 'ios  = ',ios8
          write(6,*) 'casn = ',casn
          CALL EXIT(5)
        ENDIF
      ENDIF
#else
      CALL ASNUNIT(11,'-F syscall -N ieee',IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t asnunit 11'
        write(6,*) 'ios = ',ios
        CALL EXIT(5)
      ENDIF
      CALL ASNUNIT(12,'-F syscall -N ieee',IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t asnunit 12'
        write(6,*) 'ios = ',ios
        CALL EXIT(5)
      ENDIF
      CALL ASNUNIT(21,'-F syscall -N ieee',IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t asnunit 21'
        write(6,*) 'ios = ',ios
        CALL EXIT(5)
      ENDIF
#endif
#endif
      OPEN(UNIT=11, FILE=CFILEA, FORM='UNFORMATTED', STATUS='OLD',
     +         ACCESS='DIRECT', RECL=MRECL_A, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILEA)
        write(6,*) 'ios     = ',ios
        write(6,*) 'mrecl_a = ',mrecl_a
        CALL EXIT(3)
      ENDIF
      OPEN(UNIT=12, FILE=CFILEM, FORM='UNFORMATTED', STATUS='OLD',
     +         ACCESS='DIRECT', RECL=MRECL_M, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILEM)
        write(6,*) 'ios     = ',ios
        write(6,*) 'mrecl_m = ',mrecl_m
        CALL EXIT(3)
      ENDIF
      OPEN(UNIT=21, FILE=CFILEO, FORM='UNFORMATTED', STATUS='NEW',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILEO)
        write(6,*) 'ios     = ',ios
        write(6,*) 'nrecl   = ',nrecl
        CALL EXIT(3)
      ENDIF
C
      RICE_A = 1.0/MAXICE_A
      RICE_M = 1.0/MAXICE_M
      DO 110 K= 1,9999
        READ(11,REC=K,IOSTAT=IOS) ICE_A
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(ICE_A,IDM*JDM)
#endif
        IF     (IOS.NE.0) THEN
          IF     (K.EQ.1) THEN
            WRITE(6,*) 'can''t read ',TRIM(CFILEA)
            CALL EXIT(4)
          ELSE
            GOTO 1110
          ENDIF
        ENDIF
        READ(12,REC=K,IOSTAT=IOS) ICE_M
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(ICE_M,IDM*JDM)
#endif
        IF     (IOS.NE.0) THEN
          IF     (K.EQ.1) THEN
            WRITE(6,*) 'can''t read ',TRIM(CFILEM)
            CALL EXIT(4)
          ELSE
            GOTO 1110
          ENDIF
        ENDIF
C
        AMN =  SPVAL
        AMX = -SPVAL

        IF     (ITYPE.EQ.0 .OR. ITYPE.EQ.2) THEN
C
C         ARCc0.08 approach
C
          IF     (ITYPE.EQ.0) THEN
            EDGE_A = 0.25  ! edge of obs.  ice
            EDGE_M = 0.15  ! edge of model ice
          ELSE   !ITYPE.EQ.2
            EDGE_A = 0.25  ! edge of obs.  ice
            EDGE_M = 0.45  ! edge of model ice
          ENDIF
          DO J= 1,JDM
            DO I= 1,IDM
              IF     (ICE_M(I,J).NE.SPVALM .AND.
     &                ICE_A(I,J).NE.SPVALA      ) THEN
                S_A = ICE_A(I,J)*RICE_A
                S_M = ICE_M(I,J)*RICE_M
                ICE_M(I,J) = S_M  !default is model, i.e. do nothing
                IF     (S_A.LE.S_M) THEN
                  IF (S_A.LT.0.01) THEN
                    ICE_M(I,J) = 0.0
                  ELSEIF (S_A.LT.2.0*EDGE_A)then
                    ! --- target ice conc. is obs. when obs. <   EDGE_A
                    !                    and model when obs. > 2*EDGE_A
                    Q = MAX( 0.0, (S_A-EDGE_A)/EDGE_A )  ! 0 to 1
                    ICE_M(I,J) = S_A + Q*(S_M - S_A)
                  ENDIF
                ELSEIF (S_M.LT.2.0*EDGE_M .AND. S_A.GT.0.01) then
                  ! --- target ice conc. is obs. when model <   EDGE_M
                  !                    and model when model > 2*EDGE_M
                  Q = max( 0.0, (S_M-EDGE_M)/EDGE_M )  ! 0 to 1
                  ICE_M(I,J) = S_A + Q*(S_M - S_A)
                ENDIF
                AMN = MIN( AMN, ICE_M(I,J) )
                AMX = MAX( AMX, ICE_M(I,J) )
              ELSE
                ICE_M(I,J) = SPVAL
              ENDIF
            ENDDO !i
          ENDDO !j
        ELSE
C
C         GLBb0.08 approach
C
          edge_om = 0.20  ! nominal ice edge zone
          diff_om = 0.10  ! allowed model vs obs difference
          DO J= 1,JDM
            DO I= 1,IDM
              IF     (ICE_M(I,J).NE.SPVALM .AND.
     &                ICE_A(I,J).NE.SPVALA      ) THEN
                aice_o = ICE_A(I,J)*RICE_A  !obs.  ice concentration
                aice_m = ICE_M(I,J)*RICE_M  !model ice concentration
                ICE_M(I,J) = aice_m  !default is model, i.e. do nothing
                if     (aice_o.gt.0.01 .and.
     &                  abs(aice_o-aice_m).le.0.01) then
                   ! model and obs are very close - do nothing
                elseif (min(aice_o,aice_m).ge.edge_om .and.
     &                  abs(aice_o-aice_m).le.diff_om) then
                   ! model and obs are close enough - do nothing
                elseif (aice_o.eq.aice_m) then
                   !  model and obs are identical - do nothing
                elseif (aice_o.lt.aice_m) then
                  if (aice_o.lt.0.01)then
                    ! --- remove all ice ---
                    ICE_M(I,J) = 0.0
                  else !aice_o.ge.0.01
                    if     (aice_o.lt.edge_om) then
                       ! --- target ice conc. is obs.
                       ICE_M(I,J) = aice_o
                    else !aice_m-aice_o.gt.diff_om
                       ! --- target ice conc. is obs.+diff_om
                       ICE_M(I,J) = aice_o + diff_om
                    endif
                  endif
                elseif (aice_o.gt.0.01) then  ! .and. aice_o.gt.aice_m
                  if     (aice_m.lt.edge_om) then
                     ! --- target ice conc. is obs.
                     ICE_M(I,J) = aice_o
                  else !aice_o-aice_m.gt.diff_om
                     ! --- target ice conc. is obs.-diff_om
                     ICE_M(I,J) = aice_o - diff_om
                  endif
                endif !aice_o vs aice_m
                AMN = MIN( AMN, ICE_M(I,J) )
                AMX = MAX( AMX, ICE_M(I,J) )
              ELSE
                ICE_M(I,J) = SPVAL
              ENDIF
            ENDDO !i
          ENDDO !j
       ENDIF !itype
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(ICE_M,IDM*JDM)
#endif
        IF     (NPAD.EQ.0) THEN
          WRITE(21,REC=K,IOSTAT=IOS) ICE_M
        ELSE
          WRITE(21,REC=K,IOSTAT=IOS) ICE_M,PAD
        ENDIF
        WRITE(6,'(a,2f8.4)')
     &     'iceblnd: min, max = ',AMN,AMX
  110 CONTINUE
 1110 CONTINUE
      WRITE(6,*) 
      WRITE(6,*) K-1,' FIELDS PROCESSED'
      WRITE(6,*) 
C
      CLOSE(11)
      CLOSE(12)
      CLOSE(21)
C
      RETURN
      END

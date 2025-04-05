      PROGRAM HYCOM_MASS_CORR
      IMPLICIT NONE
C
C  hycom_mass_corr - Usage:  hycom_mass_corr fin.a fout.a mean [[grid.a] mask.a]
C  hycom_mass_bal  - Usage:  hycom_mass_bal  fin.a fout.a mean [[grid.a] mask.a]
C
C                 Outputs all fields in fin.a corrected to the specified mean
C                 _corr via a region-wide constant offset
C                 _bal  by scaling down +ve or -ve anomalies
C
C                 grid.a is a hycom grid file, default regional.grid.a.
C                 Note that the corresponding grid.b must also exist.
C
C                 mask.a contains an optional mask array, that turns
C                 off sampling anywhere it is 2.0**100 (data void).
C                 Note that grid.a must be specified if mask.a is used.
C
C                 idm,jdm are taken from grid.a and the input array
C                 is assumed to be on the p-grid (i.e. plat from grid.a
C                 represents the field's latitude).
C
C  fin.a is assumed to contain idm*jdm 32-bit IEEE real values for
C   each array, in standard f77 element order, followed by padding
C   to a multiple of 4096 32-bit words, but otherwise with no control
C   bytes/words, and input values of 2.0**100 indicating a data void.
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  January 2001.
C
      REAL*4, ALLOCATABLE :: A(:,:),AMSK(:,:),PSC2(:,:)
      REAL*4              :: PAD(4096)
      INTEGER       IOS,L
      INTEGER       IARGC
      INTEGER       NARG
      CHARACTER*240 CARG
C
      LOGICAL       LBAL
      INTEGER       IDM,JDM,NPAD
      REAL*4        SMEAN
      CHARACTER*6   CVARIN
      CHARACTER*240 CFILE,CFILEG,CFILEM,CFILEB,CFILEO
C
C     READ ARGUMENTS.
C
      CALL GETARG(0,CARG)
      L = LEN_TRIM(CARG)
      LBAL = CARG(L-3:L).EQ.'_bal'
C
      NARG = IARGC()
C
      IF     (NARG.EQ.3) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CFILEO)
        CALL GETARG(3,CARG)
        READ(CARG,*) SMEAN
        CFILEG = 'regional.grid.a'
        CFILEM = 'NONE'
      ELSEIF (NARG.EQ.4) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CFILEO)
        CALL GETARG(3,CARG)
        READ(CARG,*) SMEAN
        CALL GETARG(4,CFILEG)
        CFILEM = 'NONE'
      ELSEIF (NARG.EQ.5) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CFILEO)
        CALL GETARG(3,CARG)
        READ(CARG,*) SMEAN
        CALL GETARG(4,CFILEG)
        CALL GETARG(5,CFILEM)
      ELSEIF (LBAL) THEN
        WRITE(6,*) 
     +   'Usage:  hycom_mass_bal fin.a fout.a mean [grid.a] [mask.a]'
        CALL EXIT(1)
      ELSE
        WRITE(6,*) 
     +   'Usage:  hycom_mass_corr fin.a fout.a mean [grid.a] [mask.a]'
        CALL EXIT(1)
      ENDIF
C
C     GET IDM,JDM FROM grid.b.
C
      CFILEB = CFILEG(1:LEN_TRIM(CFILEG)-1) // 'b'
C
      OPEN(UNIT=11,FILE=CFILEB,FORM='FORMATTED',
     &     STATUS='OLD',ACTION='READ')
C
      READ( 11,*) IDM,CVARIN
      IF (CVARIN.NE.'idm   ') THEN
        IF      (LBAL) THEN
          WRITE(6,*) 'hycom_mass_bal: bad header file ',
     &               TRIM(CFILEB)
        ELSE
          WRITE(6,*) 'hycom_mass_corr: bad header file ',
     &               TRIM(CFILEB)
        ENDIF
        CALL EXIT(2)
      ENDIF
      READ( 11,*) JDM,CVARIN
      IF (CVARIN.NE.'jdm   ') THEN
        IF      (LBAL) THEN
          WRITE(6,*) 'hycom_mass_bal: bad header file ',
     &               TRIM(CFILEB)
        ELSE
          WRITE(6,*) 'hycom_mass_corr: bad header file ',
     &               TRIM(CFILEB)
        ENDIF
        CALL EXIT(2)
      ENDIF
C
      CLOSE(UNIT=11)
C
      NPAD = 4096 - MOD(IDM*JDM,4096)
      IF     (NPAD.EQ.4096) THEN
        NPAD = 0
      ENDIF
C
      ALLOCATE( A(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_zonal: could not allocate ',
     +             IDM*JDM,' words for A'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( AMSK(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_zonal: could not allocate ',
     +             IDM*JDM,' words for AMSK'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( PSC2(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_zonal: could not allocate ',
     +             IDM*JDM,' words for PSC2'
        CALL EXIT(2)
      ENDIF
C
      CALL MASS(A,AMSK,PSC2,IDM,JDM,PAD,NPAD, LBAL,
     +          SMEAN, CFILE,CFILEO,CFILEG,CFILEM)
      CALL EXIT(0)
 5000 FORMAT(I4)
      END
      SUBROUTINE MASS(A,AMSK,PSC2,IDM,JDM, PAD,NPAD, LBAL,
     +                SMEAN, CFILE,CFILEO,CFILEG,CFILEM)
      IMPLICIT NONE
C
      REAL*4     SPVAL
      PARAMETER (SPVAL=2.0**100)
C
      CHARACTER*240 CFILE,CFILEO,CFILEG,CFILEM
      LOGICAL      LBAL
      INTEGER      IDM,JDM,NPAD
      REAL*4       A(IDM,JDM),AMSK(IDM,JDM),PSC2(IDM,JDM),PAD(NPAD),
     +             SMEAN
C
C     MOST OF WORK IS DONE HERE.
C
#ifdef sun
      INTEGER      IR_ISNAN
C
#endif
      CHARACTER*18 CASN
      INTEGER      I,J,K,IOS,NRECL
      REAL*4       AMN,AMX,AVE,AVEM,AVEP,Q
      REAL*8       SUMA,SUMP,SUMT,SUMTR(2)
#ifdef CRAY
      INTEGER*8    IU8,IOS8
#endif
C
C     INPUT MASK ARRAY.
C
      INQUIRE( IOLENGTH=NRECL) A,PAD
C
      IF     (CFILEM.EQ.'NONE') THEN
        AMSK(:,:) = 0.0
      ELSE
        OPEN(UNIT=11, FILE=CFILEM, FORM='UNFORMATTED', STATUS='OLD',
     +           ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
        IF     (IOS.NE.0) THEN
          write(6,*) 'Error: can''t open ',CFILEM(1:LEN_TRIM(CFILEM))
          write(6,*) 'ios   = ',ios
          write(6,*) 'nrecl = ',nrecl
          CALL EXIT(3)
        ENDIF
C
        READ(11,REC=1,IOSTAT=IOS) AMSK
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(AMSK,IDM*JDM)
#endif
        IF     (IOS.NE.0) THEN
          WRITE(6,*) 'can''t read ',CFILEM(1:LEN_TRIM(CFILEM))
          CALL EXIT(4)
        ENDIF
      ENDIF
C
      CLOSE(UNIT=11)
C
C     INPUT GRID ARRAYS.
C
      OPEN(UNIT=11, FILE=CFILEG, FORM='UNFORMATTED', STATUS='OLD',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',CFILEG(1:LEN_TRIM(CFILEG))
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
C
      READ(11,REC=10,IOSTAT=IOS) A     ! pscx
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(A,   IDM*JDM)
#endif
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'can''t read ',CFILEG(1:LEN_TRIM(CFILEG))
        CALL EXIT(4)
      ENDIF
C
      READ(11,REC=11,IOSTAT=IOS) PSC2  ! pscy
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(PSC2,IDM*JDM)
#endif
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'can''t read ',CFILEG(1:LEN_TRIM(CFILEG))
        CALL EXIT(4)
      ENDIF
C
      PSC2(:,:) = PSC2(:,:)*A(:,:)  ! psc2=pscx*pscy
C
      CLOSE(UNIT=11)
C
C     OPEN INPUT FILE.
C
#ifdef CRAY
#ifdef t3e
      IF     (MOD(NRECL,4096).EQ.0) THEN
        WRITE(CASN,8000) NRECL/4096
 8000   FORMAT('-F cachea:',I4.4,':1:0')
        IU8 = 11
        CALL ASNUNIT(IU8,CASN,IOS8)
        IF     (IOS8.NE.0) THEN
          write(6,*) 'Error: can''t asnunit 11'
          write(6,*) 'ios  = ',ios8
          write(6,*) 'casn = ',casn
          CALL EXIT(5)
        ENDIF
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
      CALL ASNUNIT(21,'-F syscall -N ieee',IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t asnunit 21'
        write(6,*) 'ios = ',ios
        CALL EXIT(5)
      ENDIF
#endif
#endif
      OPEN(UNIT=11, FILE=CFILE, FORM='UNFORMATTED', STATUS='OLD',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',CFILE(1:LEN_TRIM(CFILE))
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
      OPEN(UNIT=21, FILE=CFILEO, FORM='UNFORMATTED', STATUS='NEW',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',CFILEO(1:LEN_TRIM(CFILEO))
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
C
      IF     (NPAD.GT.0) THEN
        PAD(:) = SPVAL
      ENDIF
C
C     LOOP OVER ALL FIELDS.
C
      DO K= 1,9999
        READ(11,REC=K,IOSTAT=IOS) A
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
        IF     (IOS.NE.0) THEN
          IF     (K.EQ.1) THEN
            WRITE(6,*) 'can''t read ',CFILE(1:LEN_TRIM(CFILE))
            CALL EXIT(4)
          ELSE
            EXIT !k
          ENDIF
        ENDIF
C
        AMN  =  SPVAL
        AMX  = -SPVAL
        SUMT = 0.0
        SUMP = 0.0
        SUMA = 0.0
        DO J= 1,JDM
          DO I= 1,IDM
            IF     (   A(I,J).NE.SPVAL .AND.
     +              AMSK(I,J).NE.SPVAL      ) THEN
              AMN  = MIN( AMN, A(I,J) )
              AMX  = MAX( AMX, A(I,J) )
              IF     (LBAL .AND. A(I,J)-SMEAN.GE.0.0) THEN
                SUMP = SUMP + A(I,J)*PSC2(I,J)
              ENDIF
              SUMA = SUMA + A(I,J)*PSC2(I,J)
              SUMT = SUMT +        PSC2(I,J)
            ENDIF
          ENDDO
        ENDDO
        IF     (K.EQ.1) THEN
          SUMTR(1) = SUMT
          SUMTR(2) = SUMT
        ELSE  !k>1
          SUMTR(1) = MIN( SUMT, SUMTR(1) )
          SUMTR(2) = MAX( SUMT, SUMTR(2) )
        ENDIF !k==1:else
C
        IF     (SUMT.NE.0.0) THEN
          IF     (LBAL) THEN
            AVE  = SUMA/SUMT - SMEAN
            AVEP = SUMP/SUMT
            AVEM = AVE - AVEP  !ave = avep + avem
            WRITE(6,'(a,1p3g16.6)')
     &      'ave,avep,avem =',AVE,AVEP,AVEM
            IF     (AVE.GT.0.0) THEN
C             scale down the positive anomalies
              AMN =  SPVAL
              AMX = -SPVAL
              Q   = -AVEM/AVEP
              DO J= 1,JDM
                DO I= 1,IDM
                  IF     (   A(I,J).NE.SPVAL .AND.
     +                    AMSK(I,J).NE.SPVAL      ) THEN
                    IF     (A(I,J).GT.SMEAN) THEN
                      A(I,J) = SMEAN + Q*(A(I,J) - SMEAN)
                    ENDIF
                    AMN  = MIN( AMN, A(I,J) )
                    AMX  = MAX( AMX, A(I,J) )
                  ENDIF
                ENDDO
              ENDDO
            ELSEIF (AVE.LT.0.0) THEN
C             scale down the negative anomalies
              AMN =  SPVAL
              AMX = -SPVAL
              Q   = -AVEP/AVEM
              DO J= 1,JDM
                DO I= 1,IDM
                  IF     (   A(I,J).NE.SPVAL .AND.
     +                    AMSK(I,J).NE.SPVAL      ) THEN
                    IF     (A(I,J).LT.SMEAN) THEN
                      A(I,J) = SMEAN + Q*(A(I,J) - SMEAN)
                    ENDIF
                    AMN  = MIN( AMN, A(I,J) )
                    AMX  = MAX( AMX, A(I,J) )
                  ENDIF
                ENDDO
              ENDDO
            ELSE !ave.eq.0
              Q = 1.0 ! no scaling
            ENDIF
          ELSE !.not.LBAL
            AMN  =  SPVAL
            AMX  = -SPVAL
            AVE  = SUMA/SUMT - SMEAN
            DO J= 1,JDM
              DO I= 1,IDM
                IF     (   A(I,J).NE.SPVAL .AND.
     +                  AMSK(I,J).NE.SPVAL      ) THEN
                  A(I,J) = A(I,J) - AVE
                  AMN  = MIN( AMN, A(I,J) )
                  AMX  = MAX( AMX, A(I,J) )
                ENDIF
              ENDDO
            ENDDO
          ENDIF  !lbal:else
        ELSE
          Q   = 1.0
          AVE = 0.0
        ENDIF
        IF     (NPAD.EQ.0) THEN
          WRITE(21,REC=K,IOSTAT=IOS) A
        ELSE
          WRITE(21,REC=K,IOSTAT=IOS) A,PAD
        ENDIF
        IF     (LBAL) THEN
          WRITE(6,'(a,1p4g16.6)')
     &      'min,mn,scl,max=',AMN,SMEAN,Q,AMX
        ELSE
          WRITE(6,'(a,1p4g16.6)')
     &      'min,mn,off,max=',AMN,SMEAN,-AVE,AMX
        ENDIF
      ENDDO !k
C
      IF     (SUMTR(1).EQ.SUMTR(2)) THEN
        WRITE(6,*)
        WRITE(6,*) K-1,' FIELDS PROCESSED (AREA = ',SUMT,' M^2)'
        WRITE(6,*)
      ELSE
        WRITE(6,*)
        WRITE(6,*) K-1,' FIELDS PROCESSED (AREA = ',SUMTR(1),
     &                 ' TO ',SUMTR(2),' M^2)'
        WRITE(6,*)
      ENDIF
C
      CLOSE(UNIT=11)
      CLOSE(UNIT=21)
C
      RETURN
      END

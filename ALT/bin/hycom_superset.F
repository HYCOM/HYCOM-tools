      PROGRAM HYCOM_SUPERSET
      IMPLICIT NONE
C
C  hycom_superset - Usage:  hycom_superset fin.a [fmask.a] idm jdm i1 j1 idms jdms [spval] fout.a
C  hycom_mergeset - Usage:  hycom_superset fin.a [fmask.a] idm jdm i1 j1 idms jdms [spval] fout.a
C
C                 Outputs a superset of the input array, or
C                 merges the input array into the larger output array.
C
C  fin*.a is assumed to contain idm*jdm 32-bit IEEE real values for
C   each array, in standard f77 element order, followed by padding
C   to a multiple of 4096 32-bit words, but otherwise with no control
C   bytes/words, and input values of 2.0**100 indicating a data void.
C
C  Use hycom_superset_raw or hycom_mergeset_raw if fin.a is a raw file
C   with no control words and no padding.  Note that fmask.a and fout.a
C   are still in HYCOM format.  In these cases only, data voids in fin.a
C   are optionally indicated by the input spval (default 2.0**100).
C   If spval is present, then fmask.a must also be present.
C
C  If (i1:i1+idm-1:j1:j1+jdm-1) isn't inside (1:idms,1:jdms), the output
C  fields are assumed to be p-grid global with an arctic bi-polar patch.
C
C  For hycom_superset, fout.a must not exist and the output field will 
C   be data voids outside the specified input region, and a copy of 
C   fin.a inside the input region.
C   If fmask.a is provided fout.a will additionally have data voids 
C   in the specified input region at all locations that are voids 
C   in the 1st field of fmask.a.
C  For hycom_mergeset, fout.a must exist and will already contain
C   the result outside the specified input region, and a copy of 
C   fin.a inside the input region.
C   If fmask.a is provided fout.a will additionally be unchanged
C   in the specified input region at all locations that are voids 
C   in the 1st field of fmask.a.
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  January 2001.
C
      REAL*4, ALLOCATABLE :: A(:,:),AM(:,:),A2(:,:)
      REAL*4              :: PAD(4096),PAD2(4096)
      INTEGER       IOS,IOS2
      INTEGER       IARGC
      INTEGER       NARG
      CHARACTER*240 CARG
C
      LOGICAL       LDEBUG,LMERGE,LRAW,LSPVAL
      REAL*4        SPVAL
      INTEGER       IDM,JDM,I1,J1,IDMS,JDMS,L,NPAD,NPAD2
      CHARACTER*240 CFILE,CFILEM,CFILEO
C
C     READ ARGUMENTS.
C
      CALL GETARG(0,CARG)
      L = LEN_TRIM(CARG)
*     WRITE(6,"(4a)") TRIM(CARG),'"',CARG(L-8:L),'"'
      IF     (CARG(L-5:L).EQ.'_debug') THEN  !undocumented debug mode
        LDEBUG = .TRUE.
        L = L - 6
      ELSE
        LDEBUG = .FALSE.
      ENDIF
      IF     (CARG(L-8:L).EQ.'_superset') THEN
        LMERGE = .FALSE.
        LRAW   = .FALSE.
      ELSEIF (CARG(L-8:L).EQ.'_mergeset') THEN
        LMERGE = .TRUE.
        LRAW   = .FALSE.
      ELSEIF (CARG(L-12:L).EQ.'_superset_raw') THEN
        LMERGE = .FALSE.
        LRAW   = .TRUE.
      ELSEIF (CARG(L-12:L).EQ.'_mergeset_raw') THEN
        LMERGE = .TRUE.
        LRAW   = .TRUE.
      ELSE
        WRITE(6,'(2a)')
     &    'Usage:  ',
     &    'hycom_superset[_raw] or hycom_mergeset[_raw]'
        write(6,'(a)') CARG(L-8:L)
        write(6,'(a)') CARG(L-12:L)
        CALL EXIT(1)
      ENDIF
C
      NARG = IARGC()
C
      IF     (NARG.EQ.8) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        CFILEM = ' '
        READ(CARG,*) IDM
        CALL GETARG(3,CARG)
        READ(CARG,*) JDM
        CALL GETARG(4,CARG)
        READ(CARG,*) I1
        CALL GETARG(5,CARG)
        READ(CARG,*) J1
        CALL GETARG(6,CARG)
        READ(CARG,*) IDMS
        CALL GETARG(7,CARG)
        READ(CARG,*) JDMS
        CALL GETARG(8,CFILEO)
        LSPVAL = .FALSE.
      ELSEIF (NARG.EQ.9) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CFILEM)
        CALL GETARG(3,CARG)
        READ(CARG,*) IDM
        CALL GETARG(4,CARG)
        READ(CARG,*) JDM
        CALL GETARG(5,CARG)
        READ(CARG,*) I1
        CALL GETARG(6,CARG)
        READ(CARG,*) J1
        CALL GETARG(7,CARG)
        READ(CARG,*) IDMS
        CALL GETARG(8,CARG)
        READ(CARG,*) JDMS
        CALL GETARG(9,CFILEO)
        LSPVAL = .FALSE.
      ELSEIF (NARG.EQ.10 .AND. LRAW) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CFILEM)
        CALL GETARG(3,CARG)
        READ(CARG,*) IDM
        CALL GETARG(4,CARG)
        READ(CARG,*) JDM
        CALL GETARG(5,CARG)
        READ(CARG,*) I1
        CALL GETARG(6,CARG)
        READ(CARG,*) J1
        CALL GETARG(7,CARG)
        READ(CARG,*) IDMS
        CALL GETARG(8,CARG)
        READ(CARG,*) JDMS
        CALL GETARG(9,CARG)
        READ(CARG,*) SPVAL
        LSPVAL = .TRUE.
        CALL GETARG(10,CFILEO)
      ELSEIF (LRAW) THEN
        WRITE(6,"(2a)")
     &    'Usage: hycom_{super,merge}set_raw fin.a [fmask.a] idm jdm',
     &    ' i1 j1 idms jdms [spval] fout.a'
        CALL EXIT(1)
      ELSE
        WRITE(6,"(2a)")
     &    'Usage: hycom_{super,merge}set fin.a [fmask.a] idm jdm',
     &    ' i1 j1 idms jdms fout.a'
        CALL EXIT(1)
      ENDIF
C
      WRITE(6,'(a,2i6)')
     &     'I1,J1 = ',I1,      J1
      WRITE(6,'(a,2i6)')
     &     'IL,JL = ',I1-1+IDM,J1-1+JDM
C
      NPAD = 4096 - MOD(IDM*JDM,4096)
      IF     (NPAD.EQ.4096) THEN
        NPAD = 0
      ENDIF
      NPAD2 = 4096 - MOD(IDMS*JDMS,4096)
      IF     (NPAD2.EQ.4096) THEN
        NPAD2 = 0
      ENDIF
C
      ALLOCATE( A(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_superset: could not allocate 1st ',
     +             IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( AM(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_superset: could not allocate 1st ',
     +             IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( A2(IDMS,JDMS), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_superset: could not allocate last ',
     +             IDMS*JDMS,' words'
        CALL EXIT(2)
      ENDIF
C
      CALL SUPERSET(A,AM,IDM,JDM,A2,IDMS,JDMS,I1,J1,
     &              SPVAL,LSPVAL, LMERGE,LRAW,LDEBUG,
     &              PAD,NPAD,PAD2,NPAD2, CFILE,CFILEM,CFILEO)
      CALL EXIT(0)
      END
      SUBROUTINE SUPERSET(A,AM,IDM,JDM,A2,IDMS,JDMS,I1,J1,
     &                    SPVAL,LSPVAL, LMERGE,LRAW,LDEBUG,
     &                    PAD,NPAD,PAD2,NPAD2, CFILE,CFILEM,CFILEO)
      IMPLICIT NONE
C
      REAL*4     SPVALH
      PARAMETER (SPVALH=2.0**100)
C
      CHARACTER*240 CFILE,CFILEM,CFILEO
      LOGICAL      LSPVAL,LMERGE,LRAW,LDEBUG
      INTEGER      IDM,JDM,NPAD,IDMS,JDMS,NPAD2,I1,J1
      REAL*4       SPVAL
      REAL*4       A(IDM,JDM),AM(IDM,JDM),PAD(NPAD)
      REAL*4       A2(IDMS,JDMS),PAD2(NPAD2)
C
C     MOST OF WORK IS DONE HERE.
C
#ifdef CRAY
      INTEGER*8    IU8,IOS8
#endif
      CHARACTER*18 CASN
      INTEGER      I,II,J,JDMA,JJ,K,IOS,NRECL,NRECLM,NRECL2
      REAL*4       AMN,AMX
C
      IF     (.NOT.LSPVAL) THEN
        SPVAL = SPVALH
      ENDIF
C
      IF     (NPAD.EQ.0 .OR. LRAW) THEN
        INQUIRE(IOLENGTH=NRECL)  A
      ELSE
        INQUIRE(IOLENGTH=NRECL)  A,PAD
      ENDIF
      IF     (NPAD.EQ.0) THEN
        INQUIRE(IOLENGTH=NRECLM) AM
      ELSE
        INQUIRE(IOLENGTH=NRECLM) AM,PAD
      ENDIF
      IF     (NPAD2.EQ.0) THEN
        INQUIRE(IOLENGTH=NRECL2) A2
      ELSE
        INQUIRE(IOLENGTH=NRECL2) A2,PAD2
        PAD2(:) = SPVALH
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(PAD2,NPAD2)
#endif
      ENDIF
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
      ENDIF
      IF     (CFILEM.NE." " .AND. MOD(NRECLM,4096).EQ.0) THEN
        WRITE(CASN,8000) NRECLM/4096
 8000   FORMAT('-F cachea:',I4.4,':1:0')
        IU8 = 12
        CALL ASNUNIT(IU8,CASN,IOS8)
        IF     (IOS8.NE.0) THEN
          write(6,*) 'Error: can''t asnunit 12'
          write(6,*) 'ios  = ',ios8
          write(6,*) 'casn = ',casn
          CALL EXIT(5)
        ENDIF
      ENDIF
      IF     (MOD(NRECL2,4096).EQ.0) THEN
        WRITE(CASN,8000) NRECL2/4096
 8000   FORMAT('-F cachea:',I4.4,':1:0')
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
      IF     (CFILEM.NE." ") THEN
        CALL ASNUNIT(12,'-F syscall -N ieee',IOS)
        IF     (IOS.NE.0) THEN
          write(6,*) 'Error: can''t asnunit 12'
          write(6,*) 'ios = ',ios
          CALL EXIT(5)
        ENDIF
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
        write(6,*) 'Error: can''t open ',TRIM(CFILE)
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
      IF     (.NOT.LMERGE) THEN
        OPEN(UNIT=21, FILE=CFILEO, FORM='UNFORMATTED', STATUS='NEW',
     +           ACCESS='DIRECT', RECL=NRECL2, IOSTAT=IOS)
      ELSE
        OPEN(UNIT=21, FILE=CFILEO, FORM='UNFORMATTED', STATUS='OLD',
     +           ACCESS='DIRECT', RECL=NRECL2, IOSTAT=IOS)
      ENDIF
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILEO)
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
C
      IF     (CFILEM.EQ." ") THEN
        AM(:,:) = 0.0  !no masking
      ELSE
        OPEN(UNIT=12, FILE=CFILEM, FORM='UNFORMATTED', STATUS='OLD',
     +           ACCESS='DIRECT', RECL=NRECLM, IOSTAT=IOS)
        IF     (IOS.NE.0) THEN
          write(6,*) 'Error: can''t open ',TRIM(CFILE)
          write(6,*) 'ios   = ',ios
          write(6,*) 'nrecl = ',nrecl
          CALL EXIT(3)
        ENDIF
        READ(12,REC=1,IOSTAT=IOS) AM
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(AM,IDM*JDM)
#endif
        IF     (IOS.NE.0) THEN
          WRITE(6,*) 'can''t read ',TRIM(CFILEM)
          CALL EXIT(4)
        ENDIF
        CLOSE(12)
      ENDIF
C
      JDMA = JDMS - J1 + 1   !location of top of target array
C
      DO 110 K= 1,9999
        READ(11,REC=K,IOSTAT=IOS) A
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
        IF     (IOS.NE.0) THEN
          IF     (K.EQ.1) THEN
            WRITE(6,*) 'can''t read ',TRIM(CFILE)
            CALL EXIT(4)
          ELSE
            GOTO 1110
          ENDIF
        ENDIF
C
        IF     (LMERGE) THEN
          READ(21,REC=K,IOSTAT=IOS) A2
#ifdef ENDIAN_IO
          CALL ENDIAN_SWAP(A2,IDM*JDM)
#endif
          IF     (IOS.NE.0) THEN
            WRITE(6,*) 'can''t read ',TRIM(CFILEO),
     &                 ' at record ',K
            CALL EXIT(4)
          ENDIF
        ELSE
          DO J= 1,JDMS
            DO I= 1,IDMS
              A2(I,J) = SPVALH
            ENDDO
          ENDDO
        ENDIF !lmerge:else
C
        IF     (LDEBUG) THEN
          DO JJ= 1,MIN(JDM,JDMA)
            J = J1+JJ-1
            DO II= 1,IDM
              IF     ( A(II,JJ).NE.SPVAL  .AND.
     &                AM(II,JJ).NE.SPVALH      ) THEN
                I = MOD(I1+II-2+9*IDMS,IDMS) + 1  !assumed periodic
*               IF     (A2(I,J).NE.SPVALH) THEN
*               IF     (A2(I,J).NE.A(II,JJ)) THEN
                IF     (J.GE.JDMS-1) THEN
                  WRITE(6,'(a,2i6,1p2g16.6,2i6)')
     &              'overwrite, a2(i,j) with a(ii,jj):',
     &              I,J,A2(I,J),A(II,JJ),II,JJ
                ENDIF
                A2(I,J) = A(II,JJ)
              ENDIF !AM
            ENDDO !ii
          ENDDO !jj
          IF     (JDMA.LT.JDM) THEN
            DO JJ= JDMA-1,JDM
              J = J1+JJ-1
              J = 2*JDMS-J-1  !assumed arctic patch
              DO II= 1,IDM
                IF     ( A(II,JJ).NE.SPVAL  .AND.
     &                  AM(II,JJ).NE.SPVALH      ) THEN
                  I = MOD(I1+II-2+9*IDMS,IDMS) + 1  !assumed periodic
                  I = IDMS - MOD(I-1,IDMS)          !assumed arctic patch
*                 IF     (A2(I,J).NE.SPVALH) THEN
*                 IF     (A2(I,J).NE.A(II,JJ)) THEN
                  IF     (J.GE.JDMS-1) THEN
                    WRITE(6,'(a,2i6,1p2g16.6,2i6)')
     &                'overwrite, a2(i,j) WITH a(ii,jj):',
     &                I,J,A2(I,J),A(II,JJ),II,JJ
                  ENDIF
                  A2(I,J) = A(II,JJ)
                ENDIF !AM
              ENDDO !ii
            ENDDO !jj
          ENDIF !jdma
        ELSE
          DO JJ= 1,MIN(JDM,JDMA)
            J = J1+JJ-1
            DO II= 1,IDM
              IF     ( A(II,JJ).NE.SPVAL  .AND.
     &                AM(II,JJ).NE.SPVALH      ) THEN
                I = MOD(I1+II-2+9*IDMS,IDMS) + 1  !assumed periodic
                A2(I,J) = A(II,JJ)
              ENDIF !AM
            ENDDO !ii
          ENDDO !jj
          IF     (JDMA.LT.JDM) THEN
            DO JJ= JDMA-1,JDM
              J = J1+JJ-1
              J = 2*JDMS-J-1  !assumed arctic patch
              DO II= 1,IDM
                IF     ( A(II,JJ).NE.SPVAL  .AND.
     &                  AM(II,JJ).NE.SPVALH      ) THEN
                  I = MOD(I1+II-2+9*IDMS,IDMS) + 1  !assumed periodic
                  I = IDMS - MOD(I-1,IDMS)          !assumed arctic patch
                  A2(I,J) = A(II,JJ)
                ENDIF !AM
              ENDDO !ii
            ENDDO !jj
          ENDIF !jdma
        ENDIF
C
        AMN =  SPVALH
        AMX = -SPVALH
        DO J= 1,JDMS
          DO I= 1,IDMS
            IF     (A2(I,J).NE.SPVALH) THEN
              AMX = MAX( AMX, A2(I,J) )
              AMN = MIN( AMN, A2(I,J) )
            ENDIF
          ENDDO
        ENDDO
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(A2,IDMS*JDMS)
#endif
        IF     (NPAD2.EQ.0) THEN
          WRITE(21,REC=K,IOSTAT=IOS) A2
        ELSE
          WRITE(21,REC=K,IOSTAT=IOS) A2,PAD2
        ENDIF
        WRITE(6,'(a,1p2g16.8)')
     &     'min, max = ',AMN,AMX
  110 CONTINUE
 1110 CONTINUE
      WRITE(6,*) 
      WRITE(6,*) K-1,' FIELDS PROCESSED (IDMS,JDMS = ',IDMS,JDMS,')'
      WRITE(6,*) 
      CLOSE(21)
      RETURN
      END

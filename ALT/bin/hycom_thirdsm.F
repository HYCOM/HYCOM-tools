      PROGRAM PTHIRDSM
      IMPLICIT NONE
C
C  hycom_thirdbox  - Usage:  hycom_thirdbox  fin.a idm jdm [i1 j1 idm2 jdm2] fout.a [narctic]
C  hycom_thirdboxa - Usage:  hycom_thirdboxa fin.a idm jdm [i1 j1 idm2 jdm2] fout.a [narctic] grid.a
C  hycom_thirdsm   - Usage:  hycom_thirdsm   fin.a idm jdm [i1 j1 idm2 jdm2] fout.a [ns [narctic]]
C  hycom_thirdsma  - Usage:  hycom_thirdsma  fin.a idm jdm [i1 j1 idm2 jdm2] fout.a [ns [narctic]] grid.a
C
C                 Outputs a 1/3x as fine grid using subsampling
C                  after averaging with:
C                  (hycom_thirdsm)  a 9-point smoothing operator, or
C                  (hycom_thirdbox) a 9-point boxcar (cell average)
C                 The "a" versions allow for p-grid cell area
C
C                 i1,  j1   = 1,1 by default
C                 idm2,jdm2 = (idm+2)/3,(jdm+2)/3 by default
C                 ns is the number of smoothing passes, 0 or 1 (default)
C                 ns=0 selects the cell center (every 3rd point)
C                 narctic selects the output type (default 0=standard)
C                   0=standard; 1=p-scalar arctic; 11=p-vector arctic;
C
C  fin*.a is assumed to contain idm*jdm 32-bit IEEE real values for
C   each array, in standard f77 element order, followed by padding
C   to a multiple of 4096 32-bit words, but otherwise with no control
C   bytes/words, and input values of 2.0**100 indicating a data void.
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  January 2001.
C
      REAL*4, ALLOCATABLE :: A(:,:),B(:,:),AREA(:,:),A2(:,:)
      REAL*4              :: PAD(4096),PAD2(4096)
      LOGICAL      LBOX,LAREA
      INTEGER      IOS,IOS2,L
      INTEGER      IARGC
      INTEGER      NARG
      CHARACTER*240 CARG
C
      INTEGER      IDM,JDM,I1,J1,IDM2,JDM2,NPAD,NPAD2,NS,NARCTIC
      CHARACTER*240 CFILE,CFILEO,CFILEG
C
C     sm or box
C
      CALL GETARG(0,CARG)
      L = LEN_TRIM(CARG)
*     WRITE(6,"(4a)") TRIM(CARG),'"',CARG(L-2:L),'"'
      IF     (CARG(L-2:L).EQ.'box') THEN
        LBOX  = .TRUE.
        LAREA = .FALSE.
*       write(6,*) 'hycom_thirdbox ...'
      ELSEIF (CARG(L-3:L).EQ.'boxa') THEN
        LBOX  = .TRUE.
        LAREA = .TRUE.
*       write(6,*) 'hycom_thirdboxa ...'
      ELSEIF (CARG(L-1:L).EQ.'sm') THEN
        LBOX  = .FALSE.
        LAREA = .FALSE.
*       write(6,*) 'hycom_thirdsm ...'
      ELSEIF (CARG(L-2:L).EQ.'sma') THEN
        LBOX  = .FALSE.
        LAREA = .TRUE.
*       write(6,*) 'hycom_thirdsma ...'
      ELSE
        WRITE(6,'(2a)')
     &    'Usage:  ',
     &  'hycom_thirdsm[a] or hycom_thirdbox[a] ...'
        CALL EXIT(1)
      ENDIF
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (LAREA) THEN
        IF     (NARG.EQ.5) THEN
          CALL GETARG(1,CFILE)
          CALL GETARG(2,CARG)
          READ(CARG,*) IDM
          CALL GETARG(3,CARG)
          READ(CARG,*) JDM
          I1   = 1
          J1   = 1
          IDM2 = (IDM+2)/3
          JDM2 = (JDM+2)/3
          CALL GETARG(4,CFILEO)
          NS      = 1
          NARCTIC = 0
          CALL GETARG(5,CFILEG)
        ELSEIF (NARG.EQ.6 .AND. .NOT.LBOX) THEN
          CALL GETARG(1,CFILE)
          CALL GETARG(2,CARG)
          READ(CARG,*) IDM
          CALL GETARG(3,CARG)
          READ(CARG,*) JDM
          I1   = 1
          J1   = 1
          IDM2 = (IDM+2)/3
          JDM2 = (JDM+2)/3
          CALL GETARG(4,CFILEO)
          CALL GETARG(5,CARG)
          READ(CARG,*) NS
          NARCTIC = 0
          CALL GETARG(6,CFILEG)
        ELSEIF (NARG.EQ.6 .AND. LBOX) THEN
          CALL GETARG(1,CFILE)
          CALL GETARG(2,CARG)
          READ(CARG,*) IDM
          CALL GETARG(3,CARG)
          READ(CARG,*) JDM
          I1   = 1
          J1   = 1
          IDM2 = (IDM+2)/3
          JDM2 = (JDM+2)/3
          CALL GETARG(4,CFILEO)
          CALL GETARG(5,CARG)
          READ(CARG,*) NARCTIC
          NS = 0
          CALL GETARG(6,CFILEG)
        ELSEIF (NARG.EQ.9) THEN
          CALL GETARG(1,CFILE)
          CALL GETARG(2,CARG)
          READ(CARG,*) IDM
          CALL GETARG(3,CARG)
          READ(CARG,*) JDM
          CALL GETARG(4,CARG)
          READ(CARG,*) I1
          CALL GETARG(5,CARG)
          READ(CARG,*) J1
          CALL GETARG(6,CARG)
          READ(CARG,*) IDM2
          CALL GETARG(7,CARG)
          READ(CARG,*) JDM2
          CALL GETARG(8,CFILEO)
          NS      = 1
          NARCTIC = 0
          CALL GETARG(9,CFILEG)
        ELSEIF (NARG.EQ.10 .AND. .NOT.LBOX) THEN
          CALL GETARG(1,CFILE)
          CALL GETARG(2,CARG)
          READ(CARG,*) IDM
          CALL GETARG(3,CARG)
          READ(CARG,*) JDM
          CALL GETARG(4,CARG)
          READ(CARG,*) I1
          CALL GETARG(5,CARG)
          READ(CARG,*) J1
          CALL GETARG(6,CARG)
          READ(CARG,*) IDM2
          CALL GETARG(7,CARG)
          READ(CARG,*) JDM2
          CALL GETARG(8,CFILEO)
          CALL GETARG(9,CARG)
          READ(CARG,*) NS
          NARCTIC = 0
          CALL GETARG(10,CFILEG)
        ELSEIF (NARG.EQ.10 .AND. LBOX) THEN
          CALL GETARG(1,CFILE)
          CALL GETARG(2,CARG)
          READ(CARG,*) IDM
          CALL GETARG(3,CARG)
          READ(CARG,*) JDM
          CALL GETARG(4,CARG)
          READ(CARG,*) I1
          CALL GETARG(5,CARG)
          READ(CARG,*) J1
          CALL GETARG(6,CARG)
          READ(CARG,*) IDM2
          CALL GETARG(7,CARG)
          READ(CARG,*) JDM2
          CALL GETARG(8,CFILEO)
          CALL GETARG(9,CARG)
          READ(CARG,*) NARCTIC
          NS = 0
          CALL GETARG(10,CFILEG)
        ELSEIF (NARG.EQ.11 .AND. .NOT.LBOX) THEN
          CALL GETARG(1,CFILE)
          CALL GETARG(2,CARG)
          READ(CARG,*) IDM
          CALL GETARG(3,CARG)
          READ(CARG,*) JDM
          CALL GETARG(4,CARG)
          READ(CARG,*) I1
          CALL GETARG(5,CARG)
          READ(CARG,*) J1
          CALL GETARG(6,CARG)
          READ(CARG,*) IDM2
          CALL GETARG(7,CARG)
          READ(CARG,*) JDM2
          CALL GETARG(8,CFILEO)
          CALL GETARG(9,CARG)
          READ(CARG,*) NS
          CALL GETARG(10,CARG)
          READ(CARG,*) NARCTIC
          CALL GETARG(11,CFILEG)
        ELSEIF (LBOX) THEN
          WRITE(6,*)
     &      'Usage: hycom_thirdboxa fin.a idm jdm [i1 j1 idm2 jdm2]'//
     &     ' fout.a [narctic] grid.a'
          CALL EXIT(1)
        ELSE
          WRITE(6,*)
     &      'Usage: hycom_thirdsma fin.a idm jdm [i1 j1 idm2 jdm2]'//
     &      ' fout.a [ns [narctic]] grid.a'
          CALL EXIT(1)
        ENDIF !narg
      ELSE !original versions
        IF     (NARG.EQ.4) THEN
          CALL GETARG(1,CFILE)
          CALL GETARG(2,CARG)
          READ(CARG,*) IDM
          CALL GETARG(3,CARG)
          READ(CARG,*) JDM
          I1   = 1
          J1   = 1
          IDM2 = (IDM+2)/3
          JDM2 = (JDM+2)/3
          CALL GETARG(4,CFILEO)
          NS      = 1
          NARCTIC = 0
        ELSEIF (NARG.EQ.5 .AND. .NOT.LBOX) THEN
          CALL GETARG(1,CFILE)
          CALL GETARG(2,CARG)
          READ(CARG,*) IDM
          CALL GETARG(3,CARG)
          READ(CARG,*) JDM
          I1   = 1
          J1   = 1
          IDM2 = (IDM+2)/3
          JDM2 = (JDM+2)/3
          CALL GETARG(4,CFILEO)
          CALL GETARG(5,CARG)
          READ(CARG,*) NS
          NARCTIC = 0
        ELSEIF (NARG.EQ.5 .AND. LBOX) THEN
          CALL GETARG(1,CFILE)
          CALL GETARG(2,CARG)
          READ(CARG,*) IDM
          CALL GETARG(3,CARG)
          READ(CARG,*) JDM
          I1   = 1
          J1   = 1
          IDM2 = (IDM+2)/3
          JDM2 = (JDM+2)/3
          CALL GETARG(4,CFILEO)
          CALL GETARG(5,CARG)
          READ(CARG,*) NARCTIC
          NS = 0
        ELSEIF (NARG.EQ.8) THEN
          CALL GETARG(1,CFILE)
          CALL GETARG(2,CARG)
          READ(CARG,*) IDM
          CALL GETARG(3,CARG)
          READ(CARG,*) JDM
          CALL GETARG(4,CARG)
          READ(CARG,*) I1
          CALL GETARG(5,CARG)
          READ(CARG,*) J1
          CALL GETARG(6,CARG)
          READ(CARG,*) IDM2
          CALL GETARG(7,CARG)
          READ(CARG,*) JDM2
          CALL GETARG(8,CFILEO)
          NS      = 1
          NARCTIC = 0
        ELSEIF (NARG.EQ.9 .AND. .NOT.LBOX) THEN
          CALL GETARG(1,CFILE)
          CALL GETARG(2,CARG)
          READ(CARG,*) IDM
          CALL GETARG(3,CARG)
          READ(CARG,*) JDM
          CALL GETARG(4,CARG)
          READ(CARG,*) I1
          CALL GETARG(5,CARG)
          READ(CARG,*) J1
          CALL GETARG(6,CARG)
          READ(CARG,*) IDM2
          CALL GETARG(7,CARG)
          READ(CARG,*) JDM2
          CALL GETARG(8,CFILEO)
          CALL GETARG(9,CARG)
          READ(CARG,*) NS
          NARCTIC = 0
        ELSEIF (NARG.EQ.9 .AND. LBOX) THEN
          CALL GETARG(1,CFILE)
          CALL GETARG(2,CARG)
          READ(CARG,*) IDM
          CALL GETARG(3,CARG)
          READ(CARG,*) JDM
          CALL GETARG(4,CARG)
          READ(CARG,*) I1
          CALL GETARG(5,CARG)
          READ(CARG,*) J1
          CALL GETARG(6,CARG)
          READ(CARG,*) IDM2
          CALL GETARG(7,CARG)
          READ(CARG,*) JDM2
          CALL GETARG(8,CFILEO)
          CALL GETARG(9,CARG)
          READ(CARG,*) NARCTIC
          NS = 0
        ELSEIF (NARG.EQ.10 .AND. .NOT.LBOX) THEN
          CALL GETARG(1,CFILE)
          CALL GETARG(2,CARG)
          READ(CARG,*) IDM
          CALL GETARG(3,CARG)
          READ(CARG,*) JDM
          CALL GETARG(4,CARG)
          READ(CARG,*) I1
          CALL GETARG(5,CARG)
          READ(CARG,*) J1
          CALL GETARG(6,CARG)
          READ(CARG,*) IDM2
          CALL GETARG(7,CARG)
          READ(CARG,*) JDM2
          CALL GETARG(8,CFILEO)
          CALL GETARG(9,CARG)
          READ(CARG,*) NS
          CALL GETARG(10,CARG)
          READ(CARG,*) NARCTIC
        ELSEIF (LBOX) THEN
          WRITE(6,*)
     &      'Usage: hycom_thirdbox fin.a idm jdm [i1 j1 idm2 jdm2]'//
     &     ' fout.a [narctic]'
          CALL EXIT(1)
        ELSE
          WRITE(6,*)
     &      'Usage: hycom_thirdsm fin.a idm jdm [i1 j1 idm2 jdm2]'//
     &      ' fout.a [ns [narctic]]'
          CALL EXIT(1)
        ENDIF !narg
      ENDIF  !larea:else
C
      IF     (I1+3*IDM2-3.GT.IDM+2 .OR.
     &        J1+3*JDM2-3.GT.JDM+2     ) THEN
        WRITE(6,*) 'Error in hycom_thirdsm: idm2 or jdm2 too large'
        CALL EXIT(1)
      ENDIF
C
      IF     (NS.LT.0 .OR. NS.GT.1) THEN
        WRITE(6,*) 'Error in hycom_thirdsm: ns must be 0 or 1'
        CALL EXIT(1)
      ENDIF
C
      IF     (NARCTIC.NE.0  .AND.
     &        NARCTIC.NE.1  .AND.
     &        NARCTIC.NE.11      ) THEN
        WRITE(6,*) 'Error in hycom_thirdsm: narctic must be 0,1,11'
        CALL EXIT(1)
      ENDIF
C
      NPAD = 4096 - MOD(IDM*JDM,4096)
      IF     (NPAD.EQ.4096) THEN
        NPAD = 0
      ENDIF
      NPAD2 = 4096 - MOD(IDM2*JDM2,4096)
      IF     (NPAD2.EQ.4096) THEN
        NPAD2 = 0
      ENDIF
C
      ALLOCATE( A(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_thirdsm: could not allocate ',
     +             IDM*JDM,' words for A'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( AREA(0:IDM+3,0:JDM+3), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_thirdsm: could not allocate ',
     +             (IDM+4)*(JDM+4),' words for AREA'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( B(0:IDM+3,0:JDM+3), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_thirdsm: could not allocate  ',
     +             (IDM+4)*(JDM+4),' words for B'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( A2(IDM2,JDM2), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_thirdsm: could not allocate ',
     +             IDM2*JDM2,' words for A2'
        CALL EXIT(2)
      ENDIF
C
      CALL THIRDSM(A,AREA,B,IDM,JDM,A2,IDM2,JDM2,I1,J1,
     &             PAD,NPAD,PAD2,NPAD2, LBOX,LAREA,NS,NARCTIC,
     &             CFILE,CFILEO,CFILEG)
      CALL EXIT(0)
      END
      SUBROUTINE THIRDSM(A,AREA,B,IDM,JDM,A2,IDM2,JDM2,I1,J1,
     &                   PAD,NPAD,PAD2,NPAD2, LBOX,LAREA,NS,NARCTIC,
     &                   CFILE,CFILEO,CFILEG)
      IMPLICIT NONE
C
      REAL*4     SPVAL
      PARAMETER (SPVAL=2.0**100)
C
      CHARACTER*240 CFILE,CFILEO,CFILEG
      LOGICAL      LBOX,LAREA
      INTEGER      IDM,JDM,NPAD,IDM2,JDM2,NPAD2,NS,NARCTIC,I1,J1
      REAL*4       A(IDM,JDM),PAD(NPAD)
      REAL*4       B(0:IDM+3,0:JDM+3),AREA(0:IDM+3,0:JDM+3)
      REAL*4       A2(IDM2,JDM2),PAD2(NPAD2)
C
C     MOST OF WORK IS DONE HERE.
C
#ifdef CRAY
      INTEGER*8    IU8,IOS8
#endif
      CHARACTER*18 CASN
      LOGICAL      LPERIOD,LARCTIC
      INTEGER      LEN_TRIM
      INTEGER      I,II,IQ,J,JJ,JQ,K,IOS,NRECL,NRECL2
      REAL*4       AMN,AMX
      REAL*8       RC,RS
C
      REAL*8       C(-1:1,-1:1)
      SAVE         C
      DATA         C /  1.0, 2.0, 1.0,
     +                  2.0, 4.0, 2.0,
     +                  1.0, 2.0, 1.0 /
C
      IF     (LBOX) THEN
        C(:,:) = 1.0
      ENDIF
C
      INQUIRE( IOLENGTH=NRECL)  A, PAD
      INQUIRE( IOLENGTH=NRECL2) A2,PAD2
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
      CALL ASNUNIT(21,'-F syscall -N ieee',IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t asnunit 21'
        write(6,*) 'ios = ',ios
        CALL EXIT(5)
      ENDIF
#endif
#endif
      IF     (LAREA) THEN
        OPEN(UNIT=11, FILE=CFILEG, FORM='UNFORMATTED', STATUS='OLD',
     +           ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
        IF     (IOS.NE.0) THEN
          write(6,*) 'Error: can''t open ',TRIM(CFILEG)
          write(6,*) 'ios   = ',ios
          write(6,*) 'nrecl = ',nrecl
          CALL EXIT(3)
        ENDIF
        READ(11,REC=10,IOSTAT=IOS) A  !pscx
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
        IF     (IOS.NE.0) THEN
          WRITE(6,*) 'can''t read ',TRIM(CFILEG)
          CALL EXIT(4)
        ENDIF
        AREA(1:IDM,1:JDM) = A(:,:)
        READ(11,REC=11,IOSTAT=IOS) A  !pscy
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
        IF     (IOS.NE.0) THEN
          WRITE(6,*) 'can''t read ',TRIM(CFILEG)
          CALL EXIT(4)
        ENDIF
        CLOSE(11)
C
        AREA(1:IDM,1:JDM) = AREA(1:IDM,1:JDM)*A(:,:)  !pscx*pscy
      ELSE
        AREA(:,:) = 1.0
      ENDIF !larea
C
      OPEN(UNIT=11, FILE=CFILE, FORM='UNFORMATTED', STATUS='OLD',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',CFILE(1:LEN_TRIM(CFILE))
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
      OPEN(UNIT=21, FILE=CFILEO, FORM='UNFORMATTED', STATUS='NEW',
     +         ACCESS='DIRECT', RECL=NRECL2, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',CFILEO(1:LEN_TRIM(CFILEO))
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
C
      DO 110 K= 1,9999
        READ(11,REC=K,IOSTAT=IOS) A
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
        IF     (IOS.NE.0) THEN
          IF     (K.EQ.1) THEN
            WRITE(6,*) 'can''t read ',CFILE(1:LEN_TRIM(CFILE))
            CALL EXIT(4)
          ELSE
            GOTO 1110
          ENDIF
        ENDIF
C
        LPERIOD = MINVAL(A(IDM,:)).NE.SPVAL
        LARCTIC = MINVAL(A(:,JDM)).NE.SPVAL
*       write(6,*) 'larctic = ',larctic,narctic
*       write(6,*) 'lperiod = ',lperiod
C
        IF     (K.EQ.1 .AND. LARCTIC) THEN
          DO I= 1,IDM
            AREA(I,0) = 1.0
          ENDDO
          DO J= JDM+1,JDM+3
            JJ = JDM-1-(J-JDM)
            DO I= 1,IDM
              II = IDM-MOD(I-1,IDM)
              AREA(I,J) = AREA(II,JJ)
            ENDDO !i
          ENDDO !j
        ELSE
          DO I= 1,IDM
            AREA(I,    0) = 1.0
            AREA(I,JDM+1) = 1.0
            AREA(I,JDM+2) = 1.0
            AREA(I,JDM+3) = 1.0
          ENDDO
        ENDIF
        IF     (LPERIOD) THEN
          DO J= 0,JDM+3
            AREA(    0,J) = AREA(IDM,J)
            AREA(IDM+1,J) = AREA(  1,J)
            AREA(IDM+2,J) = AREA(  2,J)
            AREA(IDM+3,J) = AREA(  3,J)
          ENDDO
        ELSE
          DO J= 0,JDM+3
            AREA(    0,J) = 1.0
            AREA(IDM+1,J) = 1.0
            AREA(IDM+2,J) = 1.0
            AREA(IDM+3,J) = 1.0
          ENDDO
        ENDIF !larctic, area
C
        DO J= 1,JDM
          DO I= 1,IDM
            B(I,J) = A(I,J)
          ENDDO
        ENDDO
        IF     (LARCTIC) THEN
          DO I= 1,IDM
            B(I,0) = SPVAL
          ENDDO
          DO J= JDM+1,JDM+3
            JJ = JDM-1-(J-JDM)
            DO I= 1,IDM
              II = IDM-MOD(I-1,IDM)
              B(I,J) = B(II,JJ)
            ENDDO !i
          ENDDO !j
        ELSE
          DO I= 1,IDM
            B(I,    0) = SPVAL
            B(I,JDM+1) = SPVAL
            B(I,JDM+2) = SPVAL
            B(I,JDM+3) = SPVAL
          ENDDO
        ENDIF
        IF     (LPERIOD) THEN
          DO J= 0,JDM+3
            B(    0,J) = B(IDM,J)
            B(IDM+1,J) = B(  1,J)
            B(IDM+2,J) = B(  2,J)
            B(IDM+3,J) = B(  3,J)
          ENDDO
        ELSE
          DO J= 0,JDM+3
            B(    0,J) = SPVAL
            B(IDM+1,J) = SPVAL
            B(IDM+2,J) = SPVAL
            B(IDM+3,J) = SPVAL
          ENDDO
        ENDIF
        AMN =  SPVAL
        AMX = -SPVAL
        DO J= 1,JDM2
          JJ =  J1+3*J-3
          DO I= 1,IDM2
            II =  I1+3*I-3
            RS = 0.0
            RC = 0.0
            DO JQ= -NS,NS
              DO IQ= -NS,NS  !a trick: 0,0 or -1,1
                IF     (B(II+IQ,JJ+JQ).NE.SPVAL) THEN
                  RS = RS + AREA(II+IQ,JJ+JQ)*C(IQ,JQ)*B(II+IQ,JJ+JQ)
                  RC = RC + AREA(II+IQ,JJ+JQ)*C(IQ,JQ)
                ENDIF
              ENDDO
            ENDDO
*               IF     (I.EQ.15 .AND. J.EQ.3325) THEN
*                 write(6,*) ' i,ii = ',i,ii
*                 write(6,*) ' j,jj = ',j,jj
*                 write(6,*) 'rc,rs = ',rc,rs
*               ENDIF
            IF     (RC.NE.0.0) THEN
              A2(I,J) = RS/RC
              AMX = MAX( AMX, A2(I,J) )
              AMN = MIN( AMN, A2(I,J) )
            ELSE
              A2(I,J) = SPVAL
            ENDIF
          ENDDO
        ENDDO
        IF     (NARCTIC.EQ.1)  THEN  !p-scalar
          J  = JDM2
          JJ = JDM2-1
          DO I= 1,IDM2
            II = IDM2-MOD(I-1,IDM2)
            A2(I,J) = A2(II,JJ)
          ENDDO !i
        ELSEIF (NARCTIC.EQ.11) THEN  !p-vector
          J  = JDM2
          JJ = JDM2-1
          DO I= 1,IDM2
            II = IDM2-MOD(I-1,IDM2)
            IF     (A2(II,JJ).NE.SPVAL) THEN
              A2(I,J) = -A2(II,JJ)
            ELSE
              A2(I,J) = SPVAL
            ENDIF
          ENDDO !i
        ENDIF
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(A2,IDM2*JDM2)
#endif
        WRITE(21,REC=K,IOSTAT=IOS) A2
        WRITE(6,'(a,1p2g16.8)')
     &     'min, max = ',AMN,AMX
  110 CONTINUE
 1110 CONTINUE
      WRITE(6,*) 
      WRITE(6,*) K-1,' FIELDS PROCESSED (IDM2,JDM2 = ',IDM2,JDM2,')'
      WRITE(6,*) 
C
      CLOSE(11)
      CLOSE(21)
C
      RETURN
      END

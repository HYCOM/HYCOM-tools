      PROGRAM ISOLAY
      IMPLICIT NONE
C
C  hycom_1st_isopyc - Usage:  hycom_1st_isopyc archv.a [itype] il.a
C
C                 generates a 1st isopycnal surface from a HYCOM archive file.
C
C   archv.a is assumed to be an HYCOM archive data file, with companion
C   header file archv.b.  Both standard and mean archive files are allowed.
C
C   itype is the selection method (default 1)
C           = 0; mixed-layer based
C           = 1; target-density based
C
C   If itype is 0,  then the 1st isopycnal surface is the interface
C   below the mixed layer.
C   If itype is 1,  then the 1st isopycnal surface is the interface
C   above the shallowest near-isopycnal layer (based on a comparison
C   to the layer's target isopycnal density).
C
C   il.a and il.b will contain the resulting 1st isopycnal surface field.
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  April 2003.
C
      REAL*4     QONEM,SPVAL
      PARAMETER (QONEM=1.0/9806.0, SPVAL=2.0**100)
C
      REAL*4, ALLOCATABLE :: RK(:,:),DP(:,:),PK(:,:),PM(:,:),ZM(:,:)
      REAL*4              :: PAD(4096)
C
      INTEGER       IARGC
      INTEGER       NARG
      CHARACTER*240 CARG
C
      INTEGER       IDM,JDM,KDM,NSURF,NLAY,IEXPT,YRFLAG
      INTEGER       NPAD,ITYPE,ITEST,JTEST
      REAL          THBASE,SIGMA(99),TIME
      CHARACTER*30  CMTYPE
      CHARACTER*240 CFILEA,CFILEB,CFILEM
C
      CHARACTER*18  CASN
      INTEGER       I,J,K,KREC,KREC0,IOS,NRECL
      REAL          ZMIN,ZMAX
#ifdef CRAY
      INTEGER*8     IU8,IOS8
#endif
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.2) THEN
        CALL GETARG(1,CFILEA)
        CALL GETARG(2,CFILEM)
        ITYPE = 1
        ITEST = 0
        JTEST = 0
      ELSEIF (NARG.EQ.3) THEN
        CALL GETARG(1,CFILEA)
        CALL GETARG(2,CARG)
        READ(CARG,*) ITYPE
        CALL GETARG(3,CFILEM)
        ITEST = 0
        JTEST = 0
      ELSEIF (NARG.EQ.5) THEN  !undocumented, for debugging
        CALL GETARG(1,CFILEA)
        CALL GETARG(2,CARG)
        READ(CARG,*) ITYPE
        CALL GETARG(3,CFILEM)
        CALL GETARG(4,CARG)
        READ(CARG,*) ITEST
        CALL GETARG(5,CARG)
        READ(CARG,*) JTEST
      ELSE
        WRITE(6,*) 
     +    'Usage: hycom_1st_isopyc archv.a [itype] il.a'
        CALL EXIT(1)
      ENDIF
C
C     EXTRACT MODEL PARAMETERS FROM ".b" FILE.
C
      CFILEB = CFILEA(1:LEN_TRIM(CFILEA)-1) // 'b'
      CALL READ_B(CFILEB,
     +            IEXPT,YRFLAG,IDM,JDM,KDM,NSURF,NLAY,
     +            THBASE,SIGMA,TIME)
C
C     OPEN ".a" FILE.
C
      NPAD = 4096 - MOD(IDM*JDM,4096)
      IF     (NPAD.EQ.4096) THEN
        NPAD = 0
      ENDIF
C
      ALLOCATE( RK(IDM,JDM),
     +          DP(IDM,JDM),
     +          PK(IDM,JDM),
     +          PM(IDM,JDM),
     +          ZM(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_1st_isopyc: could not allocate ',
     +             8*IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
C
      IF     (NPAD.EQ.0) THEN
        INQUIRE( IOLENGTH=NRECL) RK
      ELSE
        INQUIRE( IOLENGTH=NRECL) RK,PAD(1:NPAD)
      ENDIF
*     write(6,*) 'nrecl = ',nrecl
#ifdef CRAY
#ifdef t3e
      IF     (MOD(NRECL,4096).EQ.0) THEN
        WRITE(CASN,8000) NRECL/4096
 8000   FORMAT('-F cachea:',I4.4,':1:0')
        IU8 = 11
        CALL ASNUNIT(IU8,CASN,IOS8)
        IF     (IOS8.NE.0) THEN
          WRITE(6,*) 'Error: can''t asnunit 11'
          WRITE(6,*) 'ios  = ',ios8
          WRITE(6,*) 'casn = ',casn
          CALL EXIT(5)
        ENDIF
        IU8 = 21
        CALL ASNUNIT(IU8,CASN,IOS8)
        IF     (IOS8.NE.0) THEN
          WRITE(6,*) 'Error: can''t asnunit 21'
          WRITE(6,*) 'ios  = ',ios8
          WRITE(6,*) 'casn = ',casn
          CALL EXIT(5)
        ENDIF
      ENDIF
#else
      CALL ASNUNIT(11,'-F syscall -N ieee',IOS)
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error: can''t asnunit 11'
        WRITE(6,*) 'ios = ',ios
        CALL EXIT(5)
      ENDIF
      CALL ASNUNIT(21,'-F syscall -N ieee',IOS)
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error: can''t asnunit 21'
        WRITE(6,*) 'ios = ',ios
        CALL EXIT(5)
      ENDIF
#endif
#endif
      OPEN(UNIT=11, FILE=CFILEA, FORM='UNFORMATTED', STATUS='OLD',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error: can''t open ',CFILEA(1:LEN_TRIM(CFILEA))
        WRITE(6,*) 'ios   = ',ios
        WRITE(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
C
C     OPEN OUTPUT UNITS (20 AND 21).
C
      OPEN(UNIT=21, FILE=CFILEM, FORM='UNFORMATTED', STATUS='NEW',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error: can''t open ',CFILEM(1:LEN_TRIM(CFILEM))
        WRITE(6,*) 'ios   = ',ios
        WRITE(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
      CFILEB = CFILEM(1:LEN_TRIM(CFILEM)-1) // 'b'
      OPEN(UNIT=20, FILE=CFILEB, FORM='FORMATTED', STATUS='NEW',
     +         IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error: can''t open ',CFILEB(1:LEN_TRIM(CFILEB))
        WRITE(6,*) 'ios   = ',ios
        CALL EXIT(3)
      ENDIF
C
C --- MIXED-LAYER
C
      CALL DAREAD(PM,IDM,JDM, 6, CFILEA)
C
      DO J= 1,JDM
        DO I= 1,IDM
          if     (i.eq.itest .and. j.eq.jtest) then
            write(6,'(a,f14.7)') 'mix = ',pm(i,j)*qonem
          endif
          IF     (PM(I,J).NE.SPVAL) THEN
            PK(I,J) = 0.0
            ZM(I,J) = 0.0
          ELSE
            ZM(I,J) = SPVAL
          ENDIF
        ENDDO
      ENDDO
C
C     ALL LAYERS
C
      DO K= 1,KDM
        KREC0 = NSURF+NLAY*(K-1)
        CALL DAREAD(DP,IDM,JDM, KREC0+3, CFILEA)
        CALL DAREAD(RK,IDM,JDM, KREC0+6, CFILEA)
        DO J= 1,JDM
          DO I= 1,IDM
            IF     (ZM(I,J).EQ.0.0) THEN
              PK(I,J) = PK(I,J) + DP(I,J)
              if     (i.eq.itest .and. j.eq.jtest) then
                write(6,'(a,i3,f14.7)') 'k,dr = ',
     +                                   k,rk(i,j)+thbase-sigma(k)
                write(6,'(a,i3,f14.7)') 'k,p  = ',k,pk(i,j)*qonem
              endif
              IF     (ITYPE.EQ.0) THEN
                IF     (PK(I,J).GE.PM(I,J)) THEN
                  ZM(I,J) = QONEM*PK(I,J)
                  if     (i.eq.itest .and. j.eq.jtest) then
                    write(6,'(a,f14.7)') 'iso = ',zm(i,j)
                  endif
                ENDIF
              ELSE
*               IF     (ABS(RK(I,J)+THBASE-SIGMA(K)).LT.0.005) THEN
                IF     (ABS(RK(I,J)+THBASE-SIGMA(K)).LT.0.01) THEN
                  ZM(I,J) = QONEM*(PK(I,J)-DP(I,J)) !p.k-1
                  if     (i.eq.itest .and. j.eq.jtest) then
                    write(6,'(a,f14.7)') 'iso = ',zm(i,j)
                  endif
                ENDIF
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      CLOSE(11)
C
      ZMIN = 1.0E10
      ZMAX = 0.0
      DO J= 1,JDM
        DO I= 1,IDM
          IF     (ZM(I,J).EQ.0.0) THEN
            ZM(I,J) = QONEM*PK(I,J)  ! no iso-surface above the bottom
          ENDIF
          IF     (ZM(I,J).NE.SPVAL) THEN
            ZMIN = MIN( ZMIN, ZM(I,J) )
            ZMAX = MAX( ZMAX, ZM(I,J) )
          ENDIF
        ENDDO
      ENDDO
C
C     OUTPUT THE INTERFACE.
C
      IF     (ITYPE.EQ.0) THEN
        WRITE(20,'(A,F12.2,2F10.2)') 
     +    'i.face below mix.lay.: day,min,max =',
     +    TIME,ZMIN,ZMAX
      ELSE
        WRITE(20,'(A,F12.2,2F10.2)') 
     +    '1st isopycnal i.face: day,min,max =',
     +    TIME,ZMIN,ZMAX
      ENDIF
      IF     (NPAD.EQ.0) THEN
        WRITE(21,REC=1) ZM
      ELSE
        PAD(1:NPAD) = SPVAL
        WRITE(21,REC=1) ZM,PAD(1:NPAD)
      ENDIF
      CLOSE(20)
      CLOSE(21)
      END
      SUBROUTINE READ_B(CFILEB,
     &                  IEXPT,YRFLAG,IDM,JDM,KDM,NSURF,NLAY,
     &                  THBASE,SIGMA,TIME)
      IMPLICIT NONE
C
      INTEGER       IEXPT,YRFLAG,IDM,JDM,KDM,NSURF,NLAY
      REAL          THBASE,SIGMA(99),TIME
      CHARACTER*240 CFILEB
C
C     EXTRACT NEEDED MODEL PARAMETERS FROM ARCHIVE .b FILE.
C
      INTEGER      IDUM,IOS,K,L,NSTEP
      REAL         THBASE_IN
      CHARACTER*6  CVARIN*6
      CHARACTER*240 CLINE
C
      OPEN(UNIT=12, FILE=CFILEB, FORM='FORMATTED', STATUS='OLD',
     +         IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error: can''t open ',CFILEB(1:LEN_TRIM(CFILEB))
        WRITE(6,*) 'ios   = ',ios
        CALL EXIT(3)
      ENDIF
      READ(12,*)  ! skip title(1)
      READ(12,*)  ! skip title(2)
      READ(12,*)  ! skip title(3)
      READ(12,*)  ! skip title(4)
      READ(12,*)  ! skip iversn
      READ(12,*) IEXPT,CVARIN
      IF     (CVARIN.NE.'iexpt ') THEN
        WRITE(6,*) 'Error in hycom_profile: bad .b file'
        WRITE(6,*) 'filename: ',CFILEB(1:LEN_TRIM(CFILEB))
        CALL EXIT(4)
      endif
      READ(12,*) YRFLAG
      READ(12,*) IDM
      READ(12,*) JDM
C
C     FIND KDM.
C
      NSURF = 13  ! number of surface arrays
      NLAY  =  6  ! number of arrays per layer
C
      DO K= 1,10
        READ(12,'(a)') CLINE
      ENDDO
*     write(6,*) trim(cline)
      IF     (CLINE(1:8).EQ.'thmix   ') THEN
        READ(CLINE(36:42),*) THBASE_IN
        IF     (THBASE_IN.NE.0.0) THEN
          THBASE = THBASE_IN
        ELSE
          THBASE = 25.0
        ENDIF
      ELSE
        WRITE(6,*) 
        WRITE(6,*) 'Expected thmix but got:'
        WRITE(6,*) CLINE(1:LEN_TRIM(CLINE))
        WRITE(6,*) 
        CALL EXIT(2)
      ENDIF
      DO K= 11,13
        READ(12,'(a)') CLINE
      ENDDO
*     write(6,*) trim(cline)
      IF     (CLINE(1:8).EQ.'kemix   ') THEN
        NLAY  = 7  ! mean archive
        NSURF = NSURF + 1
        READ(12,'(a)') CLINE
      ENDIF
*     write(6,*) trim(cline)
      IF     (CLINE(1:8).EQ.'covice  ') THEN
        NSURF = NSURF + 3
        READ(12,'(a)') CLINE
        READ(12,'(a)') CLINE
        READ(12,'(a)') CLINE
      ENDIF
*     write(6,*) trim(cline)
      READ(12,'(a)') CLINE
      IF     (NLAY.EQ.7) THEN
        NSURF = NSURF + 1
        READ(12,'(a)') CLINE  ! kebtrop
      ENDIF
*     write(6,*) trim(cline)
C
      L = INDEX(CLINE,'=')
      READ(CLINE(L+1:),*) NSTEP,TIME
C
C     FIND NLAY (ALLOWING FOR TRACERS)
C
      READ(12,'(a)') CLINE
      DO L= 2,99
        READ(12,'(a)',IOSTAT=IOS) CLINE
        IF     (IOS.NE.0) THEN
          EXIT
        ELSEIF (CLINE(1:8).EQ.'u-vel.  ') THEN
*         write(6,*) trim(cline)
          EXIT
        ENDIF
      ENDDO
      NLAY = L-1
C
      REWIND(UNIT=12)
      DO K= 1,NSURF+10
        READ(12,'(a)') CLINE
*       write(6,*) trim(cline)
      ENDDO
*     write(6,*) '-----------------------------------------------'
C
      DO K= 1,999
        READ(12,'(a)',IOSTAT=IOS) CLINE
        IF     (IOS.NE.0) THEN
          EXIT
        ELSEIF (CLINE(1:8).NE.'u-vel.  ') THEN
*         write(6,*) trim(cline)
          EXIT
        ENDIF
        L = INDEX(CLINE,'=')
        READ(CLINE(L+1:),*) NSTEP,TIME,IDUM,SIGMA(K)
*       write(6,*) trim(cline),"  (1)",sigma(k)
*       write(6,*) "k,sigma",k,sigma(k)
C
        DO L= 2,NLAY
          READ(12,'(a)') CLINE
        ENDDO
      ENDDO
      KDM = K-1
*     write(6,*) 'kdm = ',kdm
      CLOSE(UNIT=12)
      RETURN
      END
      SUBROUTINE DAREAD(A,IDM,JDM, KREC, CFILEA)
      IMPLICIT NONE
C
      CHARACTER*240 CFILEA
      INTEGER       IDM,JDM,KREC
      REAL*4        A(IDM,JDM)
C
C --- READ ONE RECORD ON UNIT 11
C
      INTEGER IOS
C
      READ(11,REC=KREC,IOSTAT=IOS) A
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'can''t read record ',KREC,
     +             ' from file ',TRIM(CFILEA)
        CALL EXIT(4)
        STOP
      ENDIF
      END

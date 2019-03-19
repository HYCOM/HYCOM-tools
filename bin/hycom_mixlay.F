      PROGRAM MIXLAY
      IMPLICIT NONE
C
C  hycom_mixlay - Usage:  hycom_mixlay archv.a mljmp [[ptype] mixtype] ml.a
C
C                 generates a mixed layer from a HYCOM archive file.
C
C   archv.a is assumed to be an HYCOM archive data file, with companion
C   header file archv.b.  Both standard and mean archive files are allowed.
C
C   mixed layer set to interpolated depth where the temperature or density
C   difference from the surface is mljmp.  In the default case, mljmp is
C   input as a temperature difference but converted to the equivalent
C   density difference at the surface T&S.
C
C   ptype is the profile method (default 1)
C           =-1; interpolation between cell centers
C           = 1; piecewise linear method (PLM)
C   mixtype is the mixed layer criteria (default 0)
C           = 0; mljmp is temperature but converted to density at (SST,SSS)
C           = 1; mljmp is temperature
C           = 2; mljmp is density
C
C   if mixtype is input then ptype must also be specified.
C
C   ml.a and ml.b will contain the resulting mixed layer field.
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  April 2003.
C
      REAL*4     QONEM,SPVAL
      PARAMETER (QONEM=1.0/9806.0, SPVAL=2.0**100)
C
      REAL*4, ALLOCATABLE :: RM(:,:),RK(:,:,:)
      REAL*4, ALLOCATABLE ::         PO(:,:),PN(:,:)
      REAL*4, ALLOCATABLE :: ZM(:,:),ZO(:,:),ZN(:,:)
      REAL*4              :: PAD(4096)
C
      INTEGER       IARGC
      INTEGER       NARG
      CHARACTER*240 CARG
C
      INTEGER       IDM,JDM,KDM,NSURF,NLAY,NTRC,IEXPT,YRFLAG
      INTEGER       NPAD,ITYPE,MTYPE,ITEST,JTEST
      REAL          TMLJMP,THBASE,TIME
      CHARACTER*30  CMTYPE
      CHARACTER*240 CFILEA,CFILEB,CFILEM
C
      CHARACTER*18  CASN
      INTEGER       I,J,K,KREC,KREC0,IOS,NRECL
      REAL          PK,RIA,RIB,Q,ZMIN,ZMAX
#ifdef CRAY
      INTEGER*8     IU8,IOS8
#endif
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.3) THEN
        CALL GETARG(1,CFILEA)
        CALL GETARG(2,CARG)
        READ(CARG,*) TMLJMP
        CALL GETARG(3,CFILEM)
        ITYPE = 1
        MTYPE = 0
        ITEST = 0
        JTEST = 0
      ELSEIF (NARG.EQ.4) THEN
        CALL GETARG(1,CFILEA)
        CALL GETARG(2,CARG)
        READ(CARG,*) TMLJMP
        CALL GETARG(3,CARG)
        READ(CARG,*) ITYPE
        CALL GETARG(4,CFILEM)
        MTYPE = 0
        ITEST = 0
        JTEST = 0
      ELSEIF (NARG.EQ.5) THEN
        CALL GETARG(1,CFILEA)
        CALL GETARG(2,CARG)
        READ(CARG,*) TMLJMP
        CALL GETARG(3,CARG)
        READ(CARG,*) ITYPE
        CALL GETARG(4,CARG)
        READ(CARG,*) MTYPE
        CALL GETARG(5,CFILEM)
        ITEST = 0
        JTEST = 0
      ELSEIF (NARG.EQ.7) THEN  !undocumented, for debugging
        CALL GETARG(1,CFILEA)
        CALL GETARG(2,CARG)
        READ(CARG,*) TMLJMP
        CALL GETARG(3,CARG)
        READ(CARG,*) ITYPE
        CALL GETARG(4,CARG)
        READ(CARG,*) MTYPE
        CALL GETARG(5,CFILEM)
        CALL GETARG(6,CARG)
        READ(CARG,*) ITEST
        CALL GETARG(7,CARG)
        READ(CARG,*) JTEST
      ELSE
        WRITE(6,*) 
     +    'Usage: hycom_mixlay archv.a mljmp [[ptype] mixtype] ml.a'
        CALL EXIT(1)
      ENDIF
C
C     EXTRACT MODEL PARAMETERS FROM ".b" FILE.
C
      CFILEB = CFILEA(1:LEN_TRIM(CFILEA)-1) // 'b'
      CALL READ_B(CFILEB,
     +            IEXPT,YRFLAG,IDM,JDM,KDM,NSURF,NLAY,NTRC,
     +            THBASE,TIME)
C
C     OPEN ".a" FILE.
C
      NPAD = 4096 - MOD(IDM*JDM,4096)
      IF     (NPAD.EQ.4096) THEN
        NPAD = 0
      ENDIF
C
      ALLOCATE( RM(IDM,JDM), RK(IDM,JDM,-2:1),
     +                       PO(IDM,JDM), PN(IDM,JDM),
     +          ZM(IDM,JDM), ZO(IDM,JDM), ZN(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_mixlay: could not allocate ',
     +             8*IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
C
      RK(:,:,:) = 0.0
C
      IF     (NPAD.EQ.0) THEN
        INQUIRE( IOLENGTH=NRECL) RM
      ELSE
        INQUIRE( IOLENGTH=NRECL) RM,PAD(1:NPAD)
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
C     FIRST LAYER.
C
      KREC0 = NSURF - NTRC
      CALL DAREAD(PN,       IDM,JDM, KREC0+NLAY-3,    CFILEA)
      CALL DAREAD(ZN,       IDM,JDM, KREC0+NLAY-2,    CFILEA)  !TN
      CALL DAREAD(ZO,       IDM,JDM, KREC0+NLAY-1,    CFILEA)  !SN
      IF     (MTYPE.EQ.0) THEN
        IF     (THBASE.LT.30.0) THEN
          CALL RMLJMP0(RM, ZN,ZO, IDM,JDM, TMLJMP) !sigma0
        ELSE
          CALL RMLJMP2(RM, ZN,ZO, IDM,JDM, TMLJMP) !sigma2
        ENDIF
        CALL DAREAD(RK(1,1,0),IDM,JDM, KREC0+NLAY,        CFILEA)
        CALL DAREAD(RK(1,1,1),IDM,JDM, KREC0+NLAY+NLAY,   CFILEA)
      ELSEIF (MTYPE.EQ.1) THEN
        RM(:,:)   = TMLJMP  !temperature
        CALL DAREAD(RK(1,1,1),IDM,JDM, KREC0+NLAY+NLAY-2, CFILEA)
        DO J= 1,JDM
          DO I= 1,IDM
            IF     (RK(I,J,1).NE.SPVAL) THEN
              RK(I,J,0) = -ZN(I,J)
              RK(I,J,1) = -RK(I,J,1)
            ELSE
              RK(I,J,0) = SPVAL
            ENDIF
          ENDDO
        ENDDO
      ELSE   !MTYPE.EQ.2
        RM(:,:)   = TMLJMP  !density
        CALL DAREAD(RK(1,1,0),IDM,JDM, KREC0+NLAY,        CFILEA)
        CALL DAREAD(RK(1,1,1),IDM,JDM, KREC0+NLAY+NLAY,   CFILEA)
      ENDIF
      if     (itest.ne.0) then
        write(6,'(a,i3,f14.7)')  'k,t = ',1,zn(itest,jtest)
        write(6,'(a,i3,f14.7)')  'k,s = ',1,zo(itest,jtest)
        write(6,'(a,i3,2f14.7)') 'k,r = ',1,rm(itest,jtest)
      endif
C
      DO J= 1,JDM
        DO I= 1,IDM
          IF     (RK(I,J,0).NE.SPVAL) THEN
            RM(I,J) = RK(I,J,0) + MAX(RM(I,J),0.001)
            ZN(I,J) = 0.5*QONEM*PN(I,J)
            ZM(I,J) = 0.0
          ELSE
            RM(I,J) = SPVAL
            ZN(I,J) = SPVAL
            ZM(I,J) = SPVAL
          ENDIF
        ENDDO
      ENDDO
      if     (itest.ne.0) then
        write(6,'(a,i3,2f14.7)') 'k,r = ',1,rk(itest,jtest,0),
     +                                      rm(itest,jtest)
        write(6,'(a,i3,f14.7)')  'k,z = ',1,zn(itest,jtest)
      endif
C
C     REST OF THE LAYERS
C
      DO K= 2,KDM
        DO J= 1,JDM
          DO I= 1,IDM
            IF     (ZM(I,J).EQ.0.0) THEN
              RK(I,J,-2) = RK(I,J,-1)
              RK(I,J,-1) = RK(I,J, 0)
              RK(I,J, 0) = RK(I,J, 1)
              PO(I,J)    = PN(I,J)
              ZO(I,J)    = ZN(I,J)
            ENDIF
          ENDDO
        ENDDO
        KREC0 = NSURF-NTRC+NLAY*(K-1)
        CALL DAREAD(PN,       IDM,JDM, KREC0+NLAY-3,    CFILEA)
        IF     (K.NE.KDM) THEN
          IF     (MTYPE.NE.1) THEN  !density
            CALL DAREAD(RK(1,1,1),IDM,JDM, KREC0+NLAY+NLAY,   CFILEA)
          ELSE !temperature
            CALL DAREAD(RK(1,1,1),IDM,JDM, KREC0+NLAY+NLAY-2, CFILEA)
            DO J= 1,JDM
              DO I= 1,IDM
                IF     (RK(I,J,1).NE.SPVAL) THEN
                  RK(I,J,1) = -RK(I,J,1)
                ENDIF
              ENDDO
            ENDDO
          ENDIF
        ENDIF
        DO J= 1,JDM
          DO I= 1,IDM
            IF     (ZM(I,J).EQ.0.0) THEN
              PN(I,J) = PO(I,J) + PN(I,J)
              ZN(I,J) = 0.5*QONEM*(PO(I,J) + PN(I,J))
              if     (i.eq.itest .and. j.eq.jtest) then
                write(6,'(a,i3,f14.7)') 'k,r = ',k,rk(i,j,0)
                write(6,'(a,i3,f14.7)') 'k,z = ',k,zn(i,j)
              endif
              IF     (RK(I,J,0).GE.RM(I,J)) THEN
                IF     (ITYPE.EQ.-1) THEN
C
C                 CELL CENTER METHOD: MLD IS BETWEEN ZO AND ZN.
C
                  Q = (RK(I,J,0)-RM(I,J))/(RK(I,J,0)-RK(I,J,-1))
                  ZM(I,J) = Q*ZO(I,J) + (1.0-Q)*ZN(I,J)
                  if     (i.eq.itest .and. j.eq.jtest) then
                    write(6,'(a,i3,f14.7)') 'k,q = ',k,q
                    write(6,'(a,i3,f14.7)') 'k,z = ',k,zm(i,j)
                  endif
                ELSE  !plm
C
C                 FIND THE TWO DENSITIES ON INTERFACE BETWEEN LAYERS K-1 AND K,
C                 USING PLM BUT ASSUMING LAYERS K-2 AND K+1 ARE PCM.
C
                  IF     (K.EQ.2 .OR. RK(I,J,-2).GE.RK(I,J,-1)) THEN
                    RIA = RK(I,J,-1)
                  ELSE
                    RIA = RK(I,J,-1) + MIN( RK(I,J,-1)-RK(I,J,-2),
     +                                      RK(I,J, 0)-RK(I,J,-1) )
                  ENDIF
                  IF     (K.EQ.KDM .OR. RK(I,J,+1).LE.RK(I,J,0)) THEN
                    RIB = RK(I,J,0)
                  ELSE
                    RIB = RK(I,J,0) - MIN( RK(I,J,+1)-RK(I,J, 0),
     +                                     RK(I,J, 0)-RK(I,J,-1) )
                  ENDIF
                  IF     (RIA.GT.RIB) THEN  !favor the upper layer
                    RIB = RIA
                  ENDIF
C
                  PK = QONEM*PO(I,J)
C
                  IF     (RIA.GT.RM(I,J)) THEN
C
C                   MLD IS BETWEEN ZO AND P(K)
C
                    Q       = (RIA-RM(I,J))/(RIA-RK(I,J,-1))  !R(K-1)<RIA
                    ZM(I,J) = Q*ZO(I,J) + (1.0-Q)*PK
                    if     (i.eq.itest .and. j.eq.jtest) then
                      write(6,'(a,i3,f14.7)') '-k,q = ',k,q
                      write(6,'(a,i3,f14.7)') '-k,z = ',k,zm(i,j)
                    endif
                  ELSEIF (RIB.LT.RM(I,J)) THEN
C
C                   MLD IS BETWEEN P(K) AND ZN
C
                    Q       = (RK(I,J,0)-RM(I,J))/(RK(I,J,0)-RIB)  !RIB<R(K)
                    ZM(I,J) = Q*PK + (1.0-Q)*ZN(I,J)
                    if     (i.eq.itest .and. j.eq.jtest) then
                      write(6,'(a,i3,f14.7)') '+k,q = ',k,q
                      write(6,'(a,i3,f14.7)') '+k,z = ',k,zm(i,j)
                    endif
                  ELSE
C
C                   MLD IS AT P(K)
C
                    ZM(I,J) = PK
                    if     (i.eq.itest .and. j.eq.jtest) then
                      write(6,'(a,i3,f14.7)') ' k,z = ',k,zm(i,j)
                    endif
                  ENDIF
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
            ZM(I,J) = QONEM*PN(I,J)  ! mixed layer at bottom
          ENDIF
          IF     (ZM(I,J).NE.SPVAL) THEN
            ZMIN = MIN( ZMIN, ZM(I,J) )
            ZMAX = MAX( ZMAX, ZM(I,J) )
          ENDIF
        ENDDO
      ENDDO
C
C     OUTPUT THE MIXED LAYER.
C
      IF     (MTYPE.EQ.0) THEN
        CMTYPE = 'hycom mixed layer'
      ELSEIF (MTYPE.EQ.1) THEN
        CMTYPE = 'temperature mixed layer'
      ELSE   !MTYPE.EQ.2
        CMTYPE = 'density mixed layer'
      ENDIF
      IF     (ITYPE.EQ.-1) THEN
        WRITE(20,'(A,A,F6.3,F12.2,2F10.2)') 
     +    TRIM(CMTYPE),' (cen): mljmp,day,min,max =',
     +    TMLJMP,TIME,ZMIN,ZMAX
      ELSE
        WRITE(20,'(A,A,F6.3,F12.2,2F10.2)') 
     +    TRIM(CMTYPE),' (plm): mljmp,day,min,max =',
     +    TMLJMP,TIME,ZMIN,ZMAX
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
     &                  IEXPT,YRFLAG,IDM,JDM,KDM,NSURF,NLAY,NTRC,
     &                  THBASE,TIME)
      IMPLICIT NONE
C
      INTEGER      IEXPT,YRFLAG,IDM,JDM,KDM,NSURF,NLAY,NTRC
      REAL         THBASE,TIME
      CHARACTER*240 CFILEB
C
C     EXTRACT NEEDED MODEL PARAMETERS FROM ARCHIVE .b FILE.
C
      INTEGER      IOS,K,KTR,NSTEP
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
*     write(6,*) cline(1:len_trim(cline))
      IF     (CLINE(1:8).EQ.'thmix   ') THEN
        READ(CLINE(36:42),*) THBASE_IN
        IF     (THBASE_IN.NE.0.0) THEN
          THBASE = THBASE_IN
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
*     write(6,*) cline(1:len_trim(cline))
      IF     (CLINE(1:8).EQ.'kemix   ') THEN
        NLAY  = 7  ! mean archive
        NSURF = NSURF + 1
        READ(12,'(a)') CLINE
      ENDIF
*     write(6,*) cline(1:len_trim(cline))
      IF     (CLINE(1:8).EQ.'covice  ') THEN
        NSURF = NSURF + 3
        READ(12,'(a)') CLINE
        READ(12,'(a)') CLINE
        READ(12,'(a)') CLINE
      ENDIF
*     write(6,*) cline(1:len_trim(cline))
      READ(12,'(a)') CLINE
      IF     (NLAY.EQ.7) THEN
        NSURF = NSURF + 1
        READ(12,'(a)') CLINE  ! kebtrop
      ENDIF
*     write(6,*) cline(1:len_trim(cline))
C
      K = INDEX(CLINE,'=')
      READ(CLINE(K+1:),*) NSTEP,TIME
C
C     FIRST LAYER (FIND NUMBER OF LAYERS)
C
      READ(12,'(a)') CLINE
      READ(12,'(a)') CLINE
      READ(12,'(a)') CLINE
      READ(12,'(a)') CLINE
      READ(12,'(a)') CLINE
      READ(12,'(a)') CLINE
      IF     (NLAY.EQ.7) THEN
        READ(12,'(a)') CLINE
      ENDIF
      DO KTR= 1,9999
        READ(12,'(a)',IOSTAT=IOS) CLINE
        IF     (IOS.NE.0) THEN
C
C         EXACTLY ONE LAYER
C
          NTRC = KTR-1  !includes viscty,t-diff,s-diff
          NLAY = NLAY + NTRC
          KDM  = 1
*         write(6,*) 'kdm,ndif,ntrc = ',kdm,ndif,ntrc
          CLOSE(UNIT=12)
          RETURN
        ELSEIF (CLINE(1:8).EQ.'viscty  ') THEN
          CYCLE
        ELSEIF (CLINE(1:8).EQ.'t-diff  ') THEN
          CYCLE
        ELSEIF (CLINE(1:8).EQ.'s-diff  ') THEN
          CYCLE
        ELSEIF (CLINE(1:8).NE.'tracer  ') THEN
          EXIT
        ENDIF
      ENDDO
      NTRC = KTR-1  !includes viscty,t-diff,s-diff
      NLAY = NLAY + NTRC
C
C     MORE THAN 1 LAYER
C
      DO KTR= 1,NLAY-1
        READ(12,'(a)') CLINE
      ENDDO
C
      DO K= 3,999
        READ(12,'(a)',IOSTAT=IOS) CLINE
        IF     (IOS.NE.0) THEN
          EXIT
        ELSEIF (CLINE(1:8).NE.'u-vel.  ') THEN
*         write(6,*) cline(1:len_trim(cline))
          EXIT
        ENDIF
*       write(6,*) cline(1:len_trim(cline)),"  (1)"
C
        DO KTR= 1,NLAY-1
          READ(12,'(a)') CLINE
        ENDDO
      ENDDO !k
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
      SUBROUTINE RMLJMP0(RMLJMP, SST,SSS, IDM,JDM, TMLJMP)
      IMPLICIT NONE
C
      INTEGER       IDM,JDM
      REAL*4        RMLJMP(IDM,JDM),SST(IDM,JDM),SSS(IDM,JDM),TMLJMP
C
C --- CONVERT A TEMPERATURE JUMP INTO A DENSITY JUMP
C
      REAL*4     SPVAL
      PARAMETER (SPVAL=2.0**100)
C
      INTEGER I,J
C
c-----------------------------------------------------------------------------
      real*4 dsigdt
      real*4 s,t
c
c --- coefficients for sigma-0 (based on Brydon & Sun fit)
      real*4     c1,c2,c3,c4,c5,c6,c7
      parameter (c1=-1.36471E-01, c2= 4.68181E-02, c3= 8.07004E-01,
     &           c4=-7.45353E-03, c5=-2.94418E-03,
     &           c6= 3.43570E-05, c7= 3.48658E-05)
c
c --- sigma-theta as a function of temp (deg c) and salinity (mil)
c --- (friedrich-levitus 3rd degree polynomial fit)
c
c     sig(t,s)=(c1+c3*s+t*(c2+c5*s+t*(c4+c7*s+c6*t)))
c
c --- d(sig)/dt
      dsigdt(t,s)=(c2+c5*s+2.*t*(c4+c7*s+1.5*c6*t))
c-----------------------------------------------------------------------------
C
      DO J= 1,JDM
        DO I= 1,IDM
          IF     (SST(I,J).NE.SPVAL) THEN
            RMLJMP(I,J) = -TMLJMP*DSIGDT(SST(I,J),SSS(I,J))
            RMLJMP(I,J) = MAX(RMLJMP(I,J),TMLJMP*0.1) !cold-water fix
          ELSE
            RMLJMP(I,J) = 0.0
          ENDIF
        ENDDO
      ENDDO
      END
      SUBROUTINE RMLJMP2(RMLJMP, SST,SSS, IDM,JDM, TMLJMP)
      IMPLICIT NONE
C
      INTEGER       IDM,JDM
      REAL*4        RMLJMP(IDM,JDM),SST(IDM,JDM),SSS(IDM,JDM),TMLJMP
C
C --- CONVERT A TEMPERATURE JUMP INTO A DENSITY JUMP
C
      REAL*4     SPVAL
      PARAMETER (SPVAL=2.0**100)
C
      INTEGER I,J
C
c-----------------------------------------------------------------------------
      real*4 dsigdt
      real*4 s,t
c
c --- coefficients for sigma-2 (based on Brydon & Sun fit)
      real*4     c1,c2,c3,c4,c5,c6,c7
      parameter (c1= 9.77093E+00, c2=-2.26493E-02, c3= 7.89879E-01,
     &           c4=-6.43205E-03, c5=-2.62983E-03,
     &           c6= 2.75835E-05, c7= 3.15235E-05)
c
c --- sigma-theta as a function of temp (deg c) and salinity (mil)
c --- (friedrich-levitus 3rd degree polynomial fit)
c
c     sig(t,s)=(c1+c3*s+t*(c2+c5*s+t*(c4+c7*s+c6*t)))
c
c --- d(sig)/dt
      dsigdt(t,s)=(c2+c5*s+2.*t*(c4+c7*s+1.5*c6*t))
c-----------------------------------------------------------------------------
C
      DO J= 1,JDM
        DO I= 1,IDM
          IF     (SST(I,J).NE.SPVAL) THEN
            RMLJMP(I,J) = -TMLJMP*DSIGDT(SST(I,J),SSS(I,J))
            RMLJMP(I,J) = MAX(RMLJMP(I,J),TMLJMP*0.1) !cold-water fix
          ELSE
            RMLJMP(I,J) = 0.0
          ENDIF
        ENDDO
      ENDDO
      END

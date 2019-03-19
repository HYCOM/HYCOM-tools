      PROGRAM HYCOM_SAMPLE_XY
      IMPLICIT NONE
C
C  hycom_sample_xy - Usage:  hycom_sample_xy file.a idm jdm cfmt n list.txt [void]
C
C                 prints the n points (xp,xp), one per line in list.txt,
C                 from all (1:idm,1:jdm) arrays in file.a
C                 cfmt     - output format, e.g. "(f6.2)" or "(*)"
C                 n        - number of xp,yp points in list.txt
C                 list.txt - list of xp,yp points, one per line
C                 void     - data void value (default 2.0**100)
C
C  file.a is assumed to contain idm*jdm 32-bit IEEE real values for
C   each array, in standard f77 element order, followed by padding
C   to a multiple of 4096 32-bit words, but otherwise with no control
C   bytes/words, and input values of 2.0**100 indicating a data void.
C
C  Each of the n lines in list.txt should contain two values: x,y.
C  The x,y values are the p-grid location of the point,  note that
C  hycom/ALL/bin/hycom_lonlat2xy will convert lon,lat to x,y.
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  August 2010.
C
      REAL*4     SPVAL
      PARAMETER (SPVAL=2.0**100)
C
      REAL*8,       ALLOCATABLE :: XP(:),YP(:)
      REAL*4,       ALLOCATABLE :: AP(:)
      REAL*4,       ALLOCATABLE :: A(:,:)
      REAL*4                    :: PAD(4096),VOID
      INTEGER       IOS
      INTEGER       IARGC
      INTEGER       NARG
      CHARACTER*240 CARG
C
      INTEGER       IDM,JDM,N,NPAD
      CHARACTER*240 CFILE,CFILEO,CFMT,CLIST
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.6) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,*) IDM
        CALL GETARG(3,CARG)
        READ(CARG,*) JDM
        CALL GETARG(4,CFMT)
        CALL GETARG(5,CARG)
        READ(CARG,*) N
        CALL GETARG(6,CLIST)
        VOID = SPVAL
      ELSEIF (NARG.EQ.7) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,*) IDM
        CALL GETARG(3,CARG)
        READ(CARG,*) JDM
        CALL GETARG(4,CFMT)
        CALL GETARG(5,CARG)
        READ(CARG,*) N
        CALL GETARG(6,CLIST)
        CALL GETARG(7,CARG)
        READ(CARG,*) VOID
      ELSE
        WRITE(6,*) 
     +   'Usage:  hycom_sample_xy ' //
     +            'file.a idm jdm cfmt n list.txt [void]'
        CALL EXIT(1)
      ENDIF
C
      NPAD = 4096 - MOD(IDM*JDM,4096)
      IF     (NPAD.EQ.4096) THEN
        NPAD = 0
      ENDIF
C
      ALLOCATE( A(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_sample_xy: could not allocate ',
     +             IDM*JDM,' words for A'
        CALL EXIT(2)
      ENDIF
C
      ALLOCATE( AP(N),
     +          XP(N),
     +          YP(N), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_profile_xy: could not allocate ',
     +             '(1:n) arrays'
        CALL EXIT(2)
      ENDIF

C
      CALL HSAMPL(A,IDM,JDM,PAD,NPAD,
     &            XP,YP,AP,N, VOID, CFILE,CFMT,CLIST)
      CALL EXIT(0)
      END
      SUBROUTINE HSAMPL(A,IDM,JDM,PAD,NPAD,
     &                  XP,YP,AP,N, VOID, CFILE,CFMT,CLIST)
      IMPLICIT NONE
C
      REAL*4     SPVAL
      PARAMETER (SPVAL=2.0**100)
C
      CHARACTER*240 CFILE,CFMT,CLIST
      INTEGER       IDM,JDM,NPAD,N
      REAL*8        XP(N),YP(N)
      REAL*4        AP(N),VOID
      REAL*4        A(IDM,JDM),PAD(NPAD)
C
C     MOST OF WORK IS DONE HERE.
C
#ifdef sun
      INTEGER      IR_ISNAN
C
#endif
      CHARACTER*18 CASN
      INTEGER      LEN_TRIM
      INTEGER      I,II,IP,IPP1,IPP2,J,JJ,JP,JPP1,JPP2,LANDF,LANDL,
     &             IOS,K,NRECL
      REAL*4       DX,DY
#ifdef CRAY
      INTEGER*8    IU8,IOS8
#endif
C
C     DEFINE LIST OF POINTS
C
      OPEN(UNIT=99, FILE=CLIST, FORM='FORMATTED', STATUS='OLD',
     +         IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error: can''t open ',TRIM(CLIST)
        WRITE(6,*) 'ios   = ',ios
        CALL EXIT(3)
      ENDIF
C
      DO I= 1,N
        READ(99,*) XP(I),YP(I)
      ENDDO !i
C
      INQUIRE( IOLENGTH=NRECL) A,PAD
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
#else
      CALL ASNUNIT(11,'-F syscall -N ieee',IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t asnunit 11'
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
C
      WRITE(6,6000) 
     +  TRIM(CFILE),
     +  IDM,JDM,
     +  TRIM(CFMT),
     +  N,
     +  TRIM(CLIST),
     +  VOID
      DO I= 1,N
        WRITE(6,6050) XP(I),YP(I)
      ENDDO
C
      DO K= 1,999999
        READ(11,REC=K,IOSTAT=IOS) A
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
        IF     (IOS.NE.0) THEN
          EXIT
        ENDIF
C
        DO I= 1,N
          IP = XP(I)
          JP = YP(I)
          IF     (IP.EQ.0 .AND. JP.EQ.0) THEN
C
C           STANDARD MARKER FOR POINTS OUTSIDE THE GRID.
C
            AP(I) = VOID
            CYCLE
          ENDIF
C
          DX = XP(I) - IP
          DY = YP(I) - JP
          IF     (DX.EQ.0.0) THEN !1-d interpolation
            IPP1 = IP
          ELSEIF (DX.EQ.1.0) THEN !1-d interpolation
            IPP1 = IP+1
            IP   = IPP1
          ELSEIF (IP.EQ.IDM) THEN
            IPP1 = 1
          ELSE
            IPP1 = IP+1
          ENDIF
          IF     (DY.EQ.0.0) THEN !1-d interpolation
            JPP1 = JP
          ELSEIF (DY.EQ.1.0) THEN !1-d interpolation
            JPP1 = JP+1
            JP   = JPP1
          ELSE
            JPP1 = JP+1
          ENDIF
          IF     (K.EQ.1) THEN
C
C ---       DETECT NEARBY LAND.
C ---       ONLY ONCE, I.E. ASSUME ALL FIELDS HAVE THE SAME DATA VOIDS.
C
            IF     (MAX( A(IP,  JP),
     +                   A(IP,  JPP1),
     +                   A(IPP1,JP),
     +                   A(IPP1,JPP1) ).GE.SPVAL) THEN
C
C             NEAR LAND, SO USE NEAREST POINT.
C
              IF     (DX.LE.0.5 .AND. DY.LE.0.5) THEN
                IP = IP
                JP = JP
              ELSEIF (DX.LE.0.5 .AND. DY.GT.0.5) THEN
                IP = IP
                JP = JPP1
              ELSEIF (DX.GT.0.5 .AND. DY.LE.0.5) THEN
                IP = IPP1
                JP = JP
              ELSE ! (DX.GT.0.5 .AND. DY.GT.0.5) THEN
                IP = IPP1
                JP = JPP1
              ENDIF
              XP(I) = IP  !for all subsequent interpolations
              YP(I) = JP  !for all subsequent interpolations
              DX    = 0.0
              DY    = 0.0
              IPP1  = IP
              JPP1  = JP
            ENDIF
          ENDIF !K=1
C
          IF     (A(IP,JP).GE.SPVAL) THEN
C
C           NEAREST POINT IS LAND.
C
            AP(I) = VOID
          ELSE  !all referenced points (1, 2 or 4 unique points) are ok
            AP(I) = (1.0-DX)*(1.0-DY)*A(IP,  JP)   +
     +              (1.0-DX)*     DY *A(IP,  JPP1) +
     +                   DX *(1.0-DY)*A(IPP1,JP)   +
     +                   DX *     DY *A(IPP1,JPP1)
          ENDIF
        ENDDO !i
C
        WRITE(6,CFMT) AP(1:N)
      ENDDO
      RETURN
 6000 FORMAT(
     +   '# hycom_sample_xy file.a idm jdm cfmt n list.txt void' /
     +   '# hycom_sample_xy ',A     /
     +   '#                 ',I6,I6 /
     +   '#                 ',A     /
     +   '#                 ',I12   /
     +   '#                 ',A     /
     +   '#                 ',D12.5  )
 6050 FORMAT(
     +   '#                 ',2D12.5 )
      END

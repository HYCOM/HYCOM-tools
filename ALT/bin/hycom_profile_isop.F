      PROGRAM HYCOM_PROFILE_ISOP
      IMPLICIT NONE
C
C  hycom_profile_isop - Usage:  hycom_profile_isop syn_prof_xtnd.isp nprof archo.txt [type]
C
C                 convert profile number nprof T&S from ISOP to archive
C
C   syn_prof_xtnd.isp is an ISOP syn_prof_xtnd file
C   nprof             is the profile number
C   archo.txt         will be the output text profile file
C   type              type: isop=1, background=2, climo=3 (default 1)
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  December 2015.
C
      INTEGER       IARGC
      INTEGER       NARG
      CHARACTER*240 CARG
C
      CHARACTER*240 CFILEI,CFILEO,CFORMAT
      CHARACTER*240 CLINE
      REAL          THK,DEPTH,FLAG,PLAT,PLON,DUM
      INTEGER       IOS,K,KDM,KI,KK,KP,SIGVER,NPROF,ITYPE
C
      REAL          SI(99,5),Z(99),P(100)
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.3) THEN
        CALL GETARG(1,CFILEI)
        CALL GETARG(2,CLINE)
        READ(CLINE,*) NPROF
        CALL GETARG(3,CFILEO)
        ITYPE = 1
      ELSEIF (NARG.EQ.4) THEN
        CALL GETARG(1,CFILEI)
        CALL GETARG(2,CLINE)
        READ(CLINE,*) NPROF
        CALL GETARG(3,CFILEO)
        CALL GETARG(4,CLINE)
        READ(CLINE,*) ITYPE
      ELSE
        WRITE(6,"(2a)")
     +   'Usage:  hycom_profile_isop syn_prof_xtnd.isp nprof archo.txt',
     +   ' [type]'
        CALL EXIT(1)
      ENDIF
*     write(6,*) "arguments read"
C
C     OPEN ALL FILES.
C
      OPEN(UNIT=11, FILE=CFILEI, FORM='FORMATTED', STATUS='OLD',
     +     IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error: can''t open ',TRIM(CFILEI)
        WRITE(6,*) 'ios   = ',ios
        CALL EXIT(3)
      ENDIF
      OPEN(UNIT=21, FILE=CFILEO, FORM='FORMATTED', STATUS='NEW',
     +     IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error: can''t open ',TRIM(CFILEO)
        WRITE(6,*) 'ios   = ',ios
        CALL EXIT(5)
      ENDIF
*     write(6,*) "files opened"
C
C     READ THE ISOP FILE
C
        READ(11,'(a)',IOSTAT=IOS) CLINE  !header
C
      KDM = 0
      DO
        READ(11,'(a)',IOSTAT=IOS) CLINE
        IF     (IOS.NE.0) THEN
          EXIT
        ENDIF
        READ(CLINE,*) kp
        if     (KP.EQ.NPROF) THEN
          KDM = KDM + 1
*         write(6,*) "kdm =",kdm
*         write(6,*) "cline =",trim(cline)
          IF     (ITYPE.EQ.1) THEN
            READ(CLINE,*) kp,plat,plon,z(kdm),si(kdm,3),si(kdm,4)
          ELSEIF (ITYPE.EQ.2) THEN
            READ(CLINE,*) kp,plat,plon,z(kdm),dum,dum,
     &                                        si(kdm,3),si(kdm,4)
          ELSE   !ITYPE.EQ.3
            READ(CLINE,*) kp,plat,plon,z(kdm),dum,dum,dum,dum,
     &                                        si(kdm,3),si(kdm,4)
          ENDIF
        ENDIF
      ENDDO
      IF     (KDM.EQ.0) THEN
        WRITE(6,*) 'Error: profile not in file; nprof =',nprof
        CALL EXIT(6)
      ENDIF
C
      SIGVER = 6  !17-term sigma-2
      CALL SIG_I(SIGVER)
C
      Z(KDM+1) = Z(KDM) + (Z(KDM)-Z(KDM-1))
      P(1) =  0.0
      DO K= 1,KDM
        P( K+1) = 0.5*(Z(K)+Z(K+1))
        SI(K,1) = 0.0
        SI(K,2) = 0.0
        CALL SIG_P(SI(K,3),SI(K,4), SI(K,5))
      ENDDO
C
C     OUTPUT
C
            WRITE(21,'(3a/a,6i7,2f8.2,i7)')
     &        '##   expt    idm    jdm    kdm',
     &        '   iloc   jloc  lonloc  latloc',
     &        ' yrflag',
     &        '##',    0,     0,     0,   KDM,
     &            0,    0,PLON,PLAT,
     &        3
          WRITE(21,'(3a/a,f11.4,f8.2,f8.1,2f9.3,3f8.3,6f8.2)')
     &      '## model-day  srfhgt  surflx',
     &      '     dpbl   dpmixl    tmix    smix   thmix',
     &      '    umix    vmix   ubavg   vbavg  steric  nsterc',
     &      '#',0.0,
     &      0.0,                      ! cm
     &      0.0,                      ! W/m**2
     &      0.0,                      ! m
     &      0.0,                      ! m
     &      0.0,                      ! degC
     &      0.0,                      ! psu
     &      0.0,                      ! SigmaT
     &      0.0,                      ! cm/s
     &      0.0,                      ! cm/s
     &      0.0,                      ! cm/s
     &      0.0,                      ! cm/s
     &      0.0,                      ! cm
     &      0.0                       ! cm
C

        WRITE(CFORMAT,'(a)')
     &    '(3a)'
        WRITE(21,CFORMAT)
     &      '#  k',
     &      '    utot    vtot  p.temp    saln  p.dens',
     &      '    thkns      dpth'
C
          WRITE(CFORMAT,'(a)')
     &      '(i4,2f8.2,3f8.4,f9.3,f10.3)'
C
        DO K= 1,KDM
          THK = P(K+1) - P(K)
          WRITE(21,CFORMAT)
     &      K,(SI(K,KK),KK=1,5),THK,0.5*(P(K)+P(K+1))
          write(6,*) p(k+1)
        ENDDO !k
      CLOSE(21)
      END

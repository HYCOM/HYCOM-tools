      PROGRAM HYCOM_SIGMA
      IMPLICIT NONE
C
C --- BUILD A TABLE HYCOM SIGMA TRANSFORMATION.
C
      INTEGER IVERSN,KDM,NSIGMA
      REAL    DP00,DP00X,DP00F,DS00,DS00X,DS00F,
     &        DP0K(99),DS0K(99)
C
      INTEGER ICOUNT,ID,K,KK
      REAL    D,DP0KF,DS0KF,DPMS,DSMS,DPSIG,DZ,QSIGMA,Z(0:99),ZNS
C
C --- INPUT, FROM STDIN
C
C --- 'iversn' = hycom version number x10 (20 or 21)
C --- 'kdm   ' = number of layers
C --- 'nsigma' = number of sigma  levels (kdm-nsigma always z levels)
C --- 'dp00'   = deep    z-level spacing minimum thickness (m)
C --- 'dp00x'  = deep    z-level spacing maximum thickness (m)
C --- 'dp00f'  = deep    z-level spacing stretching factor (1.0=const.z)
C --- 'ds00'   = shallow z-level spacing minimum thickness (m)
C --- 'ds00x'  = shallow z-level spacing maximum thickness (m)
C --- 'ds00f'  = shallow z-level spacing stretching factor (1.0=const.z)
C
C --- if iversn<21: ds00=ds00x=dp00s and ds00f=1.0.
C
C --- the above describe a system that is isopycnal or:
C ---     z in    deep water, based on dp00,dp00x,dp00f
C ---     z in shallow water, based on ds00,ds00x,ds00f and nsigma
C ---     sigma between them, based on ds00,ds00x,ds00f and nsigma
C
C --- for z-only set nsigma=0 (and ds00,ds00x,ds00f=dp00,dp00x,dp00f)
C --- for sigma-z (shallow-deep) use a very small ds00
C ---  (pure sigma-z also has ds00f=dp00f and ds00x=dp00x*ds00/dp00)
C --- for z-sigma (shallow-deep) use a very large dp00 (not recommended)
C --- for sigma-only set nsigma=kdm, dp00 large, and ds00 small
C
C --- or, in place of 'dp00','dp00x','dp00f','ds00','ds00x','ds00f' specify:
C --- 'dp0k  ' = layer k deep    z-level spacing minimum thickness (m)
C ---              k=1,kdm; dp0k must be zero for k>nhybrd
C --- 'ds0k  ' = layer k shallow z-level spacing minimum thickness (m)
C ---              k=1,nsigma
C
      CALL BLKINI(IVERSN,'iversn')
      IF     (IVERSN.LE.20) THEN
        CALL BLKINI(KDM,   'kdm   ')
        CALL BLKINI(NSIGMA,'nsigma')
        CALL BLKINR(DS00,  'dp00s ','("# ",a6," =",f10.4," m")')
        CALL BLKINR(DP00,  'dp00  ','("# ",a6," =",f10.4," m")')
        CALL BLKINR(DP00X, 'dp00x ','("# ",a6," =",f10.4," m")')
        CALL BLKINR(DP00F, 'dp00f ','("# ",a6," =",f10.4,"  ")')
        DS00  = MAX(DS00,0.01)
        DS00X = DS00
        DS00F = 1.0
      ELSE
        CALL BLKINI(KDM,   'kdm   ')
        CALL BLKINI(NSIGMA,'nsigma')
        CALL BLKINR2(DP00,K,
     &                     'dp00  ','("# ",a6," =",f10.4," m")',
     &                     'dp0k  ','("# ",a6," =",f10.4," m")' )
        IF     (K.EQ.1) THEN
          CALL BLKINR(DP00X, 'dp00x ','("# ",a6," =",f10.4," m")')
          CALL BLKINR(DP00F, 'dp00f ','("# ",a6," =",f10.4,"  ")')
          CALL BLKINR(DS00,  'ds00  ','("# ",a6," =",f10.4," m")')
          CALL BLKINR(DS00X, 'ds00x ','("# ",a6," =",f10.4," m")')
          CALL BLKINR(DS00F, 'ds00f ','("# ",a6," =",f10.4,"  ")')
        ELSE
          DP0K(1) = DP00
          DO K=2,KDM
            CALL BLKINR(DP0K(K), 'dp0k  ','("# ",a6," =",f10.4," m")')
          ENDDO !k
          DO K=1,NSIGMA
            CALL BLKINR(DS0K(K), 'ds0k  ','("# ",a6," =",f10.4," m")')
          ENDDO !k
          DP00    = DP0K(  1)
          DP00X   = DP0K(KDM)
          DP00F   = -1.0       !signal that dp00 is not input
          DS00    = DS0K(  1)
          DS00X   = DS0K(NSIGMA)
          DS00F   = -1.0
        ENDIF
      ENDIF
C
C --- THREE OUTPUT FILES (SIGMA AND NON-SIGMA LAYERS, CROSS-OVERS).
C
      OPEN(UNIT=11,FILE='hycom_sigma_1.tbl',STATUS='UNKNOWN')
      OPEN(UNIT=12,FILE='hycom_sigma_2.tbl',STATUS='UNKNOWN')
      OPEN(UNIT=13,FILE='hycom_sigma_3.tbl',STATUS='UNKNOWN')
C
      WRITE(11,5000) KDM,NSIGMA,DP00,DP00X,DP00F,DS00,DS00X,DS00F
      WRITE(12,5000) KDM,NSIGMA,DP00,DP00X,DP00F,DS00,DS00X,DS00F
      WRITE(13,5000) KDM,NSIGMA,DP00,DP00X,DP00F,DS00,DS00X,DS00F
 5000 FORMAT(
     &     '# kdm=',I3,
     &   ' nsigma=',I3,
     &    '  dp00=',F7.2,
     &    ' dp00x=',F7.2,
     &    ' dp00f=',F6.3 /
     &     '#        ',
     &   '           ',
     &    '  ds00=',F7.2,
     &    ' ds00x=',F7.2,
     &    ' ds00f=',F6.3)
C
      WRITE(13,'(a)') '#    DEPTH     Z-DEPTH'
C
C --- Z-LEVELS
C
      IF     (DP00F.NE.-1.0) THEN
        DS0K(1)=DS00
        DSMS=DS0K(1)
        DS0KF=1.0
        DO K=2,NSIGMA
          DS0KF=DS0KF*DS00F
          DS0K(K)=MIN(DS00*DS0KF,DS00X)
          DSMS=DSMS+DS0K(K)
        ENDDO
        DPSIG =    DSMS
        QSIGMA=1.0/DSMS
      ELSE
        DSMS=DS0K(1)
        DO K=2,NSIGMA
          DSMS=DSMS+DS0K(K)
        ENDDO
        DPSIG =    DSMS
        QSIGMA=1.0/DSMS
      ENDIF
C
      DPMS =0.0
      DSMS =0.0
      DP0KF=1.0/DP00F
      DO K=1,KDM
        IF     (DP00F.NE.-1.0) THEN
          DP0KF=DP0KF*DP00F
          DP0K(K)=MIN(DP00*DP0KF,DP00X)
        ENDIF
        DPMS=DPMS+DP0K(K)
        IF     (K.LE.NSIGMA) THEN
          DSMS=DSMS+DS0K(K)
        ELSE
          DS0K(K)=0.0
        ENDIF
        WRITE(11,6000) K,DP0K(K),DPMS,DS0K(K),DSMS
        WRITE(12,6000) K,DP0K(K),DPMS,DS0K(K),DSMS
 6000   FORMAT('# k =',I3,
     &         '  Z-thk,dpth =',F7.2,' m,',F9.2,' m',
     &         '  S-thk,dpth =',F7.2,' m,',F9.2,' m')
        IF     (K.LE.NSIGMA .AND. DP0K(K).GT.DS0K(K)) THEN
          IF     (DS0K(K).GT.0.0) THEN
C
C           WHERE THIS SIGMA LAYER REACHES DS0K(K) AND BECOMES Z AGAIN
C
            D   = 1.0/QSIGMA
            ZNS = 0.0
            DO KK= 1,K
              ZNS   = ZNS + DS0K(KK)
            ENDDO
            IF     (ZNS.LT.D) THEN
              WRITE(13,'(f10.3,2x,f10.3)') D,ZNS
            ENDIF
          ENDIF
C
C         WHERE THIS LAYER FIRST BECOMES SIGMA
C
          D   = DP0K(K)/(DS0K(K)*QSIGMA)
          ZNS = 0.0
          DO KK= 1,K
            DPSIG = D*DS0K(KK)*QSIGMA
            ZNS   = ZNS + MIN( DP0K(KK), MAX( DS0K(KK), DPSIG ) )
          ENDDO
          IF     (ZNS.LT.D .AND. D.LE.1000.0) THEN
            WRITE(13,'(f10.3,2x,f10.3)') D,ZNS
          ENDIF
        ENDIF
      ENDDO !k
      WRITE(11,'(a,i4,98i7)') '# DEPTH, K=',(K,K=1,NSIGMA),KDM
      WRITE(12,'(a,i4,98i7)') '# DEPTH, K=',(K,K=NSIGMA+1,KDM)
C
C --- TABLE OF INTERFACE DEPTHS.
C
      ICOUNT= 0
      ZNS   = -1.0
      DO ID= 0,1000
        D     = ID
        Z(0)  = 0.0
        DO K= 1,NSIGMA
          DPSIG = D*DS0K(K)*QSIGMA
          DZ    = MIN( DP0K(K), MAX( DS0K(K), DPSIG ) ) 
          Z(K)  = MIN( Z(K-1) + DZ, D )
        ENDDO
        DO K= NSIGMA+1,KDM-1
          Z(K)  = Z(K-1) +      DP0K(K)
          Z(K)  = MIN( Z(K), D )
        ENDDO
        DO K= KDM,KDM
          Z(K)  =            D
        ENDDO
        WRITE(11,'(f9.2,2x,99f9.2)') D,(Z(K),K=1,NSIGMA),Z(KDM)
        WRITE(12,'(f9.2,2x,99f9.2)') D,(Z(K),K=NSIGMA+1,KDM)
        IF     (ZNS.EQ.Z(KDM-1)) THEN
          ICOUNT = ICOUNT + 1
          IF     (ICOUNT.GT.10) THEN
            EXIT
          ENDIF
        ELSE
          ZNS = Z(KDM-1)
        ENDIF
      ENDDO
      END

      subroutine blkinr(rvar,cvar,cfmt)
      implicit none
c
      real      rvar
      character cvar*6,cfmt*(*)
c
c     read in one real value
c
      character*6 cvarin
c
      read(5,*) rvar,cvarin
*     write(6,cfmt) cvarin,rvar
c
      if     (cvar.ne.cvarin) then
        write(6,*) 
        write(6,*) 'error in blkinr - input ',cvarin,
     &                      ' but should be ',cvar
        write(6,*) 
        stop
      endif
      return
      end
      subroutine blkinr2(rvar,nvar,cvar1,cfmt1,cvar2,cfmt2)
      implicit none
c
      real      rvar
      integer   nvar
      character cvar1*6,cvar2*6,cfmt1*(*),cfmt2*(*)
c
c     read in one real value from stdin,
c     identified as either cvar1 (return nvar=1) or cvar2 (return nvar=2)
c
      character*6 cvarin
c
      read(5,*) rvar,cvarin
c
      if     (cvar1.eq.cvarin) then
        nvar = 1
        write(6,cfmt1) cvarin,rvar
      elseif (cvar2.eq.cvarin) then
        nvar = 2
        write(6,cfmt2) cvarin,rvar
      else
        write(6,*) 
        write(6,*) 'error in blkinr2 - input ',cvarin,
     +                      ' but should be ',cvar1,' or ',cvar2
        write(6,*) 
        stop
      endif
      return
      end
      subroutine blkini(ivar,cvar)
      implicit none
c
      integer     ivar
      character*6 cvar
c
c     read in one integer value
c
      character*6 cvarin
c
      read(5,*) ivar,cvarin
*     write(6,6000) cvarin,ivar
c
      if     (cvar.ne.cvarin) then
        write(6,*) 
        write(6,*) 'error in blkini - input ',cvarin,
     &                      ' but should be ',cvar
        write(6,*) 
        stop
      endif
      return
 6000 format('# ',a6,' =',i6)
      end

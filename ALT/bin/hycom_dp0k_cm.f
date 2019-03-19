      PROGRAM HYCOM_DP0K_CM
      IMPLICIT NONE
C
C --- Given a dp00-based z-level setup,
C ---  printout the corresponding dp0k's, rounded to 1 cm.
C
C --- 'kdm   ' = number of layers (-ve read in k01 and z00)
C --- 'k01'    = starting layer, default 1
C --- 'z00'    = starting depth, default 0 (m)
C --- 'dp00'   = z-level spacing minimum thickness (m)
C --- 'dp00x'  = z-level spacing maximum thickness (m)
C --- 'dp00f'  = z-level spacing stretching factor (1.0=const.spacing)
C
      INTEGER K,K01,KDM
      REAL    DP00,DP00F,DPK,DPS,DP00X,Z00
      REAL    CMK,CMO,CMS
C
      READ(5,*) KDM
      IF     (KDM.LT.0) THEN
        KDM = -KDM
        READ(5,*) K01
        READ(5,*) Z00
      ELSE
        K01 = 1
        Z00 = 0.0
      ENDIF
      READ(5,*) DP00
      READ(5,*) DP00X
      READ(5,*) DP00F
C
      DPS = Z00
      CMS = NINT(DPS*100.0)/100.0
      DO K= 1,KDM
        CMO = CMS
        DPK = MIN(DP00X,DP00*DP00F**(K-1))
        DPS = DPS + DPK
        CMS = NINT(DPS*100.0)/100.0
        CMK = CMS-CMO
        WRITE(6,'(F7.2,A,I4,2A,F8.2)')
     &    CMK,
     &    "   'dp0k  ' = layer",
     &    K+K01-1,
     &    " deep    z-level spacing minimum thickness (m)",
     &    "              TOTAL:",
     &    CMS
      ENDDO !k
      END

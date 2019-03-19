      PROGRAM DEPTH
C
C     PRINTOUT HYCOM Z-DEPTHS FOR FIVE DIFFERENT DP00F CHOICES
C
C --- 'dp00f'  = z-level spacing stretching factor (1.0=const.spacing)
C --- 'dp00'   = z-level spacing minimum thickness (m)
C --- 'dp00x'  = z-level spacing maximum thickness (m)
C
      REAL DP00(5),DP00F(5),DPK(5),DPS(5),DP00X
      DATA DPS   /  0.000, 0.000, 0.000, 0.000, 0.000 /
C
      WRITE(6,'(a)',ADVANCE='NO') 'DP00F; DP00 = '
      READ( 5,*) DP00F(1),DP00(1)
      WRITE(6,'(a)',ADVANCE='NO') 'DP00F; DP00 = '
      READ( 5,*) DP00F(2),DP00(2)
      WRITE(6,'(a)',ADVANCE='NO') 'DP00F; DP00 = '
      READ( 5,*) DP00F(3),DP00(3)
      WRITE(6,'(a)',ADVANCE='NO') 'DP00F; DP00 = '
      READ( 5,*) DP00F(4),DP00(4)
      WRITE(6,'(a)',ADVANCE='NO') 'DP00F; DP00 = '
      READ( 5,*) DP00F(5),DP00(5)
C
      WRITE(6,'(a)',ADVANCE='NO') '            DP00X = '
      READ( 5,*) DP00X
C
      WRITE(6,*)
      WRITE(6,'(A,5(A,F6.4))')
     +  '# K  ',
     +  ' DP00F=',DP00F(1),
     + '  DP00F=',DP00F(2),
     + '  DP00F=',DP00F(3),
     + '  DP00F=',DP00F(4),
     + '  DP00F=',DP00F(5)
      DO K= 1,40
        DO L= 1,5
          DPK(L) = MIN(DP00X,DP00(L)*DP00F(L)**(K-1))
          DPS(L) = DPS(L) + DPK(L)
        ENDDO
        WRITE(6,'(I3,5(2F7.1))') K,(DPK(L),DPS(L),L=1,5)
      ENDDO
      END

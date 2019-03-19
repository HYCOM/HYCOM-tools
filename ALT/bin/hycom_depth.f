      PROGRAM DEPTH
C
C     PRINTOUT HYCOM Z-DEPTHS FOR FIVE DIFFERENT DP00F CHOICES
C     THESE DP00F VALUES WERE CHOSEN BECAUSE THEIR 6-TH POWERS
C     ARE (APPROXIMATELY) ROUND NUMBERS:
C     1.07**6 = 1.5,  1.125**6 = 2,  1.2**6 = 3,  1.26**6 = 4.
C
C --- 'dp00f'  = z-level spacing stretching factor (1.0=const.spacing)
C --- 'dp00'   = z-level spacing minimum thickness (m)
C --- 'dp00x'  = z-level spacing maximum thickness (m)
C
      REAL DP00(5),DP00F(5),DPK(5),DPS(5),DP00X
      DATA DP00F /  1.000, 1.070, 1.125, 1.200, 1.260 /
      DATA DPS   /  0.000, 0.000, 0.000, 0.000, 0.000 /
C
      WRITE(6,'(a)',ADVANCE='NO') 'DP00F=1.000; DP00 = '
      READ( 5,*) DP00(1)
      WRITE(6,'(a)',ADVANCE='NO') 'DP00F=1.070; DP00 = '
      READ( 5,*) DP00(2)
      WRITE(6,'(a)',ADVANCE='NO') 'DP00F=1.125; DP00 = '
      READ( 5,*) DP00(3)
      WRITE(6,'(a)',ADVANCE='NO') 'DP00F=1.200; DP00 = '
      READ( 5,*) DP00(4)
      WRITE(6,'(a)',ADVANCE='NO') 'DP00F=1.260; DP00 = '
      READ( 5,*) DP00(5)
C
      WRITE(6,'(a)',ADVANCE='NO') '            DP00X = '
      READ( 5,*) DP00X
C
      WRITE(6,*)
      WRITE(6,'(6A)')
     + '# K  ',
     + ' DP00F=1.000  ',
     + ' DP00F=1.070  ',
     + ' DP00F=1.125  ',
     + ' DP00F=1.200  ',
     + ' DP00F=1.260  '
      DO K= 1,25
        DO L= 1,5
          DPK(L) = MIN(DP00X,DP00(L)*DP00F(L)**(K-1))
          DPS(L) = DPS(L) + DPK(L)
        ENDDO
        WRITE(6,'(I3,5(2F7.1))') K,(DPK(L),DPS(L),L=1,5)
      ENDDO
      END

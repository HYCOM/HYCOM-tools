      PROGRAM TOP_DIFF_SUB_PPMX
      use mod_za  ! HYCOM array I/O interface
      implicit none
c
      INTEGER*1     RAWPAL(0:255,3),PPMPAL(3,0:255)
      INTEGER*8     IH,JH,K,NC
      CHARACTER*5   CIH,CJH
      CHARACTER*1   CNWLIN
c
      real      hmaxa,hmaxb,hmina,hminb
      character preambl(5)*79,cline*80
c
      logical   lexist
      integer   i,j
      integer*8 i1,j1,idm8,jdm8,ksub
c
c --- This program reads in two standard HYCOM depth files, 
c --- subregions them and writes out a 24-bit PPM image file
c --- of the subregion with the differences marked.
c
c --- the subsample factor, ksub, is read from stdin
c --- the subregion, i1,j1,idms,jdms, is read from stdin
c
      character*1, allocatable :: ib(:)
      integer,     allocatable :: iop(:,:)
      real,        allocatable :: depths_in(:,:)
      real,        allocatable :: depths_a(:,:),depths_b(:,:)
c
      read(5,*) ksub,i1,j1,idm8,jdm8
      write(6,'(a, i8)') 'ksub      =',ksub
      write(6,'(a,2i8)') 'i1,  j1   =',i1,  j1
      write(6,'(a,2i8)') 'idms,jdms =',idm8,jdm8
      write(6,'(a,2i8)') 'iend,jend =',i1+idm8-1,j1+jdm8-1
      call zhflsh(6)
c
      call xcspmd  !input idm,jdm
      call zaiost
c
      if     (i1.lt.1 .or. i1+idm8-1.gt.idm .or.
     &        j1.lt.1 .or. j1+jdm8-1.gt.jdm     ) then
        write(6,'(/ a / a,2i8 /)')
     &    'error - subregion outside full region',
     &    '        idm,jdm = ',idm,jdm
        call zhflsh(6)
        stop
      endif
c
      allocate( ib( 19+3*((idm8+ksub-1)/ksub)*((jdm8+ksub-1)/ksub)) )
      allocate(       iop(idm,jdm)  )
      allocate( depths_in(idm,jdm)  )
      allocate( depths_a(idm8,jdm8) )
      allocate( depths_b(idm8,jdm8) )
      ib(:) = char(0)
c
c --- acquire 1st basin depths from unit 51.
c
      call zhopen(51, 'formatted', 'old', 0)
      read (51,'(a79)') preambl
      read (51,'(a)')   cline
      close(unit=51)
      write(6,'(/(a))') preambl,cline(1:len_trim(cline))
c
      i = index(cline,'=')
      read (cline(i+1:),*)   hminb,hmaxb
c
      call zaiopn('old', 51)
      call zaiord(depths_in,iop,.false., hmina,hmaxa, 51)
      call zaiocl(51)
c
      if     (abs(hmina-hminb).gt.abs(hminb)*1.e-4 .or.
     &        abs(hmaxa-hmaxb).gt.abs(hmaxb)*1.e-4     ) then
        write(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - .a and .b 1st topography files not consistent:',
     &    '.a,.b min = ',hmina,hminb,hmina-hminb,
     &    '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
        call zhflsh(6)
        stop
      endif
c
      depths_a(:,:) = depths_in(i1:i1+idm8-1,j1:j1+jdm8-1)
c
      call zhopen(61, 'formatted', 'old', 0)
      read (61,'(a79)') preambl
      read (61,'(a)')   cline
      close(unit=61)
      write(6,'(/(a))') preambl,cline(1:len_trim(cline))
c
      i = index(cline,'=')
      read (cline(i+1:),*)   hminb,hmaxb
c
      call zaiopn('old', 61)
      call zaiord(depths_in,iop,.false., hmina,hmaxa, 61)
      call zaiocl(61)
c
      if     (abs(hmina-hminb).gt.abs(hminb)*1.e-4 .or.
     &        abs(hmaxa-hmaxb).gt.abs(hmaxb)*1.e-4     ) then
        write(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - .a and .b 2nd topography files not consistent:',
     &    '.a,.b min = ',hmina,hminb,hmina-hminb,
     &    '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
        call zhflsh(6)
        stop
      endif
c
      depths_b(:,:) = depths_in(i1:i1+idm8-1,j1:j1+jdm8-1)
c
c --- zero land
c
      do j= 1,jdm8
        do i= 1,idm8
          if     (depths_a(i,j).gt.2.0**99) then
            depths_a(i,j) = 0.0
          endif
          if     (depths_b(i,j).gt.2.0**99) then
            depths_b(i,j) = 0.0
          endif
        enddo
      enddo
C
C     READ IN PSUDO-COLOR PALLETTE.
C
      OPEN(UNIT=41,FILE='xbathy.pal',STATUS='OLD',
     +     FORM='UNFORMATTED',ACCESS='DIRECT',RECL=768)
      READ( 41,REC=1) RAWPAL
      CLOSE(41)
      DO 110 K= 1,3
        DO 112 I= 0,255
          PPMPAL(K,I) = RAWPAL(I,K)
  112   CONTINUE
  110 CONTINUE
C
C     CONVERT TO IMAGE.
C
      IF     (MAX((IDM8+KSUB-1)/KSUB,(JDM8+KSUB-1)/KSUB).LE. 999) THEN
        NC = 15 + 3*((IDM8+KSUB-1)/KSUB)*((JDM8+KSUB-1)/KSUB)
      ELSEIF (MAX((IDM8+KSUB-1)/KSUB,(JDM8+KSUB-1)/KSUB).LE.9999) THEN
        NC = 17 + 3*((IDM8+KSUB-1)/KSUB)*((JDM8+KSUB-1)/KSUB)
      ELSE
        NC = 19 + 3*((IDM8+KSUB-1)/KSUB)*((JDM8+KSUB-1)/KSUB)
      ENDIF
C
      WRITE(6,*) 
      WRITE(6,*) "NC = ",NC,NC-3*((IDM8+KSUB-1)/KSUB)*
     &                           ((JDM8+KSUB-1)/KSUB)+1
      CALL ZHFLSH(6)
C
      CALL TOPIMG(IB(NC-3*((IDM8+KSUB-1)/KSUB)*
     &                    ((JDM8+KSUB-1)/KSUB)+1),char(PPMPAL),
     &            DEPTHS_A,DEPTHS_B,IDM8,JDM8, KSUB)
C
C     OUTPUT AS PPM IMAGE FILE.
C
      IH = (IDM8+KSUB-1)/KSUB
      JH = (JDM8+KSUB-1)/KSUB
C
      CNWLIN = CHAR(10)
      WRITE(CIH,4000) IH
      WRITE(CJH,4000) JH
C
      IB(1) = 'P'
      IB(2) = '6'
      IB(3) = CNWLIN
      IF     (MAX(IH,JH).LE. 999) THEN
        IB( 4) = CIH(3:3)
        IB( 5) = CIH(4:4)

        IB( 6) = CIH(5:5)
        IB( 7) = ' '
        IB( 8) = CJH(3:3)
        IB( 9) = CJH(4:4)
        IB(10) = CJH(5:5)
        IB(11) = CNWLIN
        IB(12) = '2'
        IB(13) = '5'
        IB(14) = '5'
        IB(15) = CNWLIN
      ELSEIF (MAX(IH,JH).LE.9999) THEN
        IB( 4) = CIH(2:2)
        IB( 5) = CIH(3:3)
        IB( 6) = CIH(4:4)
        IB( 7) = CIH(5:5)
        IB( 8) = ' '
        IB( 9) = CJH(2:2)
        IB(10) = CJH(3:3)
        IB(11) = CJH(4:4)
        IB(12) = CJH(5:5)
        IB(13) = CNWLIN
        IB(14) = '2'
        IB(15) = '5'
        IB(16) = '5'
        IB(17) = CNWLIN
      ELSE
        IB( 4) = CIH(1:1)
        IB( 5) = CIH(2:2)
        IB( 6) = CIH(3:3)
        IB( 7) = CIH(4:4)
        IB( 8) = CIH(5:5)
        IB( 9) = ' '
        IB(10) = CJH(1:1)
        IB(11) = CJH(2:2)
        IB(12) = CJH(3:3)
        IB(13) = CJH(4:4)
        IB(14) = CJH(5:5)
        IB(15) = CNWLIN
        IB(16) = '2'
        IB(17) = '5'
        IB(18) = '5'
        IB(19) = CNWLIN
      ENDIF
C
      WRITE(6,*) 'WRITING IMAGE TO fort.31'
      CALL ZHFLSH(6)
      OPEN(UNIT=31,FILE='fort.31',STATUS='NEW',
     +     FORM='UNFORMATTED',ACCESS='DIRECT',RECL=NC)
      WRITE(31,REC=1) IB(1:NC)
      CLOSE(31)
C
 4000 FORMAT( I5 )
C     END OF TOPPPM.
      END
      SUBROUTINE TOPIMG(IB,IP,HD1,HD2,IH,JH,KSUB)
      implicit none
C
      INTEGER*8     IH,JH,KSUB
      REAL          HD1(IH,JH),HD2(IH,JH)
      CHARACTER*(1) IB(3,(IH+KSUB-1)/KSUB,(JH+KSUB-1)/KSUB),IP(3,0:255)
C
C**********
C*
C 1)  CONVERT TOPOGRAPHY TO AN 8-BIT IMAGE.
C*
C**********
C
      LOGICAL   LLAND
      INTEGER   I,II,J,JJ,JB,K,M1,M2
      INTEGER*8 N0,N1,N2,N5
      REAL      D2
C
      WRITE(6,*) 
      WRITE(6,*) 'IH,JH,KSUB = ',IH,JH,KSUB
      WRITE(6,*) 
      WRITE(6,*) 'IP(:,47) = ',(ICHAR(IP(K,47)),K=1,3)
      WRITE(6,*) 'IP(:,48) = ',(ICHAR(IP(K,48)),K=1,3)
      WRITE(6,*) 'IP(:,49) = ',(ICHAR(IP(K,49)),K=1,3)
      WRITE(6,*) 'IP(:,50) = ',(ICHAR(IP(K,50)),K=1,3)
      WRITE(6,*) 
      CALL ZHFLSH(6)
      JB = (JH+KSUB-1)/KSUB
      N0 = 0
      N1 = 0
      N2 = 0
      N5 = 0
      D2 =  1.0/50.0
      WRITE(6,*) 'JB,D2 = ',JB,D2
      WRITE(6,*) 'HD.1,1 = ',HD1(1,1),HD2(1,1)
      WRITE(6,*) 
      CALL ZHFLSH(6)
      DO J= 1,JH,KSUB
        DO I= 1,IH,KSUB
          M1 = 0
          M2 = 0
          DO JJ= MAX(1,J-KSUB/2),MIN(JH,J+KSUB/2)
            DO II= MAX(1,I-KSUB/2),MIN(IH,I+KSUB/2)
              IF (HD1(II,JJ).LE.0.0) THEN
                M1 = M1 + 1
              ENDIF
              IF (HD2(II,JJ).LE.0.0) THEN
                M2 = M2 + 1
              ENDIF
            ENDDO
          ENDDO
          IF     (M1.GT.0 .AND. M1.EQ.M2) THEN  !same land (blue)
            IB(1,(I+KSUB-1)/KSUB,JB+1-(J+KSUB-1)/KSUB) = IP(1,48)
            IB(2,(I+KSUB-1)/KSUB,JB+1-(J+KSUB-1)/KSUB) = IP(2,48)
            IB(3,(I+KSUB-1)/KSUB,JB+1-(J+KSUB-1)/KSUB) = IP(3,48)
            N0 = N0 + 1
          ELSEIF (M1.GT.M2) THEN  !M1 more land (cyan)
            IB(1,(I+KSUB-1)/KSUB,JB+1-(J+KSUB-1)/KSUB) = IP(1,49)
            IB(2,(I+KSUB-1)/KSUB,JB+1-(J+KSUB-1)/KSUB) = IP(2,49)
            IB(3,(I+KSUB-1)/KSUB,JB+1-(J+KSUB-1)/KSUB) = IP(3,49)
            N1 = N1 + 1
            IF     (KSUB.EQ.1 .AND. MOD(N1,MIN(IH,JH)).EQ.1) THEN
              WRITE(6,'(a,2i8,f10.2)') 'LAND 1:',I,J,HD2(I,J)
            ENDIF
          ELSEIF (M2.GT.M1) THEN  !M2 more land (Magenta)
            IB(1,(I+KSUB-1)/KSUB,JB+1-(J+KSUB-1)/KSUB) = IP(1,50)
            IB(2,(I+KSUB-1)/KSUB,JB+1-(J+KSUB-1)/KSUB) = IP(2,50)
            IB(3,(I+KSUB-1)/KSUB,JB+1-(J+KSUB-1)/KSUB) = IP(3,50)
            N2 = N2 + 1
            IF     (KSUB.EQ.1 .AND. MOD(N2,MIN(IH,JH)).EQ.1) THEN
              WRITE(6,'(a,2i8,f10.2)') 'LAND 2:',I,J,HD1(I,J)
            ENDIF
          ELSE  !both ocean (grayscale)
            II = MIN(248, 100 + NINT((HD1(I,J)      )*D2) )
            IB(1,(I+KSUB-1)/KSUB,JB+1-(J+KSUB-1)/KSUB) = IP(1,II)
            IB(2,(I+KSUB-1)/KSUB,JB+1-(J+KSUB-1)/KSUB) = IP(2,II)
            IB(3,(I+KSUB-1)/KSUB,JB+1-(J+KSUB-1)/KSUB) = IP(3,II)
            N5 = N5 + 1
          ENDIF
        ENDDO
      ENDDO
      WRITE(6,*)
      WRITE(6,*) 'TOPIMG STATISTICS:'
      WRITE(6,*) '       LAND 1+2 POINTS = ',N0
      WRITE(6,*) '       LAND 1   POINTS = ',N1
      WRITE(6,*) '       LAND   2 POINTS = ',N2
      WRITE(6,*) '            SEA POINTS = ',N5
      WRITE(6,*)
      CALL ZHFLSH(6)
      RETURN
C     END OF TOPIMG.
      END

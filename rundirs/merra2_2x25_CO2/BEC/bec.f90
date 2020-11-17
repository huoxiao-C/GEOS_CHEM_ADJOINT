      program main

      integer i, j, k, M, mpb, mpe  
      character(len=255) str
      CHARACTER(LEN=255)                  :: DATE, TIME, ZONE
      REAL*8                              :: CO2_ORI(144, 91, 47, 1), CO2_INC(144, 91, 47, 1)
      REAL*8                              :: X(144, 91, 47), BinvX(144, 91, 47)
      REAL*8                              :: Sx(144, 144), Sy(91, 91), Sz(47, 47)
      INTEGER                             :: VALUES(8)
      LOGICAL                             :: ALIVE, DFLAG
      INTEGER                             :: L, FLAG
!variable parameters
      !time length
      call date_and_time( date , time, zone, values )
      print*, 'PROGRAM START TIME ', date(1:4)//'-'//date(5:6)// &
     &        '-'//date(7:8)//'  '//                             &
     &         time(1:2)//':'//time(3:4)//':'//time(5:10)

     
     DO WHILE(.TRUE.)
        INQUIRE( FILE = 'EXIT_FLAG', EXIST = ALIVE )
        IF ( ALIVE ) THEN
           OPEN(UNIT=1995, FILE='EXIT_FLAG')
           READ(1995, *) FLAG
           CLOSE(1995)
           IF (FLAG==1) THEN
              EXIT
           ENDIF
        ENDIF
        INQUIRE( FILE = 'gctm.gc.01', EXIST = ALIVE )
        IF ( .NOT. ALIVE  ) THEN 
           CYCLE
        ELSE 
           OPEN(UNIT=1996, FILE='gctm.gc.01')
           READ(1996, *) CO2_ORI
           CLOSE(1996)
        ENDIF
        INQUIRE( FILE = 'gctm.gc.cur', EXIST = ALIVE )
        IF ( .NOT. ALIVE ) THEN
           CYCLE
        ELSE
           OPEN(UNIT=1997, FILE='gctm.gc.cur')
           READ(1997, *) CO2_INC
           CLOSE(1997)
        ENDIF
        EXIT
     ENDDO
!     print*, 'huoxiao_debug' 
     ! read co2
!     INQUIRE( FILE = 'gctm.gc.01', EXIST = ALIVE )
!     IF ( ALIVE ) THEN
!     OPEN(UNIT=1996, FILE='gctm.gc.01')
!     READ(1996, *) CO2_ORI
!     CLOSE(1996)
!     ENDIF

!     INQUIRE( FILE = 'gctm.gc.cur', EXIST = ALIVE )
!     IF ( ALIVE ) THEN
!     OPEN(UNIT=1997, FILE='gctm.gc.cur')
!     READ(1997, *) CO2_INC
!     CLOSE(1997)
!     ENDIF
     !read Sx, Sy, Sz
     OPEN(UNIT=1998, FILE='Sx')
     OPEN(UNIT=1999, FILE='Sy')
     OPEN(UNIT=2000, FILE='Sz')
     READ(1998, *) Sx
     READ(1999, *) Sy
     READ(2000, *) Sz
     CLOSE(1998)
     CLOSE(1999)
     CLOSE(2000)
     X(:, :, :) = CO2_INC(:, :, :, 1)-CO2_ORI(:, :, :, 1) 
     print*, 'huoxiao_debug maxval X', MAXVAL(X)   
     mpb = 1
     mpe = 30

     BinvX = 0.0
!$OMP PARALLEL DO &
!$OMP DEFAULT( SHARED ) &
!$OMP PRIVATE( I, J, L) 
     DO L = 1, 47
           DO J = 1, 91
                    DO I = mpb, mpe
                       BinvX(I, J, L) = CALC_ROW(I, J, L, Sx, Sy, Sz, X)   
!                       IF (BinvX(I, J, L)>10.0) THEN
!                          print*, I, J, L
!                          STOP
!                       ENDIF
     ENDDO
     ENDDO
     ENDDO
!$OMP END PARALLEL DO
      OPEN(UNIT=2001, FILE='BINVX_01')
      WRITE(2001, *) BinvX
      CLOSE(2001)
      OPEN(UNIT=2002, FILE='FLAG_01')
      WRITE(2002, *) 1
      CLOSE(2002)

      call date_and_time( date , time, zone, values )
      print*, 'PROGRAM START TIME ', date(1:4)//'-'//date(5:6)//  &
     &        '-'//date(7:8)//'  '//                              &
     &         time(1:2)//':'//time(3:4)//':'//time(5:10)

     CONTAINS
     FUNCTION CALC_ROW(LON, LAT, LEV, Sx, Sy, Sz, X) RESULT(VALUES)

     INTEGER, INTENT(IN) :: LON
     INTEGER, INTENT(IN) :: LAT
     INTEGER, INTENT(IN) :: LEV
      REAL*8, INTENT(IN) :: Sx(144, 144)
      REAL*8, INTENT(IN) :: Sy(91, 91)
      REAL*8, INTENT(IN) :: Sz(47, 47)
      REAL*8, INTENT(IN) :: X(144, 91, 47)
      REAL*8             :: VALUES
    
     INTEGER             :: I, J, L
     
     VALUES = 0.0
     DO L =1, 47
        DO J = 1, 91
           DO I = 1, 144
              VALUES = VALUES+Sx(LON, I)*Sy(LAT, J)*Sz(LEV, L)*X(I, J, L)
           ENDDO
        ENDDO
     ENDDO 
    END FUNCTION CALC_ROW
end program

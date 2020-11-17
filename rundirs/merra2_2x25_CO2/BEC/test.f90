program huo
      integer Lon, Lat,Lev, mpb, mpe
      real*8, pointer        :: BinvX(:, :, :, :), BinvX1(:, :, :, :), &
       TEMP(:, :, :, :)   
      real*8                :: Sx(144, 144), Sy(91, 91), Sz(47, 47), MAXV
      Real*8        :: X(144, 91, 47), VALUES1, MM(2, 2), NN(2,2)
      CHARACTER(LEN=255)                  :: DATE, TIME, ZONE
      INTEGER                             :: VALUES(8), rc

 
      ALLOCATE(BinvX(144, 91, 47, 1), stat=rc)
     ALLOCATE(BinvX1(144, 91, 47, 1), stat=rc)
      ALLOCATE(temp(144, 91, 47, 1), stat=rc)
      OPEN(UNIT=1995, FILE='../gctm.gc.01')
      OPEN(UNIT=1996, FILE='../gctm.gc.02')
       READ(1995, *) BINVX
       READ(1996,*) BINVX1
      CLOSE(1995)
      CLOSE(1996)
       TEMP = BINVX-BINVX1
       MAXV = 0.0
       MM=1
       NN=2
       print*, MAXVAL(ABS(TEMP)), MM+NN
       DO LON = 1,144
          DO LAT = 1, 91
             DO LEV = 1, 47
                IF (abs(BINVX1(LON,LAT,LEV,1)-BINVX(LON,LAT,LEV,1)) >maxv) THEN
                  maxv = abs(BINVX1(LON,LAT,LEV,1)-BINVX(LON,LAT,LEV,1))
                 ENDIF
              ENDDO
           ENDDO
       ENDDO
       PRINT*, 'MAXVAL', MAXV
!      BINVX = BINVX-BINVX1
!      print*, MAXVAL(BINVX(:, :, :, 1))
!     mpb = 1
!     mpe = 30
!     Sx = 0.0 
!     Sy = 0.0
!     Sz = 0.0
!     X = 0.0
!      call date_and_time( date , time, zone, values )
!      print*, 'PROGRAM START TIME ', date(1:4)//'-'//date(5:6)// &
!     &        '-'//date(7:8)//'  '//                             &
!     &         time(1:2)//':'//time(3:4)//':'//time(5:10)

!!$OMP PARALLEL DO &
!!$OMP DEFAULT( SHARED ) &
!!$OMP PRIVATE( I, J, L, Lon, Lat, Lev) 
!     DO Lon = 1, 10
!           DO Lat = 1, 91
!                    DO Lev = 1, 47
!                            VALUES1 = 0.0
!                            DO L =1, 47
!                               DO J = 1, 91
!                                  DO I = 1, 144
!                                     VALUES1 = VALUES1
!                                  ENDDO
!                               ENDDO
!                            ENDDO
!                       BinvX(Lon, lat, lev) = VALUES1 
    
!     ENDDO
!     ENDDO
!     ENDDO
!!$OMP END PARALLEL DO
!      call date_and_time( date , time, zone, values )
!      print*, 'PROGRAM START TIME ', date(1:4)//'-'//date(5:6)// &
!     &        '-'//date(7:8)//'  '//                             &
!     &         time(1:2)//':'//time(3:4)//':'//time(5:10)

!     CONTAINS
!     FUNCTION CALC_ROW(LON, LAT, LEV, Sx, Sy, Sz, X) RESULT(VALUES)

!     INTEGER, INTENT(IN) :: LON
!     INTEGER, INTENT(IN) :: LAT
!     INTEGER, INTENT(IN) :: LEV
!      REAL*8, INTENT(IN) :: Sx(144, 144)
!      REAL*8, INTENT(IN) :: Sy(91, 91)
!      REAL*8, INTENT(IN) :: Sz(47, 47)
!      REAL*8, INTENT(IN) :: X(144, 91, 47)
!      REAL*8             :: VALUES

!     INTEGER             :: I, J, L

!     VALUES = 0.0
!     DO L =1, 47
!        DO J = 1, 91
!           DO I = 1, 144
!              VALUES = VALUES+Sx(LON, I)*Sy(LAT, J)*Sz(LEV, L)*X(I, J, L)
!           ENDDO
!        ENDDO
!     ENDDO
!    END FUNCTION CALC_ROW
end

      Module BEC_MOD

      USE PRECISION_MOD

      Integer,       Parameter   :: NVNum = 15
      CONTAINS

      FUNCTION  Calc_Val (Vo, Va, Vl, Vnum, Eonum, Eanum,&
                          Vec, x, y, z, BIndex_Array, Debug &
                          ) RESULT( RValue )

      Integer,       INTENT(IN)  :: Vnum
      Integer,       INTENT(IN)  :: Eonum, Eanum
      Integer,       INTENT(IN)  :: BIndex_Array(Vnum)
      Integer,       INTENT(IN)  :: x, y, z
      Real(fp),      INTENT(IN)  :: Vo(Vnum), Va(Vnum), Vl(z)
      Real(fp),      INTENT(IN)  :: Vec(x, y, z)
      Real(fp)                   :: RValue
      Integer,       INTENT(IN)  :: DEBUG
      !local variables
      Integer                    :: I, J, K, CIndex
      Integer                    :: Cx, Cy, Cz
      Real(fp)                   :: Temp
     
    
      RValue = 0.0_8      
      DO I = 1, Eonum
         Cx = INT(BIndex_Array(I)/(y*z))
         Cy = INT((BIndex_Array(I)-Cx*y*z)/z)
         Cz = 1
         Cx = Cx+1
         Cy = Cy+1
         DO J = 1, Eanum
            Temp = 0.0_8
            DO K = 1, z
!                IF ( DEBUG ==1 ) THEN
!                   print*, 'Cx, Cy, Cz',Cx, Cy, K, Vo(I), Va(J), Vl(k)
!                ENDIF
                Temp = Temp+Vo(I)*Va(J)*Vl(K)*Vec(Cx, Cy, K)
            ENDDO
            Cx = Cx+INT((Cy)/y)
            Cy = Cy+1
            Cy = Mod(Cy, y)
            IF ( Cy ==0 ) Cy = y
            RValue = RValue+Temp
         ENDDO
      ENDDO

      END FUNCTION Calc_Val
      
        
      FUNCTION  UTw( Sx, Sy, Sz, D, w, x, y, z ) RESULT( Vec )
      ! **********************************************
      ! Background error covariance multiply w
      ! Reference:
      ! A three-dimensional variational data assimilation system for 
      ! multiple aerosol species with WRF/Chem and an application to PM2.5 prediction
      ! *************************************************
      ! x: lon, y: lat, z: vertical layer
      ! Sx: square root of covariance matrix in lon direction
      ! Sy: square root of covariance matrix in lat direction
      ! Sz: square root of covariance matrix in vertical direction
      ! D: diagonal matrix
      ! calculate D*square_root(C)*w in this subroutine, here C is lower multiangular matrix
      Integer,       INTENT(IN)  :: x
      Integer,       INTENT(IN)  :: y
      Integer,       INTENT(IN)  :: z
      Real(fp),      INTENT(IN)  :: Sx(x, y)
      Real(fp),      INTENT(IN)  :: Sy(y, y)
      Real(fp),      INTENT(IN)  :: Sz(z, z)
      Real(fp),      INTENT(IN)  :: D(x,y,z)
      Real(fp),      INTENT(IN)  :: w(x, y, z)
      Real(fp)                   :: Vec(x, y, z)
 
      !local variable
      Real(fp)                   :: Vo(NVNum), Va(NVNum), Vl(z)
      Integer                    :: BIndex_Array(NVNum), EIdenx
      Integer                    :: Vnum
      Integer                    :: Lon, Lat, Lev, TIndex
      Integer                    :: Eonum, Eanum, Bonum, Banum
      Integer                    :: I
      !initialization
      Vo = 0.0_8
      Va = 0.0_8
      Vl = 0.0_8
      BIndex_Array = 0.0_8
      !complexity O(x*y*z*z)
      DO Lon = 1, x
         DO Lat = 1, y
            DO Lev = 1, z
               IF ( Lat <= NVNum ) THEN
                  DO TIndex = 1, Lat
                     Va(TIndex) = Sy(Lat, TIndex)
                  ENDDO
                  Eanum = Lat
                  Banum = 0
               ELSE
                  DO TIndex = 1, NVNum
                     Va(TIndex) = Sy(Lat, Lat+TIndex-NVNum)
                  ENDDO
                  Eanum = NVNum
                  Banum = Lat-NVNum
               ENDIF

               IF ( Lon <= NVNum) THEN 
                  DO TIndex = 1, Lon
                     Vo(TIndex) = Sx(Lon, TIndex)
                  ENDDO
                  Eonum = Lon
                  Bonum = 0
               ELSE
                  DO TIndex = 1, NVNum
                     Vo(TIndex) = Sx(Lon, Lon+TINdex-NVNum)
                  ENDDO
                  Eonum = NVNum
                  Bonum = Lon-NVNum
               ENDIF
               Vl = Sz(Lev, :)      
               DO I = 1, Eonum
                  BIndex_Array(I) = (Bonum+I-1)*y*z+Banum*z+1
               ENDDO         
               Vec(Lon, Lat, Lev) = Calc_Val(Vo, Va, Vl, NVNum,  &
                                      Eonum, Eanum, w, x, y, z, &
                                      BIndex_Array, 0)*D(Lon, Lat, Lev)
             ENDDO !z
         ENDDO !y
      ENDDO !x
      END FUNCTION UTw
!--------------------------------------------------------------

      FUNCTION  Uw( Sx, Sy, Sz, D, w, x, y, z ) RESULT( Vec )
      ! **********************************************
      ! Background error covariance multiply w
      ! Reference:
      ! A three-dimensional variational data assimilation system for
      ! multiple aerosol species with WRF/Chem and an application to PM2.5 prediction
      ! *************************************************
      ! x: lon, y: lat, z: vertical layer
      ! Sx: square root of covariance matrix in lon direction
      ! Sy: square root of covariance matrix in lat direction
      ! Sz: square root of covariance matrix in vertical direction
      ! D: diagonal matrix
      ! calculate D*square_root(C)*w in this subroutine, here C is upper multiangular matrix 
      Integer,       INTENT(IN)  :: x
      Integer,       INTENT(IN)  :: y
      Integer,       INTENT(IN)  :: z
      Real(fp),      INTENT(IN)  :: Sx(x, y)
      Real(fp),      INTENT(IN)  :: Sy(y, y)
      Real(fp),      INTENT(IN)  :: Sz(z, z)
      Real(fp),      INTENT(IN)  :: D(x,y,z)
      Real(fp),      INTENT(IN)  :: w(z, y, z)
      Real(fp)                   :: Vec(x, y, z)

      !local variable
      Real(fp)                   :: Vo(NVNum), Va(NVNum), Vl(z)
      Integer                    :: BIndex_Array(NVNum), EIdenx
      Integer                    :: Vnum
      Integer                    :: Lon, Lat, Lev, TIndex
      Integer                    :: Eonum, Eanum, Bonum, Banum
      Integer                    :: I

      !initialization
      Vo = 0.0_8
      Va = 0.0_8
      Vl = 0.0_8
      BIndex_Array = 0.0_8
      !complexity O(x*y*z*z)
      DO Lon = 1, x
         DO Lat = 1, y
            DO Lev = 1, z
                IF ( Lon >= x-NVNum+1 ) THEN
                  DO TIndex = 1, x-Lon+1
                     Vo(TIndex) = Sx(Lon, Lon+TIndex-1)
                  ENDDO
                  Eonum = x-Lon+1
                ELSE
                  DO TIndex = 1, NVNum
                     Vo(TIndex) = Sx(Lon, Lon+TIndex-1)
                  ENDDO
                  Eonum = NVNum
                ENDIF
                Bonum = Lon-1

                !lat direction
                IF ( Lat >= y-NVNum+1 ) THEN
                  DO TIndex = 1, y-lat+1
                     Va(TIndex) = Sy(Lat, Lat+TIndex-1)
                  ENDDO
                  Eanum = y-Lat+1
                ELSE
                  DO TIndex = 1, NVNum
                     Va(TIndex) = Sy(Lat, Lat+TIndex-1)
                  ENDDO
                  Eanum = NVNum
                ENDIF
                Banum = Lat-1

                Vl = Sz(Lev, :)
                DO I = 1, Eonum
                   BIndex_Array(I) = (Bonum+I-1)*y*z+Banum*z+1
                ENDDO
                Vec(Lon, Lat, Lev) = Calc_Val(Vo, Va, Vl, NVNum,     &
                                      Eonum, Eanum, w, x, y, z, &
                                      BIndex_Array, 0)*D(Lon, Lat, Lev)
!                if (Lon ==20 .AND. Lat == 20 .AND. Lev == 1) THEN
!                    Vec(Lon, Lat, Lev) = Calc_Val(Vo, Va, Vl, NVNum,     &
!                                      Eonum, Eanum, w, x, y, z, &
!                                      BIndex_Array, 1)*D(Lon, Lat, Lev)
!                   print*, 'Vo', Vo
!                   print*, 'Va', Va
!                   print*,  'Eo, Ea, D', Eonum, Eanum, D(Lon, Lat, Lev)
!                   print*,  'BIndex_Array', BIndex_Array
!                ENDIF
             ENDDO !z
         ENDDO !y
      ENDDO !x
      END FUNCTION Uw
      
      END Module BEC_MOD

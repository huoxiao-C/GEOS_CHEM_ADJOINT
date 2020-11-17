      MODULE OCO2_CO2_MOD
 
      USE TIME_MOD
      USE PRECISION_MOD
      IMPLICIT NONE


      !=================================================================
      ! MODULE VARIABLES
      !=================================================================

      ! Parameters
      INTEGER, PARAMETER           :: MAXLEV = 20
      INTEGER, PARAMETER           :: MAXOCO2 = 1000
      INTEGER, PARAMETER           :: LOCO2 = 20
      ! Record to store data from each GOS obs
      TYPE OCO2_CO2_OBS
         INTEGER                      :: LOCO2(1)
         REAL*8                       :: LAT(1)
         REAL*8                       :: LON(1)
         REAL*8                       :: PSURF(1)
         REAL*8                       :: CO2(MAXLEV)
         REAL*8                       :: PRES(MAXLEV)
         REAL*8                       :: PRIOR(MAXLEV)
         REAL*8                       :: AVG_KERNEL(MAXLEV,MAXLEV)
         REAL*8                       :: S_OER(MAXLEV,MAXLEV)
         REAL*8                       :: S_OER_INV(MAXLEV,MAXLEV)
      ENDTYPE OCO2_CO2_OBS

      TYPE(OCO2_CO2_OBS)                          :: OCO2(MAXOCO2)

      INTEGER, PARAMETER           :: IDTCO2   = 1
      ! Same thing for TCVV(IDTCO2)
      REAL*8,  PARAMETER           :: TCVV_CO2 = 28.97d0  / 44d0

      CONTAINS

      SUBROUTINE READ_OCO2_CO2_OBS( FILENAME, NOCO2 ) 
!
!******************************************************************************
!  Subroutine READ_GOS_CO2_OBS reads the file and passes back info contained
!  therein. (dkh, 10/12/10)
!
!  Based on READ_TES_NH3 OBS (dkh, 04/26/10)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) YYYYMMDD    INTEGER : Current year-month-day
!
!  Arguments as Output:
!  ============================================================================
!  (1 ) NGOS      (INTEGER) : Number of OCO2 retrievals for current hour!

!  Module variable as Output:
!  ============================================================================
!  (1 ) GOS    (GOS_CO2_OB         OCO2(N)%LON              =    lon_varrd_1dr(N)
!S) : CO2 retrieval for current day
!
!  NOTES:
!  (1 ) Add calculation of S_OER_INV, though eventually we probably want to
!        do this offline. (dkh, 05/04/10)
!******************************************************************************
!
      ! Reference to f90 modules

      USE NCDF_MOD
      USE m_netcdf_io_checks,   ONLY : Ncdoes_Var_Exist
      USE m_netcdf_io_read,     ONLY : NcRd
      USE m_netcdf_io_close,      ONLY : Nccl

      ! Record to store data from each GOS obs

      INTEGER, INTENT(INOUT)           :: NOCO2
      CHARACTER(LEN=*), INTENT(IN)   :: FILENAME
      ! local variables
      INTEGER                          :: FID
      INTEGER                          :: N, INC, I

      INTEGER                          :: RC
 
      REAL*8                           :: TEMP(LOCO2)

!variable parameters
      !time length
      INTEGER                          :: time_varrd_1di(1)
      INTEGER                          :: time_str1d(1)
      INTEGER                          :: time_cnt1d(1)
      
      !latitude
      REAL*8, ALLOCATABLE              :: lat_varrd_1dr(:)
      INTEGER                          :: lat_str1d(1)
      INTEGER                          :: lat_cnt1d(1)         

      !longitude
      REAL*8, ALLOCATABLE              :: lon_varrd_1dr(:)
      INTEGER                          :: lon_str1d(1)
      INTEGER                          :: lon_cnt1d(1)  

      !surface pressure
      REAL*8, ALLOCATABLE              :: psurf_varrd_1dr(:)
      INTEGER                          :: psurf_str1d(1)
      INTEGER                          :: psurf_cnt1d(1)

      !prior
      REAL*8, ALLOCATABLE              :: pri_varrd_2dr(:, :)
      INTEGER                          :: pri_str2d(2)
      INTEGER                          :: pri_cnt2d(2)  

      !pressure
      REAL*8, ALLOCATABLE              :: pre_varrd_2dr(:, :)
      INTEGER                          :: pre_str2d(2)
      INTEGER                          :: pre_cnt2d(2)  

      !co2_profile 
      REAL*8, ALLOCATABLE              :: co2_varrd_2dr(:, :)
      INTEGER                          :: co2_str2d(2)
      INTEGER                          :: co2_cnt2d(2)  
        
      !observation error covariance inverse
      REAL*8, ALLOCATABLE              :: Sinv_varrd_2dr(:, :)
      INTEGER                          :: Sinv_str2d(2)
      INTEGER                          :: Sinv_cnt2d(2)  

      !average kernel 
      REAL*8, ALLOCATABLE              :: avk_varrd_3dr(:, :, :)
      INTEGER                          :: avk_str3d(3)
      INTEGER                          :: avk_cnt3d(3)  

      LOGICAL                          :: err_stop
      INTEGER                          :: stat
      CHARACTER(LEN=255)               :: varname 
      ! filename root

      WRITE(*,"(A50,/)") '-READ_OCO2_CO2_OBS: reading file: ', TRIM(FILENAME)

      CALL NC_OPEN( TRIM(FILENAME), fID )

      ! time length
      time_str1d(1) = 1
      time_cnt1d(1) = 1
      varname = 'length'
      CALL Ncrd(time_varrd_1di, FID, TRIM(varname), time_str1d, time_cnt1d, &
     &                        err_stop, stat)
      
      NOCO2 = time_varrd_1di(1)
      ! allocate arrays and read nc file
      !latitude
      ALLOCATE( lat_varrd_1dr(time_varrd_1di(1)), STAT = RC)
      lat_str1d(1) = 1
      lat_cnt1d(1) = time_varrd_1di(1)
      varname = 'latitude'
      CALL Ncrd(lat_varrd_1dr, FID, TRIM(varname), lat_str1d, lat_cnt1d,  &
     &                        err_stop, stat)

      !longitude
      ALLOCATE( lon_varrd_1dr(time_varrd_1di(1)), STAT = RC)
      lon_str1d(1) = 1
      lon_cnt1d(1) = time_varrd_1di(1)
      varname = 'longitude'
      CALL Ncrd(lon_varrd_1dr, FID, TRIM(varname), lon_str1d, lon_cnt1d,  &
     &                        err_stop, stat)

      !surface pressure
      ALLOCATE( psurf_varrd_1dr(time_varrd_1di(1)), STAT = RC)
      psurf_str1d(1) = 1
      psurf_cnt1d(1) = time_varrd_1di(1)
      varname = 'surface_pressure'
      CALL Ncrd(psurf_varrd_1dr, FID, TRIM(varname), psurf_str1d, psurf_cnt1d,  &
     &                        err_stop, stat)

      !pressure
      ALLOCATE( pre_varrd_2dr(LOCO2, time_varrd_1di(1)), STAT = RC)
      pre_str2d(1) = 1
      pre_str2d(2) = 1
      pre_cnt2d(1) = LOCO2
      pre_cnt2d(2) = time_varrd_1di(1)
      varname = 'pressure'
      CALL Ncrd(pre_varrd_2dr, FID, TRIM(varname), pre_str2d, pre_cnt2d)

      !prior
      ALLOCATE( pri_varrd_2dr(LOCO2, time_varrd_1di(1)), STAT = RC)
      pri_str2d(1) = 1
      pri_str2d(2) = 1 
      pri_cnt2d(1) = LOCO2  
      pri_cnt2d(2) = time_varrd_1di(1)
      varname = 'prior'
      CALL Ncrd(pri_varrd_2dr, FID, TRIM(varname), pri_str2d, pri_cnt2d)



      !co2_profile
      ALLOCATE( co2_varrd_2dr(LOCO2, time_varrd_1di(1)), STAT = RC)
      co2_str2d(1) = 1 
      co2_str2d(2) = 1
      co2_cnt2d(1) = LOCO2
      co2_cnt2d(2) = time_varrd_1di(1)
      varname = 'co2_profile'
      CALL NcRd(co2_varrd_2dr, FID, TRIM(varname), co2_str2d, co2_cnt2d)
 

      !average kernel
      ALLOCATE( avk_varrd_3dr(LOCO2, LOCO2, time_varrd_1di(1)), STAT = RC)
      avk_str3d(1) = 1
      avk_str3d(2) = 1
      avk_str3d(3) = 1
      avk_cnt3d(1) = LOCO2 
      avk_cnt3d(2) = LOCO2
      avk_cnt3d(3) = time_varrd_1di(1) 
      varname = 'avg_kernel'  
      CALL NcRd(avk_varrd_3dr, FID, TRIM(varname), avk_str3d, avk_cnt3d)    

      !observation error covariance inverse
      ALLOCATE( Sinv_varrd_2dr(LOCO2, time_varrd_1di(1)), STAT = RC)
 !     Sinv = 0
      ! call subroutine to generate observation error covariance matrix
 !    CALL OBS_ERR_COVARIANCE(CHK_STT(:,:,:,IDTCO2)*TCVV_CO2, &
 ! 			      co2_varrd_2dr(LOCO2:1:-1, N) , NOCO2, Sinv)
      Sinv_str2d(1) = 1
      Sinv_str2d(2) = 1
      Sinv_cnt2d(1) = LOCO2 
      Sinv_cnt2d(2) = time_varrd_1di(1) 

      varname = 'co2_profile_uncert'
      CALL NcRd(Sinv_varrd_2dr, FID, TRIM(varname), Sinv_str2d, Sinv_cnt2d)

      CALL Nccl(FID)
      DO N = 1, NOCO2
         ! observation error covariance
	 OCO2(N)%S_OER_INV = 0
         DO I = 1, LOCO2
            OCO2(N)%S_OER_INV(I, I) = 1/Sinv_varrd_2dr(LOCO2-I+1, N)
         ENDDO
         OCO2(N)%LOCO2            =    LOCO2
         OCO2(N)%LAT              =    lat_varrd_1dr(N)
         OCO2(N)%LON              =    lon_varrd_1dr(N)    
         OCO2(N)%PSURF            =    psurf_varrd_1dr(N)
         OCO2(N)%CO2(1:LOCO2)     =    co2_varrd_2dr(LOCO2:1:-1, N)  
         OCO2(N)%PRES(1:LOCO2)    =    pre_varrd_2dr(LOCO2:1:-1, N)
         ! convert to ppmv
         OCO2(N)%PRIOR(1:LOCO2)   =    pri_varrd_2dr(LOCO2:1:-1, N)
!         print*, 'oco2', N, OCO2(N)%PRIOR(20)
 !        print*, 'oco2', N , time_varrd_1di(1)
         OCO2(N)%AVG_KERNEL(:, :) =    TRANSPOSE(avk_varrd_3dr(:, :, N))
!         OCO2(N)%S_OER_INV(:, :)  =    TRANSPOSE(Sinv_varrd_2dr(:, :, N))
!!       OCO2(N)%S_OER_INV(:, :)  =    Sinv(:, :, N)
!       WRITE(*,*) 'S_OER_INV', OCO2(N)%S_OER_INV
        !rotate matrix 180 deg
        DO INC = 0, 9
        ! averaging kernel
           TEMP(:)                           = OCO2(N)%AVG_KERNEL(LOCO2-INC, :)
           OCO2(N)%AVG_KERNEL(LOCO2-INC, :)  = OCO2(N)%AVG_KERNEL(INC+1, :)
           OCO2(N)%AVG_KERNEL(INC+1, :)      = TEMP(:)
           ! observation error covariance
    !       TEMP(:)                           = OCO2(N)%S_OER_INV(LOCO2-INC, :)
    !       OCO2(N)%S_OER_INV(LOCO2-INC, :)   = OCO2(N)%S_OER_INV(INC+1, :)
    !       OCO2(N)%S_OER_INV(INC+1, :)       = TEMP(:)

        END DO
      END DO
      
      IF ( ALLOCATED( lat_varrd_1dr        ) ) DEALLOCATE( lat_varrd_1dr        )
      IF ( ALLOCATED( lon_varrd_1dr        ) ) DEALLOCATE( lon_varrd_1dr        )
      IF ( ALLOCATED( pre_varrd_2dr        ) ) DEALLOCATE( pre_varrd_2dr        )
      IF ( ALLOCATED( pri_varrd_2dr        ) ) DEALLOCATE( pri_varrd_2dr        )
      IF ( ALLOCATED( co2_varrd_2dr        ) ) DEALLOCATE( co2_varrd_2dr        )
      IF ( ALLOCATED( avk_varrd_3dr        ) ) DEALLOCATE( avk_varrd_3dr        )
      IF ( ALLOCATED( Sinv_varrd_2dr        ) ) DEALLOCATE( Sinv_varrd_2dr       )
      !debug

      END SUBROUTINE READ_OCO2_CO2_OBS

!------------------------------------------------------------------------------

      SUBROUTINE CALC_OCO2_CO2_FORCE( COST_FUNC, State_Chm_Adj, State_Grid_Adj, &
                                      State_Met_Adj )
!
!******************************************************************************
!  Subroutine CALC_GOS_CO2_FORCE calculates the adjoint forcing from the GOSAT
!  CO2 observations and updates the cost function. (dkh, 10/12/10)
!
!
!  Arguments as Input/Output:
!  ============================================================================
!  (1 ) COST_FUNC (REAL*8) : Cost funciton                        [unitless]
!
!
!  NOTES:
!
!******************************************************************************
!
      ! Reference to f90 modules
      USE ADJ_ARRAYS_MOD,     ONLY : N_CALC
      USE ADJ_ARRAYS_MOD,     ONLY : INDEX_ARRAY
 
      USE LOGICAL_ADJ_MOD,    ONLY : FINITE_DIFF

      USE LOGICAL_ADJ_MOD,    ONLY : LTRAN_FINITE_DIFF

      USE LOGICAL_ADJ_MOD,    ONLY : LCONV_FINITE_DIFF

!      USE PRESSURE_MOD,       ONLY : GET_PCENTER, GET_PEDGE
      USE ADJ_ARRAYS_MOD,     ONLY : EXIST_OBS

      USE State_Chm_Adj_Mod, ONLY : ChmStateAdj

      USE State_Grid_Adj_Mod, ONLY : GrdStateAdj

      USE State_Met_Adj_Mod, ONLY : MetStateAdj
      ! Arguments
      REAL*8, INTENT(INOUT)       :: COST_FUNC

      !input parameters
      TYPE(GrdStateAdj), INTENT(IN) :: State_Grid_Adj
      TYPE(ChmStateAdj), INTENT(IN) :: State_Chm_Adj
      TYPE(MetStateAdj), INTENT(IN) :: State_Met_Adj

      ! Local variables
      INTEGER                     :: NTSTART, NTSTOP, NT
      INTEGER                     :: IIJJ(2), I,      J
      INTEGER                     :: L,       LL
      REAL*8                      :: GC_CO2_NATIVE(State_Grid_Adj%NZ)
      REAL*8                      :: GC_PRES(State_Grid_Adj%NZ)
      REAL*8                      :: GC_CO2(MAXLEV)
      REAL*8                      :: GC_PSURF
      REAL*8                      :: MAP(State_Grid_Adj%NZ,MAXLEV)
      REAL*8                      :: CO2_HAT(MAXLEV)
      REAL*8                      :: CO2_PERT(MAXLEV)
      REAL*8                      :: FORCE(MAXLEV)
      REAL*8                      :: DIFF(MAXLEV)
      REAL*8                      :: NEW_COST(MAXOCO2)
      REAL*8                      :: OLD_COST
      INTEGER                     :: NOCO2

      REAL*8                      :: GC_CO2_NATIVE_ADJ(State_Grid_Adj%NZ)
      REAL*8                      :: CO2_HAT_ADJ(MAXLEV)
      REAL*8                      :: CO2_PERT_ADJ(MAXLEV)
      REAL*8                      :: GC_CO2_ADJ(MAXLEV)
      REAL*8                      :: DIFF_ADJ(MAXLEV)

      LOGICAL, SAVE               :: FIRST = .TRUE.
      INTEGER                     :: IOS
   
      LOGICAL                     :: ALIVE
      CHARACTER(LEN=255)          :: EDGE_FILENAME
      CHARACTER(LEN=255)          :: OCO2_FILENAME
      CHARACTER(LEN=255)           :: PRESSURE_FILENAME
      CHARACTER(LEN=20)           :: CHA_NYMD
      CHARACTER(LEN=20)           :: CHA_HOUR
      CHARACTER(LEN=20)           :: CHA_MINUTE

      INTEGER                     :: NHMS
      INTEGER                     :: NYMD
      INTEGER                     :: HOUR 
      INTEGER                     :: MINUTE      
     
      ! center pressure && edge pressure
      REAL*8                      :: P_CENTER(State_Grid_Adj%NX, State_Grid_Adj%NY, State_Grid_Adj%NZ)
      REAL*8                      :: P_EDGE(State_Grid_Adj%NX, State_Grid_Adj%NY, State_Grid_Adj%NZ)   

      !
      REAL*8, ALLOCATABLE              :: Sinv(:, :, :)
      REAL*8, ALLOCATABLE              :: yo(:, :)         
      REAL*8, ALLOCATABLE              :: yM(:, :)  
!                  State_Chm_Adj%SpeciesAdj )
      INTEGER                          :: RC 

      INTEGER                          :: OBS_EXIST(State_Grid_Adj%NX, State_Grid_Adj%NY) 

      INTEGER                          :: Cur_X
      !=================================================================
      ! CALC_GOS_CO2_FORCE begins here!
      !=================================================================

      WRITE(*,'(A30)') 'ENTER CALC_OCO2_CO2_FORCE'
      print '(A30,/)', '     - CALC_OCO2_CO2_FORCE '
      ! Reset
      NEW_COST = 0D0
      OBS_EXIST = 0

      ! Save a value of the cost function first
      OLD_COST = COST_FUNC
      NHMS        = GET_NHMS()
      NYMD        = GET_NYMD()
      HOUR        = NHMS/10000
      MINUTE      = (NHMS-HOUR*10000)/100
      WRITE(CHA_NYMD, *) NYMD
      WRITE(CHA_HOUR, '(I2.2)') HOUR 
      WRITE(CHA_MINUTE, '(I2.2)') MINUTE
      OCO2_FILENAME = "../OCO2_CO2/oco2_"//TRIM(ADJUSTL(CHA_NYMD)) &
			//TRIM(ADJUSTL(CHA_HOUR)) &
                        //TRIM(ADJUSTL(CHA_MINUTE))//"00.nc"
      INQUIRE( FILE = TRIM(OCO2_FILENAME), EXIST = ALIVE )
      IF ( .NOT. ALIVE ) THEN
         print*, ' No matching OCO2 CO2 obs for this hour: ',OCO2_FILENAME
	 WRITE(*,*) 'ALIVE:', ALIVE
         EXIST_OBS = .FALSE.
        RETURN
      ENDIF

     ! read center pressure && edge pressure 
          ! read center pressure 
!      CALL READ_FILE( PRESSURE_FILENAME, P_CENTER, State_Grid_ADJ%NX, State_Grid_ADJ%NY, State_Grid_ADJ%NZ )

     ! read edge pressure 
!      CALL READ_FILE( PRESSURE_FILENAME, P_EDGE, State_Grid_ADJ%NX, State_Grid_ADJ%NY, State_Grid_ADJ%NZ )        

! need to update this in order to do i/o with this loop parallel
!!      ! Now do a parallel loop for analyzing data
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( NT, MAP, LGOS, IIJJ,  I, J,  L,   LL, JLOOP )
!!$OMP+PRIVATE( GC_PRES,       GC_PSURF, GC_CO2,  DIFF  )
!!$OMP+PRIVATE( GC_CO2_NATIVE, CO2_PERT,   CO2_HAT, FORCE )
!!$OMP+PRIVATE( GC_CO2_NATIVE_ADJ,        GC_CO2_ADJ     )
!!$OMP+PRIVATE( CO2_PERT_ADJ,             CO2_HAT_ADJ    )
!!$OMP+PRIVATE( DIFF_ADJ                                )  
      ! read OCO2   
      CALL READ_OCO2_CO2_OBS( TRIM(OCO2_FILENAME), NOCO2 ) 

!      ALLOCATE( yo(LOCO2, NOCO2), STAT = RC)
!      ALLOCATE( yM(LOCO2, NOCO2), STAT = RC)
!      ALLOCATE( Sinv(LOCO2, LOCO2, NOCO2), STAT = RC) 

      ! get yM
!      CALL GENERATE_yM(P_CENTER, P_EDGE, yM, NOCO2, State_Chm_Adj, State_Grid_Adj)

      !generate yo
!      DO NT = 1, NOCO2
!         yo(:, NT) = OCO2(NT)%CO2(:)
!      END DO

      ! calc obs error covariance
!      CALL OBS_ERR_COVARIANCE(yM, yo, NOCO2, Sinv)

      ! output to OCO2
!      DO NT = 1, NOCO2
!         OCO2(NT)%S_OER_INV(:, :) = Sinv(:, :, NT)
!      END DO 

       ! adjoint output
!      OPEN(UNIT=1995, FILE='COST_DELTA', ACCESS='APPEND')
!      OPEN(UNIT=1996, FILE='COST', ACCESS='APPEND')
!      OPEN(UNIT=1997, FILE='ADJOINT', ACCESS='APPEND')
!      OPEN(UNIT=1998, FILE='index_i', ACCESS='APPEND')
!      OPEN(UNIT=1999, FILE='index_j', ACCESS='APPEND')
      DO I = 1, 144
      DO J = 1, 91
      DO L = 1, 47
         IF(State_Chm_Adj%Species(I, J, L, 1)>90.0) THEN

!           print*, 'huoxiao_debug_test', I, J, L
 
         ENDIF
      ENDDO
      ENDDO
      ENDDO

      DO NT  = 1, NOCO2

         ! For safety, initialize these up to LLGOS
         GC_CO2(:)       = 0d0
         MAP(:,:)        = 0d0
         CO2_HAT_ADJ(:)  = 0d0
         FORCE(:)        = 0d0


         ! Copy LGOS to make coding a bit cleaner

         ! Get grid box of current record
         IIJJ  = GET_IJ( REAL(OCO2(NT)%LON(1),4), REAL(OCO2(NT)%LAT(1),4), State_Grid_ADJ)
         I     = IIJJ(1)
         J     = IIJJ(2)
         print*, I, J
         ! obs exist, exit current loop
         IF ( OBS_EXIST(I, J) == 1) CYCLE

         OBS_EXIST(I, J) = 1
         print*, 'huoxiao_debug_ij', I, J 
         ! adjoint verification
!         IF ( LCONV_FINITE_DIFF .OR. LTRAN_FINITE_DIFF ) THEN

!           print*, 'huoxiao debug append'
!           CALL APPEND_IJ( I, J )

!         END IF
         IF (FINITE_DIFF) THEN
            OPEN(UNIT=2000, FILE='cur_x')
            READ(2000, *) Cur_X
            if( I == Index_Array(1, Cur_X) .AND. J == Index_Array(2, Cur_X) ) then
            State_Chm_Adj%Species(I, J, 1, 1) = State_Chm_Adj%Species(I, J, 1, 1)+0.00000001

            endif
            CLOSE(2000)
         ENDIF
!         WRITE(1998, *) I
!         WRITE(1999, *) J
         ! Get GC pressure levels (bar)
         DO L = 1, State_Grid_Adj%NZ
	    ! huoxiao debug
            GC_PRES(L) = State_Met_Adj%PMID(I,J,L)*100
         ENDDO

         ! Get GC surface pressure (bar)
         GC_PSURF = State_Met_Adj%PSC2WET(I,J)*100
         ! Calculate the interpolation weight matrix
         MAP(1:State_Grid_Adj%NZ,1:LOCO2)                                      &
           = GET_INTMAP( State_Grid_Adj%NZ, GC_PRES(:),           GC_PSURF,   &
                         LOCO2,  OCO2(NT)%PRES(1:LOCO2), OCO2(NT)%PSURF(1)  )


         ! Get CO2 values at native model resolution
         GC_CO2_NATIVE(:) = State_Chm_Adj%Species(I,J,:,IDTCO2)
         ! Convert from kg/kg to v/v
         !GC_CO2_NATIVE(:) = GC_CO2_NATIVE(:) * TCVV_CO2
         ! Interpolate GC CO2 column to TES grid
         DO LL = 1, LOCO2
            GC_CO2(LL) = 0d0
            DO L = 1, State_Grid_Adj%NZ
               GC_CO2(LL) = GC_CO2(LL)       &
                          + MAP(L,LL) * GC_CO2_NATIVE(L)
            ENDDO
         ENDDO
         !--------------------------------------------------------------
         ! Apply OCO2 observation operator
         !
         !   x_hat = x_a + A_k ( x_m - x_a )
         !
         !  where
         !    x_hat = GC modeled column as seen by TES [vmr]
         !    x_a   = OCO2 apriori column               [vmr]
         !    x_m   = GC modeled column                [vmr]
         !    A_k   = OCO2 averaging kernel
         !--------------------------------------------------------------
         ! x_m - x_a
         DO L = 1, LOCO2
           CO2_PERT(L) = GC_CO2(L) - OCO2(NT)%PRIOR(L)
         ENDDO
         ! x_a + A_k * ( x_m - x_a )
         DO L = 1, LOCO2
            CO2_HAT(L)    = 0d0
            DO LL = 1, LOCO2
               CO2_HAT(L) = CO2_HAT(L)               &
                          + OCO2(NT)%AVG_KERNEL(L,LL) * CO2_PERT(LL)
            ENDDO
            CO2_HAT(L)    = CO2_HAT(L) + OCO2(NT)%PRIOR(L)

         ENDDO
         !--------------------------------------------------------------
         ! Calculate cost function, given S is error in vmr
         ! J = 1/2 [ model - obs ]^T S_{obs}^{-1} [ model - obs ]
         !--------------------------------------------------------------
         ! Calculate difference between modeled and observed profile
         DO L = 1, LOCO2
            DIFF(L) = CO2_HAT(L) - OCO2(NT)%CO2(L)
         ENDDO
         ! Calculate 1/2 * DIFF^T * S_{obs}^{-1} * DIFF
         DO L = 1, LOCO2
            FORCE(L)     = 0d0
            DO LL = 1, LOCO2
               FORCE(L)  = FORCE(L) + OCO2(NT)%S_OER_INV(L,LL) * DIFF(LL)
            ENDDO
            NEW_COST(NT) = NEW_COST(NT) + 0.5d0 * DIFF(L) * FORCE(L)
         ENDDO
         !--------------------------------------------------------------
         ! Begin adjoint calculations
         !--------------------------------------------------------------
         ! The adjoint forcing  is S_{obs}^{-1} * DIFF = FORCE
         DIFF_ADJ(:) = FORCE(:)
  
         ! Adjoint of difference
         DO L = 1, LOCO2
            CO2_HAT_ADJ(L) =  DIFF_ADJ(L)
         ENDDO
         ! adjoint of TES operator
         DO L  = 1, LOCO2
            CO2_PERT_ADJ(L)    = 0d0
            DO LL = 1, LOCO2
               CO2_PERT_ADJ(L) = CO2_PERT_ADJ(L)           &
                               + OCO2(NT)%AVG_KERNEL(LL,L) &
                               * CO2_HAT_ADJ(LL)
           ENDDO
         ENDDO
         ! huoxiao debug
         ! Adjoint of x_m - x_a
         DO L = 1, LOCO2
           ! adj code:
           GC_CO2_ADJ(L) = CO2_PERT_ADJ(L)
         ENDDO

         ! adjoint of interpolation
         DO L  = 1, State_Grid_Adj%NZ
            GC_CO2_NATIVE_ADJ(L) = 0d0
            DO LL = 1, LOCO2
               GC_CO2_NATIVE_ADJ(L) = GC_CO2_NATIVE_ADJ(L)       &
                                   + MAP(L,LL) * GC_CO2_ADJ(LL)
            ENDDO
         ENDDO
         ! huoxiao debug
         ! Adjoint of unit conversion
         GC_CO2_NATIVE_ADJ(:) = GC_CO2_NATIVE_ADJ(:) * TCVV_CO2
        ! Pass adjoint back to adjoint tracer array
         State_Chm_Adj%SpeciesAdj(I,J,:,IDTCO2)  = State_Chm_Adj%SpeciesAdj(I,J,:,IDTCO2) &        
                                                 + GC_CO2_NATIVE_ADJ(:)

         ! huoxiao debug
        ! finite difference
!         IF ( FINITE_DIFF ) THEN 
!            print*, 'finite diff' 
!            WRITE(1995, *) NEW_COST(NT)
          
!         ELSE

!            WRITE(1996, *) NEW_COST(NT) 
!            WRITE(1997, *) State_Chm_Adj%SpeciesAdj(I,J,1,IDTCO2)
!             WRITE(1995, *) NEW_COST(NT) 
!        ENDIF

      ENDDO  ! NT
!!$OMP END PARALLEL DO

      ! Update cost function
       COST_FUNC = COST_FUNC + SUM(NEW_COST(1:NOCO2))
  
!      print*, ' Updated value of COST_FUNC = ', COST_FUNC
      print*, ' OCO2 contribution           = ', SUM(NEW_COST(1:NOCO2))
!      CLOSE(1995)
!      CLOSE(1996)
!      CLOSE(1997)
!      CLOSE(1998)
!      CLOSE(1999)
      ! Return to calling program
      END SUBROUTINE CALC_OCO2_CO2_FORCE

!------------------------------------------------------------------------------
    
      FUNCTION GET_INTMAP( LGC_TOP, GC_PRESC, GC_SURFP,     &
     &                     LTM_TOP, TM_PRESC, TM_SURFP  )   &
     &         RESULT      ( HINTERPZ )
!
!******************************************************************************
!  Function GET_INTMAP linearly interpolates column quatities
!   based upon the centered (average) pressue levels.
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) LGC_TOP (TYPE) : Description                          [unit]
!  (2 ) GC_PRES (TYPE) : Description                          [unit]
!  (3 ) GC_SURFP(TYPE) : Description                          [unit]
!  (4 ) LTM_TOP (TYPE) : Description                          [unit]
!  (5 ) TM_PRES (TYPE) : Description                          [unit]
!  (6 ) TM_SURFP(TYPE) : Description                          [unit]
!
!  Arguments as Output:
!  ============================================================================
!  (1 ) HINTERPZ (TYPE) : Description                          [unit]
!
!  NOTES:
!  (1 ) Based on the GET_HINTERPZ_2 routine I wrote for read_sciano2_mod.
!
!******************************************************************************
!

      ! Arguments
      INTEGER            :: LGC_TOP, LTM_TOP
      REAL*8             :: GC_PRESC(LGC_TOP)
      REAL*8             :: TM_PRESC(LTM_TOP)
      REAL*8             :: GC_SURFP
      REAL*8             :: TM_SURFP

      ! Return value
      REAL*8             :: HINTERPZ(LGC_TOP, LTM_TOP)

      ! Local variables
      INTEGER  :: LGC, LTM
      REAL*8   :: DIFF, DELTA_SURFP
      REAL*8   :: LOW, HI

      !=================================================================
      ! GET_HINTERPZ_2 begins here!
      !=================================================================

      HINTERPZ(:,:) = 0D0

!      ! Rescale GC grid according to TM surface pressure
!!         p1_A =     (a1 + b1 (ps_A - PTOP))
!!         p2_A =     (a2 + b2 (ps_A - PTOP))
!!         p1_B =     (a + b (ps_B - PTOP))
!!         p2_B =    *(a + b (ps_B - PTOP))
!!         pc_A = 0.5(a1+a2 +(b1+b2)*(ps_A - PTOP))
!!         pc_B = 0.5(a1+a2 +(b1+b2)*(ps_B - PTOP))
!!         pc_B - pc_A = 0.5(b1_b2)(ps_B-ps_A)
!!         pc_B = 0.5(b1_b2)(ps_B-ps_A) + pc_A
!      DELTA_SURFP   = 0.5d0 * ( TM_SURFP -GC_SURFP )
!
!      DO LGC = 1, LGC_TOP
!         GC_PRESC(LGC) = ( GET_BP(LGC) + GET_BP(LGC+1))
!     &               * DELTA_SURFP + GC_PRESC(LGC)
!         IF (GC_PRESC(LGC) < 0) THEN
!            CALL ERROR_STOP( 'highly unlikey',
!     &                       'read_sciano2_mod.f')
!         ENDIF
!
!      ENDDO


      ! Loop over each pressure level of TM grid
      DO LTM = 1, LTM_TOP

         ! Find the levels from GC that bracket level LTM
         DO LGC = 1, LGC_TOP - 1

            LOW = GC_PRESC(LGC+1)
            HI  = GC_PRESC(LGC)
            IF (LGC == 0) HI = TM_SURFP

            ! Linearly interpolate value on the LTM grid
            IF ( TM_PRESC(LTM) <= HI .and.             &
     &           TM_PRESC(LTM)  > LOW) THEN

               DIFF                = HI - LOW
               HINTERPZ(LGC+1,LTM) = ( HI - TM_PRESC(LTM)  ) / DIFF
               HINTERPZ(LGC  ,LTM) = ( TM_PRESC(LTM) - LOW ) / DIFF


            ENDIF

            ! dkh debug
            !print*, 'LGC,LTM,HINT', LGC, LTM, HINTERPZ(LGC,LTM)

          ENDDO
       ENDDO
!       WRITE(*,*) 'HINTERPZ',HINTERPZ

       ! Bug fix:  a more general version allows for multiples TES pressure
       ! levels to exist below the lowest GC pressure.  (dm, dkh, 09/30/10)
       ! OLD code:
       !IF ( TM_PRESC(1) > GC_PRESC(1) ) THEN
       !   HINTERPZ(1,1)         = 1D0
       !   HINTERPZ(2:LGC_TOP,1) = 0D0
       !ENDIF
       ! New code:
       ! Loop over each pressure level of TM grid
       DO LTM = 1, LTM_TOP
          IF ( TM_PRESC(LTM) > GC_PRESC(1) ) THEN
             HINTERPZ(1,LTM)         = 1D0
             HINTERPZ(2:LGC_TOP,LTM) = 0D0
          ENDIF
       ENDDO

 
      ! Return to calling program
      END FUNCTION GET_INTMAP
!---------------------------------------------------------------------------
       !
      SUBROUTINE OBS_ERR_COVARIANCE(yM, yo, length, Sinv)

      ! calculate observation error covariance matrix 
      ! reference paper sec-3.1: 
      !  Comparative inverse analysis of satellite (MOPITT) and aircraft
      !  (TRACE-P) observations to estimate Asian sources of carbon monoxide
      ! yM: geos-chem output concentration 
      ! yo: oco2 observation concentration 
      ! length: observation numbers

      ! input parameters
      INTEGER,INTENT(IN)         :: length
      REAL*8,INTENT(IN)          :: yM(LOCO2, length)
      REAL*8,INTENT(IN)          :: yo(LOCO2, length)

      ! output parameters
      REAL*8,INTENT(OUT)         :: Sinv(LOCO2, LOCO2, length)
      ! local variables
      
      INTEGER                    :: LEVEL, NUM
      REAL*8                     :: b(LOCO2)
 
      REAL*8                     :: Resid(LOCO2, length)
      ! standard deriviation of the residual error       
      REAL*8       	   	 :: EPSON(LOCO2, length)
      
      Resid = yM(:, :)-yo(:, :)
      ! cal mean(yM-y)
      DO LEVEL  = 1, LOCO2
         b(LEVEL) = SUM(Resid(LEVEL, :))/length
      END DO
      ! epson = yM-yo-b
      EPSON(:, :) = yM(:, :)-yo(:, :)
      ! observation error covariance
      DO NUM = 1, length
         EPSON(:, NUM) = EPSON(:, NUM)-b
      END DO
       
      ! generate observation error covariance matrix
      DO LEVEL = 1 ,LOCO2
         DO NUM = 1, length
            Sinv(LEVEL, LEVEL, NUM) = EPSON(LEVEL, NUM)
         END DO
      END DO 
      
      END SUBROUTINE OBS_ERR_COVARIANCE


      SUBROUTINE GENERATE_yM(P_CENTER, P_EDGE, yM, NOCO2, State_Chm_Adj, State_Grid_Adj)

 !
!******************************************************************************
!  Subroutine CALC_GOS_CO2_FORCE calculates the adjoint forcing from the GOSAT
!  CO2 observations and updates the cost function. (dkh, 10/12/10)
!
!
!  Arguments as Input/Output:
!  ============================================================================
!  (1 ) COST_FUNC (REAL*8) : Cost funciton                        [unitless]
!
!
!  NOTES:
!
!******************************************************************************
!
      ! Reference to f90 modules

!      USE GC_GRID_MOD,        ONLY : GET_IJ
!      USE PRESSURE_MOD,       ONLY : GET_PCENTER, GET_PEDGE
      USE ADJ_ARRAYS_MOD,     ONLY : EXIST_OBS
      USE State_Grid_Adj_Mod, ONLY : GrdStateAdj
      USE State_Chm_Adj_Mod,  ONLY : ChmStateAdj
      !input parameters
      TYPE(GrdStateAdj), INTENT(IN) :: State_Grid_Adj          
      TYPE(ChmStateAdj), INTENT(IN) :: State_Chm_Adj
      ! Local variables
      INTEGER                       :: NTSTART, NTSTOP, NT
      INTEGER                       :: IIJJ(2), I,      J
      INTEGER                       :: L,       LL
      REAL*8                        :: GC_PRES(State_Grid_Adj%NZ)
      REAL*8                        :: GC_CO2_NATIVE(State_Grid_Adj%NZ)
      REAL*8                        :: GC_CO2(MAXLEV)
      REAL*8                        :: GC_PSURF
      REAL*8                        :: MAP(State_Grid_Adj%NZ,MAXLEV)
      REAL*8                        :: CO2_HAT(MAXLEV)
      REAL*8                        :: CO2_PERT(MAXLEV)


      LOGICAL, SAVE                 :: FIRST = .TRUE.
      INTEGER                       :: IOS
   
      INTEGER                       :: VAR
      CHARACTER(LEN=255)            :: FILENAME
      CHARACTER(LEN=20)             :: CHA_NYMD
      CHARACTER(LEN=20)             :: CHA_HOUR
      CHARACTER(LEN=20)             :: CHA_MINUTE

      INTEGER                       :: NHMS
      INTEGER                       :: NYMD
      INTEGER                       :: HOUR 
      INTEGER                       :: MINUTE      
      ! center pressure && edge pressure
      REAL*8, INTENT(IN)            :: P_CENTER(State_Grid_Adj%NX, State_Grid_Adj%NY, State_Grid_Adj%NZ)
      REAL*8, INTENT(IN)            :: P_EDGE(State_Grid_Adj%NX, State_Grid_Adj%NY, State_Grid_Adj%NZ)
      INTEGER, INTENT(IN)           :: NOCO2
      REAL*8, INTENT(OUT)           :: yM(LOCO2, NOCO2)



! need to update this in order to do i/o with this loop parallel
!!      ! Now do a parallel loop for analyzing data
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( NT, MAP, LGOS, IIJJ,  I, J,  L,   LL, JLOOP )
!!$OMP+PRIVATE( GC_PRES,       GC_PSURF, GC_CO2,  DIFF  )
!!$OMP+PRIVATE( GC_CO2_NATIVE, CO2_PERT,   CO2_HAT, FORCE )
!!$OMP+PRIVATE( GC_CO2_NATIVE_ADJ,        GC_CO2_ADJ     )
!!$OMP+PRIVATE( CO2_PERT_ADJ,             CO2_HAT_ADJ    )
!!$OMP+PRIVATE( DIFF_ADJ                                )
      DO NT  = 1, NOCO2
!	  print*, 'oco2', OCO2(NT)%PRIOR(:)
 !        print*, '     - CALC_OCO2_CO2_FORCE: analyzing record '

         ! For safety, initialize these up to LLGOS
         GC_CO2(:)       = 0d0
         MAP(:,:)        = 0d0


         ! Copy LGOS to make coding a bit cleaner
   !      LOCO2 = OCO2(NT)%LOCO2(1)

         ! Get grid box of current record
         IIJJ  = GET_IJ( REAL(OCO2(NT)%LON(1),4), REAL(OCO2(NT)%LAT(1),4), State_Grid_ADJ)
         I     = IIJJ(1)
         J     = IIJJ(2)


         ! Get GC pressure levels (bar)
         DO L = 1, State_Grid_Adj%NZ
	    ! huoxiao debug
            GC_PRES(L) = P_CENTER(I,J,L)*100
         ENDDO
	
         ! Get GC surface pressure (bar)
         GC_PSURF = P_EDGE(I,J,1)*100
         
         ! Calculate the interpolation weight matrix
         MAP(1:State_Grid_Adj%NZ,1:LOCO2)                                      &
           = GET_INTMAP( State_Grid_Adj%NZ, GC_PRES(:),           GC_PSURF,   &
                         LOCO2,  OCO2(NT)%PRES(1:LOCO2), GC_PSURF  )


         ! Get CO2 values at native model resolution
         GC_CO2_NATIVE(:) = State_Chm_Adj%Species(I,J,:,IDTCO2)

         ! Convert from kg/kg to v/v
         GC_CO2_NATIVE(:) = GC_CO2_NATIVE(:) * TCVV_CO2

         ! Interpolate GC CO2 column to TES grid
         DO LL = 1, LOCO2
            GC_CO2(LL) = 0d0
            DO L = 1, State_Grid_Adj%NZ
               GC_CO2(LL) = GC_CO2(LL)       &
     &                    + MAP(L,LL) * GC_CO2_NATIVE(L)
            ENDDO
         ENDDO
         !--------------------------------------------------------------
         ! Apply GOS observation operator
         !
         !   x_hat = x_a + A_k ( x_m - x_a )
         !
         !  where
         !    x_hat = GC modeled column as seen by TES [vmr]
         !    x_a   = GOS apriori column               [vmr]
         !    x_m   = GC modeled column                [vmr]
         !    A_k   = GOS averaging kernel
         !--------------------------------------------------------------

         ! x_m - x_a
         DO L = 1, LOCO2
           CO2_PERT(L) = GC_CO2(L) - OCO2(NT)%PRIOR(L)
         ENDDO

	! WRITE(*,*) 'OCO2(NT)%PRIOR', OCO2(NT)%PRIOR(1:LOCO2)
 !        WRITE(*,*) ' OCO2(NT)%PRIOR',  OCO2(NT)%PRIOR
         ! x_a + A_k * ( x_m - x_a )
         DO L = 1, LOCO2
            CO2_HAT(L)    = 0d0
            DO LL = 1, LOCO2
               CO2_HAT(L) = CO2_HAT(L)               &
     &                    + OCO2(NT)%AVG_KERNEL(L,LL) * CO2_PERT(LL)
            ENDDO
            CO2_HAT(L)    = CO2_HAT(L) + OCO2(NT)%PRIOR(L)

         ENDDO
         yM(:, NT) = CO2_HAT
      END DO 
      
      END SUBROUTINE GENERATE_yM

      !

      FUNCTION GET_IJ( LON, LAT, State_Grid_Adj ) RESULT ( IIJJ )

      USE State_Grid_Adj_Mod, ONLY : GrdStateAdj
!
! !INPUT PARAMETERS:
!
      REAL*4,         INTENT(IN)  :: LAT, LON
      TYPE(GrdStateAdj), INTENT(IN)  :: State_Grid_Adj  ! Grid State object
!
! !RETURN VALUE:
!
      INTEGER                     :: IIJJ(2)

! !LOCAL VARIABLES:
!
      REAL(fp) :: TLON, TLAT

      TLON = INT( ( LON + 180e+0_fp ) / State_Grid_Adj%DX + 1.5e+0_fp )
      TLAT = INT( ( LAT +  90e+0_fp ) / State_Grid_ADJ%DY + 1.5e+0_fp )

      !this is not nested grid
      IF ( TLON > State_Grid_Adj%NX ) TLON = TLON - State_Grid_Adj%NX
      IIJJ(1) = TLON
      IIJJ(2) = TLAT

      END FUNCTION GET_IJ
 
      SUBROUTINE APPEND_IJ(I, J)

      INTEGER,        INTENT(IN) :: I, J

      OPEN( UNIT = 2019, FILE = 'location_i', ACCESS='APPEND' )
      OPEN( UNIT = 2020, FILE = 'location_j', ACCESS='APPEND' )

      WRITE(2019, *) I
      WRITE(2020, *) J

      CLOSE(2019)
      CLOSE(2020)
      END SUBROUTINE



    
      END MODULE OCO2_CO2_MOD


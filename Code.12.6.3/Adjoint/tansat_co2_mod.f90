      MODULE TANSAT_CO2_MOD
 
      USE TIME_MOD
      USE PRECISION_MOD
      IMPLICIT NONE


      !=================================================================
      ! MODULE VARIABLES
      !=================================================================

      ! Parameters
      INTEGER, PARAMETER           :: MAXLEV = 20
      INTEGER, PARAMETER           :: MAXTANSAT = 1000
      INTEGER, PARAMETER           :: LTANSAT = 20
      ! Record to store data from each GOS obs
      TYPE TANSAT_CO2_OBS
         INTEGER                      :: LTANSAT(1)
         REAL*8                       :: LAT(1)
         REAL*8                       :: LON(1)

         REAL*8                       :: PSURF(1)
         REAL*8                       :: PRES(MAXLEV)
         REAL*8                       :: PWF(MAXLEV)

         REAL*8                       :: XCO2
         REAL*8                       :: XCO2A
         REAL*8                       :: CO2A(MAXLEV)
         REAL*8                       :: AVG_KERNEL(MAXLEV)
         REAL*8                       :: S_OER(MAXLEV)
         REAL*8                       :: S_OER_INV
      ENDTYPE TANSAT_CO2_OBS

      TYPE(TANSAT_CO2_OBS)                          :: TANSAT(MAXTANSAT)

      INTEGER, PARAMETER           :: IDTCO2   = 1
      ! Same thing for TCVV(IDTCO2)
      REAL*8,  PARAMETER           :: TCVV_CO2 = 44d0/28.97d0
      REAL*8,  PARAMETER           :: g = 9.8
      CONTAINS

      FUNCTION PRESSURE_WEIGHTING_FUNCTION(model_pressure, surface_pressure,&
                                          q,  level)
!     level = 1 is close to surface

      USE PARAMETER_ADJ_MOD,        ONLY : Mdry

      INTEGER, INTENT(IN)               :: level
      REAL(fp), INTENT(IN)              :: surface_pressure
!      REAL(fp), INTENT(IN)              :: hyai(level+1)
!      REAL(fp), INTENT(IN)              :: hybi(level+1)
      REAL(fp), INTENT(IN)              :: model_pressure(level)
      REAL(fp), INTENT(IN)              :: q(level)
      !local
!      REAL(fp)                          :: out_pressure(level+1)
      REAL(fp)                          :: delta_pressure(level-1)
      REAL(fp)                          :: hdot(level-1)
      REAL(fp)                          :: air_sum
      REAL(fp)                          :: c(level-1)
      REAL(fp)                          :: fi, fs
      REAL(fp)                          :: pressure_weighting_function(level)
      INTEGER                           :: i
!
!******************************************************************************
!

!     Reference paper: 
!     The ACOS CO2 retrieval algorithm -Part1:Description and validation against synthetic observations
!      do i = 1, level+1
!         out_pressure(i) = hyai(i)+hybi(i)*surface_pressure
!      enddo
     
!assume  
      do i = 1, level-1       
         if (i == 1) then
            delta_pressure(i) = surface_pressure-model_pressure(i+1)
         else
            delta_pressure(i) = model_pressure(i)-model_pressure(i+1)
         endif
      enddo
      air_sum = 0.0
      DO i = 1, level-1
         c(i) = (1-q(i+1))/(g*Mdry)
         air_sum = air_sum+c(i)*delta_pressure(i)
      enddo

      do i = 1, level-1
         hdot(i) = c(i)*delta_pressure(i)/air_sum
      enddo
      fs = (surface_pressure-model_pressure(2))/(model_pressure(1)-model_pressure(2))
      fi = 1/2

      do i = 1, level
         if (i == 1) then
            pressure_weighting_function(i) = fs*fi*hdot(i)
         elseif (i == 2) then
            pressure_weighting_function(i) = fi*hdot(i)+(1-fs*fi)*hdot(i-1)
         elseif (i == level) then
            pressure_weighting_function(i) = (1-fi)*hdot(i-1)
         else 
            pressure_weighting_function(i) = fi*hdot(i)+(1-fi)*hdot(i-1)
         endif
      enddo
      END FUNCTION PRESSURE_WEIGHTING_FUNCTION
      SUBROUTINE READ_OCO2_CO2_OBS( FILENAME, NTANSAT ) 
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

      INTEGER, INTENT(INOUT)           :: NTANSAT
      CHARACTER(LEN=*), INTENT(IN)   :: FILENAME
      ! local variables
      INTEGER                          :: FID
      INTEGER                          :: N, INC, I

      INTEGER                          :: RC
 
      REAL*8                           :: TEMP(LTANSAT)

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


      !pressure
      REAL*8, ALLOCATABLE              :: pre_varrd_2dr(:, :)
      INTEGER                          :: pre_str2d(2)
      INTEGER                          :: pre_cnt2d(2)

       !specific humudity

      REAL*8, ALLOCATABLE              :: pwf_varrd_2dr(:, :)
      INTEGER                          :: pwf_str2d(2)
      INTEGER                          :: pwf_cnt2d(2)

!      REAL*8, ALLOCATABLE              :: psurf_varrd_1dr(:)
!      INTEGER                          :: psurf_str1d(1)
!      INTEGER                          :: psurf_cnt1d(1)

      !prior
      REAL*8, ALLOCATABLE              :: pri_varrd_2dr(:, :)
      INTEGER                          :: pri_str2d(2)
      INTEGER                          :: pri_cnt2d(2)  


      !pressure
!      REAL*8, ALLOCATABLE              :: pre_varrd_2dr(:, :)
!      INTEGER                          :: pre_str2d(2)
!      INTEGER                          :: pre_cnt2d(2)  

      !co2_profile 
!      REAL*8, ALLOCATABLE              :: co2_varrd_2dr(:, :)
!      INTEGER                          :: co2_str2d(2)
!      INTEGER                          :: co2_cnt2d(2)  

      !hyai
!      REAL*8, ALLOCATABLE              :: hyai_varrd_2dr(:,:)
!      INTEGER                          :: hyai_str1d(2)
!      INTEGER                          :: hyai_cnt1d(2)

      !hybi
!      REAL*8, ALLOCATABLE              :: hybi_varrd_2dr(:,:)
!      INTEGER                          :: hybi_str1d(2)
!      INTEGER                          :: hybi_cnt1d(2)


      !xco2
      REAL*8, ALLOCATABLE              :: xco2_varrd_1dr(:)
      INTEGER                          :: xco2_str1d(1)
      INTEGER                          :: xco2_cnt1d(1)

      REAL*8, ALLOCATABLE              :: xco2_apri_varrd_1dr(:)        
      INTEGER                          :: xco2_apri_str1d(1)
      INTEGER                          :: xco2_apri_cnt1d(1)

      !observation error covariance inverse
      REAL*8, ALLOCATABLE              :: Sinv_varrd_1dr(:)
      INTEGER                          :: Sinv_str1d(1)
      INTEGER                          :: Sinv_cnt1d(1)  

      !average kernel 
      REAL*8, ALLOCATABLE              :: avk_varrd_2dr(:, :)
      INTEGER                          :: avk_str2d(2)
      INTEGER                          :: avk_cnt2d(2)  

      LOGICAL                          :: err_stop
      INTEGER                          :: stat
      CHARACTER(LEN=255)               :: varname 
      ! filename root


!      WRITE(*,"(A50,/)") '-READ_OCO2_CO2_OBS: reading file: ',  &
!                         TRIM(FILENAME)
!      print*, 'open the file'
      CALL NC_OPEN( TRIM(FILENAME), fID )
!      print*, 'read length'
      ! time length
      time_str1d(1) = 1
      time_cnt1d(1) = 1
      varname = 'length'
      CALL Ncrd(time_varrd_1di, FID, TRIM(varname), time_str1d, time_cnt1d, &
     &                        err_stop, stat)
      
      NTANSAT = time_varrd_1di(1)
      ! allocate arrays and read nc file
!      print*, 'read latitude'
      !latitude
      ALLOCATE( lat_varrd_1dr(time_varrd_1di(1)), STAT = RC)
      lat_str1d(1) = 1
      lat_cnt1d(1) = time_varrd_1di(1)
      varname = 'latitude'
      CALL Ncrd(lat_varrd_1dr, FID, TRIM(varname), lat_str1d, lat_cnt1d,  &
     &                        err_stop, stat)

!       print*, 'read longitude'
      !longitude
      ALLOCATE( lon_varrd_1dr(time_varrd_1di(1)), STAT = RC)
      lon_str1d(1) = 1
      lon_cnt1d(1) = time_varrd_1di(1)
      varname = 'longitude'
      CALL Ncrd(lon_varrd_1dr, FID, TRIM(varname), lon_str1d, lon_cnt1d,  &
     &                        err_stop, stat)

      !surface pressure
!      ALLOCATE( psurf_varrd_1dr(time_varrd_1di(1)), STAT = RC)
!      psurf_str1d(1) = 1
!      psurf_cnt1d(1) = time_varrd_1di(1)
!      varname = 'surface_pressure'
!      CALL Ncrd(psurf_varrd_1dr, FID, TRIM(varname), psurf_str1d, psurf_cnt1d,  &
!     &                        err_stop, stat)

      !pressure
!      ALLOCATE( pre_varrd_2dr(LOCO2, time_varrd_1di(1)), STAT = RC)
!      pre_str2d(1) = 1
!      pre_str2d(2) = 1
!      pre_cnt2d(1) = LOCO2
!      pre_cnt2d(2) = time_varrd_1di(1)
!      varname = 'pressure'
!      CALL Ncrd(pre_varrd_2dr, FID, TRIM(varname), pre_str2d, pre_cnt2d)

      !prior
!       print*, 'read prior'
      ALLOCATE( pri_varrd_2dr(LTANSAT, time_varrd_1di(1)), STAT = RC)
      pri_str2d(1) = 1
      pri_str2d(2) = 1 
      pri_cnt2d(1) = LTANSAT
      pri_cnt2d(2) = time_varrd_1di(1)
      varname = 'co2_profile_apriori'
      CALL Ncrd(pri_varrd_2dr, FID, TRIM(varname), pri_str2d, pri_cnt2d)


      !surface pressure
      ALLOCATE( psurf_varrd_1dr(time_varrd_1di(1)), STAT = RC)
      psurf_str1d(1) = 1
      psurf_cnt1d(1) = time_varrd_1di(1)
      varname = 'surface_pressure'
      CALL Ncrd(psurf_varrd_1dr, FID, TRIM(varname), psurf_str1d, psurf_cnt1d,  &
     &                        err_stop, stat)

      !pressure
      ALLOCATE( pre_varrd_2dr(LTANSAT, time_varrd_1di(1)), STAT = RC)
      pre_str2d(1) = 1
      pre_str2d(2) = 1
      pre_cnt2d(1) = LTANSAT
      pre_cnt2d(2) = time_varrd_1di(1)
      varname = 'pressure'
      CALL Ncrd(pre_varrd_2dr, FID, TRIM(varname), pre_str2d, pre_cnt2d)

      !specific humudity
      ALLOCATE( pwf_varrd_2dr(LTANSAT, time_varrd_1di(1)), STAT = RC)
      pwf_str2d(1) = 1
      pwf_str2d(2) = 1
      pwf_cnt2d(1) = LTANSAT
      pwf_cnt2d(2) = time_varrd_1di(1)
      varname = 'xco2_pressure_weighting_function'
      CALL Ncrd(pwf_varrd_2dr, FID, TRIM(varname), pwf_str2d, pwf_cnt2d)


      !xco2
      ALLOCATE( xco2_varrd_1dr(time_varrd_1di(1)), STAT = RC)
      xco2_str1d(1) = 1 
      xco2_cnt1d(1) = time_varrd_1di(1)
      varname = 'xco2'
      CALL NcRd(xco2_varrd_1dr, FID, TRIM(varname), xco2_str1d, xco2_cnt1d)
 

      !xco2 aprior
      ALLOCATE( xco2_apri_varrd_1dr(time_varrd_1di(1)), STAT = RC)
      xco2_apri_str1d(1) = 1
      xco2_apri_cnt1d(1) = time_varrd_1di(1)
      varname = 'xco2_apriori'
      CALL NcRd(xco2_apri_varrd_1dr, FID, TRIM(varname), xco2_apri_str1d, xco2_apri_cnt1d)

      ALLOCATE( avk_varrd_2dr(LTANSAT, time_varrd_1di(1)), STAT = RC)
      avk_str2d(1) = 1
      avk_str2d(2) = 1
      avk_cnt2d(1) = LTANSAT 
      avk_cnt2d(2) = time_varrd_1di(1) 
      varname = 'xco2_avg_kernel'  
      CALL NcRd(avk_varrd_2dr, FID, TRIM(varname), avk_str2d, avk_cnt2d)    

      !observation error covariance inverse
      ALLOCATE( Sinv_varrd_1dr(time_varrd_1di(1)), STAT = RC)
      Sinv_str1d(1) = 1
      Sinv_cnt1d(1) = time_varrd_1di(1) 

      varname = 'xco2_uncert'
      CALL NcRd(Sinv_varrd_1dr, FID, TRIM(varname), Sinv_str1d, Sinv_cnt1d)

      CALL Nccl(FID)
!      print*, 'close file'
      DO N = 1, NTANSAT
         ! observation error covariance
       
         TANSAT(N)%S_OER_INV        = 1/(Sinv_varrd_1dr(N)**2)
         TANSAT(N)%LTANSAT          =    LTANSAT
         TANSAT(N)%LAT              =    lat_varrd_1dr(N)
         TANSAT(N)%LON              =    lon_varrd_1dr(N)   

         TANSAT(N)%PSURF            =    psurf_varrd_1dr(N)
         TANSAT(N)%PRES(1:LTANSAT)    =    pre_varrd_2dr(:, N)
 
         TANSAT(N)%PWF(1:LTANSAT)    =    &
                                         pwf_varrd_2dr(:, N)

         TANSAT(N)%XCO2             =    xco2_varrd_1dr(N)  
         TANSAT(N)%XCO2A            =    xco2_apri_varrd_1dr(N)
            ! convert to ppmv
         TANSAT(N)%CO2A(1:LTANSAT)   =    pri_varrd_2dr(:, N)
!         print*, 'oco2', N, OCO2(N)%PRIOR(20)
         TANSAT(N)%AVG_KERNEL(1:LTANSAT) =    avk_varrd_2dr(:, N)
!         OCO2(N)%S_OER_INV(:, :)  =    TRANSPOSE(Sinv_varrd_2dr(:, :, N))
!!       OCO2(N)%S_OER_INV(:, :)  =    Sinv(:, :, N)
!       WRITE(*,*) 'S_OER_INV', OCO2(N)%S_OER_INV
        !rotate matrix 180 deg
!        DO INC = 0, 9
        ! averaging kernel
!           TEMP(:)                           = OCO2(N)%AVG_KERNEL(LOCO2-INC, :)
!           OCO2(N)%AVG_KERNEL(LOCO2-INC, :)  = OCO2(N)%AVG_KERNEL(INC+1, :)
!           OCO2(N)%AVG_KERNEL(INC+1, :)      = TEMP(:)
           ! observation error covariance
!           TEMP(:)                           = OCO2(N)%S_OER_INV(LOCO2-INC, :)
!           OCO2(N)%S_OER_INV(LOCO2-INC, :)   = OCO2(N)%S_OER_INV(INC+1, :)
!           OCO2(N)%S_OER_INV(INC+1, :)       = TEMP(:)

!        END DO
      END DO
      
      IF ( ALLOCATED( lat_varrd_1dr        ) ) DEALLOCATE( lat_varrd_1dr        )
      IF ( ALLOCATED( lon_varrd_1dr        ) ) DEALLOCATE( lon_varrd_1dr        )
      IF ( ALLOCATED( pwf_varrd_2dr        ))  DEALLOCATE( pwf_varrd_2dr )
      IF ( ALLOCATED( psurf_varrd_1dr        ) ) DEALLOCATE( psurf_varrd_1dr        )
      IF ( ALLOCATED( pre_varrd_2dr        ) ) DEALLOCATE( pre_varrd_2dr        )
      IF ( ALLOCATED( pri_varrd_2dr        ) ) DEALLOCATE( pri_varrd_2dr        )
      IF ( ALLOCATED( xco2_varrd_1dr        ) ) DEALLOCATE( xco2_varrd_1dr        )
      IF ( ALLOCATED( avk_varrd_2dr        ) ) DEALLOCATE( avk_varrd_2dr        )
      IF ( ALLOCATED( Sinv_varrd_1dr        ) ) DEALLOCATE( Sinv_varrd_1dr       )
      !debug

      END SUBROUTINE READ_OCO2_CO2_OBS

!------------------------------------------------------------------------------

      SUBROUTINE CALC_TANSAT_CO2_FORCE( COST_FUNC, State_Chm_Adj, State_Grid_Adj, &
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

      USE LOGICAL_ADJ_MOD,    ONLY : FINITE_DIFF

      USE LOGICAL_ADJ_MOD,    ONLY : LTRAN_FINITE_DIFF

      USE LOGICAL_ADJ_MOD,    ONLY : LCONV_FINITE_DIFF

!      USE PRESSURE_MOD,       ONLY : GET_PCENTER, GET_PEDGE
      USE ADJ_ARRAYS_MOD,     ONLY : EXIST_OBS
      USE ADJ_ARRAYS_MOD,     ONLY : OBS_EXIST 
      USE State_Chm_Adj_Mod, ONLY : ChmStateAdj

      USE State_Grid_Adj_Mod, ONLY : GrdStateAdj

      USE State_Met_Adj_Mod, ONLY : MetStateAdj
      ! Arguments
      REAL*8, INTENT(INOUT)       :: COST_FUNC

      !input parameters
      TYPE(GrdStateAdj), INTENT(IN) :: State_Grid_Adj
      TYPE(ChmStateAdj), INTENT(INOUT) :: State_Chm_Adj
      TYPE(MetStateAdj), INTENT(IN) :: State_Met_Adj

      ! Local variables
      INTEGER                     :: NTSTART, NTSTOP, NT
      INTEGER                     :: IIJJ(2), I,      J
      INTEGER                     :: L,       LL
      REAL*8                      :: GC_CO2_NATIVE(State_Grid_Adj%NZ)
      REAL*8                      :: GC_CO2(MAXLEV)
      REAL*8                      :: CO2_PERT(State_Grid_Adj%NZ)
      REAL*8                      :: FORCE
      REAL*8                      :: DIFF
      REAL*8                      :: NEW_COST(MAXTANSAT)
      REAL*8                      :: OLD_COST
      INTEGER                     :: NTANSAT

      REAL*8                      :: GC_CO2_NATIVE_ADJ(State_Grid_Adj%NZ)
      REAL*8                      :: GC_CO2_ADJ(MAXLEV)
      REAL*8                      :: DIFF_ADJ
      REAL*8                      :: GC_XCO2
      REAL*8                      :: XCO2_PERT

      REAL*8                      :: XCO2_PERT_ADJ
      REAL*8                      :: h(LTANSAT)


      REAL*8                      :: MAP(State_Grid_Adj%NZ,MAXLEV)

      LOGICAL, SAVE               :: FIRST = .TRUE.
      INTEGER                     :: IOS
   
      LOGICAL                     :: ALIVE
      CHARACTER(LEN=255)          :: EDGE_FILENAME
      CHARACTER(LEN=255)          :: TANSAT_FILENAME
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
      REAL*8                      :: GC_PRES(State_Grid_Adj%NZ)
      REAL*8                      :: GC_PSURF
      !
      REAL*8, ALLOCATABLE              :: Sinv(:, :, :)
      REAL*8, ALLOCATABLE              :: yo(:, :)         
      REAL*8, ALLOCATABLE              :: yM(:, :)  
      INTEGER                          :: RC 

      !=================================================================
      ! CALC_GOS_CO2_FORCE begins here!
      !=================================================================

      WRITE(*,'(A30)') 'ENTER CALC_OCO2_CO2_FORCE'
      print '(A30,/)', '     - CALC_OCO2_CO2_FORCE '

      ! Reset
      NEW_COST = 0D0

      ! Save a value of the cost function first
      OLD_COST = COST_FUNC
      NHMS        = GET_NHMS()
      NYMD        = GET_NYMD()
      HOUR        = NHMS/10000
      MINUTE      = (NHMS-HOUR*10000)/100
      WRITE(CHA_NYMD, *) NYMD
      WRITE(CHA_HOUR, '(I2.2)') HOUR 
      WRITE(CHA_MINUTE, '(I2.2)') MINUTE
      TANSAT_FILENAME =                                          & 
      "../OCO2_CO2_Litev9_ND_Observation_Error_Covariance/oco2_" &
!      "../temp/oco2_"   &
      //TRIM(ADJUSTL(CHA_NYMD)) &
      //TRIM(ADJUSTL(CHA_HOUR)) &
      //TRIM(ADJUSTL(CHA_MINUTE))//"00.nc"
      INQUIRE( FILE = TRIM(TANSAT_FILENAME), EXIST = ALIVE )
      WRITE(*,*) 'HUOXIAO_FILENAME', TRIM(TANSAT_FILENAME)
      IF ( .NOT. ALIVE ) THEN
         print*, ' No matching OCO2 CO2 obs for this hour'
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
      CALL READ_OCO2_CO2_OBS( TRIM(TANSAT_FILENAME), NTANSAT ) 
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
!      OPEN(UNIT=1995, FILE='OBS_COST_PLUS', ACCESS='APPEND')
!      OPEN(UNIT=1996, FILE='OBS_COST_MINUS', ACCESS='APPEND')
!      OPEN(UNIT=1997, FILE='OBS_ADJOINT', ACCESS='APPEND')


      DO NT  = 1, NTANSAT

         ! For safety, initialize these up to LLGOS
         FORCE         = 0d0

         GC_CO2(:)       = 0d0
         MAP(:,:)        = 0d0

         ! Copy LGOS to make coding a bit cleaner

         ! Get grid box of current record
         IIJJ  = GET_IJ( REAL(TANSAT(NT)%LON(1),4), REAL(TANSAT(NT)%LAT(1),4), State_Grid_ADJ)
         I     = IIJJ(1)
         J     = IIJJ(2)
         print*, I, J
         ! obs exist, exit current loop
         IF ( OBS_EXIST(I, J) == 1) CYCLE

         OBS_EXIST(I, J) = 1
   
         ! adjoint verification
!         IF ( LCONV_FINITE_DIFF .OR. LTRAN_FINITE_DIFF ) THEN

!           print*, 'huoxiao debug append'
!           CALL APPEND_IJ( I, J )

!         END IF
!         IF (FINITE_DIFF) THEN
!            State_Chm_Adj%Species(I, J, 1, 1) = State_Chm_Adj%Species(I, J, 1, 1)-0.000001
!         ENDIF
         
         ! Get GC pressure levels (bar)
!         DO L = 1, State_Grid_Adj%NZ
! huoxiao debug
!            GC_PRES(L) = State_Met_Adj%PMID(I,J,L)*100
!         ENDDO
	
         ! Get GC surface pressure (bar)
!         GC_PSURF = State_Met_Adj%PSC2WET(I,J)*100
         
         ! Calculate the interpolation weight matrix
!         MAP(1:State_Grid_Adj%NZ,1:LOCO2)                                      &
!           = GET_INTMAP( State_Grid_Adj%NZ, GC_PRES(:),           GC_PSURF,   &
!                         LOCO2,  OCO2(NT)%PRES(1:LOCO2), OCO2(NT)%PSURF(1)  )


         ! Get CO2 values at native model resolution
         GC_CO2_NATIVE(:) = State_Chm_Adj%Species(I,J,:,IDTCO2)


         ! Get GC pressure levels (bar)
         DO L = 1, State_Grid_Adj%NZ
            ! huoxiao debug
            GC_PRES(L) = State_Met_Adj%PMID(I,J,L)*100
         ENDDO
         ! Get GC surface pressure (bar)
         GC_PSURF = State_Met_Adj%PSC2WET(I,J)*100

         ! Calculate the interpolation weight matrix
         MAP(1:State_Grid_Adj%NZ,1:LTANSAT)                                      &
           = GET_INTMAP( State_Grid_Adj%NZ, GC_PRES(:),           GC_PSURF,   &
                         LTANSAT,  TANSAT(NT)%PRES(1:LTANSAT), TANSAT(NT)%PSURF(1)  )

         !intepolation
         DO LL = 1, LTANSAT
            GC_CO2(LL) = 0d0
            DO L = 1, State_Grid_Adj%NZ
               GC_CO2(LL) = GC_CO2(LL)       &
                          + MAP(L,LL) * GC_CO2_NATIVE(L)
            ENDDO
         ENDDO
         ! pressure weighting function
!         h = pressure_weighting_function(TANSAT(NT)%PRES(:), &
!                                         TANSAT(NT)%PSURF(1), &
!                                         TANSAT(NT)%SPECIFIC_HUMUDITY(:), &
!                                         LTANSAT)
         h = TANSAT(NT)%PWF(:)
         XCO2_PERT = 0
         !sum(h_j*a_j*(x_m-x_a)_j)
         DO L = 1, LTANSAT

            XCO2_PERT = XCO2_PERT + h(L)*TANSAT(NT)%AVG_KERNEL(L)*&
                        (GC_CO2(L)-TANSAT(NT)%CO2A(L))
         END DO   
         GC_XCO2 = TANSAT(NT)%XCO2A+XCO2_PERT
         !FORCE
         DIFF = GC_XCO2-TANSAT(NT)%XCO2
         FORCE = TANSAT(NT)%S_OER_INV*DIFF
         NEW_COST(NT) = 0.5*FORCE*DIFF
         ! Get GC pressure levels (bar)
         ! Convert from kg/kg to v/v
         !GC_CO2_NATIVE(:) = GC_CO2_NATIVE(:) * TCVV_CO2
         ! Interpolate GC CO2 column to TES grid
!         DO LL = 1, LOCO2
!            GC_CO2(LL) = 0d0
!            DO L = 1, State_Grid_Adj%NZ
!               GC_CO2(LL) = GC_CO2(LL)       &
!                          + MAP(L,LL) * GC_CO2_NATIVE(L)
!            ENDDO
!         ENDDO
         !--------------------------------------------------------------
         ! Apply OCO2 observation operator
         !
         !   x_hat = x_a + h*A_k ( x_m - x_a )
         !
         !  where
         !    x_hat = GC modeled column as seen by TES [vmr]
         !    x_a   = OCO2 apriori column               [vmr]
         !    x_m   = GC modeled column                [vmr]
         !    A_k   = OCO2 averaging kernel
         !    h     = pressure_weighting_function
         !--------------------------------------------------------------
         ! x_m - x_a
!         DO L = 1, LOCO2
!           CO2_PERT(L) = GC_CO2(L) - OCO2(NT)%PRIOR(L)
!         ENDDO
         ! x_a + A_k * ( x_m - x_a )
!         DO L = 1, LOCO2
!            CO2_HAT(L)    = 0d0
!            DO LL = 1, LOCO2
!               CO2_HAT(L) = CO2_HAT(L)               &
!                          + OCO2(NT)%AVG_KERNEL(L,LL) * CO2_PERT(LL)
!            ENDDO
!            CO2_HAT(L)    = CO2_HAT(L) + OCO2(NT)%PRIOR(L)

!         ENDDO
         !--------------------------------------------------------------
         ! Calculate cost function, given S is error in vmr
         ! J = 1/2 [ model - obs ]^T S_{obs}^{-1} [ model - obs ]
         !--------------------------------------------------------------
         ! Calculate difference between modeled and observed profile
!         DO L = 1, LOCO2
!            DIFF(L) = CO2_HAT(L) - OCO2(NT)%CO2(L)
!         ENDDO
         ! Calculate 1/2 * DIFF^T * S_{obs}^{-1} * DIFF
!         DO L = 1, LOCO2
!            FORCE(L)     = 0d0
!            DO LL = 1, LOCO2
!               FORCE(L)  = FORCE(L) + OCO2(NT)%S_OER_INV(L,LL) * DIFF(LL)
!            ENDDO
!            NEW_COST(NT) = NEW_COST(NT) + 0.5d0 * DIFF(L) * FORCE(L)
!         ENDDO
         !--------------------------------------------------------------
         ! Begin adjoint calculations
         !--------------------------------------------------------------
         ! The adjoint forcing  is S_{obs}^{-1} * DIFF = FORCE
         DIFF_ADJ = FORCE
         XCO2_PERT_ADJ = DIFF_ADJ
         GC_CO2_ADJ = 0.d0
         DO L = 1, LTANSAT

            GC_CO2_ADJ(L) = h(L)*TANSAT(NT)%AVG_KERNEL(L)*XCO2_PERT_ADJ &
                                   +GC_CO2_ADJ(L)   
         END DO
        
         ! adjoint of interpolation
         DO L  = 1, State_Grid_Adj%NZ
            GC_CO2_NATIVE_ADJ(L) = 0d0
            DO LL = 1, LTANSAT
               GC_CO2_NATIVE_ADJ(L) = GC_CO2_NATIVE_ADJ(L)       &
                                   + MAP(L,LL) * GC_CO2_ADJ(LL)
            ENDDO
         ENDDO

         ! Adjoint of difference
!         DO L = 1, LOCO2
!            CO2_HAT_ADJ(L) =  DIFF_ADJ(L)
!         ENDDO
         ! adjoint of TES operator
!         DO L  = 1, LOCO2
!            CO2_PERT_ADJ(L)    = 0d0
!            DO LL = 1, LOCO2
!               CO2_PERT_ADJ(L) = CO2_PERT_ADJ(L)           &
!                               + OCO2(NT)%AVG_KERNEL(LL,L) &
!                               * CO2_HAT_ADJ(LL)
!           ENDDO
!         ENDDO
         ! huoxiao debug
         ! Adjoint of x_m - x_a
!         DO L = 1, LOCO2
           ! adj code:
!           GC_CO2_ADJ(L) = CO2_PERT_ADJ(L)
!         ENDDO

         ! adjoint of interpolation
!         DO L  = 1, State_Grid_Adj%NZ
!            GC_CO2_NATIVE_ADJ(L) = 0d0
!            DO LL = 1, LOCO2
!               GC_CO2_NATIVE_ADJ(L) = GC_CO2_NATIVE_ADJ(L)       &
!                                   + MAP(L,LL) * GC_CO2_ADJ(LL)
!            ENDDO
         ! huoxiao debug
         ! Adjoint of unit conversion
         GC_CO2_NATIVE_ADJ(:) = GC_CO2_NATIVE_ADJ(:) / TCVV_CO2
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
 
!        ENDIF

      ENDDO  ! NT
!!$OMP END PARALLEL DO

      ! Update cost function
       COST_FUNC = COST_FUNC + SUM(NEW_COST(1:NTANSAT))
  
!      print*, ' Updated value of COST_FUNC = ', COST_FUNC
      print*, ' TANSAT contribution           = ', SUM(NEW_COST(1:NTANSAT))
      print*, 'TANSAT', NEW_COST(1:NTANSAT)
!      print*, 'huoxiao_debug jo_gradient', MAXVAL(State_Chm_Adj%SpeciesAdj(:,:,:,IDTCO2))
!      CLOSE(1995)
!      CLOSE(1996)
!      CLOSE(1997)

      ! Return to calling program
      END SUBROUTINE CALC_TANSAT_CO2_FORCE

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


!------------------------------------------------------------------------------

!------------------------------------------------------------------------------

      FUNCTION GET_INTMAP( LGC_TOP, GC_PRESC, GC_SURFP,     &
     &                     LTM_TOP, TM_PRESC, TM_SURFP  )   &
     &         RESULT      ( HINTERPZ )

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
 

       ! Loop over each pressure level of TM grid
       DO LTM = 1, LTM_TOP
          IF ( TM_PRESC(LTM) > GC_PRESC(1) ) THEN
             HINTERPZ(1,LTM)         = 1D0
             HINTERPZ(2:LGC_TOP,LTM) = 0D0
          ENDIF
       ENDDO


      ! Return to calling program
      END FUNCTION GET_INTMAP

!-----------------------------------------------------------------------
      SUBROUTINE CALC_OBSERVATION_OPERATOR_D &
        (  State_Chm_Adj, State_Grid_Adj, State_Met_Adj )
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
      USE ADJ_ARRAYS_MOD,     ONLY : N_CALC, COST_FUNC, OBS_PERTD

      USE LOGICAL_ADJ_MOD,    ONLY : FINITE_DIFF

      USE LOGICAL_ADJ_MOD,    ONLY : LTRAN_FINITE_DIFF

      USE LOGICAL_ADJ_MOD,    ONLY : LCONV_FINITE_DIFF

!      USE PRESSURE_MOD,       ONLY : GET_PCENTER, GET_PEDGE
      USE ADJ_ARRAYS_MOD,     ONLY : EXIST_OBS
      USE ADJ_ARRAYS_MOD,     ONLY : OBS_EXIST
      USE State_Chm_Adj_Mod, ONLY : ChmStateAdj

      USE State_Grid_Adj_Mod, ONLY : GrdStateAdj

      USE State_Met_Adj_Mod, ONLY : MetStateAdj

      TYPE(GRDSTATEADJ), INTENT(IN) :: state_grid_adj
      TYPE(CHMSTATEADJ), INTENT(INOUT) :: state_chm_adj
      TYPE(METSTATEADJ), INTENT(IN) :: state_met_adj
      CHARACTER(len=20) :: cha_nymd
      CHARACTER(len=20) :: cha_hour
      CHARACTER(len=20) :: cha_minute
      INTEGER :: nhms
      INTEGER :: nymd
      INTEGER :: hour
      INTEGER :: minute
      INTEGER :: ntansat
      INTEGER :: nt
      INTEGER :: iijj(2), i, j
      INTEGER :: l, ll
      REAL*8 :: result1
      REAL*8 :: gc_pres(state_grid_adj%nz)
      REAL*8 :: gc_psurf
      REAL*8 :: map(state_grid_adj%nz, maxlev)
      REAL*8 :: gc_co2_natived(state_grid_adj%nz)
      REAL*8 :: gc_co2d(maxlev)
      REAL*8 :: xco2_pertd
      REAL*8 :: h(ltansat)
      REAL*8 :: temp
      REAL*8 :: XCO2_PERT
      REAL*8 :: FORCE
      REAL*8 :: DIFF_ADJ
      REAL*8 :: DIFF
      REAL*8 :: GC_XCO2
      REAL*8 :: XCO2_PERT_ADJ
      REAL*8 :: GC_CO2(LTANSAT)
      REAL*8 :: GC_CO2_ADJ(LTANSAT)
      REAL*8 :: GC_CO2_NATIVE(State_Grid_Adj%NZ)
      REAL*8 :: GC_CO2_NATIVE_ADJ(State_Grid_Adj%NZ)
      CHARACTER(len=255) :: oco2_filename
      LOGICAL :: alive

      REAL*8 :: TMP_COST_FUNC
      REAL*8 :: TMP_XCO2_PERD

      TMP_XCO2_PERD = 0.0
      TMP_COST_FUNC=0.0
      nhms = GET_NHMS()
      nymd = GET_NYMD()
      hour = nhms/10000
      minute = (nhms-hour*10000)/100
      WRITE(cha_nymd, *) nymd
      WRITE(cha_hour, '(I2.2)') hour
      WRITE(cha_minute, '(I2.2)') minute
      oco2_filename = &
     &   '../OCO2_CO2_Litev9_ND_Observation_Error_Covariance/oco2_'//TRIM(&
!     &   '../temp/oco2_'//TRIM(&
     &   ADJUSTL(cha_nymd))//TRIM(ADJUSTL(cha_hour))//TRIM(ADJUSTL(cha_minute&
     &   ))//'00.nc'
     INQUIRE(file=trim(oco2_filename), exist=alive)
      IF (.NOT.alive) THEN
          PRINT*, ' No matching OCO2 CO2 obs for this hour'
          RETURN 
      ENDIF
         CALL READ_OCO2_CO2_OBS(TRIM(oco2_filename), ntansat)
         !tangent validation
!         OPEN(UNIT=1995, FILE='OBS_TANGENT', ACCESS='APPEND')
!         OPEN(UNIT=1996, FILE='OBS_ADJOINT', ACCESS='APPEND')
         ! end
      print*, 'NTANSAT', NTANSAT
         DO nt=1,ntansat
            iijj = GET_IJ(REAL(TANSAT(nt)%lon(1), 4), &
                   REAL(TANSAT(nt)%lat(1), 4), state_grid_adj)
            i = iijj(1)
            j = iijj(2)
            result1 = OBS_EXIST(i, j)
            print*, 'Observation, i, j', I, J
            IF (result1 .NE. 1) THEN
               obs_exist(i, j) = 1
               gc_co2_natived(:) = State_Chm_Adj%SpeciesAdj(i, j, :, idtco2) / &
                                   TCVV_CO2
!               !tangent validation 
!               gc_co2_natived(:) = 0.0
!               gc_co2_natived(:) = 1/TCVV_CO2
!               !end tangent
               DO l=1,state_grid_adj%nz
                  gc_pres(l) = state_met_adj%pmid(i, j, l)*100
               END DO
               gc_psurf = state_met_adj%psc2wet(i, j)*100
               map(1:state_grid_adj%nz, 1:ltansat) = GET_INTMAP(state_grid_adj%&
&              nz, gc_pres(:), gc_psurf, ltansat, TANSAT(nt)%pres(1:ltansat)&
&              , TANSAT(nt)%psurf(1))
               gc_co2d = 0.0_8
               DO ll=1,ltansat
                  gc_co2d(ll) = 0.0_8
                  DO l=1,state_grid_adj%nz
                  gc_co2d(ll) = gc_co2d(ll) + map(l, ll)*gc_co2_natived(l)
                  END DO
               END DO
               h = TANSAT(nt)%pwf(:)
               xco2_pertd = 0.0_8
               !sum(h_j*a_j*(x_m-x_a)_j)
               DO l=1,ltansat
                  temp = h(l)*TANSAT(nt)%avg_kernel(l)
                  xco2_pertd = xco2_pertd + temp*gc_co2d(l)
               END DO
!               !tangent validation 
!               WRITE(1995, *) xco2_pertd
!               !tangent end

               ! di = HM(xb)-yi
               DIFF = 0.0d0
!               IF (  ALIVE ) THEN
               ! Get CO2 values at native model resolution
               GC_CO2_NATIVE(:) = State_Chm_Adj%SpeciestB(I,J,:,IDTCO2)
               !intepolation
               DO LL = 1, LTANSAT
                  GC_CO2(LL) = 0.0d0
               DO L = 1, State_Grid_Adj%NZ
                  GC_CO2(LL) = GC_CO2(LL)       &
                             + MAP(L,LL) * GC_CO2_NATIVE(L)
               ENDDO
               ENDDO
               XCO2_PERT = 0
               DO L = 1, LTANSAT

                  XCO2_PERT = XCO2_PERT + h(L)*TANSAT(NT)%AVG_KERNEL(L)*&
                              (GC_CO2(L)-TANSAT(NT)%CO2A(L))
               END DO
               GC_XCO2 = TANSAT(NT)%XCO2A+XCO2_PERT
               DIFF = GC_XCO2-TANSAT(NT)%XCO2
!               ENDIF
               !S_obs{-1}*(HMUw+di)
               FORCE = TANSAT(NT)%S_OER_INV*(DIFF+xco2_pertd)
!               print*, 'huoxiao_debug diff, xco2_pertd, GC_XCO2,'//&
!               ' TANSAT(NT)%XCO2'&
!                        , DIFF, xco2_pertd, GC_XCO2, TANSAT(NT)%XCO2
               COST_FUNC = COST_FUNC+0.5*(DIFF+xco2_pertd)*FORCE
               OBS_PERTD(I, J) = xco2_pertd
               tmp_xco2_perd = xco2_pertd+TMP_XCO2_PERD
               TMP_COST_FUNC = TMP_COST_FUNC+0.5*DIFF*TANSAT(NT)%S_OER_INV*DIFF
!               print*, 'DIFF, xco2_pertd, xco2_perd+diff', DIFF, xco2_pertd, xco2_pertd+diff
               ! adjoint H^T*S_obs{-1}*(HMUw+di) 
               DIFF_ADJ = FORCE
               XCO2_PERT_ADJ = DIFF_ADJ
               GC_CO2_ADJ = 0.d0
!               !adjoint validation
!               XCO2_PERT_ADJ = 1.0
!               !adjoint end
               DO L = 1, LTANSAT

                  GC_CO2_ADJ(L) = h(L)*TANSAT(NT)%AVG_KERNEL(L)*XCO2_PERT_ADJ &
                                       +GC_CO2_ADJ(L)
               END DO
                       ! adjoint of interpolation
               DO L  = 1, State_Grid_Adj%NZ
                  GC_CO2_NATIVE_ADJ(L) = 0d0
                  DO LL = 1, LTANSAT
                     GC_CO2_NATIVE_ADJ(L) = GC_CO2_NATIVE_ADJ(L)       &
                                   + MAP(L,LL) * GC_CO2_ADJ(LL)
                  ENDDO
               ENDDO
               GC_CO2_NATIVE_ADJ(:) = GC_CO2_NATIVE_ADJ(:) / TCVV_CO2
!               !adjoint validation
!               WRITE(1996, *) GC_CO2_NATIVE_ADJ(1)
!               !adjoint end
               
               State_Chm_Adj%SpeciesAdj(I,J,:,IDTCO2)  = &
                  State_Chm_Adj%SpeciesAdj(I,J,:,IDTCO2) &
                                                    + GC_CO2_NATIVE_ADJ(:)


!               print*, 'SpeciesAdj', State_Chm_Adj%SpeciesAdj(I,J,:,IDTCO2)
 
             END IF
          END DO
          PRINT*, 'TMP_COST_FUNC', TMP_COST_FUNC, TMP_XCO2_PERD
!      print*, 'INCREMENT SPECIESADJ', State_Chm_Adj%SpeciesAdj(107, 59, :, 1)
!      print*, 'INCREMENT SPECIESADJ', State_Chm_Adj%SpeciesAdj(106, 64, :, 1) 
!      print*, 'INCREMENT SPECIESADJ', State_Chm_Adj%SpeciesAdj(105, 67, :, 1)

!          print*, 'MAX SpeciesAdj in Observation Operator', &
!                  MAXVAL(State_Chm_Adj%SpeciesAdj)
!          print*, 'COST_FUNC in OBSERVATION', COST_FUNC
          !close file
!          CLOSE(1995)
!          CLOSE(1996)
      END SUBROUTINE CALC_OBSERVATION_OPERATOR_D

!---------------------------------------------------------------------------------
 
      FUNCTION OBS_EXIT() RESULT(FLAG)

      CHARACTER(len=20) :: cha_nymd
      CHARACTER(len=20) :: cha_hour
      CHARACTER(len=20) :: cha_minute
      CHARACTER(len=255) :: oco2_filename
      INTEGER :: nhms
      INTEGER :: nymd
      INTEGER :: hour
      INTEGER :: minute
      LOGICAL :: FLAG, alive
      nhms = GET_NHMS()
      nymd = GET_NYMD()
      hour = nhms/10000
      minute = (nhms-hour*10000)/100
      WRITE(cha_nymd, *) nymd
      WRITE(cha_hour, '(I2.2)') hour
      WRITE(cha_minute, '(I2.2)') minute
      oco2_filename = &
     &   '../OCO2_CO2_Litev9_ND_Observation_Error_Covariance/oco2_'//TRIM(&
!     &   '../temp/oco2_'//TRIM(&
     &   ADJUSTL(cha_nymd))//TRIM(ADJUSTL(cha_hour))//TRIM(ADJUSTL(cha_minute&
     &   ))//'00.nc'
      INQUIRE(file=trim(oco2_filename), exist=alive)
      FLAG = alive

      END FUNCTION OBS_EXIT
      END MODULE TANSAT_CO2_MOD 

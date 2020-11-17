
      MODULE TCCON_CO2_MOD

      USE TIME_MOD
      USE PRECISION_MOD
      IMPLICIT NONE


      !=================================================================
      ! MODULE VARIABLES
      !=================================================================

      ! Parameters
      INTEGER, PARAMETER           :: MAXLEV = 71
      INTEGER, PARAMETER           :: MAXTCC = 1000
      INTEGER, PARAMETER           :: LTCC = 71
      ! Record to store data from each GOS obs
      TYPE TCC_CO2_OBS
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
      ENDTYPE TCC_CO2_OBS

      TYPE(TCC_CO2_OBS)                          :: TCC(MAXTCC)

      INTEGER, PARAMETER           :: IDTCO2   = 1
      ! Same thing for TCVV(IDTCO2)
      REAL*8,  PARAMETER           :: TCVV_CO2 = 28.97d0  / 44d0
      REAL*8,  PARAMETER           :: g = 9.8
      CONTAINS

      SUBROUTINE READ_TCC_CO2_OBS( FILENAME, NTCC )
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

      INTEGER, INTENT(INOUT)           :: NTCC
      CHARACTER(LEN=*), INTENT(IN)   :: FILENAME
      ! local variables
      INTEGER                          :: FID
      INTEGER                          :: N, INC, I

      INTEGER                          :: RC

      REAL*8                           :: TEMP(LTCC)

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


      WRITE(*,"(A50,/)") '-READ_OCO2_CO2_OBS: reading file: ', TRIM(FILENAME)

      CALL NC_OPEN( TRIM(FILENAME), fID )

      ! time length
      time_str1d(1) = 1
      time_cnt1d(1) = 1
      varname = 'length'
      CALL Ncrd(time_varrd_1di, FID, TRIM(varname), time_str1d, time_cnt1d, &
                             err_stop, stat)

      NTCC = time_varrd_1di(1)
      ! allocate arrays and read nc file
      !latitude
      ALLOCATE( lat_varrd_1dr(time_varrd_1di(1)), STAT = RC)
      lat_str1d(1) = 1
      lat_cnt1d(1) = time_varrd_1di(1)
      varname = 'latitude'
      CALL Ncrd(lat_varrd_1dr, FID, TRIM(varname), lat_str1d, lat_cnt1d,  &
                             err_stop, stat)

      !longitude
      ALLOCATE( lon_varrd_1dr(time_varrd_1di(1)), STAT = RC)
      lon_str1d(1) = 1
      lon_cnt1d(1) = time_varrd_1di(1)
      varname = 'longitude'
      CALL Ncrd(lon_varrd_1dr, FID, TRIM(varname), lon_str1d, lon_cnt1d,  &
                             err_stop, stat)

     !prior
      ALLOCATE( pri_varrd_2dr(LTCC, time_varrd_1di(1)), STAT = RC)
      pri_str2d(1) = 1
      pri_str2d(2) = 1
      pri_cnt2d(1) = LTCC
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
      ALLOCATE( pre_varrd_2dr(LTCC, time_varrd_1di(1)), STAT = RC)
      pre_str2d(1) = 1
      pre_str2d(2) = 1
      pre_cnt2d(1) = LTCC
      pre_cnt2d(2) = time_varrd_1di(1)
      varname = 'pressure'
      CALL Ncrd(pre_varrd_2dr, FID, TRIM(varname), pre_str2d, pre_cnt2d)

      !specific humudity
      ALLOCATE( pwf_varrd_2dr(LTCC, time_varrd_1di(1)), STAT = RC)
      pwf_str2d(1) = 1
      pwf_str2d(2) = 1
      pwf_cnt2d(1) = LTCC
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

     ALLOCATE( avk_varrd_2dr(LTCC, time_varrd_1di(1)), STAT = RC)
      avk_str2d(1) = 1
      avk_str2d(2) = 1
      avk_cnt2d(1) = LTCC
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
      DO N = 1, NTCC
         ! observation error covariance

         TCC(N)%S_OER_INV        = 1/(Sinv_varrd_1dr(N)*8)
         TCC(N)%LTANSAT          =    LTCC
         TCC(N)%LAT              =    lat_varrd_1dr(N)
         TCC(N)%LON              =    lon_varrd_1dr(N)

         TCC(N)%PSURF            =    psurf_varrd_1dr(N)
         TCC(N)%PRES(1:LTCC)    =    pre_varrd_2dr(:, N)

         TCC(N)%PWF(1:LTCC)    =    &
                                         pwf_varrd_2dr(:, N)

         TCC(N)%XCO2             =    xco2_varrd_1dr(N)
         TCC(N)%XCO2A            =    xco2_apri_varrd_1dr(N)
            ! convert to ppmv
         TCC(N)%CO2A(1:LTCC)   =    pri_varrd_2dr(:, N)
!         print*, 'oco2', N, OCO2(N)%PRIOR(20)
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

      END SUBROUTINE READ_TCC_CO2_OBS

!-----------------------------------------------------------------------------------
      SUBROUTINE CALC_TCC_CO2_FORCE( Cost_Func, State_Chm_Adj,   &
                                     State_Grid_Adj, State_Met_Adj )

      ! Reference to f90 modules
      USE ADJ_ARRAYS_MOD,     ONLY : N_CALC

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
      REAL*8                      :: GC_CO2(MAXLEV)
      REAL*8                      :: CO2_PERT(State_Grid_Adj%NZ)
      REAL*8                      :: FORCE
      REAL*8                      :: DIFF
      REAL*8                      :: NEW_COST(MAXTCC)
      REAL*8                      :: OLD_COST
      INTEGER                     :: NTCC

      REAL*8                      :: GC_CO2_NATIVE_ADJ(State_Grid_Adj%NZ)
      REAL*8                      :: GC_CO2_ADJ(MAXLEV)
      REAL*8                      :: DIFF_ADJ
      REAL*8                      :: GC_XCO2
      REAL*8                      :: XCO2_PERT

      REAL*8                      :: XCO2_PERT_ADJ
      REAL*8                      :: h(LTCC)

      !tccon pressure
      REAL*8                      :: p(LTCC)      

      REAL*8                      :: MAP(State_Grid_Adj%NZ,MAXLEV)

      LOGICAL, SAVE               :: FIRST = .TRUE.
      INTEGER                     :: IOS

      LOGICAL                     :: ALIVE
      CHARACTER(LEN=255)          :: EDGE_FILENAME
      CHARACTER(LEN=255)          :: TCC_FILENAME
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

      INTEGER                          :: OBS_EXIST(State_Grid_Adj%NX, State_Grid_Adj%NY) 
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
      TCC_FILENAME = "../TCC_CO2/tcc_"//TRIM(ADJUSTL(CHA_NYMD)) &
                        //TRIM(ADJUSTL(CHA_HOUR)) &
                        //TRIM(ADJUSTL(CHA_MINUTE))//"00.nc"
      INQUIRE( FILE = TRIM(TCC_FILENAME), EXIST = ALIVE )
      WRITE(*,*) 'HUOXIAO_FILENAME', TRIM(TCC_FILENAME)
      IF ( .NOT. ALIVE ) THEN
         print*, ' No matching TCCON CO2 obs for this hour'
         WRITE(*,*) 'ALIVE:', ALIVE
         EXIST_OBS = .FALSE.
        RETURN
      ENDIF


      ! read OCO2   
      CALL READ_TCC_CO2_OBS( TRIM(TCC_FILENAME), NTCC )

      DO NT  = 1, NTCC

         ! For safety, initialize these up to LLGOS
         FORCE         = 0d0

         GC_CO2(:)       = 0d0
         MAP(:,:)        = 0d0

         ! Copy LGOS to make coding a bit cleaner

         ! Get grid box of current record
         IIJJ  = GET_IJ( REAL(TCC(NT)%LON(1),4), REAL(TCC(NT)%LAT(1),4), State_Grid_ADJ)
         I     = IIJJ(1)
         J     = IIJJ(2)



         ! obs exist, exit current loop
         IF ( OBS_EXIST(I, J) == 1) CYCLE

         OBS_EXIST(I, J) = 1

         ! adjoint verification
!         IF ( LCONV_FINITE_DIFF .OR. LTRAN_FINITE_DIFF ) THEN

!           print*, 'huoxiao debug append'
!           CALL APPEND_IJ( I, J )

!         END IF
         IF (FINITE_DIFF) THEN
            State_Chm_Adj%Species(I, J, 1, 1) = State_Chm_Adj%Species(I, J, 1, 1)+0.1
         ENDIF

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
         MAP(1:State_Grid_Adj%NZ,1:LTCC)                                      &
           = GET_INTMAP( State_Grid_Adj%NZ, GC_PRES(:),           GC_PSURF,   &
                         LTCC,  TCC(NT)%PRES(1:LTCC), TCC(NT)%PSURF(1)  )
         !intepolation
         DO LL = 1, LTCC
            GC_CO2(LL) = 0d0
            DO L = 1, State_Grid_Adj%NZ
               GC_CO2(LL) = GC_CO2(LL)       &
                          + MAP(L,LL) * GC_CO2_NATIVE(L)
            ENDDO
         ENDDO


         !--------------------------------------------------------------
         ! Apply TCCON observation operator
         !
         !   Xch4_m = Xch4_a + SUM_j( h_j * a_j * (x_m - x_a) )
         !
         !   Xch4_m  - model XCH4
         !   Xch4_a  - apriori XCH4 = h^T * x_a
         !   h       - pressure weighting function
         !   a       - column averaging kernel
         !   x_m     - model CH4 [v/v]
         !   x_a     - apriori CH4 [v/v]
         !
         !   The pressure weighting function is defined in Connor et al. 2008
         !     and the OCO-2 ATBD
         !--------------------------------------------------------------

         ! pressure weighting function
         p(1:LTCC)     = TCC(NT)%PRES(1:LTCC)
         L = 1
         h(L) = 1./TCC(NT)%PSURF(1) * ABS(                               &
               ( -1e0*p(L) + ( (p(L+1)-p(L))/(LOG(p(L+1)/p(L))) ) ) )
         L = LTCC
         h(L) = 1./TCC(NT)%PSURF(1) * ABS(                               &
               (  p(L) - ( (p(L)-p(L-1))/(LOG(p(L)/p(L-1))) ) ) )
         DO L=2,LTCC-1
            h(L) = 1./TCC(NT)%PSURF(1) * ABS(                            &
               ( -1e0*p(L) + ( (p(L+1)-p(L))/(LOG(p(L+1)/p(L))) ) ) +    &
               (      p(L) - ( (p(L)-p(L-1))/(LOG(p(L)/p(L-1))) ) )   )
         ENDDO


         XCO2_PERT = 0
         !sum(h_j*a_j*(x_m-x_a)_j)
         DO L = 1, LTCC

            XCO2_PERT = XCO2_PERT + h(L)*TCC(NT)%AVG_KERNEL(L)*&
                        (GC_CO2(L)-TCC(NT)%CO2A(L))
         END DO
         GC_XCO2 = TCC(NT)%XCO2A+XCO2_PERT
         !FORCE
         DIFF = GC_XCO2-TCC(NT)%XCO2
         FORCE = TCC(NT)%S_OER_INV*DIFF
         NEW_COST(NT) = NEW_COST(NT)+0.5*FORCE*DIFF

         !--------------------------------------------------------------
         ! Begin adjoint calculations
         !--------------------------------------------------------------
         ! The adjoint forcing  is S_{obs}^{-1} * DIFF = FORCE
         DIFF_ADJ = FORCE
         XCO2_PERT_ADJ = DIFF_ADJ
         GC_CO2_ADJ = 0.d0
         DO L = 1, LTCC

            GC_CO2_ADJ(L) = h(L)*TCC(NT)%AVG_KERNEL(L)*XCO2_PERT_ADJ &
                                   +GC_CO2_ADJ(L)
         END DO

         ! adjoint of interpolation
         DO L  = 1, State_Grid_Adj%NZ
            GC_CO2_NATIVE_ADJ(L) = 0d0
            DO LL = 1, LTCC
               GC_CO2_NATIVE_ADJ(L) = GC_CO2_NATIVE_ADJ(L)       &
                                   + MAP(L,LL) * GC_CO2_ADJ(LL)
            ENDDO
         ENDDO
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

!        ENDIF

      ENDDO  ! NT
!!$OMP END PARALLEL DO

      ! Update cost function
       COST_FUNC = COST_FUNC + SUM(NEW_COST(1:NTCC))

!      print*, ' Updated value of COST_FUNC = ', COST_FUNC
      print*, ' TANSAT contribution           = ', SUM(NEW_COST(1:NTCC))
!      CLOSE(1995)
!      CLOSE(1996)
!      CLOSE(1997)

      ! Return to calling program
      END SUBROUTINE CALC_TCC_CO2_FORCE


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

     END MODULE TCCON_CO2_MOD



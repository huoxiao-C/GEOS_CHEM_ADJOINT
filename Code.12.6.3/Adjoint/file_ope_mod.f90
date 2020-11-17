
Module File_Ope_Mod
       
      USE Precision_mod
      USE NCDF_MOD
      USE History_Util_Mod, ONLY : UNDEFINED_INT
      IMPLICIT NONE
      PRIVATE

      PUBLIC  :: Read_Met_Data_For_Emis_Time
      PUBLIC  :: Read_Met_Data_For_Dyn_Time
      PUBLIC  :: Read_Met_Data_For_I3_Time
      PUBLIC  :: Read_Met_Data_For_A1_Time
      PUBLIC  :: WRITE_FILE, READ_FILE
      PUBLIC  :: Read_Grid_Data
      PUBLIC  :: Read_Co2_Flux
      PUBLIC  :: READ_GDT_FILE
      PUBLIC  :: READ_CFN_FILE 
      PUBLIC  :: READ_SF_FILE
      PUBLIC  :: MAKE_GDT_FILE
      PUBLIC  :: MAKE_CFN_FILE
      PUBLIC  :: MAKE_CFNJO_FILE
      PUBLIC  :: MAKE_SF_FILE
      PUBLIC  :: GET_GRADNT_FROM_ADJ
      PUBLIC  :: WRITE_FILE_INIT
      PUBLIC  :: NC_CONTAINER_INIT
      PUBLIC  :: WRITE_FILE_CLOSE
      PUBLIC  :: WRITE_SF_TO_RSTFILE
      PUBLIC  :: WRITE_SF_TO_SPCFILE
      PUBLIC  :: WRITE_EXIT_FLAG
      PUBLIC  :: GET_SPECIES0
      PUBLIC  :: GET_BEC
      PUBLIC  :: WRITE_N_CALC
      PUBLIC  :: READ_TOT_ITER
      PUBLIC  :: READ_OUT_ITER
      PUBLIC  :: WRITE_OBS_PERTD
      PUBLIC  :: READ_REG_PARAM
      INTERFACE WRITE_FILE
         MODULE PROCEDURE WRITE_FILE_INT_4D
         MODULE PROCEDURE WRITE_FILE_INT_3D
         MODULE PROCEDURE WRITE_FILE_R8_4D
         MODULE PROCEDURE WRITE_FILE_R8_3D
      END INTERFACE

      INTERFACE READ_FILE
         MODULE PROCEDURE READ_FILE_INT_4D
         MODULE PROCEDURE READ_FILE_INT_3D
         MODULE PROCEDURE READ_FILE_R8_4D
         MODULE PROCEDURE READ_FILE_R8_3D
      END INTERFACE

      TYPE, PUBLIC :: NC_FILE_CONTAINER

         INTEGER   :: fID 
         INTEGER   :: lonID
         INTEGER   :: latID
         INTEGER   :: levID
         INTEGER   :: timeID
         INTEGER   :: VarCt
 
      END TYPE

      TYPE(NC_FILE_CONTAINER)   :: NC_FILE_OPE_VAR
      PUBLIC                    :: NC_FILE_OPE_VAR

      CONTAINS

!--------------------------------------------------------

      SUBROUTINE NC_CONTAINER_INIT()

      Nc_File_Ope_Var%fID = UNDEFINED_INT

      END SUBROUTINE NC_CONTAINER_INIT

      SUBROUTINE Read_Met_Data_For_Emis_Time( State_Grid_Adj, State_Met_Adj )
      ! read data for observation
     ! USE 
      USE State_Grid_Adj_Mod,     ONLY : GrdStateAdj
      USE State_Met_Adj_Mod,      ONLY : MetStateAdj
      USE Ncdf_Mod
      USE m_netcdf_io_open,       ONLY : Ncop_Rd
      USE m_netcdf_io_close,      ONLY : Nccl
      USE m_netcdf_io_checks,     ONLY : Ncdoes_Var_Exist
      USE m_netcdf_io_read,       ONLY : NcRd

      USE Time_Mod,               Only : Get_Nymd
      USE Time_Mod,               Only : GET_Nhms
      USE Time_Mod,               Only : GET_Hour
      USE Time_Mod,               Only : GET_Minute
      USE Time_Mod,               Only : GET_Second

      TYPE(GrdStateAdj), INTENT(IN)    :: State_Grid_Adj  ! Grid State object
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(MetStateAdj), INTENT(INOUT) :: State_Met_Adj   ! Obj for meteorology state

      INTEGER                          ::    NYMD, NHMS
      INTEGER                        HOUR, MINUTE, SEC

      ! dimension 
      INTEGER                          :: NX, NY, NZ, Time
      ! PMID
      REAL*8, ALLOCATABLE              :: pmid_varrd_4dr(:, :, :, :)
      INTEGER                          :: pmid_str4d(4)
      INTEGER                          :: pmid_cnt4d(4)
      ! psc2wet
      REAL*8, ALLOCATABLE              :: psc2wet_varrd_3dr(:, :, :)
      INTEGER                          :: psc2wet_str3d(3)
      INTEGER                          :: psc2wet_cnt3d(3)      

      LOGICAL                          :: err_stop
      INTEGER                          :: stat
      INTEGER                          :: RC
      CHARACTER(LEN=255)               :: varname
      INTEGER                          :: fID

      CHARACTER(LEN=255)              FILENAME
      CHARACTER(LEN=20)              TEMP_HOUR
      CHARACTER(LEN=20)              TEMP_MINUTE
      CHARACTER(LEN=20)              TEMP_SEC
      CHARACTER(LEN=20)              TEMP_NYMD
      CHARACTER(LEN=20)              TEMP_TIME

      print*, 'Met Data Is Reading For Emis Time'
       ! initialization 
      NX = State_Grid_Adj%NX
      NY = State_Grid_Adj%NY
      NZ = State_Grid_Adj%NZ
      Time = 1

      NYMD             =          GET_NYMD()
      HOUR             =          GET_HOUR()
      MINUTE           =          GET_MINUTE()
      SEC              =          GET_SECOND()
     !integer to string
      WRITE(TEMP_NYMD, *)  NYMD
      WRITE(TEMP_HOUR, "(I2.2)")  HOUR
      WRITE(TEMP_MINUTE, "(I2.2)") MINUTE
      WRITE(TEMP_SEC, "(I2.2)") SEC

      TEMP_TIME =        TRIM(ADJUSTL(TEMP_NYMD))  //     &
                          "_"                       //     &
                          TRIM(ADJUSTL(TEMP_HOUR))  //     &
                          TRIM(ADJUSTL(TEMP_MINUTE))//     &
                          "z.nc4"
      FileName = "./FdCo2Out/GEOSChem.StateMet."//TEMP_TIME
      CALL Ncop_Rd( fID, TRIM(FileName) )
      ! pmid
       ALLOCATE( pmid_varrd_4dr(NX, NY, NZ, Time), STAT = RC )
       varname = 'Met_PMID'
       pmid_str4d(1) = 1
       pmid_str4d(2) = 1
       pmid_str4d(3) = 1
       pmid_str4d(4) = 1

       pmid_cnt4d(1) = NX
       pmid_cnt4d(2) = NY
       pmid_cnt4d(3) = NZ
       pmid_cnt4d(4) = 1

       CALL NcRd(pmid_varrd_4dr, fID, TRIM(varname), pmid_str4d, pmid_cnt4d)

       State_Met_Adj%PMID(:, :, :) = pmid_varrd_4dr(:, :, :, Time)

       ! psc2wet
       ALLOCATE( psc2wet_varrd_3dr(NX, NY, Time), STAT = RC )
       varname = 'Met_PSC2WET'
       psc2wet_str3d(1) = 1
       psc2wet_str3d(2) = 1
       psc2wet_str3d(3) = 1

       psc2wet_cnt3d(1) = NX
       psc2wet_cnt3d(2) = NY
       psc2wet_cnt3d(3) = 1

       CALL NcRd(psc2wet_varrd_3dr, fID, TRIM(varname), psc2wet_str3d, psc2wet_cnt3d)

       State_Met_Adj%PSC2WET( :, : ) = psc2wet_varrd_3dr( :, :, Time )

      CALL Nccl(fID)
      ! deallocate 
      IF ( ALLOCATED( pmid_varrd_4dr ) )    DEALLOCATE( pmid_varrd_4dr )
      IF ( ALLOCATED( psc2wet_varrd_3dr ) ) DEALLOCATE( psc2wet_varrd_3dr )

      END SUBROUTINE Read_Met_Data_For_Emis_Time

!--------------------------------------------------------

      SUBROUTINE Read_Met_Data_For_Dyn_Time( State_Chm_Adj, State_Grid_Adj, &
                                             State_Met_Adj )

     ! USE 
      USE State_Chm_Adj_Mod,      Only : ChmStateAdj
      USE State_Grid_Adj_Mod,     ONLY : GrdStateAdj
      USE State_Met_Adj_Mod,      ONLY : MetStateAdj
      
      USE m_netcdf_io_read,       ONLY : NcRd
      USE m_netcdf_io_close,      ONLY : Nccl
      USE m_netcdf_io_open,       ONLY : Ncop_Rd

      USE Logical_Adj_Mod,        Only : LADJ_EMS
      USE Logical_Adj_Mod,        Only : LTRAN_ADJ
      USE Logical_Adj_Mod,        Only : LCONV_ADJ
      USE Logical_Adj_Mod,        Only : LPBL_ADJ

      USE Time_Mod,               Only : Get_Nymd
      USE Time_Mod,               Only : GET_Nhms
      USE Time_Mod,               Only : GET_Hour
      USE Time_Mod,               Only : GET_Minute
      USE Time_Mod,               Only : GET_Second

      TYPE(ChmStateAdj), INTENT(INOUT)       :: State_Chm_Adj
      TYPE(GrdStateAdj), INTENT(IN)    :: State_Grid_Adj  ! Grid State object
!
! !INPUT/OUTPUT PARAMETERS: 
!
      TYPE(MetStateAdj), INTENT(INOUT) :: State_Met_Adj   ! Obj for meteorology state

      INTEGER                          ::    NYMD, NHMS
      INTEGER                          ::    HOUR, MINUTE, SEC
      ! nc reading


      INTEGER                          :: RC
      CHARACTER(LEN=255)               :: varname
      INTEGER                          :: fID

      ! tmp pointer
      REAL(fp), POINTER                :: TMPAD(:, :, :, :)
      REAL(fp), POINTER                :: TMP_FPBL(:, :, :)
      REAL(fp), POINTER                :: TMP_DELP_DRY(:, :, :, :)
      REAL(fp), POINTER                :: TMP_CO2_CONV(:, :, :, :)
      REAL(fp), POINTER                :: TMP_CO2_ADV(:, :, :, :)
      REAL(fp), POINTER                :: TMP_QDELQ(:, :, :)
      REAL(fp), POINTER                :: TMP_QDELQ1(:, :, :)

      REAL(fp), POINTER                :: TMP_XMASS(:, :, :, :)
      REAL(fp), POINTER                :: TMP_YMASS(:, :, :, :)
      REAL(fp), POINTER                :: TMP_PSC2_DRY(:, :, :)
      REAL(fp), POINTER                :: TMP_PEDGE_DRY(:, :, :)
!      REAL(fp), POINTER                :: TMP_CO2_FLUX(:, :, :)
      INTEGER , POINTER                :: TMP_IMIX(:, :, :)

      ! dimension 
      INTEGER                          :: NX, NY, NZ, Time
      CHARACTER(LEN=255)              FILENAME, FILE_PREFIX
      CHARACTER(LEN=40)              TEMP_HOUR
      CHARACTER(LEN=40)              TEMP_MINUTE
      CHARACTER(LEN=40)              TEMP_SEC
      CHARACTER(LEN=40)              TEMP_NYMD
      CHARACTER(LEN=40)              TEMP_TIME

      print*, 'Met Data Is Reading is for Dyn Time'
       ! initialization 
      NX = State_Grid_Adj%NX
      NY = State_Grid_Adj%NY
      NZ = State_Grid_Adj%NZ
      !allocate
      ALLOCATE(TMPAD(NX, NY, NZ, 1), STAT = RC)
      ALLOCATE(TMP_FPBL(NX, NY, 1), STAT = RC)
      ALLOCATE(TMP_DELP_DRY(NX, NY, NZ, 1), STAT = RC)
      ALLOCATE(TMP_CO2_CONV(NX, NY, NZ, 1), STAT = RC)
      ALLOCATE(TMP_CO2_ADV(NX, NY, NZ, 1), STAT = RC)
      ALLOCATE(TMP_QDELQ(NX, NY, NZ), STAT = RC)
      ALLOCATE(TMP_QDELQ1(NX, NY, NZ), STAT = RC)
      ALLOCATE(TMP_XMASS(NX, NY, NZ, 1), STAT = RC)
      ALLOCATE(TMP_YMASS(NX, NY, NZ, 1), STAT = RC)
      ALLOCATE(TMP_PSC2_DRY(NX, NY, 1), STAT = RC)
      ALLOCATE(TMP_PEDGE_DRY(NX, NY, 1), STAT = RC)
!      ALLOCATE(TMP_CO2_FLUX(NX, NY,  1), STAT = RC)
      ALLOCATE(TMP_IMIX(NX, NY, 1), STAT = RC)


      NYMD             =          GET_NYMD()
      NHMS             =          GET_NHMS()
      HOUR             =          GET_HOUR()
      MINUTE           =          GET_MINUTE()
      SEC              =          GET_SECOND()
!integer to string
      WRITE(TEMP_NYMD, *)  NYMD
      WRITE(TEMP_HOUR, "(I2.2)")  HOUR
      WRITE(TEMP_MINUTE, "(I2.2)") MINUTE
      WRITE(TEMP_SEC, "(I2.2)") SEC
      FILENAME = 'AdjParameters/GEOSChem.AdjParameters.' //     &
                               TRIM(ADJUSTL(TEMP_NYMD))  //     &
                               '_'                       //     &
                               TRIM(ADJUSTL(TEMP_HOUR))  //     &
                               TRIM(ADJUSTL(TEMP_MINUTE))//     &
                               TRIM(ADJUSTL(TEMP_SEC))   //     &
                               '.nc4'
      print*, 'Read File: ', FILENAME
      CALL Ncop_Rd( fID, TRIM(FileName) )

     
      ! imix, fpbl, ad for pbl
      IF( LPBL_ADJ .OR. LTRAN_ADJ ) THEN

      VarName = 'ad'
      CALL READ_FILE( fID, VarName, TMPAD )
      State_Met_Adj%AD(:, :, :) = TMPAD(:, :, :, 1)

      ENDIF

      IF( LPBL_ADJ ) THEN
      VarName = 'fpbl'
      CALL READ_FILE( fID, VarName, TMP_FPBL )
      State_Met_Adj%FPBL(:, :) = TMP_FPBL(:, :, 1)

      VarName = 'imix'
      CALL READ_FILE( fID, VarName, TMP_IMIX )   
      State_Met_Adj%IMIX(:, :) = TMP_IMIX(:, :, 1)

      ENDIF

      ! delp_dry for convection
      IF( LCONV_ADJ .OR. LADJ_EMS) THEN

      VarName = 'delp_dry'
      CALL READ_FILE( fID, VarName, TMP_DELP_DRY )
      State_Met_Adj%DELP_DRY(:, :, :) = TMP_DELP_DRY(:, :, :, 1)

      VarName = 'co2_conv'
      CALL READ_FILE( fID, VarName, TMP_CO2_CONV )
      State_Chm_Adj%CO2_CONV(:, :, :, :) = TMP_CO2_CONV(:, :, :, :)

      VarName = 'qdelq_flag_nt1'
      CALL READ_FILE( fID, VarName, TMP_QDELQ )
      State_Chm_Adj%QDELQ(:, :, :, 1) = TMP_QDELQ(:, :, :)

      VarName = 'qdelq_flag_nt2'
      CALL READ_FILE( fID, VarName, TMP_QDELQ )
      State_Chm_Adj%QDELQ(:, :, :, 2) = TMP_QDELQ(:, :, :)

      VarName = 'qdelq_flag1_nt1'
      CALL READ_FILE( fID, VarName, TMP_QDELQ1 )
      State_Chm_Adj%QDELQ1(:, :, :, 1) = TMP_QDELQ1(:, :, :)

      VarName = 'qdelq_flag1_nt2'
      CALL READ_FILE( fID, VarName, TMP_QDELQ1 )
      State_Chm_Adj%QDELQ1(:, :, :, 2) = TMP_QDELQ1(:, :, :)
      ENDIF

      !xmass ymass psc2_dry pedge_dry for advection
      IF( LTRAN_ADJ ) THEN

      VarName = 'xmass'
      CALL READ_FILE( fID, VarName, TMP_XMASS )
      State_Met_Adj%XMASS(:, :, :) = TMP_XMASS(:, :, :, 1)
      VarName = 'ymass'
      CALL READ_FILE( fID, VarName, TMP_YMASS )
      State_Met_Adj%YMASS(:, :, :) = TMP_YMASS(:, :, :, 1)

      VarName = 'psc2_dry'
      CALL READ_FILE( fID, VarName, TMP_PSC2_DRY )
      State_Met_Adj%PSC2_DRY(:, :) = TMP_PSC2_DRY(:, :, 1)

      VarName = 'pedge_dry'
      CALL READ_FILE( fID, VarName, TMP_PEDGE_DRY )
      State_Met_Adj%PEDGE_DRY(:, :) = TMP_PEDGE_DRY(:, :, 1)

      VarName = 'co2_adv'
      CALL READ_FILE( fID, VarName, TMP_CO2_ADV )
      State_Chm_Adj%CO2_ADV(:, :, :, :) = TMP_CO2_ADV(:, :, :, :)

      ENDIF

      CALL Nccl(fID)

      ! deallocate 
      IF ( ASSOCIATED( TMPAD ) )        DEALLOCATE( TMPAD )
      IF ( ASSOCIATED( TMP_FPBL ) )     DEALLOCATE( TMP_FPBL )
      IF ( ASSOCIATED( TMP_IMIX ) )     DEALLOCATE( TMP_IMIX )
      IF ( ASSOCIATED( TMP_DELP_DRY ) ) DEALLOCATE( TMP_DELP_DRY )
      IF ( ASSOCIATED( TMP_CO2_CONV ) ) DEALLOCATE( TMP_CO2_CONV )
      IF ( ASSOCIATED( TMP_CO2_ADV ) ) DEALLOCATE( TMP_CO2_ADV )
      IF ( ASSOCIATED( TMP_QDELQ ) ) DEALLOCATE( TMP_QDELQ )
      IF ( ASSOCIATED( TMP_QDELQ1 ) ) DEALLOCATE( TMP_QDELQ1 )
      IF ( ASSOCIATED( TMP_XMASS ) )    DEALLOCATE( TMP_XMASS )
      IF ( ASSOCIATED( TMP_YMASS ) )    DEALLOCATE( TMP_YMASS )
      IF ( ASSOCIATED( TMP_PSC2_DRY ) ) DEALLOCATE( TMP_PSC2_DRY )
      IF ( ASSOCIATED( TMP_PEDGE_DRY ) )DEALLOCATE( TMP_PEDGE_DRY )

      TMPAD          => NULL()
      TMP_FPBL       => NULL()
      TMP_IMIX       => NULL()
      TMP_DELP_DRY   => NULL()
      TMP_CO2_CONV   => NULL()
      TMP_CO2_ADV    => NULL()
      TMP_QDELQ      => NULL()
      TMP_QDELQ1     => NULL()
      TMP_XMASS      => NULL()
      TMP_YMASS      => NULL()
      TMP_PSC2_DRY   => NULL()
      TMP_PEDGE_DRY  => NULL()
!     IF ( ASSOCIATED( TMP_CO2_FLUX ) )    DEALLOCATE( TMP_CO2_FLUX )

      END SUBROUTINE Read_Met_Data_For_Dyn_Time
!--------------------------------------------------------
!--------------------------------------------------------
      SUBROUTINE Read_Co2_Flux( State_Grid_Adj, State_Met_Adj )


      USE State_Grid_Adj_Mod,     ONLY : GrdStateAdj
      USE State_Met_Adj_Mod,      ONLY : MetStateAdj

      USE m_netcdf_io_read,       ONLY : NcRd
      USE m_netcdf_io_close,      ONLY : Nccl
      USE m_netcdf_io_open,       ONLY : Ncop_Rd

      USE Logical_Adj_Mod,        Only : LADJ_EMS
      USE Logical_Adj_Mod,        Only : LADJ_FEMISS
      USE Logical_Adj_Mod,        Only : LADJ_FOSSIL
      USE Logical_Adj_Mod,        Only : LADJ_OCEAN
      USE Logical_Adj_Mod,        Only : LADJ_BIO
      USE Logical_Adj_Mod,        Only : LADJ_BIOFUEL


      USE Time_Mod,               Only : Get_Nymd
      USE Time_Mod,               Only : GET_Nhms
      USE Time_Mod,               Only : GET_Hour
      USE Time_Mod,               Only : GET_Minute
      USE Time_Mod,               Only : GET_Second


      TYPE(GrdStateAdj), INTENT(IN)    :: State_Grid_Adj  ! Grid State object
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(MetStateAdj), INTENT(INOUT) :: State_Met_Adj   ! Obj for meteorology state

      INTEGER                          ::    NYMD, NHMS
      INTEGER                          ::    HOUR, MINUTE, SEC

      INTEGER                          :: RC
      CHARACTER(LEN=255)               :: varname
      INTEGER                          :: fID

      REAL(fp), POINTER                :: TMP_CO2_FLUX(:, :, :)

      ! dimension
      INTEGER                          :: NX, NY, NZ, Time
      CHARACTER(LEN=255)              FILENAME, FILE_PREFIX
      CHARACTER(LEN=40)              TEMP_HOUR
      CHARACTER(LEN=40)              TEMP_MINUTE
      CHARACTER(LEN=40)              TEMP_SEC
      CHARACTER(LEN=40)              TEMP_NYMD
      CHARACTER(LEN=40)              TEMP_TIME

      NX = State_Grid_Adj%NX
      NY = State_Grid_Adj%NY
      NZ = State_Grid_Adj%NZ

      ALLOCATE(TMP_CO2_FLUX(NX, NY,  1), STAT = RC)

      NYMD             =          GET_NYMD()
      HOUR             =          GET_HOUR()
      MINUTE           =          GET_MINUTE()
      SEC              =          GET_SECOND()
!integer to string
      WRITE(TEMP_NYMD, *)  NYMD
      WRITE(TEMP_HOUR, "(I2.2)")  HOUR
      WRITE(TEMP_MINUTE, "(I2.2)") MINUTE
      WRITE(TEMP_SEC, "(I2.2)") SEC
      FILENAME = 'AdjParameters/GEOSChem.AdjParameters.' //     &
                               TRIM(ADJUSTL(TEMP_NYMD))  //     &
                               '_'                       //     &
                               TRIM(ADJUSTL(TEMP_HOUR))  //     &
                               TRIM(ADJUSTL(TEMP_MINUTE))//     &
                               TRIM(ADJUSTL(TEMP_SEC))   //     &
                               '.nc4'
      CALL Ncop_Rd( fID, TRIM(FileName) )

      IF( LADJ_EMS ) THEN

        IF ( LADJ_FOSSIL .OR. LADJ_FEMISS ) THEN
        VarName = 'co2_scaled_fossil_flux'
        CALL READ_FILE( fID, VarName, TMP_CO2_FLUX )
        State_Met_Adj%CO2_FOSSIL_FLUX(:, :) = TMP_CO2_FLUX(:, :, 1)
        ENDIF

        IF ( LADJ_OCEAN .OR. LADJ_FEMISS ) THEN
        VarName = 'co2_scaled_ocean_flux'
        CALL READ_FILE( fID, VarName, TMP_CO2_FLUX )
        State_Met_Adj%CO2_OCEAN_FLUX(:, :) = TMP_CO2_FLUX(:, :, 1)
        ENDIF

        IF ( LADJ_BIO .OR. LADJ_FEMISS ) THEN
        VarName = 'co2_scaled_bio_flux'
        CALL READ_FILE( fID, VarName, TMP_CO2_FLUX )
        State_Met_Adj%CO2_BIO_FLUX(:, :) = TMP_CO2_FLUX(:, :, 1)
!        print*, 'huoxiao_debug bio flux', MAXVAL(TMP_CO2_FLUX(:, :, 1)), &
!                                          MINVAL(TMP_CO2_FLUX(:, :, 1))
        ENDIF

        IF ( LADJ_BIOFUEL .OR. LADJ_FEMISS ) THEN
        VarName = 'co2_scaled_biofuel_flux'
        CALL READ_FILE( fID, VarName, TMP_CO2_FLUX )
        State_Met_Adj%CO2_BIOFUEL_FLUX(:, :) = TMP_CO2_FLUX(:, :, 1)
        ENDIF
      ENDIF
      CALL Nccl(fID)

      IF ( ASSOCIATED( TMP_CO2_FLUX ) )    DEALLOCATE( TMP_CO2_FLUX )

     END SUBROUTINE Read_Co2_Flux
!--------------------------------------------------------

      SUBROUTINE Read_Met_Data_For_I3_Time( State_Grid_Adj, State_Met_Adj )

           ! USE 
      USE State_Grid_Adj_Mod,     ONLY : GrdStateAdj
      USE State_Met_Adj_Mod,      ONLY : MetStateAdj
      USE Ncdf_Mod
      USE m_netcdf_io_checks,     ONLY : Ncdoes_Var_Exist
      USE m_netcdf_io_read,       ONLY : NcRd
      USE m_netcdf_io_close,      ONLY : Nccl
      USE m_netcdf_io_open,       ONLY : Ncop_Rd
   
      USE Logical_Adj_Mod,        Only : LCONV_ADJ

      USE Time_Mod,               Only : Get_Nymd
      USE Time_Mod,               Only : Get_Nhms
      USE Time_Mod,               Only : GET_Nhms
      USE Time_Mod,               Only : GET_Hour
      USE Time_Mod,               Only : GET_Minute
      USE Time_Mod,               Only : GET_Second

      TYPE(GrdStateAdj), INTENT(IN)    :: State_Grid_Adj  ! Grid State object
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(MetStateAdj), INTENT(INOUT) :: State_Met_Adj   ! Obj for meteorology state

      INTEGER                          ::    NYMD, NHMS
      INTEGER                          ::    HOUR, MINUTE, SEC
      ! dimension 
      INTEGER                          :: NX, NY, NZ, Time

      ! cmfmc
      REAL*8, ALLOCATABLE              :: cmfmc_varrd_4dr(:, :, :, :)
      INTEGER                          :: cmfmc_str4d(4)
      INTEGER                          :: cmfmc_cnt4d(4)

      ! pflcu
      REAL*8, ALLOCATABLE              :: pflcu_varrd_4dr(:, :, :, :)
      INTEGER                          :: pflcu_str4d(4)
      INTEGER                          :: pflcu_cnt4d(4)

      ! pficu
      REAL*8, ALLOCATABLE              :: pficu_varrd_4dr(:, :, :, :)
      INTEGER                          :: pficu_str4d(4)
      INTEGER                          :: pficu_cnt4d(4)

      ! dtrain
      REAL*8, ALLOCATABLE              :: dtrain_varrd_4dr(:, :, :, :)
      INTEGER                          :: dtrain_str4d(4)
      INTEGER                          :: dtrain_cnt4d(4)

      ! dqrcu
      REAL*8, ALLOCATABLE              :: dqrcu_varrd_4dr(:, :, :, :)
      INTEGER                          :: dqrcu_str4d(4)
      INTEGER                          :: dqrcu_cnt4d(4)

      LOGICAL                          :: err_stop
      INTEGER                          :: stat
      INTEGER                          :: RC
      CHARACTER(LEN=255)               :: varname
      INTEGER                          :: fID

      !time 
      CHARACTER(LEN=255)               :: FILENAME, FILE_PREFIX
      CHARACTER(LEN=20)                :: TEMP_HOUR
      CHARACTER(LEN=20)                :: TEMP_MINUTE
      CHARACTER(LEN=20)                :: TEMP_SEC
      CHARACTER(LEN=20)                :: TEMP_NYMD
      CHARACTER(LEN=20)                :: TEMP_TIME
      INTEGER                          :: Species


      IF( .NOT. LCONV_ADJ ) THEN

      print*, 'no convection adjoint, therefore there is no need to read these met data'
      return

      ENDIF

      print*, 'Met Data Is Reading For I3 Time'

       ! initialization 
      NX = State_Grid_Adj%NX
      NY = State_Grid_Adj%NY
      NZ = State_Grid_Adj%NZ
      Species = 1
      NYMD             =          GET_NYMD()
      NHMS             =          GET_NHMS()
      HOUR             =          GET_HOUR()
      MINUTE           =          GET_MINUTE()
      SEC              =          GET_SECOND()

     !integer to string
      WRITE(TEMP_NYMD, *)  NYMD
      WRITE(TEMP_HOUR, "(I2.2)")  HOUR
      WRITE(TEMP_MINUTE, "(I2.2)") MINUTE
      WRITE(TEMP_SEC, "(I2.2)") SEC

      TEMP_TIME =        TRIM(ADJUSTL(TEMP_NYMD))   //     &
                          '_'                       //     &
                          TRIM(ADJUSTL(TEMP_HOUR))  //     &
                          TRIM(ADJUSTL(TEMP_MINUTE))//     &     
                          'z.nc4'
      !pficu, pflcu, cmfmc
      FileName = "./FdCo2Out/GEOSChem.LevelEdgeDiags."//TEMP_TIME

      WRITE(*,100) 'pficu, pflcu, cmfmc', NYMD, NHMS
      WRITE(*,110) FileName

      CALL Ncop_Rd( fID, TRIM(FileName) )
       ! cmfmc
       ALLOCATE( cmfmc_varrd_4dr(NX, NY, NZ+1, Species), STAT = RC )
       varname = 'Met_CMFMC'
       cmfmc_str4d(1) = 1
       cmfmc_str4d(2) = 1
       cmfmc_str4d(3) = 1
       cmfmc_str4d(4) = 1

       cmfmc_cnt4d(1) = NX
       cmfmc_cnt4d(2) = NY
       cmfmc_cnt4d(3) = NZ+1
       cmfmc_cnt4d(4) = 1

       CALL NcRd(cmfmc_varrd_4dr, fID, TRIM(varname), cmfmc_str4d, cmfmc_cnt4d)

       State_Met_Adj%CMFMC( :, :, : ) = cmfmc_varrd_4dr( :, :, :, Species )

       ! pflcu
       ALLOCATE( pflcu_varrd_4dr(NX, NY, NZ+1, Species), STAT = RC )
       varname = 'Met_PFLCU'
       pflcu_str4d(1) = 1
       pflcu_str4d(2) = 1
       pflcu_str4d(3) = 1
       pflcu_str4d(4) = 1

       pflcu_cnt4d(1) = NX
       pflcu_cnt4d(2) = NY
       pflcu_cnt4d(3) = NZ+1
       pflcu_cnt4d(4) = 1

       CALL NcRd(pflcu_varrd_4dr, fID, TRIM(varname), pflcu_str4d, pflcu_cnt4d)

       State_Met_Adj%PFLCU( :, :, : ) = pflcu_varrd_4dr( :, :, :, Species )

       ! pficu
       ALLOCATE( pficu_varrd_4dr(NX, NY, NZ+1, Species), STAT = RC )
       varname = 'Met_PFICU'
       pficu_str4d(1) = 1
       pficu_str4d(2) = 1
       pficu_str4d(3) = 1
       pficu_str4d(4) = 1

       pficu_cnt4d(1) = NX
       pficu_cnt4d(2) = NY
       pficu_cnt4d(3) = NZ+1
       pficu_cnt4d(4) = 1

       CALL NcRd(pficu_varrd_4dr, fID, TRIM(varname), pficu_str4d, pficu_cnt4d)

       State_Met_Adj%PFICU( :, :, : ) = pficu_varrd_4dr( :, :, :, Species )

      ! close file
      CALL Nccl(fID)

      !met data for drqcu, dtrain

      FileName = "./FdCo2Out/GEOSChem.StateMet."//TEMP_TIME

      WRITE(*,100) 'drqcu, dtrain', NYMD, NHMS
      WRITE(*,110) FileName

      CALL Ncop_Rd( fID, TRIM(FileName) )


       ! dtrain
       ALLOCATE( dtrain_varrd_4dr(NX, NY, NZ, Species), STAT = RC )
       varname = 'Met_DTRAIN'
       dtrain_str4d(1) = 1
       dtrain_str4d(2) = 1
       dtrain_str4d(3) = 1
       dtrain_str4d(4) = 1

       dtrain_cnt4d(1) = NX
       dtrain_cnt4d(2) = NY
       dtrain_cnt4d(3) = NZ
       dtrain_cnt4d(4) = 1

       CALL NcRd(dtrain_varrd_4dr, fID, TRIM(varname), dtrain_str4d, dtrain_cnt4d)

       State_Met_Adj%DTRAIN( :, :, : ) = dtrain_varrd_4dr( :, :, :, Species )

       ! dqrcu
       ALLOCATE( dqrcu_varrd_4dr(NX, NY, NZ, Species), STAT = RC )
       varname = 'Met_DQRCU'
       dqrcu_str4d(1) = 1
       dqrcu_str4d(2) = 1
       dqrcu_str4d(3) = 1
       dqrcu_str4d(4) = 1

       dqrcu_cnt4d(1) = NX
       dqrcu_cnt4d(2) = NY
       dqrcu_cnt4d(3) = NZ
       dqrcu_cnt4d(4) = 1

       CALL NcRd(dqrcu_varrd_4dr, fID, TRIM(varname), dqrcu_str4d, dqrcu_cnt4d)

       State_Met_Adj%DQRCU( :, :, : ) = dqrcu_varrd_4dr( :, :, :, Species )

      ! close file
      CALL Nccl(fID)
100       FORMAT( '     - Reading file for ', a, '; reference = ',i8.8,1x,i6.6)
110       FORMAT( '        with filename = ', a                                )

      ! deallocate 
      IF ( ALLOCATED( cmfmc_varrd_4dr ) )    DEALLOCATE( cmfmc_varrd_4dr ) 
      IF ( ALLOCATED( pflcu_varrd_4dr ) )    DEALLOCATE( pflcu_varrd_4dr )
      IF ( ALLOCATED( pficu_varrd_4dr ) )    DEALLOCATE( pficu_varrd_4dr )
      IF ( ALLOCATED( dtrain_varrd_4dr ) )   DEALLOCATE( dtrain_varrd_4dr )
      IF ( ALLOCATED( dqrcu_varrd_4dr ) )    DEALLOCATE( dqrcu_varrd_4dr )


      END SUBROUTINE Read_Met_Data_For_I3_Time

!--------------------------------------------------------

      SUBROUTINE Read_Met_Data_For_A1_Time( State_Grid_Adj, State_Met_Adj )

           ! USE 
      USE State_Grid_Adj_Mod,     ONLY : GrdStateAdj
      USE State_Met_Adj_Mod,      ONLY : MetStateAdj
      USE Ncdf_Mod
      USE m_netcdf_io_checks,     ONLY : Ncdoes_Var_Exist
      USE m_netcdf_io_read,       ONLY : NcRd
      USE m_netcdf_io_close,      ONLY : Nccl
      USE m_netcdf_io_open,       ONLY : Ncop_Rd

      USE Logical_Adj_Mod,        Only : LTRAN_ADJ

      USE Time_Mod,               Only : Get_Nymd
      USE Time_Mod,               Only : GET_Nhms
      USE Time_Mod,               Only : GET_Hour
      USE Time_Mod,               Only : GET_Minute
      USE Time_Mod,               Only : GET_Second

      TYPE(GrdStateAdj), INTENT(IN)    :: State_Grid_Adj  ! Grid State object
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(MetStateAdj), INTENT(INOUT) :: State_Met_Adj   ! Obj for meteorology state

      INTEGER                          ::    NYMD, NHMS
      INTEGER                          ::    HOUR, MINUTE, SEC
      ! dimension 
      INTEGER                          :: NX, NY, NZ, Time

      ! u
      REAL*8, ALLOCATABLE              :: u_varrd_4dr(:, :, :, :)
      INTEGER                          :: u_str4d(4)
      INTEGER                          :: u_cnt4d(4)

      ! u
      REAL*8, ALLOCATABLE              :: v_varrd_4dr(:, :, :, :)
      INTEGER                          :: v_str4d(4)
      INTEGER                          :: v_cnt4d(4)

      LOGICAL                          :: err_stop
      INTEGER                          :: stat
      INTEGER                          :: RC
      CHARACTER(LEN=255)               :: varname
      INTEGER                          :: fID

      !time 
      CHARACTER(LEN=255)               :: FILENAME, FILE_PREFIX
      CHARACTER(LEN=20)                :: TEMP_HOUR
      CHARACTER(LEN=20)                :: TEMP_MINUTE
      CHARACTER(LEN=20)                :: TEMP_SEC
      CHARACTER(LEN=20)                :: TEMP_NYMD
      CHARACTER(LEN=20)                :: TEMP_TIME
      INTEGER                          :: Species


      IF( .NOT. LTRAN_ADJ ) THEN

      print*, 'no transportation adjoint, therefore there is no need to read these met data'
      return

      ENDIF


      print*, 'Met Data Is Reading For A1 Time'
       ! initialization 
      NX = State_Grid_Adj%NX
      NY = State_Grid_Adj%NY
      NZ = State_Grid_Adj%NZ
      Species = 1

      NYMD             =          GET_NYMD()
      HOUR             =          GET_HOUR()
      MINUTE           =          GET_MINUTE()
      SEC              =          GET_SECOND()
     !integer to string
      WRITE(TEMP_NYMD, *)  NYMD
      WRITE(TEMP_HOUR, "(I2.2)")  HOUR
      WRITE(TEMP_MINUTE, "(I2.2)") MINUTE
      WRITE(TEMP_SEC, "(I2.2)") SEC
      TEMP_TIME =        TRIM(ADJUSTL(TEMP_NYMD))  //      &
                          '_'                      //      &
                          TRIM(ADJUSTL(TEMP_HOUR))  //     &
                          TRIM(ADJUSTL(TEMP_MINUTE))//     &
                          'z.nc4'
      !u, v
      FileName = "./FdCo2Out/GEOSChem.StateMet."//TEMP_TIME
      CALL Ncop_Rd( fID, TRIM(FileName) )

       ! u
      ALLOCATE( u_varrd_4dr(NX, NY, NZ, Species), STAT = RC )
      varname = 'Met_U'
      u_str4d(1) = 1
      u_str4d(2) = 1
      u_str4d(3) = 1
      u_str4d(4) = 1

      u_cnt4d(1) = NX
      u_cnt4d(2) = NY
      u_cnt4d(3) = NZ
      u_cnt4d(4) = 1

      CALL NcRd(u_varrd_4dr, fID, TRIM(varname), u_str4d, u_cnt4d)

      State_Met_Adj%U( :, :, : ) = u_varrd_4dr( :, :, :, Species )

       ! u
      ALLOCATE( v_varrd_4dr(NX, NY, NZ, Species), STAT = RC )
      varname = 'Met_V'
      v_str4d(1) = 1
      v_str4d(2) = 1
      v_str4d(3) = 1
      v_str4d(4) = 1

      v_cnt4d(1) = NX
      v_cnt4d(2) = NY
      v_cnt4d(3) = NZ
      v_cnt4d(4) = 1

      CALL NcRd(v_varrd_4dr, fID, TRIM(varname), v_str4d, v_cnt4d)

      State_Met_Adj%V( :, :, : ) = v_varrd_4dr( :, :, :, Species )

      ! close file
      CALL Nccl(fID)
      ! deallocate 
      IF ( ALLOCATED( u_varrd_4dr ) )    DEALLOCATE( u_varrd_4dr )
      IF ( ALLOCATED( v_varrd_4dr ) )    DEALLOCATE( v_varrd_4dr )

      END SUBROUTINE Read_Met_Data_For_A1_Time

!--------------------------------------------------------

      SUBROUTINE Read_Grid_Data( State_Grid_Adj )

                 ! USE 
      USE State_Grid_Adj_Mod,     ONLY : GrdStateAdj
      USE Ncdf_Mod
      USE m_netcdf_io_checks,     ONLY : Ncdoes_Var_Exist
      USE m_netcdf_io_read,       ONLY : NcRd
      USE m_netcdf_io_close,      ONLY : Nccl
      USE m_netcdf_io_open,       ONLY : Ncop_Rd

      TYPE(GrdStateAdj), INTENT(IN)    :: State_Grid_Adj  ! Grid State object


      ! dimension 
      INTEGER                          :: NX, NY
      LOGICAL                          :: err_stop
      INTEGER                          :: stat
      INTEGER                          :: RC
      CHARACTER(LEN=255)               :: varname
      INTEGER                          :: fID

      !time 
      CHARACTER(LEN=255)               :: FILENAME 
      ! area
      REAL*8, ALLOCATABLE              :: area_varrd_2dr(:, :)
      INTEGER                          :: area_str2d(2)
      INTEGER                          :: area_cnt2d(2)


       ! initialization 
       NX = State_Grid_Adj%NX
       NY = State_Grid_Adj%NY

       FileName = 'FdCo2Out/GEOSChem.Area.2x25.nc4'
       CALL Ncop_Rd( fID, TRIM(FileName) )

       ! area
       ALLOCATE( area_varrd_2dr(NX, NY), STAT = RC )
       varname = 'AREA'
       area_str2d(1) = 1
       area_str2d(2) = 1

       area_cnt2d(1) = NX
       area_cnt2d(2) = NY

       CALL NcRd(area_varrd_2dr, fID, TRIM(varname), area_str2d, area_cnt2d)

       State_Grid_Adj%AREA_M2( :, : ) = area_varrd_2dr( :, : )

      ! close file
      CALL Nccl(fID)

      ! deallocate 
      IF ( ALLOCATED( area_varrd_2dr ) )    DEALLOCATE( area_varrd_2dr )

     END SUBROUTINE Read_Grid_Data


!--------------------------------------------------------
     SUBROUTINE WRITE_FILE_R8_4D( VarName, VarValue )


     CHARACTER(LEN=*),  INTENT(IN)   :: VarName
     REAL(fp), POINTER               :: VarValue(:, :, :, :)

  
     CALL Nc_Var_Def( fID         = Nc_File_Ope_Var%fID,                             &     
                      lonID       = Nc_File_Ope_Var%lonID,                           & 
                      latID       = Nc_File_Ope_Var%latID,                           &
                      levID       = Nc_File_Ope_Var%levID,                           &
                      TimeId      = Nc_File_Ope_Var%timeID,                          & 
                      VarName     = VarName,                                         &
                      VarLongName = 'this is '//VarName,                             &
                      VarUnit     = ' ',                                             &
                      DataType=8,                                                    & 
                      VarCt=Nc_File_Ope_Var%VarCT )

     CALL Nc_Var_Write( fID=Nc_File_Ope_Var%fID, VarName=VarName, Arr4D=VarValue )



     END SUBROUTINE WRITE_FILE_R8_4D

!--------------------------------------------------------
     SUBROUTINE WRITE_FILE_R8_3D( VarName, VarValue )


     CHARACTER(LEN=*),  INTENT(IN)   :: VarName
     REAL(fp), POINTER               :: VarValue(:, :, :)

     CALL Nc_Var_Def( fID         = Nc_File_Ope_Var%fID,                             &    
                      lonID       = Nc_File_Ope_Var%lonID,                           &
                      latID       = Nc_File_Ope_Var%latID,                           &
                      levID       = -1,                                              &
                      TimeId      = Nc_File_Ope_Var%timeID,                          &
                      VarName     = VarName,                                         &
                      VarLongName = 'this is '//VarName,                             &
                      VarUnit     = ' ',                                             &
                      DataType=8,                                                    & 
                      VarCt=Nc_File_Ope_Var%VarCT )

     CALL Nc_Var_Write( fID=Nc_File_Ope_Var%fID, VarName=VarName, Arr3D=VarValue )

     END SUBROUTINE WRITE_FILE_R8_3D

!--------------------------------------------------------
     SUBROUTINE WRITE_FILE_INT_4D( VarName, VarValue )


     CHARACTER(LEN=*),  INTENT(IN)   :: VarName
     INTEGER, POINTER               :: VarValue(:, :, :, :)


     CALL Nc_Var_Def( fID         = Nc_File_Ope_Var%fID,                             &
                      lonID       = Nc_File_Ope_Var%lonID,                           &
                      latID       = Nc_File_Ope_Var%latID,                           &
                      levID       = Nc_File_Ope_Var%levID,                           &
                      TimeId      = Nc_File_Ope_Var%timeID,                          &
                      VarName     = VarName,                                         &
                      VarLongName = 'this is '//VarName,                             &
                      VarUnit     = ' ',                                             &
                      DataType=8,                                                    &
                      VarCt=Nc_File_Ope_Var%VarCT )

     CALL Nc_Var_Write( fID=Nc_File_Ope_Var%fID, VarName=VarName, Arr4D=VarValue )



     END SUBROUTINE WRITE_FILE_INT_4D

!-------------------------------------------------------------
     SUBROUTINE WRITE_FILE_INT_3D( VarName, VarValue )


     CHARACTER(LEN=*),  INTENT(IN)   :: VarName
     INTEGER, POINTER                :: VarValue(:, :, :)


     CALL Nc_Var_Def( fID         = Nc_File_Ope_Var%fID,                             &
                      lonID       = Nc_File_Ope_Var%lonID,                           &
                      latID       = Nc_File_Ope_Var%latID,                           &
                      levID       = -1,                                              &
                      TimeId      = Nc_File_Ope_Var%timeID,                          &
                      VarName     = VarName,                                         &
                      VarLongName = 'this is '//VarName,                             &
                      VarUnit     = ' ',                                             &
                      DataType=8,                                                    &
                      VarCt=Nc_File_Ope_Var%VarCT )

     CALL Nc_Var_Write( fID=Nc_File_Ope_Var%fID, VarName=VarName, Arr3D=VarValue )

     END SUBROUTINE WRITE_FILE_INT_3D

!--------------------------------------------------------
     SUBROUTINE READ_FILE_R8_4D( fID, VarName, VarValue )

     USE m_netcdf_io_read,       ONLY : NcRd

     CHARACTER(LEN=*),  INTENT(IN)   :: VarName
     REAL(fp), POINTER               :: VarValue(:, :, :, :)

     LOGICAL                          :: err_stop
     INTEGER                          :: stat
     INTEGER                          :: RC
     INTEGER                          :: fID

      !time 
     CHARACTER(LEN=255)               :: FILENAME
      ! dimension
     INTEGER                          :: NX, NY, NZ, NTIME
     INTEGER                          :: str4d(4)
     INTEGER                          :: cnt4d(4)


     NX    = SIZE(VarValue, DIM=1)
     NY    = SIZE(VarValue, DIM=2)
     NZ    = SIZE(VarValue, DIM=3)
     NTIME = SIZE(VarValue, DIM=4)

     str4d(1) = 1
     str4d(2) = 1
     str4d(3) = 1
     str4d(4) = 1

     cnt4d(1) = NX
     cnt4d(2) = NY
     cnt4d(3) = NZ
     cnt4d(4) = NTIME

     CALL NcRd(varvalue, fID, TRIM(varname), str4d, cnt4d)

     END SUBROUTINE READ_FILE_R8_4D

!--------------------------------------------------------

     SUBROUTINE READ_FILE_R8_3D( fID, VarName, VarValue )

     USE m_netcdf_io_read,       ONLY : NcRd

     CHARACTER(LEN=*),  INTENT(IN)   :: VarName
     REAL(fp), POINTER               :: VarValue(:, :, :)

     LOGICAL                          :: err_stop
     INTEGER                          :: stat
     INTEGER                          :: RC
     INTEGER                          :: fID

      !time 
     CHARACTER(LEN=255)               :: FILENAME
      ! dimension
     INTEGER                          :: NX, NY,  NTIME
     INTEGER                          :: str3d(3)
     INTEGER                          :: cnt3d(3)


     NX    = SIZE(VarValue, DIM=1)
     NY    = SIZE(VarValue, DIM=2)
     NTIME = SIZE(VarValue, DIM=3)

     str3d(1) = 1
     str3d(2) = 1
     str3d(3) = 1

     cnt3d(1) = NX
     cnt3d(2) = NY
     cnt3d(3) = NTIME

     CALL NcRd(varvalue, fID, TRIM(varname), str3d, cnt3d)

     END SUBROUTINE READ_FILE_R8_3D

!-------------------------------------------------------
    SUBROUTINE READ_FILE_INT_4D( fID, VarName, VarValue )

     USE m_netcdf_io_read,       ONLY : NcRd

     CHARACTER(LEN=*),  INTENT(IN)   :: VarName
     INTEGER, POINTER               :: VarValue(:, :, :, :)

     LOGICAL                          :: err_stop
     INTEGER                          :: stat
     INTEGER                          :: RC
     INTEGER                          :: fID

      !time
     CHARACTER(LEN=255)               :: FILENAME
      ! dimension
     INTEGER                          :: NX, NY, NZ, NTIME
     INTEGER                          :: str4d(4)
     INTEGER                          :: cnt4d(4)


     NX    = SIZE(VarValue, DIM=1)
     NY    = SIZE(VarValue, DIM=2)
     NZ    = SIZE(VarValue, DIM=3)
     NTIME = SIZE(VarValue, DIM=4)

     str4d(1) = 1
     str4d(2) = 1
     str4d(3) = 1
     str4d(4) = 1

     cnt4d(1) = NX
     cnt4d(2) = NY
     cnt4d(3) = NZ
     cnt4d(4) = NTIME

     CALL NcRd(varvalue, fID, TRIM(varname), str4d, cnt4d)

     END SUBROUTINE READ_FILE_INT_4D

!-------------------------------------------------------
     SUBROUTINE READ_FILE_INT_3D( fID, VarName, VarValue )

      USE m_netcdf_io_read,       ONLY : NcRd

     CHARACTER(LEN=*),  INTENT(IN)   :: VarName
     INTEGER, POINTER                :: VarValue(:, :, :)

     LOGICAL                          :: err_stop
     INTEGER                          :: stat
     INTEGER                          :: RC
     INTEGER                          :: fID

      !time 
     CHARACTER(LEN=255)               :: FILENAME
      ! dimension
     INTEGER                          :: NX, NY, NTIME
     INTEGER                          :: str3d(3)
     INTEGER                          :: cnt3d(3)


     NX    = SIZE(VarValue, DIM=1)
     NY    = SIZE(VarValue, DIM=2)
     NTIME = SIZE(VarValue, DIM=3)

     str3d(1) = 1
     str3d(2) = 1
     str3d(3) = 1

     cnt3d(1) = NX
     cnt3d(2) = NY
     cnt3d(3) = NTIME

     CALL NcRd(varvalue, fID, TRIM(varname), str3d, cnt3d)

     END SUBROUTINE READ_FILE_INT_3D


!--------------------------------------------------------
     SUBROUTINE WRITE_FILE_INIT( NX, NY, NZ )

      !
      !read delp_dry from the file 
      USE TIME_MOD,           ONLY : GET_NYMD,   GET_HOUR 
      USE TIME_MOD,           ONLY : GET_MINUTE, GET_SECOND 
      USE TIME_MOD,           ONLY : GET_NHMS
      INTEGER, INTENT(IN)            :: NX, NY, NZ
      
      ! local variables
      INTEGER                        :: I, J, L
      INTEGER                        :: NYMD
      INTEGER                        :: HOUR
      INTEGER                        :: MINUTE
      INTEGER                        :: SECOND
      INTEGER                        :: NHMS
      CHARACTER(LEN=20)              :: CHA_NYMD
      CHARACTER(LEN=20)              :: CHA_HOUR
      CHARACTER(LEN=20)              :: CHA_MINUTE
      CHARACTER(LEN=20)              :: CHA_SECOND
      CHARACTER(LEN=255)             :: FILENAME
      ! huoxiao write variables to files
      NYMD             =          GET_NYMD()
      HOUR             =          GET_HOUR()
      MINUTE           =          GET_MINUTE()
      SECOND           =          GET_SECOND()
      NHMS             =          GET_NHMS()
 
      WRITE( CHA_NYMD,        *  ) NYMD
      WRITE( CHA_HOUR,  '(I2.2)' ) HOUR
      WRITE( CHA_MINUTE, '(I2.2)') MINUTE
      WRITE( CHA_SECOND, '(I2.2)') SECOND

   
      FILENAME = 'AdjParameters/GEOSChem.AdjParameters.'      &
                                  //TRIM(ADJUSTL(CHA_NYMD))   &
                                  //'_'                       &
                                  //TRIM(ADJUSTL(CHA_HOUR))   &
                                  //TRIM(ADJUSTL(CHA_MINUTE)) &
                                  //TRIM(ADJUSTL(CHA_SECOND)) &
                                  //'.nc4'
 
      WRITE(*,100) 'Adjoint Parameters', NYMD, NHMS 
      WRITE(*,110) FILENAME
      CALL Nc_Create( NcFile = FILENAME,               Title = 'adjoint parameters',          &
                      nLon   = NX,                     nLat  = NY,                            &  
                      nLev   = NZ,                     nTime = 1,                             &
                      fID    = Nc_File_Ope_Var%fID,    lonID = Nc_File_Ope_Var%lonID,         &
                      latID  = Nc_File_Ope_Var%latID,  levID = Nc_File_Ope_Var%levID,         &
                      timeID = Nc_File_Ope_Var%timeID, VarCt = Nc_File_Ope_Var%VarCt,         &
                      Create_NC4 = .true. )


100       FORMAT( '     - Creating file for ', a, '; reference = ',i8.8,1x,i6.6)
110       FORMAT( '        with filename = ', a                                )
 
      END SUBROUTINE WRITE_FILE_INIT


!--------------------------------------------------------

      SUBROUTINE WRITE_FILE_CLOSE

      USE m_netcdf_io_close,      ONLY : Nccl
 
      CALL Nccl( Nc_File_Ope_Var%fID )

      END SUBROUTINE WRITE_FILE_CLOSE

!      SUBROUTINE WRITE_FILE_R8_2D( FILE_PREFIX, VALUE, INDEX1, INDEX2 )

      !
      !read delp_dry from the file 
!      USE TIME_MOD,           ONLY : GET_NYMD,   GET_HOUR
!      USE TIME_MOD,           ONLY : GET_MINUTE, GET_SECOND

!      CHARACTER(LEN=*), INTENT(IN)   :: FILE_PREFIX
!      INTEGER,INTENT(IN)             :: INDEX1, INDEX2
!      REAL(fp), INTENT(INOUT)        :: VALUE(:, :)


      ! local variables
!      INTEGER                        :: I, J
!      INTEGER                        :: NYMD
!      INTEGER                        :: HOUR
!      INTEGER                        :: MINUTE
!      INTEGER                        :: SECOND
!      CHARACTER(LEN=20)              :: CHA_NYMD
!      CHARACTER(LEN=20)              :: CHA_HOUR
!      CHARACTER(LEN=20)              :: CHA_MINUTE
!      CHARACTER(LEN=20)              :: CHA_SECOND
!      CHARACTER(LEN=255)             :: FILENAME
      ! huoxiao write variables to files
!      NYMD             =          GET_NYMD()
!      HOUR             =          GET_HOUR()
!      MINUTE           =          GET_MINUTE()
!      SECOND           =          GET_SECOND()

!      WRITE( CHA_NYMD,        *  ) NYMD
!      WRITE( CHA_HOUR,  '(I2.2)' ) HOUR
!      WRITE( CHA_MINUTE, '(I2.2)') MINUTE
!      WRITE( CHA_SECOND, '(I2.2)') SECOND

!      FILENAME = TRIM(FILE_PREFIX)//TRIM(ADJUSTL(CHA_NYMD))   &
!                                  //TRIM(ADJUSTL(CHA_HOUR))   &
!                                  //TRIM(ADJUSTL(CHA_MINUTE)) &
!                                  //TRIM(ADJUSTL(CHA_SECOND))
!      OPEN(UNIT=1995, FILE=TRIM(FILENAME))

!      WRITE(1995,*) VALUE

!       CLOSE(1995)

!      END SUBROUTINE WRITE_FILE_R8_2D

!--------------------------------------------------------

!      SUBROUTINE WRITE_FILE_INT_2D( FILE_PREFIX, VALUE, INDEX1, INDEX2 )

      !
      !read delp_dry from the file 
!      USE TIME_MOD,           ONLY : GET_NYMD,   GET_HOUR 
!      USE TIME_MOD,           ONLY : GET_MINUTE, GET_SECOND 

!      CHARACTER(LEN=*), INTENT(IN)   :: FILE_PREFIX
!      INTEGER,INTENT(IN)             :: INDEX1, INDEX2
!      INTEGER, INTENT(INOUT)         :: VALUE(:, :)
!
      ! local variables
!      INTEGER                        :: I, J, L
!      INTEGER                        :: NYMD
!      INTEGER                        :: HOUR
!      INTEGER                        :: MINUTE
!      INTEGER                        :: SECOND
!      CHARACTER(LEN=20)              :: CHA_NYMD
!      CHARACTER(LEN=20)              :: CHA_HOUR
!      CHARACTER(LEN=20)              :: CHA_MINUTE
!      CHARACTER(LEN=20)              :: CHA_SECOND
!      CHARACTER(LEN=255)             :: FILENAME
      ! huoxiao write variables to files
!      NYMD             =          GET_NYMD()
!      HOUR             =          GET_HOUR()
!      MINUTE           =          GET_MINUTE()
!      SECOND           =          GET_SECOND()
 
!      WRITE( CHA_NYMD,        *  ) NYMD
!      WRITE( CHA_HOUR,  '(I2.2)' ) HOUR
!      WRITE( CHA_MINUTE, '(I2.2)') MINUTE
!      WRITE( CHA_SECOND, '(I2.2)') SECOND

!      FILENAME = TRIM(FILE_PREFIX)//TRIM(ADJUSTL(CHA_NYMD))   &
!                                  //TRIM(ADJUSTL(CHA_HOUR))   &
!       			          //TRIM(ADJUSTL(CHA_MINUTE)) &
!     			          //TRIM(ADJUSTL(CHA_SECOND))
!      OPEN(UNIT=1995, FILE=TRIM(FILENAME))
       
!       WRITE(1995,*) VALUE
 
!       CLOSE(1995)          
      
!      END SUBROUTINE WRITE_FILE_INT_2D

!--------------------------------------------------------
!      SUBROUTINE READ_FILE_R8_3D( FILE_PREFIX, VALUE, INDEX1, INDEX2, INDEX3 )

      !
      !read delp_dry from the file 
!      USE TIME_MOD,           ONLY : GET_NYMD,   GET_HOUR 
!      USE TIME_MOD,           ONLY : GET_MINUTE, GET_SECOND 

!      CHARACTER(LEN=*), INTENT(IN)   :: FILE_PREFIX
!      INTEGER,INTENT(IN)             :: INDEX1, INDEX2, INDEX3
!      REAL(fp),INTENT(INOUT)         :: VALUE(INDEX1, INDEX2, INDEX3)

      ! local variables
!      INTEGER                        :: I, J, L
!      INTEGER                        :: NYMD
!      INTEGER                        :: HOUR
!      INTEGER                        :: MINUTE
!      INTEGER                        :: SECOND
!      CHARACTER(LEN=20)              :: CHA_NYMD
!      CHARACTER(LEN=20)              :: !CHA_HOUR
!      CHARACTER(LEN=20)              :: CHA_MINUTE
!      CHARACTER(LEN=20)              :: CHA_SECOND
!      CHARACTER(LEN=255)             :: FILENAME
      ! huoxiao write variables to files
!      NYMD             =          GET_NYMD()
!      HOUR             =          GET_HOUR()
!      MINUTE           =          GET_MINUTE()
!      SECOND           =          GET_SECOND()
 
!      WRITE( CHA_NYMD,        *  ) NYMD
!      WRITE( CHA_HOUR,  '(I2.2)' ) HOUR
!      WRITE( CHA_MINUTE, '(I2.2)') MINUTE
!      WRITE( CHA_SECOND, '(I2.2)') SECOND

!      FILENAME = TRIM(FILE_PREFIX)//TRIM(ADJUSTL(CHA_NYMD))   &
!                                  //TRIM(ADJUSTL(CHA_HOUR))   &
!       			          //TRIM(ADJUSTL(CHA_MINUTE)) &
!     			          //TRIM(ADJUSTL(CHA_SECOND))
!      OPEN(UNIT=1995, FILE=TRIM(FILENAME))
       
!      READ(1995,*) VALUE

!      CLOSE(1995)          
!      END SUBROUTINE READ_FILE_R8_3D
!--------------------------------------------------------
!      SUBROUTINE READ_FILE_R8_2D( FILE_PREFIX, VALUE, INDEX1, INDEX2 )

      !
      !read delp_dry from the file 
!      USE TIME_MOD,           ONLY : GET_NYMD,   GET_HOUR
!      USE TIME_MOD,           ONLY : GET_MINUTE, GET_SECOND

!      CHARACTER(LEN=*), INTENT(IN)   :: FILE_PREFIX
!      INTEGER,INTENT(IN)             :: INDEX1, INDEX2
!      REAL(fp),INTENT(INOUT)         :: VALUE(INDEX1, INDEX2)

      ! local variables
!      INTEGER                        :: I, J
!      INTEGER                        :: NYMD
!      INTEGER                        :: HOUR
!      INTEGER                        :: MINUTE
!      INTEGER                        :: SECOND
!      CHARACTER(LEN=20)              :: CHA_NYMD
!      CHARACTER(LEN=20)              :: CHA_HOUR
!      CHARACTER(LEN=20)              :: CHA_MINUTE
!      CHARACTER(LEN=20)              :: CHA_SECOND
!      CHARACTER(LEN=255)             :: FILENAME
      ! huoxiao write variables to files
!      NYMD             =          GET_NYMD()
!      HOUR             =          GET_HOUR()
!      MINUTE           =          GET_MINUTE()
!      SECOND           =          GET_SECOND()

!      WRITE( CHA_NYMD,        *  ) NYMD
!      WRITE( CHA_HOUR,  '(I2.2)' ) HOUR
!      WRITE( CHA_MINUTE, '(I2.2)') MINUTE
!      WRITE( CHA_SECOND, '(I2.2)') SECOND

!      FILENAME = TRIM(FILE_PREFIX)//TRIM(ADJUSTL(CHA_NYMD))   &
!                                  //TRIM(ADJUSTL(CHA_HOUR))   &
!                                  //TRIM(ADJUSTL(CHA_MINUTE)) &
!                                  //TRIM(ADJUSTL(CHA_SECOND))
!      OPEN(UNIT=1995, FILE=TRIM(FILENAME))

!       READ(1995,*) VALUE

!       CLOSE(1995)
!      END SUBROUTINE READ_FILE_R8_2D

!--------------------------------------------------------

!      SUBROUTINE READ_FILE_INT_2D( FILE_PREFIX, VALUE, INDEX1, INDEX2)

      !
      !read delp_dry from the file 
!      USE TIME_MOD,           ONLY : GET_NYMD,   GET_HOUR 
!      USE TIME_MOD,           ONLY : GET_MINUTE, GET_SECOND 

!      CHARACTER(LEN=*), INTENT(IN)   :: FILE_PREFIX
!      INTEGER,INTENT(IN)             :: INDEX1, INDEX2
!      INTEGER, INTENT(INOUT)         :: VALUE(:, :)
!
      ! local variables
!      INTEGER                        :: I, J
!      INTEGER                        :: NYMD
!      INTEGER                        :: HOUR
!      INTEGER                        :: MINUTE
!      INTEGER                        :: SECOND
!      CHARACTER(LEN=20)              :: CHA_NYMD
!      CHARACTER(LEN=20)              :: CHA_HOUR
!      CHARACTER(LEN=20)              :: CHA_MINUTE
!      CHARACTER(LEN=20)              :: CHA_SECOND
!      CHARACTER(LEN=255)             :: FILENAME
      ! huoxiao write variables to files
!      NYMD             =          GET_NYMD()
!      HOUR             =          GET_HOUR()
!      MINUTE           =          GET_MINUTE()
!      SECOND           =          GET_SECOND()
 
!      WRITE( CHA_NYMD,        *  ) NYMD
!      WRITE( CHA_HOUR,  '(I2.2)' ) HOUR
!      WRITE( CHA_MINUTE, '(I2.2)') MINUTE
!      WRITE( CHA_SECOND, '(I2.2)') SECOND

!      FILENAME = TRIM(FILE_PREFIX)//TRIM(ADJUSTL(CHA_NYMD))   &
!                                  //TRIM(ADJUSTL(CHA_HOUR))   &
!       			          //TRIM(ADJUSTL(CHA_MINUTE)) &
!     			          //TRIM(ADJUSTL(CHA_SECOND))
!      OPEN(UNIT=1995, FILE=TRIM(FILENAME))
       
!      READ(1995,*) VALUE
 
!       CLOSE(1995)          
      
!      END SUBROUTINE READ_FILE_INT_2D


!------------------------------------------------------------------------------

      SUBROUTINE READ_SF_FILE (State_Grid_Adj, State_Chm_Adj)
!
!******************************************************************************
!  Subroutine MAKE_SF_FILE creates a binary file of STT_IC or EMS_ICS
!  (dkh, 9/17/04)
!
!  Module Variable as Input:
!  ============================================================================
!  (1 ) N_CALC    : Current iteration number
!  (2 ) ICS_SF    : Initial conditions scaling factors
!  (3 ) EMS_SF    : Emissions scaling factors
!
!  NOTES:
!  (1 ) Just like MAKE_ADJ_FILE except
!       - write to .ics. file
!  (2 ) Add support for ACTIVE_VARS == 'EMISSIONS' case (dkh, 11/27/04)
!  (3 ) Add support for ACTIVE_VARS == 'FDTEST' case. (dkh, 02/17/05)
!  (4 ) Change UNIT to unitless and change title to Scale factors (dkh, 03/06/05)
!  (5 ) Change output for ACTIVE_VARS == 'EMISSIONS' case.
!        Now use label IJ-EMS-$, and update gamap code accordingly.
!        First write the scaling factors, in consecutive species. Temporal
!         varations in the emissions, if any, will be in the L direction.
!        Next, write out the optimized emissions themselves.
!        Finally, write out the difference between orig and optimized emissions.
!        (dkh, 03/28/05)
!  (6 )  Use EMS_orig instead of ESO4_an_orig so that we can loop over N.
!  (7 )  Update to add support for writing NOx emissions. (dkh, 08/27/06)
!  (8 )  Only write the value of the scaling facotr in locations where the
!         actual emission is greater than zero.  Also include the current
!         scale emissions themselves in every *ics* file.  (dkh, 09/22/06)
!  (9 )  Add suppport for LOG_OPT
!  (10)  Standardize units for saving emissions. (dkh, 06/16/07)
!  (11)  Add option to print prior and posterior emissions totals. (dkh, 06/16/07)
!  (12)  Change names, replace CMN_ADJ. (dkh, ks, mak, cs  06/08/09)
!  (13)  Now support strat fluxes LADJ_STRAT (hml, dkh, 02/20/12, adj32_025)
!******************************************************************************
!
      ! References to F90 modules
      USE State_Chm_Adj_Mod,         ONLY :  EMS_SF
      USE State_Grid_Adj_Mod,        Only : GrdStateAdj
      USE State_Chm_Adj_Mod,        Only : ChmStateAdj
      ! huoxiao adjoint 
      USE ADJ_ARRAYS_MOD,      ONLY : N_CALC
    
      USE LOGICAL_ADJ_MOD,   ONLY :  LADJ_EMS, LICS
      USE LOGICAL_ADJ_MOD,   ONLY : LICS_INC


      TYPE(GrdStateAdj), INTENT(IN)       :: State_Grid_Adj
      TYPE(ChmStateAdj), INTENT(IN)       :: State_Chm_Adj

      ! Local Variables
      INTEGER              :: I,  J, M, N


      CHARACTER(LEN=255)   :: FILENAME

      CHARACTER(LEN=20)    :: TMP_FILE
      CHARACTER(LEN=20)    :: CHA_ITE_NUM

      !=================================================================
      ! READ_SF_FILE begins here!
      !=================================================================

      ! Hardwire output file for now
      WRITE(CHA_ITE_NUM, '(I2.2)') N_CALC 


      !=================================================================
      ! Open the adjoint file for output -- binary punch format
      !=================================================================

      IF ( LADJ_EMS ) THEN

      TMP_FILE = 'gctm.sf.'//TRIM( CHA_ITE_NUM )

      !=================================================================
      ! Open the adjoint file for output -- binary punch format
      !=================================================================

      FILENAME = TRIM( TMP_FILE )

      print*, 'scale factor '//TMP_FILE//' is reading'
      ! Open file for output

      OPEN( UNIT = 1995, FILE=FILENAME )

         !=================================================================
         ! Write each observed quantity to the ics file
         !=================================================================

      READ(1995, *) EMS_SF(:,:,1)
      CLOSE(1995)
      ENDIF

      IF ( LICS .AND. (.NOT. LICS_INC)) THEN

      TMP_FILE = 'gctm.gc.'//TRIM( CHA_ITE_NUM )

      !=================================================================
      ! Open the adjoint file for output -- binary punch format
      !=================================================================

      FILENAME = TRIM( TMP_FILE )

      print*, 'Species concentration '//TMP_FILE//' is reading'
      ! Open file for output

      OPEN( UNIT = 1995, FILE=FILENAME )

         !=================================================================
         ! Write each observed quantity to the ics file
         !=================================================================

      READ(1995, *) State_Chm_Adj%Species(:,:,:,1)
      CLOSE(1995)
      ENDIF

      IF ( LICS_INC ) THEN
      TMP_FILE = 'gctm.w.'//TRIM( CHA_ITE_NUM )
      FILENAME = TRIM( TMP_FILE )
      OPEN( UNIT = 1995, FILE=FILENAME )
      READ(1995, *) State_Chm_Adj%cv(:,:,:,1)
      CLOSE(1995)
      ENDIF


      ! Return to calling program
      END SUBROUTINE READ_SF_FILE

!------------------------------------------------------------------------------

      SUBROUTINE READ_CFN_FILE
!
!******************************************************************************
!  Subroutine MAKE_CFN_FILE creates a cfn.NN file which stores the current
!   iteration number and cost function value. (dkh, 02/13/06)
!
!  Module Variable as Input:
!  ============================================================================
!  (1 ) N_CALC    : Current iteration number
!
!  Module Variable as Output:
!  ============================================================================
!  (1 ) COST_FUNC : Current cost function value
!
!  NOTES:
!
!******************************************************************************
!
      ! Reference to f90 modules
      USE ADJ_ARRAYS_MOD,    ONLY : COST_FUNC
      ! huoxiao adjoint 
      USE ADJ_ARRAYS_MOD,      ONLY : N_CALC

      ! Local variables
      REAL*8             :: TEMP
      CHARACTER(LEN=80)  :: OUTPUT_CFN_FILE
      CHARACTER(LEN=20)  :: CHA_ITE_NUM
      !=================================================================
      ! MAKE_CFN_FILE begins here!
      !=================================================================


      WRITE(CHA_ITE_NUM, '(I2.2)') N_CALC 
      ! Make file name
      OUTPUT_CFN_FILE = 'cfn.'//TRIM( CHA_ITE_NUM )

      !=================================================================
      ! Open the cfn file for output
      !=================================================================

      PRINT*, 'cost function'//OUTPUT_CFN_FILE//'is reading'
      ! Open file for input
      OPEN(UNIT=1998,      FILE=TRIM(OUTPUT_CFN_FILE))

      ! Write iteration number and cost function
      READ( 1998, *) TEMP
      COST_FUNC = TEMP
     
      CLOSE(1998)

      END SUBROUTINE READ_CFN_FILE  

!-----------------------------------------------------------------------------
      SUBROUTINE READ_GDT_FILE(State_Chm_Adj, State_Grid_Adj )


!
!******************************************************************************
!  Subroutine MAKE_GDT_FILE creates a binary file of ADJ_xxx
!  (dkh, 9/17/04)
!
!  Module Variable as Input:
!  ============================================================================
!  (1 ) N_CALC     : Current iteration number
!  (2 ) ICS_SF_ADJ : Array of adjoint gradients to be written
!  (3 ) EMS_SF_ADJ : Array of adjoint gradients to be written
!
!  NOTES:
!  (1 ) Just like MAKE_OBS_FILE except
!       - write to .adj. file
!  (2 ) Changed name to MAKE_GDT_FILE.  Now the .adj. files are trajectories,
!       and the .gdt. files are final gradients  (dkh, 10/03/04)
!  (3 ) Added ACTIVE_VARS == 'EMISSIONS' case.  (dkh, 11/27/04)
!  (4 ) Added ACTIVE_VARS == 'FDTEST' case. (dkh, 02/17/05)
!  (5 ) Now use CATEGORY = 'IJ-GDE-$' for 'EMISSIONS' case. (dkh, 03/29/05)
!  (6 ) No longer pass COST_FUNC in the header; use cnf.* files. (dkh, 02/13/06)
!  (7 ) Rename everything, replace CMN_ADJ, move nonessential stuff
!       to diagnostic files  (dkh, ks, mak, cs  06/07/09)
!  (8 ) Add normalized gradients IJ-GDEN$ (dkh, 05/06/10)
!  (9 ) Now support strat fluxes LADJ_STRAT (hml, dkh, 02/20/12, adj32_025)
!******************************************************************************
!
      ! References to F90 modules
      USE LOGICAL_ADJ_MOD,     ONLY :  LADJ_EMS, LICS, LICS_INC
      USE State_Chm_Adj_Mod,         Only : ChmStateAdj
      USE State_Chm_Adj_Mod,         Only : EMS_SF_ADJ
      USE State_Grid_Adj_Mod,        Only : GrdStateAdj
      ! Added for reaction rate sensitivities (tww, 05/08/12)

      ! huoxiao adjoint 
      USE ADJ_ARRAYS_MOD,      ONLY : N_CALC

      TYPE(ChmStateAdj), INTENT(IN)       :: State_Chm_Adj
      TYPE(GrdStateAdj), INTENT(IN)       :: State_Grid_Adj

      ! Local Variables
      INTEGER              :: I, J , M, N
      CHARACTER(LEN=255)   :: FILENAME


      ! Added for reaction rate sensitivity (tww, 05/08/12)

      CHARACTER(LEN=20)    :: OUTPUT_GDT_FILE

      !huoxiao adjoint
      CHARACTER(LEN=20)    :: CHA_ITE_NUM
      !=================================================================
      ! MAKE_GDT_FILE begins here!
      !=================================================================



      ! Hardwire output file for now
      WRITE(CHA_ITE_NUM, '(I2.2)') N_CALC

      IF( LADJ_EMS ) THEN   

      OUTPUT_GDT_FILE = 'gctm.gdt.'//TRIM( CHA_ITE_NUM )
      !=================================================================
      ! Open the adjoint file for output -- binary punch format
      !=================================================================
      ! Add the OPTDATA_DIR prefix to the file name
      FILENAME = TRIM( OUTPUT_GDT_FILE )

      ! Open file for output

      OPEN( UNIT = 1997, FILE=FILENAME )


         !=================================================================
         ! Write each observed quantity to the observation file
         !=================================================================

            !Temporarily store quantities in the TRACER array
      READ(1997,*)  EMS_SF_ADJ(:,:,1)

      ENDIF


      IF( LICS .AND. (.NOT. LICS_INC) ) THEN
      OUTPUT_GDT_FILE = 'gctm.gdc.'//TRIM( CHA_ITE_NUM )
      !=================================================================
      ! Open the adjoint file for output -- binary punch format
      !=================================================================
      ! Add the OPTDATA_DIR prefix to the file name
      FILENAME = TRIM( OUTPUT_GDT_FILE )

      ! Open file for output

      OPEN( UNIT = 1997, FILE=FILENAME )


         !=================================================================
         ! Write each observed quantity to the observation file
         !=================================================================

            !Temporarily store quantities in the TRACER array
      READ(1997,*)  State_Chm_Adj%SpeciesAdj(:,:,:,1)
      print*, 'huoxiao_debug max 4dvar gradient', MAXVAL(State_Chm_Adj%SpeciesAdj(:,:,:,1))
      ! Close file
      CLOSE( 1997 )
      ENDIF

      IF ( LICS_INC ) THEN
      OUTPUT_GDT_FILE = 'gctm.gdw.'//TRIM( CHA_ITE_NUM )
      !=================================================================
      ! Open the adjoint file for output -- binary punch format
      !=================================================================
      ! Add the OPTDATA_DIR prefix to the file name
      FILENAME = TRIM( OUTPUT_GDT_FILE )

      ! Open file for output

      OPEN( UNIT = 1997, FILE=FILENAME )


         !=================================================================
         ! Write each observed quantity to the observation file
         !=================================================================

            !Temporarily store quantities in the TRACER array
      READ(1997,*)  State_Chm_Adj%SpeciesAdj(:,:,:,1)
      print*, 'huoxiao_debug max lics gradient', MAXVAL(State_Chm_Adj%SpeciesAdj(:,:,:,1))

      ! Close file
      CLOSE( 1997 )
      ENDIF

      ! Return to calling program
      END SUBROUTINE READ_GDT_FILE

!-----------------------------------------------------------------------------

      SUBROUTINE MAKE_SF_FILE (State_Grid_Adj, State_Chm_Adj)
!
!******************************************************************************
!  Subroutine MAKE_SF_FILE creates a binary file of STT_IC or EMS_ICS
!  (dkh, 9/17/04)
!
!  Module Variable as Input:
!  ============================================================================
!  (1 ) N_CALC    : Current iteration number
!  (2 ) ICS_SF    : Initial conditions scaling factors
!  (3 ) EMS_SF    : Emissions scaling factors
!
!  NOTES:
!  (1 ) Just like MAKE_ADJ_FILE except
!       - write to .ics. file
!  (2 ) Add support for ACTIVE_VARS == 'EMISSIONS' case (dkh, 11/27/04)
!  (3 ) Add support for ACTIVE_VARS == 'FDTEST' case. (dkh, 02/17/05)
!  (4 ) Change UNIT to unitless and change title to Scale factors (dkh, 03/06/05)
!  (5 ) Change output for ACTIVE_VARS == 'EMISSIONS' case.
!        Now use label IJ-EMS-$, and update gamap code accordingly.
!        First write the scaling factors, in consecutive species. Temporal
!         varations in the emissions, if any, will be in the L direction.
!        Next, write out the optimized emissions themselves.
!        Finally, write out the difference between orig and optimized emissions.
!        (dkh, 03/28/05)
!  (6 )  Use EMS_orig instead of ESO4_an_orig so that we can loop over N.
!  (7 )  Update to add support for writing NOx emissions. (dkh, 08/27/06)
!  (8 )  Only write the value of the scaling facotr in locations where the
!         actual emission is greater than zero.  Also include the current
!         scale emissions themselves in every *ics* file.  (dkh, 09/22/06)
!  (9 )  Add suppport for LOG_OPT
!  (10)  Standardize units for saving emissions. (dkh, 06/16/07)
!  (11)  Add option to print prior and posterior emissions totals. (dkh, 06/16/07)
!  (12)  Change names, replace CMN_ADJ. (dkh, ks, mak, cs  06/08/09)
!  (13)  Now support strat fluxes LADJ_STRAT (hml, dkh, 02/20/12, adj32_025)
!******************************************************************************
!
      ! References to F90 modules
      USE State_Chm_Adj_Mod,         ONLY :  EMS_SF
      USE State_Grid_Adj_Mod,        Only : GrdStateAdj
      USE State_Chm_Adj_Mod,        Only : ChmStateAdj
      ! huoxiao adjoint 
      USE ADJ_ARRAYS_MOD,      ONLY : N_CALC

      USE LOGICAL_ADJ_MOD,   ONLY :  LADJ_EMS, LICS, LICS_INC


      TYPE(GrdStateAdj), INTENT(IN)       :: State_Grid_Adj

      TYPE(ChmStateAdj), INTENT(IN)       :: State_Chm_Adj

      ! Local Variables
      INTEGER              :: I,  J, M, N


      CHARACTER(LEN=255)   :: FILENAME

      CHARACTER(LEN=20)    :: OUTPUT_SF_FILE
      CHARACTER(LEN=20)    :: CHA_ITE_NUM

      REAL                 :: TCVV_CO2
      !=================================================================
      ! READ_SF_FILE begins here!
      !=================================================================
      TCVV_CO2 = 44.0/28.97
      ! Hardwire output file for now
      WRITE(CHA_ITE_NUM, '(I2.2)') N_CALC


      !=================================================================
      ! Open the adjoint file for output -- binary punch format
      !=================================================================

      IF ( LADJ_EMS ) THEN
      OUTPUT_SF_FILE = 'gctm.sf.'//TRIM( CHA_ITE_NUM )

      !=================================================================
      ! Open the adjoint file for output -- binary punch format
      !=================================================================

      FILENAME = TRIM( OUTPUT_SF_FILE )
      print*, 'scale factor '//OUTPUT_SF_FILE//' is output'
      ! Open file for output

      OPEN( UNIT = 1995, FILE=FILENAME )


         !=================================================================
         ! Write each observed quantity to the ics file
         !=================================================================
      WRITE(1995, *) EMS_SF(:,:,1)
      CLOSE(1995)
      ENDIF


      IF ( LICS .AND. (.NOT. LICS_INC)) THEN
      OUTPUT_SF_FILE = 'gctm.gc.'//TRIM( CHA_ITE_NUM )

      !=================================================================
      ! Open the adjoint file for output -- binary punch format
      !=================================================================

      FILENAME = TRIM( OUTPUT_SF_FILE )
      print*, 'Species Concentration '//OUTPUT_SF_FILE//' is output'
      ! Open file for output

      OPEN( UNIT = 1995, FILE=FILENAME )


         !=================================================================
         ! Write each observed quantity to the ics file
         !=================================================================
      WRITE(1995, *) State_Chm_Adj%Species(:,:,:,1)*TCVV_CO2
      CLOSE(1995)

      IF(N_CALC==1) THEN
      OUTPUT_SF_FILE = 'BEC_wn13/gctm.gc.01'
      PRINT*, 'WRITE_BEC', OUTPUT_SF_FILE
      FILENAME = TRIM( OUTPUT_SF_FILE )
      OPEN( UNIT = 1995, FILE=FILENAME )
      WRITE(1995, *) State_Chm_Adj%Species(:,:,:,1)*TCVV_CO2
      CLOSE(1995)
      OUTPUT_SF_FILE = 'BEC_wn13/gctm.gc.cur'
      FILENAME = TRIM( OUTPUT_SF_FILE )
      OPEN( UNIT = 1995, FILE=FILENAME )
      WRITE(1995, *) State_Chm_Adj%Species(:,:,:,1)*TCVV_CO2
      CLOSE(1995)

      OPEN( UNIT = 1995, FILE='BEC_wn13/SF_SUCCESS' )
      WRITE(1995, *) 1
      CLOSE(1995)


      ELSE
      OUTPUT_SF_FILE = 'BEC_wn13/gctm.gc.cur'
      FILENAME = TRIM( OUTPUT_SF_FILE )
      OPEN( UNIT = 1995, FILE=FILENAME )
      WRITE(1995, *) State_Chm_Adj%Species(:,:,:,1)*TCVV_CO2
      CLOSE(1995)

      OPEN( UNIT = 1995, FILE='BEC_wn13/SF_SUCCESS' )
      WRITE(1995, *) 1
      CLOSE(1995)

      ENDIF

      
      ENDIF

      IF ( LICS_INC ) THEN
      OUTPUT_SF_FILE = 'gctm.w.'//TRIM( CHA_ITE_NUM )
      FILENAME = TRIM( OUTPUT_SF_FILE )
      OPEN( UNIT = 1995, FILE=FILENAME )
      WRITE(1995, *) State_Chm_Adj%CV(:,:,:,1)
      CLOSE(1995)
      ENDIF
      ! Return to calling program
      END SUBROUTINE MAKE_SF_FILE
!-----------------------------------------------------------------------------

      SUBROUTINE MAKE_CFN_FILE
!
!******************************************************************************
!  Subroutine MAKE_CFN_FILE creates a cfn.NN file which stores the current
!   iteration number and cost function value. (dkh, 02/13/06)
!
!  Module Variable as Input:
!  ============================================================================
!  (1 ) N_CALC    : Current iteration number
!
!  Module Variable as Output:
!  ============================================================================
!  (1 ) COST_FUNC : Current cost function value
!
!  NOTES:
!
!
      ! Reference to f90 modules
      USE ADJ_ARRAYS_MOD,    ONLY : COST_FUNC
      ! huoxiao adjoint 
      USE ADJ_ARRAYS_MOD,      ONLY : N_CALC
      USE ADJ_ARRAYS_MOD,    ONLY : TOT_ITER
      USE LOGICAL_ADJ_MOD,   ONLY : LICS_INC
      ! Local variables
      CHARACTER(LEN=80)  :: OUTPUT_CFN_FILE
      CHARACTER(LEN=20)  :: CHA_ITE_NUM
      CHARACTER(LEN=20)  :: CHA_TOT_ITER
      !=================================================================
      ! MAKE_CFN_FILE begins here!
      !=================================================================


      WRITE(CHA_ITE_NUM, '(I2.2)') N_CALC
      WRITE(CHA_TOT_ITER, '(I4.4)') TOT_ITER
      ! Make file name
      OUTPUT_CFN_FILE = 'cfn.'//TRIM( CHA_ITE_NUM )

      !=================================================================
      ! Open the cfn file for output
      !=================================================================

      PRINT*, 'cost function'//OUTPUT_CFN_FILE//'is output'

      ! Open file for input
      OPEN(UNIT=1998,      FILE=TRIM(OUTPUT_CFN_FILE))

      ! Write iteration number and cost function
      WRITE( 1998, *) COST_FUNC

      CLOSE(1998)

      IF ( LICS_INC ) THEN
      !write cost function to bac
      OUTPUT_CFN_FILE = 'cost_function/cfn.'//TRIM(CHA_TOT_ITER)//'.'//TRIM( CHA_ITE_NUM )
      OPEN(UNIT=1998,      FILE=TRIM(OUTPUT_CFN_FILE))
      WRITE( 1998, *) COST_FUNC
      CLOSE(1998)
      ENDIF
      END SUBROUTINE MAKE_CFN_FILE
!-----------------------------------------------------

      SUBROUTINE MAKE_CFNJo_FILE
!
!******************************************************************************
!  Subroutine MAKE_CFN_FILE creates a cfn.NN file which stores the current
!   iteration number and cost function value. (dkh, 02/13/06)
!
!  Module Variable as Input:
!  ============================================================================
!  (1 ) N_CALC    : Current iteration number
!
!  Module Variable as Output:
!  ============================================================================
!  (1 ) COST_FUNC : Current cost function value
!
!  NOTES:
!
!******************************************************************************
!
      ! Reference to f90 modules
      USE ADJ_ARRAYS_MOD,    ONLY : COST_FUNC_Jo
      ! huoxiao adjoint
      USE ADJ_ARRAYS_MOD,    ONLY : N_CALC
      USE ADJ_ARRAYS_MOD,    ONLY : TOT_ITER
      USE LOGICAL_ADJ_MOD,   ONLY : LICS_INC
      ! Local variables
      CHARACTER(LEN=80)  :: OUTPUT_CFN_FILE
      CHARACTER(LEN=20)  :: CHA_ITE_NUM
      CHARACTER(LEN=20)  :: CHA_TOT_ITER
      !=================================================================
      ! MAKE_CFN_FILE begins here!
      !=================================================================


      WRITE(CHA_ITE_NUM, '(I2.2)') N_CALC
      WRITE(CHA_TOT_ITER, '(I4.4)') TOT_ITER
      ! Make file name
      OUTPUT_CFN_FILE = 'cfnjo.'//TRIM( CHA_ITE_NUM )

      !=================================================================
      ! Open the cfn file for output
      !=================================================================

      PRINT*, 'cost function'//OUTPUT_CFN_FILE//'is output'

      ! Open file for input
      OPEN(UNIT=1998,      FILE=TRIM(OUTPUT_CFN_FILE))

      ! Write iteration number and cost function
      WRITE( 1998, *) COST_FUNC_Jo

      CLOSE(1998)

      !write cost function to bac
      IF ( LICS_INC ) THEN
      OUTPUT_CFN_FILE = 'cost_function/cfnjo.'//TRIM(CHA_TOT_ITER)//'.'//TRIM( CHA_ITE_NUM )
      OPEN(UNIT=1998,      FILE=TRIM(OUTPUT_CFN_FILE))
      WRITE( 1998, *) COST_FUNC_Jo
      CLOSE(1998)
      ENDIF
      END SUBROUTINE MAKE_CFNJo_FILE
!---------------------------------------------------
!-----------------------------------------------------------------------------
      SUBROUTINE MAKE_GDT_FILE(State_Chm_Adj, State_Grid_Adj )


!
!******************************************************************************
!  Subroutine MAKE_GDT_FILE creates a binary file of ADJ_xxx
!  (dkh, 9/17/04)
!
!  Module Variable as Input:
!  ============================================================================
!  (1 ) N_CALC     : Current iteration number
!  (2 ) ICS_SF_ADJ : Array of adjoint gradients to be written
!  (3 ) EMS_SF_ADJ : Array of adjoint gradients to be written
!
!  NOTES:
!  (1 ) Just like MAKE_OBS_FILE except
!       - write to .adj. file
!  (2 ) Changed name to MAKE_GDT_FILE.  Now the .adj. files are trajectories,
!       and the .gdt. files are final gradients  (dkh, 10/03/04)
!  (3 ) Added ACTIVE_VARS == 'EMISSIONS' case.  (dkh, 11/27/04)
!  (4 ) Added ACTIVE_VARS == 'FDTEST' case. (dkh, 02/17/05)
!  (5 ) Now use CATEGORY = 'IJ-GDE-$' for 'EMISSIONS' case. (dkh, 03/29/05)
!  (6 ) No longer pass COST_FUNC in the header; use cnf.* files. (dkh, 02/13/06)
!  (7 ) Rename everything, replace CMN_ADJ, move nonessential stuff
!       to diagnostic files  (dkh, ks, mak, cs  06/07/09)
!  (8 ) Add normalized gradients IJ-GDEN$ (dkh, 05/06/10)
!  (9 ) Now support strat fluxes LADJ_STRAT (hml, dkh, 02/20/12, adj32_025)
!******************************************************************************
!
      ! References to F90 modules
      USE LOGICAL_ADJ_MOD,           OnlY :  LADJ_EMS, LICS
      USE State_Chm_Adj_Mod,         Only : ChmStateAdj
      Use State_Chm_Adj_Mod,         Only : EMS_SF_ADJ
      USE State_Grid_Adj_Mod,        Only : GrdStateAdj
      ! Added for reaction rate sensitivities (tww, 05/08/12)
      USE LOGICAL_ADJ_MOD,   ONLY :  LADJ_EMS, LICS, LICS_INC
      ! huoxiao adjoint 
      USE ADJ_ARRAYS_MOD,            Only : N_CALC

      TYPE(ChmStateAdj), INTENT(IN)       :: State_Chm_Adj
      TYPE(GrdStateAdj), INTENT(IN)       :: State_Grid_Adj

      ! Local Variables
      INTEGER              :: I, J , L, M, N
      CHARACTER(LEN=255)   :: FILENAME


      ! Added for reaction rate sensitivity (tww, 05/08/12)

     CHARACTER(LEN=20)    :: OUTPUT_GDT_FILE

      !huoxiao adjoint
      CHARACTER(LEN=20)    :: CHA_ITE_NUM
      !=================================================================
      ! MAKE_GDT_FILE begins here!
      !=================================================================



      ! Hardwire output file for now
      WRITE(CHA_ITE_NUM, '(I2.2)')N_CALC

      IF ( LADJ_EMS ) THEN
      OUTPUT_GDT_FILE = 'gctm.gdt.'//TRIM( CHA_ITE_NUM )
      !=================================================================
      ! Open the adjoint file for output -- binary punch format
      !=================================================================
      ! Add the OPTDATA_DIR prefix to the file name
      FILENAME = TRIM( OUTPUT_GDT_FILE )

      ! Open file for output

      OPEN( UNIT = 1997, FILE=FILENAME )


         !=================================================================
         ! Write each observed quantity to the observation file
         !=================================================================

            !Temporarily store quantities in the TRACER array
      WRITE(1997,*)  EMS_SF_ADJ(:,:,1)
      CLOSE(1997)

      ENDIF
     
      IF ( LICS .AND. (.NOT. LICS_INC) ) THEN 
      ! Close file IF ( LADJ_EMS ) THEN
      OUTPUT_GDT_FILE = 'gctm.gdc.'//TRIM( CHA_ITE_NUM )
      !=================================================================
      ! Open the adjoint file for output -- binary punch format
      !=================================================================
      ! Add the OPTDATA_DIR prefix to the file name
      FILENAME = TRIM( OUTPUT_GDT_FILE )

      ! Open file for output

      OPEN( UNIT = 1997, FILE=FILENAME )


         !=================================================================
         ! Write each observed quantity to the observation file
         !=================================================================

      WRITE(1997,*)  State_Chm_Adj%SpeciesAdj(:,:,:,1)

      CLOSE(1995)
      ENDIF


      IF ( LICS_INC ) THEN
      ! Close file IF ( LADJ_EMS ) THEN
      OUTPUT_GDT_FILE = 'gctm.gdw.'//TRIM( CHA_ITE_NUM )
      !=================================================================
      ! Open the adjoint file for output -- binary punch format
      !=================================================================
      ! Add the OPTDATA_DIR prefix to the file name
      FILENAME = TRIM( OUTPUT_GDT_FILE )

      ! Open file for output

      OPEN( UNIT = 1997, FILE=FILENAME )


         !=================================================================
         ! Write each observed quantity to the observation file
         !=================================================================

            !Temporarily store quantities in the TRACER array
!            DO L = 1, State_Grid_Adj%NZ
!            DO J = 1, State_Grid_Adj%NY
!            DO I = 1, State_Grid_Adj%NX
               WRITE(1997,*)  State_Chm_Adj%SpeciesAdj(:,:,:,1)
!            ENDDO
!            ENDDO
!            ENDDO

       CLOSE(1997)
      ENDIF


      ! Return to calling program
      END SUBROUTINE MAKE_GDT_FILE

      SUBROUTINE GET_GRADNT_FROM_ADJ( State_Chm_Adj, State_Grid_Adj )
!
!*****************************************************************************
!  Subroutine GET_GRADNT_FROM_ADJ compiles the gradient vector from the array
!   of adjoint values.  (dkh, 9/16/04)
!
!  NOTES:
!  (1 ) Added ACTIVE_VARS == 'EMISSIONS' case.  (dkh, 11/27/04)
!  (2 ) Added ACTIVE_VARS == 'FDTEST' case. (dkh, 02/17/05)
!  (3 ) Don't zero the NIT, NH4 and NO3 gradnts (dkh, 03/03/05)
!  (4 ) Now support strat fluxes LADJ_STRAT (hml, dkh, 02/20/12, adj32_025)
!
!*****************************************************************************
!
      ! Reference to f90 modules
      USE INVERSE_MOD,               ONLY : GRADNT
      USE LOGICAL_ADJ_MOD,           ONLY : LADJ_EMS, LICS, LICS_INC
      USE State_Chm_Adj_Mod,         Only : ChmStateAdj
      USE State_Chm_Adj_Mod,         Only : EMS_SF_ADJ
      USE State_Grid_Adj_Mod,        Only : GrdStateAdj
      USE ADJ_ARRAYS_MOD,            ONLY : OBS_EXIST
      ! Added for reaction rate sensitivities (tww, 05/08/12)


      TYPE(ChmStateAdj), INTENT(IN)       :: State_Chm_Adj
      TYPE(GrdStateAdj), INTENT(IN)       :: State_Grid_Adj
 


      ! Local Variables
      INTEGER :: I, J, L, M
      INTEGER :: I_DUM

         IF ( LADJ_EMS ) THEN
!$OMP PARALLEL DO &
!$OMP DEFAULT( SHARED ) &
!$OMP PRIVATE( J, I, I_DUM ) 
            DO J = 1, State_Grid_Adj%NY
            DO I = 1, State_Grid_Adj%NX

               I_DUM    = I + State_Grid_Adj%NX * ( J - 1)    
                     
               GRADNT(I_DUM) = EMS_SF_ADJ(I,J,1)

            ENDDO
            ENDDO
!$OMP END PARALLEL DO
         ENDIF



         IF ( LICS .OR. LICS_INC) THEN
!$OMP PARALLEL DO &
!$OMP DEFAULT( SHARED ) &
!$OMP PRIVATE( L, J, I, I_DUM ) 
            DO L = 1, State_Grid_Adj%NZ
            DO J = 1, State_Grid_Adj%NY
            DO I = 1, State_Grid_Adj%NX

               I_DUM    = I + State_Grid_Adj%NX * ( J - 1)&
                            + State_Grid_Adj%NX*State_Grid_Adj%NY*(L-1)

               GRADNT(I_DUM) = State_Chm_Adj%SpeciesAdj(I,J,L,1)
!               IF ( OBS_EXIST(I, J) == 1) THEN
!               PRINT*, 'I, J, GRADIENT ', I, J, GRADNT(I_DUM)
!               ENDIF
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO
!               print*, 'MIN GRAD', MINVAL(GRADNT)
         ENDIF
     END SUBROUTINE GET_GRADNT_FROM_ADJ
!--------------------------------------------------------------------------------
       SUBROUTINE WRITE_SF_TO_RSTFILE( State_Chm_Adj, State_Grid_Adj )

       USE NCDF_MOD
       USE m_netcdf_io_checks,   ONLY : Ncdoes_Var_Exist
       USE m_netcdf_io_read,     ONLY : NcRd
       USE State_Chm_Adj_Mod,    ONLY : ChmStateAdj
       USE State_Grid_Adj_Mod,   ONLY : GrdStateAdj
       USE m_netcdf_io_close,    ONLY : Nccl
       USE Logical_Adj_Mod,      ONLY : LICS_INC

       TYPE(GrdStateAdj), INTENT(IN)    :: State_Grid_Adj  ! Grid State object
       TYPE(ChmStateAdj), INTENT(INOUT) :: State_Chm_Adj   ! Chemistry State object
       CHARACTER(LEN=255)              FILENAME
      ! dimension
       INTEGER                          :: NX, NY, NZ, nSpecies, nTimes, RC
       CHARACTER(LEN=255)               :: TIMESTR, varname
      !related nc file read
       INTEGER                          :: ori_time_str1d(1)
       INTEGER                          :: ori_time_cnt1d(1)


       INTEGER                          :: ori_lev_str1d(1)
       INTEGER                          :: ori_lev_cnt1d(1)


       INTEGER                          :: ori_lon_str1d(1)
       INTEGER                          :: ori_lon_cnt1d(1)

       INTEGER                          :: ori_lat_str1d(1)
       INTEGER                          :: ori_lat_cnt1d(1)

       INTEGER                          :: fID, lonID, levID, latID, timeID, VarCT
       REAL(fp), POINTER                :: ori_time(:), ori_lat(:), ori_lon(:), ori_lev(:)
       REAL(fp), POINTER                :: Variable(:, :, :, :)
       REAL(fp)                         :: nt_time(1), nt_lat(State_Grid_Adj%NY), &
                                           nt_lon(State_Grid_Adj%NX), nt_lev(State_Grid_Adj%NZ)

       IF ( LICS_INC ) THEN
          Variable => State_Chm_Adj%Species0
          print*, 'huoxiao_debug_max_variable, species0 ',MAXVAL(Variable), MINVAL(Variable), &
          MAXVAL(State_Chm_Adj%Species0), MINVAL(State_Chm_Adj%Species0)
       ELSE
          Variable => State_Chm_Adj%Species
       ENDIF
       NX = State_Grid_Adj%NX
       NY = State_Grid_Adj%NY
       NZ = State_Grid_Adj%NZ
       ALLOCATE(ori_time(1), STAT = RC )
       ALLOCATE(ori_lat(NY), STAT = RC )

       ALLOCATE(ori_lon(NX), STAT = RC )

       ALLOCATE(ori_lev(NZ), STAT = RC )

  
       TIMESTR = GET_MODEL_BEGIN_SIMULATION_TIME()
       FILENAME = 'GEOSRestart/initial_GEOSChem_rst.2x25_CO2.nc'
       print*, 'Write Species in ', FILENAME

       CALL Nc_Open(TRIM(FILENAME), fID)

       print*, 'reading time'
       ori_time_str1d(1) = 1
       ori_time_cnt1d(1) = 1
       VarName = 'time'
       CALL NcRd(nt_time, fID, TRIM(VarName), ori_time_str1d, ori_time_cnt1d)
       ori_time = nt_time

       print*, 'reading lat'
       ori_lat_str1d(1) = 1
       ori_lat_cnt1d(1) = NY
       VarName = 'lat'
       CALL NcRd(nt_lat, fID, TRIM(VarName), ori_lat_str1d, ori_lat_cnt1d)
       ori_lat = nt_lat

       ori_lon_str1d(1) = 1
       ori_lon_cnt1d(1) = NX
       VarName = 'lon'
       CALL NcRd(nt_lon, fID, TRIM(VarName), ori_lon_str1d, ori_lon_cnt1d)
       ori_lon = nt_lon

       ori_lev_str1d(1) = 1
       ori_lev_cnt1d(1) = NZ
       VarName = 'lev'


       CALL NcRd(ori_lev, fID, TRIM(VarName), ori_lev_str1d, ori_lev_cnt1d)
       ori_lev = nt_lev
       CALL Nccl(fID)
       VarName = 'SpeciesRst_CO2'

       CALL Nc_Create( NcFile = FILENAME,               Title = 'adjoint parameters',          &
                       nLon   = NX,                     nLat  = NY,                            &
                       nLev   = NZ,                     nTime = 1,                             &
                       fID    = fID,    lonID = lonID,         &
                       latID  = latID,  levID = levID,         &
                       timeID = timeID, VarCt = VarCt,         &
		       Create_NC4 = .true. )

       CALL Nc_Var_Def( fID         = fID,                             &
                        lonID       = lonID,                           &
                        latID       = latID,                           &
                        levID       = levID,                           &
                        TimeId      = timeID,                          &
                        VarName     = VarName,                                         &
                        VarLongName = 'CO2 tracer',                             &
                        VarUnit     = 'mol/mol',                                             &
                        DataType=8,                                                    &
                        VarCt=VarCT )

       CALL Nc_Var_Write( fID=fID, VarName=VarName, Arr4D=Variable )




       VarName = 'time'


       CALL  Nc_Var_Def( fID         = fID,                             &
                        lonID       = -1,                           &
                        latID       = -1,                           &
                        levID       = -1,                           &
                        TimeId      = timeID,                          &
                        VarName     = VarName,                                         &
                        VarLongName = 'Time',                             &
                        VarUnit     = 'hours since 1985-1-1 00:00:0.0',                                             &
                        DataType=8,                                                    &
                        VarCt=VarCT )



       VarName = 'lev'


       CALL Nc_Var_Def( fID         = fID,                             &
                        lonID       = -1,                           &
                        latID       = -1,                           &
                        levID       = levID,                           &
                        TimeId      = -1,                          &
                        VarName     = VarName,                                         &
                        VarLongName = 'Eta Centers',                             &
                        VarUnit     = 'sigma_level',                                             &
                        DataType=8,                                                    &
                        VarCt=VarCT )

       CALL Nc_Var_Write( fID=fID, VarName=VarName, Arr1D=ori_lev )

    

       VarName = 'lat'


       CALL Nc_Var_Def( fID         = fID,                             &
                        lonID       = -1,                           &
                        latID       = latID,                           &
                        levID       = -1,                           &
                        TimeId      = -1,                          &
                        VarName     = VarName,                                         &
                        VarLongName = 'Latitude',                             &
                        VarUnit     = 'degrees_north',                                             &
                        DataType=8,                                                    &
                        VarCt=VarCT )
       CALL Nc_Var_Write( fID=fID, VarName=VarName, Arr1D=ori_lat )

       VarName = 'lon'


      CALL Nc_Var_Def( fID         = fID,                             &
                        lonID       = lonID,                           &
                        latID       = -1,                           &
                        levID       = -1,                           &
                        TimeId      = -1,                          &
                        VarName     = VarName,                                         &
                        VarLongName = 'Longitude',                             &
                        VarUnit     = 'degrees_east',                                             &
                        DataType=8,                                                    &
                        VarCt=VarCT )

       CALL Nc_Var_Write( fID=fID, VarName=VarName, Arr1D=ori_lon )

       CALL Nccl(fID)

!       IF ( ASSOCIATED( ori_time ) )        DEALLOCATE( ori_time )
!       IF ( ASSOCIATED( ori_lev ) )        DEALLOCATE( ori_lev )
!       IF ( ASSOCIATED( ori_lat ) )        DEALLOCATE( ori_lat )
!       IF ( ASSOCIATED( ori_lon ) )        DEALLOCATE( ori_lon )
       END SUBROUTINE WRITE_SF_TO_RSTFILE


       SUBROUTINE WRITE_SF_TO_SPCFILE( State_Chm_Adj, State_Grid_Adj )

       USE NCDF_MOD
       USE m_netcdf_io_checks,   ONLY : Ncdoes_Var_Exist
       USE m_netcdf_io_read,     ONLY : NcRd
       USE State_Chm_Adj_Mod,    ONLY : ChmStateAdj
       USE State_Grid_Adj_Mod,   ONLY : GrdStateAdj
       USE m_netcdf_io_close,      ONLY : Nccl


       TYPE(GrdStateAdj), INTENT(IN)    :: State_Grid_Adj  ! Grid State object
       TYPE(ChmStateAdj), INTENT(INOUT) :: State_Chm_Adj   ! Chemistry State object
       CHARACTER(LEN=255)               :: FILENAME, TIMESTR, VarName
      ! dimension
       INTEGER                          :: NX, NY, NZ, nSpecies, nTimes, RC

       REAL(fp), POINTER                :: co2rst_varr4d(:, :, :, :)
       INTEGER                          :: co2rst_str4d(4)
       INTEGER                          :: co2rst_cnt4d(4)

       INTEGER                          :: fID, lonID, levID, latID, timeID, VarCT


  
       NX = State_Grid_Adj%NX
       NY = State_Grid_Adj%NY
       NZ = State_Grid_Adj%NZ

       ALLOCATE(co2rst_varr4d(NX, NY, NZ, 1), STAT = RC)

       TIMESTR = GET_MODEL_BEGIN_SIMULATION_TIME()
       FILENAME = 'GEOSRestart/initial_GEOSChem_rst.2x25_CO2.nc'
       print*, 'Write Species in ', FILENAME

       CALL Nc_Open(TRIM(FILENAME), fID)

       co2rst_str4d(1) = 1
       co2rst_str4d(2) = 1
       co2rst_str4d(3) = 1
       co2rst_str4d(4) = 1

       co2rst_cnt4d(1) = NX
       co2rst_cnt4d(2) = NY
       co2rst_cnt4d(3) = NZ
       co2rst_cnt4d(4) = 1

       VarName = 'SpeciesRst_CO2'
       CALL NcRd(co2rst_varr4d, fID, TRIM(VarName), co2rst_str4d, co2rst_cnt4d)
       CALL Nccl(fID)


       FILENAME = 'FdCo2Out/GEOSChem.SpeciesConc.'//TRIM(TIMESTR(1:13))//'z.nc4'

       CALL Nc_Create( NcFile = FILENAME,               Title = 'adjoint parameters',          &
                       nLon   = NX,                     nLat  = NY,                            &
                       nLev   = NZ,                     nTime = 1,                             &
                       fID    = fID,    lonID = lonID,         &
                       latID  = latID,  levID = levID,         &
                       timeID = timeID, VarCt = VarCt,         &
                       Create_NC4 = .true. )

       VarName = 'SpeciesConc_CO2'
       CALL Nc_Var_Def( fID         = fID,                             &
                        lonID       = lonID,                           &
                        latID       = latID,                           &
                        levID       = levID,                           &
                        TimeId      = timeID,                          &
                        VarName     = VarName,                         &
                        VarLongName = 'Dry mixing ratio of species CO2', &
                        VarUnit     = 'mol mol-1 dry',                 &
                        DataType=8,                                    &
                        VarCt=VarCT )

       CALL Nc_Var_Write( fID=fID, VarName=VarName, Arr4D=co2rst_varr4d )
       CALL Nccl(fID)

       IF (ASSOCIATED(co2rst_varr4d)) DEALLOCATE(co2rst_varr4d)
       END SUBROUTINE WRITE_SF_TO_SPCFILE 
!-----------------------------------------------------------------------------------
       SUBROUTINE WRITE_OBS_PERTD

       USE NCDF_MOD
       USE m_netcdf_io_checks,   ONLY : Ncdoes_Var_Exist
       USE m_netcdf_io_read,     ONLY : NcRd
       USE m_netcdf_io_close,      ONLY : Nccl
       USE ADJ_ARRAYS_MOD,       ONLY : OBS_PERTD

       CHARACTER(LEN=255)               :: FILENAME, TIMESTR, VarName
      ! dimension
       INTEGER                          :: NX, NY, NZ, nSpecies, nTimes, RC

       REAL(fp), POINTER                :: co2rst_varr3d(:, :, :, :)
       INTEGER                          :: co2rst_str3d(3)
       INTEGER                          :: co2rst_cnt3d(3)

       INTEGER                          :: fID, lonID, levID, latID, timeID, VarCT

       NX = 144
       NY = 91
       NZ = 1
       ALLOCATE(co2rst_varr3d(NX, NY, 1, 1), STAT = RC)

       co2rst_varr3d = 0.0

       co2rst_varr3d(:, :, 1, 1) = OBS_PERTD(:, :)
       TIMESTR = GET_MODEL_BEGIN_SIMULATION_TIME()


       FILENAME = 'PERTD/GEOSChem.SpeciesConc.'//TRIM(TIMESTR(1:13))//'z.nc4'
       CALL Nc_Create( NcFile = FILENAME,               Title = 'adjoint parameters',          &
                       nLon   = NX,                     nLat  = NY,                            &
                       nLev   = NZ,                     nTime = 1,                             &
                       fID    = fID,    lonID = lonID,         &
                       latID  = latID,  levID = levID,         &
                       timeID = timeID, VarCt = VarCt,         &
                       Create_NC4 = .true. )

       VarName = 'PERTD'
       CALL Nc_Var_Def( fID         = fID,                             &
                        lonID       = lonID,                           &
                        latID       = latID,                           &
                        levID       = levID,                           &
                        TimeId      = timeID,                          &
                        VarName     = VarName,                         &
                        VarLongName = 'Dry mixing ratio of species CO2', &
                        VarUnit     = 'mol mol-1 dry',                 &
                        DataType=8,                                    &
                        VarCt=VarCT )

       CALL Nc_Var_Write( fID=fID, VarName=VarName, Arr4D=co2rst_varr3d )
       CALL Nccl(fID)

       IF (ASSOCIATED(co2rst_varr3d)) DEALLOCATE(co2rst_varr3d)
       END SUBROUTINE WRITE_OBS_PERTD

!-----------------------------------------------------------------------------------
       FUNCTION GET_MODEL_BEGIN_SIMULATION_TIME( ) RESULT(TIMESTR)
       INTEGER                          :: LINE, LINE_NUM
       CHARACTER(LEN=255)               :: LINESTR, TIMESTR

       
       OPEN(UNIT=20, FILE='input.geos')
       DO LINE = 1, 3
          READ(20, '(a)') LINESTR
          print*, 'linestr', linestr
          IF (LINESTR(1:5) == "Start") THEN
             DO LINE_NUM = 1, 100
               ! if (str(j) == ":" and str(j+1) /= " ") then
                !   write(*,*) 'no blank in input geos Start Time:'
                IF ( LINESTR(LINE_NUM:LINE_NUM) == ":" .And. &
                     LINESTR(LINE_NUM+1:LINE_NUM+1) == " ") THEN
                     TIMESTR = LINESTR(LINE_NUM+2:LINE_NUM+9)//'_'//&
                               LINESTR(LINE_NUM+11:LINE_NUM+18)
!                ELSEIF  ( LINESTR(LINE_NUM:LINE_NUM) == ":" .And. &
!                     LINESTR(LINE_NUM+1:LINE_NUM+1) /= " ") THEN
!                    PRINT*, 'ERROR, addding black in input_geos begin time'
                ENDIF
             ENDDO
          ENDIF
      ENDDO
      CLOSE(20)
      END FUNCTION GET_MODEL_BEGIN_SIMULATION_TIME


      SUBROUTINE WRITE_EXIT_FLAG(FLAG)

      INTEGER, INTENT(IN)        :: FLAG 

      OPEN(UNIT=1995, FILE='exit.flag')
      WRITE(1995,*) FLAG
      CLOSE(1995)
 
      END SUBROUTINE WRITE_EXIT_FLAG


     SUBROUTINE GET_SPECIES0( State_Chm_Adj, State_Grid_Adj )

      USE NCDF_MOD
      USE m_netcdf_io_checks,   ONLY : Ncdoes_Var_Exist
      USE m_netcdf_io_read,     ONLY : NcRd
      USE State_Chm_Adj_Mod,    ONLY : ChmStateAdj
      USE State_Grid_Adj_Mod,   ONLY : GrdStateAdj
      USE m_netcdf_io_close,      ONLY : Nccl
      USE LOGICAL_ADJ_MOD,      ONLY : LICS_INC

      CHARACTER(LEN=255)              FILENAME


      ! dimension
      INTEGER                          :: NX, NY, NZ, nSpecies, nTimes
      !related nc file read
      REAL*8, ALLOCATABLE              :: co2_varrd_4dr(:, :, :, :)
      INTEGER                          :: co2_str4d(4)
      INTEGER                          :: co2_cnt4d(4)

      INTEGER                          :: FID


      LOGICAL                          :: err_stop
      INTEGER                          :: stat
      INTEGER                          :: RC
      CHARACTER(LEN=255)               :: varname

      TYPE(GrdStateAdj), INTENT(IN)    :: State_Grid_Adj  ! Grid State object
      TYPE(ChmStateAdj), INTENT(INOUT) :: State_Chm_Adj   ! Chemistry State object


      NX = State_Grid_Adj%NX
      NY = State_Grid_Adj%NY
      NZ = State_Grid_Adj%NZ
      nTimes = 1
 
      IF ( LICS_INC ) THEN
      FILENAME= '/home/huoxiao/GEOS_FP/GEOS_FP_ADJOINT/data/GEOSCHEM_RESTARTS/v2018-11/' &
                   //'initial_GEOSChem_rst.2x25_CO2.nc'
      print*, 'READ SPECIES0', FILENAME
      ELSE
      FILENAME= '/home/huoxiao/GEOS_FP/GEOS_FP_ADJOINT/data/GEOSCHEM_RESTARTS/v2018-11/' &
                   //'initial_GEOSChem_rst.2x25_CO2_background.nc'
      print*, 'READ SPECIES0, ', FILENAME
      ENDIF

      CALL NC_OPEN( TRIM(FILENAME), fID )
      ALLOCATE( co2_varrd_4dr(NX, NY, NZ, nTimes), STAT = RC )
      co2_str4d(1) = 1
      co2_str4d(2) = 1
      co2_str4d(3) = 1
      co2_str4d(4) = 1

      co2_cnt4d(1) = State_Grid_Adj%NX
      co2_cnt4d(2) = State_Grid_Adj%NY
      co2_cnt4d(3) = State_Grid_Adj%NZ
      co2_cnt4d(4) = 1

      varname = 'SpeciesRst_CO2'
      CALL NcRd(co2_varrd_4dr, FID, TRIM(varname), co2_str4d, co2_cnt4d)

      State_Chm_Adj%Species0(:, :, :, :) = co2_varrd_4dr(:, :, :, :)
      CALL Nccl(fID)
 
      END SUBROUTINE GET_SPECIES0
!-----------------------------------------------------------------------------------------
     SUBROUTINE GET_BEC( State_Chm_Adj, State_Grid_Adj )

      USE NCDF_MOD
      USE m_netcdf_io_checks,   ONLY : Ncdoes_Var_Exist
      USE m_netcdf_io_read,     ONLY : NcRd
      USE State_Chm_Adj_Mod,    ONLY : ChmStateAdj
      USE State_Grid_Adj_Mod,   ONLY : GrdStateAdj
      USE m_netcdf_io_close,      ONLY : Nccl


      CHARACTER(LEN=255)              FILENAME


      ! dimension
      INTEGER                          :: NX, NY, NZ, nSpecies, nTimes
      !related nc file read
      REAL*8, ALLOCATABLE              :: D_varrd_4dr(:, :, :, :)
      INTEGER                          :: D_str4d(4)
      INTEGER                          :: D_cnt4d(4)

      REAL*8, ALLOCATABLE              :: D_inv_varrd_4dr(:, :, :, :)
      INTEGER                          :: D_inv_str4d(4)
      INTEGER                          :: D_inv_cnt4d(4)

      REAL*8, ALLOCATABLE              :: Sx_varrd_2dr(:, :)
      INTEGER                          :: Sx_str2d(2)
      INTEGER                          :: Sx_cnt2d(2)

      REAL*8, ALLOCATABLE              :: Sx_inv_varrd_2dr(:, :)
      INTEGER                          :: Sx_inv_str2d(2)
      INTEGER                          :: Sx_inv_cnt2d(2)

      REAL*8, ALLOCATABLE              :: Sy_varrd_2dr(:, :)
      INTEGER                          :: Sy_str2d(2)
      INTEGER                          :: Sy_cnt2d(2)

      REAL*8, ALLOCATABLE              :: Sy_inv_varrd_2dr(:, :)
      INTEGER                          :: Sy_inv_str2d(2)
      INTEGER                          :: Sy_inv_cnt2d(2)

      REAL*8, ALLOCATABLE              :: Sz_varrd_2dr(:, :)
      INTEGER                          :: Sz_str2d(2)
      INTEGER                          :: Sz_cnt2d(2)

      REAL*8, ALLOCATABLE              :: Sz_inv_varrd_2dr(:, :)
      INTEGER                          :: Sz_inv_str2d(2)
      INTEGER                          :: Sz_inv_cnt2d(2)

      INTEGER                          :: FID

      LOGICAL                          :: err_stop
      INTEGER                          :: stat
      INTEGER                          :: RC
      CHARACTER(LEN=255)               :: varname

      TYPE(GrdStateAdj), INTENT(IN)    :: State_Grid_Adj  ! Grid State object
      TYPE(ChmStateAdj), INTENT(INOUT) :: State_Chm_Adj   ! Chemistry State object


      NX = State_Grid_Adj%NX
      NY = State_Grid_Adj%NY
      NZ = State_Grid_Adj%NZ
      nTimes = 1


      FILENAME= 'BEC.nc4'

      CALL NC_OPEN( TRIM(FILENAME), fID )
      ALLOCATE( D_varrd_4dr(NX, NY, NZ, nTimes), STAT = RC )
      D_str4d(1) = 1
      D_str4d(2) = 1
      D_str4d(3) = 1
      D_str4d(4) = 1

      D_cnt4d(1) = State_Grid_Adj%NX
      D_cnt4d(2) = State_Grid_Adj%NY
      D_cnt4d(3) = State_Grid_Adj%NZ
      D_cnt4d(4) = 1

      varname = 'D'
      CALL NcRd(D_varrd_4dr, FID, TRIM(varname), D_str4d, D_cnt4d)

      State_Chm_Adj%D(:, :, :, :) = D_varrd_4dr(:, :, :, :)

      !D_inv
      CALL NcRd(D_varrd_4dr, FID, TRIM(varname), D_str4d, D_cnt4d)

      State_Chm_Adj%D(:, :, :, :) = D_varrd_4dr(:, :, :, :)

      !D_inv
      ALLOCATE( D_inv_varrd_4dr(NX, NY, NZ, nTimes), STAT = RC )
      D_inv_str4d(1) = 1
      D_inv_str4d(2) = 1
      D_inv_str4d(3) = 1
      D_inv_str4d(4) = 1

      D_inv_cnt4d(1) = State_Grid_Adj%NX
      D_inv_cnt4d(2) = State_Grid_Adj%NY
      D_inv_cnt4d(3) = State_Grid_Adj%NZ
      D_inv_cnt4d(4) = 1

      varname = 'D_inv'
      CALL NcRd(D_inv_varrd_4dr, FID, TRIM(varname), D_inv_str4d, D_inv_cnt4d)

      State_Chm_Adj%D_inv(:, :, :, :) = D_inv_varrd_4dr(:, :, :, :)
!      State_Chm_Adj%D_inv(:, :, :, :) = 100000000
      !Sx
      ALLOCATE( Sx_varrd_2dr(NX, NX), STAT = RC )
      Sx_str2d(1) = 1
      Sx_str2d(2) = 1

      Sx_cnt2d(1) = State_Grid_Adj%NX
      Sx_cnt2d(2) = State_Grid_Adj%NX

      varname = 'Sx_h'
      CALL NcRd(Sx_varrd_2dr, FID, TRIM(varname), Sx_str2d, Sx_cnt2d)
      State_Chm_Adj%Sx(:, :) = Sx_varrd_2dr(:, :)

      !Sx_inv
      ALLOCATE( Sx_inv_varrd_2dr(NX, NX), STAT = RC )
      Sx_inv_str2d(1) = 1
      Sx_inv_str2d(2) = 1

      Sx_inv_cnt2d(1) = State_Grid_Adj%NX
      Sx_inv_cnt2d(2) = State_Grid_Adj%NX

      varname = 'Sx_inv'
      CALL NcRd(Sx_inv_varrd_2dr, FID, TRIM(varname), Sx_inv_str2d, Sx_inv_cnt2d)
      State_Chm_Adj%Sx_inv(:, :) = Sx_inv_varrd_2dr(:, :)

      print*, 'huoxiao_debug write Sx_inv'
      OPEN(UNIT=1995, FILE='Sx')
      WRITE(1995, *) State_Chm_Adj%Sx_inv 
      CLOSE(1995)
      !Sy

      ALLOCATE( Sy_varrd_2dr(NY, NY), STAT = RC )
      Sy_str2d(1) = 1
      Sy_str2d(2) = 1

      Sy_cnt2d(1) = State_Grid_Adj%NY
      Sy_cnt2d(2) = State_Grid_Adj%NY

      varname = 'Sy_h'
      CALL NcRd(Sy_varrd_2dr, FID, TRIM(varname), Sy_str2d, Sy_cnt2d)
      State_Chm_Adj%Sy(:, :) = Sy_varrd_2dr(:, :)

      !Sy_inv

      ALLOCATE( Sy_inv_varrd_2dr(NY, NY), STAT = RC )
      Sy_inv_str2d(1) = 1
      Sy_inv_str2d(2) = 1

      Sy_inv_cnt2d(1) = State_Grid_Adj%NY
      Sy_inv_cnt2d(2) = State_Grid_Adj%NY

      varname = 'Sy_inv'
      CALL NcRd(Sy_inv_varrd_2dr, FID, TRIM(varname), Sy_inv_str2d, Sy_inv_cnt2d)
      State_Chm_Adj%Sy_inv(:, :) = Sy_inv_varrd_2dr(:, :)
      OPEN(UNIT=1995, FILE='Sy')
      WRITE(1995, *) State_Chm_Adj%Sy_inv
      CLOSE(1995)
 
      !Sz
      ALLOCATE( Sz_varrd_2dr(NZ, NZ), STAT = RC )
      Sz_str2d(1) = 1
      Sz_str2d(2) = 1

      Sz_cnt2d(1) = State_Grid_Adj%NZ
      Sz_cnt2d(2) = State_Grid_Adj%NZ

      varname = 'Sz_h'
      CALL NcRd(Sz_varrd_2dr, FID, TRIM(varname), Sz_str2d, Sz_cnt2d)
      State_Chm_Adj%Sz(:, :) = Sz_varrd_2dr(:, :)

      !Sz_inv
      ALLOCATE( Sz_inv_varrd_2dr(NZ, NZ), STAT = RC )
      Sz_inv_str2d(1) = 1
      Sz_inv_str2d(2) = 1

      Sz_inv_cnt2d(1) = State_Grid_Adj%NZ
      Sz_inv_cnt2d(2) = State_Grid_Adj%NZ

      varname = 'Sz_inv'
      CALL NcRd(Sz_inv_varrd_2dr, FID, TRIM(varname), Sz_inv_str2d, Sz_inv_cnt2d)
      State_Chm_Adj%Sz_inv(:, :) = Sz_inv_varrd_2dr(:, :)
      OPEN(UNIT=1995, FILE='Sz')
      WRITE(1995, *) State_Chm_Adj%Sz_inv
      CLOSE(1995)

      CALL Nccl(fID)

      END SUBROUTINE GET_BEC

      SUBROUTINE WRITE_N_CALC

      USE ADJ_ARRAYS_MOD,        ONLY : N_CALC
     
      OPEN(UNIT=1995, FILE='N_CALC')
      WRITE(1995,*) N_CALC
      CLOSE(1995)
  
      END SUBROUTINE WRITE_N_CALC

      SUBROUTINE READ_TOT_ITER

      USE ADJ_ARRAYS_MOD,        ONLY : TOT_ITER

      OPEN(UNIT=1995, FILE='TOT_ITER')
      READ(1995,*) TOT_ITER
      CLOSE(1995)

      END SUBROUTINE READ_TOT_ITER

      SUBROUTINE READ_OUT_ITER

      USE ADJ_ARRAYS_MOD,        ONLY : OUT_ITER

      OPEN(UNIT=1995, FILE='OUT_ITER')
      READ(1995,*) OUT_ITER
      CLOSE(1995)

      END SUBROUTINE READ_OUT_ITER

      SUBROUTINE READ_REG_PARAM(REG_PARAM)

      REAL(FP),    INTENT(INOUT) :: REG_PARAM

      OPEN(UNIT=1995, FILE='REG_PARAM')
      READ(1995, *) REG_PARAM
      CLOSE(1995)
      print*, 'HUOXIAO_DEBUG READ', REG_PARAM

      END SUBROUTINE READ_REG_PARAM
End Module File_Ope_Mod


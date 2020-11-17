!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical  Adjoint Transport Model                  !
!------------------------------------------------------------------------------

MODULE State_Met_Adj_Mod

  USE Precision_Mod
  USE ErrCode_Adj_Mod

  IMPLICIT NONE
  PRIVATE

!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: Init_State_Met_Adj
  PUBLIC :: Cleanup_State_Met_Adj

  TYPE, PUBLIC :: MetStateAdj
  
     !----------------------------------------------------------------------
     ! Surface fields
     !----------------------------------------------------------------------

     !----------------------------------------------------------------------
     ! 1-D Fields
     !----------------------------------------------------------------------
      REAL(fp), ALLOCATABLE :: AP(:)                  ! "A" term for hybrid grid
      REAL(fp), ALLOCATABLE :: BP(:)                  ! "B" term for hybrid grid
    
       
     !----------------------------------------------------------------------
     ! Pbl fields
     !----------------------------------------------------------------------
     INTEGER,  POINTER :: IMIX         (:,:)    !adjoint 
     REAL(fp), POINTER :: FPBL         (:,:)    !adjoint  
     REAL(fp), POINTER :: SURFACE_PRESSURE(:, :)
     REAL(fp), POINTER :: SPECIFIC_HUMIDITY(:, :, :)
     
!----------------------------------------------------------------------
     ! 3-D Fields
     !----------------------------------------------------------------------
     REAL(fp), POINTER :: CMFMC         (:,:,:) ! Cloud mass flux [kg/m2/s]
     REAL(fp), POINTER :: DQRCU         (:,:,:) ! Conv precip production rate
                                                !  [kg/kg/s] (assume per
                                                !  dry air)
     REAL(fp), POINTER :: DTRAIN        (:,:,:) ! Detrainment flux [kg/m2/s]
     REAL(fp), POINTER :: PFICU         (:,:,:) ! Dwn flux ice prec:conv
                                                !  [kg/m2/s]
     REAL(fp), POINTER :: PFLCU         (:,:,:) ! Dwn flux liq prec:conv
                                                !  [kg/m2/s]
     REAL(fp), POINTER :: AD         (:,:,:) ! adjoint mass

     REAL(fp), POINTER :: U            (:,:,:) ! wind speed
     REAL(fp), POINTER :: V            (:,:,:) ! wind speed

     REAL(fp), POINTER :: XMASS        (:,:,:) ! temporary variable for advection
     REAL(fp), POINTER :: YMASS        (:,:,:) ! temporary variable for advection
     REAL(fp), POINTER :: PSC2_DRY     (:, :) ! temporary variable for advection
     REAL(fp), POINTER :: PEDGE_DRY    (:, :) ! temporary variable for advection

     !----------------------------------------------------------------------
     ! Air quantities assigned in AIRQNT
     !----------------------------------------------------------------------
     REAL(fp), POINTER :: DELP_DRY      (:,:,:) ! Delta-P (dry) across box [hPa]
     REAL(fp), POINTER :: PMID          (:,:,:) ! Average wet air pressure [hPa]
                                                !  defined as arithmetic
                                                !  average of edge pressures
     REAL(fp), POINTER :: PSC2WET       (:,:) ! surface pressure [hPa]

     REAL(fp), POINTER :: CO2_FOSSIL_FLUX     (:, :) !co2 fossil emission value
     REAL(fp), POINTER :: CO2_OCEAN_FLUX     (:, :) !co2 ocean emission value
     REAL(fp), POINTER :: CO2_BIO_FLUX     (:, :) !co2 bio emission value
     REAL(fp), POINTER :: CO2_BIOFUEL_FLUX     (:, :) !co2 biofuel emission value

  END TYPE MetStateAdj

  CONTAINS 

  SUBROUTINE Init_State_Met_Adj( State_Grid_Adj, State_Met_Adj, RC )

    ! USE 
    USE State_Grid_Adj_Mod,     ONLY : GrdStateAdj
    TYPE(GrdStateAdj), INTENT(IN)    :: State_Grid_Adj  ! Grid State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(MetStateAdj), INTENT(INOUT) :: State_Met_Adj   ! Obj for meteorology state
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Return code

    ! Scalars
    INTEGER            :: LX, IM, JM, LM

    ! Strings
    CHARACTER(LEN=255) :: ThisLoc


    ThisLoc = ' -> Init_State_Met (in Adjoint/state_met_mod.F90)'
    RC = GC_SUCCESS
    ! Shorten grid parameters for readability
    IM = State_Grid_Adj%NX ! # latitudes
    JM = State_Grid_Adj%NY ! # longitudes
    LM = State_Grid_Adj%NZ ! # levels

    State_Met_Adj%IMIX           => NULL()
    State_Met_Adj%FPBL           => NULL()
    State_Met_Adj%CMFMC          => NULL()
    State_Met_Adj%DQRCU          => NULL()
    State_Met_Adj%DTRAIN         => NULL()
    State_Met_Adj%PFICU          => NULL()
    State_Met_Adj%PFLCU          => NULL()
    State_Met_Adj%AD             => NULL()
    State_Met_Adj%U              => NULL()
    State_Met_Adj%V              => NULL()
    State_Met_Adj%XMASS          => NULL()
    State_Met_Adj%YMASS          => NULL()
    State_Met_Adj%PSC2_DRY       => NULL()
    State_Met_Adj%PEDGE_DRY      => NULL()
    State_Met_Adj%DELP_DRY       => NULL()
    State_Met_Adj%PMID               => NULL()
    State_Met_Adj%PSC2WET            => NULL()
    State_Met_Adj%CO2_FOSSIL_FLUX       => NULL()
    State_Met_Adj%CO2_OCEAN_FLUX       => NULL()
    State_Met_Adj%CO2_BIO_FLUX       => NULL()
    State_Met_Adj%CO2_BIOFUEL_FLUX       => NULL()
    State_Met_Adj%SURFACE_PRESSURE   => NULL()
    State_Met_Adj%SPECIFIC_HUMIDITY  => NULL()


    ALLOCATE( State_Met_Adj%AP( LM+1 ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) THEN
       print*, ThisLoc//' Allocate AP Error'
    ENDIF      
    State_Met_Adj%AP = 1e+0_fp

    ALLOCATE( State_Met_Adj%BP( LM+1 ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) THEN
       print*, ThisLoc//' Allocate BP Error'
    ENDIF
    State_Met_Adj%BP = 0e+0_fp


    ALLOCATE( State_Met_Adj%IMIX( IM, JM ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) THEN
       print*, ThisLoc//' Allocate IMIX Error'
    ENDIF

    ALLOCATE( State_Met_Adj%FPBL( IM, JM ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) THEN
       print*, ThisLoc//' Allocate FPBL Error'
    ENDIF

    ALLOCATE( State_Met_Adj%CMFMC( IM, JM, LM+1 ), STAT=RC )  
    IF ( RC /= GC_SUCCESS ) THEN 
       print*, ThisLoc//' Allocate CMFMC Error'  
    ENDIF

    ALLOCATE( State_Met_Adj%DQRCU( IM, JM, LM ), STAT=RC )  
    IF ( RC /= GC_SUCCESS ) THEN 
       print*, ThisLoc//' Allocate DQRCU Error'  
    ENDIF

    ALLOCATE( State_Met_Adj%DTRAIN( IM, JM, LM ), STAT=RC )  
    IF ( RC /= GC_SUCCESS ) THEN 
       print*, ThisLoc//' Allocate DTRAIN Error'  
    ENDIF

    ALLOCATE( State_Met_Adj%PFICU( IM, JM, LM+1 ), STAT=RC )  
    IF ( RC /= GC_SUCCESS ) THEN 
       print*, ThisLoc//' Allocate PFICU Error'  
    ENDIF

    ALLOCATE( State_Met_Adj%PFLCU( IM, JM, LM+1 ), STAT=RC )  
    IF ( RC /= GC_SUCCESS ) THEN 
       print*, ThisLoc//' Allocate PFLCU Error'  
    ENDIF

    ALLOCATE( State_Met_Adj%AD( IM, JM, LM ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) THEN
       print*, ThisLoc//' Allocate AD Error'
    ENDIF

    ALLOCATE( State_Met_Adj%V( IM, JM, LM ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) THEN
       print*, ThisLoc//' Allocate V Error'
    ENDIF

    ALLOCATE( State_Met_Adj%U( IM, JM, LM ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) THEN
       print*, ThisLoc//' Allocate U Error'
    ENDIF

    ALLOCATE( State_Met_Adj%XMASS( IM, JM, LM ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) THEN
       print*, ThisLoc//' Allocate XMASS Error'
    ENDIF

    ALLOCATE( State_Met_Adj%YMASS( IM, JM, LM ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) THEN
       print*, ThisLoc//' Allocate YMASS Error'
    ENDIF

    ALLOCATE( State_Met_Adj%PSC2_DRY( IM, JM ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) THEN
       print*, ThisLoc//' Allocate PSC2_DRY Error'
    ENDIF

    ALLOCATE( State_Met_Adj%PEDGE_DRY( IM, JM), STAT=RC )
    IF ( RC /= GC_SUCCESS ) THEN
       print*, ThisLoc//' Allocate PEDGE_DRY Error'
    ENDIF

    ALLOCATE( State_Met_Adj%DELP_DRY( IM, JM, LM ), STAT=RC )  
    IF ( RC /= GC_SUCCESS ) THEN 
       print*, ThisLoc//' Allocate DELP_DRY Error'  
    ENDIF

    !-------------------------
    ! PMID [1]
    !-------------------------
    ALLOCATE( State_Met_Adj%PMID( IM, JM, LM ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) THEN
       print*, ThisLoc//' Allocate PMID Error'
    ENDIF

    ALLOCATE( State_Met_Adj%PSC2WET( IM, JM ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) THEN
       print*, ThisLoc//' Allocate PSC2WET Error'
    ENDIF

    ALLOCATE( State_Met_Adj%CO2_FOSSIL_FLUX( IM, JM ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) THEN
       print*, ThisLoc//' Allocate CO2_FOSSIL_FLUX Error'
    ENDIF

    ALLOCATE( State_Met_Adj%CO2_OCEAN_FLUX( IM, JM ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) THEN
       print*, ThisLoc//' Allocate CO2_OCEAN_FLUX Error'
    ENDIF

    ALLOCATE( State_Met_Adj%CO2_BIO_FLUX( IM, JM ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) THEN
       print*, ThisLoc//' Allocate CO2_BIO_FLUX Error'
    ENDIF

    ALLOCATE( State_Met_Adj%CO2_BIOFUEL_FLUX( IM, JM ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) THEN
       print*, ThisLoc//' Allocate CO2_BIOFUEL_FLUX Error'
    ENDIF

    ALLOCATE( State_Met_Adj%SURFACE_PRESSURE( IM, JM ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) THEN
       print*, ThisLoc//' Allocate SURFACE_PRESSURE Error'
    ENDIF

    ALLOCATE( State_Met_Adj%SPECIFIC_HUMIDITY( IM, JM, LM ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) THEN
       print*, ThisLoc//' Allocate SPECIFIC_HUMIDITY Error'
    ENDIF


    State_Met_Adj%AP = &
         (/ 0.000000d+00, 4.804826d-02, 6.593752d+00, 1.313480d+01, &
            1.961311d+01, 2.609201d+01, 3.257081d+01, 3.898201d+01, &
            4.533901d+01, 5.169611d+01, 5.805321d+01, 6.436264d+01, &
            7.062198d+01, 7.883422d+01, 8.909992d+01, 9.936521d+01, &
            1.091817d+02, 1.189586d+02, 1.286959d+02, 1.429100d+02, &
            1.562600d+02, 1.696090d+02, 1.816190d+02, 1.930970d+02, &
            2.032590d+02, 2.121500d+02, 2.187760d+02, 2.238980d+02, &
            2.243630d+02, 2.168650d+02, 2.011920d+02, 1.769300d+02, &
            1.503930d+02, 1.278370d+02, 1.086630d+02, 9.236572d+01, &
            7.851231d+01, 5.638791d+01, 4.017541d+01, 2.836781d+01, &
            1.979160d+01, 9.292942d+00, 4.076571d+00, 1.650790d+00, &
            6.167791d-01, 2.113490d-01, 6.600001d-02, 1.000000d-02 /)

      ! Bp [unitless] for 47 levels (48 edges)
    State_Met_Adj%BP = &
         (/ 1.000000d+00, 9.849520d-01, 9.634060d-01, 9.418650d-01, &
            9.203870d-01, 8.989080d-01, 8.774290d-01, 8.560180d-01, &
            8.346609d-01, 8.133039d-01, 7.919469d-01, 7.706375d-01, &
            7.493782d-01, 7.211660d-01, 6.858999d-01, 6.506349d-01, &
            6.158184d-01, 5.810415d-01, 5.463042d-01, 4.945902d-01, &
            4.437402d-01, 3.928911d-01, 3.433811d-01, 2.944031d-01, &
            2.467411d-01, 2.003501d-01, 1.562241d-01, 1.136021d-01, &
            6.372006d-02, 2.801004d-02, 6.960025d-03, 8.175413d-09, &
            0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00, &
            0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00, &
            0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00, &
            0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00 /)

  END SUBROUTINE Init_State_Met_Adj

!
  SUBROUTINE Cleanup_State_Met_Adj( State_Met_Adj, RC )
!

!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(MetStateAdj), INTENT(INOUT) :: State_Met_Adj   ! Obj for meteorology state
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Return code

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    !========================================================================
    ! Initialize
    !========================================================================
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> Cleanup_State_Met (in Headers/state_met_mod.F90)'


    IF ( ASSOCIATED( State_Met_Adj%IMIX ) ) THEN
       DEALLOCATE( State_Met_Adj%IMIX, STAT=RC )
       IF ( RC /= GC_SUCCESS ) print*, 'unable to deallocate IMIX'
       State_Met_Adj%IMIX => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met_Adj%FPBL ) ) THEN
       DEALLOCATE( State_Met_Adj%FPBL, STAT=RC )
       IF ( RC /= GC_SUCCESS ) print*, 'unable to deallocate IMIX'
       State_Met_Adj%FPBL => NULL()
    ENDIF

    !========================================================================
    ! Deallocate 3-D fields
    !========================================================================
    IF ( ASSOCIATED( State_Met_Adj%CMFMC ) ) THEN
       DEALLOCATE( State_Met_Adj%CMFMC, STAT=RC )
       IF ( RC /= GC_SUCCESS ) print*, 'unable to deallocate CMFMC'
       State_Met_Adj%CMFMC => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met_Adj%DQRCU ) ) THEN
       DEALLOCATE( State_Met_Adj%DQRCU, STAT=RC )
       IF ( RC /= GC_SUCCESS ) print*, 'unable to deallocate DQRCU'
       State_Met_Adj%DQRCU => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met_Adj%DTRAIN ) ) THEN
       DEALLOCATE( State_Met_Adj%DTRAIN, STAT=RC )
       IF ( RC /= GC_SUCCESS ) print*, 'unable to deallocate DTRAIN'
       State_Met_Adj%DTRAIN => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met_Adj%PFICU ) ) THEN
       DEALLOCATE( State_Met_Adj%PFICU, STAT=RC )
       IF ( RC /= GC_SUCCESS ) print*, 'unable to deallocate PFICU'
       State_Met_Adj%PFICU => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met_Adj%PFLCU ) ) THEN
       DEALLOCATE( State_Met_Adj%PFLCU, STAT=RC )
       IF ( RC /= GC_SUCCESS ) print*, 'unable to deallocate PFLCU'
       State_Met_Adj%PFLCU => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met_Adj%AD ) ) THEN
       DEALLOCATE( State_Met_Adj%AD, STAT=RC )
       IF ( RC /= GC_SUCCESS ) print*, 'unable to deallocate AD'
       State_Met_Adj%AD => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met_Adj%U ) ) THEN
       DEALLOCATE( State_Met_Adj%U, STAT=RC )
       IF ( RC /= GC_SUCCESS ) print*, 'unable to deallocate U'
       State_Met_Adj%U => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met_Adj%V ) ) THEN
       DEALLOCATE( State_Met_Adj%V, STAT=RC )
       IF ( RC /= GC_SUCCESS ) print*, 'unable to deallocate V'
       State_Met_Adj%V => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met_Adj%XMASS ) ) THEN
       DEALLOCATE( State_Met_Adj%XMASS, STAT=RC )
       IF ( RC /= GC_SUCCESS ) print*, 'unable to deallocate XMASS'
       State_Met_Adj%XMASS => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met_Adj%YMASS ) ) THEN
       DEALLOCATE( State_Met_Adj%YMASS, STAT=RC )
       IF ( RC /= GC_SUCCESS ) print*, 'unable to deallocate YMASS'
       State_Met_Adj%YMASS => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met_Adj%PSC2_DRY ) ) THEN
       DEALLOCATE( State_Met_Adj%PSC2_DRY, STAT=RC )
       IF ( RC /= GC_SUCCESS ) print*, 'unable to deallocate PSC2_DRY'
       State_Met_Adj%PSC2_DRY => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met_Adj%PEDGE_DRY ) ) THEN
       DEALLOCATE( State_Met_Adj%PEDGE_DRY, STAT=RC )
       IF ( RC /= GC_SUCCESS ) print*, 'unable to deallocate PEDGE_DRY'
       State_Met_Adj%PEDGE_DRY => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met_Adj%DELP_DRY ) ) THEN
       DEALLOCATE( State_Met_Adj%DELP_DRY, STAT=RC )
       IF ( RC /= GC_SUCCESS ) print*, 'unable to deallocate DELP_DRY'
       State_Met_Adj%DELP_DRY => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met_Adj%PMID ) ) THEN
       DEALLOCATE( State_Met_Adj%PMID, STAT=RC )
       IF ( RC /= GC_SUCCESS ) print*, 'unable to deallocate PMID'
       State_Met_Adj%PMID => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met_Adj%PSC2WET ) ) THEN
       DEALLOCATE( State_Met_Adj%PSC2WET, STAT=RC )
       IF ( RC /= GC_SUCCESS ) print*, 'unable to deallocate PSC2WET'
       State_Met_Adj%PSC2WET => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met_Adj%CO2_FOSSIL_FLUX ) ) THEN
       DEALLOCATE( State_Met_Adj%CO2_FOSSIL_FLUX, STAT=RC )
       IF ( RC /= GC_SUCCESS ) print*, 'unable to deallocate CO2_FOSSIL_FLUX'
       State_Met_Adj%CO2_FOSSIL_FLUX => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met_Adj%CO2_OCEAN_FLUX ) ) THEN
       DEALLOCATE( State_Met_Adj%CO2_OCEAN_FLUX, STAT=RC )
       IF ( RC /= GC_SUCCESS ) print*, 'unable to deallocate CO2_OCEAN_FLUX'
       State_Met_Adj%CO2_OCEAN_FLUX => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met_Adj%CO2_BIO_FLUX ) ) THEN
       DEALLOCATE( State_Met_Adj%CO2_BIO_FLUX, STAT=RC )
       IF ( RC /= GC_SUCCESS ) print*, 'unable to deallocate CO2_BIO_FLUX'
       State_Met_Adj%CO2_BIO_FLUX => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met_Adj%CO2_BIOFUEL_FLUX ) ) THEN
       DEALLOCATE( State_Met_Adj%CO2_BIOFUEL_FLUX, STAT=RC )
       IF ( RC /= GC_SUCCESS ) print*, 'unable to deallocate CO2_BIOFUEL_FLUX'
       State_Met_Adj%CO2_BIOFUEL_FLUX => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met_Adj%SURFACE_PRESSURE ) ) THEN
       DEALLOCATE( State_Met_Adj%SURFACE_PRESSURE, STAT=RC )
       IF ( RC /= GC_SUCCESS ) print*, 'unable to deallocate SURFACE_PRESSURE'
       State_Met_Adj%SURFACE_PRESSURE => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met_Adj%SPECIFIC_HUMIDITY ) ) THEN
       DEALLOCATE( State_Met_Adj%SPECIFIC_HUMIDITY , STAT=RC )
       IF ( RC /= GC_SUCCESS ) print*, 'unable to deallocate SPECIFIC_HUMIDITY '
       State_Met_Adj%SPECIFIC_HUMIDITY  => NULL()
    ENDIF


    IF( ALLOCATED(State_Met_Adj%AP) ) DEALLOCATE( State_Met_Adj%AP )

    IF( ALLOCATED(State_Met_Adj%BP) ) DEALLOCATE( State_Met_Adj%BP )
       

  END SUBROUTINE Cleanup_State_Met_Adj

END MODULE State_Met_Adj_Mod

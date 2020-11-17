

! !INTERFACE:
!
MODULE State_Chm_Adj_Mod

#     include "define_adj.h"

!
  USE ErrCode_Mod                        ! Error handling

  USE Precision_Mod                      ! GEOS-Chem precision types

  IMPLICIT NONE
  PRIVATE


! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: Init_State_Chm_Adj
  PUBLIC :: Cleanup_State_Chm_Adj
  REAL(fp), PUBLIC              :: EMS_SF(144,91, 1) = 1.0 ! Emission scale factor
  REAL(fp), PUBLIC              :: EMS_SF0(144,91, 1) ! Emission scale factor
  REAL(fp), PUBLIC              :: EMS_SF_ADJ(144, 91, 1) = 0.0_fp
!  REAL(fp), PUBLIC              :: f(37) = 1.0
!  REAL(fp), PUBLIC              :: g(37) = 0.0
!  REAL(fp), PUBLIC              :: gama_df = 1.0
  !=========================================================================
  ! Derived type for Chemistry State
  !=========================================================================
  TYPE, PUBLIC :: ChmStateAdj

     INTEGER                    :: nAdvect              ! # advected specie
     INTEGER                    :: nSpecies             ! #  species
     !----------------------------------------------------------------------
     ! Chemical species
     !----------------------------------------------------------------------
     REAL(fp), POINTER  :: Species       (:,:,:,:) ! Species concentration
     REAL(fp), POINTER  :: Species0      (:,:,:,:) ! Species concentration
     REAL(fp), POINTER  :: SpeciestB     (:,:,:,:) ! Background concentration integration in time t 
     REAL(fp), POINTER  :: CV            (:,:,:,:)
     REAL(fp), POINTER  :: CO2_CONV      (:,:,:,:)
     REAL(fp), POINTER  :: CO2_ADV       (:,:,:,:)
     REAL(fp), POINTER :: QDELQ      (:,:,:,:)
     REAL(fp), POINTER :: QDELQ1      (:,:,:,:)

     REAL(fp), POINTER  :: SpeciesAdj    (:,:,:,:) ! Species adjoint
     REAL(fp), POINTER  :: EMS_SF_ADJ    (:,:, :) ! Emission scale factor adjoint

     REAL(fp), POINTER  :: emMW_g        (:)


     REAL(fp), POINTER  :: ICS_ERROR     (:)
     REAL(fp), POINTER  :: REG_PARAM_ICS (:)
     REAL(fp), POINTER  :: EMS_ERROR     (:)
     REAL(fp), POINTER  :: REG_PARAM_EMS (:)

     !----------------------------------------------------------------------
     ! Background Error Covariance
     !----------------------------------------------------------------------
     REAL(fp), POINTER  :: D       (:,:,:, :)
     REAL(fp), POINTER  :: Sx      (:,:)
     REAL(fp), POINTER  :: Sy      (:,:)
     REAL(fp), POINTER  :: Sz      (:,:)

     REAL(fp), POINTER  :: D_inv   (:,:,:, :)
     REAL(fp), POINTER  :: Sx_inv  (:,:)
     REAL(fp), POINTER  :: Sy_inv  (:,:)
     REAL(fp), POINTER  :: Sz_inv  (:,:)

     !JC
     REAL(fp), POINTER  :: C      (:,:,:, :)
     REAL(fp), POINTER  :: f      (:)
     REAL(fp), POINTER  :: g      (:)
     REAL(fp), POINTER  :: gama_df(:)
  END TYPE ChmStateAdj

CONTAINS
  SUBROUTINE Init_State_Chm_Adj( State_Chm_Adj, State_Grid_Adj, RC )

    USE State_Grid_Adj_Mod,       ONLY : GrdStateAdj

    TYPE(GrdStateAdj), INTENT(IN)    :: State_Grid_Adj  ! Grid State object

    TYPE(ChmStateAdj), INTENT(INOUT) :: State_Chm_Adj   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Return code

    ! Scalars
    INTEGER                :: IM, JM, LM

    ! Strings
    CHARACTER(LEN=255)     :: ErrMsg, ThisLoc

    ! Error handling
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at Init_State_Chm (in Headers/state_chm_mod.F90)'

    ! Shorten grid parameters for readability
    IM                      =  State_Grid_Adj%NX ! # latitudes
    JM                      =  State_Grid_Adj%NY ! # longitudes
    LM                      =  State_Grid_Adj%NZ ! # levels

    ! Number of each type of species
    State_Chm_Adj%nSpecies       =  1
    State_Chm_Adj%nAdvect       =  1
#if  defined ( LOG_OPT )
    EMS_SF0 = 0.0_fp
#else
    EMS_SF0 = 1.0_fp
#endif
    ! Chemical species
    State_Chm_Adj%Species       => NULL()

    State_Chm_Adj%Species0      => NULL()
    State_Chm_Adj%CV            => NULL()
    State_Chm_Adj%SpeciestB     => NULL()
    State_Chm_Adj%SpeciesAdj    => NULL()
!    State_Chm_Adj%EMS_SF_ADJ    => NULL()
    State_Chm_Adj%emMW_g        => NULL()
    State_Chm_Adj%ICS_ERROR     => NULL()
    State_Chm_Adj%REG_PARAM_ICS => NULL()

    State_Chm_Adj%EMS_ERROR     => NULL()
    State_Chm_Adj%REG_PARAM_EMS => NULL()

!Background Error Covariance
    State_Chm_Adj%D => NULL()
    State_Chm_Adj%Sx => NULL()
    State_Chm_Adj%Sy => NULL()
    State_Chm_Adj%Sz => NULL()
    State_Chm_Adj%D_inv => NULL()
    State_Chm_Adj%Sx_inv => NULL()
    State_Chm_Adj%Sy_inv => NULL()
    State_Chm_Adj%Sz_inv => NULL()

!co2
    State_Chm_Adj%CO2_CONV  => NULL()
    State_Chm_Adj%CO2_ADV   => NULL()
    State_Chm_Adj%QDELQ     => NULL()
    State_Chm_Adj%QDELQ1     => NULL()

!Jc
    State_Chm_Adj%C => NULL()
    State_Chm_Adj%f => NULL()
    State_Chm_Adj%g => NULL()
    State_Chm_Adj%gama_df => NULL()

    ALLOCATE( State_Chm_Adj%Species( IM, JM, LM, State_Chm_Adj%nSpecies ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Chm_Adj%Species = 0.0_fp


    ALLOCATE( State_Chm_Adj%CV( IM, JM, LM, State_Chm_Adj%nSpecies ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Chm_Adj%CV = 0.0_fp

    ALLOCATE( State_Chm_Adj%Species0( IM, JM, LM, State_Chm_Adj%nSpecies ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Chm_Adj%Species0 = 0.0_fp

    ALLOCATE( State_Chm_Adj%SpeciestB( IM, JM, LM, State_Chm_Adj%nSpecies ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Chm_Adj%SpeciestB = 0.0_fp

    ALLOCATE( State_Chm_Adj%SpeciesAdj( IM, JM, LM, State_Chm_Adj%nSpecies ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Chm_Adj%SpeciesAdj = 0.0_fp

!    ALLOCATE( State_Chm_Adj%EMS_SF_ADJ( IM, JM, State_Chm_Adj%nSpecies ), STAT=RC )
!    IF ( RC /= GC_SUCCESS ) RETURN
!    State_Chm_Adj%EMS_SF_ADJ = 0.0_fp


    ALLOCATE( State_Chm_Adj%emMW_g( State_Chm_Adj%nSpecies ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Chm_Adj%emMW_g(1) = 44.0_fp

    ALLOCATE( State_Chm_Adj%ICS_ERROR( State_Chm_Adj%nSpecies ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Chm_Adj%ICS_ERROR(1) = 0.000001_fp**2

    ALLOCATE( State_Chm_Adj%REG_PARAM_ICS( State_Chm_Adj%nSpecies ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Chm_Adj%REG_PARAM_ICS(1) = 0.0008_fp

    ALLOCATE( State_Chm_Adj%EMS_ERROR( State_Chm_Adj%nSpecies ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Chm_Adj%EMS_ERROR(1) = 1.0_fp

    ALLOCATE( State_Chm_Adj%REG_PARAM_EMS( State_Chm_Adj%nSpecies ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Chm_Adj%REG_PARAM_EMS(1) = 0.01_fp


    ALLOCATE( State_Chm_Adj%D( IM, JM, LM, State_Chm_Adj%nSpecies ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Chm_Adj%D = 0.0_fp

    ALLOCATE( State_Chm_Adj%Sx( IM, IM ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Chm_Adj%Sx = 0.0_fp

    ALLOCATE( State_Chm_Adj%Sy( JM, JM ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Chm_Adj%Sy = 0.0_fp

    ALLOCATE( State_Chm_Adj%Sz( LM, LM ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Chm_Adj%Sz = 0.0_fp


    ALLOCATE( State_Chm_Adj%D_inv( IM, JM, LM, State_Chm_Adj%nSpecies ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Chm_Adj%D_inv = 0.0_fp

    ALLOCATE( State_Chm_Adj%Sx_inv( IM, IM ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Chm_Adj%Sx_inv = 0.0_fp

    ALLOCATE( State_Chm_Adj%Sy_inv( JM, JM ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Chm_Adj%Sy_inv = 0.0_fp

    ALLOCATE( State_Chm_Adj%Sz_inv( LM, LM ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Chm_Adj%Sz_inv = 0.0_fp

    ALLOCATE( State_Chm_Adj%CO2_CONV( IM, JM, LM, State_Chm_Adj%nSpecies ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Chm_Adj%CO2_CONV = 0.0_fp

    ALLOCATE( State_Chm_Adj%CO2_ADV( IM, JM, LM, State_Chm_Adj%nSpecies ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Chm_Adj%CO2_ADV = 0.0_fp

    ALLOCATE( State_Chm_Adj%QDELQ( IM, JM, LM, 2 ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Chm_Adj%QDELQ = 0.0_fp   

    ALLOCATE( State_Chm_Adj%QDELQ1( IM, JM, LM, 2 ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Chm_Adj%QDELQ1 = 0.0_fp

    ALLOCATE( State_Chm_Adj%C( IM, JM, LM, 1 ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Chm_Adj%C = 1.0_fp

    ALLOCATE( State_Chm_Adj%f(37), STAT=RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Chm_Adj%f = 0.5_fp

    ALLOCATE( State_Chm_Adj%g(37), STAT=RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Chm_Adj%g = 0.0_fp

    ALLOCATE( State_Chm_Adj%gama_df(1), STAT=RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Chm_Adj%gama_df = 1.0_fp

  END SUBROUTINE Init_State_Chm_Adj  

  SUBROUTINE Cleanup_State_Chm_Adj( State_Chm_Adj, RC )  

! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmStateAdj), INTENT(INOUT) :: State_Chm_Adj    ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC           ! Return code
    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc
    !=======================================================================
    ! Initialize
    !=======================================================================
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at Cleanup_State_Chm (in Headers/state_chm_mod.F90)'

    !=======================================================================
    ! Deallocate and nullify pointer fields of State_Chm
    !=======================================================================
    IF ( ASSOCIATED( State_Chm_Adj%Species ) ) THEN
       DEALLOCATE( State_Chm_Adj%Species, STAT=RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm_Adj%Species=> NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm_Adj%CV ) ) THEN
       DEALLOCATE( State_Chm_Adj%CV, STAT=RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm_Adj%CV=> NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm_Adj%SpeciestB ) ) THEN
       DEALLOCATE( State_Chm_Adj%SpeciestB, STAT=RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm_Adj%Species=> NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm_Adj%Species0 ) ) THEN
       DEALLOCATE( State_Chm_Adj%Species, STAT=RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm_Adj%Species0=> NULL()
    ENDIF

  IF ( ASSOCIATED( State_Chm_Adj%SpeciesAdj ) ) THEN
       DEALLOCATE( State_Chm_Adj%SpeciesAdj, STAT=RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm_Adj%SpeciesAdj=> NULL()
    ENDIF

!  IF ( ASSOCIATED( State_Chm_Adj%EMS_SF_ADJ ) ) THEN
!       DEALLOCATE( State_Chm_Adj%EMS_SF_ADJ, STAT=RC )
!       IF ( RC /= GC_SUCCESS ) RETURN
!       State_Chm_Adj%EMS_SF_ADJ=> NULL()
!    ENDIF


  IF ( ASSOCIATED( State_Chm_Adj%emMW_g ) ) THEN
       DEALLOCATE( State_Chm_Adj%emMW_g, STAT=RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm_Adj%emMW_g=> NULL()
    ENDIF


  IF ( ASSOCIATED( State_Chm_Adj%ICS_ERROR ) ) THEN
       DEALLOCATE( State_Chm_Adj%ICS_ERROR, STAT=RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm_Adj%ICS_ERROR=> NULL()
    ENDIF

  IF ( ASSOCIATED( State_Chm_Adj%REG_PARAM_ICS ) ) THEN
       DEALLOCATE( State_Chm_Adj%REG_PARAM_ICS, STAT=RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm_Adj%REG_PARAM_ICS=> NULL()
    ENDIF

  IF ( ASSOCIATED( State_Chm_Adj%EMS_ERROR ) ) THEN
       DEALLOCATE( State_Chm_Adj%EMS_ERROR, STAT=RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm_Adj%EMS_ERROR=> NULL()
    ENDIF

  IF ( ASSOCIATED( State_Chm_Adj%REG_PARAM_EMS ) ) THEN
       DEALLOCATE( State_Chm_Adj%REG_PARAM_EMS, STAT=RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm_Adj%REG_PARAM_EMS=> NULL()
    ENDIF


  IF ( ASSOCIATED( State_Chm_Adj%D ) ) THEN
       DEALLOCATE( State_Chm_Adj%D, STAT=RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm_Adj%D=> NULL()
   ENDIF

  IF ( ASSOCIATED( State_Chm_Adj%Sx ) ) THEN
       DEALLOCATE( State_Chm_Adj%Sx, STAT=RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm_Adj%Sx=> NULL()
  ENDIF

  IF ( ASSOCIATED( State_Chm_Adj%Sy ) ) THEN
       DEALLOCATE( State_Chm_Adj%Sy, STAT=RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm_Adj%Sy=> NULL()
  ENDIF

  IF ( ASSOCIATED( State_Chm_Adj%Sz ) ) THEN
       DEALLOCATE( State_Chm_Adj%Sz, STAT=RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm_Adj%Sz=> NULL()
  ENDIF

  IF ( ASSOCIATED( State_Chm_Adj%D_inv ) ) THEN
       DEALLOCATE( State_Chm_Adj%D_inv, STAT=RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm_Adj%D_inv=> NULL()
  ENDIF

  IF ( ASSOCIATED( State_Chm_Adj%Sx_inv ) ) THEN
       DEALLOCATE( State_Chm_Adj%Sx_inv, STAT=RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm_Adj%Sx_inv=> NULL()
  ENDIF

  IF ( ASSOCIATED( State_Chm_Adj%Sy_inv ) ) THEN
       DEALLOCATE( State_Chm_Adj%Sy_inv, STAT=RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm_Adj%Sy_inv=> NULL()
  ENDIF

  IF ( ASSOCIATED( State_Chm_Adj%Sz_inv ) ) THEN
       DEALLOCATE( State_Chm_Adj%Sz_inv, STAT=RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm_Adj%Sz_inv=> NULL()
  ENDIF


    IF ( ASSOCIATED( State_Chm_Adj%CO2_CONV ) ) THEN
       DEALLOCATE( State_Chm_Adj%CO2_CONV , STAT=RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm_Adj%CO2_CONV  => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm_Adj%CO2_ADV ) ) THEN
       DEALLOCATE( State_Chm_Adj%CO2_ADV , STAT=RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm_Adj%CO2_ADV  => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm_Adj%QDELQ ) ) THEN
       DEALLOCATE( State_Chm_Adj%QDELQ , STAT=RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm_Adj%QDELQ  => NULL()
    ENDIF


    IF ( ASSOCIATED( State_Chm_Adj%QDELQ1 ) ) THEN
       DEALLOCATE( State_Chm_Adj%QDELQ1 , STAT=RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm_Adj%QDELQ1  => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm_Adj%C ) ) THEN
       DEALLOCATE( State_Chm_Adj%C , STAT=RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm_Adj%C  => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm_Adj%f ) ) THEN
       DEALLOCATE( State_Chm_Adj%f , STAT=RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm_Adj%f  => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm_Adj%g ) ) THEN
       DEALLOCATE( State_Chm_Adj%g , STAT=RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm_Adj%g  => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm_Adj%gama_df ) ) THEN
       DEALLOCATE( State_Chm_Adj%gama_df , STAT=RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm_Adj%gama_df  => NULL()
    ENDIF
  END SUBROUTINE Cleanup_State_Chm_Adj

END MODULE State_Chm_Adj_Mod

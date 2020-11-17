! $Id: co2_adj_mod.f,v 1.3 2012/03/01 22:00:26 daven Exp $
      MODULE CO2_ADJ_MOD

      USE PRECISION_MOD    ! For GEOS-Chem Precision (fp)
!
!******************************************************************************
!  Module CO2_ADJ_MOD contains variables and routines used for the CO2
!  adjoint simulation. (dkh, 04/25/10)
!
!  Based on the forward module  (pns, bmy, 8/16/05, 9/27/06)
!
!  Module Variables:
!  ============================================================================
!
!  Module Procedures:
!  ============================================================================
!  (1 ) EMISSCO2_ADJ           : Adjoint of emits CO2 into individual tracers
!
!  GEOS-CHEM modules referenced by "co2_mod.f"
!  ============================================================================
!  (1 ) biomass_mod.f          : Module w/ routines for biomass burning
!  (2 ) bpch2_mod.f            : Module w/ routines for binary punch file I/O
!  (3 ) diag04_mod.f           : Module w/ routines for CO2 diagnostics
!  (4 ) directory_mod.f        : Module w/ GEOS-CHEM data & met field dirs
!  (5 ) error_mod.f            : Module w/ I/O error and NaN check routines
!  (6 ) file_mod.f             : Module w/ file unit numbers and error checks
!  (7 ) grid_mod.f             : Module w/ horizontal grid information
!  (8 ) logical_mod.f          : Module w/ GEOS-CHEM logical switches
!  (9 ) time_mod.f             : Module w/ routines for computing time & date
!  (10) tracer_mod.f           : Module w/ GEOS-CHEM tracer array STT etc.
!  (11) transfer_mod.f         : Module w/ routines to cast & resize arrays
!  (12) dao_mod.f              : Module w/ routines for working with DAO met fields
!  (13) regrid_1x1_mod.f       : Modele w/ routines for regridding to and from 1x1
!
!  NOTES:
!  (1 ) See forward model module for complete documentation
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables
      ! and routines from being seen outside "co2_mod.f"
      !=================================================================

      ! Make everything PRIVATE ...
      PRIVATE

      ! ... except these routines
      PUBLIC :: EMISSCO2_ADJ


      !=================================================================
      ! MODULE VARIABLES
      !=================================================================

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement
      !=================================================================
      CONTAINS

!-----------------------------------------------------------------------------

      SUBROUTINE EMISSCO2_ADJ( State_Chm_Adj, State_Grid_Adj, &
                               State_Met_Adj )
!
!******************************************************************************
!  Subroutine EMISSCO2_ADJ is the adjoint routine for CO2 emissions. (dkh, 04/25/10)
!
!  Based on forward model code.  (pns, bmy, 8/16/05, 9/27/06)
!
!  NOTES:
!
!******************************************************************************
!
      ! References to F90 modules
 
      ! adj_group
      USE State_Chm_Adj_Mod,      ONLY : ChmStateAdj
!      USE State_Diag_Mod,     ONLY : DgnState
      USE State_Grid_Adj_Mod,     ONLY : GrdStateAdj
      USE State_Met_Adj_Mod,      ONLY : MetStateAdj
      USE State_Chm_Adj_Mod,      ONLY : EMS_SF_ADJ 
!      USE ADJ_ARRAYS_MOD,         ONLY : IDADJ_ECO2ff, MMSCL
 
!      USE ADJ_ARRAYS_MOD,  ONLY : EMS_SF_ADJ, EXIST_OBS
!      USE ADJ_ARRAYS_MOD,  ONLY : STT_ADJ
      USE LOGICAL_ADJ_MOD, ONLY : LADJ

      USE LOGICAL_ADJ_MOD, ONLY : LADJ_BIO

!      USE LOGICAL_ADJ_MOD, ONLY : OCEAN
!      USE LOGICAL_ADJ_MOD, ONLY : BIO
!      USE LOGICAL_ADJ_MOD, ONLY : BIOFUEL
      ! time
!      USE TIME_MOD,        ONLY : GET_TS_EMIS
       USE TIME_MOD,         ONLY : GET_TS_DYN
!      USE TIME_MOD,        ONLY : GET_NYMD
!      USE TIME_MOD,        ONLY : GET_NHMS
!      USE TIME_MOD,        ONLY : GET_NYMDb
!      USE TIME_MOD,        ONLY : GET_NHMSb
!      USE TIME_MOD,        ONLY : GET_JD
!      USE TIME_MOD,        ONLY : GET_ELAPSED_SEC
!      USE JULDAY_MOD,      ONLY :  CALDATE
  
      !constant parameter
      USE PHYSCONSTANTS,   ONLY : g0_100


      ! Inout Parameter    

      ! convection adjoint for co2 
      TYPE(MetStateAdj), INTENT(IN)    :: State_Met_Adj   ! Meteorology State object
      TYPE(GrdStateAdj), INTENT(IN)    :: State_Grid_Adj  ! Grid State object
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(ChmStateAdj), INTENT(INOUT) :: State_Chm_Adj   ! Chemistry State object

      
      ! Arguments
      LOGICAL, SAVE          :: FIRST  = .TRUE.

      ! Local variables
      INTEGER                :: I,   J   
      INTEGER                :: M

      REAL*8                 :: DTSRCE, E_CO2
      REAL*8                 :: E_CO2_ADJ

!      REAL*8                :: EMFOSSCO2(IIPAR, JJPAR), EMOCECO2(IIPAR, JJPAR)
!      REAL*8                :: EMBIOCO2(IIPAR, JJPAR),  EMBIOFCO2(IIPAR, JJPAR!!)
      REAL(fp)               :: EMCO2(State_Grid_Adj%NX, State_Grid_Adj%NY)
!      REAL*8                :: DELP_DRY(IIPAR, JJPAR, LLPAR)    
      INTEGER                :: N_TRACERS 

      CHARACTER(LEN=255)     :: FILE_PREFIX
 
      INTEGER                :: NYMDb
      INTEGER                :: NHMSb 
      INTEGER                :: NYMD, NHMS
      REAL*8                 :: ELAPSED_SEC
      REAL*8                 :: JD0
      REAL*8                 :: JD1
  
      !
      !=================================================================
      ! EMISSCO2_ADJ begins here!
      !=================================================================

  
!--------------------------------------------------------------------------------------
!      ! Read monthly-mean biomass burning emissions
!      IF ( LBIOBRNCO2 .and. ITS_A_NEW_MONTH() ) THEN
!         CALL READ_MONTH_BIOBRN_CO2( MONTH, YEAR )
!      ENDIF
!  This requires a subroutine called READ_MONTH_BIOBRN_CO2
!  GFEDv2 biomass burning emissions are a better choice !Ray Nassar
!
!--------------------------------------------------------------------------------------
!  At present, biomass burning emissions are dealt with in the following way:
!
!  1) main.f calls do_emissions in emissions.f
!  2) do_emissions calls compute_biomass_emissions in biomass_mod.f
!  3a) compute_biomass_emissions calls gfed2_compute_biomass in gfed2_biomass_mod.f
!               ** OR **
!  3b) compute_biomass_emissions calls gc_read_biomass_co2 in gc_biomass_mod.f
!--------------------------------------------------------------------------------------


!-----------------------------------------------------------------------
! Fluxes with "possible" monthly variability are called below
! In some cases the annual file is just called at the start of the month
!-----------------------------------------------------------------------

      DTSRCE = GET_TS_DYN()
      N_TRACERS = State_Chm_Adj%nSpecies
      ! file time
!      NYMDb = GET_NYMDb()
!      NHMSb = GET_NHMSb()
!      ELAPSED_SEC = GET_ELAPSED_SEC()
      ! JD0: Astronomical Julian Date at start of GEOS-Chem run
!      JD0 = GET_JD( NYMDb, NHMSb )
      ! JD1: Astronomical Julian Date at current time
!      JD1 = JD0 + ( DBLE( ELAPSED_SEC-1200 ) / 86400e+0_f8 )
      ! Call CALDATE to compute the current YYYYMMDD and HHMMSS
!      CALL CALDATE( DBLE(JD1), NYMD, NHMS )
   
!open emission file to read
   
!      FILE_PREFIX= "./EMISSION/emission_fossil_"
!      CALL READ_FILE( FILE_PREFIX, EMFOSSCO2, IIPAR, JJPAR )
     
!      FILE_PREFIX= "./EMISSION/emission_ocean_"
!      CALL READ_FILE( FILE_PREFIX, EMOCECO2, IIPAR, JJPAR )

!      FILE_PREFIX= "./EMISSION/emission_bio_"
!      CALL READ_FILE( FILE_PREFIX, EMBIOCO2, IIPAR, JJPAR )

!      FILE_PREFIX= "./EMISSION/emission_biofuel_"
!      CALL READ_FILE( FILE_PREFIX, EMBIOFCO2, IIPAR, JJPAR )

 
! open state_met file to readEMFOSSCO2
!      FILE_PREFIX = "./State_Met/delp_dry_"
!      CALL READ_DELP_DRY( FILE_PREFIX, DELP_DRY )

      IF ( LADJ_BIO ) THEN
      EMCO2 = State_Met_Adj%CO2_BIO_FLUX(:, :)
      ENDIF
      !print*, 'huoxiao_debug bio flux', MAXVAL(EMCO2), MINVAL(EMCO2)
      !=================================================================
      ! Process emissions and save diagnostics
      !=================================================================


      !huoxiao adjoint now optimization for total emission 
      ! FOSSIL
!      IF ( FOSSIL ) THEN
!      print*, 'co2 fossil adj'
      DO M = 1, N_TRACERS
      ! Loop over latitudes

      DO J = 1, State_Grid_Adj%NY

      ! Loop over longitudes
      DO I = 1, State_Grid_Adj%NX

         !-------------------------------------------
         ! #1: Total CO2
         ! #2: CO2 from fossil fuel emissions
         !-------------------------------------------

            ! fwd code:
            !STT(I,J,1,1)   = STT(I,J,1,1) + E_CO2
            ! adj code:
            E_CO2_ADJ = State_Chm_Adj%SpeciesAdj(I,J,1,M)
            ! fwd code:
            !E_CO2          = E_CO2 * A_CM2 * DTSRCE / XNUMOL_CO2
            ! adj code:
!          E_CO2_ADJ      = E_CO2_ADJ * A_CM2 * DTSRCE / XNUMOL_CO2

            !huoxiao adj
            !  fwd code  kg m-2 s-1 to kg kg 
            ! E_CO2 = E_CO2*DTSRCE*( 1.0e+0_fp/ ( g0_100* State_Met%DELP_DRY) )
            ! adj code:
            E_CO2_ADJ      = E_CO2_ADJ * DTSRCE* & 
                          (1.0/(g0_100*State_Met_Adj%DELP_DRY(I, J, 1)))
             ! Fossil fuel emissions of CO2 kg/m2/s
            E_CO2          = EMCO2(I,J)

 
 

            ! adj_group: apply scaling factors (dkh, 04/25/10)

               ! fwd code:
               !E_CO2 = E_CO2 * EMS_SF(I,J,M,IDADJ_ECO2ff)
               ! adj code:
               EMS_SF_ADJ(I,J,M) = EMS_SF_ADJ(I,J,M) + E_CO2_ADJ * E_CO2
             IF (ABS(EMS_SF_ADJ(I, J, M))>0.2) THEN
             print*, 'huoxiao_debug_value', ABS(EMS_SF_ADJ(I, J, M)), &
                E_CO2_ADJ, EMCO2(I, J), State_Chm_Adj%SpeciesAdj(I,J,1,M), &
                (1.0/(g0_100*State_Met_Adj%DELP_DRY(I, J, 1)))
             ENDIF
        END DO
        END DO
        END DO
        print*, 'huoxiao debug max value of E_CO2_ADJ', &
                 MAXVAL(EMS_SF_ADJ(:,:,1))
        print*, 'huoxiao debug min value of E_CO2_ADJ', &
                 MINVAL(EMS_SF_ADJ(:,:,1))


!      ENDIF


      ! OCEAN
!      IF ( OCEAN ) THEN

!      print*, 'co2 ocean adj'

!      DO M = 1, MMSCL
      ! Loop over latitudes

!      DO J = 1, JJPAR

      ! Loop over longitudes
!      DO I = 1, IIPAR

         !-------------------------------------------
         ! #1: Total CO2
         ! #2: CO2 from ocean emissions
         !-------------------------------------------

            ! fwd code:
            !STT(I,J,1,1)   = STT(I,J,1,1) + E_CO2
            ! adj code:
!            E_CO2_ADJ = STT_ADJ(I,J,1,1)
            ! fwd code:
            !E_CO2          = E_CO2 * A_CM2 * DTSRCE / XNUMOL_CO2
            ! adj code:
!          E_CO2_ADJ      = E_CO2_ADJ * A_CM2 * DTSRCE / XNUMOL_CO2

            !huoxiao adj
            !  fwd code  kg m-2 s-1 to kg kg 
            ! E_CO2 = E_CO2*DTSRCE*( 1.0e+0_fp/ ( g0_100* State_Met%DELP_DRY) )
            ! adj code:
!            E_CO2_ADJ      = E_CO2_ADJ * DTSRCE*(1.0/(g0_100*DELP_DRY(I, J, 1)))
             ! ocean  emissions of CO2 kg/m2/s
!            E_CO2          = EMOCECO2(I,J)

            ! adj_group: apply scaling factors (dkh, 04/25/10)

               ! fwd code:
               !E_CO2 = E_CO2 * EMS_SF(I,J,M,IDADJ_ECO2ff)
               ! adj code:
!               EMS_SF_ADJ(I,J,M,IDADJ_ECO2ff) =    &
!     &            EMS_SF_ADJ(I,J,M,IDADJ_ECO2ff) + E_CO2_ADJ * E_CO2


!        END DO
!        END DO
!        END DO
!      ENDIF


      ! BIO
!      IF ( BIO ) THEN

!      print*, 'co2 bio adj'

!      DO M = 1, MMSCL
      ! Loop over latitudes

!      DO J = 1, JJPAR

      ! Loop over longitudes
!      DO I = 1, IIPAR

         !-------------------------------------------
         ! #1: Total CO2
         ! #2: CO2 from ocean emissions
         !-------------------------------------------

            ! fwd code:
            !STT(I,J,1,1)   = STT(I,J,1,1) + E_CO2
            ! adj code:
!            E_CO2_ADJ = STT_ADJ(I,J,1,1)
            ! fwd code:
            !E_CO2          = E_CO2 * A_CM2 * DTSRCE / XNUMOL_CO2
            ! adj code:
!          E_CO2_ADJ      = E_CO2_ADJ * A_CM2 * DTSRCE / XNUMOL_CO2

            !huoxiao adj
            !  fwd code  kg m-2 s-1 to kg kg 
            ! E_CO2 = E_CO2*DTSRCE*( 1.0e+0_fp/ ( g0_100* State_Met%DELP_DRY) )
            ! adj code:
!            E_CO2_ADJ      = E_CO2_ADJ * DTSRCE*(1.0/(g0_100*DELP_DRY(I, J, 1)))
             ! BIO emissions of CO2 kg/m2/s
!            E_CO2          = EMBIOCO2(I,J)

            ! adj_group: apply scaling factors (dkh, 04/25/10)

               ! fwd code:
               !E_CO2 = E_CO2 * EMS_SF(I,J,M,IDADJ_ECO2ff)
               ! adj code:
!               EMS_SF_ADJ(I,J,M,IDADJ_ECO2ff) =    &
!     &            EMS_SF_ADJ(I,J,M,IDADJ_ECO2ff) + E_CO2_ADJ * E_CO2

               
!        END DO
!        END DO
!        END DO
!      ENDIF
!      PRINT*, 'EMS_SF_ADJ(I,J,M,IDADJ_ECO2ff)', EMS_SF_ADJ(2,2,1,IDADJ_ECO2ff)
      ! BIOF
!      IF ( BIOFUEL ) THEN

!      print*, 'co2 biofuel adj'

!      DO M = 1, MMSCL
      ! Loop over latitudes  

!      DO J = 1, JJPAR

      ! Loop over longitudes
!      DO I = 1, IIPAR

         !-------------------------------------------
         ! #1: Total CO2
         ! #2: CO2 from ocean emissions
         !-------------------------------------------

            ! fwd code:
            !STT(I,J,1,1)   = STT(I,J,1,1) + E_CO2
            ! adj code:
!            E_CO2_ADJ = STT_ADJ(I,J,1,1)
            ! fwd code:
            !E_CO2          = E_CO2 * A_CM2 * DTSRCE / XNUMOL_CO2
            ! adj code:
!          E_CO2_ADJ      = E_CO2_ADJ * A_CM2 * DTSRCE / XNUMOL_CO2

            !huoxiao adj
            !  fwd code  kg m-2 s-1 to kg kg 
            ! E_CO2 = E_CO2*DTSRCE*( 1.0e+0_fp/ ( g0_100* State_Met%DELP_DRY) )
            ! adj code:
!            E_CO2_ADJ      = E_CO2_ADJ * DTSRCE*(1.0/(g0_100*DELP_DRY(I, J, 1)))
             !  BIOF emissions of CO2 kg/m2/s
!            E_CO2          = EMBIOFCO2(I,J)

            ! adj_group: apply scaling factors (dkh, 04/25/10)

               ! fwd code:
               !E_CO2 = E_CO2 * EMS_SF(I,J,M,IDADJ_ECO2ff)
               ! adj code:
!               EMS_SF_ADJ(I,J,M,IDADJ_ECO2ff) =    &
!     &            EMS_SF_ADJ(I,J,M,IDADJ_ECO2ff) + E_CO2_ADJ * E_CO2


!        END DO
!        END DO
!        END DO
!      ENDIF
      ! Return to calling program
      END SUBROUTINE EMISSCO2_ADJ

!-----------------------------------------------------------------------------


      ! End of module
      END MODULE CO2_ADJ_MOD

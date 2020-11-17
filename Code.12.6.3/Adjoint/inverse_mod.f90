MODULE INVERSE_MOD

!
    USE PRECISION_MOD

     IMPLICIT NONE

     REAL*8, ALLOCATABLE  :: X(:)
     REAL*8, ALLOCATABLE  :: GRADNT(:)
      CONTAINS

      SUBROUTINE CALC_NOPT( State_Grid_Adj )

!
!******************************************************************************
!  Subroutine CALC_NOPT calculates the number of paramteres to optimize
!
!  NOTES:
!  (1 ) Set NOPT for initial conditions to 3D: IIPAR*JJPAR*LLPAR*N_TRACERS to
!       be consistent with other parts of the code (mak, 6/18/09)
!  (2 ) Now support strat fluxes LADJ_STRAT (hml, dkh, 02/20/12, adj32_025)
!******************************************************************************
!
      ! References to F90 modules
      USE ADJ_ARRAYS_MOD,    ONLY : NOPT
      USE State_Grid_Adj_Mod,Only : GrdStateAdj

      USE LOGICAL_ADJ_MOD,   ONLY : LADJ_EMS, LICS

      TYPE(GrdStateAdj), INTENT(IN)  :: State_Grid_Adj

      !=================================================================
      ! CALC_NOPT begins here!
      !=================================================================

      ! if optimizing both initial emissions and initial conditions
 

      ! if optimizing emissions only
      IF ( LADJ_EMS ) THEN

      NOPT = State_Grid_Adj%NX * State_Grid_Adj%NY

      PRINT*, 'Max size of control vector is:', NOPT
      ENDIF

      IF ( LICS ) THEN
 
      NOPT = State_Grid_Adj%NX * State_Grid_Adj%NY* State_Grid_Adj%NZ

      PRINT*, 'Max size of control vector is:', NOPT

      ENDIF
      ! Return to calling program
      END SUBROUTINE CALC_NOPT

!------------------------------------------------------------------------------

      SUBROUTINE INIT_INVERSE
!
!******************************************************************************
!  Subroutine INIT_INVERSE initializes and zeros all allocatable arrays
!  declared in "inverse_mod.f" (dkh, 1/26/05)
!
!  NOTES:
!  (1 ) Now also allocate EMS_ICS_orig (dkh, 03/29/05)
!  (2 ) Now check for incompatible preproc. definitions and ACTIVE_VARS. (dkh, 10/17/06)
!  (3 ) Now support strat fluxes LADJ_STRAT (hml, dkh, 02/20/12, adj32_025)
!
!******************************************************************************
!
      ! References to F90 modules
      USE ADJ_ARRAYS_MOD,  ONLY : NOPT

!      USE ADJ_ARRAYS_MOD,  ONLY : MMSCL
      USE LOGICAL_ADJ_MOD, ONLY : LADJ


      ! Local variables
      LOGICAL, SAVE            :: IS_INIT = .FALSE.
      INTEGER                  :: AS

      !=================================================================
      ! INIT_INVERSE begins here!
      !=================================================================

      ! Return if we have already initialized
!      IF ( IS_INIT ) RETURN

      !fp
!      IF ( LADJ ) THEN
         !Allocate arrays
 !        ALLOCATE( GRADNT( NOPT ), STAT = AS )

  !    ENDIF

      ALLOCATE( X( NOPT ), STAT = AS )
      ALLOCATE( GRADNT( NOPT ), STAT = AS )


      END SUBROUTINE INIT_INVERSE

!-----------------------------------------------------------------------------

      SUBROUTINE SET_SF

      USE State_Chm_Adj_Mod,         Only : EMS_SF, EMS_SF0
      USE LOGICAL_ADJ_MOD,           ONLY : LADJ_EMS, LICS


      ! Add flags (hml, 02/23/12)
      IF ( LADJ_EMS )  THEN 

      EMS_SF0 (:,:,:)    = EMS_SF (:,:,:)
    
      ENDIF
      END SUBROUTINE SET_SF


!-----------------------------------------------------------------------------

      SUBROUTINE SET_LOG_SF

      USE State_Chm_Adj_Mod,         Only : EMS_SF, EMS_SF0
      USE LOGICAL_ADJ_MOD,           ONLY : LADJ_EMS, LICS

      if ( LADJ_EMS ) THEN

      
      EMS_SF             = LOG(1.0)
      EMS_SF0 (:,:,:)    = EMS_SF (:,:,:)

      ENDIF

      END SUBROUTINE SET_LOG_SF

!-----------------------------------------------------------------------------


      SUBROUTINE GET_SF_FROM_X(State_Grid_Adj, State_Chm_Adj)
!
!*****************************************************************************
!  Subroutine GET_SF_FROM_X compiles the array of scaling factors  from
!   the vector X.  (dkh, 9/16/04)
!
!  NOTES:
!  (1 ) Added ACTIVE_VARS == 'EMISSIONS' case.  (dkh, 11/27/04)
!  (2 ) Added ACTIVE_VARS == 'FDTEST' case. (dkh, 02/17/05)
!  (3 ) Rename to SET_LOG_SF, replace CMN_ADJ with adjoint_array_mod
!        (dkh, ks, mak, cs  06/07/09)
!  (4 ) Now support strat fluxes LADJ_STRAT (hml, dkh, 02/20/12, adj32_025)
!
!*****************************************************************************

      USE State_Chm_Adj_Mod,    ONLY : EMS_SF

      USE State_Grid_Adj_Mod,   Only : GrdStateAdj
      USE State_Chm_Adj_Mod,   Only : ChmStateAdj
      USE LOGICAL_ADJ_MOD,      ONLY : LADJ_EMS, LICS, LICS_INC
      USE ADJ_ARRAYS_MOD,       ONLY : OBS_EXIST

      TYPE(GrdStateAdj), INTENT(IN)  :: State_Grid_Adj
      TYPE(ChmStateAdj), INTENT(IN)  :: State_Chm_Adj

      ! Local Variables
      INTEGER                    :: I, J, L, M, N
      INTEGER                    :: I_DUM
      REAL                       :: TCVV_CO2

      TCVV_CO2 = 44/28.97

         IF( LADJ_EMS ) THEN
!$OMP PARALLEL DO &
!$OMP DEFAULT( SHARED ) &
!$OMP PRIVATE( J, I, I_DUM )
            DO J = 1, State_Grid_Adj%NY
            DO I = 1, State_Grid_Adj%NX

               I_DUM        = I + (  State_Grid_Adj%NX * ( J - 1)  )  
               ! Update the tracer concentrations from X
               EMS_SF(I,J,1) = X(I_DUM)

            ENDDO
            ENDDO
!$OMP END PARALLEL DO
         ENDIF


         IF( LICS .AND. (.NOT.LICS_INC) ) THEN
!$OMP PARALLEL DO &
!$OMP DEFAULT( SHARED ) &
!$OMP PRIVATE( L, J, I, I_DUM )
            DO L = 1, State_Grid_Adj%NZ
            DO J = 1, State_Grid_Adj%NY
            DO I = 1, State_Grid_Adj%NX

               I_DUM        = I + (  State_Grid_Adj%NX * ( J - 1)  )&
                                + State_Grid_Adj%NX*State_Grid_Adj%NY*(L-1)
               ! Update the tracer concentrations from X
               State_Chm_Adj%Species(I,J,L, 1) = X(I_DUM)/TCVV_CO2

            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO
         ENDIF

         IF( LICS_INC ) THEN
!$OMP PARALLEL DO &
!$OMP DEFAULT( SHARED ) &
!$OMP PRIVATE( L, J, I, I_DUM )
            DO L = 1, State_Grid_Adj%NZ
            DO J = 1, State_Grid_Adj%NY
            DO I = 1, State_Grid_Adj%NX
                  
               I_DUM        = I + (  State_Grid_Adj%NX * ( J - 1)  )&
                                + State_Grid_Adj%NX*State_Grid_Adj%NY*(L-1)
               ! Update the tracer concentrations from X
               State_Chm_Adj%CV(I,J,L, 1) = X(I_DUM)
!               IF ( OBS_EXIST(I, J) == 1) THEN
!                  PRINT*, 'CONTROL_W, I, J, W, I_DUM',I, J, X(I_DUM), I_DUM
!               ENDIF

!               IF (State_Chm_Adj%CV(I,J,L, 1)>=0.000000001) THEN
!                  PRINT*, 'I, J, W, I_DUM', I, J, X(I_DUM), I_DUM
!               ENDIF
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO
         ENDIF

     END SUBROUTINE GET_SF_FROM_X

!----------------------------------------------------------------------------

     SUBROUTINE GET_X_FROM_SF( State_Chm_Adj, State_Grid_Adj )

      ! Reference to f90 modules
      USE State_Chm_Adj_Mod,    ONLY : EMS_SF
      USE LOGICAL_ADJ_MOD,      ONLY : LADJ_EMS, LICS, LICS_INC

      USE State_Grid_Adj_Mod,   ONLY : GrdStateAdj
      USE State_Chm_Adj_Mod,    ONLY : ChmStateAdj



      ! Local variables
      INTEGER                    :: I, J, L, M, N
      INTEGER                    :: I_DUM

      TYPE(GrdStateAdj)                   :: State_Grid_Adj
      TYPE(ChmStateAdj)                   :: State_Chm_Adj

      IF ( LADJ_EMS ) THEN
!$OMP PARALLEL DO &
!$OMP DEFAULT( SHARED ) &
!$OMP PRIVATE( J, I, I_DUM ) 
         DO J = 1, State_Grid_Adj%NY
         DO I = 1, State_Grid_Adj%NX

            I_DUM    = I + (  State_Grid_Adj%NX * ( J - 1)  ) 

               ! Load X from active tracer concentrations
            X(I_DUM) = EMS_SF(I,J,1)

          ENDDO
          ENDDO
!$OMP END PARALLEL DO
      ENDIF
 
      IF ( LICS .AND. (.NOT. LICS_INC)) THEN
!$OMP PARALLEL DO &
!$OMP DEFAULT( SHARED ) &
!$OMP PRIVATE( L, J, I, I_DUM )
         DO L = 1, State_Grid_Adj%NZ         
         DO J = 1, State_Grid_Adj%NY
         DO I = 1, State_Grid_Adj%NX

            I_DUM    = I + State_Grid_Adj%NX * ( J - 1) &
                         + State_Grid_Adj%NX*State_Grid_Adj%NY*(L-1)

               ! Load X from active tracer concentrations
            X(I_DUM) = State_Chm_Adj%Species(I,J,L,1)
          ENDDO
          ENDDO
          ENDDO
!$OMP END PARALLEL DO
      ENDIF

      IF ( LICS_INC ) THEN
!$OMP PARALLEL DO &
!$OMP DEFAULT( SHARED ) &
!$OMP PRIVATE( L, J, I, I_DUM )
         DO L = 1, State_Grid_Adj%NZ
         DO J = 1, State_Grid_Adj%NY
         DO I = 1, State_Grid_Adj%NX
            I_DUM    = I + State_Grid_Adj%NX * ( J - 1) &
                         + State_Grid_Adj%NX*State_Grid_Adj%NY*(L-1)

               ! Load X from active tracer concentrations
            X(I_DUM) = State_Chm_Adj%CV(I,J,L,1)
          ENDDO
          ENDDO
          ENDDO
!$OMP END PARALLEL DO
      ENDIF

      ! Return to calling program
      END SUBROUTINE GET_X_FROM_SF
!--------------------------------------------------------------------------------------
      FUNCTION PRESSURE_WEIGHTING_FUNCTION(edge_pressure, surface_pressure,&
                                          q,  level)
!     level = 1 is close to surface

      USE PARAMETER_ADJ_MOD,        ONLY : Mdry, g

      INTEGER, INTENT(IN)               :: level
      REAL(fp), INTENT(IN)              :: surface_pressure
!      REAL(fp), INTENT(IN)              :: hyai(level+1)
!      REAL(fp), INTENT(IN)              :: hybi(level+1)
      REAL(fp), INTENT(IN)              :: edge_pressure(level+1)
      REAL(fp), INTENT(IN)              :: q(level)
      !local
!      REAL(fp)                          :: out_pressure(level+1)
      REAL(fp)                          :: delta_pressure(level)
      REAL(fp)                          :: hdot(level)
      REAL(fp)                          :: air_sum
      REAL(fp)                          :: c(level)
      REAL(fp)                          :: fi, fs
      REAL(fp)                          :: pressure_weighting_function(level)
      INTEGER                           :: i
!
!******************************************************************************
!

!     Reference paper:


      do i = 1, level
         if (i == 1) then
            delta_pressure(i) = surface_pressure-edge_pressure(i+1)
         else
            delta_pressure(i) = edge_pressure(i)-edge_pressure(i+1)
         endif
      enddo
      air_sum = 0.0
      DO i = 1, level
         c(i) = (1-q(i+1))/(g*Mdry)
         air_sum = air_sum+c(i)*delta_pressure(i)
      enddo

      do i = 1, level
         pressure_weighting_function(i) = c(i)*delta_pressure(i)/air_sum
      enddo
!      fs = (surface_pressure-model_pressure(2))/(model_pressure(1)-model_pressure(2))
!      fi = 1/2

!      do i = 1, level
!         if (i == 1) then
!            pressure_weighting_function(i) = fs*fi*hdot(i)
!         elseif (i == 2) then
!            pressure_weighting_function(i) = fi*hdot(i)+(1-fs*fi)*hdot(i-1)
!         elseif (i == level) then
!            pressure_weighting_function(i) = (1-fi)*hdot(i-1)
!         else
!            pressure_weighting_function(i) = fi*hdot(i)+(1-fi)*hdot(i-1)
!         endif
!      enddo
      END FUNCTION PRESSURE_WEIGHTING_FUNCTION

!--------------------------------------------------------------------------------

      SUBROUTINE SET_OPT_RANGE( State_Chm_Adj, State_Grid_Adj )


      USE State_Grid_Adj_Mod,   Only : GrdStateAdj
      USE State_Chm_Adj_Mod,   Only : ChmStateAdj
 
      USE LOGICAL_ADJ_MOD,    ONLY : LICS
      USE LOGICAL_ADJ_MOD,    ONLY : LADJ_EMS  


      TYPE(GrdStateAdj), INTENT(IN)  :: State_Grid_Adj
      TYPE(ChmStateAdj), INTENT(IN)  :: State_Chm_Adj

      INTEGER                        :: I, J, L
      IF( LICS ) THEN 
      DO L = 1, State_Grid_Adj%NZ
      DO J = 16, 75
      DO I = 1, State_Grid_Adj%NX

      State_Chm_Adj%SpeciesAdj(I, J, L, 1) = 0.d0 

      ENDDO
      ENDDO
      ENDDO
      ENDIF
  

      END SUBROUTINE SET_OPT_RANGE 

      SUBROUTINE UPDATE_X0( State_Chm_Adj, State_Grid_Adj )

      USE State_Grid_Adj_Mod,   Only : GrdStateAdj
      USE State_Chm_Adj_Mod,   Only : ChmStateAdj
      USE BEC_Mod,                   Only : Uw
      USE Adj_Arrays_Mod,            Only : OBS_EXIST

      TYPE(ChmStateAdj), INTENT(INOUT)    :: State_Chm_Adj
      TYPE(GrdStateAdj), INTENT(INOUT)    :: State_Grid_Adj

      REAL*8                              :: TMP(144,91,47), TMP_CV(47), tmp_uw(47)
      INTEGER                             :: I, J, K


!      TMP = Uw( State_Chm_Adj%Sx(:, :), &
!                                                 State_Chm_Adj%Sy(:, :), &
!                                                 State_Chm_Adj%Sz(:, :), &
!                                                 State_Chm_Adj%D(:, :, :, 1), &
!                                                 State_Chm_Adj%CV(:, :, :, 1), &
!                                                 State_Grid_Adj%NX,      &
!                                                  State_Grid_Adj%NY,      &
!                                                 State_Grid_Adj%NZ)
!      tmp_cv = 0.0
!      tmp_uw = 0.0
!      DO J = 1, 91
!      DO I = 1, 144
!         IF ( OBS_EXIST(I, J) == 1) THEN
!            print*, 'CONTROL VARIABLE I, J, W', I, J, State_Chm_Adj%CV(I, J, :, 1)
!            print*, 'CONTROL VARIABLE UW', TMP(I,J,:)
!            tmp_cv(:) = tmp_cv(:)+State_Chm_Adj%CV(I, J, :, 1)       
!            tmp_uw(:) = tmp_uw(:)+TMP(I, J, :)
!         endif
!      enddo
!      enddo     
!      print*, 'SUM_CV AND SUM_UW', tmp_cv
!      print*, 'SUM_UW', TMP_UW 
      State_Chm_Adj%Species0(:, :, :, 1) = State_Chm_Adj%Species0(:, :, :, 1)+&
                                                 Uw( State_Chm_Adj%Sx(:, :), &
                                                 State_Chm_Adj%Sy(:, :), &
                                                 State_Chm_Adj%Sz(:, :), &
                                                 State_Chm_Adj%D(:, :, :, 1), &
                                                 State_Chm_Adj%CV(:, :, :, 1), &
                                                 State_Grid_Adj%NX,      &
                                                 State_Grid_Adj%NY,      &
                                                 State_Grid_Adj%NZ)
      END SUBROUTINE UPDATE_X0
END MODULE INVERSE_MOD

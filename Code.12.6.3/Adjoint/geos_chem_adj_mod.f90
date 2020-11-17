
!huoxiao adjoint


      MODULE GEOS_CHEM_ADJ_MOD

#     include "define_adj.h" 
      USE PRECISION_MOD

      IMPLICIT NONE


!
!******************************************************************************
!
!
!     GGGGGG   CCCCCC        A     DDDDD       J   OOO   I  N   N TTTTTTT
!    G        C             A A    D    D      J  O   O  I  NN  N    T
!    G   GGG  C        ==  AAAAA   D    D      J  0   O  I  N N N    T
!    G     G  C           A     A  D    D  J   J  0   O  I  N  NN    T
!     GGGGGG   CCCCCC    A       A DDDDD    JJJ    OOO   I  N   N    T
!
!
!           for CO2 2 x 2.5 global grids 
!
!
!******************************************************************************



      CONTAINS

      SUBROUTINE DO_GEOS_CHEM_ADJ( State_Chm_Adj, State_Grid_Adj, &
                                       State_Met_Adj, Input_Opt_Adj )

      USE State_Chm_Adj_Mod,         Only : ChmStateAdj
      USE State_Grid_Adj_Mod,        Only : GrdStateAdj
      USE State_Met_Adj_Mod,         Only : MetStateAdj

      USE Time_Mod,                  Only : Get_Elapsed_Sec
      USE Time_Mod,                  Only : Get_Ts_Dyn
      USE Time_Mod,                  Only : Get_Nhms
      USE Time_Mod,                  Only : Get_Nymd
      USE Time_Mod,                  Only : Set_Current_Time
      USE Time_Mod,                  Only : Its_Time_For_Conv
      USE Time_Mod,                  Only : Its_Time_For_Dyn
      USE Time_Mod,                  Only : Set_Elapsed_Sec_Adj

      USE Adj_Arrays_Mod,            Only : Its_Time_For_Checkpt       
      USE Adj_Arrays_Mod,            Only : Its_Time_For_Obs
      USE Adj_Arrays_Mod,            Only : Its_Time_For_I3_Field_Adj
      USE Adj_Arrays_Mod,            Only : Its_Time_For_Dyn_Adj
      USE Adj_Arrays_Mod,            Only : Its_Time_For_Emis_Adj
      USE Adj_Arrays_Mod,            Only : Convert_Units
      USE Adj_Arrays_Mod,            Only : Cost_Func
      USE Adj_Arrays_Mod,            Only : Cost_Func_Jo
      USE Adj_Arrays_Mod,            Only : N_CALC
      USE Adj_Arrays_Mod,            Only : Index_Array
      USE Adj_Arrays_Mod,            Only : Obs_Exist
      USE Adj_Arrays_Mod,            Only : Scale_Adjoint_Value

      USE LOGICAL_ADJ_MOD,           Only : LADJ_EMS
      USE LOGICAL_ADJ_MOD,           Only : LICS
      USE LOGICAL_ADJ_MOD,           Only : LICS_INC
      USE LOGICAL_ADJ_MOD,           Only : LSENS
      USE LOGICAL_ADJ_MOD,           Only : LTURB
      USE LOGICAL_ADJ_MOD,           Only : LCONV_ADJ
      USE LOGICAL_ADJ_MOD,           Only : LTRAN_ADJ
      USE LOGICAL_ADJ_MOD,           Only : LPBL_ADJ

      USE LOGICAL_ADJ_MOD,           Only : LCONV_FINITE_DIFF
      USE LOGICAL_ADJ_MOD,           Only : LTRAN_FINITE_DIFF
      USE LOGICAL_ADJ_MOD,           Only : LFORWARD_RUN

      USE File_Ope_Mod,              Only : Read_Met_Data_For_Emis_Time
      USE File_Ope_Mod,              Only : Read_Met_Data_For_Dyn_Time
      USE File_Ope_Mod,              Only : Read_Met_Data_For_I3_Time
      USE File_Ope_Mod,              Only : Read_Met_Data_For_A1_Time
      USE File_Ope_Mod,              Only : Read_Grid_Data
      USE File_Ope_Mod,              Only : Write_File_Init
      USE File_Ope_Mod,              Only : Get_Species0
      USE File_Ope_Mod,              Only : Read_Co2_Flux 
      USE File_Ope_Mod,              Only : Get_Bec

      USE Convection_Adj_Mod,        Only : Do_Convection_Adj
      USE Pbl_Mix_Adj_Mod,           Only : Do_Pbl_Mix_Adj
      USE Transport_Adj_Mod,         Only : Do_Transport_Adj
      USE Tpcore_Fvdas_Adj_Mod,      Only : Exit_Tpcore

      USE Input_Adj_Mod,             Only : OptInputAdj

      USE CO2_Adj_Mod,               Only : EMISSCO2_ADJ      

      USE Checkpt_Mod,               Only : Load_Checkpt_Data 

      USE COVARIANCE_MOD,            Only : CALC_COV_ERROR
!      USE BEC_MOD,                   Only : BinvB

      TYPE(ChmStateAdj), INTENT(INOUT)    :: State_Chm_Adj
      TYPE(GrdStateAdj), INTENT(INOUT)    :: State_Grid_Adj
      TYPE(MetStateAdj), INTENT(INOUT)    :: State_Met_Adj

      TYPE(OptInputAdj), INTENT(INOUT)    :: Input_Opt_Adj
 
      INTEGER      FINAL_ELAPSED_MIN, TS_DYN, RC
      INTEGER      FINAL_RUN_TIME
      INTEGER      MIN_ADJ
      INTEGER      NHMS, NYMD

      ! finite difference
      INTEGER                   :: Obs_Tot_Num
!      INTEGER,      Pointer     :: Index_Array(:, :) => Null()
      INTEGER      Cur_Obs_Num
      INTEGER      Index_I, Index_J

      LOGICAL      ADJ_FIRST
      WRITE(*,'(/, A100, /)')  &
          "***************** ADJOINT BEGIN ********************"
     
      ! BUG FIX: need to reset EMS_SF_ADJ so that gradients do not
      ! accumulate from one iteration to the next.

!      IF ( LADJ_EMS ) State_Chm_Adj%EMS_SF_ADJ = 0D0
      ! read obs index array and obs num
      CALL Read_Grid_Data( State_Grid_Adj )
!      CALL Read_Obs_Num( Obs_Tot_Num )
      IF ( .NOT. LICS_INC ) THEN
      CALL Get_Species0( State_Chm_Adj, State_Grid_Adj )
      CALL GET_BEC ( State_Chm_Adj, State_Grid_Adj )
!      CALL BinvB ( State_Chm_Adj%Sx, State_Chm_Adj%Sy, State_Chm_Adj%Sz, &
!                   State_Chm_Adj%Sx_inv, State_Chm_Adj%Sy_inv, State_Chm_Adj%Sz_inv)
!      STOP
      ENDIF
!      CALL Init_Index_Array( Index_Array, Obs_Tot_Num )
      ADJ_FIRST = .TRUE.

!      OPEN(UNIT=2020, FILE='forcing', ACCESS='APPEND')      
      Index_I = 126
      Index_J = 37

      FINAL_ELAPSED_MIN = GET_ELAPSED_SEC()
      ! begin time
      FINAL_RUN_TIME = 0
      TS_DYN = GET_TS_DYN()
      !set time
      CALL SET_ELAPSED_SEC_ADJ(6*TS_DYN, .TRUE., &
           FINAL_ELAPSED_MIN-6*TS_DYN, .TRUE.)
      CALL SET_CURRENT_TIME

      DO MIN_ADJ = FINAL_ELAPSED_MIN-6*TS_DYN, FINAL_RUN_TIME, -TS_DYN
       
 
         NHMS        = GET_NHMS()
         NYMD        = GET_NYMD()

         print*, 'ADJOINT CURRENT ITERATION IS ', MIN_ADJ
         
!         print '(/, A25)', 'New backward integration'
!         print *, 'nymd is', GET_NYMD(), 'nhms is',GET_NHMS()
         IF ( ITS_TIME_FOR_CHECKPT() .AND. (.NOT. LICS_INC)) THEN

            !load checkponint data 
            CALL LOAD_CHECKPT_DATA( State_Chm_Adj, State_Grid_Adj, State_Met_Adj )

         ENDIF
!         IF ( MOD(FINAL_ELAPSED_MIN, 1200) == 600 .AND. ADJ_FIRST ) THEN
                     !set time for I3 and A1 Field read
!         CALL SET_ELAPSED_SEC_ADJ(1*TS_DYN, .TRUE., 0, .FALSE.)
!         CALL SET_CURRENT_TIME
!         ENDIF
       
         ! read met data for dyn calculation             
!         IF ( Its_Time_For_I3_Field_Adj() .AND. Min_Adj == Final_Elapsed_Min ) THEN

!         CALL Read_Met_Data_For_I3_Time( State_Grid_Adj, State_Met_Adj )

!         ENDIF

!         IF ( MOD(FINAL_ELAPSED_MIN,1200) == 600 .AND. ADJ_FIRST ) THEN
                     !set time for I3 and A1 Field read
!         CALL SET_ELAPSED_SEC_ADJ(1*TS_DYN, .FALSE., 0, .FALSE.)
!         CALL SET_CURRENT_TIME
!         ENDIF
     
         IF ( Its_Time_For_Obs() ) THEN

         CALL Read_Met_Data_For_Emis_Time(State_Grid_Adj, State_Met_Adj ) 

         END IF
         ! read met data for adjoint
!         CALL Read_Met_Data( NYMD, NHMS, State_Grid_Adj, State_Met_Adj )
         IF ( Its_Time_For_Obs( ) .AND. (.NOT. LICS_INC)) THEN

                        ! Update cost function and calculate adjoint forcing

            ! for sensitivity calculations...
            IF ( LSENS ) THEN

               CALL CALC_ADJ_FORCE_FOR_SENS

            ! ... for cost functions involving observations (real or pseudo)
            ELSE
               ! temp
               CALL CALC_ADJ_FORCE_FOR_OBS( State_Chm_Adj, State_Grid_Adj, &
                                            State_Met_Adj )


            ENDIF

         ENDIF
!         print*, 'max SpeciesAdj OBS', &
!                 MAXVAL(State_Chm_Adj%SpeciesAdj)

         ! 0th time only for observation
         IF( Min_Adj == 0 ) EXIT

         IF ( MOD(FINAL_ELAPSED_MIN, 1200) == 600 .AND. ADJ_FIRST ) THEN
                     !set time for I3 and A1 Field read
         CALL SET_ELAPSED_SEC_ADJ(1*TS_DYN, .TRUE., 0, .FALSE.)
         CALL SET_CURRENT_TIME
         CALL Read_Met_Data_For_I3_Time( State_Grid_Adj, State_Met_Adj )

         CALL Read_Met_Data_For_A1_Time( State_Grid_Adj, State_Met_Adj )
         CALL SET_ELAPSED_SEC_ADJ(1*TS_DYN, .FALSE., 0, .FALSE.)
         CALL SET_CURRENT_TIME


         ! now we read met data for I3 and A1 in same frequency
         ELSEIF ( Its_Time_For_I3_Field_Adj() ) THEN
     
         !set time for I3 and A1 Field read
         CALL SET_ELAPSED_SEC_ADJ(2*TS_DYN, .TRUE., 0, .FALSE.)
         CALL SET_CURRENT_TIME

         CALL Read_Met_Data_For_I3_Time( State_Grid_Adj, State_Met_Adj )
        
         CALL Read_Met_Data_For_A1_Time( State_Grid_Adj, State_Met_Adj )

         !read co2 flux
         IF ( LADJ_EMS ) THEN
         CALL READ_CO2_FLUX( State_Grid_Adj, State_Met_Adj )
         ENDIF

         ! restore time 
         CALL SET_ELAPSED_SEC_ADJ(2*TS_DYN, .FALSE., 0, .FALSE.)
         CALL SET_CURRENT_TIME


         END IF

         ! reduce one dynamic step for convection and transportation calculation
         ! restore time?
         CALL SET_ELAPSED_SEC_ADJ(TS_DYN, .TRUE., 0, .FALSE.)
         CALL SET_CURRENT_TIME

         IF ( Its_Time_For_Dyn_Adj() ) THEN

         CALL Read_Met_Data_For_Dyn_Time( State_Chm_Adj, &
                                          State_Grid_Adj, &
                                          State_Met_Adj ) 

         END IF
         IF( Its_Time_For_Conv() ) THEN
            IF( LCONV_ADJ ) THEN

              CALL DO_CONVECTION_ADJ( State_Chm_Adj, State_Grid_Adj,     &
                                State_Met_Adj, RC )
            ENDIF
!         print*, 'max SpeciesAdj after convection', &
!                 MAXVAL(State_Chm_Adj%SpeciesAdj)
            IF( LPBL_ADJ ) THEN
              
              CALL DO_PBL_MIX_ADJ( LTURB, State_Chm_Adj, State_Grid_Adj, State_Met_Adj )

            ENDIF
!         print*, 'max SpeciesAdj after pbl', &
!                 MAXVAL(State_Chm_Adj%SpeciesAdj)
         ENDIF

         IF (Its_Time_For_Emis_Adj() ) THEN

            IF( LADJ_EMS ) THEN
               CALL EMISSCO2_ADJ(State_Chm_Adj, State_Grid_Adj, State_Met_Adj )

            ENDIF
  
         ENDIF


        IF ( Its_Time_For_Dyn() ) THEN
           IF( LTRAN_ADJ ) THEN
!             print*, 'transport adjoint'
             ! scale adjoint values for continuous adjoint 
             CALL Scale_Adjoint_Value(1, State_Chm_Adj, State_Met_Adj)  
             CALL Do_Transport_Adj(State_Chm_Adj, State_Grid_Adj, &
                                  State_Met_Adj, Input_Opt_Adj)
             ! rescale adjoint value
             CALL Scale_Adjoint_Value(0, State_Chm_Adj, State_Met_Adj)

            ENDIF
        ENDIF
        ADJ_FIRST = .FALSE.
!               print*, 'max SpeciesAdj after transport', &
!                      MAXVAL(State_Chm_Adj%SpeciesAdj)
      END DO

      !cost cost function in Jo
      IF ( .NOT. LICS_INC ) THEN
!      print*, 'max SpeciesAdj OBS', &
!              MAXVAL(State_Chm_Adj%SpeciesAdj) 
      COST_FUNC_JO = COST_FUNC
    
#if  defined ( LOG_OPT )
      print*, 'huoxiao_debug_geos_chem_adj log success'
      CALL LOG_RESCALE_ADJOINT
#endif

      print*, 'Calc background'
      IF ( LADJ_EMS ) THEN
      CALL CALC_APRIORI_CO2( State_Chm_Adj, State_Met_Adj, State_Grid_Adj )
      ELSEIF ( LICS ) THEN
      CALL CALC_NONDIAG_APRIORI_CO2( State_Chm_Adj, State_Met_Adj, State_Grid_Adj )
      ENDIF

      ! finite difference test for adjoint  output
      IF( LCONV_FINITE_DIFF .OR. LTRAN_FINITE_DIFF ) THEN

        CALL Write_Adjoint( State_Chm_Adj, Index_Array, Obs_Tot_Num )

      END IF
      ! finite difference test for forward cost function output

      IF( LForward_Run ) THEN 
        CALL Write_Cost

      END IF
      
      ENDIF !LICS_INC      
      ! clean up
!      CALL Exit_Tpcore
!      IF ( Associated( Index_Array        ) ) Deallocate( Index_Array  )
      Index_Array => Null()
!      CLOSE(2020)
      END SUBROUTINE DO_GEOS_CHEM_ADJ

      SUBROUTINE LOAD_CHECKPT_DATA( NYMD, NHMS, State_Chm_Adj, State_Grid_Adj, State_Met_Adj )
      
       ! USE 
       USE State_Chm_Adj_Mod,         Only : ChmStateAdj
       USE State_Grid_Adj_Mod,        Only : GrdStateAdj
       USE State_Met_Adj_Mod,         Only : MetStateAdj
       USE Checkpt_Mod,               Only : READ_CHECKPT_FILE
       INTEGER, INTENT(IN)                 :: NYMD, NHMS
       TYPE(GrdStateAdj), INTENT(IN)       :: State_Grid_Adj
       TYPE(ChmStateAdj), INTENT(INOUT)       :: State_Chm_Adj
       TYPE(MetStateAdj), INTENT(INOUT)    :: State_Met_Adj

!******************************************************************************
!  Subroutine LOAD_CHECKPT_DATA reads in information stored during the forward
!   calculation.  Some of the data (CSPEC) needs to be rotated.
!******************************************************************************
               ! Read data from file

       CALL READ_CHECKPT_FILE ( NYMD, NHMS, State_Chm_Adj, State_Grid_Adj, State_Met_Adj )

      END SUBROUTINE LOAD_CHECKPT_DATA

!-----------------------------------------------------------------------------
      SUBROUTINE CALC_ADJ_FORCE_FOR_OBS( State_Chm_Adj, State_Grid_Adj, &
                                         State_Met_Adj )

      USE OCO2_CO2_MOD,              Only : CALC_OCO2_CO2_FORCE
      USE TANSAT_CO2_MOD,            Only : CALC_TANSAT_CO2_FORCE
      USE ADJ_ARRAYS_MOD,            Only : COST_FUNC

      USE State_Chm_Adj_Mod,         Only : ChmStateAdj
      USE State_Grid_Adj_Mod,        Only : GrdStateAdj
      USE State_Met_Adj_Mod,         Only : MetStateAdj

      TYPE(GrdStateAdj), INTENT(IN)       :: State_Grid_Adj
      TYPE(ChmStateAdj), INTENT(INOUT)    :: State_Chm_Adj
      TYPE(MetStateAdj), INTENT(IN)       :: State_Met_Adj

 !     REAL*8              :: CF_GOSCO2
    !  REAL*8              :: CF_PRIOR
!
       ! Track cost function contributions
  !     CF_PRIOR = COST_FUNC
       WRITE(*,*) 'CALC_ADJ_FORCE_FOR_OBS'
       CALL CALC_TANSAT_CO2_FORCE( COST_FUNC, State_Chm_Adj, State_Grid_Adj, &
                                 State_Met_Adj )

       ! Track cost function contributions
    !   CF_GOSCO2 = CF_GOSCO2 + COST_FUNC - CF_PRIOR

       END SUBROUTINE CALC_ADJ_FORCE_FOR_OBS

      SUBROUTINE CALC_ADJ_FORCE_FOR_SENS(  )

                 CONTINUE

       END SUBROUTINE CALC_ADJ_FORCE_FOR_SENS

!-----------------------------------------------------------------------------

!      SUBROUTINE GET_SPECIES0( State_Chm_Adj )


!      USE State_Chm_Adj_Mod,         Only : ChmStateAdj

      
!      TYPE(ChmStateAdj), INTENT(INOUT)    :: State_Chm_Adj

!      State_Chm_Adj%Species0 = State_Chm_Adj%Species

!      END SUBROUTINE GET_SPECIES0

!-----------------------------------------------------------------------------

      SUBROUTINE WRITE_ADJOINT( State_Chm_Adj, Index_Array, Obs_Tot_Num )

      USE State_Chm_Adj_Mod,         Only : ChmStateAdj
      USE Adj_Arrays_Mod,            Only : COST_FUNC

      TYPE(ChmStateAdj), INTENT(IN)         :: State_Chm_Adj
      INTEGER, INTENT(IN)                 :: Index_Array(:,:)
      INTEGER                             :: Obs_Tot_Num
      INTEGER                             :: Index_I, Index_J, I, J, L

      OPEN( FILE = 'FULL_ADJOINT', UNIT = 1995 )
      OPEN( FILE = 'FULL_COST', UNIT = 1996 )
      DO I = 1, Obs_Tot_Num

         Index_I = Index_Array(1, I)
         Index_J = Index_Array(2, I)
 
         WRITE(1995, *) State_Chm_Adj%SpeciesAdj(Index_I, Index_J, 1, 1)

      END DO
      WRITE(1996, *) COST_FUNC
  
      CLOSE(1996)

      CLOSE(1995)
     
!      OPEN(UNIT=2000, FILE='huoxiao')
!      WRITE(2000,*) State_Chm_Adj%SpeciesAdj
!      CLOSE(2000) 
      END SUBROUTINE WRITE_ADJOINT

!-----------------------------------------------------------------------------

      SUBROUTINE Read_Obs_Num( Obs_Tot_Num )

      INTEGER, INTENT(INOUT)                   :: Obs_Tot_Num


      OPEN( UNIT=1997, FILE='Obs_Tot_Num' )

      READ(1997, *)Obs_Tot_Num

      CLOSE(1997)

      END SUBROUTINE Read_Obs_Num

!-----------------------------------------------------------------------------

      SUBROUTINE Init_Index_Array( Index_Array, Obs_Tot_Num )

      INTEGER, INTENT(INOUT), Pointer     :: Index_Array(:, :)
      INTEGER, INTENT(IN)                 :: Obs_Tot_Num

      INTEGER                            :: I

      ALLOCATE( Index_Array(2, Obs_Tot_Num) )
      OPEN( UNIT=1995, FILE='index_i' )
      OPEN( UNIT=1996, FILE='index_j' )

      DO I = 1, Obs_Tot_Num

         READ(1995, *) Index_Array(1, I)
         READ(1996, *) Index_Array(2, I)

      END DO

      CLOSE(1995)
      CLOSE(1996)

      END SUBROUTINE Init_Index_Array

!-----------------------------------------------------------------------------

      SUBROUTINE Write_Cost

      USE Adj_Arrays_Mod,            Only : Cost_Func

      OPEN( UNIT=1995, FILE='FORWARD_COST_FULL',&
            ACCESS='APPEND' )
      WRITE(1995, *) Cost_Func

      END  SUBROUTINE Write_Cost


      SUBROUTINE CALC_APRIORI_CO2( State_Chm_Adj, State_Met_Adj, State_Grid_Adj )

!******************************************************************************
!  Subroutine CALC_APRIORI_CO2 computes a priori term of the cost function and
!  gradient for the CO2 simulation. (dkh, 01/09/11)
!
!  In this routine, we assume that we have specified the standard deviation
!  (error) in the EMS_ERROR array.
!
!  For linear scaling factors, EMS_ERROR = p, where p is the pertent
!   error in the emissions (as a decimal, ie 1 = 100%)
!
!  For log scaling factors, EMS_ERROR = f, where f is a fractional error. f
!   must be greater than 1.
!
!  There is also a regularization parameter that is specified for each
!  emissions inventory, REG_PARAM_EMS, in input.gcadj.
!
!  NOTES:
!  ( 1) Based on CALC_APRIORI
!
!******************************************************************************

      USE State_Chm_Adj_Mod,         Only : EMS_SF
      USE State_Chm_Adj_Mod,         Only : EMS_SF0
      USE State_Chm_Adj_Mod,         Only : EMS_SF_ADJ

      USE ADJ_ARRAYS_MOD,            ONLY : COST_FUNC
      USE INVERSE_MOD,               Only : PRESSURE_WEIGHTING_FUNCTION
      USE State_Chm_Adj_Mod,         Only : ChmStateAdj
      USE State_Met_Adj_Mod,         Only : MetStateAdj
      USE State_Grid_Adj_Mod,        Only : GrdStateAdj
      USE LOGICAL_ADJ_MOD,           Only : LADJ_EMS, LICS
      USE Adj_Arrays_Mod,            Only : Obs_Exist
      USE COVARIANCE_MOD,            Only : CALC_COV_ERROR


      TYPE(GrdStateAdj), INTENT(IN)       :: State_Grid_Adj
      TYPE(ChmStateAdj), INTENT(INOUT)    :: State_Chm_Adj
      TYPE(MetStateAdj), INTENT(INOUT)    :: State_Met_Adj

      REAL*8                   :: S2_INV
      REAL*8                   :: REG
      REAL(fp)                 :: TEMP2(State_Grid_Adj%NX,State_Grid_Adj%NY &
                                  ,State_Chm_Adj%nSpecies)
      REAL(fp)                 :: FORCE, REG_ADJ, COLUMN_CONCENTRATION_ADJ
      REAL(fp)                 :: COLUMN_CONCENTRATION, COLUMN_CONCENTRATION0
      REAL(fp)                 :: S2_INV_2D(State_Grid_Adj%NX,State_Grid_Adj%NY)
      INTEGER                  :: I, J, L, N
      INTEGER                  :: IM, JM, LM, N_TRACERS
      REAL(fp)                 :: H(State_Grid_Adj%NZ)
      REAL(fp)                 :: PRIOR_COST
      REAL(fp)                 :: PEDGE_PRESSURE(State_Grid_Adj%NZ+1)
      REAL(fp)                 :: GC_BAC_ADJ(State_Grid_Adj%NZ)
      REAL*8,  ALLOCATABLE     :: APCOST(:,:,:,:)
      REAL(fp),  PARAMETER           :: TCVV_CO2 = 28.97d0  / 44d0
      INTEGER                  :: RC
      !=================================================================
      ! CALC_APRIORI_CO2 begins here!
      !=================================================================

!      ! For the moment, hardcode the emissions errors here.  In the
!      ! future, we should define these via input files.
!#if defined ( LOG_OPT )
!         ! assume a factor of two error
!         EMS_ERROR(:) = 2d0
!#else
!         ! assume a 100% error
!         EMS_ERROR(:) = 1d0
!
!         ! Alter a few to test if it's working
!         !EMS_ERROR(IDADJ_ECO2ff)  = 1d-2
!         !EMS_ERROR(IDADJ_ECO2ocn) = 1d2
!
!#endif

      IM = State_Grid_Adj%NX
      JM = State_Grid_Adj%NY

      LM = State_Grid_Adj%NZ
      N_TRACERS = State_Chm_Adj%nSpecies
      PRIOR_COST = 0.d0

      ALLOCATE(APCOST(State_Grid_Adj%NX, State_Grid_Adj%NY, 1, 1), STAT=RC)

      IF ( LADJ_EMS ) THEN
         
         CALL CALC_COV_ERROR( APCOST, State_Grid_Adj )
         PRIOR_COST = SUM(APCOST)
      ENDIF

      IF( LICS ) THEN 
      DO N = 1, N_TRACERS
      DO J = 1, JM
      DO I = 1, IM


         S2_INV = 1d0 / ( State_Chm_Adj%EMS_ERROR(N)** 2 )

         REG = EMS_SF(I,J,N) - EMS_SF0(I,J,N)

         ! Calculate the contribution to the cost function, weighted by REG_PARAM
         TEMP2(I,J,N) = 0.5d0 * State_Chm_Adj%REG_PARAM_EMS(N) * & 
                                  S2_INV*REG**2

         ! Add this to the gradients
         EMS_SF_ADJ(I,J,N) = EMS_SF_ADJ(I,J,N)   &
                             + State_Chm_Adj%REG_PARAM_EMS(N) * S2_INV * REG

      ENDDO
      ENDDO

      ! Add total regularization penalty to cost function

      ENDDO
      print*, 'huoxiao_debug max temp2', MAXVAL(TEMP2), MINVAL(TEMP2),&
                                         MAXVAL(EMS_SF), MAXVAL(EMS_SF0), &
                                         MINVAL(EMS_SF), MINVAL(EMS_SF0)
      PRIOR_COST = SUM(TEMP2)
      ENDIF


!      IF( LADJ_EMS ) THEN
!      DO N = 1, N_TRACERS
!      DO J = 1, JM
!      DO I = 1, IM

!         COLUMN_CONCENTRATION = 0.0
!         COLUMN_CONCENTRATION0 = 0.0

!         DO L = 1, LM+1
!            PEDGE_PRESSURE(L) = State_Met_Adj%AP(L)+&
!                 State_Met_Adj%PSC2WET(I, J)*State_Met_Adj%BP(L)
!         ENDDO
!         H = PRESSURE_WEIGHTING_FUNCTION(PEDGE_PRESSURE, &
!                                State_Met_Adj%PSC2WET(I, J), &
!                                State_Met_Adj%SPECIFIC_HUMIDITY(I, J, :), &
!                                LM)
!         DO L = 1, LM
!            COLUMN_CONCENTRATION = COLUMN_CONCENTRATION+H(L)* &
!                                  State_Chm_Adj%Species(I, J, L, 1)
            
!            COLUMN_CONCENTRATION0 = COLUMN_CONCENTRATION0+H(L)*       &
!                                State_Chm_Adj%Species0(I, J, L, 1)
!         ENDDO


!         S2_INV = 1d0 /  State_Chm_Adj%ICS_ERROR(N)

!         REG = COLUMN_CONCENTRATION  - COLUMN_CONCENTRATION0

         ! Calculate the contribution to the cost function, weighted by REG_PARAM
!         FORCE = State_Chm_Adj%REG_PARAM_ICS(N) * &
!                                  S2_INV*REG
!         REG_ADJ = FORCE
!         COLUMN_CONCENTRATION_ADJ = REG_ADJ
         ! Add this to the gradients
!         DO L = 1, LM
!            GC_BAC_ADJ(L) = H(L)*COLUMN_CONCENTRATION_ADJ
!            State_Chm_Adj%SpeciesAdj(I,J,L,N) = State_Chm_Adj%SpeciesAdj(I,J,L,N)   &
!                                  +GC_BAC_ADJ(L)*TCVV_CO2
!         ENDDO

!         PRIOR_COST = PRIOR_COST+0.5d0*FORCE*REG  
!      ENDDO
!      ENDDO
!      ENDDO
!      ENDIF


      IF( LICS ) THEN
      
  
      DO N = 1, N_TRACERS
      DO L = 1, LM
      DO J = 1, JM
      DO I = 1, IM

        S2_INV = 1d0 / State_Chm_Adj%ICS_ERROR(N) 
        REG = State_Chm_Adj%Species(I, J, L, N)-State_Chm_Adj%Species0(I, J, L, N)
        !different regularization parameter
        IF( OBS_EXIST(I, J) == 1) THEN
          State_Chm_Adj%REG_PARAM_ICS(N) = 0.0008
        ELSE
          State_Chm_Adj%REG_PARAM_ICS(N) = 0.000008
        ENDIF
        FORCE = S2_INV*REG*State_Chm_Adj%REG_PARAM_ICS(N)
        State_Chm_Adj%SpeciesAdj(I,J,L,N) = State_Chm_Adj%SpeciesAdj(I,J,L,N)   &
                                          + FORCE*TCVV_CO2
        PRIOR_COST = PRIOR_COST+0.5d0*REG*FORCE
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDIF

         
      PRINT*,' huoxiao_debug COST_FUNC ', COST_FUNC
      COST_FUNC = COST_FUNC+PRIOR_COST
      PRINT*, 'huoxiao_debug APRIOR COST FUNCTION:', PRIOR_COST
      END SUBROUTINE CALC_APRIORI_CO2
!------------------------------------------------------------------------------
      SUBROUTINE CALC_NONDIAG_APRIORI_CO2( State_Chm_Adj, State_Met_Adj, State_Grid_Adj )
      USE State_Chm_Adj_Mod,         Only : ChmStateAdj
      USE State_Met_Adj_Mod,         Only : MetStateAdj
      USE State_Grid_Adj_Mod,        Only : GrdStateAdj
      USE BEC_MOD,                   ONLY : BinvW
      USE ADJ_ARRAYS_MOD,            ONLY : COST_FUNC
      USE ADJ_ARRAYS_MOD,            ONLY : OBS_EXIST
      USE FILE_OPE_MOD,              ONLY : READ_REG_PARAM

      TYPE(GrdStateAdj), INTENT(IN)       :: State_Grid_Adj
      TYPE(ChmStateAdj), INTENT(INOUT)    :: State_Chm_Adj
      TYPE(MetStateAdj), INTENT(INOUT)    :: State_Met_Adj
      REAL(fp)                            :: BinvX(State_Grid_Adj%NX, &
                                                   State_Grid_Adj%NY, &
                                                   State_Grid_Adj%NZ, &
                                                   State_Chm_Adj%nSpecies)
      REAL(fp)                            :: DeltaX(State_Grid_Adj%NX, &
                                                   State_Grid_Adj%NY, &
                                                   State_Grid_Adj%NZ, &
                                                   State_Chm_Adj%nSpecies)
      REAL(fp)                            :: TMP_BINVX(State_Grid_Adj%NX, &
                                                   State_Grid_Adj%NY, &
                                                   State_Grid_Adj%NZ)

      REAL(fp)                            :: PRIOR_COST
      INTEGER                             :: I, J, L, N
      INTEGER                             :: NX, NY, NZ, N_TRACERS
      REAL(fp)                            :: TMP_OBS, TMP_NONOBS
      
      LOGICAL                             :: ALIVE, DFLAG
      CHARACTER(LEN=20)                   :: FILE_N
      INTEGER                             :: TMP_FLAG, FN
      REAL(fp)                            :: REG_PARAM, FORCE
      REAL(fp)                            :: TCVV_CO2

      TCVV_CO2 = 44/28.97

      NX = State_Grid_Adj%NX
      NY = State_Grid_Adj%NY
      NZ = State_Grid_Adj%NZ
      N_TRACERS = State_Chm_Adj%nSpecies
      BinvX = 0.0
      PRIOR_COST = 0.0
      TMP_OBS = 0.0
      TMP_NONOBS = 0.0
      print*, 'background batch'
      DO WHILE(.TRUE.)
         DFLAG = .TRUE.
         DO FN = 1, 5
            WRITE(FILE_n, '(I2.2)') FN
            INQUIRE( FILE = './BEC_wn13/FLAG_'//FILE_N, EXIST = ALIVE )  
!            CALL SLEEP(0.5) 

            IF ( .NOT. ALIVE ) THEN
               DFLAG = .FALSE.
               EXIT
!            ELSE
!               OPEN(UNIT=1995, FILE='FLAG_'//FILE_N)
!               READ(1995, *)TMP_FLAG
!               CLOSE(1995)
!               IF ( TMP_FLAG /= 1 ) THEN
!                   DFLAG = .FALSE.
!                   EXIT
!               ENDIF
            ENDIF
         ENDDO
         IF ( DFLAG ) EXIT
       ENDDO

       DO FN = 1, 5
          WRITE(FILE_n, '(I2.2)') FN 
          OPEN(UNIT=1995, FILE='BEC_wn13/BINVX_'//FILE_N)
          READ(1995, *) TMP_BINVX
          BINVX(:, :, :, 1) = BINVX(:, :, :, 1)+TMP_BINVX(:, :, :)
          print*, 'BLOCK MAX', MAXVAL(TMP_BINVX(:, :, :))
       ENDDO

      DeltaX = (State_Chm_Adj%Species-State_Chm_Adj%Species0)*TCVV_CO2
!      IF (N_CALC==2) THEN
!      print*, DeltaX 
!      ENDIF
      State_Chm_Adj%D_inv = 7.56*10000000000000*1/(TCVV_CO2**2)
      CALL READ_REG_PARAM(REG_PARAM)
!      DO N = 1, N_TRACERS
!      DO L = 1, NZ
!      DO J = 1, NY
!      DO I = 1, NX
!         IF ( OBS_EXIST(I, J) == 1 ) THEN
!            State_Chm_Adj%D_inv(I, J, L, 1) =1
!            TMP_OBS = TMP_OBS+DeltaX(I, J, L,1)
!         ELSE 
!            State_Chm_Adj%D_inv(I, J, L, 1) = 1
!            TMP_NONOBS = TMP_NONOBS+DeltaX(I, J, L,1)
!         ENDIF
!      ENDDO
!      ENDDO
!      ENDDO
!      ENDDO
!      PRINT*, 'TMP_OBS ', TMP_OBS, TMP_NONOBS
      !
!      BinvX(:, :, :, 1) = BinvW( State_Chm_Adj%Sx_inv(:, :), &
!                                 State_Chm_Adj%Sy_inv(:, :), &
!                                 State_Chm_Adj%Sz_inv(:, :), &
!                                 State_Chm_Adj%D_inv(:, :, :, 1), &
!                                 DeltaX         (:, :, :, 1), &
!                                 State_Grid_Adj%NX,      &
!                                 State_Grid_Adj%NY,      &
!                                 State_Grid_Adj%NZ)
!      DO N = 1, N_TRACERS
!      DO L = 1, NZ
!      DO J = 1, NY
!      DO I = 1, NX
!         IF ( OBS_EXIST(I, J) == 1 ) THEN
!            print*, 'DIFF', DeltaX(I, J, L, 1), BinvX(I, J, L, 1)
!         ELSE
!            State_Chm_Adj%D_inv(I, J, L, 1) = 100000000
!         ENDIF
!      ENDDO
!      ENDDO
!      ENDDO
!      ENDDO

!      print*, 'MAXVAL D_inv', MAXVAL(State_Chm_Adj%D_inv)
      print*, 'DeltaX', MAXVAL(DeltaX)
      print*, 'MAXVAL BinvX', MAXVAL(BinvX)
      !
      DO N = 1, N_TRACERS
      DO L = 1, NZ
      DO J = 1, NY
      DO I = 1, NX
         FORCE = BinvX(I, J, L, N)*REG_PARAM*&
                 State_Chm_Adj%D_inv(I, J, L, N)
         PRIOR_COST = PRIOR_COST+0.5*FORCE*DeltaX(I, J, L, N)
         State_Chm_Adj%SpeciesAdj(I, J, L, N) = &
         State_Chm_Adj%SpeciesAdj(I, J, L, N) + FORCE
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      print*, 'PRIOR_COST', PRIOR_COST
      print*, 'COST_FUNC', COST_FUNC
      COST_FUNC = PRIOR_COST+COST_FUNC
      END SUBROUTINE CALC_NONDIAG_APRIORI_CO2

!------------------------------------------------------------------------------

      SUBROUTINE LOG_RESCALE_ADJOINT

      USE State_Chm_Adj_Mod,         Only : EMS_SF, EMS_SF0, EMS_SF_ADJ
      USE LOGICAL_ADJ_MOD,           ONLY : LADJ_EMS
       
      IF ( LADJ_EMS )  THEN

         ! Transform back to exponential scaling factors
         !EMS_SF_ADJ(:,:,:) = EMS_SF_ADJ(:,:,:) * EMS_SF(:,:,:) use scaled flux, remove here
         EMS_SF(:,:,:)     = LOG(EMS_SF(:,:,:))
         EMS_SF0(:,:,:)    = LOG(EMS_SF0(:,:,:))

      ENDIF
     
      END SUBROUTINE

      
      !module end
      END MODULE GEOS_CHEM_ADJ_MOD

      Module LICS_INC_Mod

      USE PRECISION_MOD

      CONTAINS
      
      SUBROUTINE DO_LICS_INC         ( State_Chm_Adj, State_Grid_Adj, &
                                       State_Met_Adj, Input_Opt_Adj )
      
      USE State_Chm_Adj_Mod,         Only : ChmStateAdj
      USE State_Grid_Adj_Mod,        Only : GrdStateAdj
      USE State_Met_Adj_Mod,         Only : MetStateAdj
      USE Input_Adj_Mod,             Only : OptInputAdj

      USE BEC_Mod,                   Only : Uw
      USE BEC_Mod,                   Only : UTw
      USE Logical_Adj_Mod,           Only : LICS
      USE Logical_Adj_Mod,           Only : LICS_INC
      USE Geos_Chem_Adj_Mod,         Only : Do_Geos_Chem_Adj
      USE Geos_Chem_D_Mod,           Only : Do_Geos_Chem_D
      USE GEOS_CHEM_BALANCE_MOD,     Only : CALCULATE_JC
      USE Tansat_Co2_Mod,            Only : Calc_Observation_Operator_D
      USE Tansat_Co2_Mod,            Only : OBS_EXIT
      USE Time_Mod,                  Only : Get_Ts_Dyn
      USE Time_Mod,                  Only : Get_Elapsed_Sec
      USE Time_Mod,                  Only : Set_Elapsed_Sec_Adj
      USE Time_Mod,                  Only : Set_Current_Time
      USE TIme_Mod,                  Only : Get_Elapsed_Sec
      USE TIme_Mod,                  Only : GET_NYMD
      USE TIme_Mod,                  Only : GET_NHMS
      USE File_Ope_Mod,              Only : Get_BEC
      USE File_Ope_Mod,              Only : GET_SPECIES0
      USE Checkpt_Mod,               Only : Load_Checkpt_Data
      USE Adj_Arrays_Mod,            Only : COST_FUNC_JO, COST_FUNC 
      USE Adj_Arrays_Mod,            Only : OBS_EXIST
      TYPE(ChmStateAdj), INTENT(INOUT)    :: State_Chm_Adj
      TYPE(GrdStateAdj), INTENT(INOUT)    :: State_Grid_Adj
      TYPE(MetStateAdj), INTENT(INOUT)    :: State_Met_Adj
      TYPE(OptInputAdj), INTENT(INOUT)    :: Input_Opt_Adj

      REAL*8                              :: TMP_SpeciesAdj(State_Grid_Adj%NX, &
                                                            State_Grid_Adj%NY, &
                                                            State_Grid_Adj%NZ, &
                                                            1)
      INTEGER                             :: Obs_Period, Ts_Dyn 
      INTEGER                             :: Final_Elapsed_Min, Cur_Elapsed_Min
      INTEGER                             :: Obs_ith
  
      !initialization
      TMP_SpeciesAdj = 0.0
      !Ts_Dyn = Get_Ts_Dyn()
      Final_Elapsed_Min = Get_Elapsed_Sec()
      Obs_Period = INT(Final_Elapsed_Min/3600)-1
      OBS_EXIST = 0.0
!      print*, 'stage 1'      
!      IF ( LICS .AND. (.NOT. LICS_INC) ) THEN
!         CALL Do_Geos_Chem_Adj( State_Chm_Adj,   &
!                                State_Grid_Adj,  &
!                                State_Met_Adj,   &
!                                Input_Opt_Adj )
!         PRINT*, 'HUOXIAO_DEBUG 4DVar'
!      ELSE
         ! read bec
!         print*, 'stage 2'

!      CALL Get_BEC( State_Chm_Adj, State_Grid_Adj )     
!      CALL GET_SPECIES0( State_Chm_Adj, State_Grid_Adj )
      DO Obs_ith = 0, Obs_Period
         print*, 'Obs_ith is ', Obs_ith
         !set time
         Cur_Elapsed_Min = INT(Obs_ith*3600)
         CALL Set_Elapsed_Sec_Adj( 0, .TRUE., Cur_Elapsed_Min, .TRUE.)
         CALL Set_Current_Time
!         print*, 'huoxiao_debug Obs_ith', Obs_ith
         print*, 'NYMD ', GET_NYMD(), ' NHMS ', GET_NHMS()

         ! observation exit?
         IF( .NOT. OBS_EXIT() ) THEN
             CYCLE
          ENDIF

         ! read M_{0-t}(xb)
         CALL LOAD_CHECKPT_DATA( State_Chm_Adj, State_Grid_Adj, State_Met_Adj )

!         ! calculate Uw
         State_Chm_Adj%SpeciesAdj(:, :, :, 1) = Uw( State_Chm_Adj%Sx(:, :), &
                                                    State_Chm_Adj%Sy(:, :), &
                                                    State_Chm_Adj%Sz(:, :), &
                                                    State_Chm_Adj%D(:, :, :, 1), &
                                                    State_Chm_Adj%CV(:, :, :, 1), &
                                                    State_Grid_Adj%NX,      &
                                                    State_Grid_Adj%NY,      &
                                                    State_Grid_Adj%NZ)

!         IF ( Obs_ith == 0 ) THEN
!            print*, 'huoxiao_debug_max Uw, min Uw, '//&
!                    'max CV, min CV, max SpeciestB,'//&
!                    'min SpeciestB ',         & 
!            MAXVAL(State_Chm_Adj%SpeciesAdj), &
!            MINVAL(State_Chm_Adj%SpeciesAdj), &
!            MAXVAL(State_Chm_Adj%CV),         &
!            MINVAL(State_Chm_Adj%CV),         &
!            MAXVAL(State_Chm_Adj%SpeciestB),  &
!            MINVAL(State_Chm_Adj%SpeciestB)
            !calculate xa = xb+UW
!            State_Chm_Adj%Species0 = State_Chm_Adj%SpeciestB+&
!                                     State_Chm_Adj%SpeciesAdj
!            print*, 'huoxiao_debug_species0 ', MAXVAL(State_Chm_Adj%Species0), &
!                    MINVAL(State_Chm_Adj%Species0)
!         ENDIF
!         print*, 'max SpeciesAdj U', &
!                  MAXVAL(State_Chm_Adj%SpeciesAdj)

         ! calculate MUw
!         print*, 'stage 4'
         IF ( Obs_ith /= 0 ) THEN
         CALL Do_Geos_Chem_D( State_Chm_Adj,   &
                              State_Grid_Adj,  &
                              State_Met_Adj,   &
                              Input_Opt_Adj )
         ENDIF
         ! calculate H^T*O^{-1}*(HMUw+di)
         CALL Calc_Observation_Operator_D( State_Chm_Adj,  &
                                           State_Grid_Adj, &
                                           State_Met_Adj )
         IF ( Obs_ith /= 0 ) THEN
         ! calculate M^T*H^T*O^{-1}*(HMUw+di)
         CALL Do_Geos_Chem_Adj(State_Chm_Adj,   &
                               State_Grid_Adj,  &
                               State_Met_Adj,   &
                               Input_Opt_Adj )
         ENDIF
!         print*, 'max SpeciesAdj after adjoint', &
!                  MAXVAL(State_Chm_Adj%SpeciesAdj)
         ! calculate U^T*M^T*H^T*O^{-1}*(HMUw+di)
!         State_Chm_Adj%SpeciesAdj = 0.0
!         State_Chm_Adj%SpeciesAdj(20, 20, 1, 1) = 1.0
         State_Chm_Adj%SpeciesAdj(:, :, :, 1) = UTw( Transpose(State_Chm_Adj%Sx(:, :)), &
                                                    Transpose(State_Chm_Adj%Sy(:, :)), &
                                                    Transpose(State_Chm_Adj%Sz(:, :)), &
                                                    State_Chm_Adj%D(:, :, :, 1), &
                                                    State_Chm_Adj%SpeciesAdj(:, :, :, 1), &
                                                    State_Grid_Adj%NX,      &
                                                    State_Grid_Adj%NY,      &
                                                    State_Grid_Adj%NZ)
!         print*, 'max SpeciesAdj U^T', &
!                  MAXVAL(State_Chm_Adj%SpeciesAdj)
         TMP_SpeciesAdj = State_Chm_Adj%SpeciesAdj+TMP_SpeciesAdj
      ENDDO
      State_Chm_Adj%SpeciesAdj = TMP_SpeciesAdj
      COST_FUNC_JO = COST_FUNC
      CALL  Calc_APRIOR_CV( State_Chm_Adj,  &
                            State_Grid_Adj )
!      CALL CALCULATE_JC( State_Chm_Adj, State_Grid_Adj, &
!                         State_Met_Adj, Input_Opt_Adj )

!      CALL WRITE_SF_TO_RSTFILE( State_Chm_Adj, State_Grid_Adj )      
!      ENDIF
     
      END SUBROUTINE DO_LICS_INC
!------------------------------------------------------------------------------
      SUBROUTINE Calc_APRIOR_CV( State_Chm_Adj, State_Grid_Adj )

      USE State_Chm_Adj_Mod,         Only : ChmStateAdj
      USE State_Grid_Adj_Mod,        Only : GrdStateAdj
      USE Logical_Adj_Mod,           Only : LICS_INC
      USE Adj_Arrays_Mod,            Only : COST_FUNC
      USE Adj_Arrays_Mod,            Only : OUT_ITER

      TYPE(GrdStateAdj), INTENT(IN)       :: State_Grid_Adj
      TYPE(ChmStateAdj), INTENT(INOUT)    :: State_Chm_Adj

      REAL(fp)                 :: PRIOR_COST
      REAL(fp)                 :: TMP_CV(State_Grid_Adj%NX, &
                                         State_Grid_Adj%NY, &
                                         State_Grid_Adj%NZ, &
                                         State_Chm_Adj%nSpecies)
      REAL(fp)                 :: SUM_CV(State_Grid_Adj%NX, &
                                         State_Grid_Adj%NY, &
                                         State_Grid_Adj%NZ, &
                                         State_Chm_Adj%nSpecies)
      INTEGER                  :: I, J, L, N, ITE
      INTEGER                  :: IM, JM, LM, N_TRACERS

      IM = State_Grid_Adj%NX
      JM = State_Grid_Adj%NY
      LM = State_Grid_Adj%NZ
      N_TRACERS = State_Chm_Adj%nSpecies
      TMP_CV = 0.0
      SUM_CV = State_Chm_Adj%CV
      PRIOR_COST = 0.0

      IF ( LICS_INC ) THEN 
      DO ITE = 1, OUT_ITER-1
      print*, 'huoxiao_debug last_cv ', ITE, OUT_ITER
      CALL READ_LAST_CV(TMP_CV, ITE, IM, JM, LM, N_TRACERS)
      SUM_CV = SUM_CV+TMP_CV
      ENDDO

      DO N = 1, N_TRACERS
      DO L = 1, LM
      DO J = 1, JM
      DO I = 1, IM
         PRIOR_COST = PRIOR_COST+0.5*SUM_CV(I, J, L ,N)**2
         State_Chm_Adj%SpeciesAdj(I, J, L, N) = &
            State_Chm_Adj%SpeciesAdj(I, J, L, N)+SUM_CV(I, J, L ,N)
      ENDDO
      ENDDO
      ENDDO
      ENDDO 
      ENDIF
      COST_FUNC = COST_FUNC+PRIOR_COST
      print*, 'JB_COST:', PRIOR_COST
      print*, 'TOTAL COST_FUNC', COST_FUNC
      END SUBROUTINE Calc_APRIOR_CV
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE READ_LAST_CV( CV, ITE, NX, NY, NZ, N_TRACERS)
      INTEGER, INTENT(IN)     :: ITE, NX, NY, NZ, N_TRACERS
      REAL*8, INTENT(INOUT)   :: CV(NX, NY, NZ, N_TRACERS)
      CHARACTER(LEN=20)       :: CHA_ITE_NUM

      WRITE(CHA_ITE_NUM, '(I2.2)') ITE
      print*, 'READ FILE: '//'WRF4DVAR/gctm.w.'//CHA_ITE_NUM 
      OPEN(UNIT=1995, FILE='WRF4DVAR/gctm.w.'//CHA_ITE_NUM)
      READ(1995, *) CV
      CLOSE(1995)
      
      END SUBROUTINE READ_LAST_CV
      END Module LICS_INC_Mod

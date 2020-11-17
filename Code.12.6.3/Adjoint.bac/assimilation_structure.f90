      Module Assimilation_Structure_Mod

      USE PRECISION_MOD

      CONTAINS
      
      SUBROUTINE Do_Assimilation_Structure( State_Chm_Adj, State_Grid_Adj, &
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
      USE Tansat_Co2_Mod,            Only : Calc_Observation_Operator_D
      USE Time_Mod,                  Only : Get_Ts_Dyn
      USE Time_Mod,                  Only : Get_Elapsed_Sec
      USE Time_Mod,                  Only : Set_Elapsed_Sec_Adj
      USE Time_Mod,                  Only : Set_Current_Time
      USE TIme_Mod,                  Only : Get_Elapsed_Sec
      USE File_Ope_Mod,              Only : Get_BEC

      USE Checkpt_Mod,               Only : Load_Checkpt_Data
  
      TYPE(ChmStateAdj), INTENT(INOUT)    :: State_Chm_Adj
      TYPE(GrdStateAdj), INTENT(INOUT)    :: State_Grid_Adj
      TYPE(MetStateAdj), INTENT(INOUT)    :: State_Met_Adj
      TYPE(OptInputAdj), INTENT(INOUT)    :: Input_Opt_Adj

      INTEGER                             :: Obs_Period, Ts_Dyn 
      INTEGER                             :: Final_Elapsed_Min, Cur_Elapsed_Min
      INTEGER                             :: Obs_ith


      !Ts_Dyn = Get_Ts_Dyn()
      Final_Elapsed_Min = Get_Elapsed_Sec()
      Obs_Period = INT(Final_Elapsed_Min/3600)
      print*, 'stage 1'      
      IF ( .NOT. LICS_INC ) THEN
         CALL Do_Geos_Chem_Adj( State_Chm_Adj,   &
                                State_Grid_Adj,  &
                                State_Met_Adj,   &
                                Input_Opt_Adj )
      ELSE

         ! read bec
         print*, 'stage 2'
         CALL Get_BEC( State_Chm_Adj, State_Grid_Adj )     
         DO Obs_ith = 1, Obs_Period
         !set time
         Cur_Elapsed_Min = INT(Obs_ith*3600)
         CALL Set_Elapsed_Sec_Adj( 0, .TRUE., Cur_Elapsed_Min, .TRUE.)
         CALL Set_Current_Time
         print*, 'stage 3', GET_ELAPSED_SEC()
         ! read M_{0-t}(xb)
         CALL LOAD_CHECKPT_DATA( State_Chm_Adj, State_Grid_Adj, State_Met_Adj )
!         ! calculate Uw
!         State_Chm_Adj%SpeciesAdj = 0.0
!         State_Chm_Adj%SpeciesAdj(20, 20, 1, 1) = 1.0
         State_Chm_Adj%SpeciesAdj(:, :, :, 1) = Uw( State_Chm_Adj%Sx(:, :), &
                                                    State_Chm_Adj%Sy(:, :), &
                                                    State_Chm_Adj%Sz(:, :), &
                                                    State_Chm_Adj%D(:, :, :, 1), &
                                                    State_Chm_Adj%CV(:, :, :, 1), &
                                                    State_Grid_Adj%NX,      &
                                                    State_Grid_Adj%NY,      &
                                                    State_Grid_Adj%NZ)
         print*, 'max SpeciesAdj U', &
                  MAXVAL(State_Chm_Adj%SpeciesAdj)
         ! calculate MUw
         print*, 'stage 4'
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
         print*, 'max SpeciesAdj after adjoint', &
                  MAXVAL(State_Chm_Adj%SpeciesAdj)
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
         print*, 'max SpeciesAdj U^T', &
                  MAXVAL(State_Chm_Adj%SpeciesAdj)

         ENDDO
         CALL  Calc_APRIOR_CV( State_Chm_Adj,  &
                               State_Grid_Adj )
      ENDIF
     
      END SUBROUTINE Do_Assimilation_Structure
!------------------------------------------------------------------------------
      SUBROUTINE Calc_APRIOR_CV( State_Chm_Adj, State_Grid_Adj )

      USE State_Chm_Adj_Mod,         Only : ChmStateAdj
      USE State_Grid_Adj_Mod,        Only : GrdStateAdj
      USE Logical_Adj_Mod,           Only : LICS_INC
      USE Adj_Arrays_Mod,            Only : COST_FUNC

      TYPE(GrdStateAdj), INTENT(IN)       :: State_Grid_Adj
      TYPE(ChmStateAdj), INTENT(INOUT)    :: State_Chm_Adj

      REAL(fp)                 :: PRIOR_COST
      INTEGER                  :: I, J, L, N
      INTEGER                  :: IM, JM, LM, N_TRACERS

      IM = State_Grid_Adj%NX
      JM = State_Grid_Adj%NY
      LM = State_Grid_Adj%NZ
      N_TRACERS = State_Chm_Adj%nSpecies

      IF ( LICS_INC ) THEN 

      DO N = 1, N_TRACERS
      DO L = 1, LM
      DO J = 1, JM
      DO I = 1, IM
         PRIOR_COST = PRIOR_COST+0.5*State_Chm_Adj%CV(I, J, L ,N)**2
      ENDDO
      ENDDO
      ENDDO
      ENDDO 
      ENDIF
      COST_FUNC = COST_FUNC+PRIOR_COST
      print*, 'PRIOR_COST:', PRIOR_COST
      print*, 'COST_FUNC', COST_FUNC
      END SUBROUTINE Calc_APRIOR_CV 
      END Module Assimilation_Structure_Mod

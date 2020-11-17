      MODULE GEOS_CHEM_BALANCE_MOD
      
      USE PRECISION_MOD
      IMPLICIT NONE
    
      CONTAINS
 
      SUBROUTINE CALCULATE_JC( State_Chm_Adj, State_Grid_Adj, &
                               State_Met_Adj, Input_Opt_Adj )

      USE State_Chm_Adj_Mod,         Only : ChmStateAdj
      USE State_Grid_Adj_Mod,        Only : GrdStateAdj
      USE State_Met_Adj_Mod,         Only : MetStateAdj
      USE Input_Adj_Mod,             Only : OptInputAdj

      USE BEC_Mod,                   Only : Uw
      USE BEC_Mod,                   Only : UTw

      USE Geos_Chem_Adj_Mod,         Only : Do_Geos_Chem_Adj
      USE Geos_Chem_D_Mod,           Only : Do_Geos_Chem_D

      USE Time_Mod,                  Only : Set_Elapsed_Sec_Adj
      USE Time_Mod,                  Only : Set_Current_Time

      TYPE(ChmStateAdj), INTENT(INOUT)    :: State_Chm_Adj
      TYPE(GrdStateAdj), INTENT(INOUT)    :: State_Grid_Adj
      TYPE(MetStateAdj), INTENT(INOUT)    :: State_Met_Adj
      TYPE(OptInputAdj), INTENT(INOUT)    :: Input_Opt_Adj


      REAL*8                              :: JC_FORCE(State_Grid_Adj%NX, &
                                                      State_Grid_Adj%NY, &
                                                      State_Grid_Adj%NZ, &
                                                      1)
      REAL*8                       :: SUMD_SPECIESADJ(State_Grid_Adj%NX, &
                                                      State_Grid_Adj%NY, &
                                                      State_Grid_Adj%NZ, &
                                                      1)

      REAL*8                       :: SUMAD_SPECIESADJ(State_Grid_Adj%NX, &
                                                      State_Grid_Adj%NY, &
                                                      State_Grid_Adj%NZ, &
                                                      1)
      REAL*8                       :: TMP_SPECIESADJ(State_Grid_Adj%NX, &
                                                      State_Grid_Adj%NY, &
                                                      State_Grid_Adj%NZ, &
                                                      1)
      INTEGER                             :: INTG_Period, INTG_ITH
      INTEGER                             :: Final_Elapsed_Min, Cur_Elapsed_Min

      Intg_Period = 6      
      JC_FORCE = 0.0
      SUMD_SPECIESADJ = 0.0
      SUMAD_SPECIESADJ = 0.0
      TMP_SPECIESADJ = State_Chm_Adj%SpeciesAdj
      CALL CALC_g( State_Chm_Adj )
 
      print*, 'This is geos_chem_balance.f90'
      DO INTG_ITH = 0, INTG_PERIOD
         Cur_Elapsed_Min = INT(INTG_ITH*600)
         CALL Set_Elapsed_Sec_Adj( 0, .TRUE., Cur_Elapsed_Min, .TRUE.)
         CALL Set_Current_Time
         State_Chm_Adj%SpeciesAdj(:, :, :, 1) = Uw( State_Chm_Adj%Sx(:, :), &
                                                    State_Chm_Adj%Sy(:, :), &
                                                    State_Chm_Adj%Sz(:, :), &
                                                    State_Chm_Adj%D(:, :, :, 1), &
                                                    State_Chm_Adj%CV(:, :, :, 1), &
                                                    State_Grid_Adj%NX,      &
                                                    State_Grid_Adj%NY,      &
                                                    State_Grid_Adj%NZ)

         IF ( INTG_ith /= 0 ) THEN
         ! calculate M^T*H^T*O^{-1}*(HMUw+di)
         CALL Do_Geos_Chem_D(State_Chm_Adj,   &
                               State_Grid_Adj,  &
                               State_Met_Adj,   &
                               Input_Opt_Adj )
         ENDIF
         JC_FORCE = JC_FORCE+State_Chm_Adj%g(INTG_ITH+1)*State_Chm_Adj%SpeciesAdj
         !CALL CALC_JC_COST( State_Chm_Adj, State_Grid_Adj )
         SUMD_SPECIESADJ = SUMD_SPECIESADJ+  &
                          State_Chm_Adj%f(INTG_ITH+1)*State_Chm_Adj%SpeciesAdj
      ENDDO
      print*, 'running jc cost'
      CALL CALC_JC_COST( JC_FORCE, State_Chm_Adj, State_Grid_Adj )
      print*, 'running cw'
      CALL CALC_Cw( SUMD_SPECIESADJ, State_Chm_Adj, State_Grid_Adj )
      print*, 'running adjoint'
      DO INTG_ITH = INTG_PERIOD, 0, -1
         SUMD_SPECIESADJ = State_Chm_Adj%gama_df(1)*&
                           State_Chm_Adj%f(INTG_ITH+1)*SUMD_SPECIESADJ    
         State_Chm_Adj%SpeciesAdj = SUMD_SPECIESADJ
         Cur_Elapsed_Min = INT(INTG_ITH*600)
         CALL Set_Elapsed_Sec_Adj( 0, .TRUE., Cur_Elapsed_Min, .TRUE.)
         CALL Set_Current_Time
         IF ( INTG_ITH /= 0 ) THEN
         ! calculate M^T*H^T*O^{-1}*(HMUw+di)
         CALL Do_Geos_Chem_Adj(State_Chm_Adj,   &
                               State_Grid_Adj,  &
                               State_Met_Adj,   &
                               Input_Opt_Adj )
         ENDIF
         SUMAD_SPECIESADJ = SUMAD_SPECIESADJ + &
                            State_Chm_Adj%SpeciesAdj
      ENDDO
      SUMAD_SPECIESADJ(:, :, :, 1) = UTw( State_Chm_Adj%Sx(:, :), &
                                                    State_Chm_Adj%Sy(:, :), &
                                                    State_Chm_Adj%Sz(:, :), &
                                                    State_Chm_Adj%D(:, :, :, 1), &
                                                    SUMAD_SPECIESADJ(:, :, :, 1), &
                                                    State_Grid_Adj%NX,      &
                                                    State_Grid_Adj%NY,      &
                                                    State_Grid_Adj%NZ) 
      State_Chm_Adj%SpeciesAdj = SUMAD_SPECIESADJ+TMP_SPECIESADJ
      print*, 'JC_GRAD', MAXVAL(SUMAD_SPECIESADJ), MINVAL(SUMAD_SPECIESADJ)
      print*,'JB+JO_GRAD', MAXVAL(TMP_SPECIESADJ), MINVAL(TMP_SPECIESADJ) 
      END SUBROUTINE  CALCULATE_JC

      SUBROUTINE CALC_JC_COST ( JC_FORCE, State_Chm_Adj, State_Grid_Adj )

      USE State_Chm_Adj_Mod,         Only : ChmStateAdj
      USE State_Grid_Adj_Mod,        Only : GrdStateAdj
      USE Adj_Arrays_Mod,            Only : COST_FUNC

      TYPE(GrdStateAdj), INTENT(IN)       :: State_Grid_Adj
      TYPE(ChmStateAdj), INTENT(INOUT)    :: State_Chm_Adj
      REAL(fp),          INTENT(IN)       :: JC_FORCE( State_Grid_Adj%NX,&
                                                       State_Grid_Adj%NY, &
                                                       State_Grid_Adj%NZ, &
                                                       State_Chm_Adj%nSpecies)
      INTEGER       :: I, J, L, N, NX, NY, NZ, N_TRACERS
      REAL(fp)      :: JC_COST_FUNC
      NX = State_Grid_Adj%NX
      NY = State_Grid_Adj%NY
      NZ = State_Grid_Adj%NZ
      N_TRACERS = State_Chm_Adj%nSpecies
      JC_COST_FUNC = 0.0 
      DO N = 1, N_TRACERS
      DO L = 1, NZ
      DO J = 1, NY
      DO I = 1, NX
         JC_COST_FUNC = JC_COST_FUNC+State_Chm_Adj%C(I, J, L, N)* &
                        State_Chm_Adj%gama_df(1)*JC_FORCE(I, J, L, N)**2
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      COST_FUNC = COST_FUNC+JC_COST_FUNC
      print*, 'JC_COST_FUNC', JC_COST_FUNC
      END SUBROUTINE CALC_JC_COST

      SUBROUTINE CALC_Cw( w, State_Chm_Adj, State_Grid_Adj )

      USE State_Chm_Adj_Mod,         Only : ChmStateAdj
      USE State_Grid_Adj_Mod,        Only : GrdStateAdj

      TYPE(GrdStateAdj), INTENT(IN)       :: State_Grid_Adj
      TYPE(ChmStateAdj), INTENT(INOUT)    :: State_Chm_Adj
      REAL(fp),          INTENT(INOUT)    :: w( State_Grid_Adj%NX,&
                                                       State_Grid_Adj%NY, &
                                                       State_Grid_Adj%NZ, &
                                                       State_Chm_Adj%nSpecies)
      INTEGER       :: I, J, L, N, NX, NY, NZ, N_TRACERS

      NX = State_Grid_Adj%NX
      NY = State_Grid_Adj%NY
      NZ = State_Grid_Adj%NZ
      N_TRACERS = State_Chm_Adj%nSpecies

      DO N = 1, N_TRACERS
      DO L = 1, NZ
      DO J = 1, NY
      DO I = 1, NX
         w(I, J, L, N) = w(I, J, L, N)*State_Chm_Adj%C(I, J, L, N)
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      END SUBROUTINE CALC_Cw
 
      SUBROUTINE CALC_g( State_Chm_Adj )
      USE State_Chm_Adj_Mod,         Only : ChmStateAdj
      TYPE(ChmStateAdj), INTENT(INOUT)    :: State_Chm_Adj

      INTEGER                             :: I
   
      DO I = 1, 37
      IF (INT(37/2) == I) THEN
         State_Chm_Adj%g(I) = -State_Chm_Adj%f(I)
      ELSE
         State_Chm_Adj%g(I) = 1-State_Chm_Adj%f(I)
      ENDIF
      ENDDO
      END SUBROUTINE CALC_g
       
      END MODULE GEOS_CHEM_BALANCE_MOD

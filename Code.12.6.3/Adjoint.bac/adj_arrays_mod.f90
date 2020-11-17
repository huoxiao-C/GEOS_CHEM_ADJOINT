MODULE ADJ_ARRAYS_MOD

    USE TIME_MOD,             ONLY : GET_ELAPSED_SEC
    USE Precision_Mod

    IMPLICIT NONE
    PUBLIC

      INTEGER                      :: NOPT
      LOGICAL                      :: EXIST_OBS
      INTEGER                      :: N_CALC, N_CALC_STOP
      REAL(fp)                     :: COST_FUNC
      REAL(fp)                     :: COST_FUNC_Jo
      INTEGER                      :: OBS_EXIST(144, 91)
      INTEGER, POINTER             :: INDEX_ARRAY(:, :) 
CONTAINS

!-------------------------------------------------------
      FUNCTION ITS_TIME_FOR_CHECKPT() RESULT( FLAG )

      LOGICAL                  :: FLAG

      FLAG = MOD( GET_ELAPSED_SEC(), 1200 ) == 0 

      END FUNCTION ITS_TIME_FOR_CHECKPT



!-------------------------------------------------------
      FUNCTION ITS_TIME_FOR_OBS() RESULT( FLAG )

      LOGICAL                  :: FLAG

      FLAG = MOD( GET_ELAPSED_SEC(), 1200 ) == 0 

      END FUNCTION ITS_TIME_FOR_OBS

!-------------------------------------------------------
      FUNCTION ITS_TIME_FOR_I3_FIELD_ADJ() RESULT( FLAG )

      LOGICAL                  :: FLAG

      FLAG = MOD( GET_ELAPSED_SEC(), 1200 ) == 0

      END FUNCTION ITS_TIME_FOR_I3_FIELD_ADJ

!-------------------------------------------------------
      FUNCTION ITS_TIME_FOR_DYN_ADJ() RESULT( FLAG )

      LOGICAL                  :: FLAG

      FLAG = MOD( GET_ELAPSED_SEC(), 600 ) == 0

      END FUNCTION ITS_TIME_FOR_DYN_ADJ

!-------------------------------------------------------

      SUBROUTINE CONVERT_UNITS( IFLAG, State_Chm_Adj, State_Grid_Adj,  STT )


      USE State_Chm_Adj_Mod,         Only : ChmStateAdj
      USE State_Grid_Adj_Mod,        Only : GrdStateAdj
      USE State_Met_Adj_Mod,         Only : MetStateAdj

      USE PHYSCONSTANTS,             Only : AIRMW


      TYPE(ChmStateAdj), INTENT(INOUT)    :: State_Chm_Adj
      TYPE(GrdStateAdj), INTENT(IN)       :: State_Grid_Adj

            ! Arguments
      INTEGER, INTENT(IN)    :: IFLAG
      REAL*8,  INTENT(INOUT) :: STT(State_Grid_Adj%NX,State_Grid_Adj%NY,State_Grid_Adj%NZ,State_Chm_Adj%nSpecies)

      INTEGER                  :: IM, JM, KM
      INTEGER                  :: I, J, L, N
      REAL(fp)                 :: MwRatio
      REAL(fp)                 :: MW_g
  

      IM = State_Grid_Adj%NX
      JM = State_Grid_Adj%NY
      KM = State_Grid_Adj%NZ

      SELECT CASE ( IFLAG )

         CASE ( 1 )

            DO N = 1, State_Chm_Adj%nSpecies
               
               MW_g = State_Chm_Adj%emMW_g(N)
               MwRatio = AIRMW/MW_g
            !convert from kg/kg to v/v
            DO L = 1, KM
            DO J = 1, JM
            DO I = 1, IM
               STT(I,J,L,N) = STT(I,J,L,N)*MwRatio
            ENDDO
            ENDDO
            ENDDO
            ENDDO

         !convert from v/v to kg/kg
         CASE ( 2 )

            DO N = 1, State_Chm_Adj%nSpecies

               MW_g = State_Chm_Adj%emMW_g(N)
               MwRatio = AIRMW/MW_g

            DO L = 1, KM
            DO J = 1, JM
            DO I = 1, IM
               STT(I,J,L,N) = STT(I,J,L,N)/MwRatio
            ENDDO
            ENDDO
            ENDDO
            ENDDO
      END SELECT
 
      END SUBROUTINE CONVERT_UNITS

!-------------------------------------------------------
      SUBROUTINE Scale_Adjoint_Value(FLAG, State_Chm_Adj, State_Met_Adj)

      USE State_Chm_Adj_Mod,         Only : ChmStateAdj
      USE State_Met_Adj_Mod,         Only : MetStateAdj

      TYPE(ChmStateAdj), INTENT(INOUT)       :: State_Chm_Adj
      TYPE(MetStateAdj), INTENT(INOUT)       :: State_Met_Adj
      INTEGER, INTENT(IN)                 :: FLAG

      ! current only scale CO2
      IF ( FLAG == 1 ) THEN
      State_Chm_Adj%SpeciesAdj(:, :, :, 1) =  &
                    State_Chm_Adj%SpeciesAdj(:, :, :, 1)/State_Met_Adj%AD(:, :, :)
      ELSEIF ( FLAG == 0 ) THEN
            State_Chm_Adj%SpeciesAdj(:, :, :, 1) =  &
                    State_Chm_Adj%SpeciesAdj(:, :, :, 1)*State_Met_Adj%AD(:, :, :) 
      ELSE
            print*, 'Error: invalid flag, Exit in Scale_Adjoint_Value Subroutine'
            STOP
      ENDIF
      END SUBROUTINE Scale_Adjoint_Value 
          
!-------------------------------------------------------
      FUNCTION ITS_TIME_FOR_EMIS_ADJ() RESULT( FLAG )
  
      LOGICAL                  :: FLAG
      REAL(fp)                 :: MwRation
      REAL(fp)                 :: MW_g

      
      FLAG = MOD( GET_ELAPSED_SEC(), 600 ) == 0

      END FUNCTION ITS_TIME_FOR_EMIS_ADJ
 
END MODULE ADJ_ARRAYS_MOD


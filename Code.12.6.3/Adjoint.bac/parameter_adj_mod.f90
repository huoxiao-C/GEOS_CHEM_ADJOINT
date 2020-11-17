MODULE PARAMETER_ADJ_MOD

       USE Precision_Mod
 
       REAL(fp), ALLOCATABLE           :: TCVV(:)
       REAL(fp)                        :: Mdry

       REAL(fp)                        :: g

       CONTAINS

       SUBROUTINE Init_Parameter(State_Chm_Adj)

       USE State_Chm_Adj_Mod,    ONLY : ChmStateAdj

       TYPE(ChmStateAdj), INTENT(IN)   :: State_Chm_Adj
       INTEGER                         :: RC

       ALLOCATE(TCVV(State_Chm_Adj%nSpecies), STAT=RC)
 
       ! 1 for co2
       TCVV(1) = 28.97d0  / 44d0
       
       Mdry =  29.8
       g = 9.8
       END SUBROUTINE

!-----------------------------------------------------------

       SUBROUTINE Cleanup_Parameter

       IF(ALLOCATED(TCVV))  DEALLOCATE(TCVV)

       END SUBROUTINE Cleanup_Parameter


END MODULE PARAMETER_ADJ_MOD

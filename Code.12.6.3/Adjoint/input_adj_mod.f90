      MODULE Input_Adj_Mod


      TYPE, PUBLIC :: OptInputAdj

      INTEGER      :: IORD
      INTEGER      :: JORD
      INTEGER      :: KORD

      LOGICAL      :: LFILL
      LOGICAL      :: LPRT

      END TYPE OptInputAdj


      CONTAINS

      SUBROUTINE Init_Input_Adj_Mod( Input_Opt_Adj )

      TYPE(OptInputAdj), INTENT(INOUT)      :: Input_Opt_Adj

       
      Input_Opt_Adj%Iord = 3
      Input_Opt_Adj%Jord = 3
      Input_Opt_Adj%Kord = 7
 
      Input_Opt_Adj%LFILL = .FALSE.
      Input_Opt_Adj%LPRT = .TRUE.

      END SUBROUTINE Init_Input_Adj_Mod

      END MODULE Input_Adj_Mod




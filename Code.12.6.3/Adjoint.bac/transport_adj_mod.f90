      MODULE TRANSPORT_ADJ_MOD

      USE PRECISION_MOD      ! For GEOS-Chem Precision (fp)

      IMPLICIT NONE
      PRIVATE

      PUBLIC  :: DO_TRANSPORT_ADJ
      PUBLIC  :: DO_TRANSPORT_D
      PRIVATE :: DO_GLOBAL_ADV

      REAL(fp), ALLOCATABLE           :: Ap (:)
      REAL(fp), ALLOCATABLE           :: Bp (:)
      INTEGER                         :: NG, MG
      INTEGER                         :: JFIRST, JLAST
      INTEGER                         :: N_ADJ

      CONTAINS

      SUBROUTINE DO_TRANSPORT_ADJ(State_Chm_Adj, State_Grid_Adj, &
                                  State_Met_Adj, Input_Opt_Adj)


      USE State_Chm_Adj_Mod,  ONLY : ChmStateAdj
      USE State_Grid_Adj_Mod, ONLY : GrdStateAdj
      USE State_Met_Adj_Mod,  ONLY : MetStateAdj
      USE Input_Adj_Mod,      ONLY : OptInputAdj
      USE TIME_MOD,           ONLY : GET_ELAPSED_SEC
      LOGICAL, SAVE                    :: FIRST=.TRUE.

! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(ChmStateAdj), INTENT(INOUT) :: State_Chm_Adj   ! Chemistry State object
      TYPE(GrdStateAdj), INTENT(INOUT) :: State_Grid_Adj  ! Diagnostics State object
      TYPE(MetStateAdj), INTENT(INOUT) :: State_Met_Adj   ! Meteorology State object

      TYPE(OptInputAdj), INTENT(IN)    :: Input_Opt_Adj


      CHARACTER(LEN=255) :: ErrMsg, ThisLoc
      INTEGER            :: RC

      !adjoint validation
      INTEGER            :: I,J
      ThisLoc =   &
      ' -> at Do_Transport  (in module Adjoint/transport_adj_mod.f90)'
!      print*, ThisLoc
!      print*, 'First', First
      IF ( FIRST ) THEN
         
         CALL INIT_TRANSPORT(State_Grid_Adj,State_Met_Adj, 1)
         FIRST = .FALSE.
  
      END IF

      !
      !tangent validation
!      IF (GET_ELAPSED_SEC() == 3000 ) THEN
!      OPEN(UNIT=1995, FILE='ADV_ADJOINT', ACCESS='APPEND')
!      ENDIF
!      DO I = 50, 60
!      DO J = 40, 50
!      State_Chm_Adj%SpeciesAdj = 0.0
!      State_Chm_Adj%SpeciesAdj(I, J, 1, 1) = 1.0
      CALL DO_GLOBAL_ADV(State_Chm_Adj, State_Grid_Adj, &
                         State_Met_Adj, Input_Opt_Adj, RC, 1)

!      IF (GET_ELAPSED_SEC() == 3000 ) THEN
!      WRITE(1995, *) State_Chm_Adj%SpeciesAdj(I, J, 5, 1)
!      ENDIF

!      ENDDO
!      ENDDO
!      CLOSE(1995) 
      END SUBROUTINE DO_TRANSPORT_ADJ
!--------------------------------------------------------------
      SUBROUTINE DO_TRANSPORT_D(State_Chm_Adj, State_Grid_Adj, &
                                  State_Met_Adj, Input_Opt_Adj)


      USE State_Chm_Adj_Mod,  ONLY : ChmStateAdj
      USE State_Grid_Adj_Mod, ONLY : GrdStateAdj
      USE State_Met_Adj_Mod,  ONLY : MetStateAdj
      USE Input_Adj_Mod,      ONLY : OptInputAdj
      USE TIME_MOD,           ONLY : GET_ELAPSED_SEC
      LOGICAL, SAVE                    :: DFIRST=.TRUE.

! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(ChmStateAdj), INTENT(INOUT) :: State_Chm_Adj   ! Chemistry State object
      TYPE(GrdStateAdj), INTENT(INOUT) :: State_Grid_Adj  ! Diagnostics State object
      TYPE(MetStateAdj), INTENT(INOUT) :: State_Met_Adj   ! Meteorology State object

      TYPE(OptInputAdj), INTENT(IN)    :: Input_Opt_Adj


      CHARACTER(LEN=255) :: ErrMsg, ThisLoc
      INTEGER            :: RC
       
      !tangent validation
      INTEGER            :: I, J
      IF ( DFIRST ) THEN

         CALL INIT_TRANSPORT(State_Grid_Adj,State_Met_Adj, 0)
         DFIRST = .FALSE.

      END IF

!      !tangent validation
!      IF (GET_ELAPSED_SEC() == 3000 ) THEN
!      OPEN(UNIT=1995, FILE='ADV_TANGENT', ACCESS='APPEND')
!      ENDIF 
!      DO I = 50, 60
!      DO J = 40, 50
!      State_Chm_Adj%SpeciesAdj = 0.0
!      State_Chm_Adj%SpeciesAdj(I, J, 5, 1) = 1.0
      CALL DO_GLOBAL_ADV(State_Chm_Adj, State_Grid_Adj, &
                         State_Met_Adj, Input_Opt_Adj, RC, 0) 
!      IF (GET_ELAPSED_SEC() == 3000 ) THEN
!      WRITE(1995, *) State_Chm_Adj%SpeciesAdj(I, J, 1, 1)
!      ENDIF

!      ENDDO
!      ENDDO
!      CLOSE(1995)
      END SUBROUTINE DO_TRANSPORT_D
!--------------------------------------------------------------

!--------------------------------------------------------------
      SUBROUTINE DO_GLOBAL_ADV(State_Chm_Adj, State_Grid_Adj, &
                               State_Met_Adj, Input_Opt_Adj, RC, &
                               Adj_Flag)


      USE State_Chm_Adj_Mod,  ONLY : ChmStateAdj
      USE State_Grid_Adj_Mod, ONLY : GrdStateAdj
      USE State_Met_Adj_Mod,  ONLY : MetStateAdj
      USE Input_Adj_Mod,      ONLY : OptInputAdj
      USE TIME_MOD,           ONLY : GET_TS_DYN
      USE LOGICAL_ADJ_MOD,    ONLY : LICS_INC
      USE ErrCode_Adj_Mod
      USE PhysConstants            ! Physical constants
      USE CMN_SIZE_MOD             ! Size parameters
      USE Parameter_Adj_Mod,  ONLY : TCVV

      USE TPCORE_FVDAS_ADJ_MOD, ONLY : TPCORE_FVDAS
      USE TPCORE_FVDAS_ADJ_MOD_DIFF, ONLY : TPCORE_FVDAS_D
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
      INTEGER,        INTENT(IN)    :: Adj_Flag    ! 1 tangent; 0 adjoint
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(ChmStateAdj), INTENT(INOUT) :: State_Chm_Adj   ! Chemistry State object
      TYPE(GrdStateAdj), INTENT(INOUT) :: State_Grid_Adj  ! Diagnostics State object
      TYPE(MetStateAdj), INTENT(INOUT) :: State_Met_Adj   ! Meteorology State object

      TYPE(OptInputAdj), INTENT(IN)    :: Input_Opt_Adj


      ! Scalars
      LOGICAL            :: LFILL
      LOGICAL            :: LPRT
      INTEGER            :: IORD, JORD, KORD
      INTEGER            :: I, J, L, L2, N, N_DYN, NA, nAdvect
      INTEGER            :: ND24, ND25, ND26
      REAL(fp)           :: D_DYN

      ! dimension
      INTEGER            :: IM, JM, KM, N_TRACERS
      ! Arrays
      REAL(fp), POINTER  :: P_TP1 (:, :)

      REAL(fp), POINTER  :: P_TP2 (:, :)

      REAL(fp), POINTER  :: XMASS (:, :, :)

      REAL(fp), POINTER  :: YMASS (:, :, :)

      REAL(fp), POINTER  :: UWND (:, :, :)

      REAL(fp), POINTER  :: VWND (:, :, :)

      REAL(fp), POINTER  :: p_Sec(:, :, :, :)

      REAL(fp), POINTER  :: q(:, :, :, :)
      REAL(fp), POINTER  :: A_M2(:, :)

      REAL(fp)           :: MFLEW(State_Grid_Adj%NX, State_Grid_Adj%NY, &
                                  State_Grid_Adj%NZ, State_Chm_Adj%nSpecies)

      REAL(fp)           :: MFLNS(State_Grid_Adj%NX, State_Grid_Adj%NY, &
                                  State_Grid_Adj%NZ, State_Chm_Adj%nSpecies)

      REAL(fp)           :: MFLUP(State_Grid_Adj%NX, State_Grid_Adj%NY, &
                                  State_Grid_Adj%NZ, State_Chm_Adj%nSpecies)

      REAL(fp)           :: P_TEMP(State_Grid_Adj%NX,State_Grid_Adj%NY)
    
      REAL(fp)           :: MINVALUE(State_Chm_Adj%nSpecies)

      INTEGER            :: DI, DJ, DK

      !
      ! Assume success
      RC          =  GC_SUCCESS

      ! Initialize
      LFILL       =  Input_Opt_Adj%LFILL
      LPRT        =  Input_Opt_Adj%LPRT
      IORD        =  Input_Opt_Adj%IORD
      JORD        =  Input_Opt_Adj%JORD
      KORD        =  Input_Opt_Adj%KORD
      nAdvect     =  State_Chm_Adj%nAdvect

      N_TRACERS   = State_Chm_Adj%nSpecies
      !no diagnostic output
      ND24 = 0
      ND25 = 0
      ND26 = 0 
      N_ADJ = State_Chm_Adj%nSpecies
      !p_MFLEW     => NULL()
      ! p_MFLNS     => NULL()
      ! p_MFLUP     => NULL()


      ! Dynamic timestep [s]
      N_DYN       =  GET_TS_DYN()
      D_DYN       =  DBLE( N_DYN )



      P_TP1   => State_Met_Adj%PEDGE_DRY(:, :)
      P_TP2   => State_Met_Adj%PSC2_DRY(:, :)
      ! Flip array indices in the vertical using pointer storage
      ! continuous adjoint wind reverse
      State_Met_Adj%U = -State_Met_Adj%U
      State_Met_Adj%V = -State_Met_Adj%V
      UWND    => State_Met_Adj%U          (:, :, State_Grid_Adj%NZ:1:-1)
      VWND    => State_Met_Adj%V          (:, :, State_Grid_Adj%NZ:1:-1)
      XMASS   => State_Met_Adj%XMASS      (:, :, State_Grid_Adj%NZ:1:-1)
      YMASS   => State_Met_Adj%YMASS      (:, :, State_Grid_Adj%NZ:1:-1)


      A_M2    => State_Grid_Adj%AREA_M2   (:, :) 
      IM      = State_Grid_Adj%NX
      JM      = State_Grid_Adj%NY
      KM      = State_Grid_Adj%NZ
!      p_Sec   =>  State_Chm_Adj%SpeciesAdj(:, :, State_Grid_Adj%NZ:1:-1, :)
      p_Sec => State_Chm_Adj%SpeciesAdj(:, :, State_Grid_Adj%NZ:1:-1, :)
      q     => State_Chm_Adj%CO2_ADV(:, :, State_Grid_Adj%NZ:1:-1, :)
      ! fill value
      IF( LFILL .AND. Adj_Flag==1) THEN

        DO N = 1, N_TRACERS

           MINVALUE(N)       = MINVAL(p_Sec(:, :, :, N))
           p_Sec(:, :, :, N) = p_Sec(:, :, :, N)+ABS(MINVALUE(N))*2.0
        END DO

      ENDIF
      ! Do the advection
      ! Note: the mass flux diagnostic arrays (MASSFLEW, MASSFLNS and MASSFLUP)
      ! are incremented upside-down (level 1 = top of the atmosphere).
      ! The levels order is reversed only when written out in diag3.f
      ! (ccc, 3/8/10), (adj32_022)
      IF ( Adj_Flag == 1) THEN
      CALL TPCORE_FVDAS( D_DYN,    Re,        IM,    JM,         &
                         KM,    JFIRST,    JLAST,    NG,             &
                         MG,       N_TRACERS, Ap,       Bp,          &
                         UWND,     VWND,      P_TP1,    P_TP2,       &
                         P_TEMP,   p_Sec,  & !State_Chm_Adj%SpeciesAdj          &
                        ! (:, :, State_Grid_Adj%NZ:1:-1, :),          &
                         IORD,     JORD,      KORD,     N_ADJ,       &
                         XMASS,                                      &
                         YMASS,                                      &
!     &                   MASSFLEW(:,:,LLPAR:1:-1,:),
!     &                   MASSFLNS(:,:,LLPAR:1:-1,:),
!     &                   MASSFLUP(:,:,LLPAR:1:-1,:),    A_M2,
                         MFLEW,                                   &
                         MFLNS,                                   &
                         MFLUP,    A_M2,                          &
                         TCVV,     ND24,      ND25,     ND26,        &
                         LFILL )     

      ELSEIF ( Adj_Flag == 0 ) THEN
          CALL TPCORE_FVDAS_D( D_DYN,    Re,        IM,    JM,         &
                         KM,    JFIRST,    JLAST,    NG,             &
                         MG,       N_TRACERS, Ap,       Bp,          &
                         UWND,     VWND,      P_TP1,    P_TP2,       &
                         P_TEMP,   q,         p_Sec,  & !State_Chm_Adj%SpeciesAdj          &
                        ! (:, :, State_Grid_Adj%NZ:1:-1, :),          &
                         IORD,     JORD,      KORD,     N_ADJ,       &
                         XMASS,                                      &
                         YMASS,                                      &
!     &                   MASSFLEW(:,:,LLPAR:1:-1,:),
!     &                   MASSFLNS(:,:,LLPAR:1:-1,:),
!     &                   MASSFLUP(:,:,LLPAR:1:-1,:),    A_M2,
                         MFLEW,                                   &
                         MFLNS,                                   &
                         MFLUP,    A_M2,                          &
                         TCVV,     ND24,      ND25,     ND26,        &
                         LFILL )
       ENDIF

      ! fill va
      IF( LFILL .AND. Adj_Flag==1 ) THEN

        DO N = 1, N_TRACERS

           p_Sec(:, :, :, N) = p_Sec(:, :, :, N)-ABS(MINVALUE(N))*2.0

        END DO

      ENDIF

      ! restore negative wind
      State_Met_Adj%U = -State_Met_Adj%U
      State_Met_Adj%V = -State_Met_Adj%V

      P_TP1 => NULL()
      P_TP2 => NULL()
      UWND  => NULL()
      VWND  => NULL()
      XMASS => NULL()
      YMASS => NULL()
      A_M2  => NULL()
      p_Sec => NULL()
      END SUBROUTINE DO_GLOBAL_ADV

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

      SUBROUTINE INIT_TRANSPORT(State_Grid_Adj,State_Met_Adj, Adj_Flag)

      !USE
      USE TPCORE_FVDAS_ADJ_MOD, ONLY : INIT_TPCORE
      USE TPCORE_FVDAS_ADJ_MOD_DIFF, ONLY : INIT_TPCORE_D
      USE State_Grid_Adj_Mod,   ONLY : GrdStateAdj
      USE State_Met_Adj_Mod,    ONLY : MetStateAdj
      USE Time_Mod,             ONLY : GET_TS_DYN
      USE PHYSCONSTANTS
      INTEGER, INTENT(IN)              :: Adj_Flag
      TYPE(GrdStateAdj), INTENT(INOUT) :: State_Grid_Adj  ! Diagnostics State object
      TYPE(MetStateAdj), INTENT(INOUT) :: State_Met_Adj   ! Meteorology State object

      !local variables
      REAL*8                           :: YMID_R(State_Grid_Adj%NY)
      INTEGER                          :: RC, N_DYN
      INTEGER                          :: IM, JM, KM

     
      IM      = State_Grid_Adj%NX
      JM      = State_Grid_Adj%NY
      KM      = State_Grid_Adj%NZ

      ALLOCATE( Ap( State_Grid_Adj%NZ+1 ), STAT=RC )
      ALLOCATE( Bp( State_Grid_Adj%NZ+1 ), STAT=RC )

      Ap      = State_Met_Adj%Ap(State_Grid_Adj%NZ+1:1:-1)
      Bp      = State_Met_Adj%Bp(State_Grid_Adj%NZ+1:1:-1)
     
      ! Initialize
      N_DYN = GET_TS_DYN() * 60
      N_ADJ = 0
      NG    = 0
      MG    = 0

      YMID_R = State_Grid_Adj%YMid_R(1, :)
      ! Call INIT routine from "tpcore_fvdas_mod.f"
      if ( Adj_Flag == 1 ) THEN
      CALL INIT_TPCORE( IM,  JM, KM,  JFIRST, JLAST,                 &
                        NG, MG, DBLE( N_DYN ), Re,     YMID_R )
      ELSEIF ( Adj_Flag == 0 ) THEN 
      CALL INIT_TPCORE_D( IM,  JM, KM,  JFIRST, JLAST,                 &
                        NG, MG, DBLE( N_DYN ), Re,     YMID_R )
      ENDIF

      END SUBROUTINE INIT_TRANSPORT

      ! end

      END MODULE TRANSPORT_ADJ_MOD


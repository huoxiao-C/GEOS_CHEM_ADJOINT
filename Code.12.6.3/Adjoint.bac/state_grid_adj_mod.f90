
MODULE State_Grid_Adj_Mod
!
! USES:
!
  USE ErrCode_Adj_Mod
  USE Precision_Mod

  IMPLICIT NONE
  PRIVATE

!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: Init_State_Grid_Adj
  PUBLIC :: Allocate_State_Grid_Adj
  PUBLIC :: Cleanup_State_Grid_Adj
  PUBLIC :: Calculate_Grid_Parameter

  TYPE, PUBLIC :: GrdStateAdj
     !----------------------------------------
     ! User-defined grid fields
     !----------------------------------------
     REAL(fp)           :: DX          ! Delta X         [degrees longitude]
     REAL(fp)           :: DY          ! Delta Y         [degrees latitude]
     INTEGER            :: NX          ! # of grid boxes in X-direction
     INTEGER            :: NY          ! # of grid boxes in Y-direction
     INTEGER            :: NZ          ! # of grid boxes in Z-direction
     INTEGER            :: YMinOffset  ! Y offset from global grid

     ! Arrays
     REAL(fp),  POINTER :: Area_M2  (:,:) ! Grid box area [m2]
     REAL(fp),  POINTER :: YMid     (:,:) ! latitude center in grid 
     REAL(fp),  POINTER :: YMid_R   (:,:) ! latitude center (radian) in grid 

     !logical
     LOGICAL            :: HalfPolar
     LOGICAL            :: NestedGrid


  END TYPE GrdStateAdj

CONTAINS

  SUBROUTINE Init_State_Grid_Adj( State_Grid_Adj, RC )

! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(GrdStateAdj), INTENT(INOUT) :: State_Grid_Adj   ! Obj for grid state
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Return code

    ! Assume success
    RC            =  GC_SUCCESS

    !----------------------------------------
    ! User-defined grid fields
    !----------------------------------------
    State_Grid_Adj%NX           = 144
    State_Grid_Adj%NY           = 91
    State_Grid_Adj%NZ           = 47
    State_Grid_Adj%DX           = 2.5
    State_Grid_Adj%DY           = 2.0

    State_Grid_Adj%YMinOffset   = 0
    !--------------------------------------------
    !input parameter
    !--------------------------------------------
    State_Grid_Adj%NestedGrid = .FALSE.
    State_Grid_Adj%HalfPolar  = .FALSE.
    !---------------------------------------------------------------
    ! Nullify all fields for safety's sake before allocating them
    !---------------------------------------------------------------
    State_Grid_Adj%Area_M2      => NULL()
    State_Grid_Adj%YMid         => NULL()
    State_Grid_Adj%YMid_R       => NULL()

  END SUBROUTINE Init_State_Grid_Adj


!
  SUBROUTINE Allocate_State_Grid_Adj( State_Grid_Adj, RC )

! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(GrdStateAdj), INTENT(INOUT) :: State_Grid_Adj  ! Grid State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?

    ALLOCATE( State_Grid_Adj%Area_M2( State_Grid_Adj%NX, State_Grid_Adj%NY ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Grid_Adj%Area_M2 = 0e+0_fp

    ALLOCATE( State_Grid_Adj%YMid( State_Grid_Adj%NX, State_Grid_Adj%NY ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Grid_Adj%Area_M2 = 0e+0_fp

    ALLOCATE( State_Grid_Adj%YMid_R( State_Grid_Adj%NX, State_Grid_Adj%NY ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Grid_Adj%Area_M2 = 0e+0_fp

  END SUBROUTINE Allocate_State_Grid_Adj 
  
!
  SUBROUTINE Cleanup_State_Grid_Adj( State_Grid_Adj, RC )
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(GrdStateAdj), INTENT(INOUT) :: State_Grid_Adj   ! Obj for grid state
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Return code

    ! Assume success
    RC      = GC_SUCCESS
    !=======================================================================
    ! Deallocate arrays
    !=======================================================================
    IF ( ASSOCIATED( State_Grid_Adj%Area_M2 ) ) THEN
       DEALLOCATE( State_Grid_Adj%Area_M2, STAT=RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Grid_Adj%Area_M2 => NULL()
    ENDIF

  END SUBROUTINE Cleanup_State_Grid_Adj


  SUBROUTINE Calculate_Grid_Parameter( State_Grid_Adj, RC )

  !USE
  
    USE PhysConstants 
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(GrdStateAdj), INTENT(INOUT) :: State_Grid_Adj   ! Obj for grid state
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Return code

    INTEGER                       :: I, J, JG

    DO J = 1, State_Grid_Adj%NY
    DO I = 1, State_Grid_Adj%NX


       ! Index value for user-defined grid on the global grid
       JG = J + ( State_Grid_Adj%YMinOffset - 1 )
       !--------------------------------
       ! Latitude centers [degrees]
       !--------------------------------
       State_Grid_Adj%YMid(I,J) = ( State_Grid_Adj%DY * JG ) - 90e+0_fp
       ! If using half-sized polar boxes on a global grid,
       ! multiply DY by 1/4 at poles
       IF ( State_Grid_Adj%HalfPolar .and. .not. State_Grid_Adj%NestedGrid .and. &
           J == 1 ) THEN
          ! South Pole
          State_Grid_Adj%YMid(I,J) = -90e+0_fp + ( 0.25e+0_fp * State_Grid_Adj%DY )
       ENDIF

       IF ( State_Grid_Adj%HalfPolar .and. .not. State_Grid_Adj%NestedGrid .and. &
            J == State_Grid_Adj%NY ) THEN
          ! North Pole
          State_Grid_Adj%YMid(I,J) = +90e+0_fp - ( 0.25e+0_fp * State_Grid_Adj%DY )
       ENDIF

       !--------------------------------
       ! Latitude centers [radians]
       !--------------------------------
       State_Grid_Adj%YMid_R(I,J) = ( PI_180 * State_Grid_Adj%YMid(I,J)  )

   END DO
   END DO

   END SUBROUTINE Calculate_Grid_Parameter
END MODULE State_Grid_Adj_Mod

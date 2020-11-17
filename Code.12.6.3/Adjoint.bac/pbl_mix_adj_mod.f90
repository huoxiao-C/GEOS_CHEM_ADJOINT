      MODULE PBL_MIX_ADJ_MOD

!
      IMPLICIT NONE
      ! Make everything PRIVATE ...
      PRIVATE

      ! ... except these routines
      PUBLIC :: DO_PBL_MIX_ADJ
      PUBLIC :: DO_PBL_MIX_D
      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement
      !=================================================================
      CONTAINS


!------------------------------------------------------------------------------

      SUBROUTINE DO_PBL_MIX_ADJ( DO_TURBDAY, State_Chm_Adj, State_Grid_Adj, State_Met_Adj )
!

!      USE PBL_MIX_MOD,    ONLY : COMPUTE_PBL_HEIGHT
      USE State_Chm_Adj_Mod,      Only : ChmStateAdj
      USE State_Grid_Adj_Mod,     Only : GrdStateAdj
      USE State_Met_Adj_Mod,      Only : MetStateAdj
!      USE PRESSURE_MOD,    ONLY : GET_PEDGE
      USE Parameter_Adj_Mod,      Only : TCVV
      USE Adj_Arrays_Mod,         Only : CONVERT_UNITS
 
      ! convection adjoint for co2 
      TYPE(MetStateAdj), INTENT(IN)    :: State_Met_Adj   ! Meteorology State object
      TYPE(GrdStateAdj), INTENT(IN)    :: State_Grid_Adj  ! Grid State object
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(ChmStateAdj), INTENT(INOUT) :: State_Chm_Adj   ! Chemistry State object

      ! Arguments
      LOGICAL, INTENT(IN) :: DO_TURBDAY

      ! Now recompute these rather than checkpoint. (dkh, 07/08/09)
      ! Compute PBL height and related quantities
      ! adjoint: read from the file 
!      CALL COMPUTE_PBL_HEIGHT

      CALL CONVERT_UNITS(2, State_Chm_Adj, State_Grid_Adj, &
                         State_Chm_Adj%SpeciesAdj) 
      ! Do complete mixing of tracers in the PBL (if necessary)
      IF ( DO_TURBDAY ) THEN
         CALL TURBDAY_ADJ( State_Chm_Adj%nSpecies, TCVV(1), &
                                        State_Chm_Adj, State_Grid_Adj, State_Met_Adj )
      ENDIF

      CALL CONVERT_UNITS(1, State_Chm_Adj, State_Grid_Adj, &
                         State_Chm_Adj%SpeciesAdj)

      ! Return to calling program
      END SUBROUTINE DO_PBL_MIX_ADJ

      SUBROUTINE DO_PBL_MIX_D( DO_TURBDAY, State_Chm_Adj, State_Grid_Adj, State_Met_Adj )
!

!      USE PBL_MIX_MOD,    ONLY : COMPUTE_PBL_HEIGHT
      USE State_Chm_Adj_Mod,      Only : ChmStateAdj
      USE State_Grid_Adj_Mod,     Only : GrdStateAdj
      USE State_Met_Adj_Mod,      Only : MetStateAdj
!      USE PRESSURE_MOD,    ONLY : GET_PEDGE
      USE Parameter_Adj_Mod,      Only : TCVV
      USE Adj_Arrays_Mod,         Only : CONVERT_UNITS
      ! convection adjoint for co2
      TYPE(MetStateAdj), INTENT(IN)    :: State_Met_Adj   ! Meteorology State object
      TYPE(GrdStateAdj), INTENT(IN)    :: State_Grid_Adj  ! Grid State object
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(ChmStateAdj), INTENT(INOUT) :: State_Chm_Adj   ! Chemistry State object

      ! Arguments
      LOGICAL, INTENT(IN) :: DO_TURBDAY

      CALL CONVERT_UNITS(1, State_Chm_Adj, State_Grid_Adj, &
                         State_Chm_Adj%SpeciesAdj)
      ! Do complete mixing of tracers in the PBL (if necessary)
      IF ( DO_TURBDAY ) THEN
         CALL TURBDAY_D( State_Chm_Adj%nSpecies, TCVV(1), &
                                        State_Chm_Adj, State_Grid_Adj, State_Met_Adj )
      ENDIF
      CALL CONVERT_UNITS(2, State_Chm_Adj, State_Grid_Adj, &
                         State_Chm_Adj%SpeciesAdj)

      ! Return to calling program
      END SUBROUTINE DO_PBL_MIX_D
!------------------------------------------------------------
      SUBROUTINE TURBDAY_D(NTRC, TCVV, State_Chm_Adj, State_Grid_Adj, State_Met_Adj)
      ! References to F90 modules
      USE TIME_MOD,               ONLY : GET_TS_CONV
      USE State_Chm_Adj_Mod,      ONLY : ChmStateAdj
      USE State_Grid_Adj_Mod,     ONLY : GrdStateAdj
      USE State_Met_Adj_Mod,      ONLY : MetStateAdj
      USE TIME_MOD,               ONLY : Get_Elapsed_Sec

      ! convection adjoint for co2
      TYPE(MetStateAdj), INTENT(IN)    :: State_Met_Adj   ! Meteorology State object
      TYPE(GrdStateAdj), INTENT(IN)    :: State_Grid_Adj  ! Grid State object
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(ChmStateAdj), INTENT(INOUT) :: State_Chm_Adj   ! Chemistry State object


      ! Argument variables
      INTEGER,  INTENT(IN)     :: NTRC
      REAL*8,   INTENT(IN)     :: TCVV

      ! Local variables
      INTEGER     :: I, J, L, LTOP, N
      REAL*8      :: AA,  CC, CC_AA,   BLTOP
      REAL*8      :: PW,  PS, AREA_M2, DTCONV
      REAL*8      :: P(0:State_Grid_Adj%NZ)
      REAL*8      :: A(State_Grid_Adj%NX, State_Grid_Adj%NY)
      REAL*8      :: DTC(State_Grid_Adj%NX, State_Grid_Adj%NY, State_Grid_Adj%NZ,NTRC)
      REAL :: result1
      REAL :: result2
      REAL*8 :: temp
      ! tangent variables
      REAL*8 :: ccd, cc_aad
      REAL*8 :: dtcd(state_grid_adj%nx, state_grid_adj%ny, state_grid_adj%nz, &
                  ntrc)
      REAL*8, POINTER :: tcd(:, :, :, :)

      ! dimension
      INTEGER                  :: NX, NY, NZ

      A(:,:) = 1d0

      NX = State_Grid_Adj%NX
      NY = State_Grid_Adj%NY
      NZ = State_Grid_Adj%NZ
      dtcd = 0.0
      tcd => State_Chm_Adj%SpeciesAdj(:, :, :, :)
      ! Loop over Lat/Long grid boxes (I,J)
      DO J = 1, NY
      DO I = 1, NX
      !tangent validation 
!      IF ( GET_ELAPSED_SEC() == 3000) THEN
!      OPEN(UNIT = 1995, FILE='PBL_TANGENT_5', ACCESS='APPEND')
!      tcd = 0.0
!      tcd(I, J, 5, 1) = 1.0
!      ENDIF
      !tangent end

         ! Calculate air mass within PBL at grid box (I,J,L)
         AA = 0.d0
         DO L = 1, State_Met_Adj%IMIX(I,J)-1
            AA = AA + State_Met_Adj%AD(I,J,L)
         ENDDO

         L  = State_Met_Adj%IMIX(I,J)
         AA = AA + State_Met_Adj%AD(I,J,L) * State_Met_Adj%FPBL(I,J)
         ! Loop over tracers
         DO N = 1, NTRC
            ccd = 0.0
            DO l=1,State_Met_Adj%IMIX(i,j)-1
               ccd = ccd + State_Met_Adj%ad(i, j, l)*tcd(i, j, l, n)
            ENDDO
            l = State_Met_Adj%IMIX(i,j)
            result1 = State_Met_Adj%FPBL(i, j)
            temp = State_Met_Adj%ad(i, j, l)*result1
            ccd = ccd + temp*tcd(i, j, l, n)
! For grid boxes (I,J,L) which lie below the PBL top
            cc_aad = ccd/aa  
            DO L = 1, State_Met_Adj%IMIX(I,J)-1
               temp = a(i, j)*State_Met_Adj%ad(i, j, l)
               dtcd(i, j, l, n) = State_Met_Adj%ad(i, j, l)*a(i, j)*&
                                  cc_aad - temp*tcd(i, j, l, n)
               tcd(i, j, l, n) = tcd(i, j, l, n) + dtcd(i, j, l, n) &
                                 /State_Met_Adj%ad(i, j,l)
            ENDDO
            l = State_Met_Adj%IMIX(i,j)            
            result1 = State_Met_Adj%FPBL(i, j)
            result2 = State_Met_Adj%FPBL(i, j)
            temp = a(i, j)*result2*State_Met_Adj%ad(i, j, l)
            dtcd(i, j, l, n) = State_Met_Adj%ad(i, j, l)*a(i, j)*&
                               result1*cc_aad - temp*tcd(i, j, l, n)
            tcd(i, j, l, n) = tcd(i, j, l, n) + dtcd(i, j, l, n)/ &
            State_Met_Adj%ad(i, j, l)
         ENDDO    !N
      !tangent validation
!      IF ( I >50 .AND. I < 60 .AND. J>40 .AND. J<50 .and.&
!           GET_ELAPSED_SEC() == 3000) THEN
!      WRITE(1995, *) tcd(I, J, 1, 1)
!      ENDIF
      !tangent end
      ENDDO       !I
      ENDDO       !J

      !tangent end
      tcd => NULL()
      END SUBROUTINE TURBDAY_D
!------------------------------------------------------------
      SUBROUTINE TURBDAY_ADJ(NTRC, TCVV, State_Chm_Adj, State_Grid_Adj, State_Met_Adj)
!
!******************************************************************************
!  Subroutine TURBDAY_ADJ executes the adjoint of the GEOS-CTM dry convection
!  / boundary layer mixing algorithm from TURBDAY. It is a combination of the
!  forward code TURBDAY with TAMC generated adjoint of the loop over N.
!  See notes in turbday.f for info about the original forward code, and
!  below for notes on modifications made for the adjoint version.
!  (dkh, 10/30/04)
!
!  Arguments as Input:
!  ===========================================================================
!  (1 ) NTRC    : Number of tracers used in computation  [1 to NNPAR]
!  (2 ) TCVV    : mol. wt. air / mol. wt. tracer
!
!  Modue variable Input / Output
!  ======================================================================
!  (1 ) STT_ADJ :  Adjoint tracer array
!
!  NOTES:
!  (1 ) Rather than save / write / read the info from the forward run of TURBDAY,
!        we will just recompute most of it, hence most of the original code for
!        TURBDAY is part of ADJ_TURBDAY.  However, some alterations were made to
!        the forward code.
!        Changes to forward code:
!               -  argument list just NTRC and TCVV, which are passed
!                   the values of NADJ and ADJ_TCVV, respectfully.
!               -  add reference to CMN_ADJ for ADJ_STT
!               -  no ND15 diagnostic update, so get rid of USE DAO_MOD
!               -  get rid of XTRA2, LTURB
!               -  get rid of initial print out
!  (2 )  TAMC (and modified TAMC) code is lower case.
!        Changes to TAMC code:
!               -  The varialbes TC_IN and TC_OUT were used to construct the
!                   adjoint, but they are not needed here.  Just replace them
!                   with ADJ_STT
!               -  Replade multiple do loops with ":" operator (so no longer
!                   need integers ip1,ip2,ip3,ip4)
!               -  Force variables explicitly to double precision using .d0
!               -  Initialize and update global adjoint variables (adtc, addtc)
!                   before and after the PARALLEL DO loop
!  (3 ) Updated for v8 adjoint (dkh, 07/14/09)
!******************************************************************************
!
      ! References to F90 modules
      USE TIME_MOD,        ONLY : GET_TS_CONV
      USE State_Chm_Adj_Mod,      ONLY : ChmStateAdj
      USE State_Grid_Adj_Mod,     ONLY : GrdStateAdj
      USE State_Met_Adj_Mod,      ONLY : MetStateAdj
      USE TIME_MOD,               ONLY : Get_Elapsed_Sec 
      ! convection adjoint for co2 
      TYPE(MetStateAdj), INTENT(IN)    :: State_Met_Adj   ! Meteorology State object
      TYPE(GrdStateAdj), INTENT(IN)    :: State_Grid_Adj  ! Grid State object
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(ChmStateAdj), INTENT(INOUT) :: State_Chm_Adj   ! Chemistry State object

      ! dkh debug


      ! Argument variables
      INTEGER,  INTENT(IN)     :: NTRC
      REAL*8,   INTENT(IN)     :: TCVV

      ! Local variables
      INTEGER     :: I, J, L, LTOP, N
      REAL*8      :: AA,  CC, CC_AA,   BLTOP
      REAL*8      :: PW,  PS, AREA_M2, DTCONV
      REAL*8      :: P(0:State_Grid_Adj%NZ)
      REAL*8      :: A(State_Grid_Adj%NX, State_Grid_Adj%NY)
      REAL*8      :: DTC(State_Grid_Adj%NX, State_Grid_Adj%NY, State_Grid_Adj%NZ,NTRC)

      ! Adjoint variables
      real*8 adcc
      real*8 adcc_aa
      real*8 adtc(State_Grid_Adj%NX, State_Grid_Adj%NY, State_Grid_Adj%NZ,ntrc)
      real*8 adtc_in(State_Grid_Adj%NX, State_Grid_Adj%NY, State_Grid_Adj%NZ,ntrc)
      real*8 adtc_out(State_Grid_Adj%NX, State_Grid_Adj%NY, State_Grid_Adj%NZ,ntrc)
      ! dimension
      INTEGER                  :: NX, NY, NZ
      !=================================================================
      ! TURBDAY_ADJ begins here!
      !=================================================================

      ! Echo some input to the screen
!      WRITE( 6, '(a)'   ) '       -- TURBDAY_ADJ'

      ! Don't need DTCONV for adjoint calculation
      ! Convection timestep [s]
      !DTCONV = GET_TS_CONV() * 60d0

      ! We assume full mixing in the boundary layer, so the A
      ! coefficients are 1 everywhere, day & night (bmy, 2/11/03)
      A(:,:) = 1d0
      NX = State_Grid_Adj%NX
      NY = State_Grid_Adj%NY
      NZ = State_Grid_Adj%NZ 
      !----------------------------------------------
      ! SET GLOBAL ADJOINT VARIABLES
     !----------------------------------------------
      adtc(:,:,:,:)      = 0.d0
      adtc_in(:,:,:,:)   = 0.d0
      adtc_out(:,:,:,:)  = State_Chm_Adj%SpeciesAdj(:, :, :, :)


      ! Loop over Lat/Long grid boxes (I,J)
      DO J = 1, NY
      DO I = 1, NX

      !adjoint validation
!      IF ( GET_ELAPSED_SEC() == 3000) THEN
!      OPEN(UNIT = 1995, FILE='PBL_ADJOINT_5', ACCESS='APPEND')
!      adtc_out = 0.0
!      adtc_out(i, j, 1, 1) = 1.0
!      ENDIF
      !adjoint end
         ! Calculate air mass within PBL at grid box (I,J,L)
         AA = 0.d0
         DO L = 1, State_Met_Adj%IMIX(I,J)-1
            AA = AA + State_Met_Adj%AD(I,J,L)
         ENDDO

         L  = State_Met_Adj%IMIX(I,J)
         AA = AA + State_Met_Adj%AD(I,J,L) * State_Met_Adj%FPBL(I,J)

         ! Loop over tracers
         DO N = 1, NTRC

            !----------------------------------------------
            ! RESET LOCAL ADJOINT VARIABLES
            !----------------------------------------------
            adcc =             0.d0
            adcc_aa =          0.d0

            !----------------------------------------------
            ! ADJOINT ROUTINE BODY
            !----------------------------------------------
            l = State_Met_Adj%IMIX(i,j)
            adtc(i,j,:,n) = adtc(i,j,:,n)+adtc_out(i,j,:,n)
            adtc_out(i,j,:,n) = 0.d0
            adcc_aa = adcc_aa+adtc(i,j,l,n)*a(i,j)*State_Met_Adj%FPBL(i,j)
            adtc(i,j,l,n) = adtc(i,j,l,n)*(1.d0-a(i,j)*State_Met_Adj%FPBL(i,j))
            do l = 1, State_Met_Adj%IMIX(i,j)-1
              adcc_aa = adcc_aa+adtc(i,j,l,n)*a(i,j)
              adtc(i,j,l,n) = adtc(i,j,l,n)*(1.d0-a(i,j))
            end do
            adcc = adcc+adcc_aa/aa
            adcc_aa = 0.d0
            l = State_Met_Adj%IMIX(i,j)
            adtc(i,j,l,n) = adtc(i,j,l,n)+adcc*State_Met_Adj%ad(i,j,l)*State_Met_Adj%FPBL(i,j)
            do l = 1, State_Met_Adj%IMIX(i,j)-1
              adtc(i,j,l,n) = adtc(i,j,l,n)+adcc*State_Met_Adj%ad(i,j,l)
            end do
            adtc_in(i,j,:,n) = adtc_in(i,j,:,n)+adtc(i,j,:,n)
            adtc(i,j,:,n) = 0.d0

         ENDDO    !N
      !adjoint validation
!      IF ( I >50 .AND. I < 60 .AND. J>40 .AND. J<50 .AND. &
!          GET_ELAPSED_SEC() == 3000) THEN
!      WRITE(1995, *) adtc_in(I, J, 1, 1)
!       WRITE(1995, *) adtc_in(I, J, 5, 1)
!      ENDIF
      !adjoint end
      ENDDO       !I
      ENDDO       !J

      ! Update global adjoint variables
      State_Chm_Adj%SpeciesAdj(:,:,:,:) = adtc_in(:,:,:,:)

      adtc_in(:,:,:,:) = 0d0

      !  Return to calling program
      END SUBROUTINE TURBDAY_ADJ
      END MODULE PBL_MIX_ADJ_MOD

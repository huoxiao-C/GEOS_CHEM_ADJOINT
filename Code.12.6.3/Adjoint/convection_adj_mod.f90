!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model For Co2                  !
!------------------------------------------------------------------------------
! huoxiao adjoint for co2
! convection for geos-fp and merra2
      MODULE CONVECTION_ADJ_MOD

      USE PRECISION_MOD    ! For GEOS-Chem Precision (fp)

      IMPLICIT NONE
      PRIVATE

!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: DO_CONVECTION_ADJ
      PUBLIC  :: DO_CONVECTION_D
!
! !PRIVATE MEMBER FUNCTIONS:
!
      PRIVATE :: DO_CLOUD_CONVECTION_ADJ

      CONTAINS

      SUBROUTINE DO_CONVECTION_ADJ( State_Chm_Adj, State_Grid_Adj,     &
                                State_Met_Adj, RC )


      USE TIME_MOD,        ONLY : GET_TS_DYN
      USE TIME_MOD,        ONLY : GET_ELAPSED_SEC
      USE State_Chm_Adj_Mod,      ONLY : ChmStateAdj
!      USE State_Diag_Mod,     ONLY : DgnState
      USE State_Grid_Adj_Mod,     ONLY : GrdStateAdj
      USE State_Met_Adj_Mod,      ONLY : MetStateAdj

      ! convection adjoint for co2 
      TYPE(MetStateAdj), INTENT(IN)    :: State_Met_Adj   ! Meteorology State object
      TYPE(GrdStateAdj), INTENT(IN)    :: State_Grid_Adj  ! Grid State object
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(ChmStateAdj), INTENT(INOUT) :: State_Chm_Adj   ! Chemistry State object

! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
      ! Scalars
      INTEGER            :: NA, nAdvect
      INTEGER            :: I, J, TS_DYN
      REAL(fp)           :: AREA_M2, DT
      !arrays
      REAL(fp)           :: F(State_Grid_Adj%NZ,State_Chm_Adj%nAdvect)



      ! initialization
      TS_DYN     = GET_TS_DYN()                           ! Dyn timestep [sec]
      DT         = DBLE( TS_DYN )                         ! Dyn timestep [sec]
!      IF ( GET_ELAPSED_SEC() == 3000) THEN
      !adjoint validation
!      OPEN(UNIT = 1995, FILE='CONV_ADJOINT', ACCESS='APPEND')
!      ENDIF
      !end
      ! Number of advected species
      nAdvect = State_Chm_Adj%nAdvect
!$OMP PARALLEL DO &
!$OMP DEFAULT( SHARED ) &
!$OMP PRIVATE( NA, J,      I,      AREA_M2, RC, F) &
!$OMP SCHEDULE( DYNAMIC ) 
      DO NA = 1, nAdvect

      DO J = 1, State_Grid_Adj%NY
      DO I = 1, State_Grid_Adj%NX
         !tangent validation
!         State_Chm_Adj%SpeciesAdj = 0.0
!         State_Chm_Adj%SpeciesAdj(I, J, 1, 1) = 1.0
         !tangent end

         ! F equal to zero for co2 
         F = 0.0
         ! Grid box surface area [m2]
         AREA_M2 =  State_Grid_Adj%Area_M2(I,J)

         !--------------------------
         ! Do the cloud convection
         !--------------------------
         CALL DO_CLOUD_CONVECTION_ADJ( State_Chm_Adj  = State_Chm_Adj,              &
                                   State_Grid_Adj = State_Grid_Adj,             &
                                   State_Met_Adj  = State_Met_Adj,              &
                                   I              = I,                          &
                                   J              = J,                          &
                                   AREA_M2        = AREA_M2,                    &
                                   F              = F,                          &
                                   TS_DYN         = DT,                         &
                                      RC          = RC          )
!      IF ( I >40 .AND. I < 70 .AND. J>30 .AND. J<60 .and.&
!           GET_ELAPSED_SEC() == 3000) THEN
!          WRITE(1995, *) State_Chm_Adj%SpeciesAdj(I, J, 1, 1)
!      ENDIF
      ENDDO
      ENDDO
   
      ENDDO !nAdvect
!$OMP END PARALLEL DO
      !adjoint
!      IF ( GET_ELAPSED_SEC() == 3000) THEN
!      CLOSE(1995)
!      ENDIF
      END SUBROUTINE DO_CONVECTION_ADJ
                          
      SUBROUTINE DO_CONVECTION_D( State_Chm_Adj, State_Grid_Adj,     &
                                State_Met_Adj, RC )


      USE TIME_MOD,        ONLY : GET_TS_DYN
      USE TIME_MOD,        ONLY : GET_ELAPSED_SEC
      USE State_Chm_Adj_Mod,      ONLY : ChmStateAdj
!      USE State_Diag_Mod,     ONLY : DgnState
      USE State_Grid_Adj_Mod,     ONLY : GrdStateAdj
      USE State_Met_Adj_Mod,      ONLY : MetStateAdj

      ! convection adjoint for co2
      TYPE(MetStateAdj), INTENT(IN)    :: State_Met_Adj   ! Meteorology State object
      TYPE(GrdStateAdj), INTENT(IN)    :: State_Grid_Adj  ! Grid State object
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(ChmStateAdj), INTENT(INOUT) :: State_Chm_Adj   ! Chemistry State object

! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
      ! Scalars
      INTEGER            :: NA, nAdvect
      INTEGER            :: I, J, TS_DYN
      REAL(fp)           :: AREA_M2, DT
      !arrays
      REAL(fp)           :: F(State_Grid_Adj%NZ,State_Chm_Adj%nAdvect)



      ! initialization
      TS_DYN     = GET_TS_DYN()                           ! Dyn timestep [sec]
      DT         = DBLE( TS_DYN )                         ! Dyn timestep [sec]
!      IF ( GET_ELAPSED_SEC() == 3000) THEN
      !tangent validation 
!      OPEN(UNIT = 1995, FILE='CONV_TANGENT_TE', ACCESS='APPEND')
!      OPEN(UNIT = 1996, FILE='CONV_TANGENT_FTE', ACCESS='APPEND')
!      ENDIF
      !end
      ! Number of advected species
      nAdvect = State_Chm_Adj%nAdvect
!$OMP PARALLEL DO &
!$OMP DEFAULT( SHARED ) &
!$OMP PRIVATE( NA, J,      I,      AREA_M2, RC, F) &
!$OMP SCHEDULE( DYNAMIC ) 
      DO NA = 1, nAdvect

      DO J = 1, State_Grid_Adj%NY
      DO I = 1, State_Grid_Adj%NX
         !tangent validation
!         State_Chm_Adj%SpeciesAdj = 0.0
!         State_Chm_Adj%SpeciesAdj(I, J, 1, 1) = 1.0
!         State_Chm_Adj%Species = 0.0
!         State_Chm_Adj%Species(I, J, 1, 1) = 1.0
         !tangent end

         ! F equal to zero for co2
         F = 0.0
         ! Grid box surface area [m2]
         AREA_M2 =  State_Grid_Adj%Area_M2(I,J)
         !--------------------------
         ! Do the cloud convection
         !--------------------------
         CALL DO_CLOUD_CONVECTION_D( State_Chm_Adj  = State_Chm_Adj,              &
                                   State_Grid_Adj = State_Grid_Adj,             &
                                   State_Met_Adj  = State_Met_Adj,              &
                                   I              = I,                          &
                                   J              = J,                          &
                                   AREA_M2        = AREA_M2,                    &
                                   F              = F,                          &
                                   TS_DYN         = DT,                         &
                                      RC          = RC          )
!      IF ( I >40 .AND. I < 70 .AND. J>30 .AND. J<60 .and.&
!           GET_ELAPSED_SEC() == 3000) THEN
!         WRITE(1995, *) State_Chm_Adj%SpeciesAdj(I, J, 1, 1)
!         WRITE(1996, *) State_Chm_Adj%Species(I, J, 1, 1)
!      ENDIF
      ENDDO
      ENDDO

      ENDDO !nAdvect
!$OMP END PARALLEL DO
      !tangent
!      IF ( GET_ELAPSED_SEC() == 3000) THEN
!      CLOSE(1995)
!      ENDIF
      END SUBROUTINE DO_CONVECTION_D
     
      SUBROUTINE DO_CLOUD_CONVECTION_D( State_Chm_Adj,                        &
                                      State_Grid_Adj,                           &
                                      State_Met_Adj,                            &
                                      I,                                        &
                                      J,                                        &
                                      AREA_M2,                                  &
                                      F,                                        &
                                      TS_DYN,                                   &
                                      RC          )


      USE ErrCode_Adj_Mod
!      USE ERROR_MOD,          ONLY : IT_IS_NAN
!      USE ERROR_MOD,          ONLY : IT_IS_FINITE
!      USE Input_Opt_Mod,      ONLY : OptInput
      USE PhysConstants
      USE State_Chm_Adj_Mod,      ONLY : ChmStateAdj
!      USE State_Diag_Mod,     ONLY : DgnState
      USE State_Grid_Adj_Mod,     ONLY : GrdStateAdj
      USE State_Met_Adj_Mod,      ONLY : MetStateAdj
      USE TIME_MOD,        ONLY : GET_ELAPSED_SEC


      TYPE(ChmStateAdj), INTENT(INOUT) :: State_Chm_Adj  ! Grid State object
      TYPE(GrdStateAdj), INTENT(IN)    :: State_Grid_Adj  ! Grid State object
      TYPE(MetStateAdj), INTENT(IN)    :: State_Met_Adj   ! Meteorology State object
      INTEGER,        INTENT(IN)    :: I, J        ! Lon & lat indices
      REAL(fp),       INTENT(IN)    :: AREA_M2     ! Surface area [m2]
      REAL(fp),       INTENT(IN)    :: F(:,:)      ! Fraction of soluble species
                                                   !  for updraft scavenging
                                                   !  [unitless].  Computed by
                                                   !  routine  COMPUTE_F.
      REAL(fp),       INTENT(IN)    :: TS_DYN      ! Dynamic timestep [sec]
!      LOGICAL,        INTENT(IN)    :: USE_DIAG14  ! Archive DIAG14?
!      LOGICAL,        INTENT(IN)    :: USE_DIAG38  ! Archive DIAG38?
!
! !INPUT/OUTPUT PARAMETERS:
!
!      TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
!      REAL(fp),       INTENT(OUT)   :: DIAG14(:,:) ! Array for ND14 diagnostic
!      REAL(fp),       INTENT(OUT)   :: DIAG38(:,:) ! Array for ND38 diagnostic
      INTEGER,        INTENT(OUT)   :: RC          ! Return code

      REAL(fp), PARAMETER    :: TINYNUM = 1e-14_fp
      ! Scalars
      INTEGER                :: IC,          ISTEP,     K
      INTEGER                :: KTOP,        NC,        NDT
      INTEGER                :: NLAY,        NS,        CLDBASE
      INTEGER                :: NA,          nAdvect, NW
      REAL(fp)               :: CMFMC_BELOW, ALPHA,     ALPHA2
      REAL(fp)               :: CMOUT
      REAL(fp)               :: DNS,         ENTRN
      REAL(fp)               :: SDT,         T0_SUM
      REAL(fp)               :: QDOWN,       DT
      REAL(fp)               :: MB,          DELP_DRY_NUM
      REAL(fp)               :: temp
      REAL(fp)               :: temp0
      ! Strings
      CHARACTER(LEN=255)     :: ErrMsg, ThisLoc
      ! Arrays
      REAL(fp)               :: BMASS    (State_Grid_Adj%NZ)
      REAL(fp)               :: PDOWN    (State_Grid_Adj%NZ)

      ! Pointers
!      REAL(fp),      POINTER :: BXHEIGHT (:        )
      REAL(fp),      POINTER :: CMFMC    (:        )
      REAL(fp),      POINTER :: DQRCU    (:        )
      REAL(fp),      POINTER :: DTRAIN   (:        )
      REAL(fp),      POINTER :: PFICU    (:        )
      REAL(fp),      POINTER :: PFLCU    (:        )
!      REAL(fp),      POINTER :: REEVAPCN (:        )
      REAL(fp),      POINTER :: DELP_DRY (:        )
!      REAL(fp),      POINTER :: T        (:        )
!      REAL(fp),      POINTER :: H2O2s    (:        )
!      REAL(fp),      POINTER :: SO2s     (:        )
      REAL(fp),      POINTER :: Q        (:,:      )
      REAL(fp),      POINTER :: QD       (:,:      )
      REAL(fp),      POINTER :: QDELQ    (:,:      )
      REAL(fp),      POINTER :: QDELQ1   (:,:      )

      !tangent
      REAL(fp)               :: qc, qb_num, qb, qc_pres, tsum
      REAL(fp)               :: delq
      REAL(fp)               :: t1, t2, t3, t4

      REAL(fp)               :: delqd, qcd, qc_presd
      REAL(fp)               :: t1d, t2d, t3d, t4d
      REAL(fp)               :: tsumd, qbd, qb_numd

      rc = gc_success
      errmsg = ''
      ! Point to columns of derived-type object fields

      CMFMC    => State_Met_Adj%CMFMC   (I,J,2:State_Grid_Adj%NZ+1) ! Cloud mass flux
                                                            ! [kg/m2/s]
      DQRCU    => State_Met_Adj%DQRCU   (I,J,:        ) ! Precip production rate:
      DTRAIN   => State_Met_Adj%DTRAIN  (I,J,:        ) ! Detrainment flux [kg/m2/s]
!      REEVAPCN => State_Met%REEVAPCN(I,J,:        ) ! Evap of precip'ing conv.
      DELP_DRY => State_Met_Adj%DELP_DRY(I,J,:        ) ! Edge dry P diff [hPa]
!      T        => State_Met_Adj%T       (I,J,:        ) ! Air temperature [K]
!      H2O2s    => State_Chm%H2O2AfterChem(I,J,:   ) ! H2O2s from sulfate_mod
!      SO2s     => State_Chm%SO2AfterChem (I,J,:   ) ! SO2s from sulfate_mod
      qd        => State_Chm_Adj%SpeciesAdj (I,J,:,:      ) ! Chemical species
      Q        => State_Chm_Adj%Species(I, J, :, :)
      ! PFICU and PFLCU are on level edges
      PFICU    => State_Met_Adj%PFICU   (I,J,2:State_Grid_Adj%NZ+1) ! Dwnwd flx of conv
                                                            !  ice precip
                                                            !  [kg/m2/s]
      PFLCU    => State_Met_Adj%PFLCU   (I,J,2:State_Grid_Adj%NZ+1) ! Dwnwd flux of conv
                                                            !  liquid precip
      QDELQ    => State_Chm_Adj%QDELQ   (I,J, :, :)
      QDELQ1   => State_Chm_Adj%QDELQ1  (I,J, :, :)
      ! # of levels and # of species
      NLAY     = State_Grid_Adj%NZ
      NC       = State_Chm_Adj%nAdvect
      ! Top level for convection
      KTOP     = NLAY - 1

      ! Convection timestep [s]
      NDT      = TS_DYN

      ! Internal time step for convective mixing is 300 sec.
      ! Doug Rotman (LLNL) says that 450 sec works just as well.
      NS       = NDT / 300                ! Num internal timesteps (int)
      NS       = MAX( NS, 1 )             ! Set lower bound to 1
      DNS      = DBLE( NS )               ! Num internal timesteps (real)
      SDT      = DBLE( NDT ) / DBLE( NS ) ! seconds in internal timestep

      ! Minimum value of cloud base is the surface level
      CLDBASE = 1
      ! Find the cloud base
      DO K = 1, NLAY
         IF ( DQRCU(K) > 0e+0_fp ) THEN
            CLDBASE = K
            EXIT
         ENDIF
      ENDDO
      !-----------------------------------------------------------------
      ! Compute PDOWN and BMASS
      !-----------------------------------------------------------------
      DO K = 1, NLAY

         ! PDOWN is the convective precipitation leaving each
         ! box [cm3 H2O/cm2 air/s]. This accounts for the
         ! contribution from both liquid & ice precip.
         ! PFLCU and PFICU are converted from kg/m2/s to m3/m2/s
         ! using water and ice densities, respectively.
         ! m3/m2/s becomes cm3/cm2/s using a factor of 100.
         PDOWN(K) = ( ( PFLCU(K) / 1000e+0_fp )                  &
                 +   ( PFICU(K) /  917e+0_fp ) ) * 100e+0_fp

         ! BMASS is the dry air mass per unit area for the grid box
         ! bounded by level K and K+1 [kg/m2]
         ! BMASS is equivalent to deltaP (dry) * 100 / g
         ! This is done to keep BMASS in the same units as CMFMC * SDT
         BMASS(K) = DELP_DRY(K) * G0_100

      ENDDO
      MB = 0e+0_fp
      DO K = 1, CLDBASE-1
         MB = MB + BMASS(K)
      ENDDO

!      IF ( I >40 .AND. I < 70 .AND. J>30 .AND. J<60 .and. &
!           GET_ELAPSED_SEC() == 3000) THEN
!            print*, 'MB, CMFMC, DQRCU, DTRAIN, DELP_DRY, I, J', &
!         MB, SUM(CMFMC), SUM(DQRCU), SUM(DTRAIN),SUM(DELP_DRY), I, J
!      ENDIF

      DO na=1,nc
         ic = state_chm_adj%nadvect
         DO istep=1,ns
            qc = 0e+0_fp
            t0_sum = 0e+0_fp
            IF (cldbase .GT. 1) THEN
!
               IF (cmfmc(cldbase-1) .GT. tinynum) THEN
                  qb_num = 0e+0_fp
                  delp_dry_num = 0e+0_fp
                  qb_numd = 0.0
                  DO k=1,cldbase-1
                     qb_numd = qb_numd + delp_dry(k)*qd(k, ic)
                     qb_num = qb_num + q(k, ic)*delp_dry(k)
                     delp_dry_num = delp_dry_num + delp_dry(k)
                  END DO
                  qbd = qb_numd/delp_dry_num
                  qb = qb_num/delp_dry_num

                  temp = mb + cmfmc(cldbase-1)*sdt
                  temp0 = cmfmc(cldbase-1)*sdt
                  qcd = (mb*qbd+temp0*qd(cldbase, ic))/temp
                  qc = (mb*qb+temp0*q(cldbase, ic))/temp
                  qd(1:cldbase-1, ic) = qcd
                  q(1:cldbase-1, ic) = qc
               ELSE
                  qcd = qd(cldbase, ic)
                  qc = q(cldbase, ic)
               END IF
           ELSE

                  qcd = qd(cldbase, ic)
                  qc = q(cldbase, ic)

           END IF

           DO k=cldbase,ktop  
              cmout = 0e+0_fp
              entrn = 0e+0_fp
              qc_pres = 0e+0_fp

              IF (k .EQ. 1) THEN
                 cmfmc_below = 0e+0_fp
              ELSE
                 cmfmc_below = cmfmc(k-1)
              END IF
              IF (cmfmc_below .GT. tinynum) THEN
              ! Air mass flowing into cloud at grid box (K) [kg/m2/s]
              cmout = cmfmc(k) + dtrain(k)
!
              ! Amount of QC preserved against scavenging [kg/kg]
              entrn = cmout - cmfmc_below
!
!                  QC_PRES = QC * ( 1e+0_fp - F(K,IC) )
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Update QC taking entrainment into account [kg/kg]
! Prevent div by zero condition
              qc_presd = (1e+0_fp-f(k, na))*qcd
              qc_pres = qc*(1e+0_fp-f(k, na))
              IF (entrn .GE. 0e+0_fp .AND. cmout .GT. 0e+0_fp) THEN
                 qcd = (cmfmc_below*qc_presd+entrn*qd(k, ic))/cmout
                 qc = (cmfmc_below*qc_pres+entrn*q(k, ic))/cmout
              END IF
!
!
              t1d = cmfmc_below*qc_presd
              t1 = cmfmc_below*qc_pres
              t2d = -(cmfmc(k)*qcd)
              t2 = -(cmfmc(k)*qc)
              t3d = cmfmc(k)*qd(k+1, ic)
              t3 = cmfmc(k)*q(k+1, ic)
              t4d = -(cmfmc_below*qd(k, ic))
              t4 = -(cmfmc_below*q(k, ic))
              tsumd = t1d + t2d + t3d + t4d
              tsum = t1 + t2 + t3 + t4
!
              ! change in [kg/kg]
              ! If DELQ > Q then do not make Q negative!!!
              delqd = sdt*tsumd/bmass(k)
              delq = sdt/bmass(k)*tsum
            !huoxiao validation 
!            IF ( I == 48   .AND. J==59 .and.  &
!     &                GET_ELAPSED_SEC() == 3000) THEN
!            print*, 'QC_PRES, QC, Q(K+1), Q(K),  CMFMC_BELOW, CMFMC(K),' &
!                    //'I, J, K', QC_PRES, QC, &
!                    Q(K+1, IC), Q(K, IC),  CMFMC_BELOW, CMFMC(K), I, J, K, CLDBASE
!            print*, 'T1, T2, T3, T4', T1, T2, T3, T4
!            ENDIF
            !end

              IF (q(k, ic) + delq .LT. 0) THEN
                 delqd = -qd(k, ic)
                 delq = -q(k, ic)
              END IF
              qd(k, ic) = qd(k, ic) + delqd
              q(k, ic) = q(k, ic) + delq
 
           ELSE
              qcd = qd(k, ic)
              qc = q(k, ic)
!
              IF (cmfmc(k) .GT. tinynum) THEN
! [kg/m2/s * kg species/kg dry air]
                 t2d = -(cmfmc(k)*qcd)
                 t2 = -(cmfmc(k)*qc)
!
! Change in species concentration [kg/kg]
                 t3d = cmfmc(k)*qd(k+1, ic)
                 t3 = cmfmc(k)*q(k+1, ic)
!
! If DELQ > Q then do not make Q negative!!!
                 delqd = sdt*(t2d+t3d)/bmass(k)
                 delq = sdt/bmass(k)*(t2+t3)
! Add change in species to Q array [kg/kg]
!
              IF (q(k, ic) + delq .LT. 0.0e+0_fp) THEN
                 delqd = -qd(k, ic)
                 delq = -q(k, ic)
!
              END IF
!
              qd(k, ic) = qd(k, ic) + delqd
              q(k, ic) = q(k, ic) + delq
!
            END IF
          END IF !cmfmc_below .GT. tinynum
      END DO

      END DO
      ENDDO
      END SUBROUTINE DO_CLOUD_CONVECTION_D



      SUBROUTINE DO_CLOUD_CONVECTION_ADJ( State_Chm_Adj,                        &
                                      State_Grid_Adj,                           &
                                      State_Met_Adj,                            &
                                      I,                                        &
                                      J,                                        &
                                      AREA_M2,                                  &
                                      F,                                        &
                                      TS_DYN,                                   &
                                      RC          )

!
! !USES:
!
!       USE CMN_SIZE_MOD_ADJ
!      USE DEPO_MERCURY_MOD,   ONLY : ADD_Hg2_SNOWPACK
!      USE DEPO_MERCURY_MOD,   ONLY : ADD_Hg2_WD
!      USE DEPO_MERCURY_MOD,   ONLY : ADD_HgP_WD
      USE ErrCode_Adj_Mod
!      USE ERROR_MOD,          ONLY : IT_IS_NAN
!      USE ERROR_MOD,          ONLY : IT_IS_FINITE
!      USE Input_Opt_Mod,      ONLY : OptInput
      USE PhysConstants
      USE State_Chm_Adj_Mod,      ONLY : ChmStateAdj
!      USE State_Diag_Mod,     ONLY : DgnState
      USE State_Grid_Adj_Mod,     ONLY : GrdStateAdj
      USE State_Met_Adj_Mod,      ONLY : MetStateAdj
!      USE Species_Mod,        ONLY : Species
!      USE WETSCAV_MOD,        ONLY : WASHOUT
!      USE WETSCAV_MOD,        ONLY : LS_K_RAIN
!     USE WETSCAV_MOD,        ONLY : LS_F_PRIME
!
! !INPUT PARAMETERS:
!
!       LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
!      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
      TYPE(ChmStateAdj), INTENT(INOUT) :: State_Chm_Adj  ! Grid State object
      TYPE(GrdStateAdj), INTENT(IN)    :: State_Grid_Adj  ! Grid State object
      TYPE(MetStateAdj), INTENT(IN)    :: State_Met_Adj   ! Meteorology State object
      INTEGER,        INTENT(IN)    :: I, J        ! Lon & lat indices
      REAL(fp),       INTENT(IN)    :: AREA_M2     ! Surface area [m2]
      REAL(fp),       INTENT(IN)    :: F(:,:)      ! Fraction of soluble species
                                                   !  for updraft scavenging
                                                   !  [unitless].  Computed by
                                                   !  routine  COMPUTE_F.
      REAL(fp),       INTENT(IN)    :: TS_DYN      ! Dynamic timestep [sec]
!      LOGICAL,        INTENT(IN)    :: USE_DIAG14  ! Archive DIAG14?
!      LOGICAL,        INTENT(IN)    :: USE_DIAG38  ! Archive DIAG38?
!
! !INPUT/OUTPUT PARAMETERS:
!
!      TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
!      REAL(fp),       INTENT(OUT)   :: DIAG14(:,:) ! Array for ND14 diagnostic
!      REAL(fp),       INTENT(OUT)   :: DIAG38(:,:) ! Array for ND38 diagnostic
      INTEGER,        INTENT(OUT)   :: RC          ! Return code
!
! !REMARKS:
!  Reference:
!  ============================================================================
!  Lin, SJ.  "Description of the parameterization of cumulus transport
!     in the 3D Goddard Chemistry Transport Model, NASA/GSFC, 1996.
!                                                                             .
!  Unit conversion for BMASS:
!
!      Ps - Pt (mb)| P2 - P1 | 100 Pa |  s^2  | 1  |  1 kg        kg
!     -------------+---------+--------+-------+----+--------  =  -----
!                  | Ps - Pt |   mb   | 9.8 m | Pa | m^2 s^2      m^2
!
!                                                                             .
!  NOTE: We are passing I & J down to this routine so that it can call the
!  proper code from "mercury_mod.f".  Normally, we wouldn't pass I & J as
!  arguments to columnized code.  This prevents rewriting the mercury_mod.f
!  routines ADD_Hg2_
!
! !REVISION HISTORY:
!  15 Jul 2009 - R. Yantosca - Columnized and cleaned up.
!                            - CLDMAS renamed to CMFMC and DTRN renamed
!                              to DTRAIN for consistency w/ GEOS-5.
!  17 Jul 2009 - R. Yantosca - Now do unit conversion of Q array from
!                              [kg] --> [v/v] and vice versa internally
!  14 Dec 2009 - R. Yantosca - Now remove internal unit conversion, since
!                              Q now comes in as [mol/mol] (=[v/v]) from the
!                              calling routine.
!  14 Dec 2009 - R. Yantosca - Remove COEF from the argument list
!  06 May 2010 - R. Yantosca - Now add IDENT via the argument list
!  29 Sep 2010 - R. Yantosca - Modified for MERRA met fields
!  05 Oct 2010 - R. Yantosca - Now pass COEF via the argument list
!  05 Oct 2010 - R. Yantosca - Attach ND14 and ND38 diagnostics
!  15 Oct 2010 - H. Amos     - Added BXHEIGHT and T as arguments
!  15 Oct 2010 - R. Yantosca - Added I, J, H2O2s and SO2s as arguments
!  15 Oct 2010 - H. Amos     - Added scavenging below cloud base
!  06 Apr 2011 - M.Fu, H.Amos- Bug fix: make sure washout adheres to the same
!                              algorithm as in the wet deposition code.
!  27 Jul 2011 - R. Yantosca - Declare CLDBASE as INTEGER to avoid PGI errors
!  16 Aug 2011 - J. Fisher   - Bug fix: use IS_Hg2() and IS_HgP to test if
!                              a tracer is Hg2 or HgP (for tagged species)
!  16 Aug 2011 - J. Fisher   - Now use WETLOSS instead of T0_SUM in the ND38
!                              diagnostic below the cloud.  Using T0_SUM leads
!                              us to over-count the tracer scavenged out of
!                              the column.
!  22 Oct 2012 - R. Yantosca - Now reference Headers/gigc_errcode_mod.F90
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  31 May 2013 - R. Yantosca - Now pass State_Chm to WASHOUT
!  05 Sep 2013 - R. Yantosca - Bug fix: DT is apparently undefined, but still
!                              passed to WASHOUT.  Use SDT instead.  This
!                              avoids a floating-point error.
!  18 Apr 2014 - R. Yantosca - Now point to 3-D arrays internally
!  18 Apr 2014 - R. Yantosca - Now also pass N_TRACERS (to facilitate APM)
!  18 Apr 2014 - R. Yantosca - Remove code that we don't need anymore
!  04 Feb 2015 - M. Sulprizio- Fix calculation of WETLOSS for non-aerosol
!                              tracers below the cloud base (C. Friedman)
!  20 Apr 2015 - E. Lundgren - Use DELP*100/g instead of AD/area as grid box
!                              moist mass per area and remove AD from routine
!  20 May 2015 - M. Sulprizio- Apply bug fixes provided by Viral Shah:
!                              -- Remove F(K,IC) > 0 condition that prevents
!                                 archiving of deposited mass in DIAG38
!                              -- Add statement that subtracts the wet deposited
!                                 amount from the atmospheric mass
!                              -- Fix inconsistency in units when T0_SUM is used
!  04 Jun 2015 - E. Lundgren - Adapt Viral Shah bug fixes to moist mixing ratio
!  09 Jun 2015 - R. Yantosca - Now deposit Hg2, HgP to snowpack regardless of
!                              whether the dynamic ocean is used
!  15 Jun 2015 - E. Lundgren - Now use kg/kg total air as tracer units not v/v
!  22 Jun 2015 - E. Lundgren - Move QB_NUM calculation to within timestep loop
!  12 Aug 2015 - R. Yantosca - Treat MERRA2 in same way as we do for GEOS-FP
!  14 Sep 2015 - E. Lundgren - Apply bug fixes provided by Viral Shah:
!                              -- Prevent ALPHA > 1 in washout of aerosols
!                              -- Add tracer GAINED to Q before WETLOSS
!                                 calculation in aerosol washout
!  22 Apr 2016 - R. Yantosca - Now get Is_Hg2 & Is_HgP from species database
!  25 Apr 2016 - R. Yantosca - Now pass Hg category # to ADD_Hg2_* functions
!  28 Apr 2016 - R. Yantosca - Rewrite Is_Hg block to avoid unassociated
!                              pointer seg faults
!  29 Apr 2016 - R. Yantosca - Don't initialize pointers in declaration stmts
!  30 Jun 2016 - R. Yantosca - Remove instances of STT.  Now get the advected
!                              species ID from State_Chm%Map_Advect.
!  01 Jul 2016 - R. Yantosca - Now rename species DB object ThisSpc to SpcInfo
!  05 Jul 2016 - R. Yantosca - Now replace N_TRACERS argument with the
!                              State_Chm%nAdvect field
!  06 Jul 2016 - E. Lundgren - Now use kg/kg dry air as spc units, requiring
!                              use of DELP_DRY instead of DELP and PEDGE_DRY
!                              for avg mixing ratio
!  07 Jul 2016 - R. Yantosca - Bug fix: F and ISOL need to be indexed with
!                              the advected species index NA
!  07 Jul 2016 - R. Yantosca - DIAG14 and DIAG38 now are dimensioned with
!                              of size State_Chm%nWetDep instead of N_TRACERS
!  20 Sep 2016 - R. Yantosca - Rewrote IF test to avoid Gfortran compiler error
!  24 Aug 2017 - M. Sulprizio- Rename routine from DO_MERRA_CONVECTION to
!                              DO_CLOUD_CONVECTION
!  06 Nov 2017 - R. Yantosca - Bug fix: ND14 should use 1..nAdvect species
!  09 Nov 2017 - R. Yantosca - Return error condition to calling program
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      REAL(fp), PARAMETER    :: TINYNUM = 1e-14_fp
!
! !LOCAL VARIABLES:
!
      ! Scalars
      INTEGER                :: IC,          ISTEP,     K
      INTEGER                :: KTOP,        NC,        NDT
      INTEGER                :: NLAY,        NS,        CLDBASE
      INTEGER                :: NA,          nAdvect, NW
      REAL(fp)               :: CMFMC_BELOW, ALPHA,     ALPHA2
      REAL(fp)               :: CMOUT,       ADDELQ,    ADDQ
      REAL(fp)               :: DNS,         ENTRN,     ADQC
      REAL(fp)               :: ADQC_PRES,     ADQC_SCAV,   SDT
      REAL(fp)               :: ADT0,          T0_SUM,    ADT1
      REAL(fp)               :: ADT2,          ADT3,        ADT4
      REAL(fp)               :: ADTSUM,        LOST,      GAINED
      REAL(fp)               :: QDOWN,       DT
      REAL(fp)               :: MB,          ADQB
      REAL(fp)               :: ADQB_NUM,      DELP_DRY_NUM

      ! Strings
      CHARACTER(LEN=255)     :: ErrMsg, ThisLoc

      ! Arrays
      REAL(fp)               :: BMASS    (State_Grid_Adj%NZ)
      REAL(fp)               :: PDOWN    (State_Grid_Adj%NZ)

      ! Pointers
!      REAL(fp),      POINTER :: BXHEIGHT (:        )
      REAL(fp),      POINTER :: CMFMC    (:        )
      REAL(fp),      POINTER :: DQRCU    (:        )
      REAL(fp),      POINTER :: DTRAIN   (:        )
      REAL(fp),      POINTER :: PFICU    (:        )
      REAL(fp),      POINTER :: PFLCU    (:        )
!      REAL(fp),      POINTER :: REEVAPCN (:        )
      REAL(fp),      POINTER :: DELP_DRY (:        )
!      REAL(fp),      POINTER :: T        (:        )
!      REAL(fp),      POINTER :: H2O2s    (:        )
!      REAL(fp),      POINTER :: SO2s     (:        )
      REAL(fp),      POINTER :: ADQ        (:,:      )
      REAL(fp),      POINTER :: Q        (:,:      )
      REAL(fp),      POINTER :: QDELQ    (:,:      )
      REAL(fp),      POINTER :: QDELQ1   (:,:      )
!      TYPE(Species), POINTER :: SpcInfo     

     !========================================================================
      ! (0)  I n i t i a l i z a t i o n
      !========================================================================

      ! Initialize
      RC       = GC_SUCCESS
      ErrMsg   = ''
      ThisLoc  =     &
      ' -> at Do_Merra_Convection (in module GeosCore/convection_mod.F)'
      ! Point to columns of derived-type object fields
!      BXHEIGHT => State_Met%BXHEIGHT(I,J,:        ) ! Box height [m]
      CMFMC    => State_Met_Adj%CMFMC   (I,J,2:State_Grid_Adj%NZ+1) ! Cloud mass flux
                                                            ! [kg/m2/s]
      DQRCU    => State_Met_Adj%DQRCU   (I,J,:        ) ! Precip production rate:
      DTRAIN   => State_Met_Adj%DTRAIN  (I,J,:        ) ! Detrainment flux [kg/m2/s]
!      REEVAPCN => State_Met%REEVAPCN(I,J,:        ) ! Evap of precip'ing conv.
      DELP_DRY => State_Met_Adj%DELP_DRY(I,J,:        ) ! Edge dry P diff [hPa]
!      T        => State_Met_Adj%T       (I,J,:        ) ! Air temperature [K]
!      H2O2s    => State_Chm%H2O2AfterChem(I,J,:   ) ! H2O2s from sulfate_mod
!      SO2s     => State_Chm%SO2AfterChem (I,J,:   ) ! SO2s from sulfate_mod
      ADQ        => State_Chm_Adj%SpeciesAdj (I,J,:,:      ) ! Chemical species
      Q        => State_Chm_Adj%Species(I, J, :, :)
!      SpcInfo  => NULL()                            ! Species database entry

      ! PFICU and PFLCU are on level edges
      PFICU    => State_Met_Adj%PFICU   (I,J,2:State_Grid_Adj%NZ+1) ! Dwnwd flx of conv
                                                            !  ice precip
                                                            !  [kg/m2/s]
      PFLCU    => State_Met_Adj%PFLCU   (I,J,2:State_Grid_Adj%NZ+1) ! Dwnwd flux of conv
                                                            !  liquid precip
      QDELQ    => State_Chm_Adj%QDELQ   (I,J, :, :)
      QDELQ1   => State_Chm_Adj%QDELQ1  (I,J, :, :)
      ! # of levels and # of species
      NLAY     = State_Grid_Adj%NZ
      NC       = State_Chm_Adj%nAdvect

      ! Top level for convection
      KTOP     = NLAY - 1

      ! Convection timestep [s]
      NDT      = TS_DYN

      ! Internal time step for convective mixing is 300 sec.
      ! Doug Rotman (LLNL) says that 450 sec works just as well.
      NS       = NDT / 300                ! Num internal timesteps (int)
      NS       = MAX( NS, 1 )             ! Set lower bound to 1
      DNS      = DBLE( NS )               ! Num internal timesteps (real)
      SDT      = DBLE( NDT ) / DBLE( NS ) ! seconds in internal timestep

      !-----------------------------------------------------------------
      ! Determine location of the cloud base, which is the level where
      ! we start to have non-zero convective precipitation formation
      !-----------------------------------------------------------------

      ! Minimum value of cloud base is the surface level
      CLDBASE = 1
      ! Find the cloud base
      DO K = 1, NLAY
         IF ( DQRCU(K) > 0e+0_fp ) THEN
            CLDBASE = K
            EXIT
         ENDIF
      ENDDO
      !-----------------------------------------------------------------
      ! Compute PDOWN and BMASS
      !-----------------------------------------------------------------
      DO K = 1, NLAY

         ! PDOWN is the convective precipitation leaving each
         ! box [cm3 H2O/cm2 air/s]. This accounts for the
         ! contribution from both liquid & ice precip.
         ! PFLCU and PFICU are converted from kg/m2/s to m3/m2/s
         ! using water and ice densities, respectively.
         ! m3/m2/s becomes cm3/cm2/s using a factor of 100.
         PDOWN(K) = ( ( PFLCU(K) / 1000e+0_fp )                  &
                 +   ( PFICU(K) /  917e+0_fp ) ) * 100e+0_fp

         ! BMASS is the dry air mass per unit area for the grid box
         ! bounded by level K and K+1 [kg/m2]
         ! BMASS is equivalent to deltaP (dry) * 100 / g
         ! This is done to keep BMASS in the same units as CMFMC * SDT
         BMASS(K) = DELP_DRY(K) * G0_100

      ENDDO

      !-----------------------------------------------------------------
      ! Compute MB, the mass per unit area of dry air below the cloud
      ! base [kg/m2]. Calculate MB by looping over levels below the
      ! cloud base.
      !-----------------------------------------------------------------
      MB = 0e+0_fp
      DO K = 1, CLDBASE-1
         MB = MB + BMASS(K)
      ENDDO

      !========================================================================
      ! (1)  A d v e c t e d   S p e c i e s   L o o p
      !========================================================================

      ! Loop over only the advected species
      DO NA = 1, NC

         ! Get the species ID (modelID) from the advected species ID
         IC       =  State_Chm_Adj%nAdvect
         ! Look up the corresponding entry in the species database
!         SpcInfo  => State_Chm%SpcData(IC)%Info

         ! Also get the corresponding wetdep ID
!         NW       =  SpcInfo%WetDepId

         ! Zero the DIAG14 diagnostic array
!         DIAG14(:,NA) = 0.0_fp

         ! Zero the DIAG38 diagnostic array
!         IF ( NW > 0 ) THEN
!            DIAG38(:,NW) = 0.0_fp
!         ENDIF

         !=====================================================================
         ! (2)  I n t e r n a l   T i m e   S t e p   L o o p
         !=====================================================================
         ADQC_PRES = 0.0
         ADQC_SCAV = 0.0
         ADQC = 0.0
         ADT1 = 0.0
         ADT2 = 0.0
         ADT3 = 0.0
         ADT4 = 0.0
         ADTSUM = 0.0
         ADDELQ = 0.0
         ADQB = 0.0
         ADQB_NUM = 0.0
         DO ISTEP = NS, 1, -1
            ! Initialize
            ADQC = 0.0
            !----------------------------------------------------------
            ! B e l o w   C l o u d   B a s e   (K < CLDBASE)
            !
            ! QB is the "weighted avg" mixing ratio below the cloud
            ! base [kg/kg dry air].
            ! QC is the mixing ratio of the air that moved in cumulus
            ! transport up to the next level [kg/kg dry air].
            ! MB is the dry mass of air below the cloud base per
            ! unit area [kg/m2] (see calculation before loop).
            !-----------------------------------------------------------

            ! We need to make this a nested IF statement so that we don't
            ! get an out-of-bounds error when CLDBASE=1 (bmy, 11/18/10)
          
            !==================================================================
            ! (3)  A b o v e   C l o u d   B a s e
            !==================================================================
            DO K = KTOP, CLDBASE, -1
               ! Initialize
               ALPHA   = 0e+0_fp
               ALPHA2  = 0e+0_fp
               CMOUT   = 0e+0_fp
               ENTRN   = 0e+0_fp

               ! CMFMC_BELOW is the air mass [kg/m2/s] coming into the
               ! grid box (K) from the box immediately below (K-1).
               IF ( K == 1 ) THEN
                  CMFMC_BELOW = 0e+0_fp
               ELSE
                  CMFMC_BELOW = CMFMC(K-1)
               ENDIF

               IF ( CMFMC_BELOW > TINYNUM ) THEN

                  CMOUT   = CMFMC(K) + DTRAIN(K)
                  ENTRN   = CMOUT - CMFMC_BELOW
                  ! forward code 
                  !             Q(K,IC) = Q(K,IC) + DELQ    
                  ADDELQ = ADDELQ+ADQ(K, IC)

                  IF ( QDELQ(K, ISTEP ) == 1 ) THEN
                     print*, 'HUOXIAO_DEBUG I'
                     !ADQ(K, IC) = 0.0
                     ADQ(K, IC) = -ADDELQ+ADQ(K, IC)
                  ENDIF
                  ADTSUM = ADTSUM+ADDELQ*(SDT/BMASS(K))
                  ADDELQ = 0.

                  ! forward code
                  !             TSUM    = T1 + T2 + T3 + T4
                  ADT1 = ADT1+ADTSUM
                  ADT2 = ADT2+ADTSUM
                  ADT3 = ADT3+ADTSUM
                  ADT4 = ADT4+ADTSUM
                  ADTSUM = 0.

                  ! forward code
                  !           T0      =  CMFMC_BELOW * QC_SCAV
                  !           T1      =  CMFMC_BELOW * QC_PRES
                  !           T2      = -CMFMC(K  )  * QC
                  !           T3      =  CMFMC(K  )  * Q(K+1,IC)
                  !           T4      = -CMFMC_BELOW * Q(K,  IC)         
                  !ADQ(K, IC) = 0.0 
                  ADQ(K, IC) = ADQ(K, IC)-ADT4*CMFMC_BELOW
                  ADT4 = 0.

                  !ADQ(K+1, IC) = 0.0
                  ADQ(K+1, IC) = ADQ(K+1, IC)+ADT3*CMFMC(K)
                  ADT3 = 0.

                  !ADQC = 0.0 
                  ADQC = ADQC-ADT2*CMFMC(K)
                  ADT2 = 0. 

                  !ADQC_PRES = 0.0 
                  ADQC_PRES = ADQC_PRES+ADT1*CMFMC_BELOW
                  ADT1 = 0.  

                  IF ( ENTRN >= 0e+0_fp .and. CMOUT > 0e+0_fp ) THEN
                       !ADQ(K, IC) = 0.0
                       ADQ(K, IC) = ADQ(K, IC)+ADQC*(ENTRN/CMOUT)
                       !ADQC_PRES = 0.0
                       ADQC_PRES = ADQC_PRES+ADQC*(CMFMC_BELOW/CMOUT)
                       ADQC = 0.
                  ENDIF   
                  ! forward code 
                  !            QC_PRES = QC * ( 1e+0_fp - F(K,NA) ) for co2 adjoint F=0
                  !ADQC = 0.0
                  ADQC = ADQC+ADQC_PRES*(1.D0-F(K, IC))
                  ADQC_PRES = 0.

                ELSE
              
                  IF ( CMFMC(K) > TINYNUM ) THEN
                  ! forward code 
                  !            Q(K,IC) = Q(K,IC) + DELQ  WARNING
                     ADDELQ = ADQ(K, IC)+ADDELQ
                  
                  IF ( QDELQ1(K, ISTEP) == 1 ) THEN
                  !   ADQ(K,IC) = 0.0
                     ADQ = -ADDELQ+ADQ(K,IC)
                  ENDIF
                  ! forward code
                  !            DELQ = ( SDT / BMASS(I,J,K) ) * (T2 + T3)
                     ADT3 = ( SDT / BMASS(K) ) * ADDELQ+ADT3
                     ADT2 = ( SDT / BMASS(K) ) * ADDELQ+ADT2
                     ADDELQ = 0.
                  ! forward code 
                  !            T2   = -CMFMC(K) * QC
                  !            T3   = T3   =  CMFMC(K) * Q(K+1,IC)
                  !   ADQ(K+1, IC) = 0.0
                     ADQ(K+1, IC) = ADQ(K+1, IC) + CMFMC(K) * ADT3
                     ADT3 = 0.

                  !   ADQC = 0.0
                     ADQC = ADQC - CMFMC(K) * ADT2
                     ADT2 = 0.

                  ENDIF

                  ! forward code
                  !             QC = Q(K,IC)
                  !ADQ(K ,IC) = 0.0
                  ADQ(K, IC) = ADQ(K, IC)+ADQC          
                  ADQC = 0.

                ENDIF
            
            ENDDO


            !debug 
            IF ( CLDBASE > 1 ) THEN
               IF ( CMFMC(CLDBASE-1) > TINYNUM ) THEN 
                  
                  DELP_DRY_NUM = 0.0
                  DO K  = 1, CLDBASE-1
                     DELP_DRY_NUM = DELP_DRY_NUM + DELP_DRY(K)
                  ENDDO     
          
                  DO K = 1, CLDBASE-1
                     ADQC = ADQC+ADQ(K, IC)
                     ADQ(K, IC) = 0.
                  ENDDO


                  ! forward code 
                  !                               QC = ( MB*QB + CMFMC(CLDBASE-1) *
                  !            &                           Q(CLDBASE,IC)    * SDT  ) /
                  !            &                 ( MB    + CMFMC(CLDBASE-1) * SDT  )
                !  ADQ(CLDBASE, IC) = 0.0
                  ADQ(CLDBASE, IC) = ADQ(CLDBASE, IC)+ADQC*   &
                        (CMFMC(CLDBASE-1)*SDT/(MB+CMFMC(CLDBASE-1)*SDT))
                  ADQB = ADQB+ADQC*(MB/(MB+CMFMC(CLDBASE-1)*SDT))
                  ADQC = 0.0       
                  ! forward code 
                  !              QB = QB_NUM / DELP_DRY_NUM
                  ADQB_NUM = ADQB_NUM+ADQB/DELP_DRY_NUM
                  ADQB = 0.0

                  DO K  = 1, CLDBASE-1
                  !   ADQ(K, IC) = 0.0
                     ADQ(K, IC) = ADQ(K, IC) + DELP_DRY(K)*ADQB_NUM
                  ENDDO 
               
                  !forward code QB_NUM  = 0e+0_fp
                  ADQB_NUM = 0.0
 
                 ELSE 
                 ! ADQ(CLDBASE, IC) = 0.0 
                  ADQ(CLDBASE, IC) = ADQC+ADQ(CLDBASE, IC)
                  ADQC = 0.0
                
                 ENDIF
 
             ELSE
                 ! ADQ(CLDBASE, IC) = 0.0 
                  ADQ(CLDBASE, IC) = ADQC+ADQ(CLDBASE, IC)
                  ADQC = 0.0
                 
             ENDIF 
            !debug 
         ENDDO               ! End internal timestep loop

         ! Free pointer
      ENDDO                  ! End loop over advected species

      !================================================================
      ! Succesful return!
      !================================================================

      ! Nullify pointers
      NULLIFY( CMFMC    )
      NULLIFY( DQRCU    )
      NULLIFY( DTRAIN   )
      NULLIFY( PFICU    )
      NULLIFY( PFLCU    )
      NULLIFY( DELP_DRY )
      NULLIFY( ADQ        )

      ! Set error code to success
      RC                      = GC_SUCCESS

      END SUBROUTINE DO_CLOUD_CONVECTION_ADJ
!      END
      END MODULE CONVECTION_ADJ_MOD


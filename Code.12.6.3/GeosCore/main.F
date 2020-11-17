      PROGRAM INVERSE
      
#     include  "../Adjoint/define_adj.h" 

      USE GEOS_CHEM_MOD,        ONLY : DO_GEOS_CHEM
      USE NCDF_MOD
      USE m_netcdf_io_checks,   ONLY : Ncdoes_Var_Exist
      USE m_netcdf_io_read,     ONLY : NcRd
      USE GEOS_CHEM_ADJ_MOD,    ONLY : DO_GEOS_CHEM_ADJ
      USE State_Chm_Adj_Mod,    ONLY : EMS_SF
!      USE ADJ_ARRAYS_MOD,       ONLY : ITE_NUM, LLPAR
!      USE ADJ_ARRAYS_MOD,       ONLY : INIT_ADJ_ARRAYS
      USE ADJ_ARRAYS_MOD,       ONLY : N_CALC, N_CALC_STOP
      USE Adj_Arrays_Mod,       ONLY : NOPT
      USE Adj_Arrays_Mod,       ONLY : OBS_EXIST
!      USE ADJ_ARRAYS_MOD,       ONLY : ITE_NUM_STOP
      USE ADJ_ARRAYS_MOD,       ONLY : COST_FUNC
      USE ADJ_ARRAYS_MOD,       ONLY : OUT_ITER
!      USE ADJ_ARRAYS_MOD,       ONLY : ITE_NUM, LLPAR
!      USE ADJ_ARRAYS_MOD,       ONLY : SET_TIME_FOR_ADJ
!      USE ADJ_ARRAYS_MOD,       ONLY : EMS_SF_ADJ
!      USE ADJ_ARRAYS_MOD,       ONLY : STT_ADJ_FD, COST_FUNC_SAV
!      USE ADJ_ARRAYS_MOD,       ONLY : IFD,JFD,LFD,NFD
!      USE ADJ_ARRAYS_MOD,       ONLY : INIT_CF_REGION
!      USE ADJ_ARRAYS_MOD,       ONLY : CLEAN_ADJ_ARRAYS
!      USE INVERSE_MOD,          ONLY : GET_GRADNT_FROM_ADJ 
      USE File_Ope_Mod,          ONLY : READ_SF_FILE
      USE File_Ope_Mod,          ONLY : READ_GDT_FILE
      USE File_Ope_Mod,          ONLY : READ_CFN_FILE
      USE File_Ope_Mod,          ONLY : MAKE_SF_FILE
      USE File_Ope_Mod,          ONLY : MAKE_GDT_FILE
      USE File_Ope_Mod,          ONLY : MAKE_CFN_FILE
      USE File_Ope_Mod,          ONLY : MAKE_CFNJo_FILE
      USE File_Ope_Mod,          ONLY : GET_GRADNT_FROM_ADJ
      USE File_Ope_Mod,          ONLY : WRITE_SF_TO_RSTFILE
      USE File_Ope_Mod,          ONLY : WRITE_SF_TO_SPCFILE
      USE File_Ope_Mod,          ONLY : WRITE_EXIT_FLAG
      USE File_Ope_Mod,          ONLY : READ_TOT_ITER
      USE File_Ope_Mod,          ONLY : READ_OUT_ITER
      USE File_Ope_Mod,          ONLY : WRITE_OBS_PERTD
      USE File_Ope_Mod,          ONLY : GET_SPECIES0
      USE INVERSE_MOD,           ONLY : GET_X_FROM_SF
      USE INVERSE_MOD,           ONLY : GET_SF_FROM_X
      USE INVERSE_MOD,           ONLY : INIT_INVERSE 
!      USE INVERSE_MOD,          ONLY : DISPLAY_STUFF 
      USE INVERSE_MOD,           ONLY : SET_SF
      USE INVERSE_MOD,           ONLY : SET_LOG_SF
      USE INVERSE_MOD,           ONLY : GRADNT
      USE INVERSE_MOD,           ONLY : X
      USE INVERSE_MOD,           ONLY : SET_OPT_RANGE
      USE INVERSE_MOD,           ONLY : UPDATE_X0
!      USE INVERSE_MOD,          ONLY : CLEANUP_INVERSE
!      USE INVERSE_MOD,          ONLY : SET_SF_FORFD
      USE LOGICAL_ADJ_MOD,       ONLY : LADJ, LADJ_EMS, LICS
      USE LOGICAL_ADJ_MOD,       ONLY : LICS_INC
      USE State_Grid_Adj_Mod,        Only : Init_State_Grid_Adj
      USE State_Grid_Adj_Mod,        Only : Allocate_State_Grid_Adj
      USE State_Grid_Adj_Mod,        Only : Calculate_Grid_Parameter
      USE State_Grid_Adj_Mod,        Only : Cleanup_State_Grid_Adj
      USE State_Grid_Adj_Mod,        Only : GrdStateAdj
      USE State_Chm_Adj_Mod,         Only : Init_State_Chm_Adj
      USE State_Chm_Adj_Mod,         Only : Cleanup_State_Chm_Adj
      USE State_Chm_Adj_Mod,         Only : ChmStateAdj
      USE State_Chm_Adj_Mod,         Only : EMS_SF_ADJ
      USE State_Met_Adj_Mod,         Only : Init_State_Met_Adj
      USE State_Met_Adj_Mod,         Only : Cleanup_State_Met_Adj
      USE State_Met_Adj_Mod,         Only : MetStateAdj

      USE Parameter_Adj_Mod,         Only : Init_Parameter
      USE Parameter_Adj_Mod,         Only : Cleanup_Parameter

      USE Input_Adj_Mod,             Only : Init_Input_Adj_Mod
      USE Input_Adj_Mod,             Only : OptInputAdj
      USE LICS_INC_Mod,              ONLY : DO_LICS_INC
      USE Time_Mod,                  Only : Get_Elapsed_Sec
      USE Time_Mod,                  Only : Set_Elapsed_Sec_Adj
      USE Time_Mod,                  Only : Set_Current_Time

      USE BEC_MOD,                   Only : BinvB
      USE FILE_OPE_MOD,              ONLY : GET_BEC
      ! Reference to f90 modules
!     USE LBFGSB

      IMPLICIT NONE

      TYPE(ChmStateAdj)                   :: State_Chm_Adj
      TYPE(GrdStateAdj)                   :: State_Grid_Adj
      TYPE(MetStateAdj)                   :: State_Met_Adj

      TYPE(OptInputAdj)                   :: Input_Opt_Adj

      ! Program variables
      LOGICAL                  :: ITERATE = .TRUE.

      INTEGER                  :: RC
      ! Variables and parameters for optimization -- see setulb.f for
      ! definitions of these.
      INTEGER                  ::   iprint, isave(44)
      CHARACTER*60             ::   task, csave
      DOUBLE PRECISION         ::   factr, pgtol, dsave(29)
      DOUBLE PRECISION         ::   F
      LOGICAL                  ::   lsave(4)
      INTEGER, PARAMETER       ::   MMAX = 17
      INTEGER                  ::   LENWA
      INTEGER, ALLOCATABLE     ::   nbd(:)
      INTEGER, ALLOCATABLE     ::   iwa(:)
      REAL*8,  ALLOCATABLE     ::   llim(:)
      REAL*8,  ALLOCATABLE     ::   u(:)
      REAL*8,  ALLOCATABLE     ::   wa(:)
      INTEGER                  ::   IOPT,MM
 
      !local
      INTEGER                  ::   ILI, TOT_TIME_SEC
      INTEGER                  ::   FGEOS, MAXITE
      REAL*8                   ::   FCOST_FUNC, LCOST_FUNC
      !time variable
      REAL*8                              :: START_TIME, END_TIME
      CHARACTER(LEN=255)                  :: DATE, TIME, ZONE
      INTEGER                             :: VALUES(8)

!variable parameters
      !time length
      call date_and_time( date , time, zone, values )
      print*, 'PROGRAM START TIME ', date(1:4)//'-'//date(5:6)//
     &        '-'//date(7:8)//'  '//
     &         time(1:2)//':'//time(3:4)//':'//time(5:10)
       ! initialization 
      CALL Init_State_Grid_Adj( State_Grid_Adj, RC )
      CALL Allocate_State_Grid_Adj( State_Grid_Adj, RC )
      CALL Init_State_Chm_Adj(State_Chm_Adj, State_Grid_Adj, RC )
      CALL Init_State_Met_Adj( State_Grid_Adj, State_Met_Adj, RC )
      CALL Calculate_Grid_Parameter( State_Grid_Adj, RC )
      CALL Init_Parameter(State_Chm_Adj)
      CALL Init_Input_Adj_Mod( Input_Opt_Adj )
      CALL Write_Exit_Flag(0)

!      CALL GET_BEC ( State_Chm_Adj, State_Grid_Adj )
!      CALL BinvB ( State_Chm_Adj%Sx, State_Chm_Adj%Sy, 
!     &      State_Chm_Adj%Sz, State_Chm_Adj%Sx_inv, 
!     &      State_Chm_Adj%Sy_inv, State_Chm_Adj%Sz_inv)
!      STOP
!      ! time length
!      time_str1d(1) = 1
!      time_cnt1d(1) = 1
!      varname = 'length'
!      CALL NcRd(time_varrd_1di, FID, TRIM(varname), time_str1d, 
!     &                              time_cnt1d, err_stop, stat)

!      NOCO2 = time_varrd_1di(1)
!      print*, 'NOCO2', NOCO2

!      end do
      !read flag for geos-chem FGEOS=1 only running geos-chem
      CALL READ_FLAG_GEOS(FGEOS)
      ! iteration number 

      CALL READ_ITE_NUM(N_CALC_STOP)

      WRITE(*,*) 'N_CALC_STOP', N_CALC_STOP

      OBS_EXIST = 0
      FCOST_FUNC = 0.0
      MAXITE = 10
      ! Number of previous iterations used in the Hessian approximation
      ! Here it was set to the total number of iteration
      ! =>Memory problem for big control vector so set back to MM = 5
      MM        = 5
      task      = 'START'
      IPRINT    = 1
      FACTR     = 1.0D01
      PGTOL     = 1.0D-05
      N_CALC    = 0

      ! Initialize arrays for optimization
      IF (LADJ)  THEN

      CALL INIT_SETULB( State_Grid_Adj )

      CALL INIT_INVERSE

!      CALL INIT_ADJ_ARRAYS
 
!      CALL INIT_CF_REGION
  
      ENDIF

      DO IOPT = 1, NOPT
         NBD(IOPT)    = 1
        LLIM(IOPT)   = 0.0D0
      ENDDO  

      ! Mare array of scaling factors into a vector for optimization
      IF ( ITERATE) THEN 

      CALL GET_X_FROM_SF( State_Chm_Adj, State_Grid_Adj )

      ENDIF 

!       Beginning of the loop
      IF ( .NOT. LICS_INC .AND. FGEOS==0 ) THEN
      DO WHILE( ITERATE )


!            IF( N_CALC== 1) THEN
!              print*, 'huoxiao_debug', GRADNT
!            ENDIF
          print*, 'main before max gdt, max x',MAXVAL(GRADNT), MAXVAL(X)
             ! Call the L-BFGS-code
            
            CALL SETULB( NOPT,   MM,     X,      LLIM,   U,     NBD, 
     &		      F,      GRADNT, FACTR,  PGTOL,  WA,    IWA,  
     &                TASK,   IPRINT, CSAVE,  LSAVE,  ISAVE, DSAVE )
 
          print*, 'main after max gdt, max x',MAXVAL(GRADNT), MAXVAL(X)
 
          print*, ' done setulb ' , TASK
       
         IF ( TASK(1:2) == 'FG' ) THEN

            N_CALC = N_CALC+1
            COST_FUNC = 0.0
            
            IF ( N_CALC < N_CALC_STOP ) THEN 
                
!               write(*,*) 'read file'
!               ! Read scaling factors
               CALL READ_SF_FILE( State_Grid_Adj, State_Chm_Adj )
               !get x from sf
               CALL GET_X_FROM_SF( State_Chm_Adj, State_Grid_Adj )
!               ! read gdt from the file 
               CALL READ_GDT_FILE( State_Chm_Adj, State_Grid_Adj )
!               ! Read cost function from disk
               CALL READ_CFN_FILE
!               ! Put adjoints into GRADNT vector
               CALL GET_GRADNT_FROM_ADJ( State_Chm_Adj, State_Grid_Adj )
!               print*, 'huoxiao_debug', GRADNT              
!               ! finite difference test
!               IF ( LFDTEST .AND. LADJ_EMS ) THEN
    
!               STT_ADJ_FD(N_CALC) = EMS_SF_ADJ(IFD,JFD,LFD,NFD)

!               END IF
               ! Copy value of COST_FUNC to the optimization variable F
                F = COST_FUNC

!               ! Save current cost function for finite difference test
!               COST_FUNC_SAV(N_CALC) = COST_FUNC


            ELSEIF( N_CALC == N_CALC_STOP ) THEN
               ! reset emission scale factor adjoint
               State_Chm_Adj%SpeciesAdj = 0.0
               EMS_SF_ADJ = 0.0

!               IF ( LFDTEST .AND. N_CALC > 1 ) THEN

!                  CALL SET_SF_FORFD

               IF ( N_CALC == 1 ) THEN 
#if    defined ( LOG_OPT )
               CALL SET_LOG_SF
#else
               CALL SET_SF
#endif
               ELSE           
               ! Update the scaling factors to the current X
               CALL GET_SF_FROM_X( State_Grid_Adj, State_Chm_Adj )
!               print*, 'MAXVAL GET HUOXIAO_DEBUG', 
!     &               MAXVAL(State_Chm_Adj%Species(:, :, :, :))

               IF ( LICS  ) THEN
               CALL WRITE_SF_TO_RSTFILE(State_Chm_Adj, State_Grid_Adj)
               CALL WRITE_SF_TO_SPCFILE(State_Chm_Adj, State_Grid_Adj)
               ENDIF
               ENDIF
               IF ( N_CALC == 1 ) THEN
                  CALL GET_SPECIES0( State_Chm_Adj, State_Grid_Adj )
                  State_Chm_Adj%SPECIES = State_Chm_Adj%SPECIES0
!                  print*, 'MAIN MAX', 
!     &            MAXVAL(State_Chm_Adj%SPECIES(:, :, :, :))
               ENDIF
               ! Write scaling factors
               CALL MAKE_SF_FILE( State_Grid_Adj, State_Chm_Adj )
                CALL DO_GEOS_CHEM
!               CALL SET_TIME_FOR_ADJ

!               CALL INIT_TPCORE_ADJ( LLPAR ) 

               ! finite difference test
!               IF ( .not. LFDTEST ) THEN

                CALL DO_GEOS_CHEM_ADJ(          State_Chm_Adj, 
     &                                          State_Grid_Adj, 
     &                                          State_Met_Adj, 
     &                                          Input_Opt_Adj )
               ! For finite difference test, we may or may not do adjoint
!               ELSE

!                   CALL MAYBE_DO_GEOS_CHEM_ADJ

!               ENDIF
!               WRITE(*, *) 'COST_FUNCTION', COST_FUNC
!               CALL SET_OPT_RANGE( State_Chm_Adj, State_Grid_Adj )

               ! Write gradients
               CALL MAKE_GDT_FILE( State_Chm_Adj, State_Grid_Adj )
              
!               ! Write scaling factors
!               CALL MAKE_SF_FILE( State_Grid_Adj, State_Chm_Adj )

               ! Write cost function
               CALL MAKE_CFN_FILE

                !cost function jo
               CALL MAKE_CFNJo_FILE

           

               ITERATE = .FALSE.
 
               ! Write results to screen
!               print*, 'DISPLAY_STUFF'
!               CALL DISPLAY_STUFF( 1 )
!             ENDIF
            ENDIF  ! N_CALC < N_CALC_STOP

       

!        ELSEIF ( TASK(1:5) == 'NEW_X' ) THEN
!            CONTINUE
    
!        ELSEIF ( TASK(1:4) == 'STOP' ) THEN
!             Exit loop
!            WRITE(*,*) 'force quit'
!            STOP
     
!       ENDIF
       ELSEIF( TASK(1:5) == 'NEW_X' ) THEN
   
         PRINT*, 'now task is:', TASK

       ELSE 
         CALL Write_Exit_Flag(1)
         print*, 'iterate exit task is:', TASK
         task = 'STOP' 
         ITERATE = .FALSE.
       ENDIF  !TASK(1:2) == 'FG' 

      ENDDO
      ENDIF

      IF ( LICS_INC .AND. FGEOS==0 ) THEN

         CALL READ_OUT_ITER
         print*, 'OUTER LOOP IS ', OUT_ITER 
         N_CALC = 0
         CALL READ_TOT_ITER
         DO ILI = 1, MAXITE
            print*, ' INCREMENT 4DVAR ASSIMILATION, '
     &               //'CURRENT ITERATION IS ', ILI

            ! Call the L-BFGS-code
            CALL SETULB( NOPT,   MM,     X,      LLIM,   U,     NBD,
     &                 F,      GRADNT, FACTR,  PGTOL,  WA,    IWA,
     &                 TASK,   IPRINT, CSAVE,  LSAVE,  ISAVE, DSAVE )
            
            COST_FUNC = 0.0
            !set 0
            State_Chm_Adj%SpeciesAdj = 0.0
 
            IF ( TASK(1:2) == 'FG' ) THEN
               N_CALC = N_CALC+1
               IF ( ILI /= 1 ) THEN
               ! Update the scaling factors to the current X
               CALL GET_SF_FROM_X( State_Grid_Adj, State_Chm_Adj )
               ENDIF
               ! forward for 0th iteration
               IF ( ILI == 1 ) THEN
               
                  CALL DO_GEOS_CHEM 
                  TOT_TIME_SEC = Get_Elapsed_Sec()
               ENDIF
               
               CALL Set_Elapsed_Sec_Adj( 0, .TRUE., 
     &               TOT_TIME_SEC, .TRUE.)
               CALL Set_Current_Time
               
               CALL DO_LICS_INC( State_Chm_Adj,
     &                           State_Grid_Adj,
     &                           State_Met_Adj,
     &                           Input_Opt_Adj )
               ! Write gradients
               CALL MAKE_GDT_FILE( State_Chm_Adj, State_Grid_Adj )

               ! Write scaling factors
               CALL MAKE_SF_FILE( State_Grid_Adj, State_Chm_Adj )

               ! Write cost function
               CALL MAKE_CFN_FILE

                !cost function jo
               CALL MAKE_CFNJo_FILE

               !write N_CALC
!               CALL WRITE_N_CALC
!               !get x from sf
               CALL GET_X_FROM_SF( State_Chm_Adj, State_Grid_Adj )
               ! get gradient from adj
               CALL GET_GRADNT_FROM_ADJ( State_Chm_Adj, State_Grid_Adj )
               F = COST_FUNC
               
               IF( (ABS(COST_FUNC-FCOST_FUNC)<0.0001) .OR. 
     &             (ILI==MAXITE) ) THEN
                CALL UPDATE_X0( State_Chm_Adj, State_Grid_Adj )
                CALL WRITE_SF_TO_RSTFILE(State_Chm_Adj, State_Grid_Adj)
                CALL WRITE_SF_TO_SPCFILE(State_Chm_Adj, State_Grid_Adj)
                CALL WRITE_OBS_PERTD
                 EXIT
               ENDIF
               FCOST_FUNC = COST_FUNC
!               ! Read cost function from disk
!            CALL READ_CFN_FILE
          ELSEIF( TASK(1:5) == 'NEW_X' ) THEN

               IF ( ILI==MAXITE ) THEN
                 CALL UPDATE_X0( State_Chm_Adj, State_Grid_Adj )
                 CALL WRITE_SF_TO_RSTFILE(State_Chm_Adj, State_Grid_Adj)
                 CALL WRITE_SF_TO_SPCFILE(State_Chm_Adj, State_Grid_Adj)
                 CALL WRITE_OBS_PERTD
               ENDIF
               PRINT*, 'now task is:', TASK

          ELSE
               CALL UPDATE_X0( State_Chm_Adj, State_Grid_Adj )
               CALL WRITE_SF_TO_RSTFILE(State_Chm_Adj, State_Grid_Adj)
               CALL WRITE_SF_TO_SPCFILE(State_Chm_Adj, State_Grid_Adj)
               CALL WRITE_OBS_PERTD
               CALL Write_Exit_Flag(1)
               print*, 'iterate exit task is:', TASK
!               CALL DO_GEOS_CHEM               
               EXIT
         !ITERATE = .FALSE.
          ENDIF

        ENDDO
      ENDIF
      IF ( FGEOS == 1 ) THEN
         CALL DO_GEOS_CHEM

      ENDIF
       call date_and_time( date , time, zone, values )
      print*, 'PROGRAM START TIME ', date(1:4)//'-'//date(5:6)//
     &         '-'//date(7:8)//'  '//
     &        time(1:2)//':'//time(3:4)//':'//time(5:10) 
      ! clean up
      CALL Cleanup_State_Grid_Adj( State_Grid_Adj, RC )
      CALL Cleanup_State_Chm_Adj( State_Chm_Adj, RC )
      CALL Cleanup_State_Met_Adj( State_Met_Adj, RC )
      CALL Cleanup_Parameter

        
      CONTAINS
!------------------------------------------------------------------------------

      SUBROUTINE INIT_SETULB( State_Grid_Adj )

!******************************************************************************
!  Subroutine INIT_SETULB initializes arrays used by the optimization routine,
!  setulb, whose size depends upon the model simulation type and resolution.
!  (dkh, 06/07/09)
!
!  NOTES:
!
!******************************************************************************
!
!      USE ERROR_MOD,          ONLY : ALLOC_ERR
      USE INVERSE_MOD,        ONLY : CALC_NOPT
      USE State_Grid_Adj_Mod, Only : GrdStateAdj

      TYPE(GrdStateAdj),INTENT(IN) :: State_Grid_Adj

      ! Local variables
      INTEGER                :: AS

      !=================================================================
      ! INIT_SETULB begins here!
      !=================================================================

      ! Calculate the maximum number of control parameters that could
      ! be optimized, NOPT
      CALL CALC_NOPT( State_Grid_Adj )
      LENWA = 2 * MM * NOPT + 5 * NOPT + 11 * MM * MM + 8 * MM
      ALLOCATE( NBD( NOPT ), STAT=AS )

      ALLOCATE( IWA( 3*NOPT ), STAT=AS )
      IWA = 0

      ALLOCATE( LLIM( NOPT ), STAT=AS )

      ALLOCATE( U( NOPT ), STAT=AS )

      ALLOCATE( WA( LENWA ), STAT=AS )
      WA = 0

      ! Return to calling program
      END SUBROUTINE INIT_SETULB
!------------------------------------------------------------------------------

      SUBROUTINE READ_ITE_NUM(num)
!
!******************************************************************************
!
!  Read iteration number
!
!******************************************************************************

      INTEGER, INTENT(INOUT)       :: num
!      INTEGER, INTENT(INOUT)       :: ifd
      INTEGER                      :: TEMP
      LOGICAL                      :: ALIVE
      CHARACTER(LEN=60)            :: FILENAME
    
      FILENAME = './ITER'   
      ! if the file not exist, num = -1, return
      INQUIRE( FILE = FILENAME, EXIST = ALIVE)
      IF ( .not. ALIVE ) THEN
     
      num = -1
     
      return
      ENDIF
       
      ! open the iter file to read iteration number
      OPEN( UNIT = 1995, FILE = FILENAME)
      READ(1995,*) TEMP
      num = TEMP
      CLOSE(1995)

!      OPEN( UNIT = 1995, FILE = './IFD')
!      READ(1995,*) IFD
!      CLOSE(1995)

      END SUBROUTINE READ_ITE_NUM 


      SUBROUTINE READ_FLAG_GEOS(FGEOS)

      INTEGER,           INTENT(INOUT)  :: FGEOS
      OPEN(UNIT=1995, FILE='FGEOS')
      READ(1995, *) FGEOS
      CLOSE(1995)
     
      END SUBROUTINE READ_FLAG_GEOS
!      SUBROUTINE READ_GDT_FILE

!
!******************************************************************************
!  Subroutine MAKE_GDT_FILE creates a binary file of ADJ_xxx
!  (dkh, 9/17/04)
!
!  Module Variable as Input:
!  ============================================================================
!  (1 ) N_CALC     : Current iteration number
!  (2 ) ICS_SF_ADJ : Array of adjoint gradients to be written
!  (3 ) EMS_SF_ADJ : Array of adjoint gradients to be written
!
!  NOTES:
!  (1 ) Just like MAKE_OBS_FILE except
!       - write to .adj. file
!  (2 ) Changed name to MAKE_GDT_FILE.  Now the .adj. files are trajectories,
!       and the .gdt. files are final gradients  (dkh, 10/03/04)
!  (3 ) Added ACTIVE_VARS == 'EMISSIONS' case.  (dkh, 11/27/04)
!  (4 ) Added ACTIVE_VARS == 'FDTEST' case. (dkh, 02/17/05)
!  (5 ) Now use CATEGORY = 'IJ-GDE-$' for 'EMISSIONS' case. (dkh, 03/29/05)
!  (6 ) No longer pass COST_FUNC in the header; use cnf.* files. (dkh, 02/13/06)
!  (7 ) Rename everything, replace CMN_ADJ, move nonessential stuff
!       to diagnostic files  (dkh, ks, mak, cs  06/07/09)
!  (8 ) Add normalized gradients IJ-GDEN$ (dkh, 05/06/10)
!  (9 ) Now support strat fluxes LADJ_STRAT (hml, dkh, 02/20/12, adj32_025)
!******************************************************************************
!
      ! References to F90 modules
!      USE ADJ_ARRAYS_MOD,      ONLY : EMS_SF_ADJ
!      USE ADJ_ARRAYS_MOD,      ONLY : NNEMS, MMSCL
!      USE ADJ_ARRAYS_MOD,      ONLY : IIPAR, JJPAR
!      USE LOGICAL_ADJ_MOD,     ONLY :  LADJ_EMS


      ! Added for reaction rate sensitivities (tww, 05/08/12)

      ! huoxiao adjoint 
!      USE ADJ_ARRAYS_MOD,      ONLY : N_CALC

      ! Local Variables
!      INTEGER              :: I, J , M, N
!      CHARACTER(LEN=255)   :: FILENAME


      ! Added for reaction rate sensitivity (tww, 05/08/12)

!      CHARACTER(LEN=20)    :: OUTPUT_GDT_FILE

      !huoxiao adjoint
!      CHARACTER(LEN=20)    :: CHA_ITE_NUM
      !=================================================================
      ! MAKE_GDT_FILE begins here!
      !=================================================================



      ! Hardwire output file for now
!      WRITE(CHA_ITE_NUM, '(I2.2)') N_CALC  
!      OUTPUT_GDT_FILE = 'gctm.gdt.'//TRIM( CHA_ITE_NUM )
      !=================================================================
      ! Open the adjoint file for output -- binary punch format
      !=================================================================
      ! Add the OPTDATA_DIR prefix to the file name
!      FILENAME = TRIM( OUTPUT_GDT_FILE )

      ! Open file for output

!      OPEN( UNIT = 1997, FILE=FILENAME )
!      IF ( LADJ_EMS ) THEN


         !=================================================================
         ! Write each observed quantity to the observation file
         !=================================================================
!         DO N = 1, NNEMS

            !Temporarily store quantities in the TRACER array

!            DO M = 1, MMSCL
!            DO J = 1, JJPAR
!            DO I = 1, IIPAR
!               READ(1997,*)  EMS_SF_ADJ(I,J,M,N)
!            ENDDO
!            ENDDO
!            ENDDO
!        ENDDO

!      ENDIF

      ! Close file
!      CLOSE( 1997 )

      ! Return to calling program
!      END SUBROUTINE READ_GDT_FILE

      SUBROUTINE WRITE_TEST


       USE NCDF_MOD,             ONLY : NC_OPEN, NC_VAR_WRITE
       USE m_netcdf_io_checks,   ONLY : Ncdoes_Var_Exist
       USE m_netcdf_io_write,     ONLY : NcWr

       USE m_netcdf_io_close,      ONLY : Nccl

      USE m_netcdf_io_read,       ONLY : NcRd


       CHARACTER(LEN=255)              FILENAME
      ! dimension
       INTEGER                          :: NX, NY, NZ, nSpecies, nTimes
       CHARACTER(LEN=255)               :: TIMESTR, varname
      !related nc file read
       INTEGER                          :: co2_str4d(4)
       INTEGER                          :: co2_cnt4d(4)
       INTEGER                          :: fID, rc
       REAL*4, POINTER                  :: co2(:, :, :, :)
      
       ALLOCATE(co2(144, 91, 47, 1), stat=rc)

       FILENAME = 'GEOSChem.AdjParameters.'//'20171001_000000'//'.nc4'
       print*, 'Write Species in ', FILENAME
       CALL NC_OPEN( TRIM(FILENAME), fID )
       co2_str4d(1) = 1
       co2_str4d(2) = 1
       co2_str4d(3) = 1
       co2_str4d(4) = 1

       co2_cnt4d(1) = State_Grid_Adj%NX
       co2_cnt4d(2) = State_Grid_Adj%NY
       co2_cnt4d(3) = State_Grid_Adj%NZ
       co2_cnt4d(4) = 1
       print*, 'test', State_Grid_Adj%NX
       varname = 'delp_dry'
       CALL Nc_Var_Write(FID, TRIM(varname), co2)
       END SUBROUTINE WRITE_TEST

      END PROGRAM INVERSE
     
      

      


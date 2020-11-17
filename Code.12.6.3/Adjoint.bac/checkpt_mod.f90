      MODULE CHECKPT_MOD

      USE PRECISION_MOD
      
!
!******************************************************************************
      IMPLICIT NONE


      CONTAINS

      SUBROUTINE LOAD_CHECKPT_DATA( State_Chm_Adj, State_Grid_Adj, State_Met_Adj )

       ! USE 
       USE State_Chm_Adj_Mod,         Only : ChmStateAdj
       USE State_Grid_Adj_Mod,        Only : GrdStateAdj
       USE State_Met_Adj_Mod,         Only : MetStateAdj
       USE TIME_MOD,                  Only : GET_NYMD
       USE TIME_MOD,                  Only : GET_NHMS
       TYPE(GrdStateAdj), INTENT(IN)       :: State_Grid_Adj
       TYPE(ChmStateAdj), INTENT(INOUT)       :: State_Chm_Adj
       TYPE(MetStateAdj), INTENT(INOUT)    :: State_Met_Adj
       INTEGER                             :: NYMD, NHMS
!******************************************************************************
!  Subroutine LOAD_CHECKPT_DATA reads in information stored during the forward
!   calculation.  Some of the data (CSPEC) needs to be rotated.
!******************************************************************************
               ! Read data from file
       NYMD = GET_NYMD()
       NHMS = GET_NHMS()
       CALL READ_CHECKPT_FILE ( NYMD, NHMS, State_Chm_Adj, State_Grid_Adj, State_Met_Adj )

      END SUBROUTINE LOAD_CHECKPT_DATA
!------------------------------------------------------------------------------------------
      SUBROUTINE READ_CHECKPT_FILE( NYMD, NHMS, State_Chm_Adj, &
                                   State_Grid_Adj, State_Met_Adj )

       USE NCDF_MOD
       USE m_netcdf_io_checks,   ONLY : Ncdoes_Var_Exist
       USE m_netcdf_io_read,     ONLY : NcRd
       USE State_Chm_Adj_Mod,    ONLY : ChmStateAdj
       USE State_Grid_Adj_Mod,   ONLY : GrdStateAdj
       USE State_Met_Adj_Mod,   ONLY : MetStateAdj
       USE m_netcdf_io_close,      ONLY : Nccl
       USE Adj_Arrays_Mod,       Only : N_CALC
       USE Logical_Adj_Mod,      ONLY : LICS_INC
!       USE ADJ_ARRAYS_MOD,     ONLY:   IGLOB
!       USE ADJ_ARRAYS_MOD,     ONLY:   JGLOB
!       USE ADJ_ARRAYS_MOD,     ONLY:   LGLOB
!       USE ADJ_ARRAYS_MOD,     ONLY:   CHK_STT
!
!******************************************************************************
!  Subroutine READ_CHECKPT_FILE initializes GEOS-CHEM tracer concentrations
!  from a checkpoint file (binary punch file format)
!  (dkh, 8/30/04)
!
!  Arguments as input:
!  ============================================================================
       INTEGER, INTENT(IN)      ::    NYMD, NHMS
       TYPE(GrdStateAdj), INTENT(IN)    :: State_Grid_Adj  ! Grid State object
       TYPE(ChmStateAdj), INTENT(INOUT) :: State_Chm_Adj   ! Chemistry State object
       TYPE(MetStateAdj), INTENT(INOUT) :: State_Met_Adj   
       INTEGER                        HOUR, MINUTE, SEC
       CHARACTER(LEN=255)              FILENAME
       CHARACTER(LEN=20)              TEMP_HOUR
       CHARACTER(LEN=20)              TEMP_MINUTE
       CHARACTER(LEN=20)              TEMP_SEC
       CHARACTER(LEN=20)              TEMP_NYMD
       CHARACTER(LEN=20)              TEMP_TIME

      ! dimension 
       INTEGER                          :: NX, NY, NZ, nSpecies, nTimes
      !related nc file read
       REAL*8, ALLOCATABLE              :: co2_varrd_4dr(:, :, :, :)
       INTEGER                          :: co2_str4d(4)
       INTEGER                          :: co2_cnt4d(4)
       !surface pressure
       REAL*8, ALLOCATABLE              :: sp_varrd_3dr(:, :, :)
       INTEGER                          :: sp_str3d(3)
       INTEGER                          :: sp_cnt3d(3)
       !specific humudity
       REAL*8, ALLOCATABLE              :: sh_varrd_4dr(:, :, :, :)
       INTEGER                          :: sh_str4d(4)
       INTEGER                          :: sh_cnt4d(4)

       INTEGER                          :: FID
      
       LOGICAL                          :: err_stop
       INTEGER                          :: stat
       INTEGER                          :: RC
       CHARACTER(LEN=255)               :: varname

        
       WRITE(*,*) "READ_CHECKPT_FILE"
       ! initialization 
       NX = State_Grid_Adj%NX
       NY = State_Grid_Adj%NY
       NZ = State_Grid_Adj%NZ
       nTimes = 1

       HOUR = NHMS/10000
       MINUTE = (NHMS-HOUR*10000)/1000
       SEC = NHMS-HOUR*10000-MINUTE*1000
!integer to string
       WRITE(TEMP_NYMD, *)  NYMD
       WRITE(TEMP_HOUR, "(I2.2)")  HOUR
       WRITE(TEMP_MINUTE, "(I2.2)") MINUTE*10
       WRITE(TEMP_SEC, "(I2.2)") SEC
       TEMP_TIME =        TRIM(ADJUSTL(TEMP_NYMD))  //     &
                          "_"                       //     &
                          TRIM(ADJUSTL(TEMP_HOUR))  //     &
                          TRIM(ADJUSTL(TEMP_MINUTE))//     &
                          "z.nc4"
       FILENAME = "./FdCo2Out/GEOSChem.SpeciesConc."//TEMP_TIME
!       status = getcwd( dirname )
!       WRITE(*,*) 'getcwd( dirname )', dirname
       CALL NC_OPEN( TRIM(FILENAME), fID )
       ALLOCATE( co2_varrd_4dr(NX, NY, NZ, nTimes), STAT = RC )
       co2_str4d(1) = 1
       co2_str4d(2) = 1
       co2_str4d(3) = 1
       co2_str4d(4) = 1

       co2_cnt4d(1) = State_Grid_Adj%NX
       co2_cnt4d(2) = State_Grid_Adj%NY
       co2_cnt4d(3) = State_Grid_Adj%NZ
       co2_cnt4d(4) = 1

       varname = 'SpeciesConc_CO2'
       CALL NcRd(co2_varrd_4dr, FID, TRIM(varname), co2_str4d, co2_cnt4d)
       !convert to ppmv
       IF ( LICS_INC ) THEN
          State_Chm_Adj%SpeciestB(:, :, :, :) = co2_varrd_4dr(:, :, :, :) 
          RETURN
       ELSE
          State_Chm_Adj%Species(:, :, :, :) = co2_varrd_4dr(:, :, :, :)
       ENDIF
       !background specie concentration
       CALL Nccl(fID)

       FILENAME = "./FdCo2Out/GEOSChem.StateMet."//TEMP_TIME

       CALL NC_OPEN( TRIM(FILENAME), fID )
       !surface pressure
       ALLOCATE( sp_varrd_3dr(NX, NY,  nTimes), STAT = RC )
       sp_str3d(1) = 1
       sp_str3d(2) = 1
       sp_str3d(3) = 1

       sp_cnt3d(1) = State_Grid_Adj%NX
       sp_cnt3d(2) = State_Grid_Adj%NY
       sp_cnt3d(3) = 1

       varname = 'Met_PSC2WET'
       CALL NcRd(sp_varrd_3dr, FID, TRIM(varname), sp_str3d, sp_cnt3d)
       State_Met_Adj%SURFACE_PRESSURE(:, :) = sp_varrd_3dr(:, :, 1)

       !specific humudity
       ALLOCATE( sh_varrd_4dr(NX, NY, NZ, nTimes), STAT = RC )
       sh_str4d(1) = 1
       sh_str4d(2) = 1
       sh_str4d(3) = 1
       sh_str4d(4) = 1

       sh_cnt4d(1) = State_Grid_Adj%NX
       sh_cnt4d(2) = State_Grid_Adj%NY
       sh_cnt4d(3) = State_Grid_Adj%NZ
       sh_cnt4d(4) = 1

       varname = 'Met_SPHU'
       CALL NcRd(sh_varrd_4dr, FID, TRIM(varname), sh_str4d, sh_cnt4d)
       !convert to ppmv 
       State_Met_Adj%SPECIFIC_HUMIDITY(:, :, :) = sh_varrd_4dr(:, :, :, 1)

       CALL Nccl( FID ) 
      
      IF ( ALLOCATED( co2_varrd_4dr ) ) DEALLOCATE( co2_varrd_4dr )
      IF ( ALLOCATED( sh_varrd_4dr ) ) DEALLOCATE( sh_varrd_4dr )
      IF ( ALLOCATED( sp_varrd_3dr ) ) DEALLOCATE( sp_varrd_3dr )
END SUBROUTINE READ_CHECKPT_FILE
      

      END MODULE CHECKPT_MOD

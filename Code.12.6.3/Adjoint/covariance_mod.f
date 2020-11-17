
      MODULE COVARIANCE_MOD

!=======================================================================
! Module COVARIANCE_MOD contains routines to perform calculations
! involving non-diagonal error covariance matrices.  (nb,yd,dkh, 02/11/13)
!
! Based on Singh et al. (2011):
! Construction of non-diagonal background error covariance matrices
! for global chemical data assimilation, Geophysical Model Development 4,299-316
!
!
!
!  Module Routines
!  ============================================================================
!  (1 ) CALC_COV_ERROR      : Computes prior term in cost function (x-xb)^TB^-1(x-xb)
!                             for non-diagonal covariance matrices
!  (2 ) BVECT               : Performs covariance matrix-vector operations
!
!=======================================================================

      IMPLICIT NONE

      ! Header files
!#     include "CMN_SIZE"          ! Size parameters
!#     include "CMN_DIAG"          ! Diagnostic switches, NJDAY
!#     include "CMN_GCTM"          ! Physical constants
#     include "define_adj.h"      ! Obs operators

      CONTAINS


!-----------------------------------------------------------------------

      SUBROUTINE CALC_COV_ERROR ( APCOST, State_Grid_Adj )

!***********************************************************************
! Subroutine CALC_COV_ERROR calculates the prior term of the
! cost function when the covariance matrix of prior errors
! contains off-diagonal terms (nb,yd, dkh, 02/11/13)
!
!  Module Variables as inputs:
!  ====================================================================
! (1 ) EMS_ERROR      (REAL*8) : Diag error in scaling factor [none]
! (2 ) COV_ERROR_LX   (REAL*8) : Correlation length in x [km]
! (3 ) COV_ERROR_LY   (REAL*8) : Correlation length in y [km]
!
!  Module Variables as output:
!  ====================================================================
! (1 ) EMS_SF_ADJ     (REAL*8) : Emissions scaling factor adjoint [J]
! (2 ) TEMP2          (REAL*8) : Temp array for adjoint forcing
!
! Notes:
!
!***********************************************************************
      !huoxiao
      USE State_Grid_Adj_Mod,        Only : GrdStateAdj
!      USE ERROR_MOD,       ONLY : ALLOC_ERR, ERROR_STOP
      USE ADJ_ARRAYS_MOD,  ONLY :  COST_FUNC !COST_ARRAY,    COST_FUNC
      USE State_Chm_Adj_Mod,   ONLY : EMS_SF0, EMS_SF
      USE ADJ_ARRAYS_MOD,  ONLY : EMS_ERROR,COV_ERROR_LX,COV_ERROR_LY
      USE ADJ_ARRAYS_MOD,  ONLY : IIPAR, JJPAR
!      USE ADJ_ARRAYS_MOD,  ONLY : REG_PARAM_EMS
!      USE ADJ_ARRAYS_MOD,  ONLY : MMSCL,         NNEMS
      USE State_Chm_Adj_MOd,   ONLY : EMS_SF_ADJ !, TEMP2
      USE LOGICAL_ADJ_MOD, ONLY : LADJ_EMS  !, LBKCOV
!      USE GRID_MOD       , ONLY : GET_XMID, GET_YMID

!      USE netcdf

!#     include "CMN_SIZE"
      TYPE(GrdStateAdj), INTENT(IN)    :: State_Grid_Adj
            
      !huoxiao
      REAL*8                   :: REG_PARAM_EMS(1)
      REAL*8                   :: TEMP2(144, 91, 1, 1)
      INTEGER                  :: MMSCL,         NNEMS
      !
      REAL*8                   :: S2_INV_2D(IIPAR,JJPAR)
      REAL*8                   :: REG_4D(IIPAR, JJPAR,1, 1)
      REAL*8                   :: S2_INV
      REAL*8                   :: REG
      REAL*8,  ALLOCATABLE     :: APCOST(:,:,:,:)
      REAL                     :: TEMP(IIPAR,JJPAR)
      INTEGER                  :: I, J, M, N, STATUS, NCID, VARID
      CHARACTER(255)           :: SCALEFN

      ! Off-diagonal terms of B

      REAL*8                   :: SIGMAX, SIGMAY,CORR_LX, CORR_LY
      REAL*8                   :: MATINVX(JJPAR,IIPAR,IIPAR)
      REAL*8                   :: MATY(JJPAR,JJPAR),MATXX(IIPAR,IIPAR)
      REAL*8                   :: B(IIPAR),C(JJPAR),UIN(IIPAR,JJPAR)
      REAL*8                   :: UOUT(IIPAR,JJPAR)
      INTEGER                  :: ISTATUS, LFLAG, INFO

      !huoxiao
      REAL*8                   :: TMP_XMID(1), TMP_YMID(1)


      WRITE(*,*) 'Starting Off-Diagonal Term Calculation'
      REG_PARAM_EMS = 0.01
      MMSCL = 1
      NNEMS = 1
      TMP_XMID(1) = State_Grid_Adj%XMID(IIPAR, 1)
      TMP_YMID(1) = State_Grid_Adj%YMID(1, JJPAR)
!!!$OMP PARALLEL DO
!!!$OMP+DEFAULT( SHARED )
!!!$OMP+PRIVATE( I, J, M, N, REG_4D, S2_INV_2D )
      DO N = 1, NNEMS

         ! Put CALL BVECT here since correlations are different for each NNEMS

         !IF ( N .EQ. 1 ) THEN

            CORR_LX = COV_ERROR_LX
            CORR_LY = COV_ERROR_LY
            SIGMAX  = 1d0
            SIGMAY  = 1d0

            CALL BVECT(IIPAR,  JJPAR,  TMP_XMID, 
     &                 TMP_YMID,
     &                 SIGMAX, SIGMAY, CORR_LX,         CORR_LY,
     &                 -1,     MATINVX, MATY, 
     &                 State_Grid_Adj%DX_COV, State_Grid_Adj%DY_COV   )

         !ENDIF

         DO M = 1, MMSCL
            DO J = 1, JJPAR
            DO I = 1, IIPAR
#if defined ( LOG_OPT )
               ! inverse of error covariance (assume diagonal)
               S2_INV_2D(I,J) = 1d0 / ( EMS_ERROR(N)/EMS_SF0(I,J,1) )
     &                          ** 2
#else
               S2_INV_2D(I,J) = 1d0 / ( EMS_ERROR(N) ) ! 100% error
#endif

               REG_4D(I,J,M,N) = ( EMS_SF(I,J,1) - EMS_SF0(I,J,1) )
     &                          * S2_INV_2D(I,J)

!               IF ( REG_4D(I,J,M,N) .NE. 0.d0 ) THEN
!                  WRITE(*,*) 'REG_4DULAR=', REG_4D(I,J,M,N)
!               ENDIF
        
            ENDDO
            ENDDO

            ! Couple in x
            DO J = 1 , JJPAR
               MATXX(1:IIPAR,1:IIPAR) = MATINVX(J,1:IIPAR,1:IIPAR)

               B(1:IIPAR) = REG_4D(1:IIPAR,J, M, N)

               CALL DPOTRS( 'L', IIPAR, 1, MATXX, IIPAR, B, IIPAR, INFO)

               UOUT(1:IIPAR,J) = B(1:IIPAR)
            ENDDO
            WRITE(*, *) 'MAX UOUT IN X', MAXVAL(UOUT(:, :))
            ! COUPLE IN Y
            DO I = 1, IIPAR
               C(1:JJPAR) = UOUT(I,1:JJPAR)

               CALL DPOTRS( 'L', JJPAR, 1, MATY, JJPAR, C, JJPAR, INFO )

               UOUT(I,1:JJPAR) = C(1:JJPAR)

            ENDDO
            WRITE(*, *) 'MAX UOUT IN Y', MAXVAL(UOUT(:, :))
            DO J = 1, JJPAR
            DO I = 1, IIPAR

               !   IF(J.GE.16.AND.J.LE.60)THEN
               TEMP2(I, J, M, N)   = REG_4D(I,J,M,N)*UOUT(I,J) / 2d0
               EMS_SF_ADJ(I,J,1) = EMS_SF_ADJ(I, J,1) + 
     &                             REG_PARAM_EMS(N)*UOUT(I,J)
     &                             * S2_INV_2D(I,J)

               APCOST(I,J,M,N)     = 0.5d0 * (EMS_SF(I,J,1) -
     &                               EMS_SF0(I,J,1)) *
     &                               UOUT(I,J) * S2_INV_2D(I,J)
            ENDDO
            ENDDO
  
         ENDDO
         ! add regularization parameter from input.gcadj
         APCOST(:,:,:,N) = APCOST(:,:,:,N) * REG_PARAM_EMS(N)
         WRITE(*,*) 'SUM of UOUT =', SUM ( UOUT(:,:) ) , N
         !huoxiao
         WRITE(*,*) 'MAX EMS_SF', MAXVAL( EMS_SF(:, :, :))
      ENDDO
!!!$OMP END PARALLEL DO



      END SUBROUTINE CALC_COV_ERROR

!------------------------------------------------------------------------------

      SUBROUTINE BVECT(N,M,X,Y,SIGMAX,SIGMAY,LX,LY,DIRECTION,
     &                 MATINVX,MATY, DX_COV, DY_COV)

!*******************************************************************************
!
!       C O V A R I A N C E    M A T R I X - V E C T O R    O P E R A T I O N S
!       ( with periodicity )
!
!       INPUTS:
!         N         -  size of x vector ( dimension in x direction )
!         M         -  size of y vector ( dimension in y direction )
!         x         -  longitude coordinates of length N, degrees
!         y         -  latitude  coordinates of length M, degrees
!         sigmax    -  Background standard deviation in x direction
!         sigmay    -  Background standard deviation in y direction
!         lx        -  Correlation Lenghts in x direction (scalar)
!         ly        -  Correlation Lenghts in y direction (scalar)
!         direction -  Correlation Lenghts in x direction (scalar)
!         uin       -  Vector to be multiplied by matrix, dimension 1 x NM
!
!       OUTPUT:
!         if      (direction == 1)
!            uout   =  Covariance Matrix multiplied by uin
!         else if(direction == -1)
!            uout   =  Inverse Covariance Matrix multiplied by uin
!         else if(direction == 2)
!            uout   =  SQRT(Covariance Matrix) multiplied by uin
!         else if(direction == -2)
!            uout   =  Inverse SQRT(Covariance Matrix) multiplied by uin
!
!       DESCRIPTION:
!          Calculates the covariance matrix and it's inverse times a vector
!          without
!          explicitly contructing the covariance matrix
!
!       AUTHORS:
!          Mjardak, Kumaresh Singh
!       DATE:
!          03/03/2009
!*******************************************************************************

      ! References to F90 modules
!      USE GRID_MOD, ONLY : DX_COV, DY_COV

      IMPLICIT NONE

      ! Global variables
      
      INTEGER,INTENT(IN)               :: N,M
      DOUBLE PRECISION,INTENT(IN)      :: SIGMAX,SIGMAY,LX,LY
      DOUBLE PRECISION,INTENT(IN)      :: X(N),Y(M)
      INTEGER, INTENT(IN)              :: DIRECTION
      DOUBLE PRECISION, INTENT(INOUT)  :: MATINVX(M,N,N)
      DOUBLE PRECISION, INTENT(INOUT)  :: MATY(M,M)
      !huoxiao
      REAL*8, INTENT(IN)           :: DX_COV(M)
      REAL*8, INTENT(IN)           :: DY_COV(M)

      ! Local variables
      INTEGER :: I,J,K,L,INFO,II
      DOUBLE PRECISION :: TMP, SUM
      DOUBLE PRECISION :: ALPHA,BETA,DGRIDPTS,DKM
      DOUBLE PRECISION :: TESTVEC(N*M)
      DOUBLE PRECISION :: MATX(M,N,N),MATXX(N,N)
      
      ALPHA = 0.2  !PERTURBATION ON DIAGONAL
      BETA  = 0.2  !PERTURBATION ON DIAGONAL

!------------------------------------------------
!       One dimensional Covariance Matices
!------------------------------------------------
!        X Direction
!------------------------------------------------

      MATX(:,:,:) = 0D0
      DO K = 1 , M
      DO I = 1 , N
         MATX(K,I,I) = SIGMAX * 1.0D0
         DO J = 1 , I-1
            DGRIDPTS =  MIN( ABS( I - J ), N - ABS( I - J ) ) ! DISTANCE IN GRIDPOINTS BETWEEN (I,J)
            DKM = DGRIDPTS * DX_COV(K)             ! DISTANCE IN KM BETWEEN (I,J)
            IF ( DKM <= 3 * LX  ) THEN              !

               MATX(K,I,J) =
     &         EXP( -( DKM / LX )**2 )/( 1.D0 + ALPHA )
            ELSE

               MATX(K,I,J) = 0.D0

            ENDIF
            MATX(K,I,J) = SIGMAX * MATX(K,I,J)
            MATX(K,J,I) = MATX(K,I,J)
         ENDDO

      ENDDO
      ENDDO

!------------------------------------------------
!        Y Direction
!------------------------------------------------
      MATY = 0D0
      DO I = 1 , M

         MATY(I,I) = SIGMAY * 1D0

         DO J = 1 , I - 1

            DKM = ABS(I-J) * DY_COV(M/2)
            IF ( DKM <= 3 * LY ) THEN

               MATY(I,J) = EXP( -( DKM/LY )**2 ) / ( 1.D0 + BETA )

            ELSE

               MATY(I,J) = 0.D0

            ENDIF

            MATY(I,J) = SIGMAY * MATY(I,J)
            MATY(J,I) = MATY(I,J)

         ENDDO
      ENDDO

!-------------------------------------------------------
!       IF DIRECTION IS -1, GET INVERSE OF MATX AND MATY
!-------------------------------------------------------

C$$$	IF(DIRECTION == 1)THEN
C$$$
C$$$           ! COUPLE IN Y
C$$$	   DO K = 1,N
C$$$	      DO L = 1,M
C$$$		 UOUT(K,L) = 0D0
C$$$		 DO I = 1,M
C$$$		    UOUT(K,L) = UOUT(K,L) + MATY(L,I)*UIN(K,I)
C$$$		 END DO
C$$$	      END DO
C$$$	   END DO
C$$$
C$$$	   UMID(:,:) = UOUT(:,:)
C$$$           ! COUPLE IN X
C$$$	   DO K = 1,N
C$$$	      DO L = 1,M
C$$$		 UOUT(K,L) = 0D0
C$$$		 DO I = 1,N
C$$$		    UOUT(K,L) = UOUT(K,L) + MATX(L,K,I)*UMID(I,L)
C$$$		 END DO
C$$$	      END DO
C$$$	   END DO
C$$$
C$$$	ELSE IF(DIRECTION == 2)THEN
C$$$
C$$$	   CALL DSYEVD( 'V', 'L', M, MATY, M, SY, WORKY, 3*M*M, IWORKY,
C$$$     $                   6*M, INFO )
C$$$
C$$$	   DO L = 1 , M
C$$$	      MATXX(:,:) = MATX(L,:,:)
C$$$
C$$$	      CALL DSYEVD( 'V', 'L', N, MATXX, N, SX, WORKX, 3*N*N, IWORKX,
C$$$     $                   6*N, INFO )
C$$$
C$$$	     DO I = 1,N
C$$$	        DO J = 1,N
C$$$	           UX(L,I,J) = 0D0
C$$$	           DO K = 1,N
C$$$	               UX(L,I,J) = UX(L,I,J)+(MATXX(I,K)*SQRT(SX(K))*
C$$$     $                                        MATXX(J,K))
C$$$	           END DO
C$$$	        END DO
C$$$	     END DO
C$$$           END DO
C$$$!
C$$$	   DO I = 1,M
C$$$     	      DO J = 1,M
C$$$	         UY(I,J) = 0D0
C$$$	         DO K = 1,M
C$$$	            UY(I,J) = UY(I,J) + (MATY(I,K)*SQRT(SY(K))*MATY(J,K))
C$$$	         END DO
C$$$	      END DO
C$$$	   END DO
C$$$
C$$$           ! COUPLE IN Y
C$$$	   DO K = 1,N
C$$$	      DO L = 1,M
C$$$		 UOUT(K,L) = 0D0
C$$$		 DO I = 1,M
C$$$		    UOUT(K,L) = UOUT(K,L) + UY(L,I)*UIN(K,I)
C$$$		 END DO
C$$$	      END DO
C$$$	   END DO
C$$$
C$$$	   UMID(:,:) = UOUT(:,:)
C$$$           ! COUPLE IN X
C$$$	   DO K = 1,N
C$$$	      DO L = 1,M
C$$$		 UOUT(K,L) = 0D0
C$$$		 DO I = 1,N
C$$$		    UOUT(K,L) = UOUT(K,L) + UX(L,K,I)*UMID(I,L)
C$$$		 END DO
C$$$	      END DO
C$$$	   END DO
C$$$
C$$$	ELSE IF (DIRECTION == -1) THEN

      DO L = 1 , M
         MATXX(:,:) = MATX(L,:,:)
         ! Cholesky decomposition of x-corr
         CALL DPOTRF( 'L', N, MATXX, N, INFO )
         MATINVX(L,:,:) = MATXX(:,:)
      ENDDO

      ! Cholesky decomposition of y-corr
      CALL DPOTRF( 'L', M, MATY, M, INFO )

C$$$	ELSE IF(DIRECTION == -2)THEN
C$$$
C$$$	   CALL DSYEVD( 'V', 'L', M, MATY, M, SY, WORKY, 3*M*M, IWORKY,
C$$$     $                   6*M, INFO )
C$$$
C$$$
C$$$           DO L = 1 , M
C$$$              MATXX(:,:) = MATX(L,:,:)
C$$$      	     CALL DSYEVD('V','L',N,MATXX,N,SX,WORKX,3*N*N,IWORKX,
C$$$     $                   6*N, INFO )
C$$$
C$$$	     DO I = 1,N
C$$$	        DO J = 1,N
C$$$	           UX(L,I,J) = 0D0
C$$$	           DO K = 1,N
C$$$	               UX(L,I,J) = UX(L,I,J)+(MATXX(I,K)*MATXX(J,K))
C$$$     $                                        /SQRT(SX(K))
C$$$	           END DO
C$$$	        END DO
C$$$	     END DO
C$$$           END DO
C$$$
C$$$	   DO I = 1,M
C$$$     	      DO J = 1,M
C$$$	         UY(I,J) = 0D0
C$$$	         DO K = 1,M
C$$$	            UY(I,J) = UY(I,J) + (MATY(I,K)*MATY(J,K))
C$$$     $                                   /SQRT(SY(K))
C$$$	         END DO
C$$$	      END DO
C$$$	   END DO
C$$$
C$$$           ! COUPLE IN X
C$$$	   DO K = 1,N
C$$$	      DO L = 1,M
C$$$		 UOUT(K,L) = 0D0
C$$$		 DO I = 1,N
C$$$		    UOUT(K,L) = UOUT(K,L) + UX(L,K,I)*UIN(I,L)
C$$$		 END DO
C$$$	      END DO
C$$$	   END DO
C$$$
C$$$	   UMID(:,:) = UOUT(:,:)
C$$$           ! COUPLE IN Y
C$$$	   DO K = 1,N
C$$$	      DO L = 1,M
C$$$		 UOUT(K,L) = 0D0
C$$$		 DO I = 1,M
C$$$		    UOUT(K,L) = UOUT(K,L) + UY(L,I)*UMID(K,I)
C$$$		 END DO
C$$$	      END DO
C$$$	   END DO
C$$$
C$$$	ENDIF

!--------------------------------------------------------
      RETURN

      END SUBROUTINE BVECT

!========================================================================

      ! End of program
      END MODULE COVARIANCE_MOD


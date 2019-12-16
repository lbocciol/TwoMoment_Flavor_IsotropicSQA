MODULE CollisionsModule
  
  USE KindModule, ONLY: &
    DP, Zero, Half, One
  USE UnitsModule, ONLY: &
    Second
  USE ProgramHeaderModule, ONLY: &
    nX, nDOFX, iX_B0
  USE InitializationModule, ONLY: &
    nX_G, nE_G, nM, nF, &
    fMatrixOld, fMatrixOsc, &
    fImaginaryPart
  USE RadiationFieldsModule, ONLY: &
    iNuEX, iNuEX_Bar, nSpecies
  USE IntegrationModule, ONLY: &
    TimeStepping_Decoherence, &
    kappa_avg, eta_avg, &
    CollisionTermNu4
  USE OscillationsModule, ONLY: &
    Packf

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: CollisionsDriver
  PUBLIC :: EvaluateImpact

CONTAINS

  SUBROUTINE CollisionsDriver(TimeBlockEnd, dt, dt_coll)

    REAL(DP), INTENT(IN)    :: TimeBlockEnd, dt
    REAL(DP), INTENT(INOUT) :: dt_coll

    REAL(DP) :: TimeBlock, dt_loc
    LOGICAL  :: BlockDone

    BlockDone = .FALSE.
    dt_loc = dt
    TimeBlock = Zero

    DO WHILE( BlockDone .EQV. .FALSE. )

      IF( TimeBlock + dt_loc >= TimeBlockEnd ) THEN

        dt_loc = TimeBlockEnd - TimeBlock
        BlockDone = .TRUE.

      END IF

      CALL TimeStepping_Decoherence( dt_loc )
      
      ! --- HANDLE IMAGINARY PART OF OFF-DIAGONAL --- !
      ! ----------- ELEMENTS SEPARATELY ------------- !

      CALL SolveImaginaryPart( dt_loc )

      DEALLOCATE( kappa_avg, eta_avg, CollisionTermNu4 )
      
      TimeBlock = TimeBlock + dt_loc

      IF( .NOT. BlockDone ) dt_coll = dt_loc

    END DO

  END SUBROUTINE CollisionsDriver

  SUBROUTINE SolveImaginaryPart( dt_loc )

    REAL(DP), INTENT(IN) :: dt_loc
    INTEGER  :: iX, iS, iE

    !$OMP PARALLEL DO COLLAPSE(2) PRIVATE( iX, iE, iS )
    DO iX = 1,nX_G
      DO iE = 1,nE_G 
    
        iS = iNuEX
        fImaginaryPart(iX,iE,1) = &
           ( fImaginaryPart(iX,iE,1) + &
           dt_loc * CollisionTermNu4(iE,iX,iS) ) / &
           ( One + ( eta_avg(iE,iX,iS) + kappa_avg(iE,iX,iS) ) * dt_loc )
    
        iS = iNuEX_Bar
        fImaginaryPart(iX,iE,2) = &
           ( fImaginaryPart(iX,iE,2) + &
           dt_loc * CollisionTermNu4(iE,iX,iS) ) / & 
           ( One + ( eta_avg(iE,iX,iS) + kappa_avg(iE,iX,iS) ) * dt_loc )
      
      END DO
    END DO
    !$OMP END PARALLEL DO
    
  END SUBROUTINE SolveImaginaryPart

  SUBROUTINE EvaluateImpact( dt_rand, dt_block )

    REAL(DP), INTENT(IN)    :: dt_rand
    REAL(DP), INTENT(INOUT) :: dt_block

    INTEGER  :: iX, iX1, iX2, iX3, iNodeX
    INTEGER  :: m, iE, iS1, iS2
    REAL(DP) :: TargetImpact = 1.0d-3
    REAL(DP) :: Increase = 1.1
    REAL(DP) :: LengthL, Impact, CorrectedImpact
    REAL(DP) :: df(nF,nF)

    Impact = Zero
    DO iX = 1,nX_G

        iX3    = MOD( (iX-1) / ( nDOFX * nX(1) * nX(2) ), nX(3) ) + iX_B0(3)
        iX2    = MOD( (iX-1) / ( nDOFX * nX(1)         ), nX(2) ) + iX_B0(2)
        iX1    = MOD( (iX-1) / ( nDOFX                 ), nX(1) ) + iX_B0(1)
        iNodeX = MOD( (iX-1)                            , nDOFX ) + 1


      ! --- When you pack, add back Imaginary Part --- !
      CALL Packf( iX, iNodeX, iX1, iX2, iX3 )

      DO m = 1,nM
        DO iE = 1,nE_G
          
          CALL IsospinL( fMatrixOld(iX,m,iE,:,:) , Lengthl )

          df = fMatrixOld(iX,m,iE,:,:) - fMatrixOsc(m,iE,:,:)
	   
          DO iS1 = 1,nF
            DO iS2 = 1,nF
            
              Impact = MAX( Impact , ABS( df(iS1,iS2) ) / Lengthl )
            
            END DO
          END DO

        END DO
      END DO
    
    END DO

    IF( Impact > TargetImpact ) &
      WRITE(*,'(A,E12.3)') "WARNING: Impact=", Impact
    
    CorrectedImpact = Impact / ( dt_rand ) * dt_block;

    IF( CorrectedImpact < 0.1d0 * TargetImpact ) &
      dt_block = dt_block * MIN( Increase, 0.1d0 * TargetImpact / &
                                                CorrectedImpact )
    
    IF( CorrectedImpact > 0.1d0 * TargetImpact ) &
      dt_block = dt_block * 0.1d0 * TargetImpact / CorrectedImpact

  END SUBROUTINE EvaluateImpact

  SUBROUTINE IsospinL( Matrix, Length )

    COMPLEX(DP), INTENT(IN) :: Matrix(nF,nF)
    REAL(DP), INTENT(OUT)   :: Length

    REAL(DP) :: fx, fy, fz

    fx = REAL ( Matrix(1,2) )
    fy = AIMAG( Matrix(1,2) )
    fz = Half * ABS( Matrix(1,1) - &
                     Matrix(2,2) )

    Length = SQRT( fx**2 + fy**2 + fz**2 ) 

  END SUBROUTINE IsospinL

END MODULE CollisionsModule

MODULE FlavorOpacitiesModule
    
  USE KindModule, ONLY: &
    DP, Zero, Half, One, &
    Two, Three, Four, Five, &
    Pi, TwoPi, FourPi
  USE UnitsModule, ONLY: &
    BoltzmannConstant, &
    MeV, Second
  USE PhysicalConstantsModule, ONLY: &
    SpeedOfLightCGS
  USE ProgramHeaderModule, ONLY: &
    nDOF, nDOFX, iX_E0, iX_B0 
  USE FluidFieldsModule, ONLY: &
    uAF, iAF_T, &
    iAF_Me, iAF_Mp, iAF_Mn
  USE RadiationFieldsModule, ONLY: &
    nSpecies, iNuE, iNuE_Bar, &
    iNuEE, iNuXX, iNuEX, &
    iNuEE_Bar, iNuXX_Bar, iNuEX_Bar
  USE ReadProfileModule, ONLY: &
    LinearInterpolation1D, &
    QuadraticInterpolation1D
  USE InitializationModule, ONLY: &
    nF,IdentityR
  USE KernelsNu4Module, ONLY: &
    KernelNu4Pair, &
    KernelNu4Scat

  PRIVATE
  
  PUBLIC :: ComputeEmission
  PUBLIC :: ComputeEmissionAbsorptionAvg
  PUBLIC :: ComputeNu4CollisionTerms
  PUBLIC :: PackIntoMatrix
  PUBLIC :: FlattenMatrix

CONTAINS

  SUBROUTINE ComputeEmission(E_N, Chi, kappa, eta, &
                             nE_G, nX_G)

    INTEGER,  INTENT(IN)  :: nE_G, nX_G
    REAL(DP), INTENT(IN)  :: E_N  (nE_G)
    REAL(DP), INTENT(IN)  :: Chi  (nE_G,nX_G,nF)
    REAL(DP), INTENT(OUT) :: eta  (nE_G,nX_G,nF)
    REAL(DP), INTENT(OUT) :: kappa(nE_G,nX_G,nF)

    INTEGER  :: iS, iN_E, iN_X
    INTEGER  :: iX1, iX2, iX3, iNodeX
    INTEGER  :: nX(3)
    REAL(DP) :: Mnu, kT, FD

    nX = iX_E0 - iX_B0 + 1

    DO iN_X = 1, nX_G
  
      iX3    = MOD( (iN_X-1) / ( nDOFX * nX(1) * nX(2) ), nX(3) ) + iX_B0(3)
      iX2    = MOD( (iN_X-1) / ( nDOFX * nX(1)         ), nX(2) ) + iX_B0(2)
      iX1    = MOD( (iN_X-1) / ( nDOFX                 ), nX(1) ) + iX_B0(1) 
      iNodeX = MOD( (iN_X-1)                            , nDOFX ) + 1    
    
      DO iS = 1, nF
      
        IF ( iS .EQ. iNuE ) THEN
  
          Mnu = + ( uAF(iNodeX,iX1,iX2,iX3,iAF_Me) &
                      + uAF(iNodeX,iX1,iX2,iX3,iAF_Mp) &
                      - uAF(iNodeX,iX1,iX2,iX3,iAF_Mn) )
          
        ELSE IF ( iS .EQ. iNuE_Bar ) THEN
            
          Mnu = - ( uAF(iNodeX,iX1,iX2,iX3,iAF_Me) &
                      + uAF(iNodeX,iX1,iX2,iX3,iAF_Mp) &
                      - uAF(iNodeX,iX1,iX2,iX3,iAF_Mn) )
        ELSE
  
          Mnu = Zero
  
        END IF
          
        kT = BoltzmannConstant &
             * uAF(iNodeX,iX1,iX2,iX3,iAF_T)
  
        DO iN_E = 1, nE_G
  
          FD = MAX( One / ( EXP( (E_N(iN_E)-Mnu)/kT ) + One ), 1.0d-99 )
          
          eta(iN_E,iN_X,iS) = Chi(iN_E,iN_X,iS) * FD 
          kappa(iN_E,iN_X,iS) = Chi(iN_E,iN_X,iS) - eta(iN_E,iN_X,iS)
          
        END DO

      END DO
    
    END DO

  END SUBROUTINE ComputeEmission

  SUBROUTINE ComputeEmissionABSorptionAvg ( &
               kappa_avg, eta_avg, kappa, eta,   &
               nE_G, nX_G)
    
    INTEGER,  INTENT(IN)  :: nE_G, nX_G
    REAL(DP), INTENT(IN)  :: kappa  (nE_G,nX_G,nF)
    REAL(DP), INTENT(IN)  :: eta    (nE_G,nX_G,nF)
    REAL(DP), INTENT(OUT) :: kappa_avg(nE_G,nX_G,nSpecies)
    REAL(DP), INTENT(OUT) :: eta_avg(nE_G,nX_G,nSpecies)

    INTEGER :: iS

    DO iS = 1,nSpecies
      
      IF ( iS .EQ. iNuEE ) THEN

        kappa_avg(:,:,iS) = kappa(:,:,iNuE)
        eta_avg(:,:,iS) = eta(:,:,iNuE)

      ELSE IF ( iS .EQ. iNuEE_Bar ) THEN

        kappa_avg(:,:,iS) = kappa(:,:,iNuE_Bar)
        eta_avg(:,:,iS) = eta(:,:,iNuE_Bar)

      ELSE IF ( iS .EQ. iNuEX ) THEN

        kappa_avg(:,:,iS) = Half * kappa(:,:,iNuE)
        eta_avg(:,:,iS) = Half * eta(:,:,iNuE)

      ELSE IF ( iS .EQ. iNuEX_Bar ) THEN

        kappa_avg(:,:,iS) = Half * kappa(:,:,iNuE_Bar)
        eta_avg(:,:,iS) = Half * eta(:,:,iNuE_Bar)

      ELSE

        kappa_avg(:,:,iS) = Zero
        eta_avg(:,:,iS) = Zero

      END IF

    END DO

  END SUBROUTINE ComputeEmissionAbsorptionAvg


  SUBROUTINE ComputeNu4CollisionTerms( J0, CollisionTerm, &
                                       E_N, nE_G, iEk )
    
    INTEGER,  INTENT(IN)  :: nE_G, iEk
    REAL(DP), INTENT(IN)  :: E_N(nE_G)
    REAL(DP), INTENT(IN)  :: J0           (nE_G,nSpecies)
    REAL(DP), INTENT(OUT) :: CollisionTerm(nSpecies)
    
    INTEGER  :: i1,i2,i3
    INTEGER  :: iE
    INTEGER  :: nMatter, nAntiMatter
    REAL(DP) :: J0_E(2,nE_G,nF,nF)

    nMatter     = nSpecies / 2
    nAntiMatter = nSpecies / 2 + 1

    ! 1 is for matter, 2 is for antimatter
    DO iE = 1,nE_G

      J0_E(1,iE,:,:) = PackIntoMatrix( J0(iE,1:nMatter   ), nSpecies/2, nF )
      J0_E(2,iE,:,:) = PackIntoMatrix( J0(iE,nAntiMatter:), nSpecies/2, nF )

    END DO

    CALL ComputeNu4CollisionTerm( J0_E, CollisionTerm(1:nMatter), &
        E_N, nE_G, iEk )

    !Flip them to get the Collision term for antimatter
    DO iE = 1,nE_G

      J0_E(2,iE,:,:) = PackIntoMatrix( J0(iE,1:nMatter   ), nSpecies/2, nF )
      J0_E(1,iE,:,:) = PackIntoMatrix( J0(iE,nAntiMatter:), nSpecies/2, nF )

    END DO

    CALL ComputeNu4CollisionTerm( J0_E, CollisionTerm(nAntiMatter:), &
        E_N, nE_G, iEk )

  END SUBROUTINE ComputeNu4CollisionTerms


  SUBROUTINE ComputeNu4CollisionTerm( J0_E, CollisionTerm, &
                                      E_N, nE_G, iEk )

    INTEGER,  INTENT(IN)  :: nE_G, iEk
    REAL(DP), INTENT(IN)  :: J0_E         (2,nE_G,nF,nF)
    REAL(DP), INTENT(IN)  :: E_N(nE_G)
    REAL(DP), INTENT(OUT) :: CollisionTerm(nSpecies/2)

    REAL(DP) :: IdentityR      (nF,nF)
    REAL(DP) :: FirstTerm      (nF,nF)
    REAL(DP) :: SecondTerm     (nF,nF)
    REAL(DP) :: FirstTermG     (nF,nF)
    REAL(DP) :: SecondTermG    (nF,nF)
    REAL(DP) :: AntiCommutator1(nF,nF)
    REAL(DP) :: AntiCommutator2(nF,nF)
    REAL(DP) :: CollisionMatrix(nF,nF)
    REAL(DP) :: tmp            (nF,nF)
    REAL(DP) :: J0_E2(2,nF,nF)
    REAL(DP) :: E2

    INTEGER  :: iE, iE1, iE3
    
    IdentityR(:,:) = Zero

    FirstTerm (:,:) = Zero
    SecondTerm(:,:) = Zero

    FirstTermG (:,:) = Zero
    SecondTermG(:,:) = Zero

    DO iE1 = 1,nE_G
    DO iE3 = 1,nE_G
    
      E2 = E_N(iE1) + E_N(iE3) - E_N(iEk)
    
      IF ( (E2 >= E_N(1)) .AND. (E2 < E_N(nE_G)) ) THEN

        DO iS1 = 1,nF
          DO iS2 = 1,nF
  
            CALL QuadraticInterpolation1D( E_N, J0_E(1,:,iS1,iS2), nE_G, &
                  E2, J0_E2(1,iS1,iS2) )
  
            CALL QuadraticInterpolation1D( E_N, J0_E(2,:,iS1,iS2), nE_G, &
                  E2, J0_E2(2,iS1,iS2) )
  
          END DO
        END DO

        ! Scattering on nu, i.e. 3rd line in eq. 96 of Blaschke (2016)
        tmp = MatMul( IdentityR - J0_E(1,iE1,:,:) , J0_E2(1,:,:) )
        FirstTerm = FirstTerm + MatMul( TraceI( tmp, nF ) + tmp , &
                      IdentityR - J0_E(1,iE3,:,:) ) * &
                         KernelNu4Scat(iE3,iE1,iEk)
        
        ! Scattering on nubar, i.e. first part of 5th line in eq. 96 of Blaschke (2016)
        tmp = MatMul( J0_E2(2,:,:) , IdentityR - J0_E(2,iE1,:,:) )
        SecondTerm = SecondTerm + MatMul( TraceI( tmp, nF ) + tmp , &
                        IdentityR - J0_E(1,iE3,:,:) ) * &
                           KernelNu4Pair(iE3,iE1,iEk)

        ! Scattering on nubar, i.e. second part of 5th line in eq. 96 of Blaschke (2016)
        tmp = MatMul( IdentityR - J0_E(1,iE3,:,:) , IdentityR - J0_E(2,iE1,:,:) )
        SecondTerm = SecondTerm + &
              MatMul( TraceI( tmp, nF ) + tmp , J0_E2(2,:,:) ) * &
                KernelNu4Pair(iE3,iE1,iEk)

        ! Gain part
        ! Scattering on nu, i.e. 3rd line in eq. 96 of Blaschke (2016)
        tmp = MatMul( J0_E(1,iE1,:,:) , IdentityR - J0_E2(1,:,:) )
        FirstTermG = FirstTermG + MatMul( TraceI( tmp, nF ) + tmp , &
                        J0_E(1,iE3,:,:) ) * KernelNu4Scat(iE3,iE1,iEk)

        ! Scattering on nubar, i.e. first part of 5th line in eq. 96 of Blaschke (2016)
        tmp = MatMul( IdentityR - J0_E2(2,:,:) , J0_E(2,iE1,:,:) )
        SecondTermG = SecondTermG + MatMul( TraceI( tmp, nF ) + tmp , &
                         J0_E(1,iE3,:,:) ) * KernelNu4Pair(iE3,iE1,iEk)

        ! Scattering on nubar, i.e. second part of 5th line in eq. 96 of Blaschke (2016)
        tmp = MatMul( J0_E(1,iE3,:,:) , J0_E(2,iE1,:,:) )
        SecondTermG = SecondTermG + MatMul( TraceI( tmp, nF ) + tmp , &
                         IdentityR - J0_E2(2,:,:) ) * &
                            KernelNu4Pair(iE3,iE1,iEk)

      ENDIF

    END DO
    END DO
    
    AntiCommutator1 = MatMul( FirstTerm       , J0_E(1,iEk,:,:) ) + &
                      MatMul( J0_E(1,iEk,:,:) , FirstTerm       )

    AntiCommutator2 = MatMul( SecondTerm      , J0_E(1,iEk,:,:) ) + &
                      MatMul( J0_E(1,iEk,:,:) , SecondTerm      )

    CollisionMatrix = - ( AntiCommutator1 + AntiCommutator2 )

    !Gain part
    AntiCommutator1 = MatMul( FirstTermG       , IdentityR - J0_E(1,iEk,:,:) ) + &
                      MatMul( IdentityR - J0_E(1,iEk,:,:) , FirstTermG       )

    AntiCommutator2 = MatMul( SecondTermG      , IdentityR - J0_E(1,iEk,:,:) ) + &
                      MatMul( IdentityR - J0_E(1,iEk,:,:) , SecondTermG      )

    CollisionMatrix = CollisionMatrix + ( AntiCommutator1 + AntiCommutator2 )

    CollisionTerm   = FlattenMatrix( CollisionMatrix, nF, nSpecies / 2 )

  END SUBROUTINE ComputeNu4CollisionTerm


  FUNCTION PackIntoMatrix( Vector, dimV, dimM ) RESULT( Matrix )
    
    INTEGER,  INTENT(IN)  :: dimV, dimM
    REAL(DP), INTENT(IN)  :: Vector( dimV )
    REAL(DP)              :: Matrix( dimM,dimM )
    
    ! THIS WORKS ONLY IF THE ORIGINAL FLAVOR MATRIX WAS NUMBERED
    ! ROW-WISE : i.e. 1 2 3      1 2 
    !                   4 5  AND   3
    !                     6

    ! FILL UPPER PART
    DO i = 1,dimM
      DO j = i, dimM

        IF ( j == 3 ) THEN
          Matrix(i,j) = Vector(i+j)
        ELSE 
          Matrix(i,j) = Vector(i+j-1)
        END IF
      
      END DO
    END DO
    
    ! FILL LOWER PART
    DO i = 1,dimM-1
      DO j = i+1,dimM
        
        IF ( j == 3 ) THEN
          Matrix(j,i) = Vector(i+j)
        ELSE
          Matrix(j,i) = Vector(i+j-1)
        END IF

      END DO
    END DO

  END FUNCTION PackIntoMatrix
    
  FUNCTION FlattenMatrix( Matrix, dimM, dimV ) RESULT( Vector )

    INTEGER,  INTENT(IN)  :: dimV, dimM
    REAL(DP), INTENT(IN)  :: Matrix( dimM,dimM )
    REAL(DP)              :: Vector( dimV )

    ! THIS WORKS ONLY IF THE ORIGINAL FLAVOR MATRIX WAS NUMBERED
    ! ROW-WISE : i.e. 1 2 3      1 2 
    !                   4 5  AND   3
    !                     6

    ! FLATTEN UPPER PART, i.e. the one that is stored
    DO i = 1,dimM
      DO j = i, dimM

        IF ( j == 3 ) THEN
          Vector(i+j) = Matrix(i,j)
        ELSE
          Vector(i+j-1) = Matrix(i,j)
        END IF

      END DO
    END DO

  END FUNCTION FlattenMatrix


  FUNCTION TraceI( Matrix, dimM ) RESULT( Trace )
    
    INTEGER,  INTENT(IN)  :: dimM
    REAL(DP), INTENT(IN)  :: Matrix(dimM,dimM)
    REAL(DP)              :: Trace (dimM,dimM)

    INTEGER :: i, j

    Trace(:,:) = Zero

    DO i = 1,dimM
      DO j = 1,dimM

        Trace(i,i) = Trace(i,i) + Matrix(j,j)

      END DO
    END DO

  END FUNCTION TraceI

  FUNCTION MatrixMul( Matrix1, Matrix2, dimM) RESULT(Matrix3)

    INTEGER,  INTENT(IN) :: dimM
    REAL(DP), INTENT(IN) :: Matrix1(dimM,dimM), Matrix2(dimM,dimM)
    REAL(DP)             :: Matrix3(dimM,dimM)
    
    INTEGER :: i,j,k
    
    DO j = 1,dimM
      DO k = 1,dimM
        DO i = 1,dimM
        
          Matrix3(i,j) = Matrix1(i,j) + Matrix1(i,k) * Matrix2(k,j)

        END DO
      END DO
    END DO
  
  END FUNCTION MatrixMul

END MODULE FlavorOpacitiesModule

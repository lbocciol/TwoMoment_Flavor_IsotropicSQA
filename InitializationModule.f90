MODULE InitializationModule

  USE KindModule, ONLY: &
    DP, Zero, Half, One, Two, &
    Five, Pi
  USE UnitsModule, ONLY: &
    Gram, &
    Centimeter, &
    Kilometer,  &
    Kelvin, MeV, &
    BoltzmannConstant
  USE ProgramHeaderModule, ONLY: &       
    iX_B0, iX_E0, &                      
    iE_B0, iE_E0, &
    iZ_B0, iZ_E0, &
    nDOF, nDOFX,  &
    nDOFE, nNodesZ
  USE ReferenceElementModuleX, ONLY: &
    NodeNumberTableX
  USE ReferenceElementModule, ONLY: &  
    NodeNumbersX, &                    
    NodeNumberTable                    
  USE MeshModule, ONLY: &
    MeshE, MeshX, &
    NodeCoordinate
  USE GeometryFieldsModule, ONLY: &                
    CoordinateSystem, &                            
    uGF, iGF_Gm_dd_11, iGF_Gm_dd_22, iGF_Gm_dd_33  
  USE FluidFieldsModule, ONLY: &
    uPF, iPF_D, &
    uAF, iAF_P, iAF_T, iAF_Ye, iAF_S, iAF_E, &
    iAF_Me, iAF_Mp, iAF_Mn, iAF_Xp, iAF_Xn, &
    iAF_Xa, iAF_Xh, iAF_Gm
  USE RadiationFieldsModule, ONLY: &                 
    nSpecies, iNuEE, iNuXX, iNuEX, &
    iNuEE_Bar, iNuXX_Bar, iNuEX_Bar, &
    uPR, nPR, iPR_D, iPR_I1, iPR_I2, iPR_I3, &       
    uCR, nCR, iCR_N, iCR_G1, iCR_G2, iCR_G3          
  USE EquationOfStateModule_TABLE, ONLY: &
    ApplyEquationOfState_TABLE
  USE TwoMoment_UtilitiesModule, ONLY: &  
    ComputeConserved_TwoMoment            
  USE ReadProfileModule, ONLY: &
    Read1DChimeraProfile, &
    ReadGR1DProfile, &
    LinearInterpolation1D, &
    QuadraticInterpolation1D
  USE KernelsNu4Module, ONLY: &
    ComputeNu4Kernels

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeFields_Decoherence
  PUBLIC :: InitializeSingleZone_Decoherence
  PUBLIC :: InitializeOscillations
  PUBLIC :: CofactorMatrix
   
  ! --- OSCILLATIONS RELATED QUANTITIES --- !
  ! --- DIMENSIONS OF ARRAYS --- !
  INTEGER,            PUBLIC :: nX_G, nE_G
  INTEGER,            PUBLIC :: nF
  INTEGER, PARAMETER, PUBLIC :: nM = 2
  INTEGER, PARAMETER, PUBLIC :: nY = 6
  INTEGER, PARAMETER, PUBLIC :: nS = 2

  REAL(DP),                 PUBLIC :: TimeLocOutput
  REAL(DP),    ALLOCATABLE, PUBLIC :: Enu(:)
  REAL(DP),    ALLOCATABLE, PUBLIC :: IdentityR(:,:)
  COMPLEX(DP), ALLOCATABLE, PUBLIC :: YIdentity(:,:)
  COMPLEX(DP), ALLOCATABLE, PUBLIC :: IdentityC(:,:)
  COMPLEX(DP), ALLOCATABLE, PUBLIC :: AV(:,:,:), CV(:,:,:,:)
  COMPLEX(DP), ALLOCATABLE, PUBLIC :: Hvf(:,:,:,:) !vacuum hamiltonian in flavor basis
  COMPLEX(DP),              PUBLIC :: PMNS(2,2)
  COMPLEX(DP), ALLOCATABLE, PUBLIC :: fMatrixOsc(:,:,:,:)
  COMPLEX(DP), ALLOCATABLE, PUBLIC :: fMatrixOld(:,:,:,:,:)
  REAL(DP),    ALLOCATABLE, PUBLIC :: fImaginaryPart(:,:,:)

  REAL(DP),    ALLOCATABLE, PUBLIC :: Energies(:)

  ! --- RUNGE KUTTA PARAMETERS --- !
  INTEGER,  PARAMETER, PUBLIC :: nRK = 6
  REAL(DP), PARAMETER, PUBLIC :: nRKOrder = Five
  REAL(DP),            PUBLIC :: AA(6,5)
  REAL(DP), PARAMETER, PUBLIC :: BB5(6)=(/ 37.0d0/378.0d0 , 0.0d0 , 250.0d0/621.0d0 , &
                               125.0d0/594.0d0 , 0.0d0 , 512.0d0/1771.0d0 /)
  REAL(DP), PARAMETER, PUBLIC :: BB6(6)=(/ 2825.0d0/27648.0d0 , 0.0d0 , 18575.0d0/48384.0d0 , &
                               13525.0d0/55296.0d0 , 277.0d0/14336.0d0 , 1.0d0/4.0d0 /)
  REAL(DP), PARAMETER :: b2(5)=(/ 1.0d0/5.0d0 , 0.0d0 , 0.0d0 , 0.0d0 , 0.0d0/)
  REAL(DP), PARAMETER :: b3(5)=(/ 3.0d0/40.0d0 , 9.0d0/40.0d0 , 0.0d0 , 0.0d0 , 0.0d0/)
  REAL(DP), PARAMETER :: b4(5)=(/ 3.0d0/10.0d0 , -9.0d0/10.0d0 , 6.0d0/5.0d0 , 0.0d0 , 0.0d0 /)
  REAL(DP), PARAMETER :: b5(5)=(/ -11.0d0/54.0d0 , 5.0d0/2.0d0 , &
                            -70.0d0/27.0d0 , 35.0d0/27.0d0 ,0.0d0/)
  REAL(DP), PARAMETER :: b6(5)=(/ 1631.0d0 /55296.0d0 , 175.0d0/512.0d0 , 575.0d0/13824.0d0 , &
                              44275.0d0/110592.0d0 , 253.0d0/4096.0d0 /)

CONTAINS

  SUBROUTINE InitializeFields_Decoherence
    
    CALL InitializeFluidFields_Decoherence
    CALL InitializeRadiationFields_Decoherence
    
    CALL ComputeNu4Kernels

  END SUBROUTINE InitializeFields_Decoherence

  SUBROUTINE InitializeSingleZone_Decoherence        
                                                 
    CALL InitializeSingleZoneFields_Decoherence       
    CALL InitializeRadiationFields_Decoherence   

    CALL ComputeNu4Kernels

  END SUBROUTINE InitializeSingleZone_Decoherence    

  SUBROUTINE InitializeFluidFields_Decoherence

    INTEGER  :: iX1, iX2, iX3, iNodeX, iNodeX1
    REAL(DP) :: R, Y_New
    
    LOGICAL               :: ProfileFromGR1D
    INTEGER               :: n0
    INTEGER               :: FileNumberMax
    REAL(DP), ALLOCATABLE :: R_Ch(:), D_Ch(:), T_Ch(:), Ye_Ch(:)
    REAL(DP)              :: TimeSlice = 1d-1
    
    ProfileFromGR1D = .FALSE.

    IF ( ProfileFromGR1D ) THEN
      
      n0 = 600
      ALLOCATE(R_Ch(n0),D_Ch(n0),T_Ch(n0),Ye_Ch(n0))
      
      CALL ReadGR1DProfile(R_Ch, D_Ch, T_Ch, Ye_Ch, & 
            n0, TimeSlice)

      R_Ch = R_Ch * Centimeter
      D_Ch = D_Ch * Gram / Centimeter ** 3
      T_Ch = T_Ch * MeV

    ELSE

      FileNumberMax = 1000
      n0 = 722
      ALLOCATE(R_Ch(n0),D_Ch(n0),T_Ch(n0),Ye_Ch(n0)) 
      
      CALL Read1DChimeraProfile(R_Ch, D_Ch, T_Ch, Ye_Ch, &
            n0, TimeSlice, FileNumberMax) 

      R_Ch = R_Ch * Centimeter
      D_Ch = D_Ch * Gram / Centimeter ** 3
      T_Ch = T_Ch * Kelvin

    END IF

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)

        R = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

        ! Interpolate Chimera Data to R and set variables
        
        CALL QuadraticInterpolation1D( R_Ch, D_Ch, n0, R,  &
                    uPF(iNodeX,iX1,iX2,iX3,iPF_D) )

        CALL QuadraticInterpolation1D( R_Ch, T_Ch, n0, R,  &
                    uAF(iNodeX,iX1,iX2,iX3,iAF_T) )

        CALL QuadraticInterpolation1D( R_Ch, Ye_Ch, n0, R, &
                    uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) )
        
      END DO

      CALL ApplyEquationOfState_TABLE &
             ( uPF(:,iX1,iX2,iX3,iPF_D ), uAF(:,iX1,iX2,iX3,iAF_T ), &
               uAF(:,iX1,iX2,iX3,iAF_Ye), uAF(:,iX1,iX2,iX3,iAF_P ), &
               uAF(:,iX1,iX2,iX3,iAF_S ), uAF(:,iX1,iX2,iX3,iAF_E ), &
               uAF(:,iX1,iX2,iX3,iAF_Me), uAF(:,iX1,iX2,iX3,iAF_Mp), &
               uAF(:,iX1,iX2,iX3,iAF_Mn), uAF(:,iX1,iX2,iX3,iAF_Xp), &
               uAF(:,iX1,iX2,iX3,iAF_Xn), uAF(:,iX1,iX2,iX3,iAF_Xa), &
               uAF(:,iX1,iX2,iX3,iAF_Xh), uAF(:,iX1,iX2,iX3,iAF_Gm) )

    END DO
    END DO
    END DO

  END SUBROUTINE InitializeFluidFields_Decoherence

  SUBROUTINE InitializeSingleZoneFields_Decoherence
                                                                       
    INTEGER  :: iX1, iX2, iX3, iNodeX, iNodeX1                         
    REAL(DP) :: R, Y_New                                               
  
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)

        R = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

        uPF(iNodeX,iX1,iX2,iX3,iPF_D) = 1.0d12 * Gram / Centimeter**3
        !uPF(iNodeX,iX1,iX2,iX3,iPF_D) = 2.541656d3 * Gram / Centimeter**3

        uAF(iNodeX,iX1,iX2,iX3,iAF_T) = 10.0 * MeV

        uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) = 0.3d0                     

      END DO                                                         
      
      CALL ApplyEquationOfState_TABLE &
             ( uPF(:,iX1,iX2,iX3,iPF_D ), uAF(:,iX1,iX2,iX3,iAF_T ), &
               uAF(:,iX1,iX2,iX3,iAF_Ye), uAF(:,iX1,iX2,iX3,iAF_P ), &
               uAF(:,iX1,iX2,iX3,iAF_S ), uAF(:,iX1,iX2,iX3,iAF_E ), &
               uAF(:,iX1,iX2,iX3,iAF_Me), uAF(:,iX1,iX2,iX3,iAF_Mp), &
               uAF(:,iX1,iX2,iX3,iAF_Mn), uAF(:,iX1,iX2,iX3,iAF_Xp), &
               uAF(:,iX1,iX2,iX3,iAF_Xn), uAF(:,iX1,iX2,iX3,iAF_Xa), &
               uAF(:,iX1,iX2,iX3,iAF_Xh), uAF(:,iX1,iX2,iX3,iAF_Gm) )  
                                                                       
    END DO                  
    END DO   
    END DO

  END SUBROUTINE InitializeSingleZoneFields_Decoherence

  SUBROUTINE InitializeRadiationFields_Decoherence

    INTEGER  :: iE, iX1, iX2, iX3, iS
    INTEGER  :: iNode, iNodeE
    INTEGER  :: iN_E
    INTEGER  :: nZ(4), nX(3)
    
    REAL(DP) :: kT(nDOF)
    REAL(DP) :: Mnu(nDOF), E
    REAL(DP) :: Gm_dd_11(nDOF) 
    REAL(DP) :: Gm_dd_22(nDOF) 
    REAL(DP) :: Gm_dd_33(nDOF) 
   
    ! -- Create Array With Flattened Energies --- !
    nZ = iZ_E0 - iZ_B0 + 1
    nX = nZ(2:4)
    nE_G = nNodesZ(1) * nZ(1)
    nX_G = nDOFX * PRODUCT( nX )

    ALLOCATE( Energies( nE_G ) )

    DO iN_E = 1, nE_G

      iE     = MOD( (iN_E-1) / nDOFE, nZ(1) ) + iE_B0
      iNodeE = MOD( (iN_E-1)        , nDOFE ) + 1

      Energies(iN_E) = NodeCoordinate( MeshE, iE, iNodeE )

    END DO

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iS = 1, nSpecies
  
      Gm_dd_11 &                                        
        = uGF(NodeNumbersX,iX1,iX2,iX3,iGF_Gm_dd_11)    
                                                        
      Gm_dd_22 &                                        
        = uGF(NodeNumbersX,iX1,iX2,iX3,iGF_Gm_dd_22)    
                                                        
      Gm_dd_33 &                                        
        = uGF(NodeNumbersX,iX1,iX2,iX3,iGF_Gm_dd_33)
  
      kT = BoltzmannConstant &
             * uAF(NodeNumbersX,iX1,iX2,iX3,iAF_T)
  
      IF( iS .EQ. iNuEE )THEN
  
        Mnu = + ( uAF(NodeNumbersX,iX1,iX2,iX3,iAF_Me) &
                  + uAF(NodeNumbersX,iX1,iX2,iX3,iAF_Mp) &
                  - uAF(NodeNumbersX,iX1,iX2,iX3,iAF_Mn) )
  
      ELSEIF( iS .EQ. iNuEE_Bar )THEN
  
        Mnu = - ( uAF(NodeNumbersX,iX1,iX2,iX3,iAF_Me) &
                  + uAF(NodeNumbersX,iX1,iX2,iX3,iAF_Mp) &
                  - uAF(NodeNumbersX,iX1,iX2,iX3,iAF_Mn) )
      
      ELSE
  
        Mnu(:) = Zero
      
      END IF
  
      DO iE = iE_B0, iE_E0
        DO iNode = 1, nDOF
  
          iNodeE = NodeNumberTable(1,iNode)
  
          E = NodeCoordinate( MeshE, iE, iNodeE )
          
          IF ( (iS .EQ. iNuEE    ) .OR. (iS .EQ. iNuXX    ) .OR. &
               (iS .EQ. iNuEE_Bar) .OR. (iS .EQ. iNuXX_Bar) ) THEN
        
            uPR(iNode,iE,iX1,iX2,iX3,iPR_D,iS) &
              = MAX( One / ( EXP( (E-Mnu(iNode))/kT(iNode) ) + One ), 1.0d-99 )
            
            
          END IF
          
        END DO
      END DO

    END DO

    END DO
    END DO
    END DO


    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iE  = iE_B0, iE_E0
    DO iNode = 1, nDOF

      iNodeE = NodeNumberTable(1,iNode)

      E = NodeCoordinate( MeshE, iE, iNodeE )
      
      ! --- Now take care of non-diagonal terms --- !

      uPR(iNode,iE,iX1,iX2,iX3,iPR_D,iNuEX) &
          = SetNonDiagonalEl( &
              uPR(iNode,iE,iX1,iX2,iX3,iPR_D,iNuEE), &
              uPR(iNode,iE,iX1,iX2,iX3,iPR_D,iNuXX), &
              Mixing_Optional = One )

      uPR(iNode,iE,iX1,iX2,iX3,iPR_D,iNuEX_Bar) &
          = SetNonDiagonalEl( &
              uPR(iNode,iE,iX1,iX2,iX3,iPR_D,iNuEE_Bar), &
              uPR(iNode,iE,iX1,iX2,iX3,iPR_D,iNuXX_Bar), &
              Mixing_Optional = One )

      DO iS = 1,nSpecies

        ! --- Set First moment to zero --- !
        
        uPR(iNode,iE,iX1,iX2,iX3,iPR_I1,iS) &
           = Zero

        uPR(iNode,iE,iX1,iX2,iX3,iPR_I2,iS) &
           = Zero

        uPR(iNode,iE,iX1,iX2,iX3,iPR_I3,iS) &
           = Zero

        ! --- Compute conserved moments --- !

        CALL ComputeConserved_TwoMoment &
                 ( uPR(:,iE,iX1,iX2,iX3,iPR_D, iS), &
                   uPR(:,iE,iX1,iX2,iX3,iPR_I1,iS), &
                   uPR(:,iE,iX1,iX2,iX3,iPR_I2,iS), &
                   uPR(:,iE,iX1,iX2,iX3,iPR_I3,iS), &
                   uCR(:,iE,iX1,iX2,iX3,iCR_N, iS), &
                   uCR(:,iE,iX1,iX2,iX3,iCR_G1,iS), &
                   uCR(:,iE,iX1,iX2,iX3,iCR_G2,iS), &
                   uCR(:,iE,iX1,iX2,iX3,iCR_G3,iS), &
                   Gm_dd_11, Gm_dd_22, Gm_dd_33 )

      END DO

    END DO
    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE InitializeRadiationFields_Decoherence

  SUBROUTINE InitializeOscillations

    INTEGER :: k

    ALLOCATE( fMatrixOsc(nM,nE_G,nF,nF) )
    ALLOCATE( fMatrixOld(nX_G,nM,nE_G,nF,nF) )
    ALLOCATE( Enu(nE_G) )  
    ALLOCATE( fImaginaryPart(nX_G,nE_G,nM) )

    fImaginaryPart = Zero
    Enu = Energies / MeV
    TimeLocOutput = Zero    

    ! Initialize array for RK
    DO k = 1,nRK-1

      AA(1,k) = Zero
      AA(2,k) = b2(k)
      AA(3,k) = b3(k)
      AA(4,k) = b4(k)
      AA(5,k) = b5(k)
      AA(6,k) = b6(k)

    END DO

    CALL InitializeMatricesOscillations

  END SUBROUTINE InitializeOscillations

  SUBROUTINE InitializeMatricesOscillations
    
    INTEGER :: x
    INTEGER :: jF, iE, m

    REAL(DP), PARAMETER :: eV_to_erg = 1.60218d-12
    REAL(DP), PARAMETER :: clite = 2.99792458d10 !cm/s
    REAL(DP), PARAMETER :: m1 = 0.0d0 * eV_to_erg**2/clite**4 ! g,
    REAL(DP), PARAMETER :: dm21 = 2.43e-3 * eV_to_erg**2/clite**4 !g
    REAL(DP), PARAMETER :: theta12V = 9.0d0 * Pi/180.0d0 ! rad
    REAL(DP), PARAMETER :: c12V = COS(theta12V)
    REAL(DP), PARAMETER :: s12V = SIN(theta12V)
    REAL(DP), PARAMETER :: SIN2thetaW = 0.23122

    REAL(DP) :: kV(nE_G,nF)
    REAL(DP) :: Hvm(nF,nF)

    ALLOCATE( IdentityR(nF,nF) )
    ALLOCATE( IdentityC(nF,nF) )
    ALLOCATE( YIdentity(nS,nY) )

    ALLOCATE( Hvf(nM,nE_G,nF,nF) )

    ALLOCATE( AV(nE_G,nF,nF)    )
    ALLOCATE( CV(nE_G,nF,nF,nF) )

    DO x = 1,nS

      YIdentity(x,:) = (/ Pi/2.0d0, Pi/2.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0 /)

    END DO

    IdentityR(:,:) = Zero
    IdentityC(:,:) = Zero

    DO jF = 1,nF

      IdentityR(jF,jF) = One
      IdentityC(jF,jF) = One

    END DO

    PMNS(1,1) = c12V
    PMNS(1,2) = s12V
    PMNS(2,1) = -s12V
    PMNS(2,2) = c12V

    !Initialize vector of mass Eigenvalues
    DO iE = 1,nE_G

      kV(iE,1) = m1**2* clite**4 /(Two * Enu(iE)*1e6*eV_to_erg)
      kV(iE,2) = (kV(iE,1) + dm21) * clite**4 /(Two * Enu(iE)*1e6*ev_to_erg)

    END DO

    !Initialize Hamiltonians
    Hvm(:,:) = Zero

    DO iE = 1,nE_G

      !Vacuum Hamiltonian, constant, no need to reevaluate it
      Hvm(1,1) = kV(iE,1)
      Hvm(2,2) = kV(iE,2)
      Hvf(1,iE,:,:) = Matmul(Matmul(PMNS,Hvm(:,:)),Conjg(Transpose(PMNS)))
      Hvf(2,iE,:,:) = Conjg(Hvf(1,iE,:,:))

    END DO

    DO m = 1,nM
      DO iE = 1,nE_G

        CV(iE,:,:,:) = CofactorMatrix(Hvf(m,iE,:,:), kV(iE,:))
        AV(iE,:,:)   = Evaluate_AV    (Hvf(m,iE,:,:), kV(iE,:))

      END DO
    END DO

  END SUBROUTINE InitializeMatricesOscillations

  FUNCTION fermidirac(energy,mu,T) result(fd)

    REAL(KIND=8), INTENT(IN) :: energy,mu,T !MeV
    REAL(KIND=8) :: fd

    fd = 1.0d0/(1.0d0 + exp((energy - mu)/T))

  END FUNCTION fermidirac


  FUNCTION CofactorMatrix(M,Lambda) RESULT(CM)

    COMPLEX(DP), INTENT(IN)  :: M(2,2)
    REAL(DP),    INTENT(IN)  :: Lambda(2)

    COMPLEX(DP) :: CM(2,2,2)

    !local
    INTEGER :: j

    DO j=1,2

        CM(j,1,1) = M(2,2) - Lambda(j)
        CM(j,1,2) = -M(2,1)
        CM(j,2,1) = CONJG(CM(j,1,2))
        CM(j,2,2) = M(1,1) - Lambda(j)

    END DO

  END FUNCTION CofactorMatrix

  FUNCTION Evaluate_AV(M,Lambda) RESULT(AV_loc)

    COMPLEX(DP), INTENT(IN)  :: M(2,2)
    REAL(DP),    INTENT(IN)  :: Lambda(2)

    COMPLEX(DP) :: Av_loc(2,2)

    !local
    INTEGER :: j
    REAL(DP) :: Delta
    REAL(DP) :: re2, rmu2
    COMPLEX(DP) :: CM(2,2,2)

    CM = CofactorMatrix(M,Lambda)

    DO j=1,2

        IF(j == 1) Delta = ( lambda(2) - lambda(1) )
        IF(j == 2) Delta = ( lambda(1) - lambda(2) )

        re2  = Delta * REAL( CM(j,1,1) )
        rmu2 = Delta * REAL( CM(j,2,2) )

        IF( ABS(PMNS(1,j)) > ABS(PMNS(2,j)) ) THEN

            AV_loc(j,1)=REAL(PMNS(1,j) * SQRT(re2) / CM(j,1,1))
            AV_loc(j,2)=REAL(PMNS(1,j) * SQRT(rmu2) / CM(j,2,1))

        END IF

        IF( ABS(PMNS(2,j)) > ABS(PMNS(1,j)) ) THEN

            AV_loc(j,1)=REAL(PMNS(2,j)*SQRT(re2) / CM(j,1,2))
            AV_loc(j,2)=REAL(PMNS(2,j)*SQRT(rmu2) / CM(j,2,2))

        END IF
    END DO

  END FUNCTION Evaluate_AV


  REAL(DP) FUNCTION SetNonDiagonalEl( El_11,El_22,  &
                                     Mixing_Optional )
    
    REAL(DP), INTENT(IN) :: El_11, El_22
    REAL(DP), INTENT(IN), OPTIONAL :: Mixing_Optional

    REAL(DP) :: Mixing
    REAL(DP) :: Lmax
    REAL(DP) :: El_z, El_x, El_12
    REAL(DP) :: Trace
    
    Mixing = One
    IF( PRESENT(Mixing_Optional) ) &
      Mixing = Mixing_Optional
    
    Trace = El_11 + El_22
    Lmax = MIN( Half * Trace, One - Half * Trace )
    El_z = Half*( El_11 - El_22 )

    El_x = SQRT(Lmax**2 - El_z**2)

    SetNonDiagonalEl = El_x * Mixing

  END FUNCTION SetNonDiagonalEl

END MODULE InitializationModule

PROGRAM Decoherence

  USE KindModule, ONLY: &
    DP, Zero, One, Two, &
    Pi, TwoPi, &
    Half, Five
  USE UnitsModule, ONLY: &
    Kilometer, MeV, Millisecond, &
    Second
  USE ProgramHeaderModule, ONLY: &
    nZ, nX, nNodesZ, nDOFX, nDOFE, &
    iX_B0, iX_E0, iX_B1, iX_E1, &
    iE_B0, iE_E0, iE_B1, iE_E1, &
    iZ_B0, iZ_E0, iZ_B1, iZ_E1
  USE ProgramInitializationModule, ONLY: &
    InitializeProgram, &
    FinalizeProgram
  USE TimersModule, ONLY: &
    InitializeTimers, &
    FinalizeTimers, &
    TimersStart, &
    TimersStop, &
    Timer_InputOutput, &
    Timer_Evolve
  USE ReferenceElementModuleX, ONLY: &
    InitializeReferenceElementX, &
    FinalizeReferenceElementX
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    InitializeReferenceElementX_Lagrange, &
    FinalizeReferenceElementX_Lagrange
  USE ReferenceElementModuleE, ONLY: &
    InitializeReferenceElementE, &
    FinalizeReferenceElementE
  USE ReferenceElementModuleE_Lagrange, ONLY: &
    InitializeReferenceElementE_Lagrange, &
    FinalizeReferenceElementE_Lagrange
  USE ReferenceElementModule, ONLY: &
    InitializeReferenceElement, &
    FinalizeReferenceElement
  USE ReferenceElementModule_Lagrange, ONLY: &
    InitializeReferenceElement_Lagrange, &
    FinalizeReferenceElement_Lagrange
  USE GeometryComputationModule, ONLY: &
    ComputeGeometryX
  USE GeometryComputationModuleE, ONLY: &
    ComputeGeometryE
  USE GeometryFieldsModule, ONLY: &
    uGF
  USE GeometryFieldsModuleE, ONLY: &
    uGE
  USE FluidFieldsModule, ONLY: &
    uCF, uPF, uAF, &
    iPF_D, iAF_T, iAF_Ye
  USE RadiationFieldsModule, ONLY: &
    uCR, uPR, iPR_D
  USE InputOutputModuleHDF, ONLY: &
    WriteFieldsHDF, FileNumber
  USE EquationOfStateModule, ONLY: &
    InitializeEquationOfState, &
    FinalizeEquationOfState
  USE OpacityModule_TABLE, ONLY: &
    InitializeOpacities_TABLE, &
    FinalizeOpacities_TABLE
  USE NeutrinoOpacitiesModule, ONLY: &
    CreateNeutrinoOpacities, &
    DestroyNeutrinoOpacities
  USE NeutrinoOpacitiesComputationModule, ONLY: &
    ComputeNeutrinoOpacities  
  USE TwoMoment_ClosureModule, ONLY: &
    InitializeClosure_TwoMoment
  USE IntegrationModule, ONLY: &
    TimeStepping_Decoherence
  USE CollisionsModule, ONLY: &
    CollisionsDriver, &
    EvaluateImpact
  USE InitializationModule, ONLY: &
    InitializeFields_Decoherence, &
    InitializeSingleZone_Decoherence, &
    InitializeOscillations, &
    fMatrixOsc, fMatrixOld, &
    nF, nE_G, nX_G, &
    fImaginaryPart
  USE OscillationsModule, ONLY: &
    OscillationsDriver, &
    Packf, Flattenf
  USE InputOutputDecoherenceModule, ONLY: &
    OutputImaginaryPart, &
    InitializeRestart
  
  IMPLICIT NONE

  LOGICAL  :: wrt, Restart
  INTEGER  :: RestartNumber
  INTEGER  :: iCycle, iCycleD, iCycleW
  INTEGER  :: wrt_control
  INTEGER  :: nNodes, nSpecies
  INTEGER  :: nE, bcE, bcX(3)
  INTEGER  :: swE, swX(3)
  INTEGER  :: iX, iX1, iX2, iX3, iNodeX
  INTEGER  :: iS

  REAL(DP) :: eL, eR, xL(3), xR(3)
  REAL(DP) :: t, dt_block, t_end
  REAL(DP) :: dt_max
  REAL(DP) :: dt, dt_osc, dt_coll
  REAL(DP) :: dt_rand, Random, Exp_Rand

  Restart = .FALSE.
  RestartNumber = 493
  nNodes = 1
  nSpecies = 6
  nF = 2

  nX  = [ 2, 1, 1 ]
  xL  = [ 0.0d0 * Kilometer, 0.0_DP, 0.0_DP ]
  xR  = [ 1.0d0 * Kilometer, Pi,     TwoPi  ]
  bcX = [ 1, 0, 0 ]
  swX = [ 1, 1, 1 ]

  nE  = 50
  eL  = 1.0d0  * MeV
  eR  = 101.0d0 * MeV
  bcE = 0
  swE = 0

  t        = 0.0_DP
  t_end    = 1.0d0    * Second
  dt_block = 1.0d-13  * Second
  
  dt_osc   = dt_block
  dt_coll  = dt_block

  iCycleD = 1
  iCycleW = 1

  CALL InitializeProgram &
         ( ProgramName_Option &
             = 'Decoherence', &
           nX_Option &
             = nX, &
           swX_Option &
             = swX, &
           bcX_Option &
             = bcX, &
           xL_Option &
             = xL, &
           xR_Option &
             = xR, &
           nE_Option &
             = nE, &
           swE_Option &
             = swE, &
           bcE_Option &
             = bcE, &
           eL_Option &
             = eL, &
           eR_Option &
             = eR, &
           nNodes_Option &
             = nNodes, &
           CoordinateSystem_Option &
             = 'SPHERICAL', &
           ActivateUnits_Option &
             = .TRUE., &
           nSpecies_Option &
             = nSpecies, &
           BasicInitialization_Option &
             = .TRUE. )

  ! --- Position Space Reference Element and Geometry ---

  CALL InitializeReferenceElementX

  CALL InitializeReferenceElementX_Lagrange

  CALL ComputeGeometryX &
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF )

  ! --- Energy Space Reference Element and Geometry ---

  CALL InitializeReferenceElementE

  CALL InitializeReferenceElementE_Lagrange

  CALL ComputeGeometryE &
         ( iE_B0, iE_E0, iE_B1, iE_E1, uGE )

  ! --- Phase Space Reference Element ---

  CALL InitializeReferenceElement

  CALL InitializeReferenceElement_Lagrange

  ! --- Initialize Equation of State ---

  CALL InitializeEquationOfState &
         ( EquationOfState_Option &
             = 'TABLE', &
           EquationOfStateTableName_Option &
             = 'EquationOfStateTable.h5' )
  
! --- Initialize and create Opacities

  CALL InitializeOpacities_TABLE &                         
         ( OpacityTableName_EmAb_Option &                  
             = 'wl-Op-SFHo-15-25-50-E40-B85-AbEm.h5', &    
           OpacityTableName_Iso_Option  &                  
             = 'wl-Op-SFHo-15-25-50-E40-B85-Iso.h5', &     
!           OpacityTableName_NES_Option &                   
!             = 'wl-Op-SFHo-15-25-50-E40-B85-NES.h5', &     
!           OpacityTableName_Pair_Option &                  
!             = 'wl-Op-SFHo-15-25-50-E40-B85-Pair.h5', &    
           Verbose_Option = .TRUE. )                       

  CALL CreateNeutrinoOpacities( nZ, nNodesZ, nF )

  ! --- Initialize Moment Closure ---

  CALL InitializeClosure_TwoMoment

  ! --- Create Profile and initialize Moments ---
  
  !CALL InitializeFields_Decoherence
  
  CALL InitializeSingleZone_Decoherence
    
  CALL InitializeOscillations

  IF( Restart ) THEN
    
    CALL InitializeRestart( RestartNumber, &
        t, nX, swX, nE, swE )

    FileNumber = RestartNumber + 1

  ELSE 
    
    t = 0.0_DP

  END IF

  CALL WriteFieldsHDF &
         ( Time = t, &
           WriteGF_Option = .TRUE., &
           WriteFF_Option = .TRUE., &
           WriteRF_Option = .TRUE. )
  
  CALL OutputImaginaryPart( Time = t )

  ! --- Evolve ---

  WRITE(*,*)
  WRITE(*,'(A6,A,ES8.2E2,A8,ES8.2E2)') &
    '', 'Evolving from t = ', t / Millisecond, &
    ' to t = ', t_end / Millisecond
  WRITE(*,*)

  WRITE(*,'(A6,A5,A3,A10,A2,3A13)') '','Cycle','', 'time (s)','', &
      'dt_block (s)', 'dt_osc (s)', 'dt_coll (s)' 
 
  iCycle = 0
  wrt     = .FALSE.
  wrt_control = 1
  dt_rand = dt_block
  DO WHILE( t < t_end )
    
    iCycle = iCycle + 1
    
    IF ( MOD( iCycle, iCycleW ) == 0 ) THEN
      
      wrt = .TRUE.

    ENDIF

    IF( MOD( iCycle, iCycleD ) == 0 ) THEN
        
      WRITE(*,'(A4,I8.8,4E13.5)') '',iCycle, t / Second, dt_rand / Second, &
                  dt_osc / Second, dt_coll / Second

    END IF

    CALL TimersStart( Timer_Evolve )

    ! --- Generate Random Time Step --- !
    CALL RANDOM_NUMBER( Random )
    Exp_Rand = - LOG( Random )
    dt_rand = dt_block * MIN( Five, Exp_Rand )
    
    IF( t + dt_rand > t_end ) THEN

      dt_rand = t_end - t

    END IF

    dt = dt_osc
     
    DO iX = 1,nX_G
        
        iX3    = MOD( (iX-1) / ( nDOFX * nX(1) * nX(2) ), nX(3) ) + iX_B0(3)
        iX2    = MOD( (iX-1) / ( nDOFX * nX(1)         ), nX(2) ) + iX_B0(2)
        iX1    = MOD( (iX-1) / ( nDOFX                 ), nX(1) ) + iX_B0(1)
        iNodeX = MOD( (iX-1)                            , nDOFX ) + 1 

        ! --- When you pack, add back Imaginary Part --- !
        CALL Packf( iX, iNodeX, iX1, iX2, iX3 )
        
        CALL OscillationsDriver( dt_rand, dt, dt_osc, & 
            uPF(iNodeX,iX1,iX2,iX3,iPF_D ), &
            uAF(iNodeX,iX1,iX2,iX3,iAF_T ), &
            uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) )
        
        fMatrixOld(iX,:,:,:,:) = fMatrixOsc(:,:,:,:)
        
        fImaginaryPart(iX,:,1) = AIMAG( fMatrixOsc(1,:,1,2) )

        fImaginaryPart(iX,:,2) = AIMAG( fMatrixOsc(2,:,1,2) )
        
        CALL Flattenf( iNodeX, iX1, iX2, iX3 )

    END DO

    dt = dt_coll
    
    CALL CollisionsDriver( dt_rand, dt, dt_coll)

    CALL EvaluateImpact( dt_rand, dt_block )

    t = t + dt_rand

    CALL TimersStop( Timer_Evolve )

    IF( wrt )THEN

      CALL TimersStart( Timer_InputOutput )
      
      CALL WriteFieldsHDF &
             ( Time = t, &
               WriteGF_Option = .FALSE., &
               WriteFF_Option = .FALSE., &
               WriteRF_Option = .TRUE., &
               WriteOP_Option = .FALSE. )

      CALL OutputImaginaryPart( Time = t )

      CALL TimersStop( Timer_InputOutput )

      wrt = .FALSE.

    END IF

  END DO

  ! --- Write Final Solution ---

  CALL TimersStart( Timer_InputOutput )

  CALL WriteFieldsHDF &
         ( Time = t, &
           WriteGF_Option = .TRUE., &
           WriteFF_Option = .TRUE., &
           WriteRF_Option = .TRUE., &
           WriteOP_Option = .TRUE. )
  
  CALL OutputImaginaryPart( Time = t )

  CALL TimersStop( Timer_InputOutput )

  WRITE(*,*)
  WRITE(*,'(A6,A,I6.6,A,ES12.6E2,A)') &
    '', 'Finished ', iCycle, ' Cycles in ', Timer_Evolve, ' s'
  WRITE(*,*)

  CALL FinalizeTimers

  CALL FinalizeEquationOfState

  CALL FinalizeOpacities_TABLE
    
  CALL DestroyNeutrinoOpacities

  CALL FinalizeReferenceElementX

  CALL FinalizeReferenceElementX_Lagrange

  CALL FinalizeReferenceElementE

  CALL FinalizeReferenceElementE_Lagrange

  CALL FinalizeReferenceElement

  CALL FinalizeReferenceElement_Lagrange

  CALL FinalizeProgram

END PROGRAM Decoherence

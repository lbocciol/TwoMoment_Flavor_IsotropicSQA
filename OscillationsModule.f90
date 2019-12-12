MODULE OscillationsModule
  
  USE KindModule, ONLY: &
    DP, One, Zero, Five, &
    Two, Pi, TwoPi, FourPi
  USE UnitsModule, ONLY: &
    MeV, Erg, Gram, Centimeter, &
    Second
  USE ProgramHeaderModule, ONLY: &
    nDOFE, nE, iE_B0
  USE RadiationFieldsModule, ONLY: &
    nSpecies, uPR, iPR_D
  USE MeshModule, ONLY: &
    MeshE, NodeCoordinate
  USE ReferenceElementModuleE, ONLY: &
    WeightsE
  USE InitializationModule, ONLY: &
    BB5, BB6, AA, nRK, nRKOrder
  USE InitializationModule, ONLY: &
    YIdentity, IdentityC, Hvf, &
    Enu, fMatrixOsc, &
    fImaginaryPart, &
    nF, nM, nY, nS, nX_G, nE_G, &
    TimeLocOutput
  USE FlavorOpacitiesModule, ONLY: &
    PackIntoMatrix, FlattenMatrix
  USE OscillationsUtilsModule, ONLY: &
    Im, B, W, JInverse, &
    Eigenvalues, EigenvectorMatrix

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: OscillationsDriver
  PUBLIC :: Packf
  PUBLIC :: Flattenf

  ! --- HAMILTONIANS & CO. --- !
  COMPLEX(DP), ALLOCATABLE :: Hmsw(:,:,:,:) !MSW Hamiltonian in flavor basis
  COMPLEX(DP), ALLOCATABLE :: S(:,:,:,:) !S matrix
  COMPLEX(DP), ALLOCATABLE :: Hf(:,:,:,:)
  COMPLEX(DP), ALLOCATABLE :: U0(:,:,:,:)
  COMPLEX(DP), ALLOCATABLE :: pm0(:,:,:,:)

  REAL(DP), ALLOCATABLE :: k0(:,:,:)

  ! --- CONSTANTS USED IN OSCILLATION CALCULATIOnS --- ! 
  REAL(DP), PARAMETER :: eV_to_erg = 1.60218d-12 !erg
  REAL(DP), PARAMETER :: hbar = 1.05457266d-27 ! erg s
  REAL(DP), PARAMETER :: clite = 2.99792458d10 !cm/s
  REAL(DP), PARAMETER :: hbarc = hbar*clite
  REAL(DP), PARAMETER :: GF = 1.1663787d-5/ (1e9*1e9*eV_to_erg**2) &
                            * hbarc**3 !erg cm^3
  REAL(DP), PARAMETER :: Mp = 1.6726219e-24 ! g
  
  ! --- Local Variables --- !
  REAL(DP)              :: Rho, Temp, Ye
  REAL(DP)              :: TimeBlock, dt_loc

CONTAINS

  SUBROUTINE OscillationsDriver( TimeBlockEnd_in, dt, &
             dt_osc, Rho_in, Temp_in, Ye_in )

    REAL(DP), INTENT(IN)    :: TimeBlockEnd_in, dt
    REAL(DP), INTENT(IN)    :: Temp_in, Rho_in, Ye_in
    REAL(DP), INTENT(INOUT) :: dt_osc 

    LOGICAL  :: BlockDone 
    REAL(DP) :: TimeBlockEnd
    
    TimeBlockEnd = TimeBlockEnd_in / Second

    Rho      = Rho_in / Gram * Centimeter**3
    Temp     = Temp_in / MeV
    Ye       = Ye_in

    CALL InitializeArrays

    Hmsw(:,:,:,:) = Zero
    Hmsw(1,:,1,1) = SQRT( Two ) * GF * Rho * Ye / Mp
    Hmsw(2,:,1,1) = - CONJG( Hmsw(1,:,1,1) )

    Hf(:,:,:,:) = Hvf(:,:,:,:) + Hmsw(:,:,:,:)

    CALL Diagonalize
    
    BlockDone = .FALSE.
    TimeBlock = Zero
    dt_loc = dt / Second
    
    DO WHILE( BlockDone .EQV. .FALSE. )
      
      IF( TimeBlock + dt_loc >= TimeBlockEnd ) THEN
 
        dt_loc = TimeBlockEnd - TimeBlock
        BlockDone = .TRUE.

      END IF

      CALL EvolveOscillations( BLockDone )  
      
      IF( .NOT. BlockDone ) dt_osc = dt_loc * Second

    END DO

    !CALL OutputMatrixOsc( TimeBlockEnd )
    
    CALL FinalizeArrays

        
  END SUBROUTINE OscillationsDriver

  SUBROUTINE EvolveOscillations( BlockDone )

    LOGICAL, INTENT(IN) :: BlockDone

    LOGICAL     :: Reloop
    REAL(DP)    :: Yerror, MaxError, TimeTmp
    REAL(DP)    :: Ks(nRK,nM,nE_G,nS,nY)
    REAL(DP)    :: Ytmp(nM,nE_G,nS,nY), Y(nM,nE_G,nS,nY)
    COMPLEX(DP) :: SMSW(nF,nF), SSI(nF,nF), SNew(nF,nF)
    
    INTEGER :: k,l,m,iE,iY,x

    REAL(DP) :: Accuracy = 1d-12
    REAL(DP) :: Increase = 1.1

    Ks(:,:,:,:,:) = Zero

    DO m = 1,nM
      DO iE = 1,nE_G

        Y(m,iE,:,:) = YIdentity

      END DO
    END DO

    CALL Get_P
    Reloop = .true.

    DO WHILE(Reloop .eqv. .true.)

      DO k = 1,nRK

        Ytmp = Y
        DO l = 1,k-1
          !$OMP PARALLEL DO COLLAPSE(4) &
          !$OMP DEFAULT(SHARED) PRIVATE( m,iE,x,iY )
          DO m = 1,nM
            DO iE = 1,nE_G
              DO x = 1,nS
                DO iY = 1,nY

                  Ytmp(m,iE,x,iY) = Ytmp(m,iE,x,iY) + AA(k,l) * Ks(l,m,iE,x,iY)
        
                END DO
              END DO
            END DO
          END DO
          !$OMP END PARALLEL DO
        END DO

        CALL RK_step( Ytmp, Ks(k,:,:,:,:) )

      END DO

      ! Do the rest of the RK and evaluate error
      MaxError = 0.0d0
      Ytmp = Y

      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(m,iE,x,iY,Yerror)
      !$OMP DO COLLAPSE(2) REDUCTION(Max:MaxError)
      DO m = 1,nM
        DO iE = 1,nE_G
          DO x=1,nS
            DO iY=1,nY

              Yerror = 0.0d0

              DO k=1,nRK

                Ytmp(m,iE,x,iY) = Ytmp(m,iE,x,iY) + BB5(k) * Ks(k,m,iE,x,iY)
                Yerror = Yerror + (BB5(k) - BB6(k)) * Ks(k,m,iE,x,iY)

              END DO

              MaxError = MAX(MaxError,ABS(Yerror))

            END DO
          END DO
        END DO
      END DO
      !$OMP END DO
      !$OMP END PARALLEL
        
      TimeTmp = TimeBlock + dt_loc
      
      IF(MaxError .gt. Accuracy) THEN

        dt_loc = dt_loc * 0.9 *( (Accuracy/MaxError)**(1.0d0/(nRKOrder-1.0d0)) )
        Reloop=.true.
      
      ELSE

        dt_loc = dt_loc * Increase
        Reloop = .false.

        IF(MaxError .gt. 0.0d0) &
          dt_loc = dt_loc * MIN( 1.0, ((Accuracy/MaxError)**(1.0d0/5.0d0))/Increase)
          Y = Ytmp

      END IF
    END DO

    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(m,iE,SMSW,SSI,SNew)
    !$OMP DO COLLAPSE(2) 
    DO m = 1,nM
      DO iE = 1,nE_G

        SMSW = MATMUL(W(Y(m,iE,1,:)),B(Y(m,iE,1,:)))
        SSI  = MATMUL(W(Y(m,iE,2,:)),B(Y(m,iE,2,:)))

        !make sure that the matrices are CLOSE to diagonal
        IF(ABS(SMSW(1,1))**2 + 0.1 < ABS(SMSW(1,2))**2 .OR. &
           ABS(SSI(1,1))**2  + 0.1 < ABS(SSI(1,2))**2  .OR. &
           BlockDone ) THEN

           SNew = MATMUL(SMSW,SSI)
           S(m,iE,:,:) = MATMUL(MATMUL(U0(m,iE,:,:),SNew), &
                         CONJG(TRANSPOSE(U0(m,iE,:,:))))
           fMatrixOsc(m,iE,:,:) =  &
               MATMUL(MATMUL(S(m,iE,:,:),fMatrixOsc(m,iE,:,:)), &
                         CONJG(TRANSPOSE((S(m,iE,:,:)))))
           Y(m,iE,:,:) = YIdentity !not sure it's necessary

        ELSE

           Y(m,iE,1,3) = MOD(Y(m,iE,1,3),TwoPi)
           Y(m,iE,1,5) = MOD(Y(m,iE,1,5),1.0)
           Y(m,iE,1,6) = MOD(Y(m,iE,1,6),1.0)

           Y(m,iE,2,3) = MOD(Y(m,iE,2,3),TwoPi)
           Y(m,iE,2,5) = MOD(Y(m,iE,2,5),1.0)
           Y(m,iE,2,6) = MOD(Y(m,iE,2,6),1.0)

        END IF
      END DO
    END DO
    !$OMP END DO
    !$OMP END PARALLEL

    TimeBlock = TimeTmp
    
  END SUBROUTINE EvolveOscillations

  SUBROUTINE RK_step( Y, K )

    REAL(DP),INTENT(IN) :: Y(nM,nE_G,nS,nY)

    REAL(DP),INTENT(INOUT) :: K(nM,nE_G,nS,nY)
    !it's the integrand of the SI potential in the mass basis

    COMPLEX(DP) :: VfSI(nM,NF,NF)
    COMPLEX(DP) :: VfSIE(NF,NF) !it's the integrand of the SI potential in
    !the flavor basis 
    COMPLEX(DP) :: Sfm(NF,NF) !it transforms S to f to calculate the potential,
    !and THEN it transforms again with U0 to go to the flavor basis
    COMPLEX(DP) :: UWBW(nM,nE_G,NF,NF)
    COMPLEX(DP) :: BSI(nM,nE_G,NF,NF) !S_{si} = W(Y_{si}) * B(Y_{si}), so BSI is
    !one of the two contributions to S_{si}
    INTEGER :: m,iE,iY,l
 
    !Now variables that I'm not sure what they are
    COMPLEX(DP) :: Ha(NF,NF)
    COMPLEX(DP) :: HB0,HB1
    REAL(DP) :: dvdr(4)
    REAL(DP) :: JI(3,4)

    VfSI(:,:,:) = Zero
    
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(m,iE,Sfm,VfSIE)
    !$OMP DO COLLAPSE(2)
    DO m = 1,nM
      
      DO iE = 1,nE_G
        
        !setup transformation 
        UWBW(m,iE,:,:) = MATMUL(MATMUL(U0(m,iE,:,:)  , W(Y(m,iE,1,:))), &
                               MATMUL(B(Y(m,iE,1,:)), W(Y(m,iE,2,:))))

        !For 1=msw only 4 and 5 are relevant (i.3. B(Y(msw)) = I)
        DO iY = 1,4
        
          K(m,iE,1,iY) = 0.0d0
        
        END DO
        
        K(m,iE,1,5) = k0(m,iE,1)*dt_loc/hbar/TwoPi
        K(m,iE,1,6) = k0(m,iE,2)*dt_loc/hbar/TwoPi
        !DONE WITH MSW!

        BSI(m,iE,:,:) = B(Y(m,iE,2,:))
        
        !Evaluate Si potential
        Sfm = MATMUL(UWBW(m,iE,:,:),BSI(m,iE,:,:))
        VfSIE = MATMUL(MATMUL(Sfm,pm0(m,iE,:,:)),CONJG(TRANSPOSE(Sfm)))

        !Below is the "MINus" part of eq. 9 from Richers 2019
        IF(m .eq. 2) VfSIE = - CONJG(VfSIE)
        
        !$OMP CRITICAL        
        VfSI(1,:,:) = VfSI(1,:,:) + VfSIE(:,:)
        !$OMP END CRITICAL
      
      END DO
    
    END DO
    !$OMP END DO
    !$OMP END PARALLEL

    VfSI(2,:,:) = - CONJG(VfSI(1,:,:))
    !VfSI(:,:,:) = 0.0d0 

    !Now the part that I DOn't REALly understand
   
    !$OMP PARALLEL DEFAULT(SHARED) &
    !$OMP PRIVATE(m,iE,Ha,HB0,HB1,dvdr,JI) 
    !$OMP DO COLLAPSE(2) PRIVATE(m,iE,Ha,HB0,HB1,dvdr,JI) 
    DO m = 1,nM
      
      DO iE = 1,nE_G
        
        Ha = MATMUL(MATMUL(CONJG(TRANSPOSE(UWBW(m,iE,:,:))),VfSI(m,:,:)),UWBW(m,iE,:,:))

        K(m,iE,2,5) = dt_loc*REAL(Ha(1,1))/hbar/TwoPi
        K(m,iE,2,6) = dt_loc*REAL(Ha(2,2))/hbar/TwoPi
        
        HB0 = -Im/hbar*( Ha(1,2)*BSI(m,iE,2,1) )
        HB1 = -Im/hbar*( Ha(1,2)*BSI(m,iE,2,2) )

        dvdr(1) = REAL(HB1)
        dvdr(2) = AIMAG(HB1)
        dvdr(3) = REAL(HB0)
        dvdr(4) = AIMAG(HB0)

        JI = JInverse(Y(m,iE,2,:))

        DO iY = 1,3
          
          K(m,iE,2,iY) = Zero
          
          DO l = iY,4

            K(m,iE,2,iY) = K(m,iE,2,iY) + JI(iY,l)*dvdr(l)
          
          END DO
          
          K(m,iE,2,iY) = K(m,iE,2,iY) * dt_loc
        
        END DO

        K(m,iE,2,4) = Zero
      
      END DO
    END DO
    !$OMP END DO
    !$OMP END PARALLEL

  END SUBROUTINE RK_step


  SUBROUTINE Get_P

    REAL(DP)    :: nuHz, dnuHz
    INTEGER     :: m, iN_E, iE, iNodeE

    DO m = 1,nM 
      DO iN_E = 1,nE_G
        
        iE     = MOD( (iN_E-1) / nDOFE, nE    ) + iE_B0
        iNodeE = MOD( (iN_E-1)        , nDOFE ) + 1

        dnuHz = MeshE % Width(iE) / Erg / (TwoPi*hbar)
        nuHz = Enu(iE) * 1e6 * eV_to_erg / (TwoPi*hbar)
        
        pm0(m,iN_E,:,:) = &
            MATMUL(MATMUL(CONJG(TRANSPOSE(U0(m,iN_E,:,:))), &
            fMatrixosc(m,iN_E,:,:)), &
            U0(m,iN_E,:,:)) * &
            SQRT(2.0d0) * GF * FourPi *nuHz**2 * & 
            dnuHz * WeightsE(iNodeE) / (clite**3)
       
      END DO
    END DO

  END SUBROUTINE Get_P

  SUBROUTINE Diagonalize
    
    INTEGER :: m, iE

    DO m = 1,nM
      DO iE = 1,nE_G

         k0(m,iE,:)   = Eigenvalues       ( Hf(m,iE,:,:) )
         U0(m,iE,:,:) = EigenvectorMatrix( Hf(m,iE,:,:), k0(m,iE,:), iE )

         IF(m == 2) U0(m,iE,:,:) = CONJG( U0(m,iE,:,:) )

      END DO
    END DO


  END SUBROUTINE Diagonalize

  SUBROUTINE InitializeArrays

    ALLOCATE( Hmsw(nM,nE_G,nF,nF) )
    ALLOCATE( Hf  (nM,nE_G,nF,nF) )
    ALLOCATE( U0  (nM,nE_G,nF,nF) )
    ALLOCATE( S   (nM,nE_G,nF,nF) )
    ALLOCATE( pm0 (nM,nE_G,nF,nF) )
    ALLOCATE( k0  (nM,nE_G,nF)    )

  END SUBROUTINE InitializeArrays

  SUBROUTINE FinalizeArrays

    DEALLOCATE( Hmsw, Hf, U0, S, pm0, k0 )

  END SUBROUTINE FinalizeArrays

  SUBROUTINE Packf( iX, iNodeX, iX1, iX2, iX3 )
    
    INTEGER, INTENT(IN) :: iX, iNodeX, iX1, iX2, iX3
    INTEGER :: iN_E, iE, iNodeE
    INTEGER :: nMatter, nAntiMatter

    nMatter = nSpecies / 2
    nAntiMatter = nSpecies / 2 + 1

    DO iN_E = 1,nE_G

      iE     = MOD( (iN_E-1) / nDOFE, nE ) + iE_B0
      iNodeE = MOD( (iN_E-1)        , nDOFE ) + 1
      
      ! --- MATTER PART --- !
      fMatrixOsc(1,iN_E,:,:) = PackIntoMatrix( &
          uPR(iNodeE,iE,iX1,iX2,iX3,iPR_D,1:nMatter), &
          nSpecies/2, nF )

      fMatrixOsc(1,iN_E,1,2) = fMatrixOsc(1,iN_E,1,2) + &
          Im * fImaginaryPart(iX,iN_E,1)

      fMatrixOsc(1,iN_E,2,1) = CONJG( fMatrixOsc(1,iN_E,1,2) )
    
      ! --- ANTIMATTER PART --- !
      fMatrixOsc(2,iN_E,:,:) = PackIntoMatrix( &
          uPR(iNodeE,iE,iX1,iX2,iX3,iPR_D,nAntiMatter:), &
          nSpecies/2, nF )

      fMatrixOsc(2,iN_E,1,2) = fMatrixOsc(2,iN_E,1,2) + &
          Im * fImaginaryPart(iX,iN_E,2)

      fMatrixOsc(2,iN_E,2,1) = CONJG( fMatrixOsc(2,iN_E,1,2) )

    END DO

  END SUBROUTINE Packf
    
  SUBROUTINE Flattenf( iNodeX, iX1, iX2, iX3 )

    INTEGER, INTENT(IN) :: iNodeX, iX1, iX2, iX3
    INTEGER :: iN_E, iE, iNodeE
    INTEGER :: nMatter, nAntiMatter

    nMatter = nSpecies / 2
    nAntiMatter = nSpecies / 2 + 1

    DO iN_E = 1,nE_G

      iE     = MOD( (iN_E-1) / nDOFE, nE    ) + iE_B0
      iNodeE = MOD( (iN_E-1)        , nDOFE ) + 1

      
      uPR(iNodeE,iE,iX1,iX2,iX3,iPR_D,1:nMatter) =  &
          FlattenMatrix(REAL(fMatrixOsc(1,iN_E,:,:)), &
                         nF, nSpecies/2 )
    
      uPR(iNodeE,iE,iX1,iX2,iX3,iPR_D,nAntiMatter:) =  &
          FlattenMatrix(REAL(fMatrixOsc(2,iN_E,:,:)), &
                         nF, nSpecies/2 )

    END DO

  END SUBROUTINE Flattenf

  SUBROUTINE OutputMatrixOsc( dt_out )

    REAL(DP), INTENT(IN) :: dt_out
    CHARACTER(LEN=100) :: FileName
    INTEGER :: iE

    TimeLocOutput = TimeLocOutput + dt_out
          
    DO iE = 1,nE_G
      
      WRITE(FileName,'(A6,I2.2,A4)') 'ImPart', iE, '.dat'
      OPEN(unit=666,file=TRIM(ADJUSTL(FileName)),&
          status="unknown",form='formatted',position="appEND")

      WRITE(666,"(7E19.9)") TimeLocOutput,Real(fMatrixOsc(1,iE,1,1)),Real(fMatrixOsc(1,iE,2,2)), &
      Real(fMatrixOsc(1,iE,1,2)),Aimag(fMatrixOsc(1,iE,1,2)), &
      Real(fMatrixOsc(1,iE,2,1)),Aimag(fMatrixOsc(1,iE,2,1))

      CLOSE(666)

    END DO

  END SUBROUTINE OutputMatrixOsc

END MODULE OscillationsModule

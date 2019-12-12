MODULE KernelsNu4Module

  USE KindModule, ONLY: &
    DP, Zero, One, Two, &
    Three, Five, Half, &
    Pi, TwoPi, FourPi
  USE UnitsModule, ONLY: &
    MeV, Second
  USE PhysicalConstantsModule, ONLY: &
    SpeedOfLightCGS
  USE ProgramHeaderModule, ONLY: &
    nNodesE, nDOFE, iE_E0, iE_B0
  USE MeshModule, ONLY: &
    MeshE, NodeCoordinate
  USE ReferenceElementModuleE, ONLY: &
    WeightsE

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputeNu4Kernels

  REAL(DP), ALLOCATABLE, PUBLIC :: KernelNu4Pair(:,:,:)
  REAL(DP), ALLOCATABLE, PUBLIC :: KernelNu4Scat(:,:,:)

CONTAINS

  SUBROUTINE ComputeNu4Kernels

    REAL(DP), ALLOCATABLE :: E_N(:)
    INTEGER  :: nE_G
    
    INTEGER  :: iN_Ek, iN_E1, iN_E3, iN_E
    INTEGER  :: nE, iE, iNodeE
    INTEGER  :: iE1, iE3, iNodeE1, iNodeE3
    REAL(DP) :: qk, q1, q2, q3, E2
    REAL(DP) :: V1, V3 ! phase-space volumes (cm^-3)
    REAL(DP) :: dE
    REAL(DP) :: hbarc_mevcm = 1.97326966d-11 
    
    nE = iE_E0 - iE_B0 + 1
    nE_G = nNodesE * nE

    ALLOCATE( KernelNu4Pair(nE_G,nE_G,nE_G) )
    ALLOCATE( KernelNu4Scat(nE_G,nE_G,nE_G) )

    ALLOCATE( E_N    (nE_G) )
    
    DO iN_E = 1, nE_G

      iE     = MOD( (iN_E-1) / nDOFE, nE    ) + iE_B0
      iNodeE = MOD( (iN_E-1)        , nDOFE ) + 1

      E_N(iN_E) = NodeCoordinate( MeshE, iE, iNodeE )

      dE = MeshE % Width(iE)
      
    ENDDO

    KernelNu4Pair(:,:,:) = Zero
    KernelNu4Scat(:,:,:) = Zero

    DO iN_Ek = 1,nE_G
       
       qk = E_N(iN_Ek) / MeV
       
       DO iN_E1 = 1,nE_G
         
         q1 = E_N(iN_E1) / MeV
         
         iE1     = MOD( (iN_E1-1) / nDOFE, nE    ) + iE_B0
         iNodeE1 = MOD( (iN_E1-1)        , nDOFE ) + 1 
         dE = MeshE % Width(iE1) 
         
         V1 = FourPi * dE * WeightsE(iNodeE1) / ( TwoPi )**3 / MeV / hbarc_mevcm
        
         DO iN_E3 = 1,nE_G
           
           q3 = E_N(iN_E3) / MeV
            
           iE3     = MOD( (iN_E3-1) / nDOFE, nE    ) + iE_B0
           iNodeE3 = MOD( (iN_E3-1)        , nDOFE ) + 1
           dE = MeshE % Width(iE3)

           V3 = Fourpi * dE * WeightsE(iNodeE3) / ( TwoPi )**3 / MeV / hbarc_mevcm 

           E2 = E_N(iN_E1) + E_N(iN_E3) - E_N(iN_Ek)

           ! result is 1/cm
           IF ( E2 > E_N(1) .AND. E2 <= E_N(nE_G) ) THEN
             q2 = E2 / MeV
             
             KernelNu4Pair(iN_E3,iN_E1,iN_Ek) = Nu4Pair_Kernel_Single( qk, q1, q2, q3) * V1 * V3
             
             KernelNu4Scat(iN_E3,iN_E1,iN_Ek) = Nu4Scat_Kernel_Single( qk, q1, q2, q3) * V1 * V3
             
             !NOW FIX UNITS SINCE SHERWOOD SOLVES 1/c*df/dt
             KernelNu4Pair(iN_E3,iN_E1,iN_Ek) = KernelNu4Pair(iN_E3,iN_E1,iN_Ek) * SpeedOfLightCGS * One / Second 
             KernelNu4Scat(iN_E3,iN_E1,iN_Ek) = KernelNu4Scat(iN_E3,iN_E1,iN_Ek) * SpeedOfLightCGS * One / Second
        
             !WRITE(*,'(3I4,2ES13.5E3)') iN_E3,iN_E1,iN_Ek, KernelNu4Pair(iN_E3,iN_E1,iN_Ek), KernelNu4Scat(iN_E3,iN_E1,iN_Ek)
            IF( KernelNu4Pair(iN_E3,iN_E1,iN_Ek) < 0 ) THEN
                WRITE(*,*) E2, iN_Ek, iN_E1, iN_E3, KernelNu4Pair(iN_E3,iN_E1,iN_Ek)
                STOP "KernelNu4Pair value less than 0"
             END IF

             IF( KernelNu4Scat(iN_E3,iN_E1,iN_Ek) < 0 ) THEN
                WRITE(*,*) E2, iN_Ek, iN_E1, iN_E3, KernelNu4Scat(iN_E3,iN_E1,iN_Ek)
                STOP "KernelNu4Scat value less than 0"
             END IF
    
           ELSE
 
             KernelNu4Scat(iN_E3,iN_E1,iN_Ek) = Zero
             KernelNu4Pair(iN_E3,iN_E1,iN_Ek) = Zero
           
           ENDIF
        
         END DO
      END DO
    END DO

  DEALLOCATE( E_N )

  END SUBROUTINE ComputeNu4Kernels

  !!=============================!!
  !! nu+nu <--> nu+nu scattering !!
  !!=============================!!
  ! second line in Blaschke & Cirigliano 2016 Equation 96
  ! technically applicable only in isotropic limit
  FUNCTION Nu4Scat_Kernel_Single(qk, q1, q2, q3) result(Kernel)

    REAL(DP), INTENT(IN) :: qk, q1, q2, q3
    REAL(DP) :: Kernel
    REAL(DP) :: Gfermi = 1.16637d-11 !MeV^-2
    REAL(DP) :: hbarc_mevcm = 1.97326966d-11

    IF( ABS( q1-q2+q3-qk ) > ( q1+q2+q3+qk )*1.d-10) THEN
       WRITE(*,*) qk, q1, q2, q3
       STOP "ERROR: passed q2 Does not conserve energy"
    ENDIF

    Kernel = Gfermi**2 * Pi * hbarc_mevcm * &
      (                      nu4_D3( q1,q2,q3,qk )   &
         + q1 * q2 * q3 * qk * nu4_D1( q1,q2,q3,qk )   &
         + q2 * qk           * nu4_D2( q2,qk,q1,q3 )   &
         + q1 * q3           * nu4_D2( q1,q3,q2,qk ) ) &
       / qk**2
    
  END FUNCTION Nu4Scat_Kernel_Single

  !!=====================================!!
  !! nu+nubar <--> nu+nubar annihilation !!
  !!=====================================!!
  ! fourth line in Blaschke & Cirigliano 2016 Equation 96
  ! technically applicable only in isotropic limit
  FUNCTION Nu4Pair_Kernel_Single(qk, q1, q2, q3) result(Kernel)
    
    REAL(DP), INTENT(IN) :: qk, q1, q2, q3
    REAL(DP) :: Kernel
    REAL(DP) :: Gfermi = 1.16637d-11 !MeV^-2 
    REAL(DP) :: hbarc_mevcm = 1.97326966d-11

    IF( ABS( q1-q2+q3-qk ) > ( q1+q2+q3+qk )*1.d-10) THEN
       WRITE(*,*) qk, q1, q2, q3
       STOP "ERROR: passed q2 Does not conserve energy"
    ENDIF
    
    Kernel = Gfermi**2 * Pi * hbarc_mevcm * &
        (                      nu4_D3( q1,q2,q3,qk )   &
         + q1 * q2 * q3 * qk * nu4_D1( q1,q2,q3,qk )   &
         - q1 * qk           * nu4_D2( q1,qk,q2,q3 )   &
         - q2 * q3           * nu4_D2( q2,q3,q1,qk ) ) &
        / qk**2

  END FUNCTION Nu4Pair_Kernel_Single

  ! equation D4c in Blaschke+Cirigliano 2016
  FUNCTION nu4_D4c(q1,q2,q3,q4) RESULT(D4c) ! MeV^5
    
    REAL(DP), INTENT(IN) :: q1, q2, q3, q4
    REAL(DP) :: D4c

    D4c = One/60.0 * ( q1**5 - Five * q1**3*q2**2 + Five * q1**2*q2**3 - &
                       q2**5 - Five * q1**3*q3**2 + Five * q2**3*q3**2 + &
                               Five * q1**2*q3**3 + Five * q2**2*q3**3 - &
                       q3**5 - Five * q1**3*q4**2 + Five * q2**3*q4**2 + &
                               Five * q3**3*q4**2 + Five * q1**2*q4**3 + &
                               Five * q2**2*q4**3 + Five * q3**2*q4**3 - &
                       q4**5 )
  
  END FUNCTION nu4_D4c
  
  !!====!!
  !! D1 !!
  !!====!!
  FUNCTION nu4_D1(qin1,qin2,qin3,qin4) RESULT(D1)! MeV
  
    REAL(DP), INTENT(IN) :: qin1, qin2, qin3, qin4
    REAL(DP) :: q1, q2, q3, q4, D1

    ! all RESULTs are symmetric on 1<-->2.AND.3<-->4
    ! use the case where q1>q2.AND.q3>q4
    
    q1 = MAX(qin1,qin2)
    q2 = MIN(qin1,qin2)
    q3 = MAX(qin3,qin4)
    q4 = MIN(qin3,qin4)

    ! Case 1 (Eqns. D4)
    IF( q1+q2 >= q3+q4 .AND. q1+q4 >= q2+q3 ) THEN
      
      IF( q1 <= q2+q3+q4 ) THEN
        D1 = Half  * ( q2+q3+q4-q1 )
      ELSE
        D1 = Zero
      END IF

    ! Case 2 (Eqns. D5)
    ELSE IF( q1+q2 >= q3+q4 .AND. q1+q4 <= q2+q3 ) THEN
      
      D1 = q4

    ! Case 3 (Eqns. D6)
    ELSE IF( q1+q2 <= q3+q4 .AND. q1+q4 <= q2+q3 ) THEN
      
      IF( q3 <= q1+q2+q4 ) THEN
        D1 = Half * ( q1+q2+q4-q3 )
      ELSE
        D1 = Zero
      END IF

    ! Case 4 (Eqns. D7)
    ELSE !(q1+q2<q3+q4 .AND. q1+q4>q2+q3)
      
      D1 = q2
    
    END IF
  
  END FUNCTION nu4_D1
  
  !!====!!
  !! D2 !!
  !!====!!
  FUNCTION nu4_D2(qin1,qin2,qin3,qin4) RESULT(D2)! MeV^3
  
    REAL(DP), INTENT(IN) :: qin1, qin2, qin3, qin4
    REAL(DP) :: q1, q2, q3, q4, D2

    ! all RESULTs are symmetric on 1<-->2.AND.3<-->4
    ! use the case where q1>q2.AND.q3>q4
    q1 = MAX(qin1,qin2)
    q2 = MIN(qin1,qin2)
    q3 = MAX(qin3,qin4)
    q4 = MIN(qin3,qin4)

    ! Case 1 (Eqns. D4)
    IF(q1+q2 >= q3+q4 .AND. q1+q4>=q2+q3) THEN
      
      IF(q1<=q2+q3+q4) THEN
        D2 = One/12.0 * ( (q1-q2)**3 + Two * ( q3**3+q4**3 ) &
                          - Three * (q1-q2) * ( q3**2+q4**2) )
      ELSE
        D2 = Zero
      END IF

    ! Case 2 (Eqns. D5)
    ELSE IF ( q1+q2 >= q3+q4 .AND. q1+q4 <= q2+q3 ) THEN
      
      D2 = q4**3 / Three

    ! Case 3 (Eqns. D6)
    ELSE IF( q1+q2 <= q3+q4 .AND. q1+q4 <= q2+q3 ) THEN
    
      IF( q3<=q1+q2+q4 ) THEN
        D2 = One/12.0 * ( - (q1+q2)**3 - Two * q3**3 + &
            Two * q4**3 + Three * (q1+q2) * (q3**2+q4**2) )
      ELSE
        D2 = Zero
      END IF

    ! Case 4 (Eqns. D7)
    ELSE !(q1+q2<q3+q4 .AND. q1+q4>q2+q3)
    
      D2 = q2/6.0 * ( Three * q3**2 + Three * q4**2 - &
                      Three * q1**2 -         q2**2 )
    
    END IF

  END FUNCTION nu4_D2
  
  !!====!!
  !! D3 !!
  !!====!!
  FUNCTION nu4_D3(qin1,qin2,qin3,qin4) RESULT(D3)! MeV^5
  
    REAL(DP), INTENT(IN) :: qin1, qin2, qin3, qin4
    REAL(DP) :: q1, q2, q3, q4, D3

    ! all RESULTs are symmetric on 1<-->2.AND.3<-->4
    ! use the case where q1>q2.AND.q3>q4
    q1 = MAX(qin1,qin2)
    q2 = MIN(qin1,qin2)
    q3 = MAX(qin3,qin4)
    q4 = MIN(qin3,qin4)

    ! Case 1 (Eqns. D4)
    IF( q1+q2 >= q3+q4 .AND. q1+q4 >= q2+q3 ) THEN
      
      IF( q1 <= q2+q3+q4 ) THEN
        D3 = nu4_D4c(q1,q2,q3,q4)
      ELSE
        D3 = Zero
      END IF

    ! Case 2 (Eqns. D5)
    ELSE IF( q1+q2 >= q3+q4 .AND. q1+q4 <= q2+q3 ) THEN
      
      D3 = q4**3 / 30.0 * ( Five * q1**2 + Five * q2**2 + &
                            Five * q3**2 -        q4**2 )
    
    ! Case 3 (Eqns. D6)
    ELSE IF( q1+q2 <= q3+q4 .AND. q1+q4 <= q2+q3 ) THEN
        
      IF( q3 <= q1+q2+q4 ) THEN
        D3 = nu4_D4c(q3,q4,q1,q2)
      ELSE
        D3 = Zero
      END IF

    ! Case 4 (Eqns. D7)
    ELSE !(q1+q2<q3+q4 .AND. q1+q4>q2+q3)
      
      D3 = q2**3/30.0 * ( Five *q1**2 + Five * q3**2 + &
                          Five *q4**2 -        q2**2 )
    END IF
  
  END FUNCTION nu4_D3



END MODULE KernelsNu4Module

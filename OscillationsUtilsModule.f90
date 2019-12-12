MODULE OscillationsUtilsModule
 
  USE KindModule, ONLY: &
    DP, Zero, Half, One, &
    Two, Pi, TwoPi
  USE InitializationModule, ONLY: &
    Enu, nY, AV, CV, &
    CofactorMatrix

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: B, W
  PUBLIC :: JInverse
  PUBLIC :: Eigenvalues, EigenvectorMatrix

  COMPLEX(DP), PARAMETER, PUBLIC :: Im = (Zero,One)
  ! signs of Eigenvalues of vacuum Hamiltonian (switch them IF dm21<0)
  REAL(DP), PARAMETER :: a1 = - One
  REAL(DP), PARAMETER :: a2 = + One

CONTAINS

  FUNCTION Eigenvalues(M) RESULT(Lambda)
     
    COMPLEX(DP), INTENT(IN)  :: M(2,2)
    REAL(DP) :: Lambda(2)
    
    REAL(DP) :: Trace, Discr

    Trace = REAL(M(1,1) + M(2,2))
    Discr = 4.0d0*ABS(M(1,2))**2 + ABS(M(1,1) - M(2,2))**2
    

    Lambda(1) = Half * ( Trace + a1*SQRT(Discr) )
    Lambda(2) = Half * ( Trace + a2*SQRT(Discr) )
  
  END FUNCTION Eigenvalues

  FUNCTION EigenvectorMatrix(M,lambda,i_loc) RESULT(U)
    
    COMPLEX(DP), INTENT(IN) :: M(2,2) !i.e. Hf
    REAL(DP), INTENT(IN) :: lambda(2)

    INTEGER, INTENT(IN) :: i_loc

    COMPLEX(DP) :: U(2,2)
    
    !local variables
    COMPLEX(DP) :: C0(2,2,2) !cofactor matrices
    COMPLEX(DP) :: A0(2,2)
    REAL(DP) :: deltak,delta
    REAL(DP) :: d !discrminant
    REAL(DP) :: r2(2)
    INTEGER :: j,n
    
    C0(:,:,:) = CofactorMatrix(M,lambda) 
    
    DO j=1,2

        IF(REAL(C0(j,2,1)*CV(i_loc,j,2,1)) .lt. 0.0d0) THEN
             A0(j,1)=-AV(i_loc,j,1)
        else 
            A0(j,1)=AV(i_loc,j,1)
        END IF
        A0(j,2)=AV(i_loc,j,2)
    
    END DO   

    d = 4.0d0*ABS(M(1,2))**2 + ABS(M(1,1) - M(2,2))**2
    deltak = a1*SQRT(d)
    
    DO j=1,2

      IF(j.eq.1) delta = -deltak ! first column
      IF(j.eq.2) delta =  deltak ! second column
        
      r2(1) = REAL(C0(j,1,1))*delta
      r2(2) = REAL(C0(j,2,2))*delta
        
      IF(r2(1).ge.r2(2)) THEN

        U(1,j) = A0(j,1) * C0(j,1,1) / SQRT(r2(1))
        U(2,j) = A0(j,1) * C0(j,1,2) / SQRT(r2(1))
        
       END IF
       
       IF(r2(2).ge.r2(1)) THEN
       
         U(1,j) = A0(j,2) * C0(j,2,1) / SQRT(r2(2))
         U(2,j) = A0(j,2) * C0(j,2,2) / SQRT(r2(2))
        
       END IF
    
    END DO
  

  END FUNCTION EigenvectorMatrix

  FUNCTION W(Y) RESULT(WY)

    REAL(DP), INTENT(IN) :: Y(nY)
    
    COMPLEX(DP) :: WY(2,2)

    WY(:,:) = Zero
    WY(1,1) = EXP( -Im*TwoPi*Y(5) )
    WY(2,2) = EXP( -Im*TwoPi*Y(6) )

  END FUNCTION W

  FUNCTION B(Y) RESULT(BY)

    REAL(DP), INTENT(IN) :: Y(nY)

    COMPLEX(DP) :: BY(2,2)

    !local
    REAL(DP) :: cy1, sy1 
    REAL(DP) :: cy2, sy2
    REAL(DP) :: cy3, sy3

    cy1 = COS(Y(1)) 
    sy1 = SIN(Y(1))
    cy2 = COS(Y(2))
    sy2 = SIN(Y(2))
    cy3 = COS(Y(3))
    sy3 = SIN(Y(3))

    BY(1,2) = cy1 + im*sy1*cy2
    BY(1,1) = sy1*sy2 * (cy3 + im*sy3)
    BY(2,1) = -Y(4)*CONJG(BY(1,2))
    BY(2,2) =  Y(4)*CONJG(BY(1,1))
    
  END FUNCTION B

  FUNCTION JInverse(Y) RESULT(JI)
    
    REAL(DP), INTENT(IN) :: Y(nY)

    REAL(DP) :: JI(3,4)

    !local
    REAL(DP) :: cy1, sy1             
    REAL(DP) :: cy2, sy2
    REAL(DP) :: cy3, sy3

    cy1 = COS(Y(1))
    sy1 = SIN(Y(1))
    cy2 = COS(Y(2))
    sy2 = SIN(Y(2))
    cy3 = COS(Y(3))
    sy3 = SIN(Y(3))
    
    ! Notice that JI is a triangular superior matrix !

    JI(1,1) = -sy1
    JI(2,2) = -sy2/sy1
    JI(3,3) = -sy3/sy1
    
    JI(1,2) = cy1*cy2
    JI(2,3) = cy2*cy3/sy1
    JI(3,4) = cy3/sy1

    JI(1,3) = cy1*sy2*cy3
    JI(2,4) = cy2*sy3/sy1

    JI(1,4) = cy1*sy2*sy3

  END FUNCTION JInverse

  FUNCTION MatrixMul( Matrix1, Matrix2, dimM) RESULT(Matrix3)

    INTEGER,  INTENT(IN) :: dimM
    COMPLEX(DP), INTENT(IN) :: Matrix1(dimM,dimM), Matrix2(dimM,dimM)
    COMPLEX(DP)             :: Matrix3(dimM,dimM)

    INTEGER :: i,j,k

    DO j = 1,dimM
      DO k = 1,dimM
        DO i = 1,dimM

          Matrix3(i,j) = Matrix1(i,j) + Matrix1(i,k) * Matrix2(k,j)

        END DO
      END DO
    END DO

  END FUNCTION MatrixMul

END MODULE OscillationsUtilsModule


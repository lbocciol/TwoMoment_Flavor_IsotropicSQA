MODULE ReadProfileModule

  USE KindModule, ONLY: &
    DP
  USE wlIOModuleHDF, ONLY: &
    OpenFileHDF, &
    CloseFileHDF
  
  USE HDF5

  IMPLICIT NONE
  PRIVATE
  
  PUBLIC :: LinearInterpolation1D
  PUBLIC :: QuadraticInterpolation1D
  PUBLIC :: Read1DChimeraProfile
  PUBLIC :: ReadGR1DProfile
  
CONTAINS
    
  SUBROUTINE LinearInterpolation1D(X_Old, Y_Old, n0, &
                                   X_New, Y_new) 

    REAL(DP), INTENT(OUT)   :: Y_New
    REAL(DP), INTENT(IN)    :: X_New
    REAL(DP), INTENT(IN)    :: X_Old(n0)
    REAL(DP), INTENT(IN)    :: Y_Old(n0)
    INTEGER, INTENT(IN)     :: n0
    
    INTEGER :: i
     
    IF (X_New .le. X_Old(1)) THEN

      !Extrapolate down
      Y_New = Y_Old(1) - (X_Old(1) - X_New) * (Y_Old(2) - Y_Old(1)) &
                                            / (X_Old(2) - X_Old(1)) 

      RETURN

    ELSE IF (X_New .ge. X_Old(n0)) THEN
      
      !Exptrapolate up
      Y_New = Y_Old(n0) + (X_New - X_Old(n0)) * (Y_Old(n0) - Y_Old(n0-1)) &
                                              / (X_OLD(n0) - X_Old(n0-1))
      RETURN
    
    ELSE
      DO i = 1,n0
        
        IF(X_Old(i) .ge. X_New) THEN
            Y_New = Y_Old(i-1) + (X_New    - X_Old(i-1)) &
                               * (Y_Old(i) - Y_Old(i-1)) &
                               / (X_Old(i) - X_Old(i-1))
            RETURN
        
        END IF

      END DO
    END IF                   

  END SUBROUTINE LinearInterpolation1D

  SUBROUTINE QuadraticInterpolation1D(X_Old, Y_Old, n0, &
                                   X_New, Y_new)

    REAL(DP), INTENT(OUT)   :: Y_New
    REAL(DP), INTENT(IN)    :: X_New
    REAL(DP), INTENT(IN)    :: X_Old(n0)
    REAL(DP), INTENT(IN)    :: Y_Old(n0)
    INTEGER, INTENT(IN)     :: n0

    INTEGER  :: i
    REAL(DP) :: L2_0, L2_1, L2_2

    IF (X_New <= X_Old(1)) THEN

      !Extrapolate down
      Y_New = Y_Old(1) - (X_Old(1) - X_New) * (Y_Old(2) - Y_Old(1)) &
                                            / (X_Old(2) - X_Old(1))

      RETURN

    ELSE IF ( X_New <= X_Old(2) ) THEN

      Y_New = Y_Old(1) + (X_New    - X_Old(1)) &
                       * (Y_Old(2) - Y_Old(1)) &
                       / (X_Old(2) - X_Old(1))
    
    ELSE IF (X_New >= X_Old(n0)) THEN

      !Exptrapolate up
      Y_New = Y_Old(n0) + (X_New - X_Old(n0)) * (Y_Old(n0) - Y_Old(n0-1)) &
                                              / (X_OLD(n0) - X_Old(n0-1))
      RETURN

    ELSE IF ( X_New >= X_Old(n0-1) ) THEN

      Y_New = Y_Old(n0-1) + (X_New     - X_Old(n0-1)) &
                          * (Y_Old(n0) - Y_Old(n0-1)) &
                          / (X_Old(n0) - X_Old(n0-1))

    ELSE
  
      DO i = 1,n0
    
        IF( X_Old(i) >= X_New ) THEN
        
          L2_0 = ( (X_New      - X_Old(i)) * (X_New      - X_Old(i+1)) ) / &
                 ( (X_Old(i-1) - X_Old(i)) * (X_Old(i-1) - X_Old(i+1)) )
      
          L2_1 = ( (X_New    - X_Old(i-1)) * (X_New    - X_Old(i+1)) ) / &
                 ( (X_Old(i) - X_Old(i-1)) * (X_Old(i) - X_Old(i+1)) )
          
          L2_2 = ( (X_New      - X_Old(i-1)) * (X_New      - X_Old(i)) ) / &
                 ( (X_Old(i+1) - X_Old(i-1)) * (X_Old(i+1) - X_Old(i)) )  
          
          Y_New = Y_Old(i-1)*L2_0 + Y_Old(i)*L2_1 + Y_Old(i+1)*L2_2
          
          RETURN

        END IF

      END DO

    END IF

  END SUBROUTINE QuadraticInterpolation1D

  SUBROUTINE ReadGR1DProfile(R_GR1D, D_GR1D, T_GR1D, Ye_GR1D, n0, &
                            TimeSlice)
    
    REAL(DP), INTENT(INOUT) :: R_GR1D(n0), D_GR1D(n0)
    REAL(DP), INTENT(INOUT) :: T_GR1D(n0), Ye_GR1D(n0)
    REAL(DP), INTENT(INOUT) :: TimeSlice
    INTEGER, INTENT(IN)     :: n0
    
    INTEGER :: istate
    REAL(DP) :: Tbounce

    OPEN(UNIT = 101, FILE = "tbounce.dat", STATUS = "old", &
      FORM = 'formatted', IOSTAT = istate)
      READ(101,'(E27.18)') Tbounce
    CLOSE(101)

    TimeSlice = TimeSlice + Tbounce
    
    WRITE(*,*)

    CALL ReadGR1DProfileXg(R_GR1D, D_GR1D, n0, 4,      &
                           TimeSlice, "rho.xg")
    
    CALL ReadGR1DProfileXg(R_GR1D, T_GR1D, 600, 4,      &
                           TimeSlice, "temperature.xg")

    CALL ReadGR1DProfileXg(R_GR1D, Ye_GR1D, 600, 4,     &
                           TimeSlice, "ye.xg")
    
  END SUBROUTINE ReadGR1DProfile

  SUBROUTINE Read1DChimeraProfile(R_Ch,D_Ch,T_Ch,Ye_Ch, n0, TimeSlice, iMax)

    REAL(DP), INTENT(INOUT) :: R_Ch(n0), D_Ch(n0) 
    REAL(DP), INTENT(INOUT) :: T_Ch(n0), Ye_Ch(n0)
    REAL(DP), INTENT(IN)    :: TimeSlice              
    INTEGER,  INTENT(IN)    :: n0, iMax

    INTEGER        :: HDFERR
    INTEGER        :: iFile
    CHARACTER(5)   :: FileNumberString 
    CHARACTER(256) :: FileName
    CHARACTER(256) :: SetName

    REAL(DP)       :: Time, Tbounce


    CALL H5OPEN_F( HDFERR )

    DO iFile = 1,iMax
      
      WRITE(FileNumberString, FMT = '(I5.5)') iFile

      FileName = "CHIMERA/SFHo/E15-1DSFHo-Frames-00/chimera_" &
            // FileNumberString // "_grid_1_01.h5"

      SetName = "mesh/time"
      CALL ReadHDFScalar(FileName, SetName, Time)
    
      SetName = "mesh/t_bounce"
      CALL ReadHDFScalar(FileName, SetName, Tbounce)

      IF ( Time - Tbounce >= TimeSlice .AND. Tbounce > 0 ) &
        EXIT

    END DO

    WRITE(*,*)
    WRITE(*,'(A4,A)') '','Reading Chimera Profiles...'

    SetName = "mesh/x_cf"
    CALL ReadHDFVector( FileName, SetName, R_Ch, n0 )
    
    SetName = "fluid/rho_c"
    CALL ReadHDFVector( FileName, SetName, D_Ch, n0 )

    SetName = "fluid/t_c"
    CALL ReadHDFVector( FileName, SetName, T_Ch, n0 )

    SetName = "fluid/ye_c"
    CALL ReadHDFVector( FileName, SetName, Ye_Ch, n0 )

    CALL H5CLOSE_F( HDFERR )

  END SUBROUTINE Read1DChimeraProfile

  SUBROUTINE ReadHDFScalar(FileName, SetName, Variable)
    
    CHARACTER(256), INTENT(IN)  :: FileName
    CHARACTER(256), INTENT(IN)  :: SetName
    REAL(DP),       INTENT(OUT) :: Variable
    
    INTEGER(HID_T)              :: FILE_ID
    INTEGER(HID_T)              :: DSET_ID
    INTEGER                     :: HDFERR
    INTEGER(HSIZE_T)            :: DIMS(1)

    DIMS(1) = 1 

    ! --- Open File --- !
    CALL H5FOPEN_F( TRIM(ADJUSTL( FileName )), &
           H5F_ACC_RDONLY_F, FILE_ID, HDFERR)
    
    CALL H5DOPEN_F( FILE_ID, TRIM(ADJUSTL( SetName )), &
              DSET_ID, HDFERR)
    
    ! --- Read Data --- !
    CALL H5DREAD_F( DSET_ID, H5T_NATIVE_DOUBLE, & 
              Variable, DIMS, HDFERR )
    
    ! --- Close File --- !
    CALL H5DCLOSE_F( DSET_ID, HDFERR )
      
    CALL H5FCLOSE_F( FILE_ID, HDFERR )

  END SUBROUTINE ReadHDFScalar 


  SUBROUTINE ReadHDFVector(FileName, SetName, Variable, n0)

    CHARACTER(256), INTENT(IN)  :: FileName
    CHARACTER(256), INTENT(IN)  :: SetName
    INTEGER,        INTENT(IN)  :: n0
    REAL(DP),       INTENT(OUT) :: Variable(n0)

    INTEGER(HID_T)              :: FILE_ID
    INTEGER(HID_T)              :: DSET_ID
    INTEGER                     :: HDFERR
    INTEGER(HSIZE_T)            :: DIMS(1)

    DIMS(1) = n0

    ! --- Open File --- !
    CALL H5FOPEN_F( TRIM(ADJUSTL( FileName )), &
           H5F_ACC_RDONLY_F, FILE_ID, HDFERR)

    CALL H5DOPEN_F( FILE_ID, TRIM(ADJUSTL( SetName )), &
              DSET_ID, HDFERR)

    ! --- Read Data --- !
    CALL H5DREAD_F( DSET_ID, H5T_NATIVE_DOUBLE, &
              Variable, DIMS, HDFERR )
    
    ! --- Close File --- !
    CALL H5DCLOSE_F( DSET_ID, HDFERR )

    CALL H5FCLOSE_F( FILE_ID, HDFERR )

  END SUBROUTINE ReadHDFVector


  SUBROUTINE ReadGR1DProfileXg(R_GR1D, U_GR1D, n0, nGhost, &
                               TimeSlice, FileName)

    REAL(DP), INTENT(INOUT)       :: R_GR1D(n0)
    REAL(DP), INTENT(INOUT)       :: U_GR1D(n0)
    INTEGER, INTENT(IN)           :: n0, nGhost
    REAL(DP), INTENT(IN)          :: TimeSlice
    CHARACTER(len=*), INTENT(IN)  :: FileName

    CHARACTER(len=200) :: TimeString
    REAL(DP) :: Time 
    INTEGER  :: i, istate

    111 FORMAT(1P20E18.9)
    
    WRITE(*,'(A4,2A)') '','Reading from GR1D: ',TRIM(ADJUSTL(FileName))
    
    OPEN(UNIT=999, FILE = TRIM(ADJUSTL(FileName)), STATUS = "old", &
        FORM = 'formatted', IOSTAT = istate) 

    DO WHILE (istate .ge. 0 )
      
      READ(999,"(A8,E30.9)") TimeString, Time
      
      !Skip ghost
      DO i = 1,nGhost
        READ(999,*)
      END DO

      DO i = 1,n0 
        READ(999,111) R_GR1D(i),U_GR1D(i)
      END DO
      
      !Skip ghost and 2 blank spaces
      DO i = 1,nGhost + 2
        READ(999,*) 
      END DO
      
      IF (Time .ge. TimeSlice) &
        EXIT
    END DO
    CLOSE(999)

  END SUBROUTINE ReadGR1DProfileXg

  SUBROUTINE ReadGR1DProfileDat(TimeSlice, Var, FileName)            
                                                                      
    REAL(DP), INTENT(INOUT)       :: Var
    REAL(DP), INTENT(IN)          :: TimeSlice         
    CHARACTER(len=*), INTENT(IN)  :: FileName                           
                                                                        
    REAL(DP) :: Time                                                    
    INTEGER  :: istate                                               
                                                                        
    111 FORMAT(1P20E18.9)                                               
    WRITE(*,'(A8,A)') "Reading",TRIM(ADJUSTL(FileName))                                                                                     
    OPEN(UNIT=999, FILE = TRIM(ADJUSTL(FileName)), STATUS = "old", &    
        FORM = 'formatted', IOSTAT = istate)                            
                                                                        
    DO WHILE (istate .ge. 0 )                                           
                                                                      
      READ(999,111) Time, Var                             
                                                           
      IF (Time .ge. TimeSlice) &
        EXIT
    END DO
    CLOSE(999)

  END SUBROUTINE ReadGR1DProfileDat

END MODULE ReadProfileModule 

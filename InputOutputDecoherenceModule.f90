MODULE InputOutputDecoherenceModule

  USE KindModule, ONLY: &
    DP
  USE UnitsModule, ONLY: &
    UnitsDisplay
  USE ProgramHeaderModule, ONLY: &
    ProgramName, &
    iX_B0, iX_E0, &
    nDOF
  USE RadiationFieldsModule, ONLY: &
    uPR, iPR_D, nSpecies
  USE InitializationModule, ONLY: &
    nX_G, nE_G, nM, &
    fImaginaryPart

  USE HDF5

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: OutputImaginaryPart
  PUBLIC :: InitializeRestart
    
  CHARACTER(9),  PARAMETER :: &
    OutputDirectory = '../Output'
  CHARACTER(13), PARAMETER :: &
    ImPartSuffix = 'ImaginaryPart'
  CHARACTER(15), PARAMETER :: &
    RadiationSuffix = 'RadiationFields'
  INTEGER :: FileNumber = 0

CONTAINS

  SUBROUTINE OutputImaginaryPart( Time )
    
    REAL(DP), INTENT(IN) :: Time
    
    INTEGER(HSIZE_T) :: DATASIZE(1)
    INTEGER(HID_T)   :: FILE_ID
    INTEGER(HID_T)   :: DSPACE_ID
    INTEGER(HID_T)   :: DSET_ID
    INTEGER(HID_T)   :: DIMS1(1), DIMS3(3)

    INTEGER :: HDFERR
    INTEGER :: Rank

    CHARACTER(256) :: FileName
    CHARACTER(256) :: DatasetName
    CHARACTER(6)   :: FileNumberString

    WRITE( FileNumberString, FMT='(i6.6)') FileNumber
    
    FileName &
      = OutputDirectory // '/' // &
        TRIM( ProgramName ) // '_' // &
        ImPartSuffix // '_' // &
        FileNumberString // '.h5'
    
    ! Create File
    CALL H5OPEN_F( HDFERR )

    CALL H5FCREATE_F( TRIM( FileName ), H5F_ACC_TRUNC_F, FILE_ID, HDFERR )

    ! Create Time Dataset
    DatasetName = '/Time'
    Rank = 1
    DIMS1(1) = 1

    CALL H5SCREATE_F( H5S_SIMPLE_F, DSPACE_ID, HDFERR )

    CALL H5SSET_EXTENT_SIMPLE_F &
           ( DSPACE_ID, Rank, DIMS1, DIMS1, HDFERR )

    CALL H5DCREATE_F &
        ( FILE_ID, TRIM( DatasetName ), H5T_NATIVE_DOUBLE, &
          DSPACE_ID, DSET_ID, HDFERR )

    ASSOCIATE( U => UnitsDisplay )

    CALL H5DWRITE_F( DSET_ID, H5T_NATIVE_DOUBLE, &
          Time / U % TimeUnit, DIMS1, HDFERR )

    END ASSOCIATE

    CALL H5DCLOSE_F( DSET_ID,   HDFERR )

    CALL H5SCLOSE_F( DSPACE_ID, HDFERR )  

    ! Create ImaginaryPart Dataset

    DatasetName = '/fImaginaryPart'
    Rank = 3
    DIMS3(1) = nX_G
    DIMS3(2) = nE_G
    DIMS3(3) = nM
    
    CALL H5SCREATE_F( H5S_SIMPLE_F, DSPACE_ID, HDFERR )

    CALL H5SSET_EXTENT_SIMPLE_F &
           ( DSPACE_ID, Rank, DIMS3, DIMS3, HDFERR )

    CALL H5DCREATE_F &
        ( FILE_ID, TRIM( DatasetName ), H5T_NATIVE_DOUBLE, &
          DSPACE_ID, DSET_ID, HDFERR )

    CALL H5DWRITE_F( DSET_ID, H5T_NATIVE_DOUBLE, fImaginaryPart, DIMS3, HDFERR )

    CALL H5DCLOSE_F( DSET_ID,   HDFERR )

    CALL H5SCLOSE_F( DSPACE_ID, HDFERR )

    ! Close File

    CALL H5FCLOSE_F( FILE_ID, HDFERR )

    CALL H5CLOSE_F( HDFERR )
    
    FileNumber = FileNumber + 1 
   
  END SUBROUTINE 
  
  SUBROUTINE InitializeRestart( RestartNumber, Time, &
          nX, swX, nE, swE) 

    INTEGER, INTENT(IN) :: nX(3), swX(3)
    INTEGER, INTENT(IN) :: nE,    swE
    INTEGER, INTENT(IN) :: RestartNumber

    REAL(DP), INTENT(INOUT) :: Time
    
    REAL(DP), ALLOCATABLE :: Dummy(:,:,:,:,:)

    INTEGER(HSIZE_T) :: DATASIZE(1)
    INTEGER(HID_T)   :: FILE_ID
    INTEGER(HID_T)   :: DSPACE_ID
    INTEGER(HID_T)   :: DSET_ID
    INTEGER(HID_T)   :: DIMS1(1), DIMS3(3), DIMS5(5)

    INTEGER :: HDFERR
    INTEGER :: iS

    CHARACTER(256) :: FileName
    CHARACTER(256) :: DatasetName 
    CHARACTER(6)   :: FileNumberString
    CHARACTER(10)   :: SpeciesString
   
    FileNumber = RestartNumber
    WRITE( FileNumberString, FMT='(i6.6)') FileNumber

    CALL H5OPEN_F( HDFERR )
    
    !Open Radiation File
    FileName &
      = OutputDirectory // '/' // &
        TRIM( ProgramName ) // '_' // &
        RadiationSuffix // '_' // &
        FileNumberString // '.h5'

    CALL H5FOPEN_F( TRIM( FileName ), &
        H5F_ACC_RDONLY_F, FILE_ID, HDFERR )

    IF( HDFERR < 0 ) THEN

      WRITE(*,'(A12,A1,I6.6,A1,A9)') &
          'Restart File','',RestartNumber,'','Not Found'
      STOP

    END IF

    DIMS5(1) = nDOF
    DIMS5(2) = 2*swE + nE
    DIMS5(3) = 2*swX(1) + nX(1)
    DIMS5(4) = 2*swX(2) + nX(2)
    DIMS5(5) = 2*swX(3) + nX(3)
    
    DO iS = 1,nSpecies
      
      WRITE( SpeciesString, FMT='(A8,I2.2)') 'Species_',iS

      DataSetName = '/Radiation Fields/' // SpeciesString // '/Primitive/Lagrangian Number Density'

      CALL H5DOPEN_F( FILE_ID, DatasetName, &
          DSET_ID, HDFERR )

      CALL H5DREAD_F( DSET_ID, H5T_NATIVE_DOUBLE, &
          uPR(1:nDOF,1:nE,1:nX(1),1:nX(2),1:nX(3),iPR_D,iS), &
          DIMS5, HDFERR )
      
      CALL H5DCLOSE_F( DSET_ID, HDFERR )
    
    END DO
    
    CALL H5FCLOSE_F( FILE_ID, HDFERR )

    CALL H5CLOSE_F( HDFERR )

    ! Read Imaginary Part
    CALL H5OPEN_F( HDFERR )

    FileName &
      = OutputDirectory // '/' // &
        TRIM( ProgramName ) // '_' // &
        ImPartSuffix // '_' // &
        FileNumberString // '.h5'

    CALL H5FOPEN_F( TRIM( FileName ), &
        H5F_ACC_RDONLY_F, FILE_ID, HDFERR )

    DataSetName = '/fImaginaryPart'
    DIMS3(1) = nX_G
    DIMS3(2) = nE_G
    DIMS3(3) = nM

    CALL H5DOPEN_F( FILE_ID, DatasetName, &
        DSET_ID, HDFERR )

    CALL H5DREAD_F( DSET_ID, H5T_NATIVE_DOUBLE, &
        fImaginaryPart, DIMS3, HDFERR )

    CALL H5DCLOSE_F( DSET_ID, HDFERR )

    ! Read Time
    DataSetName = '/Time'
    DIMS1(1) = 1

    CALL H5DOPEN_F( FILE_ID, DatasetName, &
        DSET_ID, HDFERR )

    CALL H5DREAD_F( DSET_ID, H5T_NATIVE_DOUBLE, &
        Time, DIMS1, HDFERR )
    
    ASSOCIATE( U => UnitsDisplay )

      Time = Time * U % TimeUnit 

    END ASSOCIATE
    
    CALL H5DCLOSE_F( DSET_ID, HDFERR )
    
    CALL H5FCLOSE_F( FILE_ID, HDFERR )

    CALL H5CLOSE_F( HDFERR )
    
    FileNumber = FileNumber + 1

  END SUBROUTINE InitializeRestart

END MODULE InputOutputDecoherenceModule

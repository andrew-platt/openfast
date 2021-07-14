MODULE LidarSim_Subs

    USE LidarSim_Types
    USE NWTC_Library
    USE InflowWind_Subs
    USE InflowWind_Types

    IMPLICIT NONE
    PRIVATE
    
    PUBLIC  ::  LidarSim_ParsePrimaryFileInfo
    PUBLIC  ::  LidarSim_InitMeasuringPoints_Cartesian
    PUBLIC  ::  LidarSim_TransformLidarToInertial
    PUBLIC  ::  LidarSim_InitMeasuringPoints_Spherical
    PUBLIC  ::  LidarSim_InitializeWeightingGauss
    PUBLIC  ::  LidarSim_InitializeWeightingManual
    PUBLIC  ::  LidarSim_CalculateVlos
    PUBLIC  ::  LidarSim_InitializeOutputs
    PUBLIC  ::  LidarSim_SetOutputs
    PUBLIC  ::  LidarSim_CalculateIMU

CONTAINS
    

!#########################################################################################################################################################################

SUBROUTINE LidarSim_ParsePrimaryFileInfo( PriPath, InputFile, RootName, FileInfo, InputFileData, UnitEcho, ErrStat, ErrMsg)
   !  This subroutine originally opened the file as two units simultaneously. Unfortunately, the Fortran
   !  standard does not currenlty allow for opening the same file twice simultaneously. This behaviour
   !  is an extension to the standard that only the Intel Fortran compiler allows.  gfortran and other
   !  fortran compilers do not allow for this behaviour.
   !
   !  Therefore, to support an arbitrary number of comment lines, the file will be read into memory, then
   !  parsed line by line.

   implicit                                none

      ! Passed variables
   character(*),                    intent(in   )  :: PriPath           !< primary path
   CHARACTER(*),                    intent(in   )  :: InputFile         !< Name of the file containing the primary input data
   CHARACTER(*),                    intent(in   )  :: RootName          !< The rootname of the echo file, possibly opened in this routine
   type(LidarSim_InputFile),        intent(inout)  :: InputFileData     !< All the data in the AD15 primary input file
   type(FileInfoType),              intent(in   )  :: FileInfo          !< The derived type for holding the file information.
   integer(IntKi),                  intent(  out)  :: UnitEcho          !< The local unit number for this module's echo file
   integer(IntKi),                  intent(  out)  :: ErrStat           !< Error status
   CHARACTER(ErrMsgLen),            intent(  out)  :: ErrMsg            !< Error message

   character(*),                    parameter      ::  RoutineName="LidarSim_ParsePrimaryFileInfo"

   ! Local variables
   INTEGER(IntKi)                                  ::  TmpErrStat
   CHARACTER(ErrMsgLen)                            ::  TmpErrMsg           !< temporary error message
   INTEGER(IntKi)                                  ::  ErrStatIO           !< temporary error for read commands

   integer(IntKi)                                  :: i                    !< generic counter
   integer(IntKi)                                  :: CurLine              !< current entry in FileInfo%Lines array
   real(ReKi),allocatable                          :: TmpRe(:)             !< temporary 2d array for reading values in
   real(ReKi)                                      :: TmpRe2(2)            !< temporary 2 number array for reading values in
   real(ReKi)                                      :: TmpRe3(3)            !< temporary 3 number array for reading values in


   ! Initialization
   ErrStat        =  0
   ErrMsg         =  ""
   UnitEcho       = -1

   ! Allocate OutList space
   CALL AllocAry( InputFileData%OutList, 18, "LidarSim Input File's OutList", TmpErrStat, TmpErrMsg ) !Max additional output parameters = 18
        if (Failed()) return;


   !-------------------------------------------------------------------------------------------------
   ! General settings
   !-------------------------------------------------------------------------------------------------

   CurLine = 4    ! Skip the first three lines as they are known to be header lines and separators
   call ParseVar( FileInfo, CurLine, 'Echo', InputFileData%Echo, TmpErrStat, TmpErrMsg )
         if (Failed()) return;

   if ( InputFileData%Echo ) then
      CALL OpenEcho ( UnitEcho, TRIM(RootName)//'.ech', TmpErrStat, TmpErrMsg )
         if (Failed()) return;
      WRITE(UnitEcho, '(A)') 'Echo file for LidarSim input file: '//trim(InputFile)
      ! Write the first three lines into the echo file
      WRITE(UnitEcho, '(A)') FileInfo%Lines(1)
      WRITE(UnitEcho, '(A)') FileInfo%Lines(2)
      WRITE(UnitEcho, '(A)') FileInfo%Lines(3)

      CurLine = 4
      call ParseVar( FileInfo, CurLine, 'Echo', InputFileData%Echo, TmpErrStat, TmpErrMsg, UnitEcho )
            if (Failed()) return
   endif

   call ParseVar( FileInfo, CurLine, 'MAXDLLChainOutputs', InputFileData%MAXDLLChainOutputs, TmpErrStat, TmpErrMsg, UnitEcho)
         if (Failed()) return;


   !-------------------------------------------------------------------------------------------------
   ! Lidar Configuration
   !-------------------------------------------------------------------------------------------------

   ! Section break
   if ( InputFileData%Echo )   WRITE(UnitEcho, '(A)') FileInfo%Lines(CurLine)    ! Write section break to echo
   CurLine = CurLine + 1

   call ParseVar( FileInfo, CurLine, 'TrajectoryType', InputFileData%TrajectoryType, TmpErrStat, TmpErrMsg, UnitEcho)
         if (Failed()) return
   call ParseVar( FileInfo, CurLine, 'WeightingType', InputFileData%WeightingType, TmpErrStat, TmpErrMsg, UnitEcho)
         if (Failed()) return
   call ParseVar( FileInfo, CurLine, 'LidarPositionX', InputFileData%LidarPositionX, TmpErrStat, TmpErrMsg, UnitEcho)
         if (Failed()) return
   call ParseVar( FileInfo, CurLine, 'LidarPositionY', InputFileData%LidarPositionY, TmpErrStat, TmpErrMsg, UnitEcho)
         if (Failed()) return
   call ParseVar( FileInfo, CurLine, 'LidarPositionZ', InputFileData%LidarPositionZ, TmpErrStat, TmpErrMsg, UnitEcho)
         if (Failed()) return
   call ParseVar( FileInfo, CurLine, 'RollAngle', InputFileData%RollAngle, TmpErrStat, TmpErrMsg, UnitEcho)
         if (Failed()) return
   call ParseVar( FileInfo, CurLine, 'PitchAngle', InputFileData%PitchAngle, TmpErrStat, TmpErrMsg, UnitEcho)
         if (Failed()) return
   call ParseVar( FileInfo, CurLine, 'YawAngle', InputFileData%YawAngle, TmpErrStat, TmpErrMsg, UnitEcho)
         if (Failed()) return
   call ParseVar( FileInfo, CurLine, 'URef', InputFileData%URef, TmpErrStat, TmpErrMsg, UnitEcho)
         if (Failed()) return
   call ParseVar( FileInfo, CurLine, 'GatesPerBeam', InputFileData%GatesPerBeam, TmpErrStat, TmpErrMsg, UnitEcho)
         if (Failed()) return
   call ParseVar( FileInfo, CurLine, 't_measurement_interval', InputFileData%t_measurement_interval, TmpErrStat, TmpErrMsg, UnitEcho)
         if (Failed()) return


   !-------------------------------------------------------------------------------------------------
   ! Cartesian coordinates
   !-------------------------------------------------------------------------------------------------

   ! Section break
   if ( InputFileData%Echo )   WRITE(UnitEcho, '(A)') FileInfo%Lines(CurLine)    ! Write section break to echo
   CurLine = CurLine + 1

   call ParseVar( FileInfo, CurLine, 'NumberOfPoints_Cartesian', InputFileData%NumberOfPoints_Cartesian, TmpErrStat, TmpErrMsg, UnitEcho)
         if (Failed()) return

   ! Section break --  X-Tab    Y-Tab    Z-Tab
   if ( InputFileData%Echo )   WRITE(UnitEcho, '(A)') '#TABLE: '//FileInfo%Lines(CurLine)    ! Write section break to echo
   CurLine = CurLine + 1

   if (InputFileData%NumberOfPoints_Cartesian > 0) then
      CALL AllocAry( InputFileData%X_Cartesian_L, InputFileData%NumberOfPoints_Cartesian, 'X Coordinate', TmpErrStat, TmpErrMsg )
            if (Failed()) return;
      CALL AllocAry( InputFileData%Y_Cartesian_L, InputFileData%NumberOfPoints_Cartesian, 'Y Coordinate', TmpErrStat, TmpErrMsg )
            if (Failed()) return;
      CALL AllocAry( InputFileData%Z_Cartesian_L, InputFileData%NumberOfPoints_Cartesian, 'Z Coordinate', TmpErrStat, TmpErrMsg )
            if (Failed()) return;

         ! TABLE read
      do i=1,InputFileData%NumberOfPoints_Cartesian
         call ParseAry ( FileInfo, CurLine, 'Coordinates', TmpRe3, 3, TmpErrStat, TmpErrMsg, UnitEcho )
               if (Failed()) return;
         InputFileData%X_Cartesian_L(i) = TmpRe3(1)
         InputFileData%Y_Cartesian_L(i) = TmpRe3(2)
         InputFileData%Z_Cartesian_L(i) = TmpRe3(3)
      enddo
   endif


   !-------------------------------------------------------------------------------------------------
   ! Spherical coordinates
   !-------------------------------------------------------------------------------------------------

   ! Section break
   if ( InputFileData%Echo )   WRITE(UnitEcho, '(A)') FileInfo%Lines(CurLine)    ! Write section break to echo
   CurLine = CurLine + 1
   call ParseVar( FileInfo, CurLine, 'NumberOfPoints_Spherical', InputFileData%NumberOfPoints_Spherical, TmpErrStat, TmpErrMsg, UnitEcho)
         if (Failed()) return

   ! Table header -- Azimuth-Tab     Elevation-Tab   RangeGates-Tab
   if ( InputFileData%Echo )   WRITE(UnitEcho, '(A)') '#TABLE: '//FileInfo%Lines(CurLine)    ! Write section break to echo
   CurLine = CurLine + 1

   if (InputFileData%NumberOfPoints_Spherical > 0 ) then
      CALL AllocAry( InputFileData%Azimuth,   InputFileData%NumberOfPoints_Spherical, 'Azimuth', TmpErrStat, TmpErrMsg )
            if (Failed()) return;
      CALL AllocAry( InputFileData%Elevation, InputFileData%NumberOfPoints_Spherical, 'Elevation', TmpErrStat, TmpErrMsg )
            if (Failed()) return;
      CALL AllocAry( InputFileData%Range,     InputFileData%NumberOfPoints_Spherical,InputFileData%GatesPerBeam , 'Range', TmpErrStat, TmpErrMsg )
            if (Failed()) return;

      ! Temporary array for reading table
      call AllocAry( TmpRe,  2+InputFileData%GatesPerBeam, 'TmpRe', TmpErrStat, TmpErrMsg )
            if (Failed()) return;

         ! TABLE read
      do i=1,InputFileData%NumberOfPoints_Spherical
         call ParseAry ( FileInfo, CurLine, 'Coordinates', TmpRe, 2+InputFileData%GatesPerBeam, TmpErrStat, TmpErrMsg, UnitEcho )
               if (Failed()) return;
         InputFileData%Azimuth(i) = TmpRe(1)
         InputFileData%Elevation(i) = TmpRe(2)
         InputFileData%Range(i,:) = TmpRe(3:size(TmpRe))
      enddo

      deallocate(TmpRe) ! Done with this temporary array
   endif


   !-------------------------------------------------------------------------------------------------
   ! Gaussian distribution
   !-------------------------------------------------------------------------------------------------

   ! Section break -- Weighting function (gaussian distribution)
   if ( InputFileData%Echo )   WRITE(UnitEcho, '(A)') FileInfo%Lines(CurLine)    ! Write section break to echo
   CurLine = CurLine + 1

   call ParseVar( FileInfo, CurLine, 'FWHM', InputFileData%FWHM, TmpErrStat, TmpErrMsg, UnitEcho)
         if (Failed()) return
   call ParseVar( FileInfo, CurLine, 'PointsToEvaluate', InputFileData%PointsToEvaluate, TmpErrStat, TmpErrMsg, UnitEcho)
         if (Failed()) return


   !-------------------------------------------------------------------------------------------------
   ! Manual distribution
   !-------------------------------------------------------------------------------------------------

   ! Section break -- Weighting function (manual weighting)
   if ( InputFileData%Echo )   WRITE(UnitEcho, '(A)') FileInfo%Lines(CurLine)    ! Write section break to echo
   CurLine = CurLine + 1
   call ParseVar( FileInfo, CurLine, 'ManualWeightingPoints', InputFileData%ManualWeightingPoints, TmpErrStat, TmpErrMsg, UnitEcho)
         if (Failed()) return

   ! Table header -- Distance-Tab   Weighting-Tab
   if ( InputFileData%Echo )   WRITE(UnitEcho, '(A)')  '#TABLE: '//FileInfo%Lines(CurLine)    ! Write section break to echo
   CurLine = CurLine + 1

   if (InputFileData%ManualWeightingPoints > 0) then
      CALL AllocAry( InputFileData%ManualWeightingDistance, InputFileData%ManualWeightingPoints, 'ManualWeightingDistance', TmpErrStat, TmpErrMsg )
            if (Failed()) return;
      CALL AllocAry( InputFileData%ManualWeighting, InputFileData%ManualWeightingPoints, 'ManualWeighting', TmpErrStat, TmpErrMsg )
            if (Failed()) return;
         ! TABLE read
      do i=1,InputFileData%NumberOfPoints_Cartesian
         call ParseAry ( FileInfo, CurLine, 'Coordinates', TmpRe3, 3, TmpErrStat, TmpErrMsg, UnitEcho )
               if (Failed()) return;
         InputFileData%ManualWeightingDistance(i) = TmpRe2(1)
         InputFileData%ManualWeighting(i) = TmpRe2(2)
      enddo
   endif


    !-------------------------------------------------------------------------------------------------
    ! Read Outlist from FileInfo
    !-------------------------------------------------------------------------------------------------

   ! Section break -- Output list
   if ( InputFileData%Echo )   WRITE(UnitEcho, '(A)') FileInfo%Lines(CurLine)    ! Write section break to echo
   CurLine = CurLine + 1
   if ( InputFileData%Echo )   WRITE(UnitEcho, '(A)') FileInfo%Lines(CurLine)    ! Write section break to echo
   CurLine = CurLine + 1

   call ReadOutputListFromFileInfo( FileInfo, CurLine, InputFileData%OutList, &
            InputFileData%NumOuts, 'OutList', "List of user-requested output channels", TmpErrStat, TmpErrMsg, UnitEcho )
         if (Failed()) return;

   CALL Cleanup()

   RETURN


CONTAINS

   !-------------------------------------------------------------------------------------------------
   logical function Failed()
      CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
      if (Failed) call Cleanup()
   end function Failed
   !-------------------------------------------------------------------------------------------------
   SUBROUTINE Cleanup()
      if (UnitEcho  > -1_IntKi)     CLOSE( UnitEcho  )
   END SUBROUTINE Cleanup
   !-------------------------------------------------------------------------------------------------

END SUBROUTINE LidarSim_ParsePrimaryFileInfo

!#########################################################################################################################################################################
    
    SUBROUTINE LidarSim_InitMeasuringPoints_Cartesian(p, InputFileData, ErrStat, ErrMsg)

    IMPLICIT                            NONE
    CHARACTER(*),                       PARAMETER       ::  RoutineName="LidarSim_InitMeasuringPoints_Cartesian"

    TYPE(LidarSim_ParameterType),       INTENT(INOUT)   ::  p                   !parameter data (destination of the InputFileData)
    TYPE(LidarSim_InputFile),           INTENT(IN   )   ::  InputFileData       !data read from the input file
    INTEGER(IntKi),                     INTENT(  OUT)   ::  ErrStat             !< Error status of the operation
    CHARACTER(*),                       INTENT(  OUT)   ::  ErrMsg              !< Error message if ErrStat /= ErrID_None

    ! Local variables
    INTEGER(IntKi)                                      ::  LoopCounter         !< counter to run through all cartesian coordinates data
    INTEGER(IntKi)                                      ::  TmpErrStat          !< Temporary error status
    CHARACTER(ErrMsgLen)                                ::  TmpErrMsg           !< temporary error message
    
    ! Initialization
    ErrStat        =  0
    ErrMsg         =  ""
    
    CALL AllocAry( p%MeasuringPoints_L, 3, SIZE(InputFileData%X_Cartesian_L), 'MeasuringPoints_L', TmpErrStat, TmpErrMsg )     !Allocate the array size for n=CounterChannelNumbers measuringpositions
    CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
    CALL AllocAry( p%MeasuringPoints_Spherical_L, 3,  SIZE(InputFileData%X_Cartesian_L), 'MeasuringPoints_Spherical_L', TmpErrStat, TmpErrMsg )
    CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
    
    DO LoopCounter = 1,SIZE(InputFileData%X_Cartesian_L)
        p%MeasuringPoints_L(1,LoopCounter) = InputFileData%X_Cartesian_L(LoopCounter)
        p%MeasuringPoints_L(2,LoopCounter) = InputFileData%Y_Cartesian_L(LoopCounter)
        p%MeasuringPoints_L(3,LoopCounter) = InputFileData%Z_Cartesian_L(LoopCounter)
        p%MeasuringPoints_Spherical_L(:,LoopCounter) = LidarSim_Cartesian2Spherical(InputFileData%X_Cartesian_L(LoopCounter),InputFileData%Y_Cartesian_L(LoopCounter),InputFileData%Z_Cartesian_L(LoopCounter))
    END DO

    END SUBROUTINE LidarSim_InitMeasuringPoints_Cartesian
   
   
!#########################################################################################################################################################################
    
    SUBROUTINE LidarSim_InitMeasuringPoints_Spherical(p, InputFileData, ErrStat, ErrMsg)

    IMPLICIT                            NONE
    CHARACTER(*),                       PARAMETER       ::  RoutineName="LidarSim_InitMeasuringPoints_Spherical"

    TYPE(LidarSim_ParameterType),       INTENT(INOUT)   ::  p                       !parameter data (destination of the InputFileData)
    TYPE(LidarSim_InputFile),           INTENT(IN   )   ::  InputFileData           !data read from the input file
    INTEGER(IntKi),                     INTENT(  OUT)   ::  ErrStat                 !< Error status of the operation
    CHARACTER(*),                       INTENT(  OUT)   ::  ErrMsg                  !< Error message if ErrStat /= ErrID_None
    
    ! Local variables
    INTEGER(IntKi)                                      ::  OuterLoopCounter        !counter for looping through the coordinate data
    INTEGER(IntKi)                                      ::  InnerLoopCounter        !counter for looping through the multiple range gates
    INTEGER(IntKi)                                      ::  CounterChannelNumbers   !counts the amount of channels
    INTEGER(IntKi)                                      ::  TmpErrStat              !< Temporary error status
    CHARACTER(ErrMsgLen)                                ::  TmpErrMsg               !< temporary error message

    ! Initialization
    ErrStat        =  0
    ErrMsg         =  ""
    
    CALL AllocAry( p%MeasuringPoints_L, 3,InputFileData%NumberOfPoints_Spherical*InputFileData%GatesPerBeam, 'MeasuringPoints_L', TmpErrStat, TmpErrMsg )     !Allocate the array size for n=CounterChannelNumbers measuringpositions
    CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
    CALL AllocAry( p%MeasuringPoints_Spherical_L, 3,InputFileData%NumberOfPoints_Spherical*InputFileData%GatesPerBeam, 'MeasuringPoints_Spherical_L', TmpErrStat, TmpErrMsg )  
    CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)

    OuterLoopCounter = 1
    InnerLoopCounter = 1
    CounterChannelNumbers = 1    
        
    DO WHILE( OuterLoopCounter <= InputFileData%NumberOfPoints_Spherical)
        DO WHILE( InnerLoopCounter <= InputFileData%GatesPerBeam)
                IF ( InputFileData%Range(OuterLoopCounter,InnerLoopCounter) /= 0 )THEN  !Range gates mustn't be 0. => Divide by 0 !
                    p%MeasuringPoints_L(:,CounterChannelNumbers) = &   !Transformation from the spherical to cartesian coordinates
                        LidarSim_Spherical2Cartesian(InputFileData%Azimuth(OuterLoopCounter),InputFileData%Elevation(OuterLoopCounter),InputFileData%Range(OuterLoopCounter,InnerLoopCounter))
                    p%MeasuringPoints_Spherical_L(:,CounterChannelNumbers) = &
                        (/ InputFileData%Range(OuterLoopCounter,InnerLoopCounter),InputFileData%Azimuth(OuterLoopCounter),InputFileData%Elevation(OuterLoopCounter)/)
                    CounterChannelNumbers = CounterChannelNumbers + 1
                ELSE
                    CALL SetErrStat(ErrID_Fatal,"Range gates must not be 0",ErrStat,ErrMsg,RoutineName)
                ENDIF
            InnerLoopCounter = InnerLoopCounter + 1
        END DO
        InnerLoopCounter = 1
        OuterLoopCounter = OuterLoopCounter + 1
    END DO

    END SUBROUTINE LidarSim_InitMeasuringPoints_Spherical


!#########################################################################################################################################################################
    
    FUNCTION LidarSim_Spherical2Cartesian(Azimuth, Elevation, Range)

    IMPLICIT        NONE
    CHARACTER(*),   PARAMETER       ::  RoutineName="LidarSim_Spherical2Cartesian"

    REAL(ReKi),     INTENT(IN   )   ::  Azimuth                             !Azimuth angle
    REAL(ReKi),     INTENT(IN   )   ::  Elevation                           !Elevation angle
    REAL(ReKi),     INTENT(IN   )   ::  Range                               !range gate
    REAL(ReKi),     DIMENSION (3)   ::  LidarSim_Spherical2Cartesian        !Output : x,y,z coordinates
    
    LidarSim_Spherical2Cartesian(1)  =   Range*COS(Elevation)*COS(Azimuth)   !x
    LidarSim_Spherical2Cartesian(2)  =   Range*COS(Elevation)*SIN(Azimuth)   !y
    LidarSim_Spherical2Cartesian(3)  =   Range*SIN(Elevation)                !z

    END FUNCTION LidarSim_Spherical2Cartesian
    

!#########################################################################################################################################################################    
    
    FUNCTION LidarSim_Cartesian2Spherical(X, Y, Z)

    IMPLICIT        NONE
    CHARACTER(*),   PARAMETER       ::  RoutineName="LidarSim_Cartesian2Spherical"

    REAL(ReKi),     INTENT(IN   )   ::  X                             !Azimuth angle
    REAL(ReKi),     INTENT(IN   )   ::  Y                           !Elevation angle
    REAL(ReKi),     INTENT(IN   )   ::  Z                               !range gate
    REAL(ReKi),     DIMENSION (3)   ::  LidarSim_Cartesian2Spherical       !Output : x,y,z coordinates

    LidarSim_Cartesian2Spherical(1)  =   SQRT(X**2+Y**2+Z**2)   !range
    
    IF(LidarSim_Cartesian2Spherical(1) > 0 ) THEN
        LidarSim_Cartesian2Spherical(3)  =   ASIN(Z/LidarSim_Cartesian2Spherical(1))   !y
        
        LidarSim_Cartesian2Spherical(2) = ATAN2(Y,X)
    ELSE
        LidarSim_Cartesian2Spherical(2) = 0
        LidarSim_Cartesian2Spherical(3) = 0
    ENDIF

    END FUNCTION LidarSim_Cartesian2Spherical
    

!#########################################################################################################################################################################
    
    FUNCTION LidarSim_TransformLidarToInertial(LidarUnitMotion, p, MeasuringPoint_L)

    IMPLICIT        NONE
    CHARACTER(*),   PARAMETER                      ::  RoutineName="LidarSim_TransformLidarToInertial"

    REAL(ReKi)                                     ::  LidarSim_TransformLidarToInertial(3)  !Output calculated transformation from the lidar coord. sys. to the inertial system
    TYPE(MeshType),                 intent(in   )  ::  LidarUnitMotion                         !Data describing the motion of the nacelle coord. sys.
    TYPE(LidarSim_ParameterType),   intent(in   )  ::  p                                     !Parameter data 
    REAL(ReKi),                     intent(in   )  ::  MeasuringPoint_L(3)                   !point which needs to be transformed
    
      ! NOTE:  LidarOrientation follows all orientation and position changes of what it is attached to.  This mesh already includes
      !        the orientation relative to what it is attached to
    LidarSim_TransformLidarToInertial = LidarUnitMotion%Position(:,1) + LidarUnitMotion%TranslationDisp(:,1) &
                                       + MATMUL(TRANSPOSE(LidarUnitMotion%Orientation(:,:,1)),(MeasuringPoint_L ) )
 
    END FUNCTION LidarSim_TransformLidarToInertial
    

!#########################################################################################################################################################################
    
    SUBROUTINE LidarSim_InitializeWeightingManual(p, InputFileData, ErrStat, ErrMsg)
    
    IMPLICIT                        NONE
    CHARACTER(*),                   PARAMETER       ::  RoutineName="LidarSim_InitializeWeightingManual"
    
    TYPE(LidarSim_ParameterType),   INTENT(INOUT)   ::  p                   ! parameter data to write results in
    TYPE(LidarSim_InputFile),       INTENT(IN   )   ::  InputFileData       ! Inputdata from the input file
    INTEGER(IntKi),                 INTENT(  OUT)   ::  ErrStat             !< Temporary error status
    CHARACTER(*),                   INTENT(  OUT)   ::  ErrMsg              !< temporary error message
    
    ! local variables
    INTEGER(IntKi)          ::  TmpErrStat
    CHARACTER(ErrMsgLen)    ::  TmpErrMsg           !< temporary error message
    
    TmpErrStat = 0
    TmpErrMsg = ''

    CALL AllocAry( p%WeightingDistance,SIZE(InputFileData%ManualWeightingDistance), 'p%WeightingDistance', TmpErrStat, TmpErrMsg )  !Allocating the needed size for the weighting distance vector
    CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
    CALL AllocAry( p%Weighting,SIZE(InputFileData%ManualWeighting), 'p%Weighting', TmpErrStat, TmpErrMsg )                          !Allocating the needed size for the weighting vector
    CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
    p%WeightingDistance = InputFileData%ManualWeightingDistance                                                                     !writing the input distances in the parameter
    p%Weighting(:) = InputFileData%ManualWeighting(:) / SUM(InputFileData%ManualWeighting)                                          !writing the input weighting in the parameter and normalizing it (sum = 1)
    
    END SUBROUTINE LidarSim_InitializeWeightingManual
    
    
!#########################################################################################################################################################################    

    SUBROUTINE LidarSim_InitializeWeightingGauss(p, InputFileData, ErrStat, ErrMsg)
    
    IMPLICIT                        NONE
    CHARACTER(*),                   PARAMETER       ::  RoutineName="LidarSim_InitializeWeightingGauss"
    
    TYPE(LidarSim_ParameterType),   INTENT(INOUT)   ::  p                   ! parameter data to write results in
    TYPE(LidarSim_InputFile),       INTENT(IN   )   ::  InputFileData       ! Inputdata from the input file
    INTEGER(IntKi),                 INTENT(  OUT)   ::  ErrStat             !< Error status of the operation
    CHARACTER(*),                   INTENT(  OUT)   ::  ErrMsg              !< Error message if ErrStat /= ErrID_None
    
    ! local variables
    INTEGER(IntKi)          ::  Counter                                     !Loopcounter for every point to evaluate
    REAL(ReKi)              ::  Dist                                        !Distance between each evaluation point
    INTEGER(IntKi)          ::  TmpErrStat                                  !< Temporary error status
    CHARACTER(ErrMsgLen)    ::  TmpErrMsg                                   !< temporary error message
    
    ErrStat =   0
    ErrMsg  =   ''
    Dist = 2*InputFileData%FWHM/(InputFileData%PointsToEvaluate+1)          !Distance between each weighting point
    
    CALL AllocAry( p%WeightingDistance,InputFileData%PointsToEvaluate, 'p%WeightingDistance', TmpErrStat, TmpErrMsg )   !Allocating the needed size for the weighting distance vector
    CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
    CALL AllocAry( p%Weighting,InputFileData%PointsToEvaluate, 'p%Weighting', TmpErrStat, TmpErrMsg )                   !Allocating the needed size for the weighting vector
    CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
    
    DO Counter=1,InputFileData%PointsToEvaluate
        p%WeightingDistance(Counter) = (Counter) * Dist + (-InputFileData%FWHM)             !Creating the distance vector
        p%Weighting(Counter) = ( (2*SQRT(LOG(2.0)))/(InputFileData%FWHM*SQRT(Pi_D)) ) *&      !Calculation of the gaussian distribution
        EXP( -((p%WeightingDistance(Counter)**2)*4*LOG(2.0))/(InputFileData%FWHM**2) ) !&
    END DO
    p%Weighting = p%Weighting / SUM( p%Weighting )

    END SUBROUTINE LidarSim_InitializeWeightingGauss
    
    
!#########################################################################################################################################################################   
    
    SUBROUTINE LidarSim_CalculateVlos(p, UnitVector_I, Vlos, MeasuringPosition_I, LidarPosition_I,&
    Time, IfW_p, IfW_ContStates, IfW_DiscStates, IfW_ConstrStates, IfW_OtherStates, IfW_m, ErrStat, ErrMsg)
    
    IMPLICIT                                NONE
    CHARACTER(*),                           PARAMETER       ::  RoutineName="LidarSim_CalculateVlos"
    
    TYPE(LidarSim_ParameterType),           INTENT(IN   )   ::  p                           !parameter data of the lidar module
    REAL(ReKi),                             INTENT(IN   )   ::  UnitVector_I(3)             !Line of Sight Unit Vector
    REAL(ReKi),                             INTENT(INOUT)   ::  MeasuringPosition_I(3)      !Position of the measuring point
    REAL(ReKi),                             INTENT(IN   )   ::  LidarPosition_I(3)          !Position of the measuring point
    REAL(ReKi),                             INTENT(  OUT)   ::  Vlos                        !Calculated speed in los direction
    REAL(DbKi),                             INTENT(IN   )   ::  Time                        !< Current simulation time in seconds

    !Error Variables
    INTEGER(IntKi),                         INTENT(  OUT)   ::  ErrStat                     !< Error status of the operation
    CHARACTER(*),                           INTENT(  OUT)   ::  ErrMsg                      !< Error message if ErrStat /= ErrID_None
    
!FIXME: remove IfW completely from here
    !IfW Parameter
    TYPE(InflowWind_ParameterType),         INTENT(IN   )   ::  IfW_p                       !< Parameters
    TYPE(InflowWind_ContinuousStateType),   INTENT(IN   )   ::  IfW_ContStates              !< Continuous states at Time
    TYPE(InflowWind_DiscreteStateType),     INTENT(IN   )   ::  IfW_DiscStates              !< Discrete states at Time
    TYPE(InflowWind_ConstraintStateType),   INTENT(IN   )   ::  IfW_ConstrStates            !< Constraint states at Time
    TYPE(InflowWind_OtherStateType),        INTENT(IN   )   ::  IfW_OtherStates             !< Other/optimization states at Time
    TYPE(InflowWind_MiscVarType),           INTENT(INOUT)   ::  IfW_m                       !< Misc variables for optimization (not copied in glue code)
    
    !Local Variables
    TYPE(InflowWind_InputType)              ::  InputForCalculation                         ! input data field for the calculation in the InflowWind module
    TYPE(InflowWind_OutputType)             ::  OutputForCalculation                        ! data field were the calculated speed is written from the InflowWind module
    INTEGER(IntKi)                          ::  Counter                                     ! Counter for the loop for the different weightings of the point
    REAL(ReKi),DIMENSION(:), ALLOCATABLE    ::  Vlos_tmp                                    !< Array with all temporary Vlos
    
    ! Temporary variables for error handling
    INTEGER(IntKi)                                          ::  ErrStat2        
    CHARACTER(ErrMsgLen)                                    ::  ErrMsg2

    !Initialize error values
    ErrStat        =  ErrID_None
    ErrMsg         =  ""
    
    CALL AllocAry(InputForCalculation%PositionXYZ, 3,1, 'InputForCalculation%PositionXYZ',ErrStat2, ErrMsg2)        !Allocating needed space for the input
    CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
    CALL AllocAry(OutputForCalculation%VelocityUVW, 3,1, 'OutputForCalculation%VelocityUVW',ErrStat2, ErrMsg2)      !Allocating needed space for the output
    CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
    CALL AllocAry(Vlos_tmp ,SIZE(p%Weighting), 'Vlos_tmp%VelocityUVW',ErrStat2, ErrMsg2)                            !Allocating space for temporary windspeeds
    CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
     
    IF(IfW_p%WindType == 1 .OR. IfW_p%WindType == 2)Then !Uniform Wind 2 (und steady 1)
        MeasuringPosition_I(1) = MeasuringPosition_I(1)-LidarPosition_I(1)  !In the uniform wind case. the wind hits the turbine at the same time indepentend of the x shift
        DO Counter = 1, SIZE(p%Weighting)
!QUESTION: is the p%WeightingDistance about the measurement position, or from the lidar?? If the latter, there is a problem with the below calculations
            ! position of the weighted measuring point
            InputForCalculation%PositionXYZ(:,1) = MeasuringPosition_I + p%WeightingDistance(Counter) * UnitVector_I
!FIXME: Cannot call InflowWind directly like this.  This is not allowed by the framework.
            CALL CalculateOutput(Time + DBLE(-InputForCalculation%PositionXYZ(1,1)/p%Uref), InputForCalculation, IfW_p, &
                                    IfW_ContStates, IfW_DiscStates, IfW_ConstrStates, IfW_OtherStates, OutputForCalculation, &
                                    IfW_m, .FALSE., ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            Vlos_tmp(Counter) = - DOT_PRODUCT(OutputForCalculation%VelocityUVW(:,1),UnitVector_I)
        END DO
    ELSE IF(IfW_p%WindType ==  3 .OR. IfW_p%WindType == 4) THEN        !Bladed Turublent 4 ( und TurbSim 3)
        DO Counter = 1, SIZE(p%Weighting)
            ! position of the weighted measuring point
            InputForCalculation%PositionXYZ(:,1) = MeasuringPosition_I + p%WeightingDistance(Counter) * UnitVector_I
!FIXME: Cannot call InflowWind directly like this.  This is not allowed by the framework.
            CALL CalculateOutput(Time, InputForCalculation, IfW_p, IfW_ContStates, IfW_DiscStates, IfW_ConstrStates, &
                                       IfW_OtherStates, OutputForCalculation, IfW_m, .FALSE., ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            Vlos_tmp(Counter) = - DOT_PRODUCT(OutputForCalculation%VelocityUVW(:,1),UnitVector_I)
        END DO
    END IF
    Vlos = DOT_PRODUCT(Vlos_tmp, p%Weighting)           !Calculation of the weighted windspeed

    DEALLOCATE (InputForCalculation%PositionXYZ)        !Free Input Positions for the next measurement
    DEALLOCATE (OutputForCalculation%VelocityUVW)       !Free Ouput Positions for the next measurement
    DEALLOCATE (Vlos_tmp)

END SUBROUTINE LidarSim_CalculateVlos
    
    
!#########################################################################################################################################################################       
    
SUBROUTINE LidarSim_CalculateIMU(p,y,u)
!FIXME: check this against 133f6a828599671b80b8b502b7974a302a8180ec
!        There may have been some inadvertant change of the coordinate frame.  It is now returning the output in global, but local maybe be what is wanted!!!!

   IMPLICIT                                   NONE
   CHARACTER(*),                              PARAMETER            ::  RoutineName="LidarSim_LidarSim_CalculateIMU"

   TYPE(LidarSim_ParameterType),              INTENT(IN   )        ::  p
   TYPE(LidarSim_OutputType),                 INTENT(INOUT)        ::  y                       !Outputs computed at Time (IN for mesh reasons and data allocation)
   TYPE(LidarSim_InputType),                  INTENT(IN   )        ::  u
 
   ! local variables
   real(R8Ki)        :: theta(3)       ! roll, pitch, yaw
   real(R8Ki)                                                      ::  Rotation_L_I(3,3)
!    REAL(ReKi)                                                      ::  CrossProduct(3)         !Variable for the crossproduct of the rotation and the lidar position ( in the nacelle coord. system)
!    REAL(ReKi)                                                      ::  Roll
!    REAL(ReKi)                                                      ::  Pitch
!    REAL(ReKi)                                                      ::  Yaw
!    REAL(ReKi)                                                      ::  DisplacementNacelle(3)
!    REAL(ReKi)                                                      ::  DisplacementLidar(3)
!    REAL(ReKi)                                                      ::  LidarPosition_I(3)


   Rotation_L_I = MATMUL(u%LidarMesh%Orientation(:,:,1),transpose(u%LidarMesh%RefOrientation(:,:,1)))

   theta = EulerExtract( Rotation_L_I )

!FIXME: check which coordinate system the outputs from this should be in
!NOTE: if we want this in the coordinate system of the Lidar module, we will need to multiply by the u%LidarMesh%Orientation matrix
   y%IMUOutputs(1)     = theta(1)                      !Roll Angle
   y%IMUOutputs(2)     = u%LidarMesh%RotationVel(1,1)  !Roll Angle Velocity
   y%IMUOutputs(3)     = u%LidarMesh%RotationAcc(1,1)  !Roll Angle Acceleration
   y%IMUOutputs(4)     = theta(2)                      !Pitch Angle
   y%IMUOutputs(5)     = u%LidarMesh%RotationVel(2,1)  !Pitch Angle Velocity
   y%IMUOutputs(6)     = u%LidarMesh%RotationAcc(2,1)  !Pitch Angle Acceleration
   y%IMUOutputs(7)     = theta(3)                      !Yaw Angle
   y%IMUOutputs(8)     = u%LidarMesh%RotationVel(3,1)  !Yaw Angle Velocity
   y%IMUOutputs(9)     = u%LidarMesh%RotationAcc(3,1)  !Yaw Angle Acceleration

   y%IMUOutputs(10)    = u%LidarMesh%TranslationDisp(1,1)   !Displacement x 
   y%IMUOutputs(13)    = u%LidarMesh%TranslationDisp(2,1)   !Displacement y
   y%IMUOutputs(16)    = u%LidarMesh%TranslationDisp(3,1)   !Displacement z

   y%IMUOutputs(11)    = u%LidarMesh%TranslationVel(1,1)    !Velocity x
   y%IMUOutputs(14)    = u%LidarMesh%TranslationVel(2,1)    !Velocity y
   y%IMUOutputs(17)    = u%LidarMesh%TranslationVel(3,1)    !Velocity z

   y%IMUOutputs(12)    = u%LidarMesh%TranslationAcc(1,1)    !Acceleration x
   y%IMUOutputs(15)    = u%LidarMesh%TranslationAcc(2,1)    !Acceleration y
   y%IMUOutputs(18)    = u%LidarMesh%TranslationAcc(3,1)    !Acceleration z

   
!   IF(.NOT.(Rotation_L_I(3,1) == 1 .OR. Rotation_L_I(3,1) == -1)) THEN
!       Pitch = -ASIN(Rotation_L_I(3,1))
!       Roll = ATAN2(Rotation_L_I(3,2)/cos(Pitch),Rotation_L_I(3,3)/cos(Pitch))
!       Yaw = ATAN2(Rotation_L_I(2,1)/cos(Pitch), Rotation_L_I(1,1)/cos(Pitch))  
!   ELSE
!       Yaw = 0
!       IF(Rotation_L_I(3,1) == 1) THEN
!           Pitch = PiBy2_D
!           Roll = Yaw + ATAN2(Rotation_L_I(1,2),Rotation_L_I(1,3))
!       ELSE
!           Pitch = -PiBy2_D
!           Roll = -Yaw + ATAN2(-Rotation_L_I(1,2),-Rotation_L_I(1,3))
!       END IF
!   END IF
!
!    LidarPosition_I = MATMUL(TRANSPOSE(u%NacelleMotion%Orientation(:,:,1)),p%LidarPosition_N) !Rotates the position vector to the orientation of the inertial coord. system
!
!   y%IMUOutputs(1)     = Roll                              !Roll Angle
!   y%IMUOutputs(2)     = u%NacelleMotion%RotationVel(1,1)  !Roll Angle Velocity
!   y%IMUOutputs(3)     = u%NacelleMotion%RotationAcc(1,1)  !Roll Angle Acceleration
!   y%IMUOutputs(4)     = Pitch                             !Pitch Angle
!   y%IMUOutputs(5)     = u%NacelleMotion%RotationVel(2,1)  !Pitch Angle Velocity
!   y%IMUOutputs(6)     = u%NacelleMotion%RotationAcc(2,1)  !Pitch Angle Acceleration
!   y%IMUOutputs(7)     = Yaw                               !Yaw Angle
!   y%IMUOutputs(8)     = u%NacelleMotion%RotationVel(3,1)  !Yaw Angle Velocity
!   y%IMUOutputs(9)     = u%NacelleMotion%RotationAcc(3,1)  !Yaw Angle Acceleration
!   
!   y%IMUOutputs(10)    = u%NacelleMotion%TranslationDisp(1,1)  !Displacement x 
!   y%IMUOutputs(13)    = u%NacelleMotion%TranslationDisp(2,1)  !Displacement y
!   y%IMUOutputs(16)    = u%NacelleMotion%TranslationDisp(3,1)  !Displacement z
!   
!   DisplacementNacelle(1) = u%NacelleMotion%TranslationDisp(1,1)
!   DisplacementNacelle(2) = u%NacelleMotion%TranslationDisp(2,1)
!   DisplacementNacelle(3) = u%NacelleMotion%TranslationDisp(3,1)
!   
!   DisplacementLidar   = LidarPosition_I + DisplacementNacelle
!   
!   y%IMUOutputs(10)    = DisplacementLidar(1)
!   y%IMUOutputs(13)    = DisplacementLidar(2)
!   y%IMUOutputs(16)    = DisplacementLidar(3)
!   
!   CrossProduct(1)     = u%NacelleMotion%RotationVel(2,1)*LidarPosition_I(3) - u%NacelleMotion%RotationVel(3,1)*LidarPosition_I(2)
!   CrossProduct(2)     = u%NacelleMotion%RotationVel(3,1)*LidarPosition_I(1) - u%NacelleMotion%RotationVel(1,1)*LidarPosition_I(3)
!   CrossProduct(3)     = u%NacelleMotion%RotationVel(1,1)*LidarPosition_I(2) - u%NacelleMotion%RotationVel(2,1)*LidarPosition_I(1)
!   
!   y%IMUOutputs(11)    = u%NacelleMotion%TranslationVel(1,1) + CrossProduct(1)    !Velocity x
!   y%IMUOutputs(14)    = u%NacelleMotion%TranslationVel(2,1) + CrossProduct(2)    !Velocity y
!   y%IMUOutputs(17)    = u%NacelleMotion%TranslationVel(3,1) + CrossProduct(3)    !Velocity z
!   
!   CrossProduct(1)     = u%NacelleMotion%RotationAcc(2,1)*LidarPosition_I(3) - u%NacelleMotion%RotationAcc(3,1)*LidarPosition_I(2)
!   CrossProduct(2)     = u%NacelleMotion%RotationAcc(3,1)*LidarPosition_I(1) - u%NacelleMotion%RotationAcc(1,1)*LidarPosition_I(3)
!   CrossProduct(3)     = u%NacelleMotion%RotationAcc(1,1)*LidarPosition_I(2) - u%NacelleMotion%RotationAcc(2,1)*LidarPosition_I(1)    
!   
!   y%IMUOutputs(12)    = u%NacelleMotion%TranslationAcc(1,1) + CrossProduct(1)     !Acceleration x
!   y%IMUOutputs(15)    = u%NacelleMotion%TranslationAcc(2,1) + CrossProduct(2)     !Acceleration y
!   y%IMUOutputs(18)    = u%NacelleMotion%TranslationAcc(3,1) + CrossProduct(3)     !Acceleration z

END SUBROUTINE
    
    
!#########################################################################################################################################################################       
    
SUBROUTINE LidarSim_InitializeOutputs(y,p, InitOutData, InputFileData, ErrStat, ErrMsg)
    
    IMPLICIT                               NONE
    CHARACTER(*),                          PARAMETER       ::  RoutineName='LidarSim_InitializeOutputs'
    
    TYPE(LidarSim_OutputType),             INTENT(INOUT)   ::  y                   ! lidar outputs (IN since contains mesh info already)
    TYPE(LidarSim_ParameterType),          INTENT(INOUT)   ::  p                   ! Parameter data for the lidar module
    TYPE(LidarSim_InitOutputType),         INTENT(  OUT)   ::  InitOutData
    TYPE(LidarSim_InputFile),              INTENT(IN   )   ::  InputFileData
    INTEGER(IntKi),                        INTENT(  OUT)   ::  ErrStat              !< Temporary error status
    CHARACTER(ErrMsgLen),                  INTENT(  OUT)   ::  ErrMsg               !< temporary error message
    
    !Local error variables
    INTEGER(IntKi)          ::  TmpErrStat
    CHARACTER(ErrMsgLen)    ::  TmpErrMsg           !< temporary error message
    
    !Local variables
    INTEGER(IntKi)                          ::  SizeOutput
    INTEGER(IntKi)                          ::  LoopCounter
    INTEGER(IntKi)                          ::  NumberValidOutputs
    INTEGER(IntKi),DIMENSION(:),ALLOCATABLE ::  ValidOutputs
    CHARACTER(15) ,DIMENSION(:),ALLOCATABLE ::  ValidOutputNames
    CHARACTER(1024)                         ::  TmpIntegerToString
    CHARACTER(15)                           ::  OutListTmp
    
    INTEGER(IntKi), PARAMETER               ::  ROLLLI      =   1
    INTEGER(IntKi), PARAMETER               ::  ROLLDTLI    =   2
    INTEGER(IntKi), PARAMETER               ::  ROLLDTDTLI  =   3
    INTEGER(IntKi), PARAMETER               ::  PTCHLI      =   4
    INTEGER(IntKi), PARAMETER               ::  PTCHDTLI    =   5
    INTEGER(IntKi), PARAMETER               ::  PTCHDTDTLI  =   6
    INTEGER(IntKi), PARAMETER               ::  YAWLI       =   7
    INTEGER(IntKi), PARAMETER               ::  YAWDTLI     =   8
    INTEGER(IntKi), PARAMETER               ::  YAWDTDTLI   =   9
    INTEGER(IntKi), PARAMETER               ::  XLI         =   10
    INTEGER(IntKi), PARAMETER               ::  XDTLI       =   11
    INTEGER(IntKi), PARAMETER               ::  XDTDTLI     =   12
    INTEGER(IntKi), PARAMETER               ::  YLI         =   13
    INTEGER(IntKi), PARAMETER               ::  YDTLI       =   14
    INTEGER(IntKi), PARAMETER               ::  YDTDTLI     =   15
    INTEGER(IntKi), PARAMETER               ::  ZLI         =   16
    INTEGER(IntKi), PARAMETER               ::  ZDTLI       =   17
    INTEGER(IntKi), PARAMETER               ::  ZDTDTLI     =   18
    INTEGER(IntKi), PARAMETER               ::  AZIMUTHLI   =   19
    INTEGER(IntKi), PARAMETER               ::  ELEVATLI    =   20
    INTEGER(IntKi), PARAMETER               ::  MEASTIMELI  =   21
    INTEGER(IntKi), PARAMETER               ::  NEWDATALI   =   22
    INTEGER(IntKi), PARAMETER               ::  BEAMIDLI    =   23
        
    CHARACTER(ChanLen), DIMENSION(:), ALLOCATABLE  :: ValidParamAry     ! This lists the names of the allowed parameters, which must be sorted alphabetically
    INTEGER(IntKi), DIMENSION(:), ALLOCATABLE  :: ParamIndexAry     ! This lists the index into AllOuts(:) of the allowed parameters ValidParamAry(:)
    CHARACTER(ChanLen), DIMENSION(:), ALLOCATABLE  :: ParamUnitsAry     ! This lists the units corresponding to the allowed parameters
    
    !Error Variables
    TmpErrStat  = 0
    TmpErrMsg   = ''
    
    SizeOutput = InputFileData%GatesPerBeam
        
    CALL AllocAry( ValidParamAry, 23+(2*SizeOutput), 'ValidParamAry', TmpErrStat, TmpErrMsg )
    CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
    CALL AllocAry( ParamIndexAry, 23+(2*SizeOutput), 'ParamIndexAry', TmpErrStat, TmpErrMsg )
    CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
    CALL AllocAry( ParamUnitsAry, 23+(2*SizeOutput), 'ParamUnitsAry', TmpErrStat, TmpErrMsg )
    CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
    
    ! Fill ValidParamAry and ParamIndexAry in alphabetical order
    ValidParamAry( 1 : 8 ) = (/ & 
        "AZIMUTHLI ","BEAMIDLI  ","ELEVATLI  ","MEASTIMELI","NEWDATALI ","PTCHDTDTLI", & 
        "PTCHDTLI  ","PTCHLI    "/)
    ParamIndexAry( 1 : 8 ) = (/ & 
        AZIMUTHLI   ,BEAMIDLI    ,ELEVATLI    ,MEASTIMELI  ,NEWDATALI   ,PTCHDTDTLI  , & 
        PTCHDTLI    ,PTCHLI      /)
    DO LoopCounter = 1,SizeOutput
        WRITE(UNIT=TmpIntegerToString,FMT='(I2.2)') LoopCounter
        ValidParamAry( 8 + LoopCounter ) = "RANGE"//TRIM(ADJUSTL(TmpIntegerToString))//"LI"
        ParamIndexAry( 8 + LoopCounter ) = 23 + LoopCounter
    ENDDO
    ValidParamAry( (9 + SizeOutput) : (11 + SizeOutput) ) = (/ & 
        "ROLLDTDTLI","ROLLDTLI  ","ROLLLI    " /)
    ParamIndexAry( (9 + SizeOutput) : (11 + SizeOutput) ) = (/ & 
        ROLLDTDTLI  ,ROLLDTLI    ,ROLLLI      /)
    DO LoopCounter = 1,SizeOutput
        WRITE(UNIT=TmpIntegerToString,FMT='(I2.2)') LoopCounter
        ValidParamAry( 11 + SizeOutput + LoopCounter ) = "VLOS"//TRIM(ADJUSTL(TmpIntegerToString))//"LI"
        ParamIndexAry( 11 + SizeOutput + LoopCounter ) = 23 + SizeOutput + LoopCounter
    ENDDO
    ValidParamAry( (12 + (2*SizeOutput)) : (23 + (2*SizeOutput)) ) = (/ & 
        "XDTDTLI   ","XDTLI     ","XLI       ","YAWDTDTLI ","YAWDTLI   ","YAWLI     ", & 
        "YDTDTLI   ","YDTLI     ","YLI       ","ZDTDTLI   ","ZDTLI     ","ZLI       "/)
    ParamIndexAry( (12 + (2*SizeOutput)) : (23 + (2*SizeOutput)) ) = (/ & 
        XDTDTLI     ,XDTLI       ,XLI         ,YAWDTDTLI   ,YAWDTLI     ,YAWLI       , & 
        YDTDTLI     ,YDTLI       ,YLI         ,ZDTDTLI     ,ZDTLI       ,ZLI         /)
    
    ! Fill ParamUnitsAry according to parameter order
    ParamUnitsAry( 1 : 23 ) = (/ & 
        "(rad)     ","(rad/s)   ","(rad/s^2) ","(rad)     ","(rad/s)   ","(rad/s^2) ", & 
        "(rad)     ","(rad/s)   ","(rad/s^2) ","(m)       ","(m/s)     ","(m/s^2)   ", & 
        "(m)       ","(m/s)     ","(m/s^2)   ","(m)       ","(m/s)     ","(m/s^2)   ", & 
        "(rad)     ","(rad)     ","(s)       ","()        ","()        "/)
    DO LoopCounter = 1,SizeOutput
        ParamUnitsAry( 23 + LoopCounter ) = "(m)       "
        ParamUnitsAry( 23 + SizeOutput + LoopCounter ) = "(m/s)     "
    ENDDO
    
    CALL AllocAry( y%IMUOutputs, 18, 'IMUOutputs', TmpErrStat, TmpErrMsg )                              !Allocate array for all IMU data ( rotation and position displacement, velocity, acceleration)
    CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
    CALL AllocAry( ValidOutputs, 23+(2*SizeOutput), 'ValidOutputs', TmpErrStat, TmpErrMsg )             !Allocate the maximal additional outputchannels
    CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
    CALL AllocAry( ValidOutputNames, 23+(2*SizeOutput), 'ValidOutputNames', TmpErrStat, TmpErrMsg )     !Allocate the maximal additional outputchannels
    CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
    
    NumberValidOutputs=0
    DO LoopCounter = 1, InputFileData%NumOuts
        OutListTmp = InputFileData%OutList(LoopCounter)
        CALL Conv2UC(OutListTmp)
        ValidOutputs(NumberValidOutputs+1) = IndexCharAry(OutListTmp,ValidParamAry)
        IF(ValidOutputs(NumberValidOutputs+1) > 0) THEN
            ValidOutputs(NumberValidOutputs+1)      =   ParamIndexAry(ValidOutputs(NumberValidOutputs+1))
            ValidOutputNames(NumberValidOutputs+1)  =   InputFileData%OutList(LoopCounter)
            NumberValidOutputs = NumberValidOutputs + 1            
        ENDIF
    ENDDO
    
    IF(NumberValidOutputs>0) THEN
        CALL AllocAry( p%ValidOutputs, NumberValidOutputs, 'p%ValidOutputs', TmpErrStat, TmpErrMsg )                     !Allocate the fitting amount of outputchannels
        CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
        p%ValidOutputs=ValidOutputs(1:NumberValidOutputs)   
    ENDIF
    
    CALL AllocAry( y%WriteOutput, NumberValidOutputs, 'WriteOutput', TmpErrStat, TmpErrMsg )                     !Allocate the actual output array
    CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
    CALL AllocAry( InitOutData%WriteOutputHdr, NumberValidOutputs, 'WriteOutputHdr', TmpErrStat, TmpErrMsg )     !Name of the data output channels   (Size of the WriteOutputHdr array defines the number of outputs)
    CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
    CALL AllocAry( InitOutData%WriteOutputUnt, NumberValidOutputs, 'WriteOutputUnt', TmpErrStat, TmpErrMsg )     !units of the output channels
    CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
    
    DO LoopCounter = 1, NumberValidOutputs
         y%WriteOutput( LoopCounter ) = 0
         InitOutData%WriteOutputHdr( LoopCounter ) = ValidOutputNames(LoopCounter)
         InitOutData%WriteOutputUnt( LoopCounter ) = ParamUnitsAry(p%ValidOutputs(LoopCounter))
    ENDDO
    
    DEALLOCATE(ValidOutputNames)
    DEALLOCATE(ValidOutputs)
    
    ! Initialize SwapOutputs array
    CALL AllocAry( y%SwapOutputs, 2+SizeOutput+6, 'SwapOutputs', TmpErrStat, TmpErrMsg )
    CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
    
    y%SwapOutputs(1) = 0                        !NewData
    y%SwapOutputs(2) = 0                        !BeamID
    
    DO LoopCounter = 1,SizeOutput
         y%SwapOutputs(2 + LoopCounter) = 0     !Vlos
    ENDDO

    y%SwapOutputs(2 + SizeOutput + 1)   = 0     ! LdrRoll
    y%SwapOutputs(2 + SizeOutput + 2)   = 0     ! LdrPitch
    y%SwapOutputs(2 + SizeOutput + 3)   = 0     ! LdrYaw
    y%SwapOutputs(2 + SizeOutput + 4)   = 0     ! LdrXd
    y%SwapOutputs(2 + SizeOutput + 5)   = 0     ! LdrYd 
    y%SwapOutputs(2 + SizeOutput + 6)   = 0     ! LdrZd 
    
    ! Initialize AllOutputs array
    CALL AllocAry( y%AllOutputs, 23+(2*SizeOutput), 'AllOutputs', TmpErrStat, TmpErrMsg )
    CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
    
    y%AllOutputs = 0
    
END SUBROUTINE LidarSim_InitializeOutputs
    
    
!#########################################################################################################################################################################  
    
SUBROUTINE LidarSim_SetOutputs(y,p,m,Vlos,UnitVector,LoopGatesPerBeam,Time)
    
    IMPLICIT                                    NONE
    CHARACTER(*),                               PARAMETER       ::  RoutineName='LidarSim_SetOutputs'
    
    TYPE(LidarSim_OutputType),                  INTENT(INOUT)   ::  y
    TYPE(LidarSim_ParameterType),               INTENT(IN   )   ::  p
    TYPE(LidarSim_MiscVarType),                 INTENT(IN   )   ::  m
    REAL(ReKi),                                 INTENT(IN   )   ::  Vlos
    REAL(ReKi),                                 INTENT(IN   )   ::  UnitVector(3)
    INTEGER(IntKi),                             INTENT(IN   )   ::  LoopGatesPerBeam                    !from 0 to p%GatesPerBeam-1
    REAL(DbKi),                                 INTENT(IN   )   ::  Time
    
    !Local variables
    INTEGER(IntKi)                                              ::  LoopCounter
    REAL(ReKi)                                                  ::  Dot_LidarPosition_I(3)
    
!FIXME: change to refer by name???
    Dot_LidarPosition_I(1) = y%IMUOutputs(11)
    Dot_LidarPosition_I(2) = y%IMUOutputs(14)
    Dot_LidarPosition_I(3) = y%IMUOutputs(17)
    
    y%AllOutputs( 1 : 18 )  = y%IMUOutputs ( 1 : 18 )
    if (size(p%MeasuringPoints_Spherical_L,2) > 0) then
       y%AllOutputs( 19 )      = p%MeasuringPoints_Spherical_L(2,m%LastMeasuringPoint+LoopGatesPerBeam)        !Azimuth
       y%AllOutputs( 20 )      = p%MeasuringPoints_Spherical_L(3,m%LastMeasuringPoint+LoopGatesPerBeam)        !Elevation
    endif
    y%AllOutputs( 21 )      = Time                                                                          !MeasTime
    y%AllOutputs( 22 )      = 1                                                                             !NewData
    y%AllOutputs( 23 )      = REAL(m%NextBeamID)                                                            !BeamID

    if (size(p%MeasuringPoints_Spherical_L,2) > 0) then
       y%AllOutputs( 24 + LoopGatesPerBeam ) = p%MeasuringPoints_Spherical_L(1,m%LastMeasuringPoint+LoopGatesPerBeam)   !Rangegates
    endif
    y%AllOutputs( 24 + p%GatesPerBeam + LoopGatesPerBeam ) = Vlos + (  DOT_PRODUCT(Dot_LidarPosition_I,UnitVector))  !Output the measured V_los. Consiting of the windspeed and the movement of the measuring system itself
    
    DO LoopCounter = 1,SIZE(p%ValidOutputs)
        y%WriteOutput( LoopCounter ) = y%AllOutputs( p%ValidOutputs(LoopCounter) )
    END DO
    
    y%SwapOutputs( 1                        )   = y%AllOutputs( 22 )                                        !NewData
    y%SwapOutputs( 2                        )   = y%AllOutputs( 23 )                                        !BeamID
    y%SwapOutputs( 3 + LoopGatesPerBeam     )   = y%AllOutputs( 24 + p%GatesPerBeam + LoopGatesPerBeam )    !Vlos

    y%SwapOutputs( 2 + p%GatesPerBeam   + 1 )   = y%IMUOutputs(  1 )    ! LdrRoll
    y%SwapOutputs( 2 + p%GatesPerBeam   + 2 )   = y%IMUOutputs(  4 )    ! LdrPitch
    y%SwapOutputs( 2 + p%GatesPerBeam   + 3 )   = y%IMUOutputs(  7 )    ! LdrYaw
    y%SwapOutputs( 2 + p%GatesPerBeam   + 4 )   = y%IMUOutputs( 11 )    ! LdrXd
    y%SwapOutputs( 2 + p%GatesPerBeam   + 5 )   = y%IMUOutputs( 14 )    ! LdrYd
    y%SwapOutputs( 2 + p%GatesPerBeam   + 6 )   = y%IMUOutputs( 17 )    ! LdrZd
    
END SUBROUTINE
    
    
!#########################################################################################################################################################################  
    
END MODULE LidarSim_Subs

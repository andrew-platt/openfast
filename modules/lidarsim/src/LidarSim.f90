MODULE LidarSim

   USE LidarSim_Types
   USE LidarSim_Subs
   USE NWTC_Library

   IMPLICIT NONE
   PRIVATE

   type(ProgDesc), parameter   :: LidarSim_ProgDesc = ProgDesc( 'LidarSim', 'v1.0', '6 Oct 2021' )

   PUBLIC                      ::  LidarSim_Init
   PUBLIC                      ::  LidarSim_CalcOutput
   PUBLIC                      ::  LidarSim_End

      ! These routines satisfy the framework, but do nothing at present.
   PUBLIC :: LidarSim_UpdateStates               !< Loose coupling routine for solving for constraint states, integrating continuous states, and updating discrete states
   PUBLIC :: LidarSim_CalcConstrStateResidual    !< Tight coupling routine for returning the constraint state residual
   PUBLIC :: LidarSim_CalcContStateDeriv         !< Tight coupling routine for computing derivatives of continuous states
   PUBLIC :: LidarSim_UpdateDiscState            !< Tight coupling routine for updating discrete states


   CONTAINS

   !#########################################################################################################################################################################

SUBROUTINE LidarSim_Init(InitInp, u, p, x, xd, z, OtherState, y, m, Interval, InitOutData, ErrStat, ErrMsg )

    IMPLICIT                                NONE
    CHARACTER(*),                           PARAMETER       ::  RoutineName="LidarSim_Init"
    
    TYPE(LidarSim_InitInputType),           INTENT(IN   )   ::  InitInp             ! Input data for initialization routine
    TYPE(LidarSim_InputType),               INTENT(  OUT)   ::  u                   !< An initial guess for the input; input mesh must be defined
    TYPE(LidarSim_ParameterType),           INTENT(  OUT)   ::  p                   !< Parameters
    TYPE(LidarSim_ContinuousStateType),     INTENT(  OUT)   ::  x                   !< Initial continuous states
    TYPE(LidarSim_DiscreteStateType),       INTENT(  OUT)   ::  xd                  !< Initial discrete states
    TYPE(LidarSim_ConstraintStateType),     INTENT(  OUT)   ::  z                   !< Initial guess of the constraint states
    TYPE(LidarSim_OtherStateType),          INTENT(  OUT)   ::  OtherState          !< Initial other states
    TYPE(LidarSim_OutputType),              INTENT(  OUT)   ::  y                   !< Initial system outputs (outputs are not calculated;
    TYPE(LidarSim_MiscVarType),             INTENT(  OUT)   ::  m                   !< MiscVars
    REAL(DbKi),                             INTENT(INOUT)   ::  Interval            !< timestep OpenFAST is using
    TYPE(LidarSim_InitOutputType),          INTENT(  OUT)   ::  InitOutData         !< Data to initialize the outputs, and pass to other modules
    INTEGER(IntKi),                         INTENT(  OUT)   ::  ErrStat             !< Error status of the operation
    CHARACTER(*),                           INTENT(  OUT)   ::  ErrMsg              !< Error message if ErrStat /= ErrID_None

    !Local Variables
    type(FileInfoType)                                      :: FileInfo_In          !< The derived type for holding the full input file for parsing -- we may pass this in the future
    TYPE(LidarSim_InputFile)                                :: InputFileData        !< Structure to load the input file data into
    CHARACTER(1024)                                         :: RootFileName
    character(1024)                                         :: PriPath              !< Primary path
    integer(IntKi)                                          :: i                    !< Generic counter
    integer(IntKi)                                          :: UnEcho               ! Unit number for the echo file

    ! Temporary variables for error handling
    INTEGER(IntKi)                                          ::  TmpErrStat          !< temporary error message
    CHARACTER(ErrMsgLen)                                    ::  TmpErrMsg           
    
    
    ! Initial error values
    ErrStat        =  0
    ErrMsg         =  ""  
    
    ! Initialize the NWTC Subroutine Library
    CALL NWTC_Init( )

    ! Display the module information
    CALL DispNVD( LidarSim_ProgDesc )
    InitOutData%Ver = LidarSim_ProgDesc

    ! Filename stuff
    RootFileName  = InitInp%RootName
    IF (LEN_TRIM(RootFileName) == 0) CALL GetRoot( InitInp%InputInitFile, RootFileName )

    p%RootName  = TRIM(InitInp%RootName)
    CALL GetPath( InitInp%InputInitFile, PriPath )     ! Input files will be relative to the path where the primary input file is located.




      ! -----------------------------------------------------------------
      ! Read the primary AeroDyn input file, or copy from passed input
   if (InitInp%UsePrimaryInputFile) then
      ! Read the entire input file, minus any comment lines, into the FileInfo_In
      ! data structure in memory for further processing.
      call ProcessComFile( InitInp%InputInitFile, FileInfo_In, TmpErrStat, TmpErrMsg )
   else
      call NWTC_Library_CopyFileInfoType( InitInp%PassedFileData, FileInfo_In, MESH_NEWCOPY, TmpErrStat, TmpErrMsg )
   endif
   if (Failed()) return;


   ! For diagnostic purposes, the following can be used to display the contents
   ! of the FileInfo_In data structure.
   ! call Print_FileInfo_Struct( CU, FileInfo_In ) ! CU is the screen -- different number on different systems.

      !  Parse the FileInfo_In structure of data from the inputfile into the InitInp%InputFile structure
   CALL LidarSim_ParsePrimaryFileInfo( PriPath, InitInp%InputInitFile, p%RootName, FileInfo_In, InputFileData, UnEcho, TmpErrStat, TmpErrMsg )
      if (Failed()) return;
!FIXME:  some sanity checks on the input file would be useful


   ! Convert angles read in degrees to radians
   InputFileData%RollAngle  =  InputFileData%RollAngle  *  D2R_D
   InputFileData%PitchAngle =  InputFileData%PitchAngle *  D2R_D
   InputFileData%YawAngle   =  InputFileData%YawAngle   *  D2R_D
   do i=1,InputFileData%NumberOfPoints_Spherical
      InputFileData%Azimuth(i)   = InputFileData%Azimuth(i)   * D2R_D
      InputFileData%Elevation(i) = InputFileData%Elevation(i) * D2R_D
   enddo

   !Transfering InputFileData to the p
   p%MeasurementMaxSteps   =   CEILING(REAL(NINT(InputFileData%t_measurement_interval*100000))/REAL(NINT(InitInp%DT*100000))) !NINT to remove float precision errors. Back to REAL, otherwise the divion ignores everything behind the decima point. Ceiling to round up to next integer
   p%GatesPerBeam          =   InputFileData%GatesPerBeam
!FIXME: can we add some error checking on this to make sure it is a sane value?
   p%MAXDLLChainOutputs    =   InputFileData%MAXDLLChainOutputs
   ! Are we asking IfW to override the uniform wind and calculate x-offset (violation of Uniform Wind methods -- not very physical)?
   if (InputFileData%URef > -9999.0_ReKi) then    ! default value set at parse is -9999.0
      p%UnifWind_Xoffset_URef    = InputFileData%URef
      p%UnifWind_Xoffset_URefF   = .true.
   else
      p%UnifWind_Xoffset_URef    = 0.0_ReKi
      p%UnifWind_Xoffset_URefF   = .false.
   endif
   InitOutData%UnifWind_Xoffset_URef   = p%UnifWind_Xoffset_URef
   InitOutData%UnifWind_Xoffset_URefF  = p%UnifWind_Xoffset_URefF

    ! Create the mesh for the u%LidarMotion mesh
   call Init_LidarMountMesh( InitInp, InputFileData, y, p, TmpErrStat, TmpErrMsg )
   if (Failed())  return


   !----- Calls different Subroutines to initialize the measuring points
   IF(InputFileData%TrajectoryType == 0) THEN
      CALL LidarSim_InitMeasuringPoints_Cartesian(p, InputFileData, TmpErrStat, TmpErrMsg)    ! Calls Routine to initialize cartesian coordinate inputs
      if (Failed())  return
   ELSEIF(InputFileData%TrajectoryType == 1)THEN
       CALL LidarSim_InitMeasuringPoints_Spherical(p, InputFileData, TmpErrStat, TmpErrMsg )   ! Calls Routine to initialize spherical coordinate inputs
      if (Failed())  return
   END IF

   !----- Calls different Subroutines to initialize the weighting points
   IF(InputFileData%WeightingType == 0) THEN                                                   ! Single point
       CALL AllocAry( p%WeightingDistance,1, 'p%WeightingDistance', TmpErrStat, TmpErrMsg )
      if (Failed())  return
       CALL AllocAry( p%Weighting,1, 'p%Weighting', TmpErrStat, TmpErrMsg )
      if (Failed())  return
       p%WeightingDistance(1) = 0
       p%Weighting(1) = 1
       p%nWeightPts = 1
   ELSEIF(InputFileData%WeightingType == 1) THEN                                               ! Calls Routine to initialize the weighting with gaussian distribution
       CALL LidarSim_InitializeWeightingGauss(p, InputFileData, TmpErrStat, TmpErrMsg )
       if (Failed())  return
       p%nWeightPts = size(p%WeightingDistance)
   ELSEIF(InputFileData%WeightingType == 2) THEN
       CALL LidarSim_InitializeWeightingManual(p, InputFileData, TmpErrStat, TmpErrMsg )       ! Calls Routine to initialize with manual weighting settings
       if (Failed())  return
       p%nWeightPts = size(p%WeightingDistance)
   ENDIF

   ! create mesh for the y%LidarMeasPos
   !     This mesh is used for visualization and for retrieving wind velocity information
   call Init_LidarMeasMesh( u, p, y, TmpErrStat, TmpErrMsg )
   if (Failed())  return

   ! Initialize all outputs
   CALL LidarSim_InitializeOutputs(y,p, InitOutData, InputFileData, TmpErrStat, TmpErrMsg )
   if (Failed())  return

   !initialize variables and outputs
   m%MeasurementCurrentStep = -1                                                               !< there was no measurement yet
   m%LastMeasuringPoint = 1                                                                    !< First measurement point
   m%NextBeamID = 0

   call Cleanup()
   return
 
   contains
   logical function Failed()
      CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
      Failed = ErrStat >= AbortErrLev
      if (Failed) then
         call Cleanup
      endif
   end function Failed

   subroutine Cleanup()
      CALL LidarSim_DestroyInputFile(InputFileData, TmpErrStat, TmpErrMsg)
      CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
   end subroutine Cleanup

   subroutine Init_LidarMountMesh( InitInp, InputFileData, y, p, ErrStat2, ErrMsg2 )
      TYPE(LidarSim_InitInputType),          INTENT(IN   )  :: InitInp           ! Input data for initialization routine
      TYPE(LidarSim_InputFile),              INTENT(IN   )  :: InputFileData     !< Structure to load the input file data into
      TYPE(LidarSim_ParameterType),          INTENT(INOUT)  :: p                 !< Parameters
      TYPE(LidarSim_OutputType),             INTENT(INOUT)  :: y                 !< Initial system outputs (outputs are not calculated;
      INTEGER(IntKi),                        INTENT(  OUT)  :: ErrStat2           !< Error status of the operation
      CHARACTER(ErrMsgLen),                  INTENT(  OUT)  :: ErrMsg2            !< Error message if ErrStat /= ErrID_None

      INTEGER(IntKi)                                        :: TmpErrStat
      CHARACTER(ErrMsgLen)                                  :: TmpErrMsg
      real(ReKi)                                            :: Pos(3)            ! Position    of the lidar unit (global coords)
      real(R8Ki)                                            :: Orient(3,3)       ! Orientation of the lidar unit (global coords)
      real(R8Ki)                                            :: theta(3)          ! Euler angles input
      ! Initial error values
      ErrStat     =  ErrID_None
      ErrMsg      = ""


      !  Creates the static rotationmatrix from the lidar system to the reference system
      !  Note: reference system could be the nacelle, hub, ground, or floating platform.  Depends what point is passed in
      !        for where the Lidar is mounted.  The calling code decides this.

      !---------------------------------------
      ! Position of the lidar module in global coordinates
      Pos      =  InitInp%LidarRefPosition + (/ InputFileData%LidarPositionX, InputFileData%LidarPositionY, InputFileData%LidarPositionZ /)
      theta    = (/ InputFileData%RollAngle, InputFileData%PitchAngle, InputFileData%YawAngle /)
      Orient   = EulerConstruct(theta)
      Orient   =  MATMUL( transpose(Orient), InitInp%LidarRefOrientation )

      ! Create the input mesh for the Lidar unit
      CALL MeshCreate( BlankMesh        = u%LidarMesh       &
                     ,IOS               = COMPONENT_INPUT   &
                     ,Nnodes            = 1                 &
                     ,ErrStat           = TmpErrStat        &
                     ,ErrMess           = TmpErrMsg         &
                     ,TranslationDisp   = .TRUE.            &
                     ,Orientation       = .TRUE.            &
                     ,TranslationVel    = .TRUE.            &
                     ,RotationVel       = .TRUE.            &
                     ,TranslationAcc    = .TRUE.            &
                     ,RotationAcc       = .TRUE.)
         CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat2, ErrMsg2, '')
         if (ErrStat2 >= AbortErrLev) return

      ! Create the node on the mesh
      CALL MeshPositionNode ( u%LidarMesh,1, Pos, TmpErrStat, TmpErrMsg, Orient )
         CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat2, ErrMsg2, '')

      ! Create the mesh element
      CALL MeshConstructElement (  u%LidarMesh         &
                                  , ELEMENT_POINT      &
                                  , TmpErrStat         &
                                  , TmpErrMsg          &
                                  , 1                  )
         CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat2, ErrMsg2, '')
      CALL MeshCommit ( u%LidarMesh         &
                      , TmpErrStat          &
                      , TmpErrMsg           )
         CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat2, ErrMsg2, '')

      u%LidarMesh%TranslationDisp = 0.0_R8Ki
      u%LidarMesh%TranslationVel  = 0.0_ReKi
      u%LidarMesh%TranslationAcc  = 0.0_ReKi
      u%LidarMesh%Orientation     = u%LidarMesh%RefOrientation
      u%LidarMesh%RotationVel     = 0.0_ReKi
      u%LidarMesh%RotationAcc     = 0.0_ReKi

      u%LidarMesh%RemapFlag  = .false.
   end subroutine Init_LidarMountMesh

   !> This creates a mesh with all the measurement points.  This is used in finding the wind velocity measurements
   !! to pass in, and for visualization of the measurement locations.
   subroutine Init_LidarMeasMesh( u, p, y, ErrStat2, ErrMsg2 )
      type(LidarSim_InputType),              intent(inout)  :: u                 !< input
      type(LidarSim_ParameterType),          intent(in   )  :: p                 !< Parameters
      type(LidarSim_OutputType),             intent(inout)  :: y                 !< Initial system outputs (outputs are not calculated;
      integer(IntKi),                        intent(  out)  :: ErrStat2          !< Error status of the operation
      character(ErrMsgLen),                  intent(  out)  :: ErrMsg2           !< Error message if ErrStat /= ErrID_None

      integer(IntKi)                                        :: i,j               !< generic counter
      integer(IntKi)                                        :: iPt               !< index to mesh
      integer(IntKi)                                        :: nTot              !< total number of measurement points
      REAL(ReKi)                                            :: MeasPos_I(3)      !< Measuring Position of given weighted point in global coords
      REAL(ReKi)                                            :: MeasPosCtr_I(3)   !< Central point in the weighted measurement points
      REAL(ReKi)                                            :: LidarPos_I(3)     !< Lidar Mount point
      REAL(ReKi)                                            :: UnitVector(3)     !< Vector from Lidar mount to measurement position
      integer(IntKi)                                        :: TmpErrStat
      character(ErrMsgLen)                                  :: TmpErrMsg
      ! Initial error values
      ErrStat     =  ErrID_None
      ErrMsg      = ""

      ! Total number of measurement points
      nTot = SIZE(p%MeasuringPoints_L,2)*p%nWeightPts

      ! Create the output mesh for the Lidar measurement locations
      CALL MeshCreate( BlankMesh        = y%LidarMeasMesh   &
                     ,IOS               = COMPONENT_OUTPUT  &
                     ,Nnodes            = nTot              &     ! total number of measurement points
                     ,ErrStat           = TmpErrStat        &
                     ,ErrMess           = TmpErrMsg         &
                     ,TranslationDisp   = .TRUE.            &
                     ,Orientation       = .false.           &
                     ,TranslationVel    = .false.           &
                     ,RotationVel       = .false.           &
                     ,TranslationAcc    = .false.           &
                     ,RotationAcc       = .false.           )
         CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat2, ErrMsg2, '')
         if (ErrStat2 >= AbortErrLev) return

      ! Lidar position
      LidarPos_I =  u%LidarMesh%Position(1:3,1) + u%LidarMesh%TranslationDisp(1:3,1)

      ! Loop over all beams
      do i=1,size(p%MeasuringPoints_L,2)

         ! Get central position in inertial coordinates:
         MeasPosCtr_I = LidarSim_TransformLidarToInertial(u%LidarMesh, p, p%MeasuringPoints_L(:,i))

         ! Line of Sight unit vector
         UnitVector    =   MeasPosCtr_I - LidarPos_I          ! Calculation of the Line of Sight Vector
         UnitVector    =   UnitVector/NORM2(UnitVector)       ! =>Magnitude = 1

         ! loop over weighted measure points along beam about central point
         do j=1,p%nWeightPts
            iPt = (i-1)*p%nWeightPts+j

            MeasPos_I = MeasPosCtr_I + p%WeightingDistance(j) * UnitVector

            ! Create the node on the mesh
            CALL MeshPositionNode ( y%LidarMeasMesh,iPt, MeasPos_I, TmpErrStat, TmpErrMsg )
               CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat2, ErrMsg2, '')

            ! Create the mesh element
            CALL MeshConstructElement (  y%LidarMeasMesh     &
                                        , ELEMENT_POINT      &
                                        , TmpErrStat         &
                                        , TmpErrMsg          &
                                        , iPt                )
               CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat2, ErrMsg2, '')
         enddo
      enddo

      CALL MeshCommit ( y%LidarMeasMesh     &
                      , TmpErrStat          &
                      , TmpErrMsg           )
         CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat2, ErrMsg2, '')

      y%LidarMeasMesh%TranslationDisp = 0.0_R8Ki
      y%LidarMeasMesh%RemapFlag  = .false.

      !--------------------
      ! create mesh mapping
      call MeshMapCreate( u%LidarMesh, y%LidarMeasMesh, m%u_L_p2p_y_Lmeas, TmpErrStat, TmpErrMsg )
         call SetErrStat( TmpErrStat, TmpErrMsg, ErrStat2, ErrMsg2, '')
         if (ErrStat2 >= AbortErrLev) return
      call AllocAry( u%WindVelocity, 3, y%LidarMeasMesh%Nnodes, 'WindVelocity', TmpErrStat, TmpErrMsg )
         call SetErrStat( TmpErrStat, TmpErrMsg, ErrStat2, ErrMsg2, '')
         if (ErrStat2 >= AbortErrLev) return

   end subroutine Init_LidarMeasMesh

END SUBROUTINE LidarSim_Init

!#########################################################################################################################################################################

SUBROUTINE LidarSim_CalcOutput (Time, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg )
    IMPLICIT                                    NONE
    CHARACTER(*),                               PARAMETER           ::  RoutineName="LidarSim_CalcOutput"

    REAL(DbKi),                                 INTENT(IN   )       ::  Time                !< Current simulation time in seconds
    TYPE(LidarSim_InputType),                   INTENT(IN   )       ::  u                   !< An initial guess for the input; input mesh must be defined
    TYPE(LidarSim_ParameterType),               INTENT(IN   )       ::  p                   !< Parameters
    TYPE(LidarSim_ContinuousStateType),         INTENT(IN   )       ::  x                   !< Initial continuous states
    TYPE(LidarSim_DiscreteStateType),           INTENT(IN   )       ::  xd                  !< Initial discrete states
    TYPE(LidarSim_ConstraintStateType),         INTENT(IN   )       ::  z                   !< Initial guess of the constraint states
    TYPE(LidarSim_OtherStateType),              INTENT(IN   )       ::  OtherState          !< Initial other states
    TYPE(LidarSim_OutputType),                  INTENT(INOUT)       ::  y                   !< Initial system outputs (outputs are not calculated;
    TYPE(LidarSim_MiscVarType),                 INTENT(INOUT)       ::  m                   !< MiscVars
    INTEGER(IntKi),                             INTENT(  OUT)       ::  ErrStat                     !< Error status of the operation
    CHARACTER(*),                               INTENT(  OUT)       ::  ErrMsg                      !< Error message if ErrStat /= ErrID_None

    !Local Variables
    REAL(ReKi)                                                      ::  UnitVector(3)               !Line of Sight Unit Vector
    REAL(ReKi)                                                      ::  MeasuringPosition_I(3)      !Transformed Measuring Position
    REAL(ReKi)                                                      ::  MeasuringPositionMesh_I(3)  !Transformed Measuring Position
    REAL(ReKi)                                                      ::  LidarPosition_I(3)          !Transformed Lidar Position
    REAL(ReKi)                                                      ::  Vlos                        !Line of sight speed
    REAL(ReKi)                                                      ::  WeightPos(3,p%nWeightPts)   ! Inertial position of weighted measure point
    REAL(ReKi)                                                      ::  WindVel(3,p%nWeightPts)     ! Velocity of wind at weighted measure point
    INTEGER(IntKi)                                                  ::  LoopGatesPerBeam            ! Counter to loop through all gate points of a line
    INTEGER(IntKi)                                                  ::  LoopCounter
    integer(IntKi)                                                  ::  StrtPt                      ! Starting index for set of weighted beam locations in mesh

    ! Temporary variables for error handling
    INTEGER(IntKi)                                                  ::  TmpErrStat
    CHARACTER(ErrMsgLen)                                            ::  TmpErrMsg

    !Initialize error values
    ErrStat        =  ErrID_None
    ErrMsg         =  ""

    ! update the output mesh
    call Transfer_Point_to_Point( u%LidarMesh, y%LidarMeasMesh, m%u_L_p2p_y_Lmeas, TmpErrStat, TmpErrMsg )
    call SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)

    IF(m%MeasurementCurrentStep>=p%MeasurementMaxSteps .OR. m%MeasurementCurrentStep == -1)THEN         !Check if there must be a new measurement     !(NINT(Time*1000)-NINT(p%t_last_measurement*1000)) >= NINT(p%t_measurement_interval*1000)
        m%MeasurementCurrentStep = 0
        LidarPosition_I =  u%LidarMesh%Position(1:3,1) + u%LidarMesh%TranslationDisp(1:3,1)

        CALL LidarSim_CalculateIMU(p, y, u)
        DO LoopGatesPerBeam = 0,p%GatesPerBeam-1
            StrtPt = p%nWeightPts * (m%LastMeasuringPoint-1 + LoopGatesPerBeam) + 1
            ! lidar measurement position in inertial coordinates
!FIXME: can use the first weighted point for the unit vector
            MeasuringPosition_I = LidarSim_TransformLidarToInertial(u%LidarMesh, p, p%MeasuringPoints_L(:,m%LastMeasuringPoint+LoopGatesPerBeam)) ! Calculate the Measuringpoint coordinate in the initial system
            WeightPos(1:3,1:p%nWeightPts) = y%LidarMeasMesh%Position(       :,StrtPt:StrtPt+p%nWeightPts-1) &
                                          + y%LidarMeasMesh%TranslationDisp(:,StrtPt:StrtPt+p%nWeightPts-1)
            WindVel(1:3,1:p%nWeightPts)   = u%WindVelocity(                 :,StrtPt:StrtPt+p%nWeightPts-1)

            !Line of Sight
            UnitVector    =   MeasuringPosition_I - LidarPosition_I             !Calculation of the Line of Sight Vector
            UnitVector    =   UnitVector/NORM2(UnitVector)                      !=>Magnitude = 1

            ! Calculation of the wind speed at the calculated position
            CALL LidarSim_CalculateVlos( p, UnitVector, Vlos, WeightPos, WindVel, TmpErrStat, TmpErrMsg) !Calculation of the line of sight wind speed
            CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
            CALL LidarSim_SetOutputs(y,p,m,Vlos,UnitVector,LoopGatesPerBeam,Time)    !Set all outputs to the output variable
        ENDDO

        !Choosing which measuring point has to be calculated
        IF(m%LastMeasuringPoint+p%GatesPerBeam > SIZE(p%MeasuringPoints_L,2)) THEN
            m%LastMeasuringPoint = 1        ! already reached the last point before ? => start over from the beginning
            m%NextBeamID = 0
        ELSE
            m%LastMeasuringPoint = m%LastMeasuringPoint + p%GatesPerBeam
            m%NextBeamID = m%NextBeamID + 1
        END IF
    ELSE                                            !Set NewData signals to zero
        IF(ANY(p%ValidOutputs == 22)) THEN
            DO LoopCounter = 1,SIZE(p%ValidOutputs)
                IF(p%ValidOutputs(LoopCounter) == 22) THEN
                    y%WriteOutput(LoopCounter) = 0
                END IF
            END DO
        END IF
        y%SwapOutputs(1) = 0
    ENDIF
    m%MeasurementCurrentStep = m%MeasurementCurrentStep + 1

END SUBROUTINE LidarSim_CalcOutput

!#########################################################################################################################################################################

SUBROUTINE LidarSim_End( InputData, p, ContStates, DiscStates, ConstrStateGuess, OtherStates, y, m, ErrStat, ErrMsg )

   IMPLICIT NONE

   CHARACTER(*),              PARAMETER                     :: RoutineName="LidarSim_End"

   TYPE(LidarSim_InputType),               INTENT(INOUT)  :: InputData         !< Input data for initialization
   TYPE(LidarSim_ParameterType),           INTENT(INOUT)  :: p         !< Parameters
   TYPE(LidarSim_ContinuousStateType),     INTENT(INOUT)  :: ContStates        !< Continuous states
   TYPE(LidarSim_DiscreteStateType),       INTENT(INOUT)  :: DiscStates        !< Discrete states
   TYPE(LidarSim_ConstraintStateType),     INTENT(INOUT)  :: ConstrStateGuess  !< Guess of the constraint states
   TYPE(LidarSim_OtherStateType),          INTENT(INOUT)  :: OtherStates       !< Other/optimization states
   TYPE(LidarSim_OutputType),              INTENT(INOUT)  :: y           !< Output data
   TYPE(LidarSim_MiscVarType),             INTENT(INOUT)  :: m          !< Misc variables for optimization (not copied in glue code)


   ! Error Handling
   INTEGER( IntKi ),                        INTENT(  OUT)   :: ErrStat      !< error status
   CHARACTER(*),                            INTENT(  OUT)   :: ErrMsg       !< error message


   ErrStat = ErrID_None
   ErrMsg = ""

   CALL LidarSim_DestroyInput( InputData, ErrStat, ErrMsg )
   CALL LidarSim_DestroyParam( p, ErrStat, ErrMsg )
   CALL LidarSim_DestroyContState( ContStates, ErrStat, ErrMsg )
   CALL LidarSim_DestroyDiscState( DiscStates, ErrStat, ErrMsg )
   CALL LidarSim_DestroyConstrState( ConstrStateGuess, ErrStat, ErrMsg )
   CALL LidarSim_DestroyOtherState( OtherStates, ErrStat, ErrMsg )
   CALL LidarSim_DestroyOutput( y, ErrStat, ErrMsg )
   CALL LidarSim_DestroyMisc( m, ErrStat, ErrMsg )

END SUBROUTINE LidarSim_End

    !#########################################################################################################################################################################





!====================================================================================================
! The following routines were added to satisfy the framework, but do nothing useful.
!====================================================================================================
!> This is a loose coupling routine for solving constraint states, integrating continuous states, and updating discrete and other
!! states. Continuous, constraint, discrete, and other states are updated to values at t + Interval.
SUBROUTINE LidarSim_UpdateStates( t, n, Inputs, InputTimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )

   REAL(DbKi),                            INTENT(IN   ) :: t               !< Current simulation time in seconds
   INTEGER(IntKi),                        INTENT(IN   ) :: n               !< Current step of the simulation: t = n*Interval
   TYPE(LidarSim_InputType),              INTENT(INOUT) :: Inputs(:)       !< Inputs at InputTimes (output only for mesh record-keeping in ExtrapInterp routine)
   REAL(DbKi),                            INTENT(IN   ) :: InputTimes(:)   !< Times in seconds associated with Inputs
   TYPE(LidarSim_ParameterType),          INTENT(IN   ) :: p               !< Parameters
   TYPE(LidarSim_ContinuousStateType),    INTENT(INOUT) :: x               !< Input: Continuous states at t;
                                                                             !!    Output: Continuous states at t + Interval
   TYPE(LidarSim_DiscreteStateType),      INTENT(INOUT) :: xd              !< Input: Discrete states at t;
                                                                             !!    Output: Discrete states at t  + Interval
   TYPE(LidarSim_ConstraintStateType),    INTENT(INOUT) :: z               !< Input: Constraint states at t;
                                                                             !!   Output: Constraint states at t + Interval
   TYPE(LidarSim_OtherStateType),         INTENT(INOUT) :: OtherState      !< Other states: Other states at t;
                                                                             !!   Output: Other states at t + Interval
   TYPE(LidarSim_MiscVarType),            INTENT(INOUT) :: m               !< Misc variables for optimization (not copied in glue code)
   INTEGER(IntKi),                        INTENT(  OUT) :: ErrStat         !< Error status of the operation
   CHARACTER(*),                          INTENT(  OUT) :: ErrMsg          !< Error message if ErrStat /= ErrID_None


      ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ""

   x%DummyContState     = 0.0_ReKi
   xd%DummyDiscState    = 0.0_ReKi
   z%DummyConstrState   = 0.0_ReKi

   RETURN


END SUBROUTINE LidarSim_UpdateStates

!----------------------------------------------------------------------------------------------------------------------------------
!> Tight coupling routine for computing derivatives of continuous states
SUBROUTINE LidarSim_CalcContStateDeriv( Time, u, p, x, xd, z, OtherState, m, dxdt, ErrStat, ErrMsg )
!..................................................................................................................................

   REAL(DbKi),                               INTENT(IN   )  :: Time        !< Current simulation time in seconds
   TYPE(LidarSim_InputType),                 INTENT(IN   )  :: u           !< Inputs at Time
   TYPE(LidarSim_ParameterType),             INTENT(IN   )  :: p           !< Parameters
   TYPE(LidarSim_ContinuousStateType),       INTENT(IN   )  :: x           !< Continuous states at Time
   TYPE(LidarSim_DiscreteStateType),         INTENT(IN   )  :: xd          !< Discrete states at Time
   TYPE(LidarSim_ConstraintStateType),       INTENT(IN   )  :: z           !< Constraint states at Time
   TYPE(LidarSim_OtherStateType),            INTENT(IN   )  :: OtherState  !< Other states at Time
   TYPE(LidarSim_MiscVarType),               INTENT(INOUT)  :: m           !< Misc variables for optimization (not copied in glue code)
   TYPE(LidarSim_ContinuousStateType),       INTENT(  OUT)  :: dxdt        !< Continuous state derivatives at Time
   INTEGER(IntKi),                           INTENT(  OUT)  :: ErrStat     !< Error status of the operation
   CHARACTER(*),                             INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None


      ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ""


      ! Compute the first time derivatives of the continuous states here:

   dxdt%DummyContState = 0.0_ReKi


END SUBROUTINE LidarSim_CalcContStateDeriv

!----------------------------------------------------------------------------------------------------------------------------------
!> Tight coupling routine for updating discrete states
SUBROUTINE LidarSim_UpdateDiscState( Time, u, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )

   REAL(DbKi),                               INTENT(IN   )  :: Time        !< Current simulation time in seconds
   TYPE(LidarSim_InputType),                 INTENT(IN   )  :: u           !< Inputs at Time
   TYPE(LidarSim_ParameterType),             INTENT(IN   )  :: p           !< Parameters
   TYPE(LidarSim_ContinuousStateType),       INTENT(IN   )  :: x           !< Continuous states at Time
   TYPE(LidarSim_DiscreteStateType),         INTENT(INOUT)  :: xd          !< Input: Discrete states at Time;
                                                                             !! Output: Discrete states at Time + Interval
   TYPE(LidarSim_ConstraintStateType),       INTENT(IN   )  :: z           !< Constraint states at Time
   TYPE(LidarSim_OtherStateType),            INTENT(IN   )  :: OtherState  !< Other states at Time
   TYPE(LidarSim_MiscVarType),               INTENT(INOUT)  :: m           !< Misc variables for optimization (not copied in glue code)
   INTEGER(IntKi),                           INTENT(  OUT)  :: ErrStat     !< Error status of the operation
   CHARACTER(*),                             INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None


      ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ""


      ! Update discrete states here:

   ! StateData%DiscState =

END SUBROUTINE LidarSim_UpdateDiscState

!----------------------------------------------------------------------------------------------------------------------------------
!> Tight coupling routine for solving for the residual of the constraint state equations
SUBROUTINE LidarSim_CalcConstrStateResidual( Time, u, p, x, xd, z, OtherState, m, z_residual, ErrStat, ErrMsg )

   REAL(DbKi),                               INTENT(IN   )  :: Time        !< Current simulation time in seconds
   TYPE(LidarSim_InputType),                 INTENT(IN   )  :: u           !< Inputs at Time
   TYPE(LidarSim_ParameterType),             INTENT(IN   )  :: p           !< Parameters
   TYPE(LidarSim_ContinuousStateType),       INTENT(IN   )  :: x           !< Continuous states at Time
   TYPE(LidarSim_DiscreteStateType),         INTENT(IN   )  :: xd          !< Discrete states at Time
   TYPE(LidarSim_ConstraintStateType),       INTENT(IN   )  :: z           !< Constraint states at Time (possibly a guess)
   TYPE(LidarSim_OtherStateType),            INTENT(IN   )  :: OtherState  !< Other states at Time
   TYPE(LidarSim_MiscVarType),               INTENT(INOUT)  :: m           !< Misc variables for optimization (not copied in glue code)
   TYPE(LidarSim_ConstraintStateType),       INTENT(  OUT)  :: z_residual  !< Residual of the constraint state equations using
                                                                           !! the input values described above
   INTEGER(IntKi),                           INTENT(  OUT)  :: ErrStat     !< Error status of the operation
   CHARACTER(*),                             INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None


      ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ""


      ! Solve for the constraint states here:

   z_residual%DummyConstrState = 0

END SUBROUTINE LidarSim_CalcConstrStateResidual


END MODULE LidarSim

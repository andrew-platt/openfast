MODULE LidarSim

   USE LidarSim_Types
   USE LidarSim_Subs
   USE NWTC_Library
!FIXME: remove IfW completely from this module.
   USE InflowWind
   USE InflowWind_Subs
   USE InflowWind_Types

   IMPLICIT NONE
   PRIVATE

   TYPE(ProgDesc), PARAMETER   ::  IfW_Ver = ProgDesc( 'LidarSim', 'v0.20', '12-December-2019' )
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
    TYPE(LidarSim_InitOutputType),          INTENT(  OUT)   ::  InitOutData         !< Data to initialize the outputs
    INTEGER(IntKi),                         INTENT(  OUT)   ::  ErrStat             !< Error status of the operation
    CHARACTER(*),                           INTENT(  OUT)   ::  ErrMsg              !< Error message if ErrStat /= ErrID_None

    !Local Variables
    TYPE(LidarSim_InputFile)                                ::  InputFileData      !< Structure to load the input file data into
    CHARACTER(1024)                                         ::  RootFileName
    CHARACTER(1024)                                         ::  EchoFileName
   integer(IntKi)                                           :: i                    !< Generic counter

    ! Temporary variables for error handling
    INTEGER(IntKi)                                          ::  TmpErrStat          !< temporary error message
    CHARACTER(ErrMsgLen)                                    ::  TmpErrMsg           
    
    
    ! Initial error values
    ErrStat        =  0
    ErrMsg         =  ""  
    
    RootFileName  = InitInp%RootName
    IF (LEN_TRIM(RootFileName) == 0) CALL GetRoot( InitInp%InputInitFile, RootFileName )
    EchoFileName  = TRIM(RootFileName)//".ech"

    
    ! Reads the Config from the Input file and writes it into the LidarSim_InputFile data 
    CALL LidarSim_ReadInputFile(InitInp%InputInitFile,EchoFileName , InputFileData, TmpErrStat, TmpErrMsg)!EchoFileName
      if (Failed()) return

   ! Convert angles read in degrees to radians
   InputFileData%RollAngle_N  =  InputFileData%RollAngle_N  *  D2R_D
   InputFileData%PitchAngle_N =  InputFileData%PitchAngle_N *  D2R_D
   InputFileData%YawAngle_N   =  InputFileData%YawAngle_N   *  D2R_D
   do i=1,InputFileData%NumberOfPoints_Spherical
      InputFileData%Azimuth(i)   = InputFileData%Azimuth(i)   * D2R_D
      InputFileData%Elevation(i) = InputFileData%Elevation(i) * D2R_D
   enddo
 
    !Transfering InputFileData to the p
    p%MeasurementMaxSteps   =   CEILING(REAL(NINT(InputFileData%t_measurement_interval*100000))/REAL(NINT(InitInp%DT*100000))) !NINT to remove float precision errors. Back to REAL, otherwise the divion ignores everything behind the decima point. Ceiling to round up to next integer
    p%LidarPosition_N(1)    =   InputFileData%LidarPositionX_N
    p%LidarPosition_N(2)    =   InputFileData%LidarPositionY_N
    p%LidarPosition_N(3)    =   InputFileData%LidarPositionZ_N
    p%URef                  =   InputFileData%URef
    p%GatesPerBeam          =   InputFileData%GatesPerBeam
!FIXME: can we add some error checking on this to make sure it is a sane value?
    p%MAXDLLChainOutputs    =   InputFileData%MAXDLLChainOutputs

    ! Create the mesh for the u%LidarMotion mesh
   call Init_LidarMountMesh( InitInp, InputFileData, y, p, TmpErrStat, TmpErrMsg )


    !----- Calls different Subroutines to initialize the measuring points   
    IF(InputFileData%TrajectoryType == 0) THEN
        CALL LidarSim_InitMeasuringPoints_Cartesian(p, InputFileData, TmpErrStat, TmpErrMsg)    ! Calls Routine to initialize cartesian coordinate inputs
        CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
    ELSEIF(InputFileData%TrajectoryType == 1)THEN
        CALL LidarSim_InitMeasuringPoints_Spherical(p, InputFileData, TmpErrStat, TmpErrMsg )   ! Calls Routine to initialize spherical coordinate inputs
        CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
    END IF
    
    !----- Calls different Subroutines to initialize the weighting points
    IF(InputFileData%WeightingType == 0) THEN                                                   ! Single point
        CALL AllocAry( p%WeightingDistance,1, 'p%WeightingDistance', TmpErrStat, TmpErrMsg )
        CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
        CALL AllocAry( p%Weighting,1, 'p%Weighting', TmpErrStat, TmpErrMsg )
        CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
        p%WeightingDistance(1) = 0
        p%Weighting(1) = 1
    ELSEIF(InputFileData%WeightingType == 1) THEN                                               ! Calls Routine to initialize the weighting with gaussian distribution
        CALL LidarSim_InitializeWeightingGauss(p, InputFileData, TmpErrStat, TmpErrMsg )       
        CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
    ELSEIF(InputFileData%WeightingType == 2) THEN
        CALL LidarSim_InitializeWeightingManual(p, InputFileData, TmpErrStat, TmpErrMsg )       ! Calls Routine to initialize with manual weighting settings
        CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
    ENDIF

    CALL LidarSim_InitializeOutputs(y,p, InitOutData, InputFileData, TmpErrStat, TmpErrMsg )
    CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
    
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

   subroutine Init_LidarMountMesh( InitInp, InputFileData, y, p, ErrStat, ErrMsg )
      TYPE(LidarSim_InitInputType),          INTENT(IN   )  :: InitInp           ! Input data for initialization routine
      TYPE(LidarSim_InputFile),              INTENT(IN   )  :: InputFileData     !< Structure to load the input file data into
      TYPE(LidarSim_ParameterType),          INTENT(INOUT)  :: p                 !< Parameters
      TYPE(LidarSim_OutputType),             INTENT(INOUT)  :: y                 !< Initial system outputs (outputs are not calculated;
      INTEGER(IntKi),                        INTENT(  OUT)  :: ErrStat           !< Error status of the operation
      CHARACTER(ErrMsgLen),                  INTENT(  OUT)  :: ErrMsg            !< Error message if ErrStat /= ErrID_None

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
      Pos      =  p%LidarPosition_N + InitInp%LidarRefPosition
      theta    = (/ InputFileData%RollAngle_N, InputFileData%PitchAngle_N, InputFileData%YawAngle_N /)
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
         CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)

      ! Create the node on the mesh
      CALL MeshPositionNode ( u%LidarMesh,1, Pos, TmpErrStat, TmpErrMsg, Orient )
         CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)

      ! Create the mesh element
      CALL MeshConstructElement (  u%LidarMesh         &
                                  , ELEMENT_POINT      &
                                  , TmpErrStat         &
                                  , TmpErrMsg          &
                                  , 1                  )
         CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
      CALL MeshCommit ( u%LidarMesh         &
                      , TmpErrStat          &
                      , TmpErrMsg           )
         CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)

      u%LidarMesh%RemapFlag  = .TRUE.
   end subroutine Init_LidarMountMesh

END SUBROUTINE LidarSim_Init

    !#########################################################################################################################################################################

SUBROUTINE LidarSim_CalcOutput (Time, u, p, x, xd, z, OtherState, y, m,&
    IfW_p, IfW_ContStates, IfW_DiscStates, IfW_ConstrStates, IfW_OtherStates,  IfW_m,&
    ErrStat, ErrMsg )
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

!FIXME: Remove IfW completely from here
    !Data for CalcOutput of IfW_Subs
    TYPE(InflowWind_ParameterType),             INTENT(IN   )       ::  IfW_p                       !< Parameters
    TYPE(InflowWind_ContinuousStateType),       INTENT(IN   )       ::  IfW_ContStates              !< Continuous states at Time
    TYPE(InflowWind_DiscreteStateType),         INTENT(IN   )       ::  IfW_DiscStates              !< Discrete states at Time
    TYPE(InflowWind_ConstraintStateType),       INTENT(IN   )       ::  IfW_ConstrStates            !< Constraint states at Time
    TYPE(InflowWind_OtherStateType),            INTENT(IN   )       ::  IfW_OtherStates             !< Other/optimization states at Time
    TYPE(InflowWind_MiscVarType),               INTENT(INOUT)       ::  IfW_m                       !< Misc variables for optimization (not copied in glue code)    

    !Local Variables
!FIXME: remove IfW completely from here
    TYPE(InflowWind_InputType)                                      ::  InputForCalculation         !Data Field needed for the calculation of the windspeed
    TYPE(InflowWind_OutputType)                                     ::  OutputForCalculation        !datafield in which the calculated speed is stored
    REAL(ReKi)                                                      ::  UnitVector(3)               !Line of Sight Unit Vector
    REAL(ReKi)                                                      ::  MeasuringPosition_I(3)      !Transformed Measuring Position
    REAL(ReKi)                                                      ::  LidarPosition_I(3)          !Transformed Lidar Position
    REAL(ReKi)                                                      ::  Vlos                        !Line of sight speed
    INTEGER(IntKi)                                                  ::  LoopGatesPerBeam            !Counter to loop through all gate points of a line
    INTEGER(IntKi)                                                  ::  LoopCounter

    ! Temporary variables for error handling
    INTEGER(IntKi)                                                  ::  TmpErrStat        
    CHARACTER(ErrMsgLen)                                            ::  TmpErrMsg
    
    !Initialize error values
    ErrStat        =  0
    ErrMsg         =  ""

    
    IF(m%MeasurementCurrentStep>=p%MeasurementMaxSteps .OR. m%MeasurementCurrentStep == -1)THEN         !Check if there must be a new measurement     !(NINT(Time*1000)-NINT(p%t_last_measurement*1000)) >= NINT(p%t_measurement_interval*1000)
        m%MeasurementCurrentStep = 0
        LidarPosition_I =  u%LidarMesh%Position(1:3,1) + u%LidarMesh%TranslationDisp(1:3,1)

        CALL LidarSim_CalculateIMU(p, y, u)
        DO LoopGatesPerBeam = 0,p%GatesPerBeam-1
            !Transform measuring and lidar measurement position to the inertial system
            MeasuringPosition_I = LidarSim_TransformLidarToInertial(u%LidarMesh, p, p%MeasuringPoints_L(:,m%LastMeasuringPoint+LoopGatesPerBeam)) ! Calculate the Measuringpoint coordinate in the initial system

            !Line of Sight
            UnitVector    =   MeasuringPosition_I - LidarPosition_I             !Calculation of the Line of Sight Vector          
            UnitVector    =   UnitVector/NORM2(UnitVector)                      !=>Magnitude = 1

            !Calculation of the wind speed at the calculated position
            CALL LidarSim_CalculateVlos( p, UnitVector, Vlos, MeasuringPosition_I, LidarPosition_I, Time, IfW_p, IfW_ContStates, IfW_DiscStates, IfW_ConstrStates, IfW_OtherStates, IfW_m, TmpErrStat, TmpErrMsg) !Calculation of the line of sight wind speed
            CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)    
            CALL LidarSim_SetOutputs(y,p,m,Vlos,UnitVector,LoopGatesPerBeam,Time)    !Set all outputs to the output variable
        ENDDO                
        
        !Choosing which measuring point has to be calculated
        IF(m%LastMeasuringPoint+p%GatesPerBeam > SIZE(p%MeasuringPoints_L,2))THEN                      
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

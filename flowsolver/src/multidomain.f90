! *****************************************************************
! *****************************************************************
! *** multidomain module, contains outflo and inflow 0d models ***
! *****************************************************************
! *****************************************************************
!
      module multidomain
!
      use pointer_data, only: r0d
      use datatypes
      use debuggingTools

      implicit none
!
! --- by default define all types/subroutines/functions non-public
      private 
!
! --- public constructors
      public :: nrcrconstructor
      private :: nrcrconstructor_params, nrcrconstructor_file         
      interface nrcrconstructor
         module procedure nrcrconstructor_params, nrcrconstructor_file         
      end interface 

      public :: ntrcrconstructor
      public :: sysconstructor
      public :: multidomconstructor
      public :: hrtconstructor
      public :: controlledCoronaryConstructor
      public :: netlistLPNConstructor
!
! --- public type in order to access get/set
      public :: reducedorder
      public :: multidomaincontainer
      public :: advancedreducedorder
      public :: controlledCoronaryModel !\todo needed?
      public :: netlistLPN !\todo needed?

!
! --- declare public subroutines not bound to a type
      public :: setsimv_maxsurf
      public :: setsimv_delt
      public :: setsimv_alfi
      public :: setsimv_ntout
      public :: setsimv_nstep
      public :: setsimv_lstep
      public :: setsimv_rho
      public :: initmultidomain
      public :: multidomainstatus
      public :: startmultidomain
      private :: startmultidomain_file
      private :: startmultidomain_int
      private :: startmultidomain_int_char     
      interface startmultidomain
         module procedure startmultidomain_file,    &
                         startmultidomain_int,      &
                         startmultidomain_int_char
      end interface 
!
! --- rcr data types
      type rcr
         private
         real*8 :: rp
         real*8 :: c
         real*8 :: rd
         real*8 :: pd
      end type rcr
!
      type rcrdata
         private
         real*8 :: rp
         real*8 :: c
         real*8 :: rd   
         type(timedata) :: pd
      end type rcrdata
!
      type trcrdata
         private
         type(timedata) :: rp
         type(timedata) :: c
         type(timedata) :: rd   
         type(timedata) :: pd
      end type trcrdata
!
      type, extends(rcr) :: rcr_pbias
         private
         real*8, public :: pb
      end type rcr_pbias
!
! *******************************
! *** reduced order data type ***
! *******************************
!
      type reducedorder  
         private                            
         integer :: isactive = int(0)
         logical :: first_order_fd = .true.
         logical :: second_order_fd = .false.                
         integer :: updatepressure = int(0)
         integer :: surfnum = 0
         integer :: statenum
         integer, allocatable :: surfids(:)
         real*8, allocatable  :: surfarea(:)
         real*8, allocatable  :: flow_n(:)
         real*8, allocatable  :: flow_n1(:)
         real*8, allocatable  :: pressure_n(:)
         ! real*8, pointer  :: pressure_n(:) => null()
         type(r0d), allocatable :: pressure_n_ptr(:)



         real*8, allocatable  :: implicitcoeff(:,:)    ! at n+alfi
         real*8, allocatable  :: implicitcoeff_n1(:,:) ! at n+1
         real*8, allocatable :: flowhist(:,:)
         real*8, allocatable :: pressurehist(:,:)
         character(len=50) :: flowfile
         character(len=50) :: pressurefile         
         character(len=50), allocatable :: variablesfile(:)
         character(len=50), public :: classNameString
         contains
         private 
         procedure, public :: loadflowfile => loadflowfile
         procedure, public :: loadpressurefile => loadpressurefile
         procedure, public :: getsurfnum => getsurfnum
         procedure, public :: getsurfids => getsurfids
         procedure, public :: hassurfid => hassurfid
         procedure, public :: setarea => setarea
         procedure, public :: getarea => getarea
         procedure, public :: setflow_n => setflow_n   ! sets flow at previous time step 
         procedure, public :: setflow_n1 => setflow_n1 ! sets current flow at n+alf, updated during multicorrect
         procedure, public :: setpressure_n => setpressure_n   
         procedure, public :: updpressure_n1_withflow
         procedure, public :: updpressure_n1_withvalue
         generic, public   :: updatepressure_n1 => updpressure_n1_withflow, &
                                                   updpressure_n1_withvalue         
         procedure, public :: ispressureupdate => ispressureupdate     
         procedure, public :: getimplicitcoeff => getimplicitcoeff
         procedure, public :: setimplicitcoeff => setimplicitcoeff              
         procedure, public :: initxvars => initxvars
         procedure, public :: updxvars => updxvars         
         procedure, public :: writexvars => writexvars
         procedure, public :: setpresspntrs => setpresspntrs_ro

         procedure, public :: assign_ptrs_ext => assign_ptrs_ext
      end type reducedorder

      
!           
!
! ****************************
! *** reduced order models ***
! ****************************
!
! *** constant rcr with time varying pdistal
!
      type, extends(reducedorder) :: numericalrcr
         real*8, allocatable  :: parameters_RCR(:,:) ! RCR parameter array - for filter
         real*8               :: parameters_Pd       ! Pd parameter - for filter
         private 
         type(rcrdata), allocatable :: rcrparams(:) ! RCR parameter list
         real*8, allocatable        :: pRes_n(:)    ! Reservoir pressure
         integer                    :: init_pRes    ! Integer test
         contains
         private
         procedure :: initialise_rcr => initialise_rcr 
         procedure :: setimplicitcoeff_rcr => setimplicitcoeff_rcr
         procedure :: updxvars_rcr => updxvars_rcr         
         procedure :: writexvars_rcr => writexvars_rcr
         procedure :: assign_ptrs_ext_rcr => assign_ptrs_ext_rcr
      end type numericalrcr

      ! real*8, allocatable  :: parameters_RCR(:,:) ! RCR parameter array - for filter
      ! real*8               :: parameters_Pd       ! Pd parameter - for filter    
      ! real*8, allocatable  :: pressure_n(:)



!
! *** time varying rcr with time varying pdistal
!
      type, extends(reducedorder) :: numericaltrcr
         private 
         type(trcrdata), allocatable :: trcrparams(:)
         type(rcr), allocatable :: currparams(:)
         type(rcr), allocatable :: initparams(:)
         logical, allocatable :: isfeedback(:)
         real*8, allocatable :: feedbackcontrol(:)
         contains
         private
         procedure :: initialise_trcr => initialise_trcr 
         procedure :: setimplicitcoeff_trcr => setimplicitcoeff_trcr
         procedure, pass, public :: setfeedbacksurfs => setfeedbacksurfs
         procedure, pass, public :: setfeedbackcontrol => setfeedbackcontrol 
      end type numericaltrcr
!
! *** heart model 
!
      type, extends(reducedorder) :: numericalheart
         private
         real*8 :: patrial                      ! atrial pressure
         real*8 :: edv                          !
         real*8 :: rmv, lmv                     ! mitral valve resistance and inductor         
         real*8 :: qmv_n                        ! mitral valve flow at t=n and t=n+1         
         real*8 :: rav, lav                     ! aortic valve resistance and inductor
         real*8 :: vlv_n                        ! left ventricular volume at t=n and t=n+1
         real*8 :: plv_n                        ! left ventricular pressure at t=n and t=n+1         
         real*8 :: vulv                         ! left ventricular unstressed volume
         real*8 :: pla_n                        ! left atrial pressure at t=n and t=n+1         
         real*8 :: vula                         ! left atrial unstressed volume
         real*8 :: emax, emin                   ! maximum/minimum elastance
         real*8 :: period                       ! heart period
         real*8 :: tmax                         ! 
         real*8 :: trelax                       ! 
         real*8 :: kelv                         ! defined value for elastance varying resistance 
         real*8 :: vlv_coeff(2)                 ! vlv coefficients
         real*8 :: activationtime               ! elastance activation time
         real*8, allocatable :: vlv_hist(:)     ! vlv history
         real*8, allocatable :: plv_hist(:)     ! plv history         
         real*8, allocatable :: qmv_hist(:)     ! qmv history       
         real*8, allocatable :: elv_hist(:)     ! elv history            
         real*8, allocatable :: eval_hist(:,:)  ! eval history            
         real*8, allocatable :: stab_hist(:)    ! stabilisation pressure history
         integer :: inputHRandSP                ! Flag to indicate use of HeartRate.dat and SystolicPressure.dat input files
         type(linkedlist), pointer :: heartRateList ! prescribed heart rate data, when we want to input this from data
         type(linkedlist), pointer :: systolicPressureList ! prescribed peak systolic pressure array, for when we want to VERY APPROXIMATELY prescribe this from input data
         real*8 :: totalTimeOfAllCompletedPeriods ! provides the time at which the present beat started
         integer :: avopen = int(0)             ! aortic valve logical              
         type(timedata) :: elv_senzaki
         type(timedata) :: elv_input
         integer :: input_elastance         
         real*8 :: evalval 
         integer :: evalcount
         real*8, pointer :: sPress => null()
         integer :: backflow = int(0)             ! backflow logical
         real*8  :: m_backflow                    ! backflow magnitude
         real*8  :: s_backflow                    ! backflow steepness 
         real*8  :: c_backflow                    ! backflow closure 
         real*8  :: t_backflow                    ! time spent in backflow         
         real*8  :: max_backflow = real(39.5e-3,8)! maximum backflow time = 39.5 ms from Leyh et al. Circulation, 1999, 100, 2153-2160          
         real*8, allocatable :: act_hist(:)       ! activation history
         integer :: ibackflow = int(0)             !
         integer :: timestepsSinceLastHeartPeriodReset ! counter to store how far through the current heart beat we are
         integer :: newBeatJustStarted            ! This should only be set to 1 on the first timestep of the current beat, starting from the /second/ beat of the simulation. Otherwise, zero.
         character(len=50) :: ahistfilename     ! activation time, normalised by period [0-1]
         character(len=50) :: avarsformat       ! 
         contains
         private         
         procedure :: initialise_hrt => initialise_hrt
         procedure :: initxvars_hrt => initxvars_hrt
         procedure :: updxvars_hrt => updxvars_hrt
         procedure :: setimplicitcoeff_hrt => setimplicitcoeff_hrt        
         procedure, public :: iterate_hrt
         procedure :: writexvars_hrt => writexvars_hrt               
         procedure, public :: isavopen => isavopen         
         procedure :: getelastance_senzaki         
         procedure :: setpresspntrs_hrt => setpresspntrs_hrt
         procedure, public :: update_activationtime_hrt => update_activationtime_hrt
         procedure, public :: write_activation_history_hrt => write_activation_history_hrt
         procedure :: read_activation_history_hrt => read_activation_history_hrt
         procedure, public :: set_sPress              
         procedure :: assign_ptrs_ext_hrt => assign_ptrs_ext_hrt
      end type numericalheart

      interface
         function iterate_interface(a) result(b)
         import :: numericalheart
         class(numericalheart) :: a
         real*8 :: b(2)
         end function 
         function elastance_interface(aa) result(bb)
         import :: numericalheart
         class(numericalheart) :: aa
         real*8 :: bb
         end function                   
      end interface        

      abstract interface
         function return_scalar(a) result(b)
            real*8 :: a
            real*8 :: b
         end function

      end interface 

      abstract interface
         subroutine lv_ao_circuit(a)
         import :: numericalheart
         class(numericalheart) :: a
         end subroutine
      end interface
!         
!
! ************************************************************
! *** feedback data type, stores feedback input and output ***
! ************************************************************
!
      type feedback
         private
         logical :: active = .false.
         real*8 :: activation_tol
         integer :: surfnum
         integer, allocatable :: surfids(:)
         type(pntr), allocatable, public :: pntr(:)
         real*8, allocatable :: input(:,:)
         real*8, allocatable :: output(:,:)
         real*8, allocatable :: outputval(:)
         real*8, allocatable :: targetval(:)
         real*8, allocatable :: params(:)   
         real*8, allocatable :: tolhist(:,:)
         type(pntr) :: hr_ctrl
         procedure(feedback_interface), pointer :: calculate => null()
         procedure(feedback_interface), pointer :: update => null()
         contains
         private
         procedure :: initialise => initialise_feedback
         procedure :: isactive => isactive_feedback
         procedure :: setcntlpntrs => setcntlpntrs_feedback
      end type 

      abstract interface
         subroutine feedback_interface(a,b)
         import :: feedback
         class(feedback) :: a
         integer :: b
         end subroutine
      end interface
!      
! **************************
! *** control 
! **************************
!
      type control                                                      ! first order ode: tau dx/dt + x  = f(t)
         private                                                        ! 
         logical :: active                                              ! active logical
         logical :: first_order_fd = .true.                             ! finite difference logical 
         logical :: second_order_fd = .false.                           ! finite difference logical 
         integer :: dim                                                 ! control dimension
         real*8, pointer :: n(:) => null()                              ! control at t_{n} 
         real*8, pointer :: n1(:)  => null()                            ! control at t_{n+1}
         real*8, allocatable :: tau(:)                                  ! ode rate
         real*8, allocatable :: params(:,:)                             ! control parameters
         real*8, allocatable :: hist(:,:)                               ! 
         procedure(control_interface), pointer :: getsource => null()   ! pointer variable to procedure         
         contains
         procedure :: initialise => initialise_control
         procedure :: set => set_control
         procedure :: update => update_control
      end type

      abstract interface                  ! abstract interface
         function control_interface(a,b)  ! interface name
         import :: control, feedback      ! import derived types
         class(control) :: a              ! derived type
         type(feedback) :: b              ! 
         real*8 :: control_interface(a%dim)   !
         end function
      end interface           
!
! ********************************************
! *** extended advanced reduced order type ***
! ********************************************
!
      type, extends(reducedorder) :: advancedreducedorder
         private
         type(pntr), allocatable, public :: flowpntr(:)
         type(pntr), allocatable, public :: presspntr(:)
         type(pntr), allocatable, public :: dflowpntr(:)
         type(pntr), allocatable, public :: stabpntr(:)
         type(feedback) :: feedback  
         type(control) :: control
         integer :: xdim
         real*8, allocatable :: x_n(:)
         real*8, allocatable :: x_n1(:)
         real*8, allocatable :: x_hist(:,:)     ! x-variables
         real*8, allocatable :: p_hist(:,:)     ! pressure
         real*8, allocatable :: q_hist(:,:)     ! flow
         real*8, allocatable :: stab_hist(:,:)  ! stabilisation history
         real*8, allocatable :: r_hist(:,:,:)   ! restart
         character(len=50) :: xvarsfilename     !  
         character(len=50) :: phistfilename     !
         character(len=50) :: qhistfilename     !
         character(len=50) :: chistfilename     ! control  
         character(len=50) :: fhistfilename     ! feedback  
         character(len=50) :: xvarsformat       !
         character(len=50) :: pvarsformat       !
         character(len=50) :: qvarsformat       !
         character(len=50) :: cvarsformat       ! control  
         character(len=50) :: rvarsformat       ! restart
         character(len=50) :: fvarsformat       ! restart
         character(len=50) :: tvarsformat       ! tolerance
         logical :: feedbackactive = .false. 
         contains                
         procedure, public :: setflowpntrs => setflowpntrs
         procedure, public :: setpresspntrs => setpresspntrs              
         procedure, public :: setfeedbackpntrs => setfeedbackpntrs           
         procedure, public :: solve => solve
         procedure, public :: update => update
         procedure, public :: isfeedbackactive => isfeedbackactive
         procedure, public :: writeRestart => writeRestart
      end type advancedreducedorder
!
!
!

!
! *************************************
! *** advanced reduced order models ***
! *************************************
!
! *** systemic circuit 
!
      type, extends(advancedreducedorder) :: systemiccircuit
         private                               
         integer :: rcrsurfnum                  ! number of 3D branches         
         real*8, allocatable :: rp(:)           ! proximal resistance for each 3D branch   
         real*8, allocatable :: c(:)            ! capacitance for each 3D branch
         real*8, allocatable :: rd(:)           ! distal resistance for each 3D branch
         real*8 :: ca3, vua3, ra3               ! arterial section 3 [rc]
         real*8 :: cv1, vuv1, rv1               ! venous section 1 [rc]
         real*8 :: cv2, vuv2, rv2, lv2          ! venous section 2 [rcl]
         real*8 :: rmv, lmv, qmv_n              ! mitral valve resistance/inductor/flow 
         real*8 :: rav, lav                     ! aortic valve resistance/inductor
         real*8 :: vlv_n, vulv, plv_n           ! left ventricular volume/pressure
         real*8 :: elv_n, elv_n1                ! left ventricular elastances
         real*8 :: vla_n, vula, ela, pla_n      ! left atrial volume/elastance/pressure 
         real*8 :: emax, emin                   ! max/min elastance         
         real*8 :: tmax, trel                   ! time to max elastance/diastolic relaxation time
         real*8 :: period                       ! heart period
         real*8 :: activationtime               ! elastance activation time
         real*8 :: edv, esv                     ! end diastolic/systolic volumes
         integer :: va3indx, vv1indx, vv2indx   ! arterial/venous volume indices
         integer :: vlaindx, vlvindx            ! atrial/ventricular volume indices 
         integer :: qv2indx, qmvindx            ! venous/mitral flow indices
         integer :: paoindx                     ! aortic pressure index
         integer :: readinitialvalues           ! integer test to read initial values
         integer :: vardim                      ! x-variable array dimension         
         logical :: avopen, mvopen              ! valve states
         real*8, allocatable :: elv_hist(:)     ! elastance history
         real*8 :: paorta_n, qaorta_n           ! aortic pressure and flow
         real*8 :: vlvcoeff(2)                  ! vlv coefficients
         integer :: hrindx = 1                  ! heart rate
         integer :: emindx = 2                  ! max elastance
         integer :: rdindx = 3                  ! arterial resistance
         integer :: cpindx = 4                  ! venous compliance
         integer :: vuindx = 5                  ! venous unstressed volume
         character(len=50) :: ahistfilename     ! activation time, normalised by period [0-1]
         character(len=50) :: avarsformat       ! 
         real*8, allocatable :: act_hist(:)     ! activation history
         logical :: closedloop = .false.        ! closed loop
         real*8 :: patrial
         real*8 :: kelv = 0.0005                ! kplv
         real*8 :: tmaxratio, trelratio         ! tmax and trel ratios
         contains
         private
         procedure :: initialise_sys => initialise_sys        
         procedure :: setflowpntrs_sys => setflowpntrs_sys
         procedure :: setpresspntrs_sys => setpresspntrs_sys
         procedure :: setfeedbackpntrs_sys => setfeedbackpntrs_sys
         procedure :: initxvars_sys => initxvars_sys
         procedure :: setimplicitcoeff_sys => setimplicitcoeff_sys         
         procedure :: setsystem_sys => setsystem_sys
         procedure :: setsystem_sys_open => setsystem_sys_open
         procedure :: solve_sys => solve_sys
         procedure :: solve_sys_open => solve_sys_open
         procedure :: update_sys => update_sys
         procedure :: update_sys_open => update_sys_open
         procedure :: getelastance => getelastance_sys
         procedure, public :: writefile => writefile_sys        
         procedure :: update_activationtime => update_activationtime         
      end type systemiccircuit
!
! **********************************************************************
! *** multidomain container used by advanced reduced order data type ***
! **********************************************************************
!
!     here the flow_alfi contains the flows at the current n+alf_{i} step
!     in order to allow the reduced order models to access these value/values
!     this has been defined as a pointer rather than allocatable 
!
      type multidomaincontainer
         private
         integer :: surfnum = 0
         integer, pointer :: surfids(:) => null()
         real*8, pointer :: flow(:) => null()
         real*8, pointer :: flowderivative(:) => null()
         real*8, pointer :: pressure(:) => null()
         real*8, pointer :: areas(:) => null()
         real*8, pointer :: stabilisationpressure(:) => null()
         real*8, pointer :: stb_pres(:)  => null()
         contains
         procedure :: addsurfids => addsurfids
         procedure :: setflow => setflow_mdc
         procedure :: setflowderivative => setflowderivative_mdc
         procedure :: setarea => setarea_mdc
         procedure :: setpressure => setpressure_mdc
         procedure, public :: getsurfnum => getsurfnum_mdc
         procedure, public :: getsurfids => getsurfids_mdc
         procedure, public :: getarea => getarea_mdc    
         procedure, public :: resetstb_pres => resetstb_pres
         procedure, public :: addstb_pres => addstb_pres
         procedure, public :: sumstb_pres => sumstb_pres
         procedure, public :: getstb_pres => getstb_pres 
      end type                        

!     ***************~~~~BEGIN CORONARY MODEL TYPES~~~~***************
      type, extends(advancedreducedorder) :: controlledCoronaryModel
            private
            real*8, allocatable :: O2supplyDemandHistoryWindowLength_seconds(:) ! width of history over which to consider myocardial O2 supply/demand discrepancy
            integer, allocatable :: O2supplyDemandHistoryWindowLength_timesteps(:) ! the length in time-steps of the moving average window for O2 supply-demand averaging
            integer :: numberOfControlledCoronaryModelSurfaces
            integer, allocatable :: localToGlobalSurfaceIndexMap(:)
            real*8, allocatable :: currentMVO2(:) ! stores myocardial oxygen demand for the myocardial region perfused by this outlet, on the present time-step (number of corntrolledCoronaryModel outlets)
            real*8, allocatable :: currentMyocardialHungerSignal(:) ! stores the moving average of the O2 demand - supply (upper case H(t) in the paper, when used in this mode)
            real*8, allocatable :: hungerDelta(:) ! lower-case h(t) in the paper
            real*8, allocatable :: O2supplyDemandHistory(:,:) ! History of the MVO2 - O2 supply (data,number of controlledCoronaryModel outlets)
            real*8, allocatable :: oneOverR(:) ! 1/total resitance of each outlet LPN. "one over" is for ease of writing down the ODEs
            real*8 :: oneOverTotalResistanceOfModel
            real*8 :: arterial_O2_volume_proportion
            real*8, allocatable :: feedbackGain(:) ! metabolic feedback vasodilation
            real*8, allocatable :: betaFeedforwardGain(:) ! vasodilation on microvasculature
            real*8, allocatable :: alphaFeedforwardGain(:) ! vasoconstriction on the small arteries
            real*8, allocatable :: dampingCoefficient(:) ! Oscillation-damping in the harmonic control method
            real*8 :: intramyocardialPressureToLVScaling ! proportion of LV pressure applied to "ground terminal" of C_im capacitor
            real*8, allocatable :: P_1(:) ! pressure at circuit-side terminal of capacitor C_a (used in setimplicitcoeff_controlledCoronary)
            real*8, allocatable :: P_2(:) ! pressure at circuit-side terminal of capacitor C_im (used in setimplicitcoeff_controlledCoronary)
            real*8, allocatable :: P_a(:) ! pressure at each coronary surface (badly-dubbed "aortic pressure", despite not being in the aorta)
            real*8, allocatable :: rp(:) ! proximal resistance (one for each coronary outlet.. same goes for each other electrical component parameter below!)
            real*8, allocatable :: rd(:) ! distal resistance
            real*8, allocatable :: ra(:) ! entrance resistance at coronary outlet (arterial resistance)
            real*8, allocatable :: c_im(:)  ! intramyocardial capacitance (grounded on the cyclic intramyocardial pressure)
            real*8, allocatable :: c_a(:)  ! epicardial artery capacitance (grounded to zero)
            real*8, allocatable :: deltaMVO2(:) ! The change in MVO2 since the last time-step for each outlets perfused region
            real*8, allocatable :: MVO2(:) ! MVO2 for the region perfused by each outlet
            real*8, allocatable :: MVO2previousDt(:) ! MVO2 on the previous time step for the region perfused by each outlet
            real*8, allocatable :: MVO2computedAtEndOfBeatPreviousToPreviousBeat(:)
            real*8, allocatable :: MVO2computedAtEndOfPreviousBeat(:)
            real*8, allocatable :: deltaMVO2PerTimestepOverPreviousBeat(:) ! for feedforward; to track MVO2 changes per-dt
            real*8, allocatable :: proportionOfMyocardiumPerfusedByOutlet(:) ! used to assign proportion of MVO2 each outlet needs to provide flow to satisfy
            real*8, allocatable :: R_min(:) ! minimum total effective resistances of outlets
            real*8, allocatable :: R_max(:) ! maximum total effective resistances of outlets
            real*8 :: P_IM_mid ! mid myocardial extravascular pressure
            real*8 :: P_IM_mid_lasttimestep ! mid myocardial extravascular pressure at last timestep
            contains
            private
            procedure :: initialiseCoronaryControlModel => initialiseCoronaryControlModel
            procedure :: computeHungerAndSetHistories => computeHungerAndSetHistories_negativeHungerCapped
            procedure :: computeResistanceChanges => computeResistanceChanges_harmonic
            procedure :: updateFeeds => updateFeedforwardAndFeedbackInCoronaryControlModel
            procedure :: updateMVO2 => updateMVO2
            procedure :: setflowpntrs_coronary => setflowpntrs_coronary
            procedure :: setimplicitcoeff_controlledCoronary => setimplicitcoeff_controlledCoronary
            procedure, public :: updateLPN_coronary => updateLPN_coronary
            procedure, public :: setSurfacePressure_coronary => setSurfacePressure_coronary
            procedure, public :: getNumberOfControlledCoronaryModelSurfaces => getNumberOfControlledCoronaryModelSurfaces
            procedure, public :: getRpByIndex => getRpByIndex
            procedure, public :: getRdByIndex => getRdByIndex
            procedure, public :: getP_1ByIndex => getP_1ByIndex
            procedure, public :: getP_2ByIndex => getP_2ByIndex
            procedure, public :: getcurrentMeanMyocardialOxygenHungerByIndex => getcurrentMeanMyocardialOxygenHungerByIndex
            procedure, public :: getMVO2previousByIndex => getMVO2previousByIndex
            procedure, public :: getP_aByIndex => getP_aByIndex
            procedure :: writeO2supplyDemandHistory => writeO2supplyDemandHistory
            procedure :: writeDebugHistories => writeDebugHistories
            procedure, public :: setRpByIndex => setRpByIndex
            procedure, public :: setRdByIndex => setRdByIndex
            procedure, public :: setP_1ByIndex => setP_1ByIndex
            procedure, public :: setP_2ByIndex => setP_2ByIndex
            procedure, public :: setcurrentMeanMyocardialOxygenHungerByIndex => setcurrentMeanMyocardialOxygenHungerByIndex
            procedure, public :: setMVO2previousByIndex => setMVO2previousByIndex
            procedure, public :: loadO2supplyDemandHistory => loadO2supplyDemandHistory
            procedure :: writeRestart_controlledCoronaryModel => writeRestart_controlledCoronaryModel
      end type controlledCoronaryModel

! *************** END CORONARY MODEL TYPES ***************

      type, extends(advancedreducedorder) :: netlistLPN
         private
         real*8, allocatable :: systemMatrix(:,:,:)
         real*8, allocatable :: RHS(:,:)
         real*8, allocatable :: pressuresInLPN(:,:)                       ! Pressure at each LPN node, using the same node indexing as in the netlist
         real*8, allocatable :: historyPressuresInLPN(:,:)                ! As pressuresInLPN, but for any nodes with histories. /Most/ of the entries in this array will never be used.
         real*8, allocatable :: flowsInLPN(:,:)                           ! Flow through each component in the LPN, in the order they appear in the netlist
         real*8, allocatable :: solutionVector(:,:)
         integer, allocatable :: numberOfComponents(:)
         integer, allocatable :: numberOfPressureNodes(:)
         integer :: numberOfLPNSurfaces
         integer, allocatable :: localToGlobalSurfaceIndexMap(:)
         character, allocatable :: circuitData_componentTypes(:,:) ! the data in here will be the stripped first column of the netilst, identifying each line of circuitData as being r=resistor, c=capacitor, etc.
         real*8, allocatable :: circuitData(:,:,:)
         integer, allocatable :: nodeIndexToPressureHistoryNodeOrderingMap(:,:)
         integer, allocatable :: componentIndexToFlowHistoryComponentOrderingMap(:,:)
         integer, allocatable :: systemSize(:)
         integer, allocatable :: listOfNodesWithMultipleIncidentCurrents(:,:)
         integer, allocatable :: numberOfMultipleIncidentCurrentNodes(:)
         integer, allocatable :: listOfPrescribedPressures(:,:)         ! input data, listing node indices of prescribed pressures (e.g. LV pressure on a capacitor, venous pressure in open loop, etc.).
         integer, allocatable :: listOfPrescribedFlows(:,:)             ! input data, listing node indices of prescribed flows, listing the component indices with prescribed flows (e.g. 3D outlet flow)
         integer, allocatable :: listOfHistoryPressures(:,:)            ! generated from input data, listing pressure node indices and component flow indices where a history is needed (i.e. last time-step values for capacitors/inductors)
         integer, allocatable :: listOfHistoryFlows(:,:)
         integer, allocatable :: numberOfPrescribedPressures(:)
         real*8, allocatable :: valueOfPrescribedPressures(:,:)
         real*8, allocatable :: valueOfPrescribedFlows(:,:)
         character, allocatable :: typeOfPrescribedPressures(:,:)
         integer, allocatable :: numberOfPrescribedFlows(:)
         character, allocatable :: typeOfPrescribedFlows(:,:)
         integer, allocatable :: numberOfPrescribedPressuresAndFlows(:)           ! Just the sum of the previous two declared integers
         integer,allocatable :: numberOfHistoryPressures(:)
         integer,allocatable :: numberOfHistoryFlows(:)
         integer, allocatable :: columnMap(:,:)
         integer :: maxColumnMapSize
         real*8, allocatable :: inverseOfSystemMatrix(:,:,:)
         real*8, allocatable :: P_a(:)
         integer :: maxSystemSize
         integer :: maxComponents
         integer :: maxPressureNodes
         integer :: maxPrescribedPressureNodes
         integer :: maxPrescribedFlows
         integer, allocatable :: columnIndexOf3DInterfaceFlowInLinearSystem(:)
        
         contains
         private

         procedure :: initialiseNetlistLPN => initialiseNetlistLPN
         procedure :: setimplicitcoeff_netlistLPN => setimplicitcoeff_netlistLPN
         procedure :: setflowpntrs_netlistLPN => setflowpntrs_netlistLPN
         procedure, public :: generateLinearSystemFromPrescribedCircuit => generateLinearSystemFromPrescribedCircuit
         procedure, public ::  getListOfNodesWithMultipleIncidentCurrents => getListOfNodesWithMultipleIncidentCurrents
         procedure :: getMapOfPressHistoriesToCorrectPressNodes => getMapOfPressHistoriesToCorrectPressNodes
         procedure :: getMapOfFlowHistoriesToCorrectComponents => getMapOfFlowHistoriesToCorrectComponents
         procedure, public :: assembleRHS_netlistLPN => assembleRHS_netlistLPN
         procedure, public :: updateLPN_netlistLPN => updateLPN_netlistLPN
         procedure, public :: setSurfacePressure_netlistLPN => setSurfacePressure_netlistLPN
         procedure, public :: writeInternalLPNPressuresAndFlows => writeInternalLPNPressuresAndFlows
         procedure :: writeRestart_netlistLPNSurfaces => writeRestart_netlistLPNSurfaces

      end type netlistLPN


!
! ***********************************************************************
! *** public objects and internal interfaces contained in this module ***
! ***********************************************************************
!
!     ! data types
      type(numericalrcr), public, save :: nrcr 
      type(numericaltrcr), public, save :: ntrcr
      type(systemiccircuit), public, save :: sys
      type(multidomaincontainer), public, save :: multidom
      type(numericalheart), public, save :: hrt
      type(controlledCoronaryModel), public, save :: controlledCoronarySurfaces
      type(netlistLPN), public, save :: netlistLPNSurfaces
!
!     ! simvascular parameters
      integer :: maxsurf   ! maximum number of surfaces to correspond to the arrays used in getflow, etc.
      real*8 :: alfi       ! \alpha_{f}
      real*8 :: delt       ! \Delta{t}
      real*8 :: rho        ! \rho_{fluid}
      integer :: ntout     ! frequency of output inputted from solver.inp 
      integer :: nstep     ! number of timesteps to do (inputted from solver.inp
      integer :: lstep     ! \todo remove this probably, pass lstep in with hrtconstructor and sysconstructor (because some functions pass lstep in explicitly, if we use lstep in a function expecting it to have the right value and without passing it into that function, the compiler thinks this is ok because of this here lstep; in reality this one is never updated though and the function is then buggy!) total number of steps completed (read from numstart.dat?)
      real*8, parameter :: mmhgtodynes = real(1333.3,8)
      real*8, parameter :: pi = real(4.0d+0,8) &
                              * atan(1.0d+0)
! 
!     ! public logicals
!      
      integer, public :: multidomainactive = int(0)
      integer, public :: nrcractive = int(0)
      integer, public :: ntrcractive = int(0)
      integer, public :: sysactive = int(0)
      integer, public :: hrtactive = int(0)
      integer, public :: newCoronaryActive = int(0)
      integer, public :: netlistActive = int(0)
!
!     ! matrix/matrix and matrix/vector multiplication interface
!
      interface multiply_ab
         module procedure matrixmatrix
         module procedure matrixvector
         module procedure matrixvectorpointer
      end interface 
!
      interface solve_axb
         module procedure solvematrixvector
         module procedure solvematrixmatrix
      end interface
!
      contains

!     ***************~~~~BEGIN CORONARY MODEL CODE~~~~***************

      function controlledCoronaryConstructor(numCoronarySrfs,nSrfListCoronary,numCalcSrfs,nsrflistCalc) result(a)
         implicit none
         integer, intent(in) :: numCoronarySrfs
         integer, intent(in) :: nSrfListCoronary(0:maxsurf)
         integer, intent(in) :: numCalcSrfs
         integer, intent(in) :: nsrflistCalc(0:maxsurf)
         type(controlledCoronaryModel) :: a
         call a%initialiseCoronaryControlModel(numCoronarySrfs,nSrfListCoronary,numCalcSrfs,nsrflistCalc)
         return
      end function controlledCoronaryConstructor

      subroutine initialiseCoronaryControlModel(this,numCoronarySrfs,nSrfListCoronary,numCalcSrfs,nsrflistCalc)

         implicit none
         integer, intent(in) :: numCoronarySrfs
         integer, intent(in) :: nSrfListCoronary(0:maxsurf)
         integer, intent(in) :: numCalcSrfs
         integer, intent(in) :: nsrflistCalc(0:maxsurf)
         class(controlledCoronaryModel) :: this
         integer :: ii, jj
         integer :: numberOfLines, iostatus
         integer :: hstep
         logical :: file_exists ! for testing whether there are restart files on the disk
         real*8 readvalue ! just a temporary variable to store reads before saving them in the appropriate places


         ! A label so we can identify this class (avoids a world of pain when working with a tower of derived types)
         this%classNameString = 'controlledCoronaryModel'

         ! set surface numbers and lists   
         this%numberOfControlledCoronaryModelSurfaces = numCoronarySrfs
         this%surfnum = numCoronarySrfs !\todo check this and line above both necessary! (this line is needed)
         allocate(this%surfids(numCoronarySrfs))
         this%surfids(1:numCoronarySrfs) = nSrfListCoronary(1:numCoronarySrfs)

         allocate(this%flowpntr(this%numberOfControlledCoronaryModelSurfaces))

         allocate(this%MVO2previousDt(this%numberOfControlledCoronaryModelSurfaces))
         allocate(this%MVO2(this%numberOfControlledCoronaryModelSurfaces))
         allocate(this%currentMyocardialHungerSignal(this%numberOfControlledCoronaryModelSurfaces))
         allocate(this%hungerDelta(this%numberOfControlledCoronaryModelSurfaces))
         allocate(this%rp(this%numberOfControlledCoronaryModelSurfaces))
         allocate(this%rd(this%numberOfControlledCoronaryModelSurfaces))
         allocate(this%ra(this%numberOfControlledCoronaryModelSurfaces))
         allocate(this%c_im(this%numberOfControlledCoronaryModelSurfaces))
         allocate(this%c_a(this%numberOfControlledCoronaryModelSurfaces))
         allocate(this%deltaMVO2(this%numberOfControlledCoronaryModelSurfaces))
         allocate(this%deltaMVO2PerTimestepOverPreviousBeat(this%numberOfControlledCoronaryModelSurfaces)) ! There's some redundancy between this and deltaMVO2 at the time of writing, but the separation is likely to be useful in future modifications
         allocate(this%R_min(this%numberOfControlledCoronaryModelSurfaces))
         allocate(this%R_max(this%numberOfControlledCoronaryModelSurfaces))
         allocate(this%oneOverR(this%numberOfControlledCoronaryModelSurfaces))
         allocate(this%localToGlobalSurfaceIndexMap(this%numberOfControlledCoronaryModelSurfaces))
         allocate(this%P_1(this%numberOfControlledCoronaryModelSurfaces))
         allocate(this%P_2(this%numberOfControlledCoronaryModelSurfaces))
         allocate(this%P_a(this%numberOfControlledCoronaryModelSurfaces))
         allocate(this%MVO2computedAtEndOfBeatPreviousToPreviousBeat(this%numberOfControlledCoronaryModelSurfaces))
         allocate(this%MVO2computedAtEndOfPreviousBeat(this%numberOfControlledCoronaryModelSurfaces))
         allocate(this%proportionOfMyocardiumPerfusedByOutlet(this%numberOfControlledCoronaryModelSurfaces))
         allocate(this%feedbackGain(this%numberOfControlledCoronaryModelSurfaces))
         allocate(this%betaFeedforwardGain(this%numberOfControlledCoronaryModelSurfaces))
         allocate(this%alphaFeedforwardGain(this%numberOfControlledCoronaryModelSurfaces))
         allocate(this%dampingCoefficient(this%numberOfControlledCoronaryModelSurfaces))
         allocate(this%O2supplyDemandHistoryWindowLength_seconds(this%numberOfControlledCoronaryModelSurfaces))
         allocate(this%O2supplyDemandHistoryWindowLength_timesteps(this%numberOfControlledCoronaryModelSurfaces))

         this%deltaMVO2PerTimestepOverPreviousBeat(:) = real(0.0d0,8)

         ! these values should never be used as they are re-assigned before the first time they're read, but set zero anyway
         this%P_IM_mid = 0.0d0
         this%P_IM_mid_lasttimestep = 0.0d0

         do ii=1, this%numberOfControlledCoronaryModelSurfaces
            this%currentMyocardialHungerSignal(ii) = 0.0d0
            this%hungerDelta(ii) = 0.0d0
            this%P_a(ii) = 10000d0 !\todo think of a better way of initialising this (& restarting it!)
            this%P_1(ii) = 10000d0
            this%P_2(ii) = 10000d0
         enddo

!     ! zero variables \todo check all necessary:
         !     ! allocate and zero arrays 
         allocate(this%surfarea(this%surfnum))
         allocate(this%flow_n(this%surfnum))
         allocate(this%flow_n1(this%surfnum))
         allocate(this%pressure_n(this%surfnum))
         allocate(this%implicitcoeff(this%surfnum,2))
         allocate(this%implicitcoeff_n1(this%surfnum,2))
         this%surfarea(:) = real(0.0,8)
         this%flow_n(:) = real(0.0,8)
         this%flow_n1(:) = real(0.0,8)
         this%pressure_n(:) = real(0.0,8)
         this%implicitcoeff(:,:) = real(0.0,8)
         this%implicitcoeff_n1(:,:) = real(0.0,8)

         this%intramyocardialPressureToLVScaling = 1.0d0 ! \todo try other values here (0.4 in MATLAB)


         inquire(file="controlled_coronaries.dat", exist=file_exists)
         if (file_exists) then
            ! Do some sanity-checking on the input file:
            open(unit=73,file='controlled_coronaries.dat',status='old')
            numberOfLines = 0
            iostatus = 1
            do while (iostatus.ge.0)
               read(73,*,IOSTAT=iostatus)
            !   write(*,*) iostatus
               numberOfLines = numberOfLines + 1
            end do
            ! Check the file length was as expected:
            if ((numberOfLines-3).ne.(this%numberOfControlledCoronaryModelSurfaces*31)) then
               write(*,*) 'ERROR: controlled_coronaries.dat format is incorrect. Example file:'
               write(*,*) '# File contains the coronary model parameters for each outlet.'
               write(*,*)    '#### Coronary 1 ####'
               write(*,*)    '# R_a:'
               write(*,*)    '1.286208e5'
               write(*,*)    '# R_p:'
               write(*,*)    '1.286208e5'
               write(*,*)    '# R_d:'
               write(*,*)    '7.5e5'
               write(*,*)    '#C_a:'
               write(*,*)    '4.5e-7'
               write(*,*)    '# C_im:'
               write(*,*)    '2.7e-6'
               write(*,*)    '# R_min:'
               write(*,*)    '10'
               write(*,*)    '# R_max:'
               write(*,*)    '1e8'
               write(*,*)    '# MVO2 on previos time-step'
               write(*,*)    '0.25'
               write(*,*)    '# MVO2 on current time-step'
               write(*,*)    '0.25'
               write(*,*)    '# MVO2 for this outlets perfusion territory on previous time-step:'
               write(*,*)    '0.0272'
               write(*,*)    '# MVO2 for this outlets perfusion territory on current time-step:'
               write(*,*)    '0.0272'
               write(*,*)    '# Proportion of myocardium perfused by this outlet (decimal between 0 and 1):'
               write(*,*)    '0.1'
               write(*,*)    '# Metabolic feedback control gain:'
               write(*,*)    '2.5'
               write(*,*)    '# Alpha adrenergic feedforward control gain:'
               write(*,*)    '-0.1'
               write(*,*)    '# Beta adrenergic feedforward control gain:'
               write(*,*)    '0.25'
               write(*,*)    '# Myocardial O2 demand/supply integration window duration (seconds)'
               write(*,*)    '5'
               write(*,*)    '#### Coronary 2 ####'
               write(*,*)    '# R_a:'
               write(*,*)    '1.286208e5'
               write(*,*)    '# R_p:'
               write(*,*)    '1.286208e5'
               write(*,*)    '# R_d:'
               write(*,*)    '7.5e5'
               write(*,*)    '#C_a:'
               write(*,*)    '4.5e-7'
               write(*,*)    '# C_im:'
               write(*,*)    '2.7e-6'
               write(*,*)    '# R_min:'
               write(*,*)    '10'
               write(*,*)    '# R_max:'
               write(*,*)    '1e8'
               write(*,*)    '# MVO2 for this outlets perfusion territory on previous time-step:'
               write(*,*)    '0.0272'
               write(*,*)    '# MVO2 for this outlets perfusion territory on current time-step:'
               write(*,*)    '0.0272'
               write(*,*)    '# Proportion of myocardium perfused by this outlet (decimal between 0 and 1):'
               write(*,*)    '0.1'
               write(*,*)    '# Metabolic feedback control gain:'
               write(*,*)    '0.31'
               write(*,*)    '# Alpha adrenergic feedforward control gain:'
               write(*,*)    '-0.1'
               write(*,*)    '# Beta adrenergic feedforward control gain:'
               write(*,*)    '1.0'
               write(*,*)    '# Control damping coefficient (only used in harmonic mode, but must always be present for reading this file)'
               write(*,*)    '0.35'
               write(*,*)    '# Myocardial O2 demand/supply integration window duration (seconds)'
               write(*,*)    '5'
               write(*,*)    '# This line must not be blank (SimVascular expects it when validating this file), must only occur at the end of the file, and there should be no whitespace after this line'
               stop
            end if
            ! Done sainity-checking the input file.


            call fseek(73,1,0)
            open(unit=73,file='controlled_coronaries.dat',status='old')
            read(73,*) ! Read & ignore comment line: "# File contains the coronary model parameters for each outlet."
            do ii=1, this%numberOfControlledCoronaryModelSurfaces
               read(73,*) ! Read & ignore comment line: "# Coronary ii:"
               read(73,*) ! Read & ignore comment line: "# R_a:"
               read(73,*) this%ra(ii)
               read(73,*) ! Read & ignore comment line: "# R_p:"
               read(73,*) this%rp(ii)
               read(73,*) ! Read & ignore comment line: "# R_d:"
               read(73,*) this%rd(ii)
               read(73,*) ! Read & ignore comment line: "# C_a:"
               read(73,*) this%c_a(ii)
               read(73,*) ! Read & ignore comment line: "# C_im:"
               read(73,*) this%c_im(ii)
               read(73,*) ! Read & ignore comment line: "# R_min:"
               read(73,*) this%R_min(ii)
               read(73,*) ! Read & ignore comment line: "# R_max:"
               read(73,*) this%R_max(ii)
               read(73,*) ! Read & ignore comment line: "# MVO2 for this outlets perfusion territory on previous time-step"
               read(73,*) this%MVO2previousDt(ii)
               read(73,*) ! Read & ignore comment line: "# MVO2 for this outlets perfusion territory current time-step"
               read(73,*) readvalue
               this%MVO2(ii) = readvalue
               this%MVO2computedAtEndOfBeatPreviousToPreviousBeat(ii) = readvalue
               this%MVO2computedAtEndOfPreviousBeat(ii) = readvalue
               read(73,*) ! Read & ignore comment line: "# Proportion of myocardium perfused by this outlet (decimal between 0 and 1)"
               read(73,*) this%proportionOfMyocardiumPerfusedByOutlet(ii)
               read(73,*) ! Read & ignore comment line: "# Metabolic feedback control gain:"
               read(73,*) this%feedbackGain(ii)
               read(73,*) ! Read & ignore comment line: "# Alpha adrenergic feedforward control gain:"
               read(73,*) this%alphaFeedforwardGain(ii)
               read(73,*) ! Read & ignore comment line: "# Beta adrenergic feedforward control gain:"
               read(73,*) this%betaFeedforwardGain(ii)
               read(73,*) ! Read & ignore comment line: "# Control damping coefficient (only used in harmonic mode, but must always be present for reading this file)"
               read(73,*) this%dampingCoefficient(ii)
               read(73,*) ! Read & ignore comment line: "# Myocardial O2 demand/supply integration window duration (seconds)"
               read(73,*) this%O2supplyDemandHistoryWindowLength_seconds(ii)
            enddo
            close(73)
         else
            write(*,*) 'ERROR: controlled_coronaries.dat not found. Exiting.'
            stop
         endif

         do ii=1, this%numberOfControlledCoronaryModelSurfaces
            this%O2supplyDemandHistoryWindowLength_timesteps(ii) = nint(this%O2supplyDemandHistoryWindowLength_seconds(ii)/delt) !nint rounds to the nearest integer
         end do

         allocate(this%O2supplyDemandHistory(nstep + maxval(this%O2supplyDemandHistoryWindowLength_timesteps) + 1, this%numberOfControlledCoronaryModelSurfaces))
           
         do ii=1,  this%numberOfControlledCoronaryModelSurfaces
            ! set the history of the supply-demand discrepancy at each outlet (note that this array is long enough for the entire simulation + a fake "negative time" history; we set this fake history now:
            !\todo worry about restart issues here.
            do jj=1, maxval(this%O2supplyDemandHistoryWindowLength_timesteps)  ! \todo this needs to be set so that it doesn't overwrite a real history on a restart
               this%O2supplyDemandHistory(jj,ii) = 0d0
            enddo
            ! initialise the moving average O2 supply/demand history to zero
            this%currentMyocardialHungerSignal(ii) = 0d0 ! \todo have a way of setting the initial discrepancy history to be non-zero
            ! set oneOverR to be one over the initial total effective resistance of the LPN
            this%oneOverR(ii) = 1/(this%rp(ii) + this%rd(ii) + this%ra(ii))
         enddo

         this%arterial_O2_volume_proportion = real(7.0d0,8)/real(40d0,8) ! percent \todo this could be dynamic!

!        Populate the map which converts coronary surface indices (local) to surface indices (global) - for FlowHist etc.
         do ii=1, this%numberOfControlledCoronaryModelSurfaces
            do jj=1, numCalcSrfs
               if(this%surfids(ii).eq.nsrflistCalc(jj)) then
                  this%localToGlobalSurfaceIndexMap(ii) = jj
               endif
            enddo
         enddo

! ----------------------NEW CORONARY RESTART CODE COMPONENT----------------------
         ! Test whether this is a restart; if so, overwrite the necessary values with those from disk:
         inquire(file="coronary_restart.dat", exist=file_exists)
         if (file_exists) then
         ! \todo confirm this works with multiple coronaries - it looks likely that it does, but should make certain.
            ! Save data necessary to restart the LPN models at coronary surfaces:
              ! this is done using type-bound public "get" procedures which allow read-only
              ! access to the private variables. (we don't actually need these procedures now, as this code is now internal to the multidomain module [it initially was outside, this is copy-paste, so this could be made cleaner!])
              open(unit=73,file='coronary_restart.dat',status='old')
              read(73,*) !'# ~=Controlled Coronary Model Restart File=~'
              do ii = 1,this%numberOfControlledCoronaryModelSurfaces
                read(73,*) !'# Coronary outlet index number:', ii
                read(73,*) !'# Resistance rp:'
                read(73,*) readvalue ! \todo should probably check if the read succeeded so we don't get silent failures here!
                call this%setRpByIndex(ii,readvalue)
                read(73,*) !'# Resistance rd:'
                read(73,*) readvalue
                call this%setRdByIndex(ii,readvalue)
                read(73,*) !'# Coronary LPN internal pressure P_1:'
                read(73,*) readvalue
                call this%setP_1ByIndex(ii,readvalue)
                read(73,*) !'# Coronary LPN internal pressure P_2:'
                read(73,*) readvalue
                call this%setP_2ByIndex(ii,readvalue)
                read(73,*) !'# currentMyocardialHungerSignal:'
                read(73,*) readvalue
                call this%setcurrentMeanMyocardialOxygenHungerByIndex(ii,readvalue)
                read(73,*) !'# MVO2previousDt:'
                read(73,*) readvalue
                call this%setMVO2previousByIndex(ii,readvalue)
                read(73,*) !'# MVO2'
                read(73,*) this%MVO2(ii)
                read(73,*) !'# P_a (LPN pressure at outlet itself):'
                read(73,*) readvalue
                this%P_a(ii) = readvalue
                read(73,*) !'# deltaMVO2PerTimestepOverPreviousBeat (used for MVO2 ramping on the present beat)'
                read(73,*) readvalue
                this%deltaMVO2PerTimestepOverPreviousBeat(ii) = readvalue
                read(73,*) !'# Total myocardial oxygen demand during the last complete beat:'
                read(73,*) this%MVO2computedAtEndOfPreviousBeat(ii)
                read(73,*) !'# Total myocardial oxygen demand during the beat before the last complete beat:'
                read(73,*) this%MVO2computedAtEndOfBeatPreviousToPreviousBeat(ii)
                read(73,*) !'# Latest myocardial hunger delta ("small h(t)")'
                read(73,*) this%hungerDelta(ii)
              enddo
              read(73,*) !'# Timesteps since the start of the last heart period (for MVO2 calculation)'
              read(73,*) hrt%timestepsSinceLastHeartPeriodReset
              write(73,*) !'# Are we starting on the first dt of a new heart-beat? 1=yes, 0=no.'
              write(73,*) hrt%newBeatJustStarted
              close(73)

              ! Save O2supplyDemandHistory as a binary file:
              ! (easier just to have this done in the multidomain module, and called here)
              call this%loadO2supplyDemandHistory()

              do ii = 1,this%numberOfControlledCoronaryModelSurfaces
                this%oneOverR(ii) = 1d0/(this%rp(ii) + this%rd(ii) + this%ra(ii))
              enddo
         endif
      end subroutine initialiseCoronaryControlModel


      subroutine updateMVO2(this,lstep_passedIn)
         implicit none
         ! include 'mpif.h'
         
         class(controlledCoronaryModel) :: this
         integer, intent(in) :: lstep_passedIn

         real*8, allocatable :: plv_overLastBeat(:)
         real*8, allocatable :: vlv_overLastBeat(:)
         ! real*8 :: tempPreviousMVO2
         integer :: ii
         integer :: isystemic
         integer :: rank
         integer :: ierr

         ! real*8 :: testArray1(4) !\todo remove this testcode
         ! real*8 :: testArray2(4)
         ! real*8 :: testResult(1)

         ! testArray1(1) = 2
         ! testArray1(2) = 1
         ! testArray1(3) = 1
         ! testArray1(4) = 2

         ! testArray2(1) = 2
         ! testArray2(2) = 2
         ! testArray2(3) = 1
         ! testArray2(4) = 1

         ! testResult = computePV_Area(testArray1,testArray2)
         ! write(*,*) 'area:', testResult
         ! stop

         do ii=1, this%numberOfControlledCoronaryModelSurfaces
            this%MVO2previousDt(ii) = this%MVO2(ii)
         enddo

         if (hrt%newBeatJustStarted .eq. int(1)) then !if this is the first time-step after a heart period reset
            allocate(plv_overLastBeat(hrt%timestepsSinceLastHeartPeriodReset))
            allocate(vlv_overLastBeat(hrt%timestepsSinceLastHeartPeriodReset))

            plv_overLastBeat = hrt%plv_hist(lstep_passedIn-hrt%timestepsSinceLastHeartPeriodReset:lstep_passedIn-1) !\todo ensure the "-1"s belong here.
            vlv_overLastBeat = hrt%vlv_hist(lstep_passedIn-hrt%timestepsSinceLastHeartPeriodReset:lstep_passedIn-1) !\todo ensure the "-1"s belong here.

            do ii=1, this%numberOfControlledCoronaryModelSurfaces
               this%MVO2computedAtEndOfBeatPreviousToPreviousBeat(ii) = this%MVO2computedAtEndOfPreviousBeat(ii)
               ! units for MVO2 are Joules: (so must convert volume to m^3, and pressure to Pa in the next line)
               this%MVO2computedAtEndOfPreviousBeat(ii) = computePV_Area(vlv_overLastBeat/real(1E6,8),plv_overLastBeat/real(10.0,8))
               ! Conversion factors on the next line of code:
               ! *3 - the PVA to total energy requirements ratio
               ! /20 is to convert energy to ml O2
               ! Details of the *3 and the /20 below can be found in Kameyama et
               ! al. 1992, "Energy conversion efficiency in human left ventricle".
               ! See Table 2 for the *3 (PVA/MVO2 means of 33 and 36);
               ! see the text on page 990 for the /20.
               ! /tperiod converts to MVO2 per second rather than per period
               this%MVO2computedAtEndOfPreviousBeat(ii) = this%MVO2computedAtEndOfPreviousBeat(ii) * this%proportionOfMyocardiumPerfusedByOutlet(ii) * real(3.0,8) / real(20.0,8) / hrt%period !\todo consider whether this period should just be the current (continuously-changing, dynamic) period, or (one-over) the duration of the previous beat, or whether it doesn't matter.
               ! The value of MVO2 we use lags the PV-loop calculations by one beat. This is because we wish to
               ! linearly interpolate MVO2 between two known values of MVO2, so we use the last two /complete/ beats.
               ! This introduces an additional one-beat delay into the system. \todo consider alternatives, but this is reasonable for now.
               this%deltaMVO2PerTimestepOverPreviousBeat(ii) = (this%MVO2computedAtEndOfPreviousBeat(ii) - this%MVO2computedAtEndOfBeatPreviousToPreviousBeat(ii)) / hrt%timestepsSinceLastHeartPeriodReset
               ! if (abs(this%deltaMVO2PerTimestepOverPreviousBeat(ii))*hrt%timestepsSinceLastHeartPeriodReset .gt. this%MVO2(ii)) then
               !    this%deltaMVO2PerTimestepOverPreviousBeat(ii) = 0
               !    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
               !    if (rank .eq. int(0)) then
               !       write(*,*) '(WW) Warning: MVO2 was to change too much over this heartbeat. Setting deltaMVO2PerTimestepOverPreviousBeat to zero.'
               !       write(*,*) '(WW) If this happens during the first few beats, you probably initialised MVO2 to an unphysiological value,'
               !       write(*,*) '(WW) and the problem will resolve itself after a few beats. After that, see coronary_restart.dat for a more physiological MVO2.'
               !    end if
               ! end if
            enddo
            hrt%timestepsSinceLastHeartPeriodReset = 0
         else !this is any other time-step within the heart cycle
            ! continue counting the number of steps taken this period:
            hrt%timestepsSinceLastHeartPeriodReset = hrt%timestepsSinceLastHeartPeriodReset + 1
            ! We ramp the MVO2 each time-step by assuming that the change per-dt on this beat will be the same as on the previous beat.
            ! \todo (low priority) - consider alternate approaches based on the current period here - and whether doing something different would add unnecessary complexity
            do ii=1, this%numberOfControlledCoronaryModelSurfaces
               this%MVO2(ii) = this%MVO2(ii) + this%deltaMVO2PerTimestepOverPreviousBeat(ii)
               !\todo probably want to set some minimum/maximum MVO2 values, or at least generate a warning if MVO2 becomes unphysiological.
            end do
         end if

         ! write(*,*) 'MVO2 at outlet 1 is: ', this%MVO2(1)

      end subroutine updateMVO2


      subroutine computeHungerAndSetHistories(this,lstep_passedIn)
         ! This version of the model does not track O2 mass; it just adjusts flow
         ! according to whether the historial O2 supply was too high or too low.
         ! This model also has a "limit on memory": O2 supply discrepancy from
         ! more than a (user-set) time ago is "forgotten", which may not be physiological
         ! because oxygen debts need to be (over-)repaid by reactive hyperaemia.
         use calcFlowPressure, only: FlowHist
         implicit none

         class(controlledCoronaryModel) :: this
         integer, intent(in) :: lstep_passedIn
         integer :: ii, currentSurface

         ! update the moving average of the O2 supply - demand difference
         do ii=1, this%numberOfControlledCoronaryModelSurfaces
            ! this is computed as a "rolling mean" - subtract out the value that just fell out of the start of the averaging window, and add the new value
            currentSurface = this%localToGlobalSurfaceIndexMap(ii) ! this is needed where an array (eg (and currently only) FlowHist) is indexed by global surface number, rather than coronary surface number

            !\todo are you reading FlowHist from the right place?
            this%currentMyocardialHungerSignal(ii) = this%currentMyocardialHungerSignal(ii) - this%O2supplyDemandHistory(lstep_passedIn,ii)/this%O2supplyDemandHistoryWindowLength_timesteps(ii) + (this%MVO2(ii) - FlowHist(lstep_passedIn+1,currentSurface) * this%arterial_O2_volume_proportion)/this%O2supplyDemandHistoryWindowLength_timesteps(ii)

            ! record the O2 supply-demand difference for this timestep:
            this%O2supplyDemandHistory(lstep_passedIn + this%O2supplyDemandHistoryWindowLength_timesteps(ii),ii) = this%MVO2(ii) - FlowHist(lstep_passedIn+1,currentSurface)*this%arterial_O2_volume_proportion !\todo removed brackets here - pretty sure this is correct now, but be aware!
         enddo
      end subroutine computeHungerAndSetHistories

      subroutine computeHungerAndSetHistories_negativeHungerCapped(this,lstep_passedIn)
         ! This version of the model does oxygen mass delivery vs oxygen mass consumption tracking.
         ! It allows the hunger to get arbitrarily large, but caps the amount of O2 the myocardium
         ! can store up for a rainy day (capped negative hunger).
         use calcFlowPressure, only: FlowHist
         implicit none

         class(controlledCoronaryModel) :: this
         integer :: ii, currentSurface
         integer, intent(in) :: lstep_passedIn
         real*8 :: maxNegativeHunger

         maxNegativeHunger = -0.025 !\todo input this from config files

         ! update the moving average of the O2 supply - demand difference
         do ii=1, this%numberOfControlledCoronaryModelSurfaces
            ! this is computed as a "rolling mean" - subtract out the value that just fell out of the start of the averaging window, and add the new value
            currentSurface = this%localToGlobalSurfaceIndexMap(ii) ! this is needed where an array (eg (and currently only) FlowHist) is indexed by global surface number, rather than coronary surface number

            !\todo are you reading FlowHist from the right place?
            this%hungerDelta(ii) = (this%MVO2(ii) - FlowHist(lstep_passedIn+1,currentSurface) * this%arterial_O2_volume_proportion)
             ! this%hungerDelta(ii) = max(maxNegativeHunger, hungerDelta(ii)) !\todo remove

            ! record the current total level of myocardial O2 hunger:
            this%O2supplyDemandHistory(lstep_passedIn + this%O2supplyDemandHistoryWindowLength_timesteps(ii),ii) = max(maxNegativeHunger, this%O2supplyDemandHistory(lstep_passedIn + this%O2supplyDemandHistoryWindowLength_timesteps(ii) - 1,ii) + this%hungerDelta(ii)*delt)

            this%currentMyocardialHungerSignal(ii) = this%O2supplyDemandHistory(lstep_passedIn+1,ii)
         enddo
      end subroutine computeHungerAndSetHistories_negativeHungerCapped


      subroutine  updateFeedforwardAndFeedbackInCoronaryControlModel(this,lstep_passedIn)
         use calcFlowPressure, only: FlowHist
         implicit none
         class(controlledCoronaryModel) :: this
         integer, intent(in) :: lstep_passedIn
         integer :: ii, currentSurface

         call this%updateMVO2(lstep_passedIn)
         ! Get the latest MVO2 gradient, for use in the feedforward control:
         do ii=1, this%numberOfControlledCoronaryModelSurfaces
            this%deltaMVO2(ii) = (this%MVO2(ii) - this%MVO2previousDt(ii))/delt
         enddo

         call this%computeHungerAndSetHistories(lstep_passedIn)

         call this%computeResistanceChanges(lstep_passedIn)

         ! this is to handle the case where oneOverR goes negative:
         do ii=1, this%numberOfControlledCoronaryModelSurfaces
            if (this%oneOverR(ii) < 0) then
#if EXTRA_CONSOLE_OUTPUT == 1
               write(*,*) 'oneOverR went negative, setting to max resistance resistance, was ', 1/this%oneOverR(ii), ' vs ', this%R_max(ii) + this%rp(ii) + this%ra(ii)
#endif
               this%oneOverR(ii) = 1/(this%R_max(ii) + this%rp(ii) + this%ra(ii))
            end if
#if EXTRA_CONSOLE_OUTPUT == 1
            write(*,*) 'oneOverR ', 'ii: ', this%oneOverR(ii)
#endif
         end do


         ! Limit the total resistance to be between two (physiologically-set) values
         do ii=1, this%numberOfControlledCoronaryModelSurfaces
            if (1/this%oneOverR(ii) < this%R_min(ii) + this%rp(ii) + this%ra(ii)) then
#if EXTRA_CONSOLE_OUTPUT == 1
               write(*,*) 'Hit minimum resistance, was ', 1/this%oneOverR(ii), ' vs ', this%R_min(ii) + this%rp(ii) + this%ra(ii)
#endif
               this%oneOverR(ii) = 1/(this%R_min(ii) + this%rp(ii) + this%ra(ii))
            elseif (1/this%oneOverR(ii) > this%R_max(ii) + this%rp(ii) + this%ra(ii)) then
#if EXTRA_CONSOLE_OUTPUT == 1
               write(*,*) 'Hit maximum resistance, was ', 1/this%oneOverR(ii), ' vs ', this%R_max(ii) + this%rp(ii) + this%ra(ii)
#endif
               this%oneOverR(ii) = 1/(this%R_max(ii) + this%rp(ii) + this%ra(ii))
            endif
            ! Set rd according to the total resistance that we have found:
            this%rd(ii) = 1/this%oneOverR(ii) - this%rp(ii) - this%ra(ii)

            ! Do the alpha-feedforward control of the resistance
            this%rp(ii) = 1/(1/this%rp(ii) + delt * this%alphaFeedforwardGain(ii)*this%deltaMVO2(ii)/this%arterial_O2_volume_proportion / (100*1333.2237)) ! \todo Ensure scaling is correct here

#if EXTRA_CONSOLE_OUTPUT == 1
            write(*,*) 'The current values of rp and rd are:', this%rp(ii), this%rd(ii), 'for outlet ', ii
#endif
         enddo


      end subroutine updateFeedforwardAndFeedbackInCoronaryControlModel


      subroutine computeResistanceChanges_undamped(this, lstep_passedIn)
         ! This version of the resistance control ODE is the undamped, origina
         ! model from the CMBE 2013 conference presentation (see the published abstract)
         implicit none

         class(controlledCoronaryModel) :: this
         integer, intent(in) :: lstep_passedIn
         integer :: ii


         do ii=1, this%numberOfControlledCoronaryModelSurfaces
            this%oneOverR(ii) = this%oneOverR(ii) + delt*(this%feedbackGain(ii)*this%currentMyocardialHungerSignal(ii) + this%betaFeedforwardGain(ii)*this%deltaMVO2(ii)/this%arterial_O2_volume_proportion) / (100*1333.2237) !\todo sort out units here!
         enddo

      end subroutine computeResistanceChanges_undamped


      subroutine computeResistanceChanges_harmonic(this, lstep_passedIn)
         ! This version of the resistance control ODE is the spring-style model
         ! from the paper.
         implicit none

         class(controlledCoronaryModel) :: this
         integer, intent(in) :: lstep_passedIn
         integer :: ii


         do ii=1, this%numberOfControlledCoronaryModelSurfaces
            this%oneOverR(ii) = this%oneOverR(ii) + ( delt*this%feedbackGain(ii)*this%currentMyocardialHungerSignal(ii) +      &
                                                      delt*this%dampingCoefficient(ii)*this%hungerDelta(ii) +                  &
                                                      this%betaFeedforwardGain(ii)*this%deltaMVO2(ii)*delt                  )  &
                                             / (100*1333.2237) / this%arterial_O2_volume_proportion!\todo sort out units here!
         enddo

      end subroutine computeResistanceChanges_harmonic


      ! This function takes in a polygon as a set of points (x-coords coords in polygonNodesX, y-coords in polygonNodesY)
      ! and returns the area enclosed by the convex hull consisting of the polygoin points PLUS THE ORIGIN.
      ! The sides of the polygon should not intersect. (0,0) should not be included.
      ! polygonNodesX should be the ventricular volumes,
      ! polygonNodesY should be the ventricular pressures.
      ! This is a MODIFICATION of a standard formula; see eg. http://www.mathopenref.com/coordpolygonarea.html
      ! The modification is to include the origin in the convex hull.
      ! Note that the points MUST be numbered anticlockwise!
      function computePV_Area(polygonNodesX,polygonNodesY)
         implicit none

         real*8 :: computePV_Area
         real*8, intent(in) :: polygonNodesX(:)
         real*8, intent(in) :: polygonNodesY(:)

         real*8 :: accumulator
         real*8 :: potentialNextAddition
         integer :: ii, numberOfPolygonVertices

         numberOfPolygonVertices = size(polygonNodesX) ! this array must not be padded with zeros (or with anything else!)

         ! Form the sum at the core of the formula:
         accumulator=0
         ii=1
         do while (ii<numberOfPolygonVertices)
            potentialNextAddition = polygonNodesX(ii)*polygonNodesY(ii+1) - polygonNodesX(ii+1)*polygonNodesY(ii)
            ! It is this if-guard that allows the inclusion of the area to the origin.
            if (potentialNextAddition .gt. int(0)) then
               accumulator = accumulator + potentialNextAddition
            endif
            ii = ii + 1
         end do
         ! Add the last sum term, which doesn't work nicely with the above loop:
         potentialNextAddition = polygonNodesX(numberOfPolygonVertices)*polygonNodesY(1) - polygonNodesX(1)*polygonNodesY(numberOfPolygonVertices)
         if (potentialNextAddition .gt. int(0)) then
            accumulator = accumulator + potentialNextAddition
         endif

         ! Just need half now, to complete the formula:
         ! (the absolute value here could be taken if we didn't want the origin included, which would relax the constraint that the points must be anticlockwise.)
         computePV_Area = accumulator/2

      end function computePV_Area


!     *************** END CORONARY MODEL CODE ***************


!     *************** BEGIN NETLIST ARBITRARY-DESIGN LPN MODEL CODE ***************


      function netlistLPNConstructor(numNetlistSrfs,nSrfListNetlist,numCalcSrfs,nsrflistCalc) result(a)
         implicit none
         integer, intent(in) :: numNetlistSrfs
         integer, intent(in) :: nSrfListNetlist(0:maxsurf)
         integer, intent(in) :: numCalcSrfs
         integer, intent(in) :: nsrflistCalc(0:maxsurf)
         type(netlistLPN) :: a
         call a%initialiseNetlistLPN(numNetlistSrfs,nSrfListNetlist,numCalcSrfs,nsrflistCalc)
         return
      end function netlistLPNConstructor

      
      subroutine initialiseNetlistLPN(this,numNetlistSrfs,nSrfListNetlist,numCalcSrfs,nsrflistCalc)

         implicit none

         include "mpif.h"

         integer, intent(in) :: numNetlistSrfs
         integer, intent(in) :: nSrfListNetlist(0:maxsurf)
         integer, intent(in) :: numCalcSrfs
         integer, intent(in) :: nsrflistCalc(0:maxsurf)
         class(netlistLPN) :: this
         integer :: ii
         integer :: jj

         integer :: rank
         integer :: err

         integer :: expectedNumberOfLines
         integer :: nextReadComponentCount
         integer :: newval
         integer :: numberOfLines
         integer :: iostatus
         logical :: file_exists


         this%classNameString = 'netlistLPN'

         ! set surface numbers and lists
         this%numberOfLPNSurfaces = numNetlistSrfs
         allocate(this%numberOfComponents(this%numberOfLPNSurfaces))
         allocate(this%P_a(this%numberOfLPNSurfaces))
         allocate(this%localToGlobalSurfaceIndexMap(this%numberOfLPNSurfaces))
         this%surfnum = numNetlistSrfs !\todo check this and line above both necessary! (this line is needed)
         allocate(this%surfids(numNetlistSrfs))
         this%surfids(1:numNetlistSrfs) = nSrfListNetlist(1:numNetlistSrfs)
         allocate(this%flow_n(this%surfnum))
         allocate(this%flow_n1(this%surfnum))
         allocate(this%pressure_n(this%surfnum))
         allocate(this%implicitcoeff(this%surfnum,2))
         allocate(this%implicitcoeff_n1(this%surfnum,2))

         allocate(this%flowpntr(this%numberOfLPNSurfaces))

         allocate(this%numberOfPrescribedPressures(this%numberOfLPNSurfaces))
         allocate(this%numberOfPrescribedFlows(this%numberOfLPNSurfaces))
         allocate(this%numberOfPrescribedPressuresAndFlows(this%numberOfLPNSurfaces))

         allocate(this%numberOfPressureNodes(this%numberOfLPNSurfaces))

         call MPI_COMM_RANK(MPI_COMM_WORLD,rank,err)
         
         inquire(file='netlist_surfaces_restart.dat', exist = file_exists)
         if (file_exists) then
            open(unit=73,file='netlist_surfaces_restart.dat',status='old')
            if (rank .eq. int(0)) then
               write(*,*) '====> Loading netlist restart status from a previous simulation (netlist_surfaces_restart.dat)'
            end if
         else
            inquire(file="netlist_surfaces.dat", exist=file_exists)
            if (file_exists) then
               open(unit=73,file='netlist_surfaces.dat',status='old')
            else
               if (rank .eq. int(0)) then
                  write(*,*) 'ERROR: netlist_surfaces.dat not found. Exiting.'
               end if
               stop
            end if
         end if

         read(73,*) ! Read & ignore comment line: "# List of components in a format similar to that for netlist. All comment lines must be present, but the comment content is itself irrelevant."
         
         read(73,*) ! Read & ignore comment line: "# Maximum number of components that any of the netlist boundary LPNs have:"
         read(73,*) this%maxComponents
         allocate(this%circuitData(this%maxComponents, 3, this%numberOfLPNSurfaces))
         allocate(this%circuitData_componentTypes(this%maxComponents,this%numberOfLPNSurfaces))
         
         read(73,*) ! Read & ignore comment line: "# Maximum number of prescribed pressure nodes that the LPNs have:"
         read(73,*) this%maxPrescribedPressureNodes
         allocate(this%listOfPrescribedPressures(this%maxPrescribedPressureNodes,this%numberOfLPNSurfaces))
         allocate(this%valueOfPrescribedPressures(this%maxPrescribedPressureNodes,this%numberOfLPNSurfaces))
         allocate(this%typeOfPrescribedPressures(this%maxPrescribedPressureNodes,this%numberOfLPNSurfaces))

         read(73,*) ! Read & ignore comment line: "# Maximum number of prescribed flows that the LPNs have:"
         read(73,*) this%maxPrescribedFlows
         allocate(this%listOfPrescribedFlows(this%maxPrescribedFlows,this%numberOfLPNSurfaces))
         allocate(this%valueOfPrescribedFlows(this%maxPrescribedFlows,this%numberOfLPNSurfaces))
         allocate(this%typeOfPrescribedFlows(this%maxPrescribedFlows,this%numberOfLPNSurfaces))

         read(73,*) ! Read & ignore comment line: "# Maximum number of nodes that any LPN has:"
         read(73,*) this%maxPressureNodes
         allocate(this%pressuresInLPN(this%maxPressureNodes,this%numberOfLPNSurfaces))   ! Pressure at each LPN node, using the same node indexing as in the netlist
         allocate(this%historyPressuresInLPN(this%maxPressureNodes,this%numberOfLPNSurfaces))

         do ii=1, this%numberOfLPNSurfaces
            read(73,*) ! Read & ignore comment line: "### Begin first netlist boundary condition model"
            read(73,*) ! Read & ignore comment line: "# Number of Components"
            read(73,*) this%numberOfComponents(ii)
            do jj=1, this%numberOfComponents(ii)
               read(73,*) ! Read & ignore comment line: "# Component jj type:"
               read(73,*) this%circuitData_componentTypes(jj,ii)
               read(73,*) ! Read & ignore comment line: "# Component jj details (start-node index, end-node index, associated parameter (resistance for resistors, capacitance for capacitors):"
               read(73,*) this%circuitData(jj,1,ii)
               read(73,*) this%circuitData(jj,2,ii)
               read(73,*) this%circuitData(jj,3,ii)
            end do
            read(73,*) ! Read & ignore comment line: "# Number of prescribed pressure nodes:"
            read(73,*) this%numberOfPrescribedPressures(ii)
            if (this%numberOfPrescribedPressures(ii) .gt. int(0)) then
               read(73,*) ! Read & ignore comment line: "# Indices of nodes with prescribed pressures:" !\todo make this work with not-three columns
               do jj=1, this%numberOfPrescribedPressures(ii)
                  read(73,*) this%listOfPrescribedPressures(jj,ii)
               end do
               read(73,*) ! Read & ignore comment line: "# Prescribed pressure values / scalings (dependent on types, given by next component):"
               do jj=1, this%numberOfPrescribedPressures(ii)
                  read(73,*) this%valueOfPrescribedPressures(jj,ii)
               end do
               read(73,*) ! Read & ignore comment line: "# Prescribed pressure types (f=fixed to value given in previous line, l=left ventricular pressure, scaled by value given in previous line):"
               do jj=1, this%numberOfPrescribedPressures(ii)
                  read(73,*) this%typeOfPrescribedPressures(jj,ii)
               end do
            end if
            read(73,*) ! Read & ignore comment line: "# Number of prescribed flows:"
            read(73,*) this%numberOfPrescribedFlows(ii)
            if (this%numberOfPrescribedFlows(ii) .gt. int(0)) then
               read(73,*) ! Read & ignore comment line: "# Indices of components with prescribed flows"
               do jj=1, this%numberOfPrescribedFlows(ii)
                  read(73,*) this%listOfPrescribedFlows(jj,ii)
               end do
               read(73,*) ! Read & ignore comment line: "# Values of prescribed flows (3D interface set to -1; this value is irrelevant and unused by the code):"
               do jj=1, this%numberOfPrescribedFlows(ii)
                  read(73,*) this%valueOfPrescribedFlows(jj,ii)
               end do
               read(73,*) ! Read & ignore comment line: "# Types of prescribed flows (t=threeD domain interface)"
               do jj=1, this%numberOfPrescribedFlows(ii)
                  read(73,*) this%typeOfPrescribedFlows(jj,ii)
               end do
            end if
            read(73,*) ! Read & ignore comment line: "# Number of pressure nodes (including everything- 3D interface, zero-pressure points, internal nodes, etc.):"
            read(73,*) this%numberOfPressureNodes(ii)
            read(73,*) ! Read & ignore comment line: "# Initial pressures at the pressure nodes:"
            do jj=1, this%numberOfPressureNodes(ii)
               read(73,*) this%pressuresInLPN(jj,ii) !\todo this will break restarts - you must load these from an alternative place in that case!
            end do
            read(73,*) ! Read & ignore comment line: "### End first netlist boundary condition model"
          enddo
         close(73)

         this%numberOfPrescribedPressuresAndFlows = this%numberOfPrescribedPressures + this%numberOfPrescribedFlows ! Just the sum of the previous two declared integers

         !\todo make these dynamic
         ! allocate(this%circuitData(5,3))
         ! allocate(this%circuitData_componentTypes(5))

         ! this%circuitData = reshape((/1.0d0,2.0d0,2.0d0,4.0d0,4.0d0,  2.0d0,3.0d0,4.0d0,5.0d0,6.0d0,  1.286208d5,4.5d-7,1.286208d5,2.7d-6,7.5d5/), (/5,3/)) !column-major content, then the array shape
         ! this%circuitData_componentTypes = (/'r','c','r','c','r'/) ! the data in here will be the stripped first column of the netilst, identifying each line of circuitData as being r=resistor, c=capacitor, etc.

         ! this%circuitData(1,3) = 1.286208000000000e5
         ! this%circuitData(2,3) = real(4.5000000000e-7,8)

         ! open(unit=73,file='netlist.dat',status='old')
         ! read(73,*) this%circuitData(1,3)
         ! read(73,*) this%circuitData(2,3)
         ! read(73,*) this%circuitData(3,3)
         ! read(73,*) this%circuitData(4,3)
         ! read(73,*) this%circuitData(5,3)
         ! close(73)

         ! this%circuitData(3,3) = 1.286208000000000e5
         ! this%circuitData(4,3) = 2.7000000000e-6
         ! this%circuitData(5,3) = real(7.5000000000e5,8)

         ! this%numberOfPrescribedPressures = 3 ! \todo set this dynamically
         ! this%numberOfPrescribedFlows = 1 ! \todo set this dynamically (when inductor added) - the 1 is for the inflow face (will need to hard-code this for now!)
         
         ! allocate(this%listOfPrescribedFlows(this%numberOfPrescribedFlows))

         ! this%listOfPrescribedPressures = (/3,5,6/) !\todo make dynamic
         ! this%listOfPrescribedFlows = (/1/) !\todo needed? \todo make dynamic

         !\todo make these an input
         !\note that for non-fixed values (non-'f'), such as 'l' for LV pressure, the valueOfPrescribedPressure gives a scaling (so here we're using 1x LV pressure for the 2nd prescribed node)
         ! this%valueOfPrescribedPressures = (/0d0,1d0,0d0/)
         ! this%typeOfPrescribedPressures = (/'f','l','f'/)
         ! this%typeOfPrescribedFlows = (/'3'/) !'3' for 3D domain

         ! this%numberOfComponents = 5 ! \todo set this dynamically
         ! this%numberOfPressureNodes = 6 ! \todo set this dynamically
         allocate(this%flowsInLPN(this%maxComponents,this%numberOfLPNSurfaces))                           ! Flow through each component in the LPN, in the order they appear in the netlist

         this%flowsInLPN = 0d0

         call this%getMapOfPressHistoriesToCorrectPressNodes() !initialises this%numberOfHistoryPressures
         call this%getMapOfFlowHistoriesToCorrectComponents() !initialises this%numberOfHistoryFlows

         allocate(this%systemSize(this%numberOfLPNSurfaces))

         this%systemSize = this%numberOfPressureNodes + this%numberOfHistoryPressures + this%numberOfComponents + this%numberOfHistoryFlows

         this%maxSystemSize = maxval(this%systemSize)

         allocate(this%systemMatrix(this%maxSystemSize,this%maxSystemSize,this%numberOfLPNSurfaces))
         this%systemMatrix = 0d0
         allocate(this%RHS(this%maxSystemSize,this%numberOfLPNSurfaces))
         allocate(this%columnIndexOf3DInterfaceFlowInLinearSystem(this%numberOfLPNSurfaces))
         allocate(this%solutionVector(this%maxSystemSize,this%numberOfLPNSurfaces))

         allocate(this%inverseOfSystemMatrix(this%maxSystemSize,this%maxSystemSize,this%numberOfLPNSurfaces))

         this%maxColumnMapSize = maxval(this%numberOfHistoryPressures + this%numberOfHistoryFlows + this%numberOfPrescribedPressures + this%numberOfPrescribedFlows)
         allocate(this%columnMap(this%maxColumnMapSize,this%numberOfLPNSurfaces))

         call this%getListOfNodesWithMultipleIncidentCurrents()

         this%surfnum = this%numberOfLPNSurfaces !\todo is this necessary?

!        Populate the map which converts LPN surface indices (local) to surface indices (global) - for FlowHist etc.
         do ii=1, this%numberOfLPNSurfaces
            do jj=1, numCalcSrfs
               if(this%surfids(ii).eq.nsrflistCalc(jj)) then
                  this%localToGlobalSurfaceIndexMap(ii) = jj
               endif
            enddo
         enddo

      end subroutine initialiseNetlistLPN

      subroutine writeRestart_netlistLPNSurfaces(this)

         implicit none

         include "mpif.h"

         class(netlistLPN) :: this
         integer :: ii
         integer :: jj
         integer :: rank
         integer :: err

         call MPI_COMM_RANK(MPI_COMM_WORLD,rank,err)
         if (err .ne. int(0)) then
            write(*,*) 'MPI Error in writeRestart_netlistLPNSurfaces. Code: ', err
            stop
         end if

         if (rank .eq. int(0)) then
            open(unit=73,file='netlist_surfaces_restart.dat',status='replace')
            write(73,*) '# Netlist LPN restart file. See netlist_surfaces.dat.'
            write(73,*) '# Max components:'
            write(73,*) this%maxComponents
                     
            write(73,*) '# Max prescribed pressure nodes:'
            write(73,*) this%maxPrescribedPressureNodes

            write(73,*) '# Max prescribed flows:'
            write(73,*) this%maxPrescribedFlows

            write(73,*) '# Max pressure nodes:'
            write(73,*) this%maxPressureNodes

            do ii=1, this%numberOfLPNSurfaces
               write(73,*) '### Begin model for surface:', ii
               write(73,*) '# Number of Components'
               write(73,*) this%numberOfComponents(ii)
               do jj=1, this%numberOfComponents(ii)
                  write(73,*) '# Component type:', jj
                  write(73,*) this%circuitData_componentTypes(jj,ii)
                  write(73,*) '# Component details:'
                  write(73,*) this%circuitData(jj,1,ii)
                  write(73,*) this%circuitData(jj,2,ii)
                  write(73,*) this%circuitData(jj,3,ii)
               end do
               write(73,*) '# Number of prescribed pressure nodes:'
               write(73,*) this%numberOfPrescribedPressures(ii)
               write(73,*) '# Indices of nodes with prescribed pressures:'
               do jj=1, this%numberOfPrescribedPressures(ii)
                  write(73,*) this%listOfPrescribedPressures(jj,ii)
               end do
               write(73,*) '# Prescribed pressure values OR scalings:'
               do jj=1, this%numberOfPrescribedPressures(ii)
                  write(73,*) this%valueOfPrescribedPressures(jj,ii)
               end do
               write(73,*) '# Prescribed pressure types:'
               do jj=1, this%numberOfPrescribedPressures(ii)
                  write(73,*) this%typeOfPrescribedPressures(jj,ii)
               end do
               write(73,*) '# Number of prescribed flows:'
               write(73,*) this%numberOfPrescribedFlows(ii)
               write(73,*) '# Indices of components with prescribed flows'
               do jj=1, this%numberOfPrescribedFlows(ii)
                  write(73,*) this%listOfPrescribedFlows(jj,ii)
               end do
               write(73,*) '# Values of prescribed flows: (-1.0 @ 3D)'
               do jj=1, this%numberOfPrescribedFlows(ii)
                  write(73,*) this%valueOfPrescribedFlows(jj,ii)
               end do
               write(73,*) '# Types of prescribed flows'
               do jj=1, this%numberOfPrescribedFlows(ii)
                  write(73,*) this%typeOfPrescribedFlows(jj,ii)
               end do
               write(73,*) '# Number of pressure nodes:'
               write(73,*) this%numberOfPressureNodes(ii)
               write(73,*) '# Initial pressures at the pressure nodes:'
               do jj=1, this%numberOfPressureNodes(ii)
                  write(73,*) this%pressuresInLPN(jj,ii) !\todo this will break restarts - you must load these from an alternative place in that case!
               end do
               write(73,*) '### End model for one surface'
             enddo
            close(73)
         end if

      end subroutine writeRestart_netlistLPNSurfaces


      subroutine writeRestart_controlledCoronaryModel(this,lstep_passedIn)
         ! Save data necessary to restart the LPN models at coronary surfaces:
         ! this is done using type-bound public "get" procedures which allow read-only
         ! access to the private variables.
         implicit none

         include 'mpif.h'

         class(controlledCoronaryModel), intent(in) :: this
         integer, intent(in) :: lstep_passedIn
         integer :: ii
         integer :: ierr
         integer :: rank

         call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
         if (rank .eq. int(0)) then
            ! Save O2supplyDemandHistory as a binary file:
            ! (easier just to have this done in the multidomain module, and called here)
            call controlledCoronarySurfaces%writeO2supplyDemandHistory()

            ! Output some debug info on the dynamic cor. resistances
            call controlledCoronarySurfaces%writeDebugHistories(lstep_passedIn)

            ! save the current resistances, etc. at each coronary outlet 
            open(unit=73,file='coronary_restart.dat',status='replace')
            write(73,*) '# ~=Controlled Coronary Model Restart File=~'
            do ii = 1,this%getNumberOfControlledCoronaryModelSurfaces()
               write(73,*) '# Coronary outlet index number:', ii
               write(73,*) '# Resistance rp:'
               write(73,*) this%getRpByIndex(ii)
               write(73,*) '# Resistance rd:'
               write(73,*) this%getRdByIndex(ii)
               write(73,*) '# Coronary LPN internal pressure P_1:'
               write(73,*) this%getP_1ByIndex(ii)
               write(73,*) '# Coronary LPN internal pressure P_2:'
               write(73,*) this%getP_2ByIndex(ii)
               write(73,*) '# currentMyocardialHungerSignal:'
               write(73,*) this%getcurrentMeanMyocardialOxygenHungerByIndex(ii)
               write(73,*) '# MVO2previousDt:'
               write(73,*) this%getMVO2previousByIndex(ii)
               write(73,*) '# MVO2'
               write(73,*) this%MVO2(ii)
               write(73,*) '# P_a (LPN pressure at outlet itself):'
               write(73,*) this%getP_aByIndex(ii)
               write(73,*) '# deltaMVO2PerTimestepOverPreviousBeat (used for MVO2 ramping on present beat)'
               write(73,*) this%deltaMVO2PerTimestepOverPreviousBeat(ii)
               write(73,*) '# Total myocardial oxygen demand during the last complete beat:'
               write(73,*) this%MVO2computedAtEndOfPreviousBeat(ii)
               write(73,*) '# Total myocardial oxygen demand during the beat before the last complete beat:'
               write(73,*) this%MVO2computedAtEndOfBeatPreviousToPreviousBeat(ii)
               write(73,*) '# Latest myocardial hunger delta ("small h(t)")'
               write(73,*) this%hungerDelta(ii)
            enddo
            write(73,*) '# Timesteps since the start of the last heart period (for MVO2 calculation)'
            write(73,*) hrt%timestepsSinceLastHeartPeriodReset
            write(73,*) '# Are we starting on the first dt of a new heart-beat? 1=yes, 0=no.'
            write(73,*) hrt%newBeatJustStarted
            close(73)
         end if

      end subroutine writeRestart_controlledCoronaryModel


      subroutine generateLinearSystemFromPrescribedCircuit(this,alfi_delt)
         ! This subroutine assembles the system of (time-discretised) linear algebraic equations for the LPN.
         
          implicit none

          integer :: ll, mm
          integer :: rowsDoneSoFar
          integer :: tempIndexingShift
          integer :: tempUnknownVariableIndexWithinLinearSystem
          integer :: kk

          integer :: ierr

          class(netlistLPN) :: this
          real*8, intent(in) :: alfi_delt

          do kk=1, this%numberOfLPNSurfaces
             do ll=1, this%numberOfComponents(kk)
               if (this%circuitData_componentTypes(ll,kk) == 'r') then
                  ! insert resistor relationship into equation system
                  this%systemMatrix(ll,int(this%circuitData(ll,1,kk)),kk) = 1.0d0
                  this%systemMatrix(ll,int(this%circuitData(ll,2,kk)),kk) = -1.0d0
                  this%systemMatrix(ll,ll+this%numberOfPressureNodes(kk)+this%numberOfHistoryPressures(kk),kk) = -this%circuitData(ll,3,kk)
               elseif (this%circuitData_componentTypes(ll,kk) == 'c') then
                  ! insert capacitor relationship into equation system
                  this%systemMatrix(ll,int(this%circuitData(ll,1,kk)),kk) = 1.0d0
                  this%systemMatrix(ll,int(this%circuitData(ll,2,kk)),kk) = -1.0d0
                  this%systemMatrix(ll,ll+this%numberOfPressureNodes(kk)+this%numberOfHistoryPressures(kk),kk) = -alfi_delt/this%circuitData(ll,3,kk) ! note that this needs to change for alfi \todo
                  this%systemMatrix(ll,this%nodeIndexToPressureHistoryNodeOrderingMap(int(this%circuitData(ll,1,kk)),kk) + this%numberOfPressureNodes(kk),kk) = -1.0d0
                  this%systemMatrix(ll,this%nodeIndexToPressureHistoryNodeOrderingMap(int(this%circuitData(ll,2,kk)),kk) + this%numberOfPressureNodes(kk),kk) = 1.0d0
               elseif (this%circuitData_componentTypes(ll,kk) == 'i') then
                  ! insert inductor relationship into equation system
                  this%systemMatrix(ll,int(this%circuitData(ll,1,kk)),kk) = 1.0d0
                  this%systemMatrix(ll,int(this%circuitData(ll,2,kk)),kk) = -1.0d0
                  this%systemMatrix(ll,ll+this%numberOfPressureNodes(kk)+this%numberOfHistoryPressures(kk),kk) = -this%circuitData(ll,3,kk)/alfi_delt
                  this%systemMatrix(ll,this%componentIndexToFlowHistoryComponentOrderingMap(ll,kk) + this%numberOfPressureNodes(kk) + this%numberOfHistoryPressures(kk) + this%numberOfComponents(kk),kk) = this%circuitData(ll,3,kk)/alfi_delt
               else
                  write(*,*) 'EE: Unknown component type in netlist. Halting.'
                  stop
               end if
             end do

             ! Do the equations for the nodes with multiple incident currents:
             do mm=1, this%numberOfMultipleIncidentCurrentNodes(kk)
               do ll=1, this%numberOfComponents(kk)
                  if (int(this%circuitData(ll,2,kk)) .eq. this%listOfNodesWithMultipleIncidentCurrents(mm,kk)) then
                     this%systemMatrix(mm+this%numberOfComponents(kk), this%numberOfPressureNodes(kk) + this%numberOfHistoryPressures(kk) + ll,kk) = 1.0d0
                  end if
                  if (int(this%circuitData(ll,1,kk)) .eq. this%listOfNodesWithMultipleIncidentCurrents(mm,kk)) then
                     this%systemMatrix(mm+this%numberOfComponents(kk), this%numberOfPressureNodes(kk) + this%numberOfHistoryPressures(kk) + ll,kk) = -1.0d0
                  end if
               end do
             end do

             rowsDoneSoFar = this%numberOfComponents(kk) + this%numberOfMultipleIncidentCurrentNodes(kk)

             ! create the columnMap which tells us which system column each of the prescribed pressure, pressure-history or flow values belong to
             tempIndexingShift = 0 ! we use this with zero initially, just for code readability in what follows
             tempUnknownVariableIndexWithinLinearSystem = 0 ! we use this with zero initially, just for code readability in what follows
             do ll=1, this%numberOfPrescribedPressures(kk)
               this%columnMap(ll + tempIndexingShift,kk) = this%listOfPrescribedPressures(ll,kk) + tempUnknownVariableIndexWithinLinearSystem
             end do

             tempIndexingShift = tempIndexingShift + this%numberOfPrescribedPressures(kk) ! tempIndexingShift is zero before this line; I'm doing it like this for clarity & consistency
             tempUnknownVariableIndexWithinLinearSystem = tempUnknownVariableIndexWithinLinearSystem + this%numberOfPressureNodes(kk) ! tempUnknownVariableIndexWithinLinearSystem is zero before this line; I'm doing it like this for clarity & consistency
             do ll=1, this%numberOfHistoryPressures(kk)
               ! this%columnMap(ll + tempIndexingShift) = this%listOfHistoryPressures(ll) + tempUnknownVariableIndexWithinLinearSystem
               this%columnMap(ll + tempIndexingShift,kk) = ll + tempUnknownVariableIndexWithinLinearSystem
             end do

             tempIndexingShift = tempIndexingShift + this%numberOfHistoryPressures(kk)
             tempUnknownVariableIndexWithinLinearSystem = tempUnknownVariableIndexWithinLinearSystem + this%numberOfHistoryPressures(kk)
             do ll=1, this%numberOfPrescribedFlows(kk)
               this%columnMap(ll + tempIndexingShift,kk) = this%listOfPrescribedFlows(ll,kk) + tempUnknownVariableIndexWithinLinearSystem
             end do

             tempIndexingShift = tempIndexingShift + this%numberOfPrescribedFlows(kk)
             tempUnknownVariableIndexWithinLinearSystem = tempUnknownVariableIndexWithinLinearSystem + this%numberOfComponents(kk)
             do ll=1, this%numberOfHistoryFlows(kk)
               this%columnMap(ll + tempIndexingShift,kk) = ll + tempUnknownVariableIndexWithinLinearSystem
             end do

             ! Set the prescribed-value equations (i.e. pressure_1 (LHS) = pressure_1 (RHS) - so really just a way of setting the prescribed values within the linear system)
             do ll = 1, this%systemSize(kk) - rowsDoneSoFar 
               this%systemMatrix(rowsDoneSoFar + ll, this%columnMap(ll,kk),kk) = 1d0
             end do
         end do
      end subroutine generateLinearSystemFromPrescribedCircuit


      subroutine assembleRHS_netlistLPN(this,stepn)

         use calcFlowPressure, only: FlowHist
         implicit none

         ! include "mpif.h"

         class(netlistLPN) :: this
         integer, intent(in) :: stepn

         integer :: tempIndexingShift
         integer :: ll
         integer :: kk
         integer :: nn
         real*8 :: P_IM_mid_lasttimestep
         real*8 :: P_IM_mid

         ! integer :: err, rank

         ! call MPI_COMM_RANK(MPI_COMM_WORLD,rank,err)

         this%RHS = 0d0

         this%historyPressuresInLPN = this%pressuresInLPN

         do kk=1, this%numberOfLPNSurfaces
            ! Prescribed pressures
            tempIndexingShift = this%numberOfComponents(kk) + this%numberOfMultipleIncidentCurrentNodes(kk)
            do ll=1, this%numberOfPrescribedPressures(kk)
               ! 'f' for 'fixed'
               if (this%typeOfPrescribedPressures(ll,kk) .eq. 'f') then
                  this%RHS(ll + tempIndexingShift,kk) = this%valueOfPrescribedPressures(ll,kk)
               ! 'l' for 'leftVentricular'
               else if (this%typeOfPrescribedPressures(ll,kk) .eq. 'l') then

                  if ((stepn .eq. int(0)).or.(stepn .eq. int(1))) then ! treat case with no known IM pressure yet
                     P_IM_mid_lasttimestep = 5000d0 ! \todo find a better way of doing this; maybe input this value from file...
                     P_IM_mid = 5000d0 ! ... or set it based on the aortic valve state at simulation start
                  elseif (stepn .eq. int(2)) then ! treat case where only one IM pressure history point is known
                     P_IM_mid_lasttimestep = this%valueOfPrescribedPressures(ll,kk) * hrt%plv_hist(stepn-1)
                     P_IM_mid = this%valueOfPrescribedPressures(ll,kk) * hrt%plv_hist(stepn)
                  else ! get the previous intramyocardial pressure in the case where we have enough doata for this (see comment before start of "if" block)
                     P_IM_mid_lasttimestep = this%valueOfPrescribedPressures(ll,kk) * hrt%plv_hist(stepn-1) !\todo check these actually exist on first iteration
                     ! Get IM pressure for now (this will be adjusted in a moment if we're on, according to alpha in gen alpha method)
                     P_IM_mid = this%valueOfPrescribedPressures(ll,kk) * hrt%plv_hist(stepn) !\todo check these actually exist on first iteration
                  end if

                  this%RHS(ll + tempIndexingShift,kk) = P_IM_mid
                  do nn=1, this%numberOfHistoryPressures(kk)
                     !\todo check this if guard is correct
                     if (this%listOfHistoryPressures(nn,kk) .eq. this%listOfPrescribedPressures(ll,kk)) then
                        ! this%RHS(ll + tempIndexingShift + this%numberOfPrescribedPressures(kk), kk) = P_IM_mid_lasttimestep
                        this%historyPressuresInLPN(this%listOfHistoryPressures(nn,kk),kk) = P_IM_mid_lasttimestep
                     end if
                  end do


               end if
            end do
            ! History Pressures
            tempIndexingShift = tempIndexingShift + this%numberOfPrescribedPressures(kk)
            do ll=1,this%numberOfHistoryPressures(kk)
               this%RHS(ll + tempIndexingShift,kk) = this%historyPressuresInLPN(this%listOfHistoryPressures(ll,kk),kk)
            end do

            ! Prescribed Flows
            tempIndexingShift = tempIndexingShift + this%numberOfHistoryPressures(kk)
            do ll=1,this%numberOfPrescribedFlows(kk)
               ! 't 'three-D interface'
               if (this%typeOfPrescribedFlows(ll,kk) .eq. 't') then
                  this%columnIndexOf3DInterfaceFlowInLinearSystem(kk) = ll + tempIndexingShift
                  !\todo ensure this is the correct step's flow:
                  this%RHS(ll + tempIndexingShift,kk) = this%flowpntr(kk)%p ! \todo ensure this is set up properly over the kk in setflowpntrs
                  ! this%RHS(ll + tempIndexingShift) = FlowHist(lstep_passedIn+1,this%localToGlobalSurfaceIndexMap(1)) 
               ! 'f' for 'fixed'
               else if (this%typeOfPrescribedFlows(ll,kk) .eq. 'f') then
                  this%RHS(ll + tempIndexingShift,kk) = this%valueOfPrescribedFlows(ll,kk)
               end if
            end do

            ! History Flows
            tempIndexingShift = tempIndexingShift + this%numberOfPrescribedFlows(kk)
            do ll=1,this%numberOfHistoryFlows(kk)
               this%RHS(ll + tempIndexingShift,kk) = this%flowsInLPN(this%listOfHistoryFlows(ll,kk),kk)
            end do
         end do

         ! if (rank .eq. int(0)) then
         !    write(*,*) 'RHSprint', this%RHS
         ! end if

      end subroutine assembleRHS_netlistLPN


      subroutine getListOfNodesWithMultipleIncidentCurrents(this)
      ! Note that this function also counts pressure nodes which are just
      ! /between/ two components, eg. for the (two resistor) subcircuit:
      !        N0--[==R1==]--N1--[==R2==]--N2
      ! This would count node N1 as appearing twice, and do a "Kirchoff" current
      ! balance of the form "flow through R1 = flow through R2".
      !
      ! It also catches and deals with true Kirchoff equations wherea third (fourth, fifth,...)
      ! component is connected to N1.

         implicit none

         class(netlistLPN) :: this
         integer :: node
         integer :: nextWriteLocation
         integer :: numberOfTimesNodeSeen
         integer :: listOfNodesWithMultipleIncidentCurrents_temp(this%maxPressureNodes,this%numberOfLPNSurfaces)
         integer :: ii
         integer :: kk
         integer :: maxMultipleIncidentCurrents
         integer, allocatable :: nextWriteLocation2(:)

         allocate(this%numberOfMultipleIncidentCurrentNodes(this%numberOfLPNSurfaces))

         this%numberOfMultipleIncidentCurrentNodes = int(0)

         do kk=1, this%numberOfLPNSurfaces
            do node=1, this%numberOfPressureNodes(kk)
               numberOfTimesNodeSeen = int(0)
               do ii = 1, this%numberOfComponents(kk)
                  if ((int(this%circuitData(ii,1,kk)) .eq. node) .or. (int(this%circuitData(ii,2,kk)) .eq. node)) then
                     numberOfTimesNodeSeen = numberOfTimesNodeSeen + int(1)
                  end if
               end do
               if (numberOfTimesNodeSeen .gt. int(1)) then
                  ! this acts as a flag to the next loop, which will make the final listOfNodesWithMultipleIncidentCurrents.
                  listOfNodesWithMultipleIncidentCurrents_temp(node,kk) = int(1)
                  this%numberOfMultipleIncidentCurrentNodes(kk) = this%numberOfMultipleIncidentCurrentNodes(kk) + int(1)
               end if
            end do
         end do

         maxMultipleIncidentCurrents = maxval(this%numberOfMultipleIncidentCurrentNodes)
         allocate(this%listOfNodesWithMultipleIncidentCurrents(int(maxMultipleIncidentCurrents),this%numberOfLPNSurfaces))
         
         ! I'm using the allocatable array nextWriteLocation2 instead of the single nextWriteLocation counter (and the =int(0) resetting inside the loop below)
         ! because we get weird crashes otherwise. This seems to be related to a bug with ifort's loop collapse optimisation -http://software.intel.com/en-us/forums/topic/505605
         ! I can avoid it by just having the nextWriteLocation2 with allocatable attribute above, or with the write(*,*) below, or with no compiler optimisations, or with -O1, or
         ! with -O2 -no-simd -no-vec. I think this hack fixes it; try again when the bug is fixed in ifort to confirm (probably ifort > 14; follow the link: http://software.intel.com/en-us/forums/topic/505605)
         allocate(nextWriteLocation2(this%numberOfLPNSurfaces))
         nextWriteLocation2 = int(0)

         do kk=1, this%numberOfLPNSurfaces
            ! Write the final listOfNodesWithMultipleIncidentCurrents by collapsing the flagged locations from the previous loop into their corresponding node indices
            ! nextWriteLocation = int(0)
            do node = 1, this%numberOfPressureNodes(kk)
               if (listOfNodesWithMultipleIncidentCurrents_temp(node,kk) .eq. int(1)) then
                  nextWriteLocation2(kk) = nextWriteLocation2(kk) + int(1)
                  this%listOfNodesWithMultipleIncidentCurrents(nextWriteLocation2(kk),kk) = node
               end if
            end do
            ! write(*,*) 'nrlis',nextWriteLocation2(kk)
         end do
      end subroutine getListOfNodesWithMultipleIncidentCurrents


      subroutine getMapOfPressHistoriesToCorrectPressNodes(this)

         implicit none

         class(netlistLPN) :: this
         integer :: ii, jj
         integer :: cursor
         integer :: locationOfSmallestValueSoFar
         integer, allocatable :: listOfHistoryPressures_temp(:,:)
         integer :: tempValue
         integer :: maxHistoryPressures
         integer :: kk
         integer :: writeCounter

         allocate(listOfHistoryPressures_temp(2*this%maxComponents,this%numberOfLPNSurfaces))
         allocate(this%numberOfHistoryPressures(this%numberOfLPNSurfaces))
         listOfHistoryPressures_temp = int(0)
         do kk = 1, this%numberOfLPNSurfaces
            cursor = 1
            do ii=1, this%numberOfComponents(kk)
               ! Check for capacitor. Will want to do something similar for inductors and flow histories, when that functionality is added
               if (this%circuitData_componentTypes(ii,kk) == 'c') then
                  listOfHistoryPressures_temp(cursor:cursor+1,kk) = int(this%circuitData(ii,1:2,kk))
                  cursor = cursor + 2
               end if
            end do

            ! Sort the history nodes into ascending order (will have lots of zeros at beginning of array after sort)
            ! This is an implementation of the simple, inefficient "Selection Sort" algorithm.
            do ii=1, cursor - 1
               locationOfSmallestValueSoFar = ii
               do jj=ii+1, cursor - 1
                  if(listOfHistoryPressures_temp(locationOfSmallestValueSoFar,kk) > listOfHistoryPressures_temp(jj,kk)) then
                     locationOfSmallestValueSoFar = jj
                  end if
               end do
               tempValue = listOfHistoryPressures_temp(ii,kk)
               listOfHistoryPressures_temp(ii,kk) = listOfHistoryPressures_temp(locationOfSmallestValueSoFar,kk)
               listOfHistoryPressures_temp(locationOfSmallestValueSoFar,kk) = tempValue
            end do
            
            ! Remove duplicates from the history nodes array:
            ! Begin by counting the number of unique entries:

            ! First, ensure that there are any history pressures at all:
            if (maxval(listOfHistoryPressures_temp(:,kk)) .ne. int(0)) then
               tempValue = listOfHistoryPressures_temp(1,kk)
               this%numberOfHistoryPressures(kk) = int(1)
               do ii = 2, 2*this%numberOfComponents(kk)
                  if ((listOfHistoryPressures_temp(ii,kk) .ne. tempValue) .and. (listOfHistoryPressures_temp(ii,kk) .ne. int(0))) then
                     tempValue = listOfHistoryPressures_temp(ii,kk)
                     this%numberOfHistoryPressures(kk) = this%numberOfHistoryPressures(kk) + 1
                  else
                     listOfHistoryPressures_temp(ii,kk) = int(0) ! we just do this for ease when removing the duplicates in the next loop... read ahead and you'll see what I mean
                  end if
               end do
            else
               this%numberOfHistoryPressures(kk) = int(0)
            end if
         end do

         ! Now get the unique listOfHistoryPressures array, using the number of unique entries, and the fact that listOfHistoryPressures_temp now contains unique values, interspersed with zeros
         maxHistoryPressures = maxval(this%numberOfHistoryPressures)

         allocate(this%nodeIndexToPressureHistoryNodeOrderingMap(this%maxPressureNodes,this%numberOfLPNSurfaces))
         this%nodeIndexToPressureHistoryNodeOrderingMap = int(-1) ! this initialisation to a clearly-invalid index should help catch indexing bugs

         ! Copy the sorted values into the member array:
         allocate(this%listOfHistoryPressures(maxHistoryPressures,this%numberOfLPNSurfaces))

         do kk=1, this%numberOfLPNSurfaces
            writeCounter = 1
            do ii = 1, 2*this%numberOfComponents(kk)
               if (listOfHistoryPressures_temp(ii,kk) .ne. int(0)) then
                  this%listOfHistoryPressures(writeCounter,kk) = listOfHistoryPressures_temp(ii,kk)
                  writeCounter = writeCounter + 1
               end if
            end do

            ! Now do the actual generation of the pressure history node ordering map:
            do ii=1, this%numberOfHistoryPressures(kk)
               this%nodeIndexToPressureHistoryNodeOrderingMap(this%listOfHistoryPressures(ii,kk),kk) = ii
            end do
         end do

      end subroutine getMapOfPressHistoriesToCorrectPressNodes


      subroutine getMapOfFlowHistoriesToCorrectComponents(this)

         implicit none

         class(netlistLPN) :: this
         integer :: ii, jj
         integer :: cursor
         integer :: locationOfSmallestValueSoFar
         integer, allocatable :: listOfHistoryFlows_temp(:,:)
         integer :: tempValue
         integer :: maxHistoryFlows
         integer :: kk

         allocate(listOfHistoryFlows_temp(2*this%maxComponents,this%numberOfLPNSurfaces))
         allocate(this%numberOfHistoryFlows(this%numberOfLPNSurfaces))
         listOfHistoryFlows_temp = int(0)
         do kk = 1, this%numberOfLPNSurfaces
            cursor = int(0)
            do ii=1, this%numberOfComponents(kk)
               ! Check for inductor.
               if (this%circuitData_componentTypes(ii,kk) == 'i') then
                  cursor = cursor + int(1)
                  listOfHistoryFlows_temp(cursor,kk) = ii
               end if
            end do

            ! Sort the history flows into ascending order (will have lots of zeros at beginning of array after sort)
            ! This is an implementation of the simple, inefficient "Selection Sort" algorithm.
            do ii=1, cursor - int(1)
               locationOfSmallestValueSoFar = ii
               do jj=ii+1, cursor - 1
                  if(listOfHistoryFlows_temp(locationOfSmallestValueSoFar,kk) > listOfHistoryFlows_temp(jj,kk)) then
                     locationOfSmallestValueSoFar = jj
                  end if
               end do
               tempValue = listOfHistoryFlows_temp(ii,kk)
               listOfHistoryFlows_temp(ii,kk) = listOfHistoryFlows_temp(locationOfSmallestValueSoFar,kk)
               listOfHistoryFlows_temp(locationOfSmallestValueSoFar,kk) = tempValue
            end do

            
            ! Remove duplicates from the history flows array:
            ! Begin by counting the number of unique entries:

            ! First, ensure that there are any history flows at all:
            if (maxval(listOfHistoryFlows_temp(:,kk)) .ne. int(0)) then
               tempValue = listOfHistoryFlows_temp(1,kk)
               this%numberOfHistoryFlows(kk) = int(1)
               do ii = 1, 2*this%numberOfComponents(kk)
                  if ((listOfHistoryFlows_temp(ii,kk) .ne. tempValue) .and. (listOfHistoryFlows_temp(ii,kk) .ne. int(0))) then
                     tempValue = listOfHistoryFlows_temp(ii,kk)
                     this%numberOfHistoryFlows(kk) = this%numberOfHistoryFlows(kk) + 1
                  ! else
                  !    listOfHistoryFlows_temp(ii,kk) = int(0) ! we just do this for ease when removing the duplicates in the next loop... read ahead and you'll see what I mean
                  end if
               end do
            else
               this%numberOfHistoryFlows(kk) = int(0)
            end if
            ! Now get the unique listOfHistoryFlows array, using the number of unique entries, and the fact that listOfHistoryFlows_temp now contains unique values, interspersed with zeros
         end do

         maxHistoryFlows = maxval(this%numberOfHistoryFlows)

         allocate(this%componentIndexToFlowHistoryComponentOrderingMap(this%maxComponents,this%numberOfLPNSurfaces))
         this%componentIndexToFlowHistoryComponentOrderingMap = int(-1) ! this initialisation to a clearly-invalid index should help catch indexing bugs

         ! Copy the sorted values into the member array:
         allocate(this%listOfHistoryFlows(maxHistoryFlows,this%numberOfLPNSurfaces))
         this%listOfHistoryFlows = int(0)
         do kk=1, this%numberOfLPNSurfaces
            do ii = 1, 2*this%numberOfComponents(kk)
               if (listOfHistoryFlows_temp(ii,kk) .ne. int(0)) then
                  this%listOfHistoryFlows(ii,kk) = listOfHistoryFlows_temp(ii,kk)
               end if
            end do

            ! Now do the actual generation of the pressure history node ordering map:
            do ii=1, this%numberOfHistoryFlows(kk)
               this%componentIndexToFlowHistoryComponentOrderingMap(this%listOfHistoryFlows(ii,kk),kk) = ii
            end do
         end do

      end subroutine getMapOfFlowHistoriesToCorrectComponents
      

!     *************** END NETLIST ARBITRARY-DESIGN LPN MODEL CODE ***************


!
! *******************************************************************
! *** set key simvascular parameters to internal module variables ***
! *******************************************************************
!
      subroutine setsimv_maxsurf(max_surf_val)
      implicit none 
      integer :: max_surf_val
      maxsurf = max_surf_val
      end subroutine setsimv_maxsurf
!
      subroutine setsimv_delt(dtime_val)
      implicit none 
      real*8 :: dtime_val
      delt = dtime_val
      end subroutine setsimv_delt
!
      subroutine setsimv_alfi(alfi_val)
      implicit none 
      real*8 :: alfi_val
      alfi = alfi_val
      end subroutine setsimv_alfi
!
      subroutine setsimv_ntout(ntout_val)
      implicit none 
      integer :: ntout_val
      ntout = ntout_val
      end subroutine setsimv_ntout
!
      subroutine setsimv_nstep(nstep_val)
      implicit none 
      integer :: nstep_val
      nstep = nstep_val
      end subroutine setsimv_nstep    
!
      subroutine setsimv_lstep(lstep_val)
      implicit none 
      integer :: lstep_val
      lstep = lstep_val
      end subroutine setsimv_lstep    
!

      subroutine setsimv_rho(rho_val)
      implicit none 
      real*8 :: rho_val
      rho = rho_val
      end subroutine setsimv_rho
!      
! *******************************************************************
! *** start multidomain subroutine - activates logicals in module ***
! *******************************************************************
!
! *** start multidomain with multidomain.dat file
!
      subroutine startmultidomain_file() 
      implicit none
      integer :: ierr, rerr      
      integer :: fnum = 671
      real*8 :: testvalue
!
      open(fnum, file='multidomain.dat', status='old', iostat=ierr)   

      if (ierr .ne. int(0)) then            
         if (multidomainactive .ne. int(1)) then
            multidomainactive = int(0)
            nrcractive = int(0)
            ntrcractive = int(0)
         end if 
!            
      else  
!         
         if (multidomainactive .ne. int(1)) then   
            multidomainactive = int(1)
         end if 
!           
         ! numerical rcr on/off   
         read(fnum,*)
         read(fnum,*,iostat=rerr) testvalue
!            
         if (testvalue .gt. int(0)) then
            nrcractive = int(1)          
         else
            nrcractive = int(0)
         end if
!
         ! numerical time varying rcr on/off            
         read(fnum,*)
         read(fnum,*,iostat=rerr)  testvalue
!            
         if (testvalue .gt. int(0)) then
            ntrcractive = int(1)
         else
            ntrcractive = int(0)        
         end if

!                       
      end if
!         
      close(fnum)     
!
      end subroutine
!      
! *** start multidomain with integer 
!
      subroutine startmultidomain_int(iflag) 
      implicit none
      integer :: iflag
!
      if (iflag .gt. int(0)) then
         if (multidomainactive .ne. int(1)) then   
            multidomainactive = int(1)
         end if
      end if
!
      end subroutine
!
! *** start multidomain with integer and character string 
!      
      subroutine startmultidomain_int_char(iflag,charflag)          
      implicit none
      integer, intent(in) :: iflag
      character(len=*), intent(in) :: charflag
      character(len=*), parameter :: heartflag = 'heart'
      character(len=*), parameter :: systemicflag = 'systemic'
      character(len=*), parameter :: coronaryflag = 'coronary'
      character(len=*), parameter :: netlistflag = 'netlist'

      integer :: ierr, rerr
      integer :: fnum = 505
      integer :: gnum = 606
      integer :: testvalue  
!
      if( iflag .gt. int(0)) then
! 
!        ! systemic circuit      
         if (charflag .eq. systemicflag) then
!        
            open(gnum, file='sys.dat', status='old', iostat=ierr)
!            
               if (ierr .ne. int(0)) then            
!               
                  write(*,*) '** '
                  write(*,*) '** Error!! Cannot open file: "sys.dat"'
                  write(*,*) '** '               
                  stop
!                  
               else
!               
                  if (multidomainactive .ne. int(1)) then   
                     multidomainactive = int(1)
                  end if
!                  
                  sysactive = int(1)
                  nrcractive = int(1)
!                  
               end if      
            close(gnum)
!               
!        ! heart circuit
         elseif (iflag .gt. int(0) .and. charflag .eq. heartflag) then        
!
            if (multidomainactive .ne. int(1)) then   
               multidomainactive = int(1)
            end if
!            
            hrtactive = int(1)

         elseif (charflag .eq. coronaryflag) then
            newCoronaryActive = int(1)
         elseif (charflag .eq. netlistflag) then
            netlistActive = int(1)
         end if
!         
      end if
!
      end subroutine 
!
!
! 
      subroutine multidomainstatus()
      
      use mpi, only : MPI_COMM_RANK, MPI_COMM_WORLD
      implicit none
      character(len=20) :: format = '(5x,a20,a5,5x)'
      integer :: mpirank, mpierr
     
      call MPI_COMM_RANK(MPI_COMM_WORLD, mpirank, mpierr)        

      if (mpirank .eq. int(0)) then

         write(*,*)  
         write(*,'(a35)') ' ************************************'
         write(*,*)  
         if (multidomainactive) then
            write(*,format) 'Multidomain module: ','On'
            if (nrcractive) then
               write(*,format) 'Numerical RCR: ','On'
            else 
               write(*,format) 'Numerical RCR: ','Off'
            end if
            if (ntrcractive) then
               write(*,format) 'Numerical TRCR: ','On'
            else 
               write(*,format) 'Numerical TRCR: ','Off'
            end if
            if (hrtactive) then
               write(*,format) 'Numerical heart: ','On'
            else 
               write(*,format) 'Numerical heart: ','Off'
            end if
            if (sysactive) then
               write(*,format) 'Systemic circuit: ','On'
            else 
               write(*,format) 'Systemic circuit: ','Off'
            end if
            if (newCoronaryActive) then
               write(*,format) 'Coronary control: ','On'
            else
               write(*,format) 'Coronary control: ','Off'
            end if
            if (netlistActive) then
               write(*,format) 'Netlist surfaces: ', 'On'
            else
               write(*,format) 'Netlist surfaces: ', 'Off'
            end if
         else
            write(*,format) 'Multidomain module: ','Off'
         end if
         write(*,*)  
         write(*,'(a35)') ' ************************************'
         write(*,*)  

      end if 
      
      end subroutine
!

      subroutine initmultidomain(isys)	  
!
      implicit none
      integer :: isys ! systemic flag
      integer :: ierr, rerr
      integer :: fnum = 505
      integer :: gnum = 606
      integer :: testvalue
!
      if (isys .gt. int(0)) then
         open(gnum, file='sys.dat', status='old', iostat=ierr)
         if (ierr .ne. int(0)) then
            sysactive = int(0)
            write(*,*) '** '
            write(*,*) '** Error!! Cannot open file: "sys.dat"'
            write(*,*) '** '
            close(gnum)
            stop
         else
            ! set sysytemic circuit active
            sysactive = int(1)
            write(*,*) '** '            
            write(*,*) '** Systemic circuit: On'
            write(*,*) '** '          
            ! set multidomain active
            multidomainactive = int(1)
            ! set numerical rcr active
            nrcractive = int(1)
         end if      
      else 
         sysactive = int(0)
         write(*,*) '** '         
         write(*,*) '** Systemic circuit: Off'
         write(*,*) '** '         
      end if 
!
      return      
      end subroutine
!
! ************************************************
! *** reduced order type get and set functions ***
! ************************************************
!
! *** initialise internal variables 
!
      subroutine initxvars(ro,stepnum)
      implicit none
!
      class(reducedorder) :: ro
      integer :: stepnum

      ! if (ro%classNameString .eq. 'numericalheart') then
      !    call ro%initxvars_hrt(stepnum)
      ! else if (ro%classNameString .eq. 'systemiccircuit') then
      !    call ro%initxvars_sys(stepnum)
      ! else if (ro%classNameString .eq. 'numericalrcr') then
      ! else if (ro%classNameString .eq. 'numericaltrcr') then
      ! end if

      select type (ro)
         type is (numericalrcr) 
         type is (numericaltrcr) 
         type is (numericalheart)
            call ro%initxvars_hrt(stepnum)
         type is (systemiccircuit)
            call ro%initxvars_sys(stepnum)
      end select 
!      
      return
      end subroutine initxvars
!
! *** update internal variables 
!
      subroutine updxvars(ro,stepnum)
      implicit none
!
      class(reducedorder) :: ro
      integer :: stepnum

      ! if (ro%classNameString .eq. 'numericalheart') then
      !    call ro%updxvars_hrt(stepnum)
      ! else if (ro%classNameString .eq. 'numericalrcr') then

      ! else if (ro%classNameString .eq. 'numericaltrcr') then

      ! end if

      !
      select type (ro)
         type is (numericalrcr) 
            call ro%updxvars_rcr(stepnum)
         type is (numericaltrcr) 
         type is (numericalheart)
            call ro%updxvars_hrt(stepnum)
      end select 
!      
      return
      end subroutine
!
! *** write internal variables
!
      subroutine writexvars(ro,stepnum)
      implicit none
!
      class(reducedorder) :: ro
      integer :: stepnum
!

      ! if (ro%classNameString .eq. 'numericalheart') then
      !    call ro%writexvars_hrt(stepnum)
      ! else if (ro%classNameString .eq. 'numericalrcr') then

      ! else if (ro%classNameString .eq. 'numericaltrcr') then

      ! end if

      select type (ro)
         type is (numericalrcr) 
            call ro%writexvars_rcr(stepnum)
         type is (numericaltrcr) 
         type is (numericalheart)
            call ro%writexvars_hrt(stepnum)
      end select 
!      
      return
      end subroutine
!
! *** get implicit coefficients
!
      function getimplicitcoeff(a)
      implicit none
      class(reducedorder) :: a
      real*8 :: getimplicitcoeff(a%surfnum,2)
      getimplicitcoeff(:,:) = a%implicitcoeff(:,:)      
!      
!!      real*8 :: getimplicitcoeff(0:maxsurf,2)
!!      getimplicitcoeff(1:a%surfnum,1:2) = a%implicitcoeff(:,:)      
!      
      end function getimplicitcoeff
!
! *** set implicit coefficient, calls individual model type 
!     subroutines to set values 
!
      subroutine setimplicitcoeff(ro,stepnum)
      implicit none 
      class(reducedorder) :: ro
      integer :: stepnum
      real*8 :: coeffs(ro%surfnum,2)

      ! if (ro%classNameString .eq. 'numericalrcr') then
      !    call ro%setimplicitcoeff_rcr(stepnum,'solve')
      !    call ro%setimplicitcoeff_rcr(stepnum,'update')
      ! else if (ro%classNameString .eq. 'numericaltrcr') then
      !    coeffs = ro%setimplicitcoeff_trcr(stepnum)
      !    ro%implicitcoeff(:,:) = coeffs(:,:)
      ! else if (ro%classNameString .eq. 'numericalheart') then
      !    call ro%setimplicitcoeff_hrt(stepnum,'solve') 
      !    ! get implicit coefficients for t = t_{n+1}
      !    call ro%setimplicitcoeff_hrt(stepnum,'update')
      ! else if (ro%classNameString .eq. 'controlledcoronarymodel') then
      !    call ro%setimplicitcoeff_controlledCoronary(stepnum,'solve')!\todo check this solves/updates in all the right places
      !    call ro%setimplicitcoeff_controlledCoronary(stepnum,'update')!\todo check that these should be done one after another like this...
      ! end if
      select type (ro)
         class is (numericalrcr)
            call ro%setimplicitcoeff_rcr(stepnum,'solve')
            call ro%setimplicitcoeff_rcr(stepnum,'update')
         class is (numericaltrcr)
            call ro%setimplicitcoeff_trcr(stepnum,'solve')
            call ro%setimplicitcoeff_trcr(stepnum,'update')
            ! coeffs = ro%setimplicitcoeff_trcr(stepnum) !\todo remove
            ! ro%implicitcoeff(:,:) = coeffs(:,:) !\todo remove
         class is (numericalheart)
            call ro%setimplicitcoeff_hrt(stepnum,'solve') 
            ! get implicit coefficients for t = t_{n+1}
            call ro%setimplicitcoeff_hrt(stepnum,'update')
         class is(controlledCoronaryModel)
            call ro%setimplicitcoeff_controlledCoronary(stepnum,'solve')!\todo check this solves/updates in all the right places
            call ro%setimplicitcoeff_controlledCoronary(stepnum,'update')!\todo check that these should be done one after another like this...
         class is(netlistLPN)
            call ro%setimplicitcoeff_netlistLPN(stepnum,'solve')
            call ro%setimplicitcoeff_netlistLPN(stepnum,'update')
      end select

      end subroutine setimplicitcoeff
!
! *** get number of surfaces in reduced order model
!
      function getsurfnum(a)
      implicit none
      class(reducedorder) :: a
      integer :: getsurfnum
      getsurfnum = a%surfnum
      end function getsurfnum 
!
! *** get surface ids in reduced order model
!
      function getsurfids(a)
      implicit none
      class(reducedorder) :: a
      integer :: getsurfids(0:maxsurf)
      getsurfids(:) = 0
      getsurfids(1:a%surfnum) = a%surfids(1:a%surfnum)
      end function getsurfids
!
! *** has surface id in reduced order model 
!
      function hassurfid(a,surfid) 
      implicit none
      class(reducedorder) :: a
      integer :: surfid
      logical :: hassurfid
      integer :: i
!
!     ! loop through surfaces and check
      hassurfid = .false.
      do i = 1, a%surfnum       
         if (a%surfids(i) .eq. surfid) then
            hassurfid = .true.
            exit
         end if
      end do
!
      end function
!
! *** get surface index
!
      function getsurfidindx(a,surfid) result(surfindx)
      implicit none
      class(reducedorder) :: a
      integer :: surfid
      integer :: surfindx
      integer :: i
      do i = 1, a%surfnum
         if (a%surfids(i) .eq. surfid) then
            surfindx = i
            exit
         end if
      end do 
      end function 
!
! *** load flow file in reduced order model
!
      subroutine loadflowfile(a,torow,charvar)
      implicit none
      class(reducedorder) :: a
      integer :: ierr
      integer :: torow
      real*8 :: flowfile(torow,a%surfnum)
      integer :: flownum = 689
      real*8 :: temp     
      integer :: i
      integer :: j
      character(len=*) :: charvar
      character(len=50) :: legacyformat = 'legacy'
      character(len=50) :: multidomformat = 'multidomain'      

      if (charvar .eq. legacyformat) then
         
         ! open file, read header and values
         open(flownum, file=a%flowfile, status='old', iostat=ierr)
         read(flownum,*) temp
         

         do i = 1, torow
           read(flownum,*) flowfile(i,:)
         end do
         
         close(flownum)         

      else if (charvar .eq. multidomformat) then         

         ! open file, read header and values
         open(flownum, file=a%flowfile, status='old', iostat=ierr)
         
         ! data now organised with step number in first column         
         do i = 1, torow
            ! read(flownum,*) j, flowfile(i,1:a%surfnum)
            read(flownum,*) j, a%flowhist(i,1:a%surfnum)
         end do
         
         close(flownum)

      end if 

      ! set flows
      ! a%flow_n(1:a%surfnum) = flowfile(torow,1:a%surfnum)
      a%flow_n(1:a%surfnum) = a%flowhist(torow,1:a%surfnum)
      
      end subroutine
!
! *** load pressure file
!
      subroutine loadpressurefile(a,torow,charvar)
      implicit none
      class(reducedorder) :: a
      integer :: ierr
      integer :: torow
      real*8 :: pressfile(torow,a%surfnum)
      integer :: presnum = 374
      real*8 :: temp     
      integer :: i
      integer :: j
      character(len=*) :: charvar
      character(len=50) :: legacyformat = 'legacy'
      character(len=50) :: multidomformat = 'multidomain'      
!       

      if (charvar .eq. legacyformat) then
         
         ! open file, read header and values      
         open(presnum, file=a%pressurefile, status='old', iostat=ierr)
         read(presnum,*) temp
         
         do i = 1, torow
            read(presnum,*) pressfile(i,:)
         end do
         
         close(presnum)

      elseif (charvar .eq. multidomformat) then    
         
         ! open file, read header and values      
         open(presnum, file=a%pressurefile, status='old', iostat=ierr)
         
         ! data now organised with step number in first column                  
         do i = 1, torow
            !read(presnum,*) j, pressfile(i,1:a%surfnum)
            read(presnum,*) j, a%pressurehist(i,1:a%surfnum)      
         end do
         
         close(presnum)
      
      end if

      ! set pressure
      ! a%pressure_n(1:a%surfnum) = pressfile(torow,1:a%surfnum)
      a%pressure_n(1:a%surfnum) = a%pressurehist(torow,1:a%surfnum)      
      
!      write(*,*) a%pressure_n(1:a%surfnum)
 
      end subroutine



!     
! *** set area subroutine 
!
      subroutine setarea(a,narea,surfacearea)
      implicit none
      class(reducedorder) :: a
      integer :: narea
      real*8 :: surfacearea(0:maxsurf)
      a%surfarea(1:narea) = surfacearea(1:narea)          
      return
      end subroutine 
!
! *** get area subroutine 
!
      function getarea(a)
      implicit none
      class(reducedorder) :: a
      real*8 :: getarea(0:maxsurf)
      getarea(1:a%surfnum) = a%surfarea(1:a%surfnum)      
      return
      end function
!
! *** set flow_n subroutine 
!
      subroutine setflow_n(a,nflows,flows)
      implicit none
      class(reducedorder) :: a
      integer :: nflows
      real*8 :: flows(0:maxsurf)
      a%flow_n(1:nflows) = flows(1:nflows)     
      return
      end subroutine 
!
! *** set flow_n1 subroutine 
!
      subroutine setflow_n1(a,nflows,flows)
      implicit none
      class(reducedorder) :: a
      integer :: nflows
      real*8 :: flows(0:maxsurf)
      a%flow_n1(1:nflows) = flows(1:nflows)
      return
      end subroutine       
!
! *** set pressure_n subroutine, only called at initialisation 
!     with no pressure history
!
      subroutine setpressure_n(a,npressures,pressures)
      implicit none
      class(reducedorder) :: a
      integer :: npressures
      real*8 :: pressures(0:maxsurf)
      a%pressure_n(1:npressures) = pressures(1:npressures)           
      return
      end subroutine 
!
! *** is pressure update logical, if true pressure is updated from 
!     integrated value
!
      function ispressureupdate(ro)
      implicit none      
      class(reducedorder) :: ro
      integer :: ispressureupdate      
      ispressureupdate = ro%updatepressure      
      return
      end function      
!
! *** update pressure_n_with flow subroutine, called with final update before
!     moving onto next step, uses flow at n+1 which is now stored in flow_n
!
      subroutine updpressure_n1_withflow(a)
      implicit none
      class(reducedorder) :: a
      integer :: i
!
      do i = 1, a%surfnum
         a%pressure_n(i) = a%implicitcoeff_n1(i,1)*a%flow_n(i) &
                         + a%implicitcoeff_n1(i,2)
      end do
!
      return
      end subroutine


!     assign the external pointers
      subroutine assign_ptrs_ext(ro)

      implicit none

      class(reducedorder) :: ro
      select type (ro)
         class is (numericalrcr)
            call ro%assign_ptrs_ext_rcr()
         class is (numericalheart)
            call ro%assign_ptrs_ext_hrt()

      end select
      end subroutine


!
! *** update pressure_n_withvalue subroutine, called with final update before
!     moving onto the next step, here the value in the 3D domain is set as the 
!     pressure value
!
      subroutine updpressure_n1_withvalue(a, pressurevalue)
      implicit none
      class(reducedorder) :: a
      real*8 :: pressurevalue(0:maxsurf)
      integer :: i
!
      do i = 1, a%surfnum
         a%pressure_n(i) = pressurevalue(i)
      end do
!
      return
      end subroutine
!
! ********************************************************
! *** advanced reduced order functions and subroutines ***
! ********************************************************
!
! *** solve 
!
      subroutine solve(this,stepn)
      implicit none
      class(advancedreducedorder) :: this 
      integer :: stepn
      select type (this)
         class is (systemiccircuit) 
         if (sys%closedloop) then
            call sys%solve_sys(stepn)
         else
            call sys%solve_sys_open(stepn)
         end if
      end select
      end subroutine
!
! *** update
!
      subroutine update(this,stepn)
      implicit none
      class(advancedreducedorder) :: this 
      integer :: stepn
      select type (this)
         class is (systemiccircuit) 
         if (sys%closedloop) then
            call sys%update_sys(stepn)
         else
            call sys%update_sys_open(stepn)
         end if
      end select
      end subroutine
!
! *** is feedback active
!
      function isfeedbackactive(this)
      implicit none
      class(advancedreducedorder) :: this 
      logical :: isfeedbackactive
      isfeedbackactive = this%feedbackactive
      return
      end function isfeedbackactive

      subroutine writeRestart(this,lstep_passedIn)
         implicit none
         class(advancedreducedorder) :: this
         integer, intent(in) :: lstep_passedIn

         select type(this)
            class is (controlledCoronaryModel)
               call this%writeRestart_controlledCoronaryModel(lstep_passedIn)
            class is(netlistLPN)
               call this%writeRestart_netlistLPNSurfaces()
         end select

         ! if (this%classNameString .eq. 'controlledCoronaryModel') then
         !    call this%writeRestart_controlledCoronaryModel()
         ! else if (this%classNameString .eq. 'netlistLPN') then
         !    call this%writeRestart_netlistLPNSurfaces()
         ! end if
      end subroutine writeRestart
            
!!c
!!c *** add surfids
!!c      
!!      subroutine addsurfs(this,surfnum,surfids,varchar)
!!      implicit none
!!      class(advancedreducedorder) :: this 
!!      integer :: surfnum
!!      integer :: surfids(0:maxsurf)
!!      character(len=*) :: varchar
!!      character(len=*), parameter :: feedbackchar = 'feedback'
!!      if (varchar .eq. feedbackchar) then
!!         this%feedbacksurfs%num = surfnum
!!         allocate(this%feedbacksurfs%ids(surfnum))
!!         this%feedbacksurfs%ids(1:surfnum) = surfids(1:surfnum)
!!c        ! array of pointers
!!         allocate(this%feedbackpntr(surfnum))         
!!      end if
!!      return
!!      end subroutine
!
!
! *******************************************************
! *** multidomain container functions and subroutines ***
! *******************************************************
!
! *** constructor 
!
      function multidomconstructor() result(a)
      implicit none
      type(multidomaincontainer) :: a
      end function multidomconstructor
!
! *** add surfaces ids
!
      subroutine addsurfids(a,num,ids)
      implicit none
!
      class(multidomaincontainer) :: a      
      integer :: num
      integer :: ids(0:maxsurf)
      integer :: i
      integer :: j
      logical :: addid = .false.
!
!     ! check if number of surfaces is > 0
      if (num .gt. int(0)) then
!        ! check to see if the pointer is associated          
         if (associated(a%surfids)) then
!           ! check if id already exists
!           ! loop over new list
            do i = 1, num
!              ! loop over old list
               do j = 1, a%surfnum
                  if (ids(i) .eq. a%surfids(j)) then
                     addid = .false.
                     exit
                  else
                     addid = .true.
                  end if 
               end do 
               if (addid) then
                  ! add surface ids
                  a%surfids(a%surfnum+1) = ids(i)
                  ! update number of surfaces by one
                  !!a%surfnum = a%surfnum + num
                  a%surfnum = a%surfnum + int(1)
!                  write(*,*) 'a%surfnum:',a%surfnum
!                  write(*,*) 'a%surfid:',a%surfids(a%surfnum)
               end if
            end do
         else
!
            ! allocate surfids array, size 0:maxsurf 
            allocate(a%surfids(0:maxsurf))     
            a%surfids(0:maxsurf) = int(0)
!
            ! add surface ids
            a%surfids(1:num) = ids(1:num)
!
            ! update surface number (initialised as 0)
            a%surfnum = a%surfnum + num
!
            ! allocate area array, size 0:maxsurf
            allocate(a%areas(0:maxsurf))     
            a%areas(0:maxsurf) = real(0.0,8)
!
            ! allocate pressure_n array, size 0:maxsurf
            allocate(a%pressure(0:maxsurf))     
            a%pressure(0:maxsurf) = real(0.0,8)          
!
            ! allocate flow_nalf array, size 0:maxsurf
            allocate(a%flow(0:maxsurf))     
            a%flow(0:maxsurf) = real(0.0,8)            
!
            ! allocate flow_nalf derivative array, size 0:maxsurf
            allocate(a%flowderivative(0:maxsurf))     
            a%flowderivative(0:maxsurf) = real(0.0,8)   
!
            ! allocate stabilisation pressure array, size 0:maxsurf
            allocate(a%stabilisationpressure(0:maxsurf))     
            a%stabilisationpressure(0:maxsurf) = real(0.0,8)   

            allocate(a%stb_pres(0:maxsurf))     
            a%stb_pres(0:maxsurf) = real(0.0,8)
!            
        end if         
      end if
! #if EXTRA_CONSOLE_OUTPUT == 1
!       write(*,*) 'addsurfids: cor info', a%surfids(1:5) !\todo remove
! #endif
!
      end subroutine
!
! *** set flow pointers
!
      subroutine setflowpntrs(this, container)
      implicit none
      class(advancedreducedorder) :: this
      class(multidomaincontainer) :: container      
      select type (this)
         class is (systemiccircuit) 
            call this%setflowpntrs_sys(container)
         class is (controlledCoronaryModel)
            call this%setflowpntrs_coronary(container)
         class is (netlistLPN)
            call this%setflowpntrs_netlistLPN(container)
      end select
      end subroutine setflowpntrs
!
! *** set flow pointers
!
      subroutine setpresspntrs(this, container)
      implicit none
      class(advancedreducedorder) :: this
      class(multidomaincontainer) :: container
      select type (this)
         class is (systemiccircuit) 
         call this%setpresspntrs_sys(container)
      end select
      end subroutine setpresspntrs

      subroutine setpresspntrs_ro(this, container)
      implicit none
      class(reducedorder) :: this
      class(multidomaincontainer) :: container
      select type (this)
         class is (numericalheart) 
         call this%setpresspntrs_hrt(container)
      end select
      end subroutine setpresspntrs_ro      
!
!!      subroutine setspresspntr(this, container)
!!      implicit none
!!      class(reducedorder) :: this
!!      class(multidomaincontainer) :: container
!!      select type (this)
!!         class is (numericalheart) 
!!         call hrt%setspresspntr_hrt(container)
!!      end select
!!      end subroutine setspresspntr
            
!
! *** set feeback pointers
!
      subroutine setfeedbackpntrs(this, container)
      implicit none
      class(advancedreducedorder) :: this
      class(multidomaincontainer) :: container
      select type (this)
         class is (systemiccircuit) 
         call this%setfeedbackpntrs_sys(container)
      end select
      end subroutine setfeedbackpntrs
!
! *** get number of surfaces
!
      function getsurfnum_mdc(this) 
      implicit none
      class(multidomaincontainer) :: this
      integer :: getsurfnum_mdc
      getsurfnum_mdc = this%surfnum
      end function
!
! *** get surface ids
!
      function getsurfids_mdc(this)
      implicit none
      class(multidomaincontainer) :: this
      integer :: getsurfids_mdc(0:maxsurf)
      getsurfids_mdc(0:maxsurf) = this%surfids(0:maxsurf)
      end function
!
! *** set area 
!
      subroutine setarea_mdc(this,nareas,areas)
      implicit none
      class(multidomaincontainer) :: this
      integer :: nareas
      real*8 :: areas(0:maxsurf)
      this%areas(1:nareas) = areas(1:nareas)
      end subroutine      
!
! *** get area
!
      function getarea_mdc(this)
      implicit none
      class(multidomaincontainer) :: this
      real*8 :: getarea_mdc(0:maxsurf)
      getarea_mdc = this%areas
      end function
!
! *** set pressure at t_{n}
!
      subroutine setpressure_mdc(this, npress, press)
      implicit none
      class(multidomaincontainer) :: this     
      integer :: npress
      real*8 :: press(0:maxsurf)
      this%pressure(1:npress) = press(1:npress)          
      end subroutine
!
! *** set flow at t_{n+alf_{i}}
!
      subroutine setflow_mdc(this, nflows, flows)
      implicit none
      class(multidomaincontainer) :: this     
      integer :: nflows
      real*8 :: flows(0:maxsurf)
      this%flow(1:nflows) = flows(1:nflows)  
      ! write(*,*) this%flow(1:nflows)       
      end subroutine        
!
! *** set flow derivative at t_{n+alf_{i}}
!
      subroutine setflowderivative_mdc(this, nflows, dflows)
      implicit none
      class(multidomaincontainer) :: this     
      integer :: nflows
      real*8 :: dflows(0:maxsurf)
      this%flowderivative(1:nflows) = dflows(1:nflows)          
      end subroutine   
!
! *** reset stabilisation pressure
!      
      subroutine resetstb_pres(this)      
      implicit none
      class(multidomaincontainer) :: this
      this%stb_pres(:) = real(0.0,8)
      end subroutine
!
! *** add stabilisation pressure
!
      subroutine addstb_pres(this,nelm,pres,ibcb)
      implicit none
      include "mpif.h" 
      class(multidomaincontainer) :: this

      integer :: nelm                     ! number of element 
      real*8  :: pres(nelm)               ! stabilised pressure 
      integer :: ibcb(nelm)               ! surface ids
      integer :: i, j, ierr               !
      real*8  :: psum(this%surfnum)       ! 
      real*8  :: psum_mpi(this%surfnum)   !
      real*8  :: sendarray(this%surfnum)  !
      real*8  :: recvarray(this%surfnum)  !
      integer :: num 
!
      !       
      num = this%surfnum
!
!     ! zero pressure on all surfaces
      psum(:) = real(0.0,8)
!
!     ! loop through surfaces in container
      do i = 1, this%surfnum         
!        ! loop through elements      
         do j = 1, nelm
!           ! if this surface is in the container
!           ! add pressure
            if (this%surfids(i) .eq. ibcb(j))then
               psum(i) = psum(i) + pres(j)
            end if             
         end do
      end do
!!      
!!      sendarray(:) = psum
!!      recvarray(:) = real(0.0,8)
!!c
!!c     ! allreduce to communicate sum over all processors
!!      call MPI_ALLREDUCE(sendarray, 
!!     &                   recvarray, 
!!     &                   num,
!!     &                   MPI_DOUBLE_PRECISION, 
!!     &                   MPI_SUM, 
!!     &                   MPI_COMM_WORLD, 
!!     &                   ierr)  
!!
!!      write(*,*) 'ierr: ', recvarray
!!      psum_mpi(:) = recvarray(:)
!
!     ! add contribution of these nodes to the surface
!     ! this sum is reset in elmgmr.f before looping 
!     ! through the boundary element blocks
      !!do i = 1, this%surfnum
      !!   if
      !!   end if
      !!end do 

      this%stb_pres(1:this%surfnum) = this%stb_pres(1:this%surfnum) &
                                    + psum(1:this%surfnum)

      ! write(*,*) 'addstb_pres: ',this%stb_pres(1:this%surfnum)

      end subroutine      
!
! *** 
!
      subroutine sumstb_pres(a)
      implicit none
      include "mpif.h" 
      class(multidomaincontainer) :: a
      real*8  :: psum_mpi(a%surfnum)   !
      real*8  :: parea(a%surfnum)
      !real*8  :: sendarray(this%surfnum)  !
      !real*8  :: recvarray(this%surfnum)  !
      integer :: ierr !, num 
!   
!     ! allreduce to communicate sum over all processors
      call MPI_ALLREDUCE(a%stb_pres(1:a%surfnum), &
                         psum_mpi, &
                         a%surfnum, &
                         MPI_DOUBLE_PRECISION, &
                         MPI_SUM, &
                         MPI_COMM_WORLD, &
                         ierr)  

      !!write(*,*) 'psum_mpi: ', psum_mpi(1:a%surfnum)
      !!psum_mpi(:) = recvarray(:)
!
!     ! add contribution of these nodes to the surface
!     ! this sum is reset in elmgmr.f before looping 
!     ! through the boundary element blocks
      a%stb_pres(1:a%surfnum) = psum_mpi(:)
      !!write(*,*) 'stb_pres: ', this%stb_pres(1:this%surfnum)
      a%stb_pres(1:a%surfnum) = a%stb_pres(1:a%surfnum)/a%areas(1:a%surfnum)

      ! write(*,*) 'sumstb_pres: ',a%stb_pres(1:a%surfnum)

      end subroutine


      ! function to return stabilisation pressure for specific surface

      function getstb_pres(a, surfnum)      
      
      implicit none
      
      class(multidomaincontainer) :: a
      integer :: surfnum
      real*8  :: getstb_pres

      integer :: i

      do i = 1, a%surfnum
         if (a%surfids(i) .eq. surfnum) then
            getstb_pres = a%stb_pres(i)
!!            write(*,*) 'GETSTB_PRES = ', getstb_pres
            exit
         end if 
      end do 

      end function getstb_pres
!
!
!!c *** set stabilisation pressure at t_{n+alf_{i}}
!!c
!!      subroutine setstabpressure_mdc(this, npress, spress)
!!      implicit none
!!      class(multidomaincontainer) :: this     
!!      integer :: npress
!!      real*8 :: spress(0:maxsurf)
!!      this%stabilisationpressure(1:npress) = spress(1:npress)          
!!      end subroutine         
!
! ***************************************************
! *** feedback specific subroutines and functions ***
! ***************************************************
!
      subroutine initialise_feedback(a)
      implicit none
      class(feedback) :: a   
      a%active = .false.  ! initialise as false      
      return
      end subroutine

      function isactive_feedback(a) 
      implicit none
      class(feedback) :: a   
      logical :: isactive_feedback
      isactive_feedback = a%active
      return
      end function

      subroutine setcntlpntrs_feedback(a,b,hr_indx)
      implicit none
      class(feedback) :: a
      type(control) :: b
      integer :: hr_indx
      a%hr_ctrl%p => b%n1(hr_indx)
      end subroutine
!
! **************************************************
! *** control specific subroutines and functions ***
! **************************************************
!
      subroutine initialise_control(a,cdim,pdim)
      implicit none
      class(control) :: a   
      integer :: cdim, pdim, i
      a%dim = cdim
      allocate(a%n(cdim))
      allocate(a%n1(cdim))
      allocate(a%tau(cdim))
      allocate(a%params(cdim,pdim)) 
      do i = 1,cdim
         a%n(i) = real(1.0,8)
         a%n1(i) = real(1.0,8)
      end do
      return
      end subroutine
!
! *** set control at t = t_{n+alfi}
!
      subroutine set_control(a,b)           
      implicit none
      class(control) :: a
      type(feedback) :: b
      real*8 :: delt_tau
      real*8 :: odesource(a%dim)
      integer :: i
!
!     ! if feedback active
      if (b%active) then
!
!        ! get source
         odesource = a%getsource(b)
!         
!        ! first order finite difference      
         if (a%first_order_fd) then
!            
            do i = 1,a%dim 
               delt_tau = (alfi*delt)/a%tau(i)
               a%n1(i) = a%n(i) + delt_tau*odesource(i) 
               a%n1(i) = a%n1(i)/(real(1.0,8) + delt_tau)
            end do
!
!        ! second order finite difference
         elseif (a%second_order_fd) then
         end if
!         
      else
!
         do i = 1,a%dim      
            a%n1(i) = real(1.0,8)
         end do
!         
      end if
!
!
      return
      end subroutine
!
! *** update at t = t_{n+1}
!
      subroutine update_control(a,b)
      implicit none
      class(control) :: a
      type(feedback) :: b
      real*8 :: delt_tau
      real*8 :: odesource(a%dim)
      integer :: i
!
!
      if (b%active) then
!      
!        ! get source
         odesource = a%getsource(b)
!         
!        ! first order finite difference      
         if (a%first_order_fd) then
!            
            do i = 1,a%dim 
               delt_tau = (alfi*delt)/a%tau(i)
               a%n1(i) = a%n(i) + delt_tau*odesource(i) 
               a%n1(i) = a%n1(i)/(real(1.0,8) + delt_tau)
            end do
!            
         elseif (a%second_order_fd) then         
         end if
!         
      else
!
         do i = 1,a%dim      
            a%n1(i) = real(1.0,8) 
         end do
      end if
!
      do i = 1,a%dim
         a%n(i) = a%n1(i) !!! n+1 cycled to n 
      end do
!
!
      end subroutine
!
! ***********************************************************
! *** systemic circuit specific functions and subroutines ***
! ***********************************************************
!
! *** constructor 
!
      function sysconstructor(params) result(a)
      implicit none
      real*8 :: params(0:2)
      type(systemiccircuit) :: a
      call a%initialise_sys(params)
      return
      end function sysconstructor
!
! *** initialise 
!
      subroutine initialise_sys(this,params) 
      implicit none
      class(systemiccircuit) :: this
      real*8 :: params(0:2)
      integer :: fnum = 784
      integer :: gnum = 481
      integer :: i
      integer :: j
      integer :: ierr
      integer :: surfidsindx
      integer :: nvolumes, nflows, ninlets
      character(len=100) :: dimchar       
      type(rcrdata) :: rcrparam 
      integer :: intialvaluefile
      integer :: val 
      integer :: hstep

      ! A label so we can identify this class (avoids a world of pain when working with a tower of derived types)
      this%classNameString = 'systemiccircuit'
!
!     ! set as active
      this%isactive = int(1)
!
!     ! closed loop
      if (params(0) .gt. real(0.0,8)) then
         this%closedloop = .true.
      else
         this%closedloop = .false.
      end if      
!
!     !!! open loop preload
      !!this%patrial = params(2)
!
!     ! check if numerical rcr and numerical heart are active
      if (nrcr%isactive .and. hrt%isactive) then
      else
         write(*,*) '** WARNING numerical RCR and heart model not active'
         stop
      end if      
!      
!     ! set filenames for internal files
      this%xvarsfilename = 'sysXhist.dat'
      this%phistfilename = 'sysPhist.dat'
      this%qhistfilename = 'sysQhist.dat'
      this%chistfilename = 'sysChist.dat'
      this%fhistfilename = 'sysFhist.dat'
      this%ahistfilename = 'sysAhist.dat'
!
!     ! set filenames for restart files       
      this%flowfile = 'Qsystemic.dat'
      this%pressurefile = 'Psystemic.dat'
!
!     ! open file, should be openable as tested in initmultidomain
      open(fnum, file='sys.dat', status='old', iostat=ierr)
!
!     ! read switch header 
      read(fnum,*)
      read(fnum,*) 
!
!     ! read intial values flag
      read(fnum,*)
      read(fnum,*) this%readinitialvalues    
!
!     ! initialise feedback
      call this%feedback%initialise() 
!
!     ! initialise control
!     ! 1 = heart rate
!     ! 2 = maximum elastance
!     ! 3 = arterial resistance
!     ! 4 = venous compliance
!     ! 5 = venous unstressed volume
      call this%control%initialise(int(5),int(3))      
      this%control%getsource => getsource_sys  
!
!     ! set pointer to heart rate control
      call this%feedback%setcntlpntrs(this%control,this%hrindx)
!      
!     ! read number of feedback surfaces
      read(fnum,*)
      read(fnum,*) val     
!
      if (val .gt. int(0)) then
!
!        ! as the feedback class is not visible outside this module
         this%feedbackactive = .true. 
!         
         this%feedback%surfnum = val
!         
         allocate(this%feedback%surfids(val))
         allocate(this%feedback%pntr(val))
         allocate(this%feedback%params(2)) ! parameters set later          
!
         read(fnum,*)
         read(fnum,*) this%feedback%surfids(:)
!
!        ! set procedure pointers for feedback
         this%feedback%update => update_pressavg
         this%feedback%calculate => calculate_pressavg
!
      else
         read(fnum,*)
         read(fnum,*) 
      end if         
!
!     ! read number of rcr surfaces and add heart surface
      read(fnum,*)
      read(fnum,*) this%rcrsurfnum      
      if (this%closedloop) then
         this%surfnum = this%rcrsurfnum + int(1) ! includes heart surface
      else
         this%surfnum = this%rcrsurfnum ! without heart surface
      end if
      allocate(this%surfids(this%surfnum))   
      allocate(this%surfarea(this%surfnum))      
!
!     ! read rcr surface ids
      read(fnum,*)
      read(fnum,*) this%surfids(1:this%rcrsurfnum)
!      
!     ! add heart surface as last surface
      if (this%closedloop) then
         this%surfids(this%surfnum) = hrt%surfids(1)
      end if
!
!     ! read a3 compartment parameters
      read(fnum,*)
      read(fnum,*) this%ca3
      read(fnum,*) this%vua3 
      read(fnum,*) this%ra3
!
!     ! read v1 compartment parameters
      read(fnum,*)
      read(fnum,*) this%cv1
      read(fnum,*) this%vuv1
      read(fnum,*) this%rv1
!
!     ! read v2 compartment parameters
      read(fnum,*)
      read(fnum,*) this%cv2
      read(fnum,*) this%vuv2
      read(fnum,*) this%rv2
      read(fnum,*) this%lv2
!
!     ! read heart parameters      
!     ! read la parameters      
      read(fnum,*)
      read(fnum,*) this%ela
      read(fnum,*) this%vula
      read(fnum,*) this%rmv
      read(fnum,*) this%lmv
!         
      ! read lv parameters      
      read(fnum,*)
      read(fnum,*) this%vulv
      read(fnum,*) this%rav
      read(fnum,*) this%lav         
!
!     ! read elastance parameters      
      read(fnum,*)
      read(fnum,*) this%emax
      read(fnum,*) this%emin
      read(fnum,*) this%tmax
      read(fnum,*) this%trel
!
!     ! read heart period      
      read(fnum,*)
      read(fnum,*) this%period     
      if (this%feedbackactive) then    
         this%feedback%params(1) = this%period
      end if
!
!     ! save the tmax and trel ratios
      this%tmaxratio = this%tmax/this%period
      this%trelratio = this%trel/this%period
!
!     ! read control 
      if (this%feedbackactive) then
!
         ! activation tolerance 
         read(fnum,*)
         read(fnum,*) this%feedback%activation_tol
!
         ! nu steepness          
         read(fnum,*)
         read(fnum,*) this%feedback%params(2)

!        ! heart rate - tau, alpha, beta, gamma
         read(fnum,*)
         read(fnum,*) this%control%tau(this%hrindx)
         do i = 1,3
            read(fnum,*) this%control%params(this%hrindx,i)
         end do 
!
!        ! maximum elastance - tau, alpha, beta, gamma
         read(fnum,*)
         read(fnum,*) this%control%tau(this%emindx)
         do i = 1,3
            read(fnum,*) this%control%params(this%emindx,i)
         end do
!
!        ! arterial resistance - tau, alpha, beta, gamma
         read(fnum,*)
         read(fnum,*) this%control%tau(this%rdindx)
         do i = 1,3
            read(fnum,*) this%control%params(this%rdindx,i)
         end do 
!
!        ! venous compliance - tau, alpha, beta, gamma
         read(fnum,*)
         read(fnum,*) this%control%tau(this%cpindx)
         do i = 1,3
            read(fnum,*) this%control%params(this%cpindx,i)
         end do
!
!        ! venous unstressed volume - tau, alpha, beta, gamma
         read(fnum,*)
         read(fnum,*) this%control%tau(this%vuindx)
         do i = 1,3 
            read(fnum,*) this%control%params(this%vuindx,i)     
         end do
!         
      end if
!     
      ! close file
      close(fnum)
!             
!     ! number of volumes (capacitors) including heart chambers (va3 vv1 vv2 vla vlv)
      if (this%closedloop) then
         nvolumes = int(5)
      else
         nvolumes = int(3)
      end if
!
!     ! number of flows (inductors), including heart chamber (qv2 qmv)
      if (this%closedloop) then
         nflows = int(2) 
      else
         nflows = int(1) 
      end if
!
!     ! number of x variables
      this%xdim = int(2)*this%rcrsurfnum + nvolumes + nflows !!+ ninlets
!
!     ! allocate flow, etc arrays
      allocate(this%pressure_n(this%surfnum))             
      allocate(this%flow_n(this%surfnum))             ! flow values at t_{n}
      allocate(this%flow_n1(this%surfnum))             ! flow values at t_{n}
      allocate(this%rp(this%rcrsurfnum))              ! proximal resistance
      allocate(this%c(this%rcrsurfnum))               ! capacitance
      allocate(this%rd(this%rcrsurfnum))              ! distal resistance
      allocate(this%x_n(this%xdim))                   ! x-variables at t_{n}
      allocate(this%x_n1(this%xdim))                  ! x-variables at t_{n+alfi}
!
!     ! allocate history arrays
      if (lstep .gt. int(0)) then
         hstep = nstep + lstep
      else
         hstep = nstep
      end if
      allocate(this%x_hist(hstep+1,this%xdim))              ! x-variable history
      allocate(this%p_hist(hstep+1,this%xdim-nflows))       ! pressure history
      allocate(this%q_hist(hstep+1,this%surfnum))           ! flow history       
      allocate(this%elv_hist(hstep+1))                      ! elv history      
      allocate(this%r_hist(hstep+1,this%surfnum,2))         ! restart history
      allocate(this%control%hist(hstep+1,this%control%dim)) ! control history   
      allocate(this%act_hist(hstep+1))                    ! activation history
      allocate(this%stab_hist(hstep+1,this%surfnum))                    ! stabilisation history
!
!     ! allocate feedback arrays, here used for pressure
      if (this%feedbackactive) then
         allocate(this%feedback%input(hstep+1,this%feedback%surfnum))   ! feedback input - surface pressure
         allocate(this%feedback%output(hstep+1,this%feedback%surfnum))  ! feedback output - average pressure
         allocate(this%feedback%outputval(this%feedback%surfnum))       ! feedback value         
         allocate(this%feedback%targetval(this%feedback%surfnum))       ! target value         
         this%feedback%targetval(:) = real(0.0,8)                       ! zero target values
!!         allocate(this%avgpress(this%feedback%surfnum))
         allocate(this%feedback%tolhist(hstep+1,this%feedback%surfnum))
         this%feedback%tolhist(:,:) = real(0.0,8)
      end if
!
!     ! set up variable indices
!     !              1:rcrsurfnum = p_{rcrsurf}
!     ! rcrsurfnum+1:2*rcrsurfnum = pc_{rcrsurf}
!     !            2*rcrsurfnum+1 = va3
!     !            2*rcrsurfnum+2 = vv1
!     !            2*rcrsurfnum+3 = vv2
!     !            2*rcrsurfnum+4 = qv2
!     !            2*rcrsurfnum+5 = vla
!     !            2*rcrsurfnum+6 = qmv
!     !            2*rcrsurfnum+7 = vlv
!     !            2*rcrsurfnum+8 = part
      this%va3indx = int(2)*this%rcrsurfnum + int(1)
      this%vv1indx = int(2)*this%rcrsurfnum + int(2)
      this%vv2indx = int(2)*this%rcrsurfnum + int(3)
      this%qv2indx = int(2)*this%rcrsurfnum + int(4)
      this%vlaindx = int(2)*this%rcrsurfnum + int(5)
      this%vlvindx = int(2)*this%rcrsurfnum + int(6)      
      this%qmvindx = int(2)*this%rcrsurfnum + int(7)
!
!     ! zero arrays
      this%pressure_n(:) = real(0.0,8)
      this%flow_n(:) = real(0.0,8)
      this%flow_n1(:) = real(0.0,8)
      this%x_n(:) = real(0.0,8)
      this%x_n1(:) = real(0.0,8)      
      this%x_hist(:,:) = real(0.0,8)
      this%p_hist(:,:) = real(0.0,8)
      this%q_hist(:,:) = real(0.0,8)
      this%control%hist(:,:) = real(0.0,8)
      this%elv_hist(:) = real(0.0,8)  
      this%r_hist(:,:,:) = real(0.0,8)
!!      this%f_hist(:,:) = real(0.0,8)
!
!     ! initialise control as one
!!      this%c_n1(:) = real(1.0,8) 
!
!     ! intialise volumes 
      this%x_n(this%va3indx) = real(523.8990,8)
      this%x_n(this%vv1indx) = real(692.6520,8)
      this%x_n(this%vv2indx) = real(2329.564,8)
      if (this%closedloop) then
         this%x_n(this%vlaindx) = real(96.9333,8)
         this%x_n(this%vlvindx) = real(130,8)
      end if
!
!     ! array of flow pointers 
      allocate(this%flowpntr(this%surfnum))
      allocate(this%dflowpntr(this%surfnum))
      allocate(this%presspntr(this%surfnum))
      allocate(this%stabpntr(this%surfnum))
!      
!     ! load rcr parameters from nrcr, looping over the rcrsurfnum only
      do i = 1, this%rcrsurfnum
!!!         if (hassurfid(nrcr,this%surfids(i))) then                      
         ! is calling with the object OK? (see above ^^)
         if (nrcr%hassurfid(this%surfids(i))) then                               
            surfidsindx = getsurfidindx(nrcr,this%surfids(i))
            rcrparam = nrcr%rcrparams(surfidsindx)
            this%rp(i) = rcrparam%rp
            this%c(i) = rcrparam%c
            this%rd(i) = rcrparam%rd
         end if
      end do
!
!     ! implicit coefficients to replace nrcr and hrt values
      allocate(this%implicitcoeff(this%surfnum,2))
      this%implicitcoeff(:,:) = real(0.0,8)
!
!     ! write formats
!
!     ! x-variables
      write(dimchar,'(i10)') this%xdim
      write(this%xvarsformat,'(3(a))') &
      '(i8,',trim(adjustl(dimchar)),'(e20.10))'
!     ! flows
      write(dimchar,'(i10)') this%surfnum
      write(this%qvarsformat,'(3(a))') &
      '(i8,',trim(adjustl(dimchar)),'(e20.10))'
!     ! pressure
      write(dimchar,'(i10)') this%xdim - nflows
      write(this%pvarsformat,'(3(a))') &
      '(i8,',trim(adjustl(dimchar)),'(e20.10))'
!     ! control
      write(dimchar,'(i10)') this%control%dim
      write(this%cvarsformat,'(3(a))') &
      '(i8,',trim(adjustl(dimchar)),'(e20.10))'
!     ! restart
      write(dimchar,'(i10)') this%surfnum
      write(this%rvarsformat,'(3(a))') &
      '(',trim(adjustl(dimchar)),'(e20.10))'
!     ! feedback
      write(dimchar,'(i10)') int(3)*this%feedback%surfnum ! input(s), output(s) and target value(s)
      write(this%fvarsformat,'(3(a))') &
      '(i8,',trim(adjustl(dimchar)),'(e20.10))'  
      ! tolerance
      write(dimchar,'(i10)') this%feedback%surfnum 
      write(this%tvarsformat,'(3(a))') &
      '(i8,',trim(adjustl(dimchar)),'(e20.10))'  
!     
!     ! activation
      write(this%avarsformat,'(3(a))') '(i8,e20.10)' 
!      
      end subroutine initialise_sys
!
! *** set flow pointers to multidomain container
!
      subroutine setflowpntrs_sys(this, container)
      implicit none
!
      class(systemiccircuit) :: this
      type(multidomaincontainer) :: container
      integer :: i
      integer :: j
!
      do i = 1, this%surfnum
          do j = 1, container%surfnum
             if (this%surfids(i) .eq. container%surfids(j)) then
                this%flowpntr(i)%p => container%flow(j)            
             end if
         end do 
      end do
!
      end subroutine setflowpntrs_sys 

      subroutine setflowpntrs_coronary(this, container)
         implicit none
   !
         class(controlledCoronaryModel) :: this
         type(multidomaincontainer) :: container
         integer :: i
         integer :: j
#if EXTRA_CONSOLE_OUTPUT == 1
         write(*,*) 'numberOfControlledCoronaryModelSurfaces', this%numberOfControlledCoronaryModelSurfaces
#endif
         do i = 1, this%numberOfControlledCoronaryModelSurfaces
#if EXTRA_CONSOLE_OUTPUT == 1
         withrite(*,*) 'flowpntr', this%surfids(i) !\todo remove
#endif
             do j = 1, container%surfnum
                if (this%surfids(i) .eq. container%surfids(j)) then
                   this%flowpntr(i)%p => container%flow(j)            
                end if
            end do 
         end do
!
      end subroutine setflowpntrs_coronary

      subroutine setflowpntrs_netlistLPN(this,container)
         implicit none

         class(netlistLPN) :: this
         type(multidomaincontainer) :: container
         integer :: ii
         integer :: jj

         do ii = 1, this%numberOfLPNSurfaces
             do jj = 1, container%surfnum
                if (this%surfids(ii) .eq. container%surfids(jj)) then
                   this%flowpntr(ii)%p => container%flow(jj)            
                end if
            end do 
         end do

      end subroutine setflowpntrs_netlistLPN
!
! *** set pressure pointers to multidomain container
!
      subroutine setpresspntrs_sys(this, container)
      implicit none
!
      class(systemiccircuit) :: this
      type(multidomaincontainer) :: container
      integer :: i, j
      do i = 1, this%surfnum
          do j = 1, container%surfnum
             if (this%surfids(i) .eq. container%surfids(j)) then
                this%presspntr(i)%p => container%pressure(j)            
                this%stabpntr(i)%p => container%stb_pres(j)           
             end if
         end do 
      end do
      end subroutine setpresspntrs_sys
!
! *** set feedback pointers to pressure pointers in multidomain container
!
      subroutine setfeedbackpntrs_sys(this, container)
      implicit none
      class(systemiccircuit) :: this
      type(multidomaincontainer) :: container
      integer :: i, j
      do i = 1, this%feedback%surfnum
          do j = 1, container%surfnum
             if (this%feedback%surfids(i) .eq. container%surfids(j)) then
                this%feedback%pntr(i)%p => container%pressure(j)            
             end if
         end do 
      end do
      end subroutine setfeedbackpntrs_sys

!
! *** initialise internal variables
!
      subroutine initxvars_sys(a,stepn) 
!      
      implicit none
      class(systemiccircuit) :: a
      integer :: stepn
      integer :: i ,j 
      integer :: xnum = 387
      integer :: fnum = 652
      integer :: rnum = 743
      integer :: qnum = 573
      integer :: pnum = 983
      integer :: enum = 772
      integer :: cnum = 734
      integer :: anum = 377
      real*8 :: temp(a%xdim+1) ! includes step number
      real*8 :: dt_n, elv_n, cntl, cntlv(a%control%dim)
      integer :: ierr

!
      if (stepn .eq. int(0)) then
!
!        ! set initial values of flow for rcr surfaces
         do i = 1, a%rcrsurfnum
            a%flow_n(i) = a%flowpntr(i)%p         
         end do
!
!        ! set rcr surface pressures and capacitor pressures from flow
         do i = 1, a%rcrsurfnum
            a%x_n(i) = a%presspntr(i)%p                     ! set pressure at interface
            a%x_n(a%rcrsurfnum+i) = a%presspntr(i)%p      & ! set capacitor pressure using initial flow
                                  - a%rp(i)*a%flowpntr(i)%p
         end do
!
!        ! initialise heart surface
         if (a%closedloop) then
            a%flow_n(a%surfnum) = real(0.0,8)
            a%paorta_n = a%presspntr(a%surfnum)%p
            a%qaorta_n = a%flow_n(a%surfnum)
         end if
!
!        ! open initial values file
         if (a%readinitialvalues .gt. int(0)) then
            open(xnum, file='sysXinit.dat', status='old', iostat=ierr)         
            if (ierr .ne. int(0)) then
               write(*,*) ' ** Error: No sysXinit.dat file'
               close(xnum)
               stop
            else
!              ! read header
               read(xnum,*) 
               read(xnum,*) temp
               close(xnum)
               do i = 1,a%rcrsurfnum
                  a%x_n(i) = temp(i)
                  j = a%rcrsurfnum + i
                  a%x_n(j) = temp(j)
               end do
               a%x_n(a%va3indx) = temp(a%va3indx)
               a%x_n(a%vv1indx) = temp(a%vv1indx)
               a%x_n(a%vv2indx) = temp(a%vv2indx)
               a%x_n(a%qv2indx) = temp(a%qv2indx)
               if (a%closedloop) then
                  a%x_n(a%vlaindx) = temp(a%vlaindx)
                  a%x_n(a%qmvindx) = temp(a%qmvindx)
                  a%x_n(a%vlvindx) = temp(a%vlvindx)
               end if
            end if 
         end if 
!
!        ! set initial average pressure to value in pressure pointers
         if (a%feedbackactive) then
            do i = 1, a%feedback%surfnum
!!               a%avgpress(i) = a%feedback%pntr(i)%p
               a%feedback%outputval = a%feedback%pntr(i)%p
            end do 
         end if
!
!        ! initialise control
         do i = 1, a%control%dim
            a%control%n(i) = real(1.0,8)
         end do
!
!        ! initialise activation time to zero
         a%activationtime = real(0.0,8)
!         
      else
!
!        ! open flow restart history
         open(rnum, file=a%flowfile, status='old', iostat=ierr)               
         do i = 1, stepn
            read(rnum,*) a%r_hist(i,:,1)
         end do
         close(rnum)
         a%flow_n(:) = a%r_hist(stepn,:,1)
         a%qaorta_n = a%flow_n(a%surfnum)
!
!        ! open pressure restart history
         open(rnum, file=a%pressurefile, status='old', iostat=ierr)               
         do i = 1, stepn
            read(rnum,*) a%r_hist(i,:,2)
         end do
         close(rnum)
         a%paorta_n = a%r_hist(stepn,a%surfnum,2)
!
!        ! open x-vars history
         open(xnum, file='sysXhist.dat', status='old', iostat=ierr)               
         do i = 1, stepn
            read(xnum,*) temp
            a%x_hist(i,:) = temp(2:a%xdim+1)
         end do 
         close(xnum)
!
!       ! set x variables at t_{n}         
         do i = 1, a%rcrsurfnum
            a%x_n(i) = temp(i+1) ! step number stored in first entry
            j = a%rcrsurfnum + i 
            a%x_n(j) = temp(j+1)
         end do
         a%x_n(a%va3indx) = temp(a%va3indx+1)
         a%x_n(a%vv1indx) = temp(a%vv1indx+1)
         a%x_n(a%vv2indx) = temp(a%vv2indx+1)
         a%x_n(a%qv2indx) = temp(a%qv2indx+1)
         if (a%closedloop) then
            a%x_n(a%vlaindx) = temp(a%vlaindx+1)
            a%x_n(a%qmvindx) = temp(a%qmvindx+1)
            a%x_n(a%vlvindx) = temp(a%vlvindx+1)         
         end if
!
!        ! read flow history
         open(qnum, file=a%qhistfilename, status='old', iostat=ierr)
         do i = 1, stepn
            read(qnum,*) j, a%q_hist(i,:)
         end do
         close(qnum) 
!
!        ! read pressure history
         open(pnum, file=a%phistfilename, status='old', iostat=ierr)
         do i = 1, stepn
            read(pnum,*) j, a%p_hist(i,:)            
         end do
         close(pnum) 
!
!        ! read control history
         open(cnum, file=a%chistfilename, status='old', iostat=ierr)
         do i = 1, stepn
            read(cnum,*) j, a%control%hist(i,:)            
         end do
         close(cnum)  
         do i = 1, a%control%dim
            a%control%n(i) = a%control%hist(stepn,i)            
         end do
!
!        ! check to see if control is activated    
         if (a%feedbackactive) then         
            ! normalise by x = {1,1,1,1,1}
            cntlv = a%control%n(:)/sqrt(real(5,8))
            ! determine length of control vector
!            cntl = norm2(cntlv)   
            cntl = dot_product(cntlv,cntlv)
            cntl = sqrt(cntl)            
            ! subtract 1
            cntl = cntl - real(1.0,8)
            if ( abs(cntl) .gt. real(1.0e-3,8)) then                  
               a%feedback%active = .true.
            end if
         end if 
!
!        ! read feedback history
         if (a%feedbackactive) then            
            open(fnum, file=a%fhistfilename, status='old', iostat=ierr)
            do i = 1, stepn
               read(fnum,*) j, a%feedback%input(i,:), &
                               a%feedback%output(i,:), &
                               a%feedback%targetval(:)
            end do
            close(fnum)   
!
!           ! set feedback value (average pressure) at stepn from file            
!!            i = a%feedback%surfnum + int(1)
!!            j = int(2)*a%feedback%surfnum
!!            a%avgpress(:) = a%feedback%hist(stepn,i:j)
            a%feedback%outputval(:) = a%feedback%output(stepn,:)            

         end if
!
!        ! read elastance history
         open(enum, file='sysEhist.dat', status='old', iostat=ierr)
         do i = 1, stepn
            read(enum,*) j, a%elv_hist(i)
         end do
         close(enum) 
!
!        ! read activation history
         open(anum, file=a%ahistfilename, status='old')      
         do i = 1, stepn
            read(anum,*) j, a%act_hist(i)
         end do
         a%activationtime = a%act_hist(stepn)
         close(anum)
!         
      end if
!
      if (a%closedloop) then
!      
!        ! delta time to be at t =  t_{n} 
         dt_n = real(0.0,8)
!
!        ! elv and plv at t_{n}, uses current activation time
         elv_n = a%getelastance(dt_n,a%control%n)
         a%plv_n = elv_n*(a%x_n(a%vlvindx) - a%vulv)
!
!        ! pla at t_{n}
         a%pla_n = a%ela*(a%x_n(a%vlaindx) - a%vula)      
!
!        ! set av state
         if (a%plv_n .gt. a%paorta_n) then
            hrt%avopen = int(1)             
            a%updatepressure = int(0)
            a%qaorta_n = a%flow_n(a%surfnum)
         else if (stepn .gt. 0 .and. a%qaorta_n .lt. real(0.0,8)) then !! flow is negative !!
            hrt%avopen = int(1)
            a%updatepressure = int(0)             
         else
            hrt%avopen = int(0)             
            a%updatepressure = int(1)
            a%qaorta_n = real(0.0,8)      
         end if
!         
      end if
!
      return
      end subroutine initxvars_sys
!
! *** get elastance function wrapper
!     note that the nondimensionalised periodic time is stored in the 
!     systemic circuit class, so here the delta time to the required
!     elastance is provided
!
      function getelastance_sys(a,delta_time,ctrlpntr) result(elastance)
      implicit none
      class(systemiccircuit) :: a
      real*8 :: elastance
      real*8 :: delta_time, periodtime
      real*8 :: ctrlpntr(a%control%dim)
      real*8 :: heartrate, period, emax 
!
!     ! current heart rate in control(1)
      if (a%feedback%active) then
         heartrate = ctrlpntr(a%hrindx)/a%period         
      else
         heartrate = real(1.0,8)/a%period
      end if
      period = real(1.0,8)/heartrate
!
!     ! current max elastance in control(2)
      if (a%feedback%active) then
         emax = ctrlpntr(a%emindx)*a%emax
      else
         emax = a%emax
      end if
!
!!c     ! non-dim relaxation time based on heart rate
!!c     ! values after rideout, via ottensen 
!!      trel = real(0.29,8)*heartrate - real(0.18,8)
!
!
!     ! using the nondim activation time [0->1] add dt/period
      periodtime = a%activationtime + delta_time/period
      ! redimensionalise
      periodtime = periodtime*period
!
!     ! get elastance function
      elastance = getelastance(periodtime, period, emax, &
                               a%emin, a%tmaxratio, a%trelratio)
!      
      return
      end function getelastance_sys
!
!
! *** analytical elastance function from:
!     pope, s. r.; ellwein, l. m.; zapata, c. l.; novak, v.; kelley, c. t. & olufsen, m. s.  
!     estimation and identification of parameters in a lumped cerebrovascular model.
!     math biosci eng, 2009, 6, 93-115
!
      function getelastance(periodtime, period, emax, emin, tmaxval, trelval)
      implicit none                       
      real*8 :: periodtime, period, emax, emin, tmaxval, trelval
      real*8 :: tper, tmax, trel
      real*8 :: evalue, getelastance
!!      integer :: nper
!      
!
      tper = period
      tmax = tmaxval*tper
      trel = trelval*tper
!                             
      if ( periodtime .le. tmax ) then
         evalue = emin &
                + real(0.5,8)*(emax - emin) &
                * (real(1.0,8) - cos((periodtime*pi)/tmax))
      elseif ( periodtime .gt. tmax .and. periodtime .le. (tmax + trel) ) then
         evalue = emin &
                + real(0.5,8)*(emax-emin) &
                * (real(1.0,8) + cos((periodtime-tmax)*(pi/trel)))            
      elseif ( periodtime .gt. (tmax + trel) ) then            
         evalue = emin
      end if
!
      getelastance = evalue

      !\todo remove
#if EXTRA_CONSOLE_OUTPUT == 1      
      write(*,*) 'elastance info:', emax
      write(*,*) 'tmax', tmaxval
      write(*,*) 'trelax', trelval
#endif      
!
      return      
      end function
!
! *** senzaki elastance function
!      
      function getelastance_senzaki(a,periodtime, period, emax, emin, tmaxval, &
                                    trelval) result(evalue)
!
      class(numericalheart) :: a
      real*8 :: periodtime, period, emax, emin, tmaxval, trelval
      real*8 :: tper, tmax, trel, evalue, etime 
!
      tper = period
      tmax = tmaxval*tper
      trel = trelval*tper
!
!     ! normalised elastance time 
      etime = periodtime/tmax
!
      if ( periodtime .le. (tmax + trel) ) then         
         evalue = emin + (emax - emin)*getvalue(etime,a%elv_senzaki)
      else
         evalue = emin
      end if
!
      return
      end function
!      
      subroutine set_sPress(a, press)
      class(numericalheart) :: a
      real*8 :: press(0:maxsurf)
      a%sPress = press(1)
      end subroutine

!!      subroutine set_dsPress(a, press)
!!      class(numericalheart) :: a
!!      real*8 :: press(0:maxsurf)
!!      a%dsPress = press(1)
!!      end subroutine
!
! *** set stabilisation pressure pointer to multidomain container
!
      subroutine setpresspntrs_hrt(this, container)        
      implicit none
!
      class(numericalheart) :: this
      type(multidomaincontainer) :: container
      integer :: i, j
      do i = 1, this%surfnum
          do j = 1, container%surfnum
             if (this%surfids(i) .eq. container%surfids(j)) then
                this%sPress => container%stb_pres(j) 
!!                write(*,'(a,i5,a,i5)') '*** this surfid: ',this%surfids(i), &
!!                                      ' pointed to: ',container%surfids(j)
             end if
         end do 
      end do
      end subroutine 
!
!
! *** solve open loop systemic circuit
!
      subroutine solve_sys_open(a,stepn)
!
      implicit none
      class(systemiccircuit) :: a
      integer :: stepn
      integer :: i
      real*8 :: a_lhs(a%xdim,a%xdim)
      real*8 :: a_rhs(a%xdim,a%xdim)
      real*8 :: b_rhs(a%xdim,a%xdim)
      real*8 :: c_rhs(a%xdim)
      real*8 :: d_rhs(a%xdim,a%rcrsurfnum)            
      real*8 :: dsegv_a(a%xdim,a%xdim)
      real*8 :: dsegv_b(a%xdim)
      real*8 :: tmpvec(a%xdim)                    
      real*8 :: tmpmatx(a%xdim,a%rcrsurfnum)
      real*8 :: tmpvecx(a%xdim)  
      real*8 :: tmpsol(a%xdim)  
      real*8 :: flowsurf(a%rcrsurfnum)
      real*8 :: dt_alfi              
      real*8 :: amat(2,2), bvec(2), cvec(2), xvec(2)
      real*8 :: t_n
      integer :: nsteps, startindx, endindx, counter      
!  
!
!       ! NO CONTROL AT THE MOMENT
!!!c     ! set control at t = t_{n+alfi}
!!!      if (a%feedbackactive) then
!!!         call a%control%set(a%feedback)
!!!      end if                 
!
!     ! set system to be solved
      call a%setsystem_sys_open(a_lhs, a_rhs, b_rhs, c_rhs, d_rhs, 'solve') 
!
!     ! dsegv matrix a
      dsegv_a(:,:) = a_lhs(:,:) - a_rhs(:,:)
!
!     ! solve for state variables at n+1
!
!     ! dsegv vector b => b_rhs * xn + c_rhs
      dsegv_b = multiply_ab(b_rhs, a%x_n) 
      dsegv_b = dsegv_b + c_rhs
!
!     ! solve for the rcr non-flow implicit coefficient
      tmpvecx = solve_axb(dsegv_a, dsegv_b)
!
      a%implicitcoeff(1:a%rcrsurfnum,2) = tmpvecx(1:a%rcrsurfnum)               
!
!     ! solve for the flow implicit coefficient
      tmpmatx = solve_axb(dsegv_a, d_rhs)        
!
!     ! here each rcr flow only interacts with its own surface 
      do i = 1, a%rcrsurfnum
         flowsurf(:) = real(0.0,8)
         flowsurf(i) = real(1.0,8)
         tmpvecx = multiply_ab(tmpmatx,flowsurf)
         a%implicitcoeff(i,1) = tmpvecx(i)
      end do
!
!     ! set implicit coefficients in the nrcr object
      call a%setimplicitcoeff_sys(nrcr)    
!
      end subroutine       
!
! *** solve systemiccircuit
!
      subroutine solve_sys(a,stepn)
!
      implicit none
      class(systemiccircuit) :: a
      integer :: stepn
      integer :: i
      real*8 :: a_lhs(a%xdim,a%xdim)
      real*8 :: a_rhs(a%xdim,a%xdim)
      real*8 :: b_rhs(a%xdim,a%xdim)
      real*8 :: c_rhs(a%xdim)
      real*8 :: d_rhs(a%xdim,a%rcrsurfnum)            
      real*8 :: dsegv_a(a%xdim,a%xdim)
      real*8 :: dsegv_b(a%xdim)
      real*8 :: tmpvec(a%xdim)                    
      real*8 :: tmpmatx(a%xdim,a%rcrsurfnum)
      real*8 :: tmpvecx(a%xdim)  
      real*8 :: tmpsol(a%xdim)  
      real*8 :: flowsurf(a%rcrsurfnum)
      real*8 :: dt_alfi              
      real*8 :: amat(2,2), bvec(2), cvec(2), xvec(2)
      real*8 :: t_n
      integer :: nsteps, startindx, endindx, counter
      real*8 :: xiter(2), j(2,2), r(2), dx(2)    
      real*8 :: dtime, eval  
      real*8, parameter :: tol = real(1d-6,8)
      integer, parameter :: itermax = 10
!
!     ! set mv valve state based on pressures at t_{n}
!     ! if pla > plv, mv open
!     ! elseif qmv +ve, mv open
!     ! else mv closed
      !!if (a%x_n(a%qmvindx) .lt. real(0.0,8)) then
      !!   a%mvopen = .false.
      if (a%pla_n .gt. a%plv_n) then                           
         a%mvopen = .true.                                     
      elseif (a%x_n(a%qmvindx) .gt. real(0.0,8)) then           
         a%mvopen = .true.                                        
      else                                                     
         a%mvopen = .false.                                    
      end if                                                   
!     
!     ! set av valve state based on pressure at t_{n}, if plv > paorta then av open,
!     ! elseif qaorta -ve then av open, else av closed
      if (a%plv_n .gt. a%paorta_n) then                            
         a%avopen = .true.                                     
         hrt%avopen = int(1) 
         a%updatepressure = int(0)
      elseif (a%qaorta_n .lt. real(0.0,8)) then   
         a%avopen = .true.                                     
         hrt%avopen = int(1) 
         a%updatepressure = int(0)
      else                                                     
         a%avopen = .false.                                    
         hrt%avopen = int(0) 
         a%updatepressure = int(1)
      end if     
!
!     ! delta time to t_{n+alfi}
      dt_alfi = alfi*delt

!     ! set control at t = t_{n+alfi}
      if (a%feedbackactive) then
         call a%control%set(a%feedback)
      end if      
!      
!     ! set elastance at t = t_{n+alfi}
      a%elv_n1 = a%getelastance(dt_alfi,a%control%n1)                
!
!     ! set system to be solved
      call a%setsystem_sys(a_lhs, a_rhs, b_rhs, c_rhs, d_rhs, 'solve') 
!
!     ! dsegv matrix a
      dsegv_a(:,:) = a_lhs(:,:) - a_rhs(:,:)
!
!     ! solve for state variables at n+1
!
!     ! dsegv vector b => b_rhs * xn + c_rhs
      if (a%mvopen) then
         dsegv_b(:) = multiply_ab(b_rhs, a%x_n) 
         dsegv_b(:) = dsegv_b(:) + c_rhs(:)
      elseif (a%avopen) then
         dsegv_b(1:a%vlaindx) = multiply_ab(b_rhs(1:a%vlaindx,1:a%vlaindx), &
                                            a%x_n(1:a%vlaindx)) 
         dsegv_b(1:a%vlaindx) = dsegv_b(1:a%vlaindx) + c_rhs(1:a%vlaindx)      
      else
         dsegv_b(1:a%vlvindx) = multiply_ab(b_rhs(1:a%vlvindx,1:a%vlvindx), &
                                            a%x_n(1:a%vlvindx)) 
         dsegv_b(1:a%vlvindx) = dsegv_b(1:a%vlvindx) + c_rhs(1:a%vlvindx)
      end if
!
!     ! solve for the rcr non-flow implicit coefficient
      if (a%mvopen) then
         tmpvecx = solve_axb(dsegv_a, dsegv_b)
      else if (a%avopen) then
         tmpvecx(1:a%vlaindx) = solve_axb(dsegv_a(1:a%vlaindx,1:a%vlaindx), &
                                          dsegv_b(1:a%vlaindx))      
      else
         tmpvecx(1:a%vlvindx) = solve_axb(dsegv_a(1:a%vlvindx,1:a%vlvindx), &
                                          dsegv_b(1:a%vlvindx))
      end if
      a%implicitcoeff(1:a%rcrsurfnum,2) = tmpvecx(1:a%rcrsurfnum)               
!
!     ! solve for the flow implicit coefficient
      if (a%mvopen) then
         tmpmatx = solve_axb(dsegv_a, d_rhs)        
      else if (a%avopen) then
         tmpmatx(1:a%vlaindx,:) = solve_axb(dsegv_a(1:a%vlaindx,1:a%vlaindx), & 
                                            d_rhs(1:a%vlaindx,:))              
      else
         tmpmatx(1:a%vlvindx,:) = solve_axb(dsegv_a(1:a%vlvindx,1:a%vlvindx), &
                                            d_rhs(1:a%vlvindx,:))        
      end if
!
!     ! here each rcr flow only interacts with its own surface 
      do i = 1, a%rcrsurfnum
         flowsurf(:) = real(0.0,8)
         flowsurf(i) = real(1.0,8)
         tmpvecx = multiply_ab(tmpmatx,flowsurf)
         a%implicitcoeff(i,1) = tmpvecx(i)
      end do
!
!     ! set implicit coefficients in the nrcr object
      call a%setimplicitcoeff_sys(nrcr)
!
!!c     ! if valve is open then solve for av flow
!!      if (a%avopen) then
!!
!!         amat(1,1) = real(1.0,8)
!!         amat(1,2) = real(-1.0,8)*a%elv_n1 
!!         amat(2,1) = real(0.0,8)
!!         amat(2,2) = real(1.0,8)/(alfi*delt)
!!c
!!c        ! nb aortic flow is -ve on the inflow surface
!!         bvec(1) = (a%rav+a%lav/(alfi*delt))
!!         bvec(2) = real(1.0,8)
!!
!!         cvec(1) = (real(-1.0,8)*a%lav*a%qaorta_n)/(alfi*delt) 
!!     &           - a%elv_n1*a%vulv
!!         cvec(2) = a%x_n(a%vlvindx)/(alfi*delt) 
!!
!!         xvec = solve_axb(amat, bvec)
!!
!!         a%implicitcoeff(a%surfnum,1) = xvec(1)
!!         a%vlvcoeff(1) = xvec(2)
!!
!!         xvec = solve_axb(amat, cvec)
!!
!!         a%implicitcoeff(a%surfnum,2) = xvec(1)
!!         a%vlvcoeff(2) = xvec(2)
!!
!!         call a%setimplicitcoeff_sys(hrt)            
!!      end if
!
!     ! if valve open
      if (a%avopen) then
!      
!        ! initial pressure and volume at previous step
         xiter(1) = a%paorta_n
         xiter(2) = a%x_n(a%vlvindx)
!
!        ! solve timestep        
         dtime = alfi*delt        

!        ! jacobian
         j(:,:) = real(0.0,8)
         j(1,1) = real(1.0,8)
         j(1,2) = real(-1.0,8)*a%elv_n1*(real(1.0,8) + a%flowpntr(a%surfnum)%p) ! flow -ve
         j(2,2) = real(1.0,8)/dtime
!!         j(1,1) = real(1.0,8)
!!         j(2,2) = real(-1.0,8)*a%elv_n1
!!         j(2,2) = real(1.0,8)/dtime
!
!        ! initial value and evaluation
         eval = real(1e20,8)
!
         do i = 1,itermax
!
            r(1) = xiter(1) - a%stabpntr(a%surfnum)%p &
                 - a%elv_n1*(xiter(2) - a%vulv) &
                 *(real(1.0,8) + a%kelv*a%flowpntr(a%surfnum)%p) &! flow -ve
                 - a%rav*a%flowpntr(a%surfnum)%p &
                 + (a%lav/dtime)*(a%qaorta_n - a%flowpntr(a%surfnum)%p)

            r(2) = xiter(2)/dtime - a%x_n(a%vlvindx)/dtime &
                 - a%flowpntr(a%surfnum)%p
!
!
!!            r(1) = xiter(1) - a%elv_n1*(xiter(2) - a%vulv) 
!!     &           - a%rav*a%flowpntr(a%surfnum)%p 
!!     &           - (a%lav/dtime)*(a%flowpntr(a%surfnum)%p -a%qaorta_n)
!!            r(2) = xiter(2)/dtime - a%x_n(a%vlvindx)/dtime 
!!     &           - a%flowpntr(a%surfnum)%p
!
            dx = solve_axb(j,real(-1.0,8)*r)  
            xiter = xiter + dx   
!
!            eval = norm2(dx) 
			eval = dot_product(dx,dx)
            eval = sqrt(eval) 
!
!           ! if evaluation is less than tolerance exit
            if (eval .lt. tol) then
               exit
            end if
!         
         end do
!    
!        ! implicit coefficients
         a%implicitcoeff(a%surfnum,1) = a%rav + a%lav/dtime &
                                      + a%elv_n1*(xiter(2) - a%vulv)*a%kelv
         a%implicitcoeff(a%surfnum,2) = a%elv_n1*(xiter(2) - a%vulv) &
                                      - (a%lav/dtime)*a%qaorta_n
!
         call a%setimplicitcoeff_sys(hrt)            
!
      end if      
!      
      end subroutine             
!
! *** get source term for control ODE - here this the activity function 
!     based on the feedback at step t = t_{n} within the feedback class
!
      function getsource_sys(a,b) result(c)
      implicit none
      class(control) :: a
      type(feedback) :: b
      real*8 :: c(a%dim)
      real*8 :: diff(b%surfnum) 
      integer :: maxindx
      real*8 :: nsympathetic, nvagal
!
!     ! find maximum difference between current each branch's output and target
!     ! maximum difference used to drive the autoregulation
      diff(:) = b%outputval(:) - b%targetval(:)
!     ! absolute value taken to account for both min and max      
      diff(:) = abs(diff)
      maxindx = maxloc(diff,1)
!
!     ! parasympathetic
      !!nvagal = (b%targetval/b%outputval(maxindx))**b%params(2)      
      nvagal = (b%targetval(maxindx)/b%outputval(maxindx))**b%params(2)      
      nvagal = nvagal + real(1.0,8)
      nvagal = real(1.0,8)/nvagal
!      
!     ! sympathetic (inverse of above)
      !!nsympathetic = (b%outputval(maxindx)/b%targetval)**b%params(2)
      nsympathetic = (b%outputval(maxindx)/b%targetval(maxindx))**b%params(2)
      nsympathetic = nsympathetic + real(1.0,8)
      nsympathetic = real(1.0,8)/nsympathetic
!
!     ! source
      c(:) = a%params(:,1)*nsympathetic  &! alpha x sympathetic   
           - a%params(:,2)*nvagal        &! beta x parasympathetic
           + a%params(:,3)                ! gamma
!     
      return
      end function
!
! *** set up matrices for the simultaneous solve at t = t_{n+alfi}
!
!     a_lhs * x_{n+1} = a_rhs * x_{n+1} + b_rhs * x_{n} + c_rhs + d_rhs * q
!
!     a_lhs - x_{n+alfi} coefficients [xdmin,xdim]
!     a_rhs - x_{n+alfi} coefficients [xdim,xdim]
!     b_rhs - x_{n} coefficients [xdim,xdim]
!     c_rhs - constant coefficients [xdim]
!     d_rhs - flow_{n+alfi} coefficients [xdim,rcrflowsurfs]
!
      subroutine setsystem_sys_open(a, a_lhs, a_rhs, b_rhs, c_rhs, d_rhs, &
                                    charvar)
!      
      implicit none
      class(systemiccircuit) :: a
      real*8, intent(out) :: a_lhs(a%xdim,a%xdim)
      real*8, intent(out) :: a_rhs(a%xdim,a%xdim)
      real*8, intent(out) :: b_rhs(a%xdim,a%xdim)
      real*8, intent(out) :: c_rhs(a%xdim)
      real*8, intent(out) :: d_rhs(a%xdim,a%rcrsurfnum)            
      character(len=*) :: charvar
      character(len=*), parameter :: solve = 'solve'
      character(len=*), parameter :: update = 'update'
      real*8 :: rd_ctrl, cv_ctrl, vu_ctrl
      real*8 :: alfi_delt 
      integer :: i, j, nsrf      
!
!     ! timestep size
      if (charvar .eq. solve) then
         alfi_delt = real(alfi*delt,8)
      elseif (charvar .eq. update) then
         alfi_delt = delt
      end if
!
!     ! number of rcr surfaces
      nsrf = a%rcrsurfnum 
!
!     ! zero matrices
      a_lhs(:,:) = real(0.0,8)
      a_rhs(:,:) = real(0.0,8)
      b_rhs(:,:) = real(0.0,8)
      c_rhs(:) = real(0.0,8)
      d_rhs(:,:) = real(0.0,8)
!
!     ! SET TO UNITY
!     ! control at t_{n+alfi} or t_{n+1} depending upon solve/update
      rd_ctrl = real(1.0,8) !a%control%n1(a%rdindx)    
      cv_ctrl = real(1.0,8) !a%control%n1(a%cpindx)    
      vu_ctrl = real(1.0,8) !a%control%n1(a%vuindx)     
!
!     ! rcr surfaces
      do i = 1, nsrf 
!
!        ! pc index
         j = i + nsrf      
!
         ! p - pc = r q
         a_lhs(i,i) = real(1,8)
         a_lhs(i,j) = real(-1,8)
         d_rhs(i,i) = a%rp(i)
!        
         ! c dpc/dt = q - (pc - pa3)/rd 
         ! convert pa3 to va3 
         ! pa3 = (va3-vua3)/ca3        
         a_lhs(j,j) = a%c(i)/alfi_delt
         b_rhs(j,j) = a%c(i)/alfi_delt
!        
!        ! rd
         a_rhs(j,j) = real(-1,8)/(rd_ctrl*a%rd(i))
         a_rhs(j,a%va3indx) = real(1,8)/(rd_ctrl*a%rd(i)*a%ca3)
         c_rhs(j) = -a%vua3/(rd_ctrl*a%rd(i)*a%ca3)         
!                        
         d_rhs(j,i) = real(1,8)        
!        
        ! dva3/dt = q(1) + ... + q(nbranch)
         a_rhs(a%va3indx,j) = real(1,8)/(rd_ctrl*a%rd(i))
         a_rhs(a%va3indx,a%va3indx) = a_rhs(a%va3indx,a%va3indx) &
                                    - real(1,8)/(rd_ctrl*a%rd(i)*a%ca3)  
         c_rhs(a%va3indx) = c_rhs(a%va3indx) &
                          + a%vua3/(rd_ctrl*a%rd(i)*a%ca3)
!    
      end do
!
!     ! time derivatives - q_{i} equations are divided by l_{i} 
      do i = a%va3indx, a%qv2indx
         a_lhs(i,i) = real(1,8)/alfi_delt
         b_rhs(i,i) = real(1,8)/alfi_delt
      end do
!
!     ! va3
      a_rhs(a%va3indx,a%va3indx) = a_rhs(a%va3indx,a%va3indx) &
                                 - real(1,8)/(rd_ctrl*a%ra3*a%ca3)
      a_rhs(a%va3indx,a%vv1indx) = real(1,8)/(rd_ctrl*a%ra3*cv_ctrl*a%cv1)
      c_rhs(a%va3indx) = c_rhs(a%va3indx) &
                       + a%vua3/(rd_ctrl*a%ra3*a%ca3) &
                       - (vu_ctrl*a%vuv1)/(rd_ctrl*a%ra3*cv_ctrl*a%cv1)
!        
!     ! vv1
      a_rhs(a%vv1indx,a%va3indx) = real(1,8)/(rd_ctrl*a%ra3*a%ca3)
      a_rhs(a%vv1indx,a%vv1indx) = real(-1,8)/(rd_ctrl*a%ra3*cv_ctrl*a%cv1) &
                                 + real(-1,8)/(a%rv1*cv_ctrl*a%cv1)
      a_rhs(a%vv1indx,a%vv2indx) = real(1,8)/(a%rv1*cv_ctrl*a%cv2)
      c_rhs(a%vv1indx) = (vu_ctrl*a%vuv1)/(rd_ctrl*a%ra3*cv_ctrl*a%cv1)     &
                       + (vu_ctrl*a%vuv1)/(a%rv1*cv_ctrl*a%cv1)             &
                       - a%vua3/(rd_ctrl*a%ra3*a%ca3)                       &
                       - (vu_ctrl*a%vuv2)/(a%rv1*cv_ctrl*a%cv2) 
!        
!     ! vv2
      a_rhs(a%vv2indx,a%vv1indx) = real(1,8)/(a%rv1*cv_ctrl*a%cv1)
      a_rhs(a%vv2indx,a%vv2indx) = real(-1,8)/(a%rv1*cv_ctrl*a%cv2)
      a_rhs(a%vv2indx,a%qv2indx) = real(-1,8)
      c_rhs(a%vv2indx) = (vu_ctrl*a%vuv2)/(a%rv1*cv_ctrl*a%cv2) &
                       - (vu_ctrl*a%vuv1)/(a%rv1*cv_ctrl*a%cv1)
!    
!     ! qv2
      a_rhs(a%qv2indx,a%vv2indx) = real(1,8)/(a%lv2*cv_ctrl*a%cv2)
      a_rhs(a%qv2indx,a%qv2indx) = real(-1,8)*a%rv2/a%lv2
      c_rhs(a%qv2indx) = real(-1,8)*(vu_ctrl*a%vuv2)/(a%lv2*cv_ctrl*a%cv2) &
                       - a%patrial/a%lv2
!             
!
      return
      end subroutine setsystem_sys_open
!
! *** set up matrices for the simultaneous solve at t = t_{n+alfi}
!
!     a_lhs * x_{n+1} = a_rhs * x_{n+1} + b_rhs * x_{n} + c_rhs + d_rhs * q
!
!     a_lhs - x_{n+alfi} coefficients [xdmin,xdim]
!     a_rhs - x_{n+alfi} coefficients [xdim,xdim]
!     b_rhs - x_{n} coefficients [xdim,xdim]
!     c_rhs - constant coefficients [xdim]
!     d_rhs - flow_{n+alfi} coefficients [xdim,rcrflowsurfs]
!
      subroutine setsystem_sys(a, a_lhs, a_rhs, b_rhs, c_rhs, d_rhs, charvar)
      implicit none
      class(systemiccircuit) :: a
      real*8, intent(out) :: a_lhs(a%xdim,a%xdim)
      real*8, intent(out) :: a_rhs(a%xdim,a%xdim)
      real*8, intent(out) :: b_rhs(a%xdim,a%xdim)
      real*8, intent(out) :: c_rhs(a%xdim)
      real*8, intent(out) :: d_rhs(a%xdim,a%rcrsurfnum)            
      character(len=*) :: charvar
      character(len=*), parameter :: solve = 'solve'
      character(len=*), parameter :: update = 'update'
      real*8 :: rd_ctrl, cv_ctrl, vu_ctrl
      real*8 :: alfi_delt 
      integer :: i, j, nsrf      
!
!     ! timestep size
      if (charvar .eq. solve) then
         alfi_delt = real(alfi*delt,8)
      elseif (charvar .eq. update) then
         alfi_delt = delt
      end if
!
!     ! number of rcr surfaces
      nsrf = a%rcrsurfnum 
!
!     ! zero matrices
      a_lhs(:,:) = real(0.0,8)
      a_rhs(:,:) = real(0.0,8)
      b_rhs(:,:) = real(0.0,8)
      c_rhs(:) = real(0.0,8)
      d_rhs(:,:) = real(0.0,8)
!
!     ! control at t_{n+alfi} or t_{n+1} depending upon solve/update
      rd_ctrl = a%control%n1(a%rdindx)    
      cv_ctrl = a%control%n1(a%cpindx)    
      vu_ctrl = a%control%n1(a%vuindx)     
!
!     ! rcr surfaces
      do i = 1, nsrf 
!
!        ! pc index
         j = i + nsrf      
!
         ! p - pc = r q
         a_lhs(i,i) = real(1,8)
         a_lhs(i,j) = real(-1,8)
         d_rhs(i,i) = a%rp(i)
!        
         ! c dpc/dt = q - (pc - pa3)/rd 
         ! convert pa3 to va3 
         ! pa3 = (va3-vua3)/ca3        
         a_lhs(j,j) = a%c(i)/alfi_delt
         b_rhs(j,j) = a%c(i)/alfi_delt
!        
!        ! rd
         a_rhs(j,j) = real(-1,8)/(rd_ctrl*a%rd(i))
         a_rhs(j,a%va3indx) = real(1,8)/(rd_ctrl*a%rd(i)*a%ca3)
         c_rhs(j) = -a%vua3/(rd_ctrl*a%rd(i)*a%ca3)         
!                        
         d_rhs(j,i) = real(1,8)        
!        
        ! dva3/dt = q(1) + ... + q(nbranch)
         a_rhs(a%va3indx,j) = real(1,8)/(rd_ctrl*a%rd(i))
         a_rhs(a%va3indx,a%va3indx) = a_rhs(a%va3indx,a%va3indx) &
                                    - real(1,8)/(rd_ctrl*a%rd(i)*a%ca3)  
         c_rhs(a%va3indx) = c_rhs(a%va3indx) &
                          + a%vua3/(rd_ctrl*a%rd(i)*a%ca3)
!    
      end do
!
!     ! time derivatives - q_{i} equations are divided by l_{i} 
      do i = a%va3indx, a%qmvindx
         a_lhs(i,i) = real(1,8)/alfi_delt
         b_rhs(i,i) = real(1,8)/alfi_delt
      end do
!
!     ! va3
      a_rhs(a%va3indx,a%va3indx) = a_rhs(a%va3indx,a%va3indx) &
                                 - real(1,8)/(rd_ctrl*a%ra3*a%ca3)
      a_rhs(a%va3indx,a%vv1indx) = real(1,8)/(rd_ctrl*a%ra3*cv_ctrl*a%cv1)
      c_rhs(a%va3indx) = c_rhs(a%va3indx) &
                       + a%vua3/(rd_ctrl*a%ra3*a%ca3) &
                       - (vu_ctrl*a%vuv1)/(rd_ctrl*a%ra3*cv_ctrl*a%cv1)
!        
!     ! vv1
      a_rhs(a%vv1indx,a%va3indx) = real(1,8)/(rd_ctrl*a%ra3*a%ca3)
      a_rhs(a%vv1indx,a%vv1indx) = real(-1,8)/(rd_ctrl*a%ra3*cv_ctrl*a%cv1) &
                                 + real(-1,8)/(a%rv1*cv_ctrl*a%cv1)
      a_rhs(a%vv1indx,a%vv2indx) = real(1,8)/(a%rv1*cv_ctrl*a%cv2)
      c_rhs(a%vv1indx) = (vu_ctrl*a%vuv1)/(rd_ctrl*a%ra3*cv_ctrl*a%cv1)  &
                       + (vu_ctrl*a%vuv1)/(a%rv1*cv_ctrl*a%cv1)          &
                       - a%vua3/(rd_ctrl*a%ra3*a%ca3)                    &
                       - (vu_ctrl*a%vuv2)/(a%rv1*cv_ctrl*a%cv2) 
!        
!     ! vv2
      a_rhs(a%vv2indx,a%vv1indx) = real(1,8)/(a%rv1*cv_ctrl*a%cv1)
      a_rhs(a%vv2indx,a%vv2indx) = real(-1,8)/(a%rv1*cv_ctrl*a%cv2)
      a_rhs(a%vv2indx,a%qv2indx) = real(-1,8)
      c_rhs(a%vv2indx) = (vu_ctrl*a%vuv2)/(a%rv1*cv_ctrl*a%cv2)          &
                       - (vu_ctrl*a%vuv1)/(a%rv1*cv_ctrl*a%cv1)
!    
!     ! qv2
      a_rhs(a%qv2indx,a%vv2indx) = real(1,8)/(a%lv2*cv_ctrl*a%cv2)
      a_rhs(a%qv2indx,a%qv2indx) = real(-1,8)*a%rv2/a%lv2
      a_rhs(a%qv2indx,a%vlaindx) = real(-1,8)*(a%ela/a%lv2)
      c_rhs(a%qv2indx) = a%ela*(a%vula/a%lv2)                            &
                       - (vu_ctrl*a%vuv2)/(a%lv2*cv_ctrl*a%cv2)
!
!     ! vla
      a_rhs(a%vlaindx,a%qv2indx) = real(1,8)
      a_rhs(a%vlaindx,a%qmvindx) = real(-1,8)               
!
      ! vlv      
      a_rhs(a%vlvindx,a%qmvindx) = real(1,8)
!
      ! qmv
      a_rhs(a%qmvindx,a%vlaindx) = a%ela/a%lmv
      a_rhs(a%qmvindx,a%qmvindx) = real(-1,8)*(a%rmv/a%lmv)
      a_rhs(a%qmvindx,a%vlvindx) = real(-1,8)*(a%elv_n1/a%lmv)
      c_rhs(a%qmvindx) = a%elv_n1*(a%vulv/a%lmv) - a%ela*(a%vula/a%lmv)
!
      return
      end subroutine
!
! *** set implicit coefficients in the reduced order model
!
      subroutine setimplicitcoeff_sys(from_multdom,to_rom)
      implicit none
      class(systemiccircuit) :: from_multdom
      class(reducedorder) :: to_rom
      integer :: i, j
      do i = 1, from_multdom%surfnum
         do j = 1, to_rom%surfnum
            if (from_multdom%surfids(i) .eq. to_rom%surfids(j)) then
               to_rom%implicitcoeff(j,1:2) = from_multdom%implicitcoeff(i,1:2)
            end if
         end do
      end do
      end subroutine
!
!
!
!
! *** subroutine update xn to xn+1
!
      subroutine update_sys_open(a,stepn)
      implicit none
      class(systemiccircuit) :: a 
      integer :: stepn
      real*8 :: a_lhs(a%xdim,a%xdim)
      real*8 :: a_rhs(a%xdim,a%xdim)
      real*8 :: b_rhs(a%xdim,a%xdim)
      real*8 :: c_rhs(a%xdim)
      real*8 :: d_rhs(a%xdim,a%rcrsurfnum)      
      real*8 :: amat(a%xdim,a%xdim)
      real*8 :: bvec(a%xdim)
      real*8 :: qvec(a%xdim)
      real*8 :: c_temp(a%xdim)   
      real*8 :: tmpvec(a%xdim) 
      real*8 :: tmpmat(a%xdim,a%xdim)      
      real*8 :: tmpsol(a%xdim)       
      real*8 :: paorta
      real*8 :: dtime_n1
      integer :: i, j
!
!     ! update control to n+1
      call a%control%update(a%feedback)
!
!     ! set system to be solved - TECHNICALLY WE ARE AT N+1 NOT ALFI, NEED TO MODIFY THIS
      call a%setsystem_sys_open(a_lhs, a_rhs, b_rhs, c_rhs, d_rhs, 'update') 
!
!     ! dsegv matrix a
      amat = a_lhs - a_rhs
!
!     ! dsegv vector b => b_rhs * xn + c_rhs
      bvec = multiply_ab(b_rhs, a%x_n) 
      bvec = bvec + c_rhs
!
!     ! add latest flows => dsegv_b = dsegv_b + d_rhs * flows      
      qvec = multiply_ab(d_rhs, a%flowpntr(1:a%rcrsurfnum))
      bvec = bvec + qvec
!
      a%x_n1 = solve_axb(amat, bvec)
!      
!     ! store current x vars
      a%x_hist(stepn,:) = a%x_n1(:) 
!
!     ! store all current flow for next time step
      do i = 1, a%rcrsurfnum         
         a%flow_n(i) = a%flowpntr(i)%p 
      end do        
!
!     ! restart flow history
      a%r_hist(stepn,:,1) = a%flow_n(:)
!
!     ! restart pressure history - rcr surfaces
      do i = 1, a%rcrsurfnum
         a%r_hist(stepn,i,2) = a%implicitcoeff(i,1)*a%flow_n(i) &
                             + a%implicitcoeff(i,2)
      end do
!
!
!     ! store current flow, pressure, outlet and capacitance pressures      
!
!     ! q 
      a%q_hist(stepn,1:a%surfnum) = a%flow_n(1:a%surfnum)
!     ! p/pc      
      a%p_hist(stepn,1:2*a%rcrsurfnum) = a%x_n1(1:2*a%rcrsurfnum)
!     ! pa3
      a%p_hist(stepn,a%va3indx) = (a%x_n1(a%va3indx)-a%vua3)/a%ca3
!     ! pv1
      a%p_hist(stepn,a%vv1indx) = (a%x_n1(a%vv1indx)-a%vuv1)/a%cv1
!     ! pv2
      a%p_hist(stepn,a%vv2indx) = (a%x_n1(a%vv2indx)-a%vuv2)/a%cv2
!     ! pla - minus 1 for qv2
!
!     ! control
      do i = 1, a%control%dim
         a%control%hist(stepn,i) = a%control%n1(i)
      end do
      
!      
!     ! feedback
      if (a%feedbackactive) then
!        ! update feedback - update latest value of average pressure
         call a%feedback%update(stepn)
!
      end if 
!
!     ! update
      a%x_n(:) = a%x_n1(:) 
!
!
!
      return
      end subroutine update_sys_open      
!
! *** subroutine update xn to xn+1
!
      subroutine update_sys(a,stepn)
      implicit none
      class(systemiccircuit) :: a 
      integer :: stepn
      real*8 :: a_lhs(a%xdim,a%xdim)
      real*8 :: a_rhs(a%xdim,a%xdim)
      real*8 :: b_rhs(a%xdim,a%xdim)
      real*8 :: c_rhs(a%xdim)
      real*8 :: d_rhs(a%xdim,a%rcrsurfnum)      
      real*8 :: amat(a%xdim,a%xdim)
      real*8 :: bvec(a%xdim)
      real*8 :: qvec(a%xdim)
      real*8 :: c_temp(a%xdim)   
      real*8 :: tmpvec(a%xdim) 
      real*8 :: tmpmat(a%xdim,a%xdim)      
      real*8 :: tmpsol(a%xdim)       
      real*8 :: paorta
      real*8 :: dtime_n1
      integer :: i
      real*8 :: xiter(2), j(2,2), r(2), dx(2)    
      real*8 :: dtime, eval  
      real*8, parameter :: tol = real(1e-6,8)
      integer, parameter :: itermax = 10      
!
!     ! update control to n+1
      call a%control%update(a%feedback)
! 
!     ! update activation time, needs latest feedback
      call a%update_activationtime() 
      a%act_hist(stepn) = a%activationtime
!
!     ! delta time to t = t_{n+1}
!     ! this is zero as we are already at at n+1 and activation time is updated above
      dtime_n1 = real(0.0,8)
!
!     ! elastance at t = t_{n+1}
      a%elv_n1 = a%getelastance(dtime_n1,a%control%n1)

!     ! set system to be solved - TECHNICALLY WE ARE AT N+1 NOT ALFI, NEED TO MODIFY THIS
      call a%setsystem_sys(a_lhs, a_rhs, b_rhs, c_rhs, d_rhs, 'update') 
!
!     ! dsegv matrix a
      amat = a_lhs - a_rhs
!
!     ! dsegv vector b => b_rhs * xn + c_rhs
      bvec = multiply_ab(b_rhs, a%x_n) 
      bvec = bvec + c_rhs
!
!     ! add latest flows => dsegv_b = dsegv_b + d_rhs * flows      
      qvec = multiply_ab(d_rhs, a%flowpntr(1:a%rcrsurfnum))
      bvec = bvec + qvec
!
!     ! solve for state variables at n+1
      if (a%mvopen) then
         a%x_n1 = solve_axb(amat, bvec)
      else if (a%avopen) then
         tmpsol(1:a%vlaindx) = solve_axb(amat(1:a%vlaindx,1:a%vlaindx), &
                                         bvec(1:a%vlaindx))
         a%x_n1(1:a%vlaindx) = tmpsol(1:a%vlaindx)
         a%x_n1(a%qmvindx) = real(0.0,8)
      else
         tmpsol(1:a%vlvindx) = solve_axb(amat(1:a%vlvindx,1:a%vlvindx), &
                                         bvec(1:a%vlvindx))
         a%x_n1(1:a%vlvindx) = tmpsol(1:a%vlvindx)
         a%x_n1(a%qmvindx) = real(0.0,8)
      end if
!
!     ! if valve open
      if (a%avopen) then
!              
!        ! initial pressure and volume at previous step
         xiter(1) = a%paorta_n
         xiter(2) = a%x_n(a%vlvindx)
!
!        ! solve timestep        
         dtime = delt        

!        ! jacobian
         j(:,:) = real(0.0,8)
         j(1,1) = real(1.0,8)
         j(1,2) = real(-1.0,8)*a%elv_n1*(real(1.0,8) + a%flowpntr(a%surfnum)%p) ! flow -ve
         j(2,2) = real(1.0,8)/dtime
!!         j(1,1) = real(1.0,8)
!!         j(2,2) = real(-1.0,8)*a%elv_n1
!!         j(2,2) = real(1.0,8)/dtime         
!
!        ! initial value and evaluation
         eval = real(1e20,8)
!
         do i = 1,itermax
!
            r(1) = xiter(1) - a%stabpntr(a%surfnum)%p                       &
                 - a%elv_n1*(xiter(2) - a%vulv)*(real(1.0,8) +              &
                                            a%kelv*a%flowpntr(a%surfnum)%p) &! flow -ve
                 - a%rav*a%flowpntr(a%surfnum)%p                            &
                 + (a%lav/dtime)*(a%qaorta_n - a%flowpntr(a%surfnum)%p)
            r(2) = xiter(2)/dtime - a%x_n(a%vlvindx)/dtime                  &
                 - a%flowpntr(a%surfnum)%p
!
!!            r(1) = xiter(1) - a%elv_n1*(xiter(2) - a%vulv) 
!!     &           - a%rav*a%flowpntr(a%surfnum)%p 
!!     &           - (a%lav/dtime)*(a%flowpntr(a%surfnum)%p -a%qaorta_n)
!!            r(2) = xiter(2)/dtime - a%x_n(a%vlvindx)/dtime 
!!     &           - a%flowpntr(a%surfnum)%p
!
            dx = solve_axb(j,real(-1.0,8)*r)  
            xiter = xiter + dx   
!
            eval = dot_product(dx,dx)
            eval = sqrt(eval) 
!
!           ! if evaluation is less than tolerance exit
            if (eval .lt. tol) then
               exit
            end if
!         
         end do
!    
!        ! implicit coefficients
         a%implicitcoeff(a%surfnum,1) = a%rav + a%lav/dtime                   &
                                      + a%elv_n1*(xiter(2) - a%vulv)*a%kelv
         a%implicitcoeff(a%surfnum,2) = a%elv_n1*(xiter(2) - a%vulv)          &
                                      - (a%lav/dtime)*a%qaorta_n
!!c
!
!!         a%implicitcoeff(a%surfnum,1) = a%rav + a%lav/dtime
!!         a%implicitcoeff(a%surfnum,2) = a%elv_n1*(xiter(2) - a%vulv)
!!     &                                - (a%lav/dtime)*a%qaorta_n
!
         ! vlvcoeff         
         a%vlvcoeff(1) = dtime
         a%vlvcoeff(2) = a%x_n(a%vlvindx)
!
!        ! paorta
         a%paorta_n = a%implicitcoeff(a%surfnum,1)*a%flowpntr(a%surfnum)%p    &
                    + a%implicitcoeff(a%surfnum,2)                            &
                    - a%stabpntr(a%surfnum)%p    

         a%x_n1(a%vlvindx) = a%vlvcoeff(1)*a%flowpntr(a%surfnum)%p            &
                           + a%vlvcoeff(2)   
         
         a%plv_n = a%elv_n1*(a%x_n1(a%vlvindx) - a%vulv)
         a%plv_n = a%plv_n*(real(1.0,8) + a%kelv*a%flowpntr(a%surfnum)%p) !! flow negative    
         a%plv_n = a%plv_n - a%stabpntr(a%surfnum)%p    
!
!!         a%paorta_n = xiter(1)
!!         a%x_n1(a%vlvindx) = xiter(2)
!     
      else
         ! ?????
         a%plv_n = a%elv_n1*(a%x_n1(a%vlvindx) - a%vulv)
      end if   
!
!     ! store elv history
      a%elv_hist(stepn) = a%elv_n1
!      
!     ! store current x vars
      a%x_hist(stepn,:) = a%x_n1(:) 
!
!     ! store all current flow for next time step
      do i = 1, a%rcrsurfnum         
         a%flow_n(i) = a%flowpntr(i)%p 
      end do        
!
!     ! zero flow if av closed
      if (a%avopen) then
         a%flow_n(a%surfnum) = a%flowpntr(a%surfnum)%p 
      else
         a%flow_n(a%surfnum) = real(0.0,8)
      end if
!
!     !! set left atrium pressure for next time step !!
      a%pla_n = a%ela*(a%x_n1(a%vlaindx)-a%vula)
!
!     ! restart flow history
      a%r_hist(stepn,:,1) = a%flow_n(:)
!
!     ! restart pressure history - rcr surfaces
      do i = 1, a%rcrsurfnum
         a%r_hist(stepn,i,2) = a%implicitcoeff(i,1)*a%flow_n(i)      &
                             + a%implicitcoeff(i,2)
      end do
!
!     ! restart pressure history - heart surface
      if (a%avopen) then
         a%r_hist(stepn,a%surfnum,2) = a%paorta_n
      else
         a%r_hist(stepn,a%surfnum,2) = a%presspntr(a%surfnum)%p
      end if
!
!     ! store current flow, pressure, outlet and capacitance pressures      
!
!     ! q 
      a%q_hist(stepn,1:a%surfnum) = a%flow_n(1:a%surfnum)
!     ! p/pc      
      a%p_hist(stepn,1:2*a%rcrsurfnum) = a%x_n1(1:2*a%rcrsurfnum)
!     ! pa3
      a%p_hist(stepn,a%va3indx) = (a%x_n1(a%va3indx)-a%vua3)/a%ca3
!     ! pv1
      a%p_hist(stepn,a%vv1indx) = (a%x_n1(a%vv1indx)-a%vuv1)/a%cv1
!     ! pv2
      a%p_hist(stepn,a%vv2indx) = (a%x_n1(a%vv2indx)-a%vuv2)/a%cv2
!     ! pla - minus 1 for qv2
      a%p_hist(stepn,a%vlaindx-1) = a%ela*(a%x_n1(a%vlaindx)-a%vula)

!     ! set av valve state based on pressure at t_{n+1}
!     ! stored in plv_n but pressures at t_{n+1}
!
!!      a%plv_n = a%elv_n1*(a%x_n1(a%vlvindx)-a%vulv)
!!      a%plv_n = a%plv_n*(real(1.0,8) + a%kelv*a%flowpntr(a%surfnum)%p)      
!      
!     ! plv
      a%p_hist(stepn,a%vlvindx-1) = a%plv_n
!
      do i = 1, a%surfnum
         a%stab_hist(stepn,i) = a%stabpntr(i)%p
      end do
!
!     ! control
      do i = 1, a%control%dim
         a%control%hist(stepn,i) = a%control%n1(i)
      end do
      
!      
!     ! feedback
      if (a%feedbackactive) then
!        ! update feedback - update latest value of average pressure
         call a%feedback%update(stepn)
!
      end if 
!
!     ! update
      a%x_n(:) = a%x_n1(:) 
!
!
!     ! set aortic pressure from 3D or 0D domain (stored in the last variable)
      if (a%updatepressure) then
         a%paorta_n = a%presspntr(a%surfnum)%p                                         
      end if
!      
      if (a%plv_n .gt. a%paorta_n) then                            
         a%avopen = .true.                                     
         hrt%avopen = int(1) 
         a%updatepressure = int(0)
         a%qaorta_n = a%flowpntr(a%surfnum)%p
      elseif (a%flow_n(a%surfnum) .lt. real(0.0,8)) then   
         a%avopen = .true.                                        
         hrt%avopen = int(1) 
         a%updatepressure = int(0)                            
         a%qaorta_n = a%flowpntr(a%surfnum)%p
      elseif (a%flow_n(a%surfnum) .gt. real(0.0,8)) then ! reverse flow
         a%avopen = .false.                                     
         hrt%avopen = int(0) 
         a%updatepressure = int(1)                            
         a%qaorta_n = real(0.0,8)
      else                                                     
         a%avopen = .false.                                    
         hrt%avopen = int(0) 
         a%updatepressure = int(1)
         a%qaorta_n = real(0.0,8)
      end if  
!
!
      return
      end subroutine update_sys
!
! *** update activation time, increase 
!
      subroutine update_activationtime(a)
      implicit none
      class(systemiccircuit) :: a
      real*8 :: heartrate
      real*8 :: period
!
!     ! current heart rate in control(1)
      if (a%feedback%active) then
         heartrate = a%control%n1(a%hrindx)/a%period         
      else
         heartrate = real(1.0,8)/a%period
      end if
      period = real(1.0,8)/heartrate
!
!     ! update activation time by delt/period
      a%activationtime = a%activationtime + delt/period
!
!     ! if greater than 1, reset      
      if (a%activationtime .gt. real(1.0,8)) then
         a%activationtime = a%activationtime - real(1.0,8)
      end if 
!      
      end subroutine
!
!
!
      subroutine update_activationtime_hrt(this,stepn)
         use externalDataTools
         implicit none
         class(numericalheart) :: this
         integer :: stepn

         ! update activation time by delt/period
         this%activationtime = this%activationtime + delt/this%period !\todo where is delt coming from?
         this%act_hist(stepn) = this%activationtime
         ! write(*,*) 'in update_activationtime_hrt', this%activationtime !\todo remove

         ! if greater than 1, reset
         if (this%newBeatJustStarted .eq. int(1)) then
            if (this%inputHRandSP .eq. int(1)) then
               call writeRestartForLinkedList('SystolicPressure.dat')
               call writeRestartForLinkedList('HeartRate.dat')
            end if
            this%newBeatJustStarted = int(0)
         end if
         if (this%activationtime .gt. real(1.0,8)) then
            this%activationtime = this%activationtime - real(1.0,8)
            this%newBeatJustStarted = int(1)

            this%totalTimeOfAllCompletedPeriods = this%totalTimeOfAllCompletedPeriods + this%period

            write(*,*) 'updating heart activation times...' !\todo remove

            ! Things to do at the start of a new beat:
            if (this%inputHRandSP .eq. int(1)) then
               ! Prepare to update emax by removing the systolicPressure value from the previous beat:
               this%emax = this%emax / this%systolicPressureList%value ! (note that we dont need to take out the factor 0.9, as we'd just put it back below..)

               this%heartRateList => this%heartRateList%next
               this%systolicPressureList => this%systolicPressureList%next

               if ((this%heartRateList%value .eq. -1.0) .or. (this%systolicPressureList%value .eq. -1.0)) then
                  write(*,*) 'The imposed systolic pressure or heart rate data has run out (Input data too short, or simulation time too long).'
                  stop
               end if

               ! Complete the update of emax, begun above:
               this%emax = this%emax * this%systolicPressureList%value

               this%period = 60/this%heartRateList%value
               this%tmax = 0.3*sqrt(this%period)
               this%trelax = 0.5*this%tmax
            end if
         end if

      end subroutine update_activationtime_hrt



!
! *** update feedback with current values
!
      subroutine update_pressavg(a,stepn)
      implicit none
      class(feedback) :: a
      integer :: i, j, stepn, tstep, pstep, tmaxloc
      real*8 :: tperiod, tmax
      real*8 :: taverage(a%surfnum), tdifference(a%surfnum)
      real*8 :: outputmax
!
      real*8 :: sum_x(a%surfnum)
      real*8 :: sum_y(a%surfnum)            
      real*8 :: sum_x_y(a%surfnum)
      real*8 :: sum_x_2(a%surfnum)
      real*8 :: m(a%surfnum)
      real*8 :: c(a%surfnum)
      integer :: n
      real*8 :: amat(2,2), bvec(2), xvec(2)
!
!      
!     ! update feedback input (current pressure)
      do i = 1, a%surfnum   
         a%input(stepn,i) = a%pntr(i)%p                        
      end do
!
!     ! calculate feedback output
      call a%calculate(stepn) 
!
!     ! if feedback not activated
      if (a%active .eq. .false.) then

         tperiod = a%params(1)
         !!tstep = floor(tperiod/delt)
         ! go back 10% of the cycle
         tstep = floor(real(0.1,8)*tperiod/delt)
!
         if (stepn .gt. tstep) then

!           ! compare current cycle average to 
!           ! same average one cycle prior (pstep)
            pstep = stepn - tstep + int(1);
!         
! 
            ! loop over surfaces
            do i = 1, a%surfnum
               sum_x(i) = real(0.0,8)
               sum_y(i) = real(0.0,8)
               sum_x_y(i) = real(0.0,8)
               sum_x_2(i) = real(0.0,8)
               n = 0
               do j = pstep, stepn
                  sum_x(i) = sum_x(i) + real(j,8)*delt      
                  sum_y(i) = sum_y(i) + a%input(j,i)      
                  sum_x_y(i) = sum_x_y(i) + real(j,8)*delt*a%input(j,i)      
                  sum_x_2(i) = sum_x_2(i) + real(j,8)*real(j)*delt*delt
                  n = n + 1
               end do
               amat(1,1) = sum_x_2(i)
               amat(1,2) = sum_x(i)
               amat(2,1) = sum_x(i)
               amat(2,2) = real(n,8)
               bvec(1) = sum_x_y(i)
               bvec(2) = sum_y(i)
               xvec = solve_axb(amat,bvec)
               m(i) = xvec(1)
               c(i) = xvec(2)
               tdifference(i) = m(i)*(real(n,8)*delt)/a%input(stepn,i)               
            end do 
!            
!!c           ! difference between current average 
!!c           ! and previous cycle average
!!            tdifference(:) = a%outputval(:) - a%output(pstep,:) 
!!            tdifference(:) = abs(tdifference(:))/a%outputval(:)
!
            tmax = maxval(abs(tdifference))
            a%tolhist(stepn,:) = tdifference(:)
            tmaxloc = maxloc(tdifference,1)
!
!           ! if maximum average difference is small, then activate the feedback
            if (tmax .lt. a%activation_tol) then
               a%active = .true.
!              !!! target value is set as the largest difference at activation 
               !!a%targetval = a%outputval(tmaxloc)                 
!
               ! target value for each branch recorded               
               a%targetval(:) = a%outputval(:)            
            end if
!
         end if
!
      end if
!
      return
      end subroutine
!
! *** calculate average pressure
!
      subroutine calculate_pressavg(a,stepn)
      implicit none
      class(feedback) :: a
      integer :: stepn
      real*8 :: time_n1, period, heartrate
      integer :: i, startindx, endindx, counter, nsteps
!
!     ! period stored in param array 
      period = a%params(1)
!
!     ! if feedback active, then use heart rate 
      if (a%active) then
         heartrate = a%hr_ctrl%p/period
         period = real(1.0,8)/heartrate
      end if
!
!     ! time at current time t = t_{n+1}
      time_n1 = real(stepn,8)*delt
!
      if (time_n1 .lt. period) then
         nsteps = stepn 
      else
         nsteps = floor(period/delt) 
      end if                  
!
      a%outputval(1:a%surfnum) = real(0.0,8)
      startindx = stepn
      endindx = startindx - nsteps + int(1)
      counter = int(0)
!         
!     ! loop backwards
      do i = startindx, endindx, -1 
         a%outputval(1:a%surfnum) = a%outputval(1:a%surfnum)               &
                                  + a%input(i,1:a%surfnum)
         counter = counter + int(1)
      end do
!         
!     ! current value
      a%outputval(1:a%surfnum) = a%outputval(1:a%surfnum)/real(counter,8)
!
!     ! store current value 
      a%output(stepn,1:a%surfnum) = a%outputval(1:a%surfnum)
!         
      return
      end subroutine
!
! *** write x variable history to file
!
      subroutine writefile_sys(a,stepn)
      implicit none
!      
      include "mpif.h"                      
      class(systemiccircuit) :: a
      integer :: stepn
      integer :: rank
      integer :: ierr 
      integer :: xnum = 478
      integer :: qnum = 874
      integer :: pnum = 561
      integer :: cnum = 156
      integer :: enum = 678
      integer :: rnum = 101
      integer :: fnum = 359
      integer :: anum = 376
      integer :: tnum = 888
      integer :: i
!
      if (mod(stepn,ntout) .eq. int(0)) then
         call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)              
         if (rank .eq. int(0)) then
!
!           ! write x variable history
            open(xnum, file=a%xvarsfilename, status='replace', iostat=ierr)
            do i = 1, stepn
               write(xnum,a%xvarsformat) i, a%x_hist(i,:)            
            end do
            close(xnum) 
!
!           ! write flow history
            open(qnum, file=a%qhistfilename, status='replace', iostat=ierr)
            do i = 1, stepn
               write(qnum,a%qvarsformat) i, a%q_hist(i,:)            
            end do
            close(qnum)
!
!           ! write pressure history
            open(pnum, file=a%phistfilename, status='replace', iostat=ierr)
            do i = 1, stepn
               write(pnum,a%pvarsformat) i, a%p_hist(i,:)            
            end do
            close(pnum) 
!
!           ! write control history
            open(cnum, file=a%chistfilename, status='replace', iostat=ierr)
            do i = 1, stepn
               write(cnum,a%cvarsformat) i, a%control%hist(i,:)            
            end do
            close(cnum) 
!               
!           ! write feedback history
            if (a%feedbackactive) then                           
               open(fnum, file=a%fhistfilename, status='replace', iostat=ierr)
               do i = 1, stepn
                  write(fnum,a%fvarsformat) i,                                      &
                                      a%feedback%input(i,1:a%feedback%surfnum),     &
                                      a%feedback%output(i,1:a%feedback%surfnum),    &
                                      a%feedback%targetval(1:a%feedback%surfnum)  
               end do
               close(fnum) 
!  
!              ! write tolerance history
               open(tnum, file='sysThist.dat', status='replace', iostat=ierr)
               do i = 1, stepn
                  write(tnum,a%tvarsformat) i, a%feedback%tolhist(i,:)     
               end do
               close(tnum)   


            end if
!
!           ! write elastance history
            open(enum, file='sysEhist.dat', status='replace', iostat=ierr)
            do i = 1, stepn
               write(enum,*) i, a%elv_hist(i)
            end do
            close(enum) 
!
!           ! write activation history
            open(anum, file=a%ahistfilename, status='replace')
            do i = 1, stepn
               write(anum,a%avarsformat) i, a%act_hist(i)
            end do
            close(anum) 
!            
!           ! write flow restart history
            open(rnum, file=a%flowfile, status='replace', iostat=ierr)
            do i = 1, stepn
               write(rnum,a%rvarsformat) a%r_hist(i,:,1)
            end do
            close(rnum)             
!            
!           ! write pressure restart history
            open(rnum, file=a%pressurefile, status='replace', iostat=ierr)
            do i = 1, stepn
               write(rnum,a%rvarsformat) a%r_hist(i,:,2)
            end do
            close(rnum)             

!           ! write stabilsation history
            open(129, file='PStab.dat', status='replace', iostat=ierr)
            do i = 1, stepn
               write(129,a%pvarsformat) i, a%stab_hist(i,:)            
            end do
            close(129)             
!
         end if
      end if
!     
      end subroutine writefile_sys


      subroutine write_activation_history_hrt(this,stepn)
         implicit none

         include "mpif.h" 

         class(numericalheart) :: this
         integer :: anum = 376
         integer :: ii
         integer :: stepn
         integer :: rank
         integer :: ierr

         ! write activation history
         call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)              
         if (rank .eq. int(0)) then
            open(anum, file=this%ahistfilename, status='replace')
            do ii = 1, stepn
               write(anum,this%avarsformat) ii, this%act_hist(ii)
            end do
            close(anum) 
         end if

      end subroutine write_activation_history_hrt


      subroutine read_activation_history_hrt(this,stepn)
         implicit none
         class(numericalheart) :: this
         integer :: anum = 472
         integer :: ii
         integer :: jj
         integer :: stepn
         integer :: ierr

         ! read activation history
         open(anum, file=this%ahistfilename, status='old',iostat=ierr)      
         if (ierr .eq. int(0)) then   
            do ii = 1, stepn
               read(anum,*) jj, this%act_hist(ii)
            end do
            this%activationtime = this%act_hist(stepn)
         else 
            write(*,*) ' WARNING! No ***Ahist.dat found! Continuing ...'
         end if
         close(anum)

!
      end subroutine read_activation_history_hrt
!
! *** load x variable history from file
!
      subroutine loadxvars(a,torow)
      !\todo work out why this subroutine not called anywhere!
      implicit none
!
      class(systemiccircuit) :: a     
      integer :: torow

      integer :: ierr 
      integer :: fnum = 357
      integer :: i
      character(len=100) :: char1
      character(len=100) :: char2           
      integer :: stepn
      real*8 :: xvars(a%xdim)      
! 
      open(fnum, file=a%xvarsfilename, status='old', iostat=ierr)
      do i = 1, torow
         read(fnum,*) stepn, xvars
      end do      
      close(fnum)
!
!     ! set x_{n} from last values read
      a%x_n(:) = xvars
!
      return
      end subroutine loadxvars      
!
! *********************************************
! *** heart model functions and subroutines ***
! *********************************************
!
! *** constructor function
!     heart parameters parsed in input_fform.cxx and stored in common block
!
      function hrtconstructor(heartparameters,inputHRandSP,lstep_passedIn)
      
      implicit none
      type(numericalheart) :: hrtconstructor
      real*8 :: heartparameters(0:15)
      integer, intent(in) :: inputHRandSP
      integer, intent(in) :: lstep_passedIn

      call hrtconstructor%initialise_hrt(heartparameters,inputHRandSP,lstep_passedIn)
      
      return
      end function
!
! *** initialisation subroutine
!
      subroutine initialise_hrt(this, parameters,inputHRandSP,lstep_passedIn)

      use externalDataTools

      implicit none

      include "mpif.h"               

      class(numericalheart), intent(inout) :: this
      real*8 :: parameters(0:15)  
      integer :: hstep
      integer :: i
      integer :: hnum = 743
      integer :: elvnum = 744      
      integer :: anum = 377
      integer :: rank, ierr, nelv
      integer, intent(in) :: inputHRandSP
      integer, intent(in) :: lstep_passedIn
      real*8 :: tn(34)
      real*8 :: en(34)
      real*8  :: tval, eval      
      real*8 :: approxInitialStrokeVolume

      call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)

      ! active model
      this%isactive = int(1)

      ! Copy the common block flag, making it a member variable:
      this%inputHRandSP = inputHRandSP

      ! A label so we can identify this class (avoids a world of pain when working with a tower of derived types)
      this%classNameString = 'numericalheart'

      ! This has been refactored here from the systemic circuit      
      this%ahistfilename = 'sysAhist.dat'
      
      ! allocate history arrays
      if (lstep_passedIn .gt. int(0)) then
         hstep = nstep + lstep_passedIn
      else
         hstep = nstep
      end if

      ! \todo restart this!
      this%totalTimeOfAllCompletedPeriods = 0.0d0

      allocate(this%act_hist(hstep+1))                    ! activation history
      if (lstep_passedIn .eq. int(0)) then
         this%activationtime = real(0.0,8)
      else
         call this%read_activation_history_hrt(lstep_passedIn)
      end if

      ! activation
      write(this%avarsformat,'(3(a))') '(i8,e20.10)' 

      ! number of surfaces, set to 1
      this%surfnum = int(1)
      allocate(this%surfids(this%surfnum))
      
      ! assign parameters form array set in input_fform.cxx
      this%surfids(1:this%surfnum) = int(parameters(0))  ! aortic surface id
      this%patrial = parameters(1)                       ! atrial pressure       
      this%rmv = parameters(2)                           ! mitral valve    
      this%edv = parameters(3)                           ! end diastolic volume
      this%vulv = parameters(4)                          ! unstressed volume  
      this%emax = parameters(5)                          ! maximum elastance  
      this%tmax = parameters(6)                          ! time to maximum elastance
      this%trelax = parameters(7)                        ! time to relax - i.e. min. elastance after max. elastance
      this%period = parameters(8)                        ! period 
      this%kelv = parameters(9)                          ! ventricular resistance - kelv      
      this%rav = parameters(10)                          ! aortic valve
      this%ibackflow = int(parameters(11))               ! back flow flag
      this%m_backflow = parameters(12)                   ! back flow magnitude
      this%s_backflow = parameters(13)                   ! back flow steepness
      this%c_backflow = parameters(14)                   ! back flow closure
      this%max_backflow = parameters(15)                 ! back flow time

      ! set model variables from parameters
      this%vlv_n = this%edv            
      this%emin = this%patrial/(this%edv - this%vulv)     

      ! store tmax and trelax in non-dimensional units
      this%tmax = this%tmax/this%period            
      this%trelax = this%trelax/this%period

      ! reset values if inputHRandSP = 1
      ! period, time to tmax and trelax input in seconds, and maximum elastance
      if (this%inputHRandSP .eq. int(1)) then

         this%systolicPressureList => loadFileAsLinkedList('SystolicPressure.dat',lstep_passedIn)
         this%heartRateList => loadFileAsLinkedList('HeartRate.dat',lstep_passedIn)

         if (rank .eq. int(0)) then
            write(*,*) 'Over-writing heart period and elastance using data from HeartRate.dat and SystolicPressure.dat.'
         end if
         
         this%period = 60/this%heartRateList%value

         this%tmax = 0.3*sqrt(this%period)
         this%trelax = 0.5*this%tmax

         approxInitialStrokeVolume = 83.3
         this%emax = 0.9*this%systolicPressureList%value*(this%edv - approxInitialStrokeVolume)      
      ! else                
         ! this%period = parameters(1)

         ! this%tmax = parameters(11)
         ! this%tmax = this%tmax/this%period            
         ! this%trelax = parameters(12)
         ! this%trelax = this%trelax/this%period

         ! this%emax = parameters(2)
      end if
  
      ! write out values to file
      call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)              
      if (rank .eq. int(0)) then
         open(hnum, file='heartmodel.dat', status='replace')      
         write(hnum,*) '# heart model input parameters'
         write(hnum,*) '# '
         write(hnum,*) '# aortic surface'
         write(hnum,*) this%surfids(1:this%surfnum) 
         write(hnum,*) '# heart period [time]'
         write(hnum,*) this%period
         write(hnum,*) '# end diastolic volume [volume]'            
         write(hnum,*) this%edv
         write(hnum,*) '# unstressed volume [volume]'
         write(hnum,*) this%vulv
         write(hnum,*) '# preload [pressure]'
         write(hnum,*) this%patrial 
         write(hnum,*) '# tmax/tperiod [-]'
         write(hnum,*) this%tmax
         write(hnum,*) '# trelax/tperiod [-]'
         write(hnum,*) this%trelax
         write(hnum,*) '# emax [pressure/volume]'
         write(hnum,*) this%emax
         write(hnum,*) '# emin [pressure/volume]'
         write(hnum,*) this%emin
         write(hnum,*) '# rav [pressure/volume]'
         write(hnum,*) this%rav 
         write(hnum,*) '# kelv [1/time]'
         write(hnum,*) this%kelv         
         write(hnum,*) '# rmv [pressure/flow]'
         write(hnum,*) this%rmv
         close(hnum)
      end if

!
!     ! set logical to update inflow pressure value from face
      this%updatepressure = int(0)          
!
!     ! allocate and zero arrays 
      allocate(this%surfarea(this%surfnum))
      allocate(this%flow_n(this%surfnum))
      allocate(this%flow_n1(this%surfnum))
      allocate(this%pressure_n(this%surfnum))
      allocate(this%implicitcoeff(this%surfnum,2))
      allocate(this%implicitcoeff_n1(this%surfnum,2))      

      this%surfarea(:) = real(0.0,8)
      this%flow_n(:) = real(0.0,8)
      this%flow_n1(:) = real(0.0,8)
      this%pressure_n(:) = real(0.0,8)
      this%implicitcoeff(:,:) = real(0.0,8)            
      this%vlv_coeff(:) = real(0.0,8)      
!
!     ! allocate and zero history arrays      
      if (lstep_passedIn .gt. int(0)) then
         hstep = nstep + lstep_passedIn
      else
         hstep = nstep
      end if
!     
!     ! allocate aortic flow and pressure history arrays
      allocate(this%flowhist(hstep+1,1))       
      allocate(this%pressurehist(hstep+1,1))             
!      
!     ! write flow and pressure history file names
      write (this%flowfile,'(a)') 'Qaorta.dat'
      write (this%pressurefile,'(a)') 'Paorta.dat'
!            
!     ! allocate internal history arrays
      allocate(this%vlv_hist(hstep+1))      
      allocate(this%plv_hist(hstep+1))      
      allocate(this%qmv_hist(hstep+1))      
      allocate(this%elv_hist(hstep+1))         
      allocate(this%eval_hist(hstep+1,2))         
      allocate(this%stab_hist(hstep+1))         
!      
!     ! zero arrays
      this%flowhist(:,:) = real(0.0,8)
      this%pressurehist(:,:) = real(0.0,8)
      this%vlv_hist(:) = real(0.0,8)
      this%plv_hist(:) = real(0.0,8)
      this%qmv_hist(:) = real(0.0,8)
      this%elv_hist(:) = real(0.0,8)
      this%eval_hist(:,:) = real(0.0,8)
      this%stab_hist(:) = real(0.0,8)
!
!     ! senzaki elastance, zeroed and modified
!     !
      tn = (/  0.000000000, 0.066587859, 0.108926098, 0.159141738, 0.204260893,  &
               0.253208241, 0.342776414, 0.430810673, 0.565132010, 0.672851660,  &
               0.741853430, 0.821271212, 0.879264898, 0.925757782, 1.000000000,  &
               1.052920518, 1.078931147, 1.111102775, 1.136603451, 1.167835266,  &
               1.183370068, 1.218328693, 1.249934608, 1.278512236, 1.297685646,  &
               1.324613784, 1.366134883, 1.400884661, 1.420531527, 1.448895238,  &
               1.481469353, 1.516814245, 1.553492314, 1.601764457 /)
      en = (/  0.000000000, 0.049089405, 0.115038061, 0.203182740, 0.291300279,  &
               0.363929671, 0.441941430, 0.514037796, 0.632163282, 0.729471447,  &
               0.814021729, 0.890502803, 0.943240884, 0.981885804, 1.000000000,  &
               0.990938280, 0.972407549, 0.916584276, 0.856419017, 0.774032075,  &
               0.705920401, 0.552490548, 0.439957667, 0.357044200, 0.297462444,  &
               0.221308336, 0.143998827, 0.099877285, 0.074752855, 0.050286845,  &
               0.031380477, 0.017409095, 0.007752265, 0.000000000 /)
!
!     
      allocate(this%elv_senzaki%v(34,2))       
!
      do i = 1,34
         this%elv_senzaki%v(i,1) = tn(i)
         this%elv_senzaki%v(i,2) = en(i)
      end do

      ! load elastance from input file ELV.in 
      open(elvnum, file='ELV.in', status='old', iostat=ierr)
      
      ! if file exists
      if (ierr .eq. int(0)) then         

         ! set boolean true
         this%input_elastance = int(1)
         
         ! read number of data points and read allocate time data object
         read(elvnum,*) nelv
         allocate(this%elv_input%v(nelv,2))

         ! read data from file
         do i = 1,nelv         
            read(elvnum,*) tval,eval
            this%elv_input%v(i,1) = tval
            this%elv_input%v(i,2) = eval
         end do 

      else

         ! set boolean false      
         this%input_elastance = int(0)

      end if

      ! close file      
      close(elvnum)

      ! 
      if (lstep_passedIn .eq. int(0)) then
         ! Note that currently these are only loaded-on-restart if the controlledcoronarymodel is active (from coronary_restart.dat).
         ! If you want to use it somewhere else in the code, make sure it's loaded appropriately on restarts.
         this%timestepsSinceLastHeartPeriodReset = int(0)
         this%newBeatJustStarted = int(0)
      end if

      return
      end subroutine initialise_hrt

      ! assign pointers for filter

      subroutine assign_ptrs_ext_hrt(a)

      use iso_c_binding
      use phcommonvars, only: PhAssignPointerDP

      implicit none

      class(numericalheart) :: a

      ! set pointer to EMax, added to map of pointers in SimvascularGlobalArrayTransfer.cxx
      call PhAssignPointerDP(c_loc(a%emax), c_char_"Heart_EMax"//c_null_char)

      ! set pointer to TMax, added to map of pointers
      call PhAssignPointerDP(c_loc(a%tmax), c_char_"Heart_TMax"//c_null_char)

      ! set pointer to TRelax, added to map of pointers
      call PhAssignPointerDP(c_loc(a%trelax), c_char_"Heart_TRel"//c_null_char)
      
      ! set pointer to VLV^{t_{n}}, added to map of pointers
      call PhAssignPointerDP(c_loc(a%vlv_n), c_char_"Heart_LV_Vol"//c_null_char)
      
      end subroutine



!
! *** initialise internal variables in the heart model at t_{n}
!
      subroutine initxvars_hrt(a, stepnum)
      implicit none
!      
      class(numericalheart) :: a
      integer :: stepnum
      real*8 :: current_time, periodic_time
      integer :: periods, i, j
      real*8 :: time, elv      
      logical :: mvopen
      real*8 :: paortic, c_temp
      integer :: anum = 123
      integer :: bnum = 143
      integer :: cnum = 173
      integer :: dnum = 183
      integer :: enum = 103
!
!     ! convert current time to periodic time: 0 -> T
      current_time = delt*real(stepnum,8) 
      if (current_time .gt. a%period) then
         periods = floor(current_time/a%period)
         periodic_time = current_time - real(periods,8)*a%period
         time = periodic_time      
      else
         time = current_time
      end if
!

!     ! current elastance
      elv = getelastance(time, a%period, a%emax, a%emin, a%tmax, a%trelax)
!
!     ! if at step zero
      if (stepnum .eq. int(0)) then
!
!        ! initial aortic pressure set in pressure_n set in initreducedordermodel()
         paortic = a%pressure_n(1)  
!
!        ! lv pressure, here vlv_n set from edv 
         a%plv_n = elv*(a%vlv_n - a%vulv)
!
!        ! mitral valve check
!        ! if mv open calculate flow, assuming that dqdt = 0
         if (a%patrial .gt. a%plv_n) then
            mvopen = .true.
            a%qmv_n = (a%patrial - a%plv_n)/a%rmv
         else 
            mvopen = .false.
            a%qmv_n = real(0.0,8) 
         end if 
!
!        ! aortic valve check
!        ! if av open then determine pressure from flow relationship
         if ( a%plv_n .gt. paortic) then
            a%avopen = int(1) 
            a%updatepressure = int(0)
         else 
            a%avopen = int(0) 
            a%updatepressure = int(1)
            a%flow_n(1) = real(0.0,8)
         end if 
!
!     ! restart
      else
!
         open(cnum, file='VLV.dat', status='old')
         do i = 1, stepnum
            read(cnum,*) j, a%vlv_hist(i)
         end do
         close(cnum) 
         a%vlv_n = a%vlv_hist(stepnum)
!
         open(cnum, file='ELV.dat', status='old')
         do i = 1, stepnum
            read(cnum,*) j, a%elv_hist(i)
         end do
         close(cnum) 
!
         open(dnum, file='PLV.dat', status='old')
         do i = 1, stepnum
            read(dnum,*) j, a%plv_hist(i)
         end do
         close(dnum) 
         a%plv_n = a%plv_hist(stepnum)
!
         open(enum, file='QMV.dat', status='old')
         do i = 1, stepnum
            read(enum,*) j, a%qmv_hist(i)
         end do
         close(enum)          
         a%qmv_n = a%qmv_hist(stepnum)

!        ! initial aortic pressure set in pressure_n set in initreducedordermodel()
         paortic = a%pressure_n(1)  
!
!        ! mitral valve check
!        ! if mv open calculate flow, assuming that dqdt = 0
         if (a%patrial .gt. a%plv_n) then
            mvopen = .true.
            a%qmv_n = (a%patrial - a%plv_n)/a%rmv
         else 
            mvopen = .false.
            a%qmv_n = real(0.0,8) 
         end if
!
!        ! aortic valve check
!        ! if av open then determine pressure from flow relationship
         if ( a%plv_n .gt. paortic) then
            a%avopen = int(1) 
            a%updatepressure = int(0)
         elseif ( a%flow_n(1) .lt. real(0.0,8)) then
            a%avopen = int(1) 
            a%updatepressure = int(0)            
         else 
            a%avopen = int(0) 
            a%updatepressure = int(1)
         end if 

      end if
!
!
      return
      end subroutine initxvars_hrt
!
! *** iterate state variables for 
!
      subroutine iterate_hrt(a,stepn,varchar)
!     Note that this subroutine is only used if the heart is present
!     /without/ a systemic circuit.
      implicit none
      class(numericalheart) :: a
      integer :: stepn      
      character(len=*) :: varchar
      character(len=*), parameter :: solvechar = 'solve'
      character(len=*), parameter :: updatechar = 'update'      
!      
      real*8 :: time_n1  
      real*8 :: ptime         
      real*8 :: dtime
      real*8 :: elv_n1      
      real*8 :: x0(2)
      real*8 :: x(2)
      real*8 :: dx(2)
      real*8 :: j(2,2)
      real*8 :: r(2)
      real*8 :: eval
      real*8, parameter :: tol = real(1e-6,8)
      integer :: nperiods
      integer :: i
      integer :: counter = 0
      integer :: itermax = 10
      real*8 :: ps, dps
      real*8 :: f_backflow
!
!
!     ! solve or update
      if (varchar .eq. solvechar) then
         time_n1 = delt*(real(stepn,8)+alfi)
      elseif (varchar .eq. updatechar) then
!!         time_n1 = delt*real(stepn+1,8)
         time_n1 = delt*real(stepn,8) !! update called after lstep increase         
      end if
!
!     ! solve
      if (varchar .eq. solvechar) then
         if (a%plv_n .gt. a%pressure_n(1)) then 
            a%avopen = int(1) 
            a%updatepressure = int(0)
         elseif (a%flow_n(1) .lt. real(0.0,8)) then ! there is still forward flow
            a%avopen = int(1) 
            a%updatepressure = int(0)
         else
            a%avopen = int(0) 
            a%updatepressure = int(1)
         end if
!
         
         if (a%ibackflow) then
            if (a%backflow) then
               a%avopen = int(1)
            end if            
         end if 
!
      end if
!
!     ! periodic time

      ! TODO - These versions on ptime clash
      ! The different forms of ptime need an option flag   

      ! ptime = time_n1 - a%totalTimeOfAllCompletedPeriods
      ! ! Beacuse this subroutine looks ahead in time by dt or alfi*dt, it is possible that
      ! ! this look-ahead moves us to the next period , in which case we need to
      ! ! correct ptime:
      ! if (ptime .ge. a%period) then
      !    ptime = ptime - a%period
      ! end if
      
      if (time_n1 .gt. a%period) then
         nperiods = floor(time_n1/a%period)
         ptime = time_n1 - real(nperiods,8)*a%period
      else
         ptime = time_n1
      end if


      ! elastance
      elv_n1 = getelastance(ptime,a%period,a%emax,a%emin,a%tmax,a%trelax)
      ! write(*,*) 'periodic time:', ptime, a%activationTime * a%period !\todo remove
      ! write(*,*) varchar, a%activationTime * a%period
      
      ! finite difference time step
      if (varchar .eq. solvechar) then
         dtime = alfi*delt               
      elseif (varchar .eq. updatechar) then
         dtime = delt
      end if 

      ! jacobian
      j(:,:) = real(0.0,8)
      j(1,1) = real(1.0,8)
      j(1,2) = real(-1.0,8)*elv_n1*(real(1.0,8) + a%kelv*a%flow_n1(1)) ! flow -ve
      j(2,2) = real(1.0,8)/dtime

      ! initial state variables 
      x0(1) = a%pressure_n(1)
      x0(2) = a%vlv_n

      ! initial state value and evaluation
      x = x0
      eval = real(-1.0,8)

      ! newton solve only if valve is open    
      if (a%avopen .eq. int(1)) then

         ! zero state and counter  
         !x(:) = real(0.0,8)
         counter = int(0)
 
         do i = 1,itermax

            ! write(*,*) 'a%sPress',a%sPress
            ! write(*,*) 'x(1)',x(1)
            ! write(*,*) 'x(2)',x(2)
            ! write(*,*) 'elv_n1',elv_n1
            ! write(*,*) 'a%vulv',a%vulv
            ! write(*,*) 'a%kelv',a%kelv
            ! write(*,*) 'a%flow_n1(1)',a%flow_n1(1)
            ! write(*,*) 'a%rav',a%rav
            ! write(*,*) 'dtime',dtime  
            ! write(*,*) 'a%vlv_n',a%vlv_n

            r(1) = x(1) - a%sPress &                                            ! add stabilising pressure, this \Delta P acts to drive flow
                 - elv_n1*(x(2) - a%vulv)*(real(1.0,8) + a%kelv*a%flow_n1(1)) & ! flow -ve
                 - a%rav*a%flow_n1(1) 
            r(2) = x(2)/dtime - a%vlv_n/dtime - a%flow_n1(1)




            dx(:) = real(0.0,8)
            dx = solve_axb(j,real(-1.0,8)*r)  
            x = x + dx   

            eval = dx(1)*dx(1) + dx(2)*dx(2)
            eval = sqrt(eval)

            ! if evaluation is less than tolerance exit
            if (eval .lt. tol) then
               exit
            end if

            counter = counter + int(1)        
         
         end do
    
         ! calculate implicit coefficients for solve step
         if (varchar .eq. solvechar) then
         
            ! calculate pressure from flow, does not include stabilisation pressure as 
            ! this is already on the right hand side of the residual
            a%implicitcoeff(1,1) = a%rav &
                                 + elv_n1*(x(2) - a%vulv)*a%kelv
            a%implicitcoeff(1,2) = elv_n1*(x(2) - a%vulv)     
     
         end if     
      
         ! calculate implicit coefficients for update step
         if (varchar .eq. updatechar) then
         
            ! calculate pressure from flow, does not include stabilisation pressure as 
            ! this is already on the right hand side of the residual
            a%implicitcoeff_n1(1,1) = a%rav &
                                    + elv_n1*(x(2) - a%vulv)*a%kelv
            a%implicitcoeff_n1(1,2) = elv_n1*(x(2) - a%vulv) 
     
            ! set vlv from updated state
            a%vlv_n = x(2)

         end if 

      end if      
      
      ! store eval and counter
      a%evalval = eval
      a%evalcount = counter   

      return
      end subroutine iterate_hrt
!
! *** set implicit coefficients, although called each time step 
!     and always calculated the ibc bit-test sets if this is solved
!     or not 
!
      subroutine setimplicitcoeff_hrt(a,stepn,varchar) 
!
      implicit none
!
      class(numericalheart) :: a
      integer :: stepn      
      character(len=*) :: varchar
      character(len=*), parameter :: solvechar = 'solve'
      character(len=*), parameter :: updatechar = 'update'
      real*8 :: time_n1  
      integer :: nperiods
      real*8 :: ptime         
      real*8 :: dtime
      real*8 :: elv_n1
      integer :: i
      real*8 :: xmat(2,2)
      real*8 :: qvec(2)
      real*8 :: cvec(2)
      real*8 :: temp(2)
      real*8 :: qsol(2)
      real*8 :: csol(2)
      real*8 :: lrt

!
!     ! time at t_{n+alfi} [solve] or t_{n+1} [update]
      if (varchar .eq. solvechar) then
         time_n1 = delt*(real(stepn,8)+alfi)
      elseif (varchar .eq. updatechar) then
         time_n1 = delt*real(stepn+1,8)
      end if
!
!     ! if solve step, aortic valve check 
!     ! aortic pressure stored in pressure_n
      if (varchar .eq. solvechar) then
         if (a%plv_n .gt. a%pressure_n(1)) then 
            a%avopen = int(1) 
            a%updatepressure = int(0)
         elseif (a%flow_n(1) .lt. real(0.0,8)) then ! there is still forward flow
            a%avopen = int(1) 
            a%updatepressure = int(0)
         else
            a%avopen = int(0) 
            a%updatepressure = int(1)
         end if
      end if

!
!     ! convert time to periodic time
      if (time_n1 .gt. a%period) then
         nperiods = floor(time_n1/a%period)
         ptime = time_n1 - real(nperiods,8)*a%period         
      else
         ptime = time_n1
      end if      
!
!     ! elastance at ptime/time_n1
      elv_n1 = getelastance(ptime,a%period,a%emax,a%emin,a%tmax,a%trelax)
!!      elv_n1 = a%getelastance_senzaki(ptime, a%period, a%emax, a%emin, a%tmax,
!!     &                                a%trelax)             
!      
!     ! finite difference time step
      if (varchar .eq. solvechar) then
         dtime = alfi*delt               
      elseif (varchar .eq. updatechar) then
         dtime = delt
      end if 
!
      xmat(:,:) = real(0.0,8)
      qvec(:) = real(0.0,8)
      cvec(:) = real(0.0,8)
!
      xmat(1,1) = real(1.0,8)
      xmat(1,2) = real(-1.0,8)*elv_n1*(a%kelv*a%flow_n1(1) + real(1.0,8))
      xmat(2,2) = real(1.0,8)/dtime
!
!     ! here inflow is minus [!!]
      qvec(1) = a%rav + a%lav/dtime + elv_n1*a%kelv*a%vulv 
      qvec(2) = real(1.0,8)
!
!     ! minus flow as inflow is negative [!!]
      cvec(1) = -a%flow_n(1)*(a%lav/dtime) - elv_n1*a%vulv 
      cvec(2) = a%vlv_n/dtime
!
!     ! solve system
      qsol = solve_axb(xmat, qvec)
      csol = solve_axb(xmat, cvec)
!
!     ! set implicit coeffs
      if (varchar .eq. solvechar) then
         a%implicitcoeff(1,1) = qsol(1)
         a%implicitcoeff(1,2) = csol(1)
      elseif (varchar .eq. updatechar) then
!        ! inflow pressure coefficients for n+1
         a%implicitcoeff_n1(1,1) = qsol(1)
         a%implicitcoeff_n1(1,2) = csol(1)
!        ! ventricular volume coefficients for n+1
         a%vlv_coeff(1) = qsol(2)
         a%vlv_coeff(2) = csol(2)
      end if
!
      return
      end subroutine setimplicitcoeff_hrt     
!
!
! *** update variables in the heart model
!
      subroutine updxvars_hrt(a,lstep_passedIn)
      implicit none
      class(numericalheart) :: a
      integer :: lstep_passedIn
      logical :: mvopen
      real*8 :: time_n1
      real*8 :: current_time, periodic_time
      integer :: periods
      real*8 :: nelv, elv
      real*8 :: amat(2,2)
      real*8 :: bvec(2)
      real*8 :: xsol(2)
      real*8 :: xmat(2,2)
      real*8 :: qvec(2)
      real*8 :: cvec(2)
      real*8 :: qsol(2)
      real*8 :: csol(2)
      real*8 :: ps
!
!
!     ! check if mv is open using previous pressures at t = t_{n} 
!     ! or if there is forward flow at t = t_{n}
      if (a%patrial .gt. a%plv_n) then
         mvopen = .true.
      elseif (a%qmv_n .gt. real(0.0,8)) then
         mvopen = .true.
      else 
         mvopen = .false.
      end if
!
!     ! time at t = t_{n+1} (no longer at t = t_{n+alfi})
      time_n1 = delt*real(lstep_passedIn,8) 
!
!     ! periodic time
      ! TO DO - There are multiple definitions of periodic time - see iterate_hrt
      ! add a way to define this ???

      ! periodic_time = time_n1 - a%totalTimeOfAllCompletedPeriods

      current_time = time_n1 
      if (current_time .gt. a%period) then
         periods = floor(current_time/a%period)
         periodic_time = current_time - real(periods,8)*a%period
         ! time_n1 = periodic_time
      else
         ! time_n1 = current_time
         periodic_time = time_n1
      end if
!      
!     ! elv at t = t_{n+1}            
      elv = getelastance(periodic_time,a%period,a%emax,a%emin,a%tmax,a%trelax)

!     ! if mitral flow 
      if (mvopen) then 
!
         amat(:,:) = real(0.0,8)
         bvec(:) = real(0.0,8)
!
         amat(1,1) = 1/delt         ! update is at t = t_{n+1}
         amat(1,2) = real(-1.0,8)
         amat(2,1) = elv 
         amat(2,2) = a%rmv + a%lmv/delt
!
!        ! using previous vlv at t = t_{n}
         bvec(1) = a%vlv_n/delt
         bvec(2) = (a%lmv/delt)*a%qmv_n + a%patrial + elv*a%vulv
!
!        ! solve
         xsol = solve_axb(amat, bvec)
!
!        ! vlv and plv at t = t_{n+1}, stored in at t = t_{n} for next time step
         a%vlv_n = xsol(1)
         a%plv_n = elv*(xsol(1) - a%vulv)
!
!        ! qmv at t = t_{t+1}         
!!         if (xsol(2) .lt. real(0.0,8)) then
!!            a%qmv_n = real(0.0,8)
!!         else
            a%qmv_n = xsol(2)
!!         end if
!
!     ! else if aortic valve open
      elseif (a%avopen) then             
!        
         a%plv_n = elv*(a%vlv_n - a%vulv)         
         a%plv_n = a%plv_n*(real(1.0,8) + a%kelv*a%flow_n(1)) !! flow negative      
         a%plv_n = a%plv_n - a%sPress ! as this pressure is also driving flow forward
 
         a%qmv_n = real(0.0,8)
     
!
!        ! overwrite pressure with hydrostatic pressure 
         ! as pressure on face is 
         a%pressure_n(1) = a%implicitcoeff_n1(1,1)*a%flow_n(1)          &
                         + a%implicitcoeff_n1(1,2)                      &
                         - a%sPress ! needs to be subtracted as this is the pressure relative to PLV
!                                   ! and is included on this face

!         
      else ! both valves closed
!      
         a%vlv_n = a%vlv_n
         a%plv_n = elv*(a%vlv_n - a%vulv)
         a%qmv_n = real(0.0,8)         
!
      end if
!
!     ! store history variables
      a%vlv_hist(lstep_passedIn) = a%vlv_n
      a%elv_hist(lstep_passedIn) = elv
      a%plv_hist(lstep_passedIn) = a%plv_n
      a%qmv_hist(lstep_passedIn) = a%qmv_n
      a%flowhist(lstep_passedIn,1) = a%flow_n(1)
      a%pressurehist(lstep_passedIn,1) = a%pressure_n(1)
      a%eval_hist(lstep_passedIn,1) = a%evalcount 
      a%eval_hist(lstep_passedIn,2) = a%evalval
      a%stab_hist(lstep_passedIn) = a%sPress

!     ! check if av is open using new pressures and flows at t = t_{n+1}
      if (a%plv_n .gt. a%pressure_n(1)) then
         a%avopen = int(1) 
         a%updatepressure = int(0)
         if (a%ibackflow) then
            a%t_backflow = real(0.0,8)
         end if 
      elseif (a%flow_n(1) .lt. real(0.0,8)) then ! there is still forward flow
         a%avopen = int(1) 
         a%updatepressure = int(0)
      else
         a%avopen = int(0) 
         a%updatepressure = int(1)
      end if  
!      
      if (a%ibackflow) then
         if (a%flow_n(1) .gt. real(0.0,8)) then                     
            if (a%t_backflow < a%max_backflow) then            
               a%t_backflow = a%t_backflow + delt
               a%backflow = int(1)
               a%avopen = int(1)            
            else if (a%flow_n(1) > a%c_backflow) then            
               a%t_backflow = a%t_backflow + delt
               a%backflow = int(1)            
               a%avopen = int(1)            
            else             
               a%avopen = int(0)
               a%backflow = int(0)
            end if 
         else            
            a%backflow = int(0)
         end if          
      end if 
!
      return
      end subroutine      

   
!
! *** write internal variables to file
!
      subroutine writexvars_hrt(a,stepn)
      implicit none
!      
      include "mpif.h"                      
      class(numericalheart) :: a
      integer :: stepn
      integer :: rank
      integer :: ierr 
      integer :: anum = 619
      integer :: bnum = 740
      integer :: cnum = 478
      integer :: dnum = 874
      integer :: enum = 561
      integer :: i
      character(len=*), parameter :: format = '(i8,e20.10)'
      character(len=*), parameter :: format2 = '(i8,i8,e20.10)'      
!
      if (mod(stepn,ntout) .eq. int(0)) then
         call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)              
         if (rank .eq. int(0)) then
!
            open(anum, file='Paorta.dat', status='replace')
            do i = 1, stepn
               write(anum,format) i, a%pressurehist(i,1)
            end do
            close(anum)
!            
            open(bnum, file='Qaorta.dat', status='replace')
            do i = 1, stepn
               write(bnum,format) i, a%flowhist(i,1)
            end do
            close(bnum)
!
            open(cnum, file='VLV.dat', status='replace')
            do i = 1, stepn
               write(cnum,format) i, a%vlv_hist(i)
            end do
            close(cnum) 
!
            open(cnum, file='ELV.dat', status='replace')
            do i = 1, stepn
               write(cnum,format) i, a%elv_hist(i)
            end do
            close(cnum)             
!
            open(dnum, file='PLV.dat', status='replace')
            do i = 1, stepn
               write(dnum,format) i, a%plv_hist(i)
            end do
            close(dnum) 
!
            open(enum, file='QMV.dat', status='replace')
            do i = 1, stepn
               write(enum,format) i, a%qmv_hist(i)
            end do
            close(enum)                         
!
            open(enum, file='Eval.dat', status='replace')
            do i = 1, stepn
               write(enum,format2) i, int(a%eval_hist(i,1)), a%eval_hist(i,2)
            end do
            close(enum)                                  

            open(enum, file='PStab.dat', status='replace')
            do i = 1, stepn
               write(enum,format) i, a%stab_hist(i)
            end do
            close(enum)                                  

            
!
         end if
      end if
!     
      end subroutine 
!
! *** logical check if aortic valve is open
!
      function isavopen(a) 
      implicit none
      class(numericalheart) :: a
      integer :: isavopen
      isavopen = a%avopen    
      end function
!
!
! ********************************************************
! *** numerical rcr specific functions and subroutines ***
! ********************************************************
!
! *** nrcr constructor 
!
      function nrcrconstructor_params(surfnum,surflist,rcrval,pdmax,pdval,initrcr)    
      implicit none
      integer :: surfnum
      integer :: surflist(0:maxsurf)
      real*8  :: rcrval(surfnum,3)
      integer :: pdmax
      real*8 :: pdval(pdmax,2,surfnum) 
      integer :: initrcr

      type(numericalrcr) :: nrcrconstructor_params
!
      call nrcrconstructor_params%initialise_rcr(surfnum,          &
                                          surflist,         &
                                          rcrval,           &
                                          pdmax,            &
                                          pdval,            &
                                          initrcr)
!
      end function nrcrconstructor_params

      function nrcrconstructor_file(surfnum,surflist)    
      implicit none
      integer :: surfnum
      integer :: surflist(0:maxsurf)

      real*8  :: rcrval(surfnum,3)
      integer :: pdmax
      integer :: initrcr = 0
      integer :: numtimepoints_Pdist(surfnum)
      integer :: numDataRCR
      integer :: k, j, n
      real*8, allocatable :: pdval(:,:,:)
        
      type(numericalrcr) :: nrcrconstructor_file
       
      ! open rcrt.dat
      open(unit=818, file='rcrt.dat',status='old')
      read (818,*) pdmax

      allocate(pdval(pdmax,2,surfnum))

      do k = 1,surfnum

         read (818,*) numDataRCR
         numtimepoints_Pdist(k) = numDataRCR

         do j = 1,3
            read(818,*) rcrval(k,j) ! reads Rp,C,Rd in that order
         enddo
         
         do j = 1,numtimepoints_Pdist(k)
            read(818,*) (pdval(j,n,k), n=1,2) ! n=1 time, 2 value
         enddo

      enddo
      close(818)

      !PDist = pdval(1,2,1) ! read just the first value

      call nrcrconstructor_file%initialise_rcr(surfnum, surflist, rcrval, &
                                               pdmax, pdval, initrcr)
!
      end function nrcrconstructor_file
!
! *** initialise rcr
!
      subroutine initialise_rcr(this,surfnum,surflist,rcrcoeff,pdmax,pdval, &
                               initrcr)

      use iso_c_binding
      use datatypes, only: nrcr_states
      use phcommonvars, only: PhAssignPointerDP

      implicit none

      class(numericalrcr), intent(inout) :: this
      integer :: surfnum
      integer :: surflist(0:maxsurf)
      real*8 :: rcrcoeff(surfnum,3)
      integer :: pdmax, hstep
      real*8 :: pdval(pdmax,2,surfnum) ! this array is initialised as zero
!                                      ! then filled in, therefore values 
                                       ! that are 0 then the index is > 2 
                                       ! should be ignored  
      integer :: initrcr, ierr, fnum
      integer :: i, j ,k
!
      this%isactive = int(1)

      ! A label so we can identify this class (avoids a world of pain when working with a tower of derived types)
      this%classNameString = 'numericalrcr'
!      
!     ! set surface numbers and lists   
      this%surfnum = surfnum
      allocate(this%surfids(surfnum))
      this%surfids(1:surfnum) = surflist(1:surfnum)


      ! allocate history arrays
      if (lstep .gt. int(0)) then
         hstep = nstep + lstep
      else
         hstep = nstep
      end if      
      allocate(this%flowhist(hstep+1,surfnum))       
      allocate(this%pressurehist(hstep+1,surfnum))             

      ! zero history arrays
      this%flowhist(:,:) = real(0.0,8)
      this%pressurehist(:,:) = real(0.0,8)

      ! set flow and pressure file names
      write(this%flowfile,'(a)') 'QHistRCR.dat'
      write(this%pressurefile,'(a)') 'PHistRCR.dat'      

      ! allocate arrays for input parameters & data
      allocate(this%rcrparams(surfnum))
      allocate(this%parameters_RCR(3,surfnum)) ! testing

      ! allocate and zero other arrays 
      allocate(this%surfarea(surfnum))
      allocate(this%flow_n(surfnum))
      allocate(this%flow_n1(surfnum)) !! flow at n+1 / n+alf NOT USED
      allocate(this%pressure_n(surfnum))      
      allocate(this%implicitcoeff(surfnum,2)) 
      allocate(this%implicitcoeff_n1(surfnum,2)) 

      ! initialise reservoir pressure  
      if (initrcr) then
         ! set int
         this%init_pRes = 1
         allocate(this%pRes_n(surfnum))
         ! open rcr.x.dat
         open(fnum, file='rcrt.x.dat', status='old', iostat=ierr)         
         do i = 1, surfnum
            read(fnum,*) this%pRes_n(i)
         end do
         close(fnum)
      else 
         this%init_pRes = 0
      end if

      ! zero variables 
      this%surfarea(:) = real(0.0,8)
      this%flow_n(:) = real(0.0,8)
      this%flow_n1(:) = real(0.0,8)
      this%pressure_n(:) = real(0.0,8)
      this%implicitcoeff(:,:) = real(0.0,8)
      this%implicitcoeff_n1(:,:) = real(0.0,8)

      do i = 1, surfnum

         ! set rcr parameters
         this%rcrparams(i)%rp = rcrcoeff(i,1)
         this%rcrparams(i)%c = rcrcoeff(i,2)
         this%rcrparams(i)%rd = rcrcoeff(i,3)

         this%parameters_RCR(3,i) = rcrcoeff(i,3)
         this%parameters_RCR(1,i) = rcrcoeff(i,1)
         this%parameters_RCR(2,i) = rcrcoeff(i,2)

         ! count time points 
         k = int(1)
         do j = 2, pdmax
           if (pdval(j,1,i) .lt. pdval(j-1,1,i)) then
              exit
           else
              k = k + 1
           end if
         end do        
 
         ! allocate i'th entry with k x 2 size
         allocate(this%rcrparams(i)%pd%v(k,2))
     
         ! set values
         do j = 1, k
            this%rcrparams(i)%pd%v(j,1) = pdval(j,1,i)
            this%rcrparams(i)%pd%v(j,2) = pdval(j,2,i)
         end do

      end do

      end subroutine initialise_rcr      



     
!
! *** set implicit coefficients
!
! 
! *** get implicit coefficients
!
      subroutine setimplicitcoeff_rcr(a,stepn,varchar)
!
      use datatypes
      implicit none
!
      class(numericalrcr) :: a
      integer :: stepn
      character(len=*) :: varchar
      character(len=*), parameter :: solvechar = 'solve'
      character(len=*), parameter :: updatechar = 'update'      
!
      type(timedata) :: pd
      real*8 ::coeff(a%surfnum,2)    
      real*8 :: timen
      real*8 :: timen_1     
      real*8 :: pdistn
      real*8 :: pdistn_1
      real*8 :: rdn
      real*8 :: rdn_1
      real*8 :: rp
      real*8 :: c
      real*8 :: denom
      real*8 :: alfi_delt
!
      real*8 :: rd,pc_n1
!      
      integer :: i
!
!     ! time at stepn
      timen = delt*real(stepn,8)
!      
!     ! time and time step at solve/update
      if (varchar .eq. solvechar) then
         timen_1 = delt*(real(stepn,8)+alfi)         
         alfi_delt = alfi*delt
      elseif (varchar .eq. updatechar) then
         timen_1 = delt*real(stepn+1,8)
         alfi_delt = delt         
      end if
!
      if (a%init_pRes) then
!
         do i = 1, a%surfnum
!
            rdn_1 = a%rcrparams(i)%rd
            rp = a%rcrparams(i)%rp
            c  = a%rcrparams(i)%c

            pd = a%rcrparams(i)%pd
            pdistn = getvalue(timen,pd)
            pdistn_1 = getvalue(timen_1,pd)

            ! dirty hack for filtering
             rdn_1 = a%parameters_RCR(3,i)
             rp = a%parameters_RCR(1,i)
             c = a%parameters_RCR(2,i)
            ! rdn_1 = parameters_RCR(3,i)
            ! rp = parameters_RCR(1,i)
            ! c = parameters_RCR(2,i)    

             pdistn = a%parameters_Pd
             pdistn_1 = a%parameters_Pd

            denom = c/alfi_delt + 1/rdn_1
!
            coeff(i,1) = real(1.0,8) + rp*denom
            coeff(i,1) = coeff(i,1) / denom

            coeff(i,2) = pdistn_1/rdn_1 + (c/alfi_delt)*a%pRes_n(i)
            coeff(i,2) = coeff(i,2) / denom
!
         end do
!
      else
!
         do i = 1, a%surfnum

            ! parameters from derived type
            rdn = a%rcrparams(i)%rd
            rdn_1 = a%rcrparams(i)%rd

            rp = a%rcrparams(i)%rp
            c = a%rcrparams(i)%c

            pd = a%rcrparams(i)%pd
            pdistn = getvalue(timen,pd)
            pdistn_1 = getvalue(timen_1,pd)

            ! parameters overwritten
            ! dirty hack for filtering
             rdn_1 = a%parameters_RCR(3,i)
             rp = a%parameters_RCR(1,i)
             c = a%parameters_RCR(2,i)
            ! rdn_1 = parameters_RCR(3,i)
            ! rp = parameters_RCR(1,i)
            ! c = parameters_RCR(2,i)    

             pdistn = a%parameters_Pd
             pdistn_1 = a%parameters_Pd

            denom = real(1.0,8) + ((c*rdn_1)/alfi_delt)

            coeff(i,1) = rdn_1 + rp*(real(1.0,8) + ((c*rdn_1)/alfi_delt))
!
            coeff(i,2) = a%pressure_n(i) + pdistn_1 - pdistn - rp*a%flow_n(i)
            coeff(i,2) = ((c*rdn_1)/alfi_delt)*coeff(i,2)+ pdistn_1
!
            coeff(i,1) = coeff(i,1) / denom
            coeff(i,2) = coeff(i,2) / denom
!
         end do
!
      end if
!
!     ! set coefficients at solve/update
      if (varchar .eq. solvechar) then
         a%implicitcoeff(:,1) = coeff(:,1)
         a%implicitcoeff(:,2) = coeff(:,2)
      elseif (varchar .eq. updatechar) then
         a%implicitcoeff_n1(:,1) = coeff(:,1)
         a%implicitcoeff_n1(:,2) = coeff(:,2)         
      end if
!
      end subroutine setimplicitcoeff_rcr
!
! *** 
!
!       subroutine updxvars_rcr(a,lstep_passedIn)
! !
!       use datatypes
!       implicit none

!       class(numericalrcr) :: a      
!       integer :: i, lstep_passedIn
!       real*8 :: time_n1
!       real*8 :: r1, r2, c, pd
!       type(timedata) :: pdvals
!       real*8 :: denom, pc_n1
! !
! !     ! time at t = t_{n+1} (no longer at t = t_{n+alfi})
!       time_n1 = delt*real(lstep_passedIn,8) 
! !      
!       do i = 1, a%surfnum

!          ! parameters set from derived types
!          r1 = a%rcrparams(i)%rp
!          r2 = a%rcrparams(i)%rd            
!          c = a%rcrparams(i)%c

!          pdvals = a%rcrparams(i)%pd
!          pd = getvalue(time_n1,pdvals)

!          ! parameters overwritten
!          ! dirty hack for filtering
!           r2 = a%parameters_RCR(3,i)
!           r1 = a%parameters_RCR(1,i)
!           c = a%parameters_RCR(2,i)
!             ! r2 = parameters_RCR(3,i)
!             ! r1 = parameters_RCR(1,i)
!             ! c = parameters_RCR(2,i)    




!           pd = a%parameters_Pd

! !
!          denom = c/delt + real(1,8)/r2
!          pc_n1 = a%flow_n1(i) + a%pRes_n(i)*(c/delt) + pd/r2
!          a%pRes_n(i) = pc_n1/denom
! !         
!       end do
! !
!       end subroutine
!

      ! subroutine to update the pressure and flow history arrays
      subroutine updxvars_rcr(a,lstep_passedIn)

      implicit none

      class(numericalrcr) :: a      
      integer :: i, lstep_passedIn

      do i = 1, a%surfnum
         a%flowhist(lstep_passedIn,i) = a%flow_n(i)
         a%pressurehist(lstep_passedIn,i) = a%pressure_n(i)
      end do 
      
      end subroutine

      ! subroutine to write out the pressure and flow history arrays
      subroutine writexvars_rcr(a,stepn)
      
      use mpi, only: MPI_COMM_WORLD
      
      implicit none

      class(numericalrcr) :: a
      integer :: stepn
      integer :: rank
      integer :: ierr 
      integer :: anum = 619
      integer :: bnum = 740
      integer :: i
      character(len=100) :: dimchar             
      character(len=100) :: varsformat

      write(dimchar,'(i10)') a%surfnum
      write(varsformat,'(3(a))') '(i8,',trim(adjustl(dimchar)),'(e20.10))'

      if (mod(stepn,ntout) .eq. int(0)) then

         call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)              

         if (rank .eq. int(0)) then

            open(anum, file=a%pressurefile, status='replace')
            do i = 1, stepn
               write(anum,varsformat) i, a%pressurehist(i,:)
            end do
            close(anum)
            
            open(bnum, file=a%flowfile, status='replace')
            do i = 1, stepn
               write(bnum,varsformat) i, a%flowhist(i,1)
            end do
            close(bnum)

         end if

      end if

      end subroutine       




      subroutine assign_ptrs_ext_rcr(a)

      use iso_c_binding
      use phcommonvars, only: PhAssignPointerDP

      implicit none

      class(numericalrcr) :: a

      call PhAssignPointerDP(c_loc(a%parameters_RCR), c_char_"WindkesselRCR_Params"//c_null_char)
      call PhAssignPointerDP(c_loc(a%parameters_Pd), c_char_"WindkesselRCR_Pdist"//c_null_char)
      call PhAssignPointerDP(c_loc(a%pressure_n), c_char_"WindkesselRCR_P"//c_null_char)

      end subroutine


      subroutine setimplicitcoeff_netlistLPN(this,stepn,varchar)

         implicit none

         include "mpif.h"

         class(netlistLPN), intent(inout) :: this
         integer, intent(in) :: stepn
         character(len=*), intent(in) :: varchar
         character(len=5), parameter :: solvechar = 'solve'
         character(len=6), parameter :: updatechar = 'update'

         ! integer :: ierr
         ! integer :: rank

         real*8 :: alfi_delt
         real*8 :: coeff(this%surfnum,2)

         integer :: kk

         if (varchar.eq.solvechar) then
            alfi_delt = alfi*delt
         elseif (varchar.eq.updatechar) then
            alfi_delt = delt
         end if

         call this%generateLinearSystemFromPrescribedCircuit(alfi_delt)
         call this%assembleRHS_netlistLPN(stepn)

         ! call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
         ! if (rank .eq. int(0)) then
         !    call flush(6)
         !    write(*,*) this%systemMatrix
         !    call flush(6)
         !    stop
         ! end if
         ! call MPI_BARRIER(MPI_COMM_WORLD,ierr)
         do kk=1, this%numberOfLPNSurfaces
            this%inverseOfSystemMatrix(1:this%systemSize(kk),1:this%systemSize(kk),kk) = invertSquareMatrix(this%systemMatrix(1:this%systemSize(kk),1:this%systemSize(kk),kk))
            this%solutionVector(1:this%systemSize(kk),kk) = matmul(this%inverseOfSystemMatrix(1:this%systemSize(kk),1:this%systemSize(kk),kk), this%RHS(1:this%systemSize(kk),kk))
            ! this%solutionVector = solvematrixvector(this%systemMatrix,this%RHS)

            !\todo make this generic!
            coeff(kk,1) = this%inverseOfSystemMatrix(1,this%columnIndexOf3DInterfaceFlowInLinearSystem(kk),kk)
            coeff(kk,2) = this%solutionVector(1,kk) - this%inverseOfSystemMatrix(1,this%columnIndexOf3DInterfaceFlowInLinearSystem(kk),kk)*this%RHS(this%columnIndexOf3DInterfaceFlowInLinearSystem(kk),kk) !\todo make dynamic

            if (varchar.eq.solvechar) then
               this%implicitcoeff(kk,1) = coeff(kk,1)
               this%implicitcoeff(kk,2) = coeff(kk,2)
               ! write(*,*) 'varchar', varchar
               ! write(*,*) 'IC1', coeff(kk,1)
               ! write(*,*) 'IC2', coeff(kk,2)
               ! write(*,*) 'solv1', this%solutionVector(1)
               ! write(*,*) 'rhs', this%RHS
               ! write(*,*) 'PIM:', this%RHS(9)
               ! write(*,*) 'PIM_old:', this%RHS(14)
               ! write(*,*) 'P1:', this%RHS(11)
               ! write(*,*) 'P2:', this%RHS(13)
               ! write(*,*) 'row1 inversemat:', this%inverseOfSystemMatrix(1,:)
               ! write(*,*) 'alfidelt', alfi_delt
               ! write(*,*) 'alfi', alfi
               ! write(*,*) 'R_etc',this%circuitData
            elseif (varchar.eq.updatechar) then
               this%implicitcoeff_n1(kk,1) = coeff(kk,1)
               this%implicitcoeff_n1(kk,2) = coeff(kk,2)
               ! write(*,*) 'IC1b', coeff(kk,1)
               ! write(*,*) 'IC2b', coeff(kk,2)
            end if
         end do


      end subroutine setimplicitcoeff_netlistLPN

      subroutine updateLPN_netlistLPN(this)

         use calcFlowPressure, only: FlowHist
         implicit none

         class(netlistLPN) :: this
         integer :: ll
         integer :: kk
         integer :: indexShift

         do kk=1, this%numberOfLPNSurfaces
            this%RHS(this%columnIndexOf3DInterfaceFlowInLinearSystem(kk),kk) = this%flowpntr(kk)%p!\todo make this write to the correct entry of RHS, dynamically, and read the correct pointer when there are multiple netlist LPNs
            ! write(*,*) 'flow added: ', this%flowpntr(1)%p
            this%solutionVector(1:this%systemSize(kk),kk) = matmul(this%inverseOfSystemMatrix(1:this%systemSize(kk),1:this%systemSize(kk),kk), this%RHS(1:this%systemSize(kk),kk))

            do ll=1, this%numberOfPressureNodes(kk)
               this%pressuresInLPN(ll,kk) = this%solutionVector(ll,kk)
               ! write(*,*) 'pressure==> ', ll, ": ", this%pressuresInLPN(ll)
            end do

            indexShift = this%numberOfPressureNodes(kk) + this%numberOfHistoryPressures(kk)
            do ll=1, this%numberOfComponents(kk)
               this%flowsInLPN(ll,kk) = this%solutionVector(ll+indexShift,kk)
            end do
         end do

            ! write(*,*) 'discrepancy:', (-this%P_a(1) - this%pressuresInLPN(2))/1.2862d5 - this%flowsInLPN(1)

      end subroutine updateLPN_netlistLPN


      subroutine setimplicitcoeff_controlledCoronary(this,stepn,varchar)

         use datatypes
         implicit none

         class(controlledCoronaryModel), intent(inout) :: this
         integer, intent(in) :: stepn
         character(len=*), intent(in) :: varchar
         character(len=5), parameter :: solvechar = 'solve'
         character(len=6), parameter :: updatechar = 'update'
         real*8 ::coeff(this%surfnum,2)
         real*8 alfi_delt
         real*8 :: R_a
         real*8 :: R_p
         real*8 :: R_d
         real*8 :: C_a
         real*8 :: C_im
         real*8 :: m11                ! just a convenient constant
         real*8 :: m12                ! just a convenient constant
         real*8 :: m22                ! just a convenient constant
         real*8 :: temporary_variable
         real*8 :: P_IM_mid !\todo ensure that the numerical heart is active, to give us these pressures
         real*8 :: P_IM_mid_lasttimestep
         real*8 :: timen_1
         integer :: ii

!         write(*,*) 'set coronary implicits', this%surfnum !\todo remove
!         write(*,*) 'set coronary implicits', this%numberOfControlledCoronaryModelSurfaces !\todo remove

         ! Get the intramyocardial pressure for the previous time-step and previous-previous time-step
         ! note that these should both really be evaluated one time-step later, but to avoid the need
         ! for a simultaneous solve of the heart and coronaries, we do it this way for now
         ! \todo fix this
         if ((stepn .eq. int(0)).or.(stepn .eq. int(1))) then ! treat case with no known IM pressure yet
            P_IM_mid_lasttimestep = 5000d0 ! \todo find a better way of doing this; maybe input this value from file...
            P_IM_mid = 5000d0 ! ... or set it based on the aortic valve state at simulation start
         elseif (stepn .eq. int(2)) then ! treat case where only one IM pressure history point is known
            P_IM_mid_lasttimestep = this%intramyocardialPressureToLVScaling * hrt%plv_hist(stepn-1)
            P_IM_mid = this%intramyocardialPressureToLVScaling * hrt%plv_hist(stepn)
         else ! get the previous intramyocardial pressure in the case where we have enough doata for this (see comment before start of "if" block)
            P_IM_mid_lasttimestep = this%intramyocardialPressureToLVScaling * hrt%plv_hist(stepn-1) !\todo check these actually exist on first iteration
            ! Get IM pressure for now (this will be adjusted in a moment if we're on, according to alpha in gen alpha method)
            P_IM_mid = this%intramyocardialPressureToLVScaling * hrt%plv_hist(stepn) !\todo check these actually exist on first iteration
         end if


!        record these so they're available elsewhere (just updateLPN_coronary() at time-of-writing)
!        Note that these may want to become per-outlet at some point \todo
         this%P_IM_mid = P_IM_mid
         this%P_IM_mid_lasttimestep = P_IM_mid_lasttimestep

         ! a test:   \todo remove
#if EXTRA_CONSOLE_OUTPUT == 1
         write(*,*) 'plv_hist test:'
         write(*,*) hrt%plv_hist(1), stepn, hrt%plv_hist(stepn-1)
         write(*,*) stepn, P_IM_mid, P_IM_mid_lasttimestep
#endif

         ! time at stepn:
         if (varchar.eq.solvechar) then
            timen_1 = delt*(real(stepn,8)+alfi)
            alfi_delt = alfi*delt
!            P_IM_mid = this%intramyocardialPressureToLVScaling * hrt%plv_hist(stepn)
         elseif (varchar.eq.updatechar) then
            timen_1 = delt*real(stepn,8)
            alfi_delt = delt
!            P_IM_mid = this%intramyocardialPressureToLVScaling * hrt%plv_hist(stepn)
         end if

         do ii=1, this%surfnum

            R_a = this%ra(ii)
            R_p = this%rp(ii)
            R_d = this%rd(ii)
#if EXTRA_CONSOLE_OUTPUT == 1
            write(*,*) 'R_d from within setimplicitcoeff_controlledCoronary', R_d, ii !\todo remove
#endif
            C_a = this%c_a(ii)
            C_im = this%c_im(ii)


            ! we need these a couple of times to compute coeff, so precompute them for convenience.
            m11 = 1.0d0 + R_a*C_a/alfi_delt
            m12 = R_a * (1.0d0/R_d + C_im/alfi_delt)
            m22 = R_p * (C_im/alfi_delt + 1.0d0/R_d) + 1.0d0
            temporary_variable = m22/(m11*m22-m22+m12) !real(1.0,8) / (m11 - real(2.0,8))

            ! The actual differential equation "solve" - here we obtain the operators
            ! which describe the boundary condtion from the ODE.
            coeff(ii,1) = (m11 + m12/m22) * R_a * temporary_variable
            coeff(ii,2) = ((this%P_1(ii)*C_a + this%P_2(ii)*C_im + P_IM_mid*C_im - P_IM_mid_lasttimestep &
                           *C_im)*R_a/alfi_delt - m12/m22*R_p*C_im/alfi_delt           &
                           *(this%P_2(ii)+P_IM_mid-P_IM_mid_lasttimestep))*temporary_variable
            ! write(*,*) 'pim mid', P_IM_mid
            ! write(*,*) 'pim mid LTS', P_IM_mid_lasttimestep
            ! write(*,*) 'P1', this%P_1(ii)
            ! write(*,*) 'P2', this%P_2(ii)
#if EXTRA_CONSOLE_OUTPUT == 1
            write(*,*) '------=====>>>>> coeff', ii, coeff, m11,m12,m22,R_a, alfi_delt,alfi
            write(*,*) 'should be updated LPN pressure nodes:',this%P_1(ii)/1333,this%P_2(ii)/1333
#endif
         end do

#if EXTRA_CONSOLE_OUTPUT == 1
         write(*,*) 'av status:', hrt%avopen !\todo remove
#endif

         if (varchar.eq.solvechar) then
            this%implicitcoeff(:,1) = coeff(:,1)
            this%implicitcoeff(:,2) = coeff(:,2)
            ! write(*,*) 'varchar', varchar
            ! write(*,*) '->IC1', coeff(1,1)
            ! write(*,*) '->IC2', coeff(1,2)
            ! write(*,*) 'PIM:', P_IM_mid
            ! write(*,*) 'PIM_old:', P_IM_mid_lasttimestep
            ! write(*,*) 'P1:', this%P_1(1)
            ! write(*,*) 'P2:', this%P_2(1)
            ! write(*,*) 'ii', ii
            ! write(*,*) 'alfi', alfi
            ! write(*,*) 'alfidelt', alfi_delt
            ! write(*,*) 'R_a', R_a
            ! write(*,*) 'R_p', R_p
            ! write(*,*) 'R_d', R_d
            ! write(*,*) 'C_a', C_a
            ! write(*,*) 'C_im', C_im
         elseif (varchar.eq.updatechar) then
            this%implicitcoeff_n1(:,1) = coeff(:,1)
            this%implicitcoeff_n1(:,2) = coeff(:,2)
            ! write(*,*) '->IC1b', coeff(1,1)
            ! write(*,*) '->IC2b', coeff(1,2)
         end if

      end subroutine setimplicitcoeff_controlledCoronary

!     The coronary system-solve in setimplicitcoeff_controlledCoronary relies on knowing the 
!     values of the pressure at various points throughout the LPN model. These need updating
!     each step; this is done by the subroutine updateLPN_coronary.
      subroutine updateLPN_coronary(this,stepn)

         use datatypes
         implicit none

         class(controlledCoronaryModel), intent(inout) :: this
         integer, intent(in) :: stepn
         integer :: ii
         real*8 :: determinant
         real*8 :: P_a !\todo set this - maybe should be generated at the same time as the initial conditions?
         real*8 :: rhs(2)
         real*8 :: alfi_delt
         real*8 :: resultVector(2)
         real*8 :: inverseM(2,2)
         real*8 :: m11
         real*8 :: m12
         real*8 :: m22
         real*8 :: R_a
         real*8 :: R_p
         real*8 :: R_d
         real*8 :: C_a
         real*8 :: C_im

         alfi_delt = delt !\todo maybe it should be possible to multiply this by alfi too (i.e. if we call this subroutine on both solves and on updates)
                          ! Actually, I dont think you'd ever want to call it on a varchar=='solve' step; maybe rename this just "delt".

         do ii=1, this%surfnum

            R_a = this%ra(ii)
            R_p = this%rp(ii)
            R_d = this%rd(ii)
            C_a = this%c_a(ii)
            C_im = this%c_im(ii)

            P_a = -this%P_a(ii) ! this minus sign is correct!

            ! We're now going to solve the 2x2 system m*[P_1;P_2] = rhs.
            ! Define the LPN system matrix:
            m11 = real(1.0,8) + R_a*C_a/alfi_delt
            m12 = R_a * (1/R_d + C_im/alfi_delt)
            m22 = R_p * (C_im/alfi_delt + real(1.0,8)/R_d) + real(1.0,8)
            !the m21 entry of this matrix is just -1.
            determinant = m11*m22 + m12


            rhs(1) = P_a + (this%P_1(ii)*C_a + this%P_2(ii)*C_im + this%P_IM_mid*C_im - this%P_IM_mid_lasttimestep &
                     *C_im) * R_a/alfi_delt
            rhs(2) = R_p*C_im/alfi_delt * (this%P_2(ii)+this%P_IM_mid-this%P_IM_mid_lasttimestep)

            ! Find the inverse of this matrix
            inverseM(1,1) = m22 / determinant
            inverseM(1,2) = -m12 / determinant
            inverseM(2,1) = real(1.0,8) / determinant
            inverseM(2,2) = m11 / determinant

            ! Solve the system for the pressures within the LPN that we need to know:
            this%P_1(ii) = inverseM(1,1) * rhs(1) + inverseM(1,2) * rhs(2)
            this%P_2(ii) = inverseM(2,1) * rhs(1) + inverseM(2,2) * rhs(2)
         end do

         ! write(*,*) 'discrepancy:', (-this%P_a(1) - this%P_2(1))/1.2862e5 - this%flowpntr(1)%p
         ! this%P_1(1) = -(1.2862d5*this%flowpntr(1)%p + this%P_a(1))

      end subroutine updateLPN_coronary

      subroutine setSurfacePressure_coronary(this,pressure,coronarySurfaceIndex)

         implicit none

         class(controlledCoronaryModel) :: this
         real*8, intent(in) :: pressure
         integer, intent(in) :: coronarySurfaceIndex

         this%P_a(coronarySurfaceIndex) = pressure
         return
      end subroutine setSurfacePressure_coronary

      subroutine setSurfacePressure_netlistLPN(this,pressure,netlistSurfaceIndex)

         implicit none

         class(netlistLPN) :: this
         real*8, intent(in) :: pressure
         integer, intent(in) :: netlistSurfaceIndex

         !\todo this function currently doesn't do anything useful (the written pressuresInLPN value is never used). Maybe just remove this function.
         this%P_a(netlistSurfaceIndex) = pressure
         return
      end subroutine setSurfacePressure_netlistLPN

      ! for writing out restarts in restar.f
      integer function getNumberOfControlledCoronaryModelSurfaces(this)
         
         implicit none

         class(controlledCoronaryModel) :: this

         getNumberOfControlledCoronaryModelSurfaces = this%numberOfControlledCoronaryModelSurfaces
         
         return
      end function getNumberOfControlledCoronaryModelSurfaces

      ! for writing out restarts in restar.f
      real*8 function getRpByIndex(this,ii)

         implicit none

         class(controlledCoronaryModel) :: this
         integer :: ii

         getRpByIndex = this%rp(ii)

         return
      end function getRpByIndex

      subroutine setRpByIndex(this,ii,restart_rp)

         implicit none

         class(controlledCoronaryModel) :: this
         integer :: ii
         real*8 :: restart_rp

         this%rp(ii) = restart_rp

         return
      end subroutine setRpByIndex

      ! for writing out restarts in restar.f
      real*8 function getRdByIndex(this,ii)

         implicit none

         class(controlledCoronaryModel) :: this
         integer :: ii

         getRdByIndex = this%rd(ii)

         return
      end function getRdByIndex

      subroutine setRdByIndex(this,ii,restart_rd)

         implicit none

         class(controlledCoronaryModel) :: this
         integer :: ii
         real*8 :: restart_rd

         this%rd(ii) = restart_rd

         return
      end subroutine setRdByIndex

      ! for writing out restarts in restar.f
      real*8 function getP_1ByIndex(this,ii)

         implicit none

         class(controlledCoronaryModel) :: this
         integer :: ii

         getP_1ByIndex = this%P_1(ii)

         return
      end function getP_1ByIndex

      subroutine setP_1ByIndex(this,ii,restart_P_1)

         implicit none

         class(controlledCoronaryModel) :: this
         integer :: ii
         real*8 :: restart_P_1

         this%P_1(ii) = restart_P_1

         return
      end subroutine setP_1ByIndex

      ! for writing out restarts in restar.f
      real*8 function getP_2ByIndex(this,ii)

         implicit none

         class(controlledCoronaryModel) :: this
         integer :: ii

         getP_2ByIndex = this%P_2(ii)

         return
      end function getP_2ByIndex

      subroutine setP_2ByIndex(this,ii,restart_P_2)

         implicit none

         class(controlledCoronaryModel) :: this
         integer :: ii
         real*8 restart_P_2

         this%P_2(ii) = restart_P_2

         return
      end subroutine setP_2ByIndex

      ! for writing out restarts in restar.f
      real*8 function getcurrentMeanMyocardialOxygenHungerByIndex(this,ii)

         implicit none

         class(controlledCoronaryModel) :: this
         integer :: ii

         getcurrentMeanMyocardialOxygenHungerByIndex = this%currentMyocardialHungerSignal(ii)

         return
      end function getcurrentMeanMyocardialOxygenHungerByIndex

      subroutine setcurrentMeanMyocardialOxygenHungerByIndex(this,ii,restart_meanO2Discrepancy)

         implicit none

         class(controlledCoronaryModel) :: this
         integer :: ii
         real*8 :: restart_meanO2Discrepancy

         this%currentMyocardialHungerSignal(ii) = restart_meanO2Discrepancy

         return
      end subroutine setcurrentMeanMyocardialOxygenHungerByIndex

      ! for writing out restarts in restar.f
      real*8 function getMVO2previousByIndex(this,ii)

         implicit none

         class(controlledCoronaryModel) :: this
         integer :: ii

         getMVO2previousByIndex = this%MVO2previousDt(ii)

         return
      end function getMVO2previousByIndex

      subroutine setMVO2previousByIndex(this,ii,restart_MVO2previous)

         implicit none

         class(controlledCoronaryModel) :: this
         integer :: ii
         real*8 :: restart_MVO2previous

         this%MVO2previousDt(ii) = restart_MVO2previous

         return
      end subroutine setMVO2previousByIndex

      real*8 function getP_aByIndex(this,ii)

         implicit none

         class(controlledCoronaryModel) :: this
         integer :: ii

         getP_aByIndex = this%P_a(ii)

         return
      end function getP_aByIndex

      ! for writing out restarts, as called in restar.f
      subroutine writeO2supplyDemandHistory(this)

         implicit none

         class(controlledCoronaryModel) :: this
         integer :: recordLength

         ! compiler is set to take longword-sized chunks. These are 4 bytes each,
         ! a real*8 is 8 bytes, so we need 2 longwords for each real*8:
         ! This is adjustable with the compiler flag -assume byterecl , which
         ! would make the chunk size equal to one byte, instead of one longword.
         ! \todo worry about endianness for other machines!
         recordLength = 2*(nstep + maxval(this%O2supplyDemandHistoryWindowLength_timesteps) + 1)*(this%numberOfControlledCoronaryModelSurfaces)

         open(unit=74,file='O2supplyDemandHistory.bin',status='replace',form='unformatted',access='direct',recl=recordLength)
         write(74,rec=1) this%O2supplyDemandHistory
         close(74)

      end subroutine writeO2supplyDemandHistory

      subroutine loadO2supplyDemandHistory(this)

         implicit none

         class(controlledCoronaryModel) :: this
         integer :: recordLength

         ! compiler is set to take longword-sized chunks. These are 4 bytes each,
         ! a real*8 is 8 bytes, so we need 2 longwords for each real*8:
         ! This is adjustable with the compiler flag -assume byterecl , which
         ! would make the chunk size equal to one byte, instead of one longword.
         ! \todo worry about endianness for other machines!
         recordLength = 2*(nstep + maxval(this%O2supplyDemandHistoryWindowLength_timesteps) + 1)*(this%numberOfControlledCoronaryModelSurfaces)

         open(unit=75,file='O2supplyDemandHistory.bin',status='old',form='unformatted',access='direct',recl=recordLength)
         read(75,rec=1) this%O2supplyDemandHistory

         close(75)

      end subroutine loadO2supplyDemandHistory

      subroutine writeInternalLPNPressuresAndFlows(this,lstep_passedIn)
         use calcFlowPressure, only: FlowHist
         implicit none
         include "mpif.h"
         
         class(netlistLPN), intent(in) :: this
         integer, intent(in) :: lstep_passedIn
         integer :: ii
         integer :: rank
         integer :: ierr
         logical :: exist

         integer :: jj
         integer :: kk
         character(len=14) :: basename_pressure
         character(len=20) :: basename_flow
         character(len=10) :: ofSurfaceWord
         character(len=3) :: internalIndex
         character(len=3) :: surfaceIndex
         character(len=60) :: fullname
         character(len=4) :: dotdat

         ofSurfaceWord = '_ofSurface'
         dotdat = '.dat'

         call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)           
         if (rank .eq. int(0)) then
            basename_pressure = 'PressureAtNode'
            basename_flow = 'FlowThroughComponent'
            do kk=1, this%numberOfLPNSurfaces
               write(surfaceIndex,'(I3)') kk
               do jj=1, this%numberOfPressureNodes(kk)
                  write(internalIndex,'(I3)') jj
                  fullname = basename_pressure//trim(adjustl(internalIndex))//ofSurfaceWord//trim(adjustl(surfaceIndex))//dotdat
                  inquire(file=fullname,exist=exist)
                  if(exist) then
                     open(unit=97, file=fullname, status='old',position='append',action='write')
                  else
                     open(97, file=fullname,status='new',action='write')
                  end if
                  write(97,*) lstep_passedIn, this%pressuresInLPN(jj,kk)
                  close(97)
               end do

               ! Write flows
               do jj=1, this%numberOfComponents(kk)
                  write(internalIndex,'(I3)') jj
                  fullname = basename_flow//trim(adjustl(internalIndex))//ofSurfaceWord//trim(adjustl(surfaceIndex))//dotdat
                  inquire(file=fullname,exist=exist)
                  if(exist) then
                     open(unit=97, file=fullname, status='old',position='append',action='write')
                  else
                     open(97, file=fullname,status='new',action='write')
                  end if
                  write(97,*) lstep_passedIn, this%flowsInLPN(jj,kk)
                  close(97)
               end do
            end do


            ! inquire(file='r_p_history.dat',exist=exist)
            ! if(exist) then
            !    open(unit=97, file='r_p_history.dat',status='old',position='append',action='write')
            ! else
            !    open(97, file='r_p_history.dat',status='new',action='write')
            ! end if
            ! write(97,*) (this%rp(ii),ii=1,this%numberOfControlledCoronaryModelSurfaces) !\todo make sure this works for multiple coronaries (i.e. that it writes all the resistances, not just one. Check corresponding code for other output variables below, too.)
            ! close(97)

            ! inquire(file='r_d_history.dat',exist=exist)
            ! if(exist) then
            !    open(unit=97, file='r_d_history.dat',status='old',position='append',action='write')
            ! else
            !    open(97, file='r_d_history.dat',status='new',action='write')
            ! end if
            ! write(97,*) (this%rd(ii),ii=1,this%numberOfControlledCoronaryModelSurfaces)
            ! close(97)

            ! inquire(file='MVO2_history.dat',exist=exist)
            ! if(exist) then
            !    open(unit=97, file='MVO2_history.dat',status='old',position='append',action='write')
            ! else
            !    open(97, file='MVO2_history.dat',status='new',action='write')
            ! end if
            ! write(97,*) (this%MVO2(ii),ii=1,this%numberOfControlledCoronaryModelSurfaces)
            ! close(97)
            
            ! inquire(file='meanO2Discrepancy_history.dat',exist=exist)
            ! if(exist) then
            !    open(unit=97, file='meanO2Discrepancy_history.dat',status='old',position='append',action='write')
            ! else
            !    open(97, file='meanO2Discrepancy_history.dat',status='new',action='write')
            ! end if
            ! write(97,*) (this%currentMyocardialHungerSignal(ii),ii=1,this%numberOfControlledCoronaryModelSurfaces)
            ! close(97)

            ! inquire(file='FlowDebug_history.dat',exist=exist)
            ! if(exist) then
            !    open(unit=97, file='FlowDebug_history.dat',status='old',position='append',action='write')
            ! else
            !    open(97, file='FlowDebug_history.dat',status='new',action='write')
            ! end if
            ! write(97,*) (FlowHist(lstep_passedIn,this%localToGlobalSurfaceIndexMap(ii)),ii=1,this%numberOfControlledCoronaryModelSurfaces)
            ! close(97)

            ! inquire(file='P1_hist.dat',exist=exist)
            ! if(exist) then
            !    open(unit=97, file='P1_hist.dat',status='old',position='append',action='write')
            ! else
            !    open(97, file='P1_hist.dat',status='new',action='write')
            ! end if
            ! write(97,*) (this%P_1(ii),ii=1,this%numberOfControlledCoronaryModelSurfaces)
            ! close(97)

            ! inquire(file='P2_hist.dat',exist=exist)
            ! if(exist) then
            !    open(unit=97, file='P2_hist.dat',status='old',position='append',action='write')
            ! else
            !    open(97, file='P2_hist.dat',status='new',action='write')
            ! end if
            ! write(97,*) (this%P_2(ii),ii=1,this%numberOfControlledCoronaryModelSurfaces)
            ! close(97)
         end if

         return
      end subroutine writeInternalLPNPressuresAndFlows

      ! This writes the resistances, O2 supply discrepancies, etc.
      ! to disk for debug purposes.
      subroutine writeDebugHistories(this,lstep_passedIn)
         use calcFlowPressure, only: FlowHist
         implicit none
         include "mpif.h"
         
         class(controlledCoronaryModel), intent(in) :: this
         integer, intent(in) :: lstep_passedIn
         integer :: ii
         integer :: rank
         integer :: ierr
         logical :: exist

         call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)           
         if (rank .eq. int(0)) then
            inquire(file='r_p_history.dat',exist=exist)
            if(exist) then
               open(unit=97, file='r_p_history.dat',status='old',position='append',action='write')
            else
               open(97, file='r_p_history.dat',status='new',action='write')
            end if
            write(97,*) (this%rp(ii),ii=1,this%numberOfControlledCoronaryModelSurfaces) !\todo make sure this works for multiple coronaries (i.e. that it writes all the resistances, not just one. Check corresponding code for other output variables below, too.)
            close(97)

            inquire(file='r_d_history.dat',exist=exist)
            if(exist) then
               open(unit=97, file='r_d_history.dat',status='old',position='append',action='write')
            else
               open(97, file='r_d_history.dat',status='new',action='write')
            end if
            write(97,*) (this%rd(ii),ii=1,this%numberOfControlledCoronaryModelSurfaces)
            close(97)

            inquire(file='MVO2_history.dat',exist=exist)
            if(exist) then
               open(unit=97, file='MVO2_history.dat',status='old',position='append',action='write')
            else
               open(97, file='MVO2_history.dat',status='new',action='write')
            end if
            write(97,*) (this%MVO2(ii),ii=1,this%numberOfControlledCoronaryModelSurfaces)
            close(97)
            
            inquire(file='meanO2Discrepancy_history.dat',exist=exist)
            if(exist) then
               open(unit=97, file='meanO2Discrepancy_history.dat',status='old',position='append',action='write')
            else
               open(97, file='meanO2Discrepancy_history.dat',status='new',action='write')
            end if
            write(97,*) (this%currentMyocardialHungerSignal(ii),ii=1,this%numberOfControlledCoronaryModelSurfaces)
            close(97)

            inquire(file='FlowDebug_history.dat',exist=exist)
            if(exist) then
               open(unit=97, file='FlowDebug_history.dat',status='old',position='append',action='write')
            else
               open(97, file='FlowDebug_history.dat',status='new',action='write')
            end if
            write(97,*) (FlowHist(lstep_passedIn,this%localToGlobalSurfaceIndexMap(ii)),ii=1,this%numberOfControlledCoronaryModelSurfaces)
            close(97)

            inquire(file='P1_hist.dat',exist=exist)
            if(exist) then
               open(unit=97, file='P1_hist.dat',status='old',position='append',action='write')
            else
               open(97, file='P1_hist.dat',status='new',action='write')
            end if
            write(97,*) (this%P_1(ii),ii=1,this%numberOfControlledCoronaryModelSurfaces)
            close(97)

            inquire(file='P2_hist.dat',exist=exist)
            if(exist) then
               open(unit=97, file='P2_hist.dat',status='old',position='append',action='write')
            else
               open(97, file='P2_hist.dat',status='new',action='write')
            end if
            write(97,*) (this%P_2(ii),ii=1,this%numberOfControlledCoronaryModelSurfaces)
            close(97)
         end if

         return
      end subroutine writeDebugHistories


!!c 
!!c *** get implicit coefficients
!!c
!!      function getimplicitcoeff_rcr(a,stepn) result(coeff)
!!c
!!      use datatypes
!!      implicit none
!!c
!!      class(numericalrcr) :: a
!!      integer :: stepn
!!c
!!      type(timedata) :: pd
!!      real*8 ::coeff(a%surfnum,2)    
!!      real*8 :: timen
!!      real*8 :: timen_1     
!!      real*8 :: pdistn
!!      real*8 :: pdistn_1
!!      real*8 :: rdn
!!      real*8 :: rdn_1
!!      real*8 :: rp
!!      real*8 :: c
!!      real*8 :: denom
!!      real*8 :: alfi_delt
!!      integer :: i
!!c
!!      timen = delt*real(stepn,8)
!!      timen_1 = delt*(real(stepn,8)+alfi)
!!c
!!      alfi_delt = alfi*delt
!!c
!!      do i = 1, a%surfnum
!!c
!!         rdn = a%rcrparams(i)%rd
!!         rdn_1 = a%rcrparams(i)%rd
!!c
!!         rp = a%rcrparams(i)%rp
!!         c = a%rcrparams(i)%c
!!c
!!         pd = a%rcrparams(i)%pd
!!         pdistn = getvalue(timen,pd)
!!         pdistn_1 = getvalue(timen_1,pd)
!!c
!!         denom = real(1.0,8) + ((c*rdn_1)/alfi_delt)
!!c
!!         coeff(i,1) = rdn_1 + rp*(real(1.0,8) + ((c*rdn_1)/alfi_delt))
!!c
!!         coeff(i,2) = a%pressure_n(i) + pdistn_1 - pdistn - rp*a%flow_n(i)
!!         coeff(i,2) = ((c*rdn_1)/alfi_delt)*coeff(i,2)+ pdistn_1
!!c
!!         coeff(i,1) = coeff(i,1) / denom
!!         coeff(i,2) = coeff(i,2) / denom
!!c
!!      end do
!!c
!!      end function getimplicitcoeff_rcr
!
! *********************************************************
! *** numerical trcr specific functions and subroutines ***
! *********************************************************
!
! *** constructor
!
      function ntrcrconstructor(surfnum,surflist,nptsmax,trcrvals)
!
      implicit none
!
      type(numericaltrcr) :: ntrcrconstructor
      integer :: surfnum
      integer :: surflist(0:maxsurf)
      integer :: nptsmax
      real*8 :: trcrvals(nptsmax,5,surfnum)
!
      call ntrcrconstructor%initialise_trcr(surfnum,        &
                                            surflist,       &
                                            nptsmax,        &
                                            trcrvals)
!
      end function ntrcrconstructor

!
! *** initialise time varying rcr
!
      subroutine initialise_trcr(this,surfnum,surflist,nptsmax,trcrvals)
!
      implicit none
!
      class(numericaltrcr), intent(inout) :: this
      integer :: surfnum
      integer :: surflist(0:maxsurf)
      integer :: nptsmax
      real*8 :: trcrvals(nptsmax,5,surfnum)
      integer :: i
      integer :: j
      integer :: k
!
      this%isactive = int(1)

      ! A label so we can identify this class (avoids a world of pain when working with a tower of derived types)
      this%classNameString = 'numericaltrcr'
!      
!     ! set surface numbers and lists   
      this%surfnum = surfnum
      allocate(this%surfids(surfnum))
      this%surfids(1:surfnum) = surflist(1:surfnum)
!
!     ! allocate arrays for input parameters & data
      allocate(this%trcrparams(surfnum))
!
!     ! allocate and zero other arrays 
      allocate(this%surfarea(surfnum))
      allocate(this%flow_n(surfnum))
      allocate(this%flow_n1(surfnum))
      allocate(this%pressure_n(surfnum))
      allocate(this%implicitcoeff(surfnum,2))
      allocate(this%implicitcoeff_n1(surfnum,2))
!
!     ! zero variables 
      this%surfarea(:) = real(0.0,8)
      this%flow_n(:) = real(0.0,8)
      this%flow_n1(:) = real(0.0,8)
      this%pressure_n(:) = real(0.0,8)
      this%implicitcoeff(:,:) = real(0.0,8)
      this%implicitcoeff_n1(:,:) = real(0.0,8)
!
      do i = 1, surfnum
!
!        ! count time points 
         k = int(1)
         do j = 2, nptsmax
           if (trcrvals(j,1,i) .lt. trcrvals(j-1,1,i)) then
              exit
           else
              k = k + 1
           end if
         end do        
!
         allocate(this%trcrparams(i)%rp%v(k,2))
         allocate(this%trcrparams(i)%c%v(k,2))
         allocate(this%trcrparams(i)%rd%v(k,2))
         allocate(this%trcrparams(i)%pd%v(k,2))
!   
!        ! set values
         do j = 1, k
!
            this%trcrparams(i)%rp%v(j,1) = trcrvals(j,1,i)
            this%trcrparams(i)%rp%v(j,2) = trcrvals(j,2,i)
!
            this%trcrparams(i)%c%v(j,1) = trcrvals(j,1,i)
            this%trcrparams(i)%c%v(j,2) = trcrvals(j,3,i)
!
            this%trcrparams(i)%rd%v(j,1) = trcrvals(j,1,i)
            this%trcrparams(i)%rd%v(j,2) = trcrvals(j,4,i)
!
            this%trcrparams(i)%pd%v(j,1) = trcrvals(j,1,i)
            this%trcrparams(i)%pd%v(j,2) = trcrvals(j,5,i)
!
         end do
!
      end do
!
      write (this%flowfile,'(a)') 'QHistTRCR.dat'
      write (this%pressurefile,'(a)') 'PHistTRCR.dat'
!
!!      allocate(this%variablesfile(4))
!!      write(this%variablesfile(1),'(a)') 'ValueTRCRRp.dat'
!!      write (this%variablesfile(2),'(a)') 'ValueTRCRC.dat'
!!      write (this%variablesfile(3),'(a)') 'ValueTRCRRd.dat'
!!      write (this%variablesfile(4),'(a)') 'ValueTRCRPd.dat'
!!c      
!!      allocate(this%currparams(surfnum))
!!      this%currparams(:)%rp = real(0.0,8)
!!      this%currparams(:)%c = real(0.0,8)
!!      this%currparams(:)%rd = real(0.0,8)
!!      this%currparams(:)%pd = real(0.0,8)
!!c
!!      allocate(this%initparams(surfnum))
!!      do i = 1, this%surfnum
!!         this%initparams(i)%rp = this%trcrparams(i)%rp%v(1,2)
!!         this%initparams(:)%c = this%trcrparams(i)%c%v(1,2)
!!         this%initparams(:)%rd = this%trcrparams(i)%rd%v(1,2)
!!         this%initparams(:)%pd = this%trcrparams(i)%pd%v(1,2)
!!      end do
!
      allocate(this%isfeedback(surfnum))
      this%isfeedback(:) = .false.
!
      allocate(this%feedbackcontrol(surfnum))
!
      end subroutine
!
! *** get implicit coefficients
!
      subroutine setimplicitcoeff_trcr(a,stepn,varchar)

      use datatypes
      implicit none
!
      class(numericaltrcr) :: a
      integer :: stepn
      character(len=*) :: varchar
      character(len=*), parameter :: solvechar = 'solve'
      character(len=*), parameter :: updatechar = 'update'   
!
      type(timedata) :: rp
      type(timedata) :: c
      type(timedata) :: rd
      type(timedata) :: pd

!!      real*8 :: setimplicitcoeff_trcr(a%surfnum,2)
      real*8 :: coeff(a%surfnum,2)        
      real*8 :: timen
      real*8 :: timen_1     
      real*8 :: rpn_1
      real*8 :: cn_1
      real*8 :: rdn_1
      real*8 :: pdistn
      real*8 :: pdistn_1
      real*8 :: denom
      real*8 :: alfi_delt 
      integer :: i
!
!       timen = delt*real(stepn,8) ! \todo remove
!       timen_1 = delt*(real(stepn,8)+alfi) ! \todo remove
!
      ! alfi_delt = alfi*delt! \todo remove

      ! time at stepn
      timen = delt*real(stepn,8)

!     ! time and time step at solve/update
      if (varchar .eq. solvechar) then
         timen_1 = delt*(real(stepn,8)+alfi)         
         alfi_delt = alfi*delt
      elseif (varchar .eq. updatechar) then
         timen_1 = delt*real(stepn+1,8)
         alfi_delt = delt         
      end if

!
      do i = 1, a%surfnum
!
         rp = a%trcrparams(i)%rp
         rpn_1 = getvalue(timen_1,rp)
!
         c = a%trcrparams(i)%c
         cn_1 = getvalue(timen_1,c)
!
         rd = a%trcrparams(i)%rd
         rdn_1 = getvalue(timen_1,rd)

         if (a%isfeedback(i)) then
            rdn_1 = rdn_1*a%feedbackcontrol(i)
         else

         end if      
!
         pd = a%trcrparams(i)%pd
         pdistn = getvalue(timen,pd)
         pdistn_1 = getvalue(timen_1,pd)         
!
         denom = real(1.0,8) + ((cn_1*rdn_1)/alfi_delt)
!
         coeff(i,1) = rdn_1 + rpn_1*(real(1.0,8) + ((cn_1*rdn_1)/alfi_delt))
!
         coeff(i,2) = a%pressure_n(i) + pdistn_1 - pdistn - rpn_1*a%flow_n(i)
         coeff(i,2) = ((cn_1*rdn_1)/alfi_delt)*coeff(i,2) + pdistn_1

         coeff(i,1) = coeff(i,1) / denom
         coeff(i,2) = coeff(i,2) / denom
!
      end do

      ! set coefficients at solve/update
      if (varchar .eq. solvechar) then
         a%implicitcoeff(:,1) = coeff(:,1)
         a%implicitcoeff(:,2) = coeff(:,2)
      elseif (varchar .eq. updatechar) then
         a%implicitcoeff_n1(:,1) = coeff(:,1)
         a%implicitcoeff_n1(:,2) = coeff(:,2)
      end if

      end subroutine setimplicitcoeff_trcr
!
! *** set variables feedback control at step n+1
!
      subroutine setfeedbackcontrol(a,controlparameter)
!
      implicit none
      class(numericaltrcr) :: a
      real*8 :: controlparameter
      integer :: i
!
      do i = 1, a%surfnum
         if (a%isfeedback(i)) then         
            a%feedbackcontrol(i) = controlparameter
         end if
      end do
!
      return
      end subroutine setfeedbackcontrol
!
!!c *** load trcr variables 
!!c
!!      subroutine loadtrcrfiles(a,torow)
!!c
!!      use datatypes
!!      implicit none
!!c
!!      class(numericaltrcr) :: a
!!      integer :: torow
!!      real*8 :: datafile(torow,a%surfnum)
!!      integer :: unitnum = 450
!!      real*8 :: temp     
!!      integer :: i
!!      integer :: j
!!      integer :: ierr
!!      integer :: fnum
!!c       
!!      if (torow .gt. int(0)) then
!!c
!!         do fnum = 1, 4     
!!            open(unitnum, file=a%variablesfile(i), status='old', iostat=ierr)
!!            read(unitnum,*) temp
!!            do i = 1, torow
!!               do j = 1, a%surfnum
!!                  read(unitnum,*) datafile(i,j)
!!               end do
!!            end do
!!            close(unitnum)
!!c      
!!            select case (fnum)
!!            case (1) 
!!               a%currparams(1:a%surfnum)%rp = datafile(torow,1:a%surfnum) 
!!            case (2)
!!               a%currparams(1:a%surfnum)%c = datafile(torow,1:a%surfnum) 
!!            case (3)
!!               a%currparams(1:a%surfnum)%rd = datafile(torow,1:a%surfnum) 
!!            case (4)
!!               a%currparams(1:a%surfnum)%pd = datafile(torow,1:a%surfnum) 
!!            end select 
!!         end do
!!c      
!!      else
!!c
!!         do i = 1, a%surfnum
!!            a%currparams(i)%rp = getvalue(real(0.0,8),a%trcrparams(i)%rp)
!!            a%currparams(i)%c = getvalue(real(0.0,8),a%trcrparams(i)%c)
!!            a%currparams(i)%rd = getvalue(real(0.0,8),a%trcrparams(i)%rd)
!!            a%currparams(i)%pd = getvalue(real(0.0,8),a%trcrparams(i)%pd)
!!         end do
!!c
!!      end if
!!c
!!      return
!!      end subroutine
!
! *** add feedback surfaces
!
      subroutine setfeedbacksurfs(a,surflist)
!
      implicit none
!
      class(numericaltrcr) :: a
      integer :: surflist(0:maxsurf)
      integer :: i
      integer :: j
!
      do i = 1, a%surfnum
         do j = 1, maxsurf        
            if (a%surfids(i) .eq. surflist(j)) then
               a%isfeedback(i) = .true.

               exit
            end if
         end do
      end do
!
      return
      end subroutine

!

!
! *******************************************
! *** matrix and vector utility functions ***
! *** NB dependent upon intel mkl library ***
! *******************************************
!
!     ! matrix matrix multiply
      function matrixmatrix(amatrix,bmatrix) result(cmatrix)	
      implicit none
      real*8, dimension(:,:), intent(in) :: amatrix
      real*8, dimension(:,:), intent(in) :: bmatrix
      real*8, dimension(:,:) :: cmatrix(size(amatrix,1),size(bmatrix,2))
      integer :: m 
      integer :: n 
      integer :: p 
      real*8  :: alpha 
      real*8  :: beta 
!
!     ! check length of matrix and vector match
      if (size(amatrix,2) .ne. size(bmatrix,1)) then 
         write(*,*) ' *** WARNING size mismatch between matrix and matrix ***'
         stop
      end if 
!
      m = size(amatrix,1)
      n = size(amatrix,2)
      p = size(bmatrix,2)
      alpha = real(1.0,8)
      beta = real(0.0,8)
!
      call dgemm('N',      &
                 'N',      &
                  m,       &
                  n,       &
                  p,       &
                  alpha,	&
                  amatrix, &
                  m,       &
                  bmatrix, &
                  p,       &
                  beta,    &
                  cmatrix, &
                  m)
!
      return
      end function matrixmatrix
!
!     ! matrix vector multiply
      function matrixvector(amatrix,bvector) result(cvector)	
      implicit none
      real*8, dimension(:,:), intent(in) :: amatrix
      real*8, dimension(:), intent(in)   :: bvector
      real*8, dimension(:)   :: cvector(size(amatrix,1))
      integer :: m 
      integer :: n 
      integer :: p 
      real*8  :: alpha 
      real*8  :: beta 
!
!     ! check length of matrix and vector match
      if (size(amatrix,2) .ne. size(bvector,1)) then 
         write(*,*) ' *** WARNING size mismatch between matrix and vector ***'
         stop
      end if 
!
      m = size(amatrix,1)
      n = int(1) 
      p = size(bvector,1)
      alpha = real(1.0,8)
      beta = real(0.0,8)
!
      call dgemm('N',      &
                 'N',      &
                  m,       &
                  n,       &
                  p,       &
                  alpha,	&
                  amatrix, &
                  m,       &
                  bvector, &
                  p,       &
                  beta,    &
                  cvector, &
                  m)       
!
      return
      end function matrixvector
!
! *** wrapper for dgemm to to multiply a matrix with a vector (pointer)
! 
      function matrixvectorpointer(aamatrix,bbvector) result(ccvector)	
      implicit none
      real*8, dimension(:,:), intent(in) :: aamatrix
      type(pntr), dimension(:), intent(in) :: bbvector
      real*8, dimension(:) :: ccvector(size(aamatrix,1))
      real*8, dimension(:) :: dummyvector(size(aamatrix,1))
      integer :: m 
      integer :: n 
      integer :: p 
      real*8  :: alpha 
      real*8  :: beta 
      integer :: i
!
!     ! check length of matrix and vector match
      if (size(aamatrix,2) .ne. size(bbvector,1)) then 
         write(*,*) ' *** WARNING size mismatch between matrix and vector ***'
         stop
      end if 
!
      m = size(aamatrix,1)
      n = int(1) 
      p = size(bbvector,1)
      alpha = real(1.0,8)
      beta = real(0.0,8)
!
      do i = 1, size(bbvector,1)
         dummyvector(i) = bbvector(i)%p
!!         write(*,'(a15,e15.7)') 'dummyflow:', dummyvector(i)
      end do
!
      call dgemm('N',            &
                 'N',            &
                  m,             &
                  n,             &
                  p,             &
                  alpha,	      &
                  aamatrix,      &
                  m,             &
                  dummyvector,   &
                  p,             &
                  beta,          &
                  ccvector,      &
                  m)             
!
      return
      end function matrixvectorpointer
!
! *** wrapper for dgesv to solve matrix vector system ***
!
      function solvematrixvector(amat,bvec) result(xvec) 
      implicit none
      real*8, dimension(:,:), intent(in) :: amat
      real*8, dimension(:), intent(in)   :: bvec
      real*8, dimension(:,:) :: amatcopy(size(amat,1),size(amat,2))
      real*8, dimension(:) :: bveccopy(size(bvec,1))
      real*8, dimension(:) :: xvec(size(amat,1))
      integer :: n
      integer :: nrhs 
      integer :: lda 
      integer :: ipiv(size(amat,1)) 
      integer :: ldb       
      integer :: info
!
!     ! check length of matrix and vector match
      if (size(amat,1) .ne. size(bvec,1)) then 
         write(*,*) ' *** WARNING size mismatch between matrix and vector ***'
         stop
      elseif (size(amat,1) .ne. size(amat,2)) then
         write(*,*) ' *** WARNING non square matrix ***'
         stop
      end if 
!
      n = size(amat,1)
      nrhs = int(1)   !! vector
      lda = n 
      ldb = n
      amatcopy = amat
      bveccopy = bvec
!
      call dgesv(n,           &
                 nrhs,        &
                 amatcopy,    &
                 lda,         &
                 ipiv,        &
                 bveccopy,    &
                 ldb,         &
                 info)
!      
      if (info .ne. int(0)) then 
         write(*,'(a,i3)') ' *** WARNING dgesv error, info = ',info
         stop
      end if
!
      xvec = bveccopy
!
      return
      end function
!
!c
! *** wrapper for dgesv to solve matrix matrix system ***
!
      function solvematrixmatrix(aamat,bbmat) result(xxmat) 
      implicit none
      real*8, dimension(:,:), intent(in) :: aamat
      real*8, dimension(:,:), intent(in) :: bbmat
      real*8, dimension(:,:) :: aamatcopy(size(aamat,1),size(aamat,2))
      real*8, dimension(:,:) :: bbmatcopy(size(bbmat,1),size(bbmat,2))
      real*8, dimension(:) :: xxmat(size(aamat,1),size(bbmat,2))
      integer :: n
      integer :: nrhs 
      integer :: lda 
      integer :: ipiv(size(aamat,1)) 
      integer :: ldb       
      integer :: info
!
!     ! check length of matrix and vector match
      if (size(aamat,1) .ne. size(bbmat,1)) then 
         write(*,*) ' *** WARNING size mismatch between matrix and matrix ***'
         stop
      elseif (size(aamat,1) .ne. size(aamat,2)) then
         write(*,*) ' *** WARNING non square matrix ***'
         stop
      end if 
!
      n = size(aamat,1)
      nrhs = size(bbmat,2)   !! vector
      lda = n 
      ldb = n
      aamatcopy = aamat
      bbmatcopy = bbmat
!
      call dgesv(n,           &
                 nrhs,        &
                 aamatcopy,   &
                 lda,         &
                 ipiv,        &
                 bbmatcopy,   &
                 ldb,         &
                 info)
!      
      if (info .ne. int(0)) then 
         write(*,'(a,i3)') ' *** WARNING dgesv error, info = ',info
         stop
      end if
!
      xxmat = bbmatcopy
!
      return
      end function


      ! wrapper to call the lapack subroutine dgetri to invert a matrix      
      ! see http://www.netlib.org/lapack/double/dgetri.f for dgetri reference,
      ! also, http://www.netlib.no/netlib/lapack/double/dgetrf.f
      function invertSquareMatrix(aamat) result(aainverse_mat)

         implicit none

         include "mpif.h"

         integer :: irank, ierr

         real*8, dimension(:,:), intent(in) :: aamat
         real*8, allocatable, dimension(:,:) :: aainverse_mat

         integer :: rows, columns
         integer, allocatable :: IPIV(:) ! will hold row permutation data from dgetrf for passing to dgetri
         integer :: INFO
         real*8, allocatable :: WORK(:) ! just some memory for dgetri to use as a scratchpad

         rows = size(aamat,1)
         columns = size(aamat,2)
         allocate(aainverse_mat(rows,columns))

         allocate(IPIV(rows))
         allocate(WORK(rows*columns))
         INFO=0

         aainverse_mat = aamat

         call MPI_COMM_RANK(MPI_COMM_WORLD,irank,ierr)

         ! First, we need to get the P*L*U factorisation of aamat:
         ! The factorisation comes out in aainverse_mat, presumably with L and U occupying the two triangles of the one square matrix.
         call dgetrf(rows,columns,aainverse_mat,rows,IPIV,INFO)
         call MPI_BARRIER(MPI_COMM_WORLD,ierr)
         if (INFO .ne. int(0)) then
            write(*,*) 'Error in invertSquareMatrix: dgetrf call. Code: ', INFO
            if (irank .eq. int(0)) then
               write(*,*) 'Problematic computed factorisation is:'
               write(*,*) aainverse_mat
               write(*,*) 'and there are', rows, 'rows,'
               write(*,*) 'and', columns, 'columns.'
               call flush(6)
               write(*,*) 'matrix is'
               write(*,*) aamat
            end if
            call flush(6)
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
            stop
         end if

         ! rows*columns gives the size of the scratchpad workspace WORK.
         ! Stackoverflow says this is the optimal size for inverting the whole matrix.
         ! Note that if INFO=0 on exit, the WORK(1) will contain the optimal WORK size.
         call dgetri(columns, aainverse_mat, rows, IPIV, WORK, rows*columns, INFO)
         if (INFO .ne. int(0)) then
            write(*,*) 'Error in invertSquareMatrix: dgetri call. Code: ', INFO
            stop
         end if

      end function invertSquareMatrix
!
      end module multidomain   
!



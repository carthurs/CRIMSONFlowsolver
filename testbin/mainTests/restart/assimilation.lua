----------------------------------- GLOBAL -----------------------------------

 Delta_t_model = 1
 final_time_model = 2200
 -- Saving period.
 Nskip_save = 1100

 problem_directory = "/home/nxiao/simulations/carotid_segment_mgs_verycoarse/def_pulse/"

 output_mode = "text"
 output_directory = problem_directory .. "result/"
 output_mode_scalar = "text"

----------------------------------- MODEL ------------------------------------
simvascular_model = {
    -- Choose what goes into the reduced state
    state_reduced_has_wall_parameters = 0;
    state_reduced_has_coupled_parameters = 1;

    RCR_parameters_info = {
        estimate_resistance = 1,
        estimate_compliance = 1,
        estimate_prox_resistance = 1,
        estimate_pstates = 0,
        estimate_pout = 0,
        resistance_included = {1},
        compliance_included = {1},
        prox_resistance_included = {1},
    },

    Heart_parameters_info = {
        estimate_emax = 0,
        estimate_tmax = 0,
        estimate_trel = 0,
        estimate_vlv = 0,
    }, 

    error_statistics = {
        -- Diagonal of state error variance
        state_error_variance = {0.2, 0.2, 0.2},
    },
}

-------------------------------- OBSERVATION ---------------------------------
observation = {
   -- Do we use simulation restarts as data?
   -- Deprecated
   use_restarts = 0,

   -- Path to the stored data files.
   data_directory = problem_directory .. "truth/",

   -- Times at which the data are available
   initial_time = 0,

   final_time = 2200,

   data_period = 2200, --set to the same value as final_time if data is not periodic

   Nskip = 10,

   -- Simple nodal observations
   execute_nodal_observations_ = 0,  

   -- Distance observations
   execute_distance_observations_ = 0,

   -- Cross-sectional flow observation
   Nobservation_flow = 1,

   Nobservation_avgpressure = 1,

   csobs_origins = {0, 0, 0,
                    0, 0, 0},

   csobs_normals = {0, 0, 1,
                    0, 0, 1},

   csobs_radii = {4, 4},

   error = {
      -- diagonal entries of observation error variance
      variance = 1,
      variance_nodal = 1e-2,
      variance_dist = 1e-3,
      variance_avgpress = 1.7774e+06,
      variance_flow = 1e6,
   },
}

----------------------------------- METHOD -----------------------------------
-- Simulation with assimilation using ROUKF.
reduced_order_unscented_kalman_filter = {
   data_assimilation = {
      analyze_first_step = false,
      with_resampling = false,
      -- Indicates how R is stored: "matrix", "matrix_inverse".
      observation_error_variance = "matrix_inverse",
      -- Should the computations be carried out on observations (false) or
      -- directly on innovations (true)?
      with_innovation = true,
      with_innovation_computation = true

   },

   sigma_point = {
      -- Choice of sigma-points: "canonical", "star" or "simplex".
      type = "simplex"
   },

   display = {
      show_iteration = true,
      show_time = true
   },

   output_saver = {
       variable_list = {"state_forecast", "state_analysis"},
       file = output_directory .. "roukf-%{name}.%{extension}",
       time = "step " .. Delta_t_model * Nskip_save .. " 1.0e-6",
       mode = output_mode,
       mode_scalar = output_mode_scalar
   },

   output = {
       configuration = output_directory .. "roukf.lua",
       log = output_directory .. "roukf_%{rank}.log"
   },

   mpi = {
      algorithm = 0,
      master_process_contribution = 1.0
   }
}
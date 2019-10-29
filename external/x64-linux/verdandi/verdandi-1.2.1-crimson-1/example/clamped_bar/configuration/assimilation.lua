----------------------------------- GLOBAL -----------------------------------


Delta_t_clamped_bar = 0.01
final_time_clamped_bar = 10.0
-- Saving period.
Nskip_save = 1

observation_file = "result/truth-state_forecast.bin"

output_mode = "text"
output_directory = "result/"
output_mode_scalar = "text"


----------------------------------- MODEL ------------------------------------


dofile("configuration/clamped_bar.lua")

-- In order to demonstrate the assimilation, errors are introduced in the
-- model.
clamped_bar.physics.theta_force = {1., 1.}


-------------------------------- OBSERVATION ---------------------------------


dofile("configuration/observation.lua")


----------------------------------- METHOD -----------------------------------


-- Simulation with assimilation using optimal interpolation.
optimal_interpolation = {

   -- Computation mode for BLUE: "vector" or "matrix".
   BLUE_computation = "vector",

   data_assimilation = {

      analyze_first_step = true,

   },

   display = {

      show_iteration = false,
      show_time = true
   },

   output_saver = {

      variable_list = {"state_forecast", "state_analysis"},
      file = output_directory .. "oi-%{name}.%{extension}",
      time = "step " .. Delta_t_clamped_bar * Nskip_save .. " 1.e-6",
      mode = output_mode,
      mode_scalar = output_mode_scalar

   },

   output = {

      configuration = output_directory .. "oi.lua",
      log = output_directory .. "oi.log"

   }

}


-- Simulation with assimilation using EKF.
extended_kalman_filter = {

   -- Computation mode for BLUE: "vector" or "matrix".
   BLUE_computation = "matrix",
   -- Computation mode for covariance: "vector" or "matrix".
   covariance_computation = "vector",

   data_assimilation = {

      analyze_first_step = false,

   },

   display = {

      show_iteration = false,
      show_time = true
   },

   output_saver = {

      variable_list = {"state_forecast", "state_analysis"},
      file = output_directory .. "ekf-%{name}.%{extension}",
      time = "step " .. Delta_t_clamped_bar * Nskip_save .. " 1.e-6",
      mode = output_mode,
      mode_scalar = output_mode_scalar

   },

   output = {

     configuration = output_directory .. "ekf.lua",
     log = output_directory .. "ekf.log"

  }

}


-- Simulation with assimilation using UKF.
unscented_kalman_filter = {

   data_assimilation = {

      analyze_first_step = false

   },

   sigma_point = {

      -- Choice of sigma-points: "canonical", "star" or "simplex".
      type = "simplex"

   },

   display = {

      show_iteration = false,
      show_time = true

   },

   output_saver = {

      variable_list = {"state_forecast", "state_analysis"},
      file = output_directory .. "ukf-%{name}.%{extension}",
      time = "step " .. Delta_t_clamped_bar * Nskip_save .. " 1.e-6",
      mode = output_mode,
      mode_scalar = output_mode_scalar

   },

   output = {

     configuration = output_directory .. "ukf.lua",
     log = output_directory .. "ukf.log"

  }

}


-- Simulation with assimilation using ROUKF.
reduced_order_unscented_kalman_filter = {

   data_assimilation = {

      analyze_first_step = false,
      with_resampling = false,
      -- Indicates how R is stored: "matrix", "matrix_inverse".
      observation_error_variance = "matrix_inverse",
      -- Should the computations be carried out on observations (false) or
      -- directly on innovations (true)?
      with_innovation = true

   },

   sigma_point = {

      -- Choice of sigma-points: "canonical", "star" or "simplex".
      type = "simplex"

   },

   display = {

      show_iteration = false,
      show_time = true

   },

   output_saver = {

      variable_list = {"state_forecast", "state_analysis"},
      file = output_directory .. "roukf-%{name}.%{extension}",
      time = "step " .. Delta_t_clamped_bar * Nskip_save .. " 1.e-6",
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


-- Simulation with assimilation using ROEKF.
reduced_order_extended_kalman_filter = {

   data_assimilation = {

      analyze_first_step = false,

   },

   display = {

      show_iteration = false,
      show_time = true

   },

   output_saver = {

      variable_list = {"state_forecast", "state_analysis"},
      file = output_directory .. "roekf-%{name}.%{extension}",
      time = "step " .. Delta_t_clamped_bar * Nskip_save .. " 1.e-6",
      mode = output_mode,
      mode_scalar = output_mode_scalar

   },

   output = {

     configuration = output_directory .. "roekf.lua",
     log = output_directory .. "roekf_%{rank}.log"

   },

}


-- Simulation with assimilation using 4D-VAR.
four_dimensional_variational = {

   nlopt = {

      -- Optimization algorithm (LD_VAR1, LD_LBFGS, LD_SLSQP, LD_MMA,
      -- LD_TNEWTON ...).
      algorithm = "LD_VAR1",

      -- If you do not want to use a particular tolerance termination, you can
      -- just set that tolerance to zero and it will be ignored.
      -- Relative tolerance on the optimization parameters.
      parameter_tolerance = 1.e-4,
      -- Relative tolerance on the cost function.
      cost_function_tolerance = 1.e-4,
      -- Maximum number of function evaluations. Criterion is disabled
      -- if the value is non-positive.
      Niteration_max = -1

   },

   trajectory_manager = {

      -- Checkpoint recording mode (memory or disk).
      checkpoint_recording_mode = "disk",
      -- Checkpoint recording file
      checkpoint_recording_file = "checkpoint.bin",
      -- Trajectory recording mode.
      trajectory_recording_mode = "memory",
      -- Trajectory recording file.
      trajectory_recording_file = "trajectory.bin",
      -- Recording period.
      Nskip_step = 10,


   },

   display = {

      show_optimization_iteration = true,
      show_optimized_parameter = false,
      show_iteration = false,
      show_time = true

   },

   output_saver = {

      variable_list = {"state_forecast", "state_analysis"},
      file = output_directory .. "4dvar-%{name}.%{extension}",
      time = "step " .. Delta_t_clamped_bar * Nskip_save .. " 1.e-6",
      mode = output_mode,
      mode_scalar = output_mode_scalar

   },

   output = {

     configuration = output_directory .. "4dvar.lua",
     log = output_directory .. "4dvar.log"

   },

}


-- Forward simulation.
forward = {

   display = {

      show_iteration = false,
      show_time = true

   },

   output_saver = {

      variable_list = {"state_forecast"},
      file = output_directory .. "forward-%{name}.%{extension}",
      time = "step " .. Delta_t_clamped_bar * Nskip_save .. " 1.e-6",
      mode = output_mode,
      mode_scalar = output_mode_scalar

   },

   output = {

      configuration = output_directory .. "forward.lua",
      log = output_directory .. "forward.log"

   }

}

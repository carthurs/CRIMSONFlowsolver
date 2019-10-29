----------------------------------- GLOBAL -----------------------------------


Delta_t_petsc_clamped_bar = 0.01
final_time_petsc_clamped_bar = 10
-- Saving period.
Nskip_save = 1

observation_file = "result/truth-state_forecast.bin"

output_mode = "text"
output_directory = "result/"
output_mode_scalar = "text"


----------------------------------- MODEL ------------------------------------


dofile("configuration/petsc_clamped_bar.lua")

-- In order to demonstrate the assimilation, errors are introduced in the
-- model.
petsc_clamped_bar.physics.theta_force = {1., 1.}


-------------------------------- OBSERVATION ---------------------------------


dofile("configuration/observation.lua")


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
      time = "step " .. Delta_t_petsc_clamped_bar * Nskip_save .. " 1.e-6",
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


-- Forward simulation.
forward = {

   display = {

      show_iteration = false,
      show_time = true

   },

   output_saver = {

      variable_list = {"state_forecast"},
      file = output_directory .. "forward-%{name}.%{extension}",
      time = "step " .. Delta_t_petsc_clamped_bar * Nskip_save .. " 1.e-6",
      mode = output_mode,
      mode_scalar = output_mode_scalar

   },

   output = {

      configuration = output_directory .. "forward.lua",
      log = output_directory .. "forward.log"

   }

}



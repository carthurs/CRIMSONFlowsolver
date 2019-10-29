----------------------------------- GLOBAL -----------------------------------


Delta_t_petsc_clamped_bar = 0.01
final_time_petsc_clamped_bar = 10
-- Saving period.
Nskip_save = 1

output_directory = "result/"
output_mode = "binary"
output_mode_scalar = "text"


----------------------------------- MODEL ------------------------------------


dofile("configuration/petsc_clamped_bar.lua")


----------------------------------- METHOD -----------------------------------


forward = {

   output_saver = {

      variable_list = {"state_forecast"},
      file = output_directory .. "truth-%{name}.%{extension}",
      time = "step " .. Delta_t_petsc_clamped_bar * Nskip_save .. " 1.e-6",
      mode = output_mode,
      mode_scalar = output_mode_scalar

   },

   display = {

      show_iteration = false,
      show_time = true

   },

   output = {

      configuration = output_directory .. "truth.lua",
      log = output_directory .. "truth.log"

   }

}

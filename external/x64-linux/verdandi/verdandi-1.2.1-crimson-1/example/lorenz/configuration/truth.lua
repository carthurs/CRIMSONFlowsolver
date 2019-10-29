----------------------------------- GLOBAL -----------------------------------


output_directory = "result/"


----------------------------------- MODEL ------------------------------------


dofile("configuration/lorenz.lua")


----------------------------------- METHOD -----------------------------------


-- Forward simulation.
forward = {

   display = {

      show_iteration = false,
      show_time = false

   },

   output_saver = {

      variable_list = {"state_forecast"},
      file = output_directory .. "truth-%{name}.%{extension}",
      mode = "binary",
      mode_scalar = "text"

   },

   output = {

      configuration = output_directory .. "truth.lua",
      log = output_directory .. "truth.log"

   }

}

----------------------------------- GLOBAL -----------------------------------


Delta_t_model = 0.0015
Nskip_save = 100

output_directory = "result/"


----------------------------------- MODEL ------------------------------------


dofile("configuration/quadratic_model.lua")


python_model = {

   module = "QuadraticModel",
   directory = "../../model/",
   class_name = "QuadraticModel"

}


----------------------------------- METHOD -----------------------------------


-- Forward simulation.
forward = {

   display = {

      show_iteration = false,
      show_time = true

   },

   output_saver = {

      variable_list = {"state_forecast"},
      file = output_directory .. "truth-%{name}.%{extension}",
      time = "step " .. Delta_t_model * Nskip_save .. " 1.e-6"

   },

   output = {

      configuration = output_directory .. "truth.lua",
      log = output_directory .. "truth.log"

   }

}

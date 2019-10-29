perturbation_manager = {

   newran = {

      -- Choose between
      -- * "time" for machine time,
      -- * "number" for a fixed seed value between 0 and 1, and
      -- * "directory" to use a Newran directory which stores the seeds.
      seed_type = "number",

      -- If 'seed_type' is set to "number", choose a number in ]0, 1[.
      seed_number = 0.4,

      -- If 'seed_type' is set to "directory", put the path to the Newran
      -- directory. It must end with a slash.
      seed_directory = os.getenv("HOME") .. "/.newran/"

   },

   output = {

      configuration = output_directory .. "perturbation.lua",
      log = output_directory .. "perturbation.log"

   }

}

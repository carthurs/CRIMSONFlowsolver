-------------------------------- OBSERVATION ---------------------------------


observation = {

   -- Path to the file storing the observations.
   file = observation_file,
   -- How are defined the observations? If the type is "observation, only
   -- observations are stored in the file. If the type is "state", the whole
   -- model state is stored.
   type = "state",
   -- Is the period with which observations are available constant?
   Delta_t_constant = true,
   -- If the period with which observations are available non constant
   -- one should define the observation time file.
   observation_time_file = "result/time.dat",
   -- Period with which observations are available.
   Delta_t = Delta_t_model * Nskip_save,
   -- Period with which available observations are actually assimilated.
   Nskip = 1,
   -- Duration during which observations are assimilated.
   final_time = 2.3,

   aggregator = {

      type = "step",
      width_left = Delta_t_model * Nskip_save / 2.,
      width_right = Delta_t_model * Nskip_save / 2.,

      -- If the type is "triangle", the triangles widths may be the same for
      -- all observations ("constant") or not ("per-observation").
      width_property = "constant",

      -- If the triangles widths are not constant, or in the case of
      -- "interpolation", one should define an observation interval. It is
      -- assumed that the observations outside this interval have no
      -- contribution.
      width_left_upper_bound = 1.,
      width_right_upper_bound = 1.,

      -- If the value is true, each observation can be used only one time.
      discard_observation = false

   },

   error = {

      -- Variance of observational errors.
      variance = 10.

   },

   operator = {

      -- Is the operator a scaled identity matrix?
      scaled_identity = true,
      -- If so, put the diagonal value:
      diagonal_value = 1.,
      -- Otherwise, the operator value (file name or table):
      value = {}

   },

}

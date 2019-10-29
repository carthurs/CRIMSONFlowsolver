-------------------------------- OBSERVATION ---------------------------------


observation = {

   -- Path to the file storing the observations.
   file = observation_file,
   -- How are defined the observations? If the type is "observation, only
   -- observations are stored in the file. If the type is "state", the whole
   -- model state is stored.
   type = "state",
   -- Period with which observations are available.
   Delta_t = Delta_t_petsc_clamped_bar * Nskip_save,
   -- Period with which available observations are actually assimilated.
   Nskip = 1,
   -- Duration during which observations are assimilated.
   final_time = final_time_petsc_clamped_bar,

   -- In case of triangles widths defined in a file.
   width_file = "configuration/width.bin",

   aggregator = {

      -- The interpolation type may be "step", "triangle" or "interpolation".
      type = "step",
      width_left = 0.005,
      width_right = 0.005,

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
      discard_observation = true

   },

   error = {

      -- Variance of observational errors.
      variance = 1.e-4

   },

   operator = {

      -- The first "Nobservation" values of the state vector are observed.
      -- "Nobservation", must be lower than the size of the local state
      -- on processor 0.
      Nobservation = 5,
      -- If so, put the diagonal value:
      diagonal_value = 1.,
      -- Otherwise, the operator value (file name or table):
      value = {}


   },

   location = {

      observation_location = {1, 0}

   }

}

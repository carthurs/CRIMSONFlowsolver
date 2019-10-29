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
   Delta_t = Delta_t_shallow_water * Nskip_save,
   -- Period with which available observations are actually assimilated.
   Nskip = 1,
   -- Duration during which observations are assimilated.
   final_time = final_time_shallow_water,

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
      variance = 100.

   },

   operator = {

      -- Is the operator a scaled identity matrix?
      scaled_identity = false,
      -- If so, put the diagonal value:
      diagonal_value = 1.,
      -- Otherwise, the operator value (file name or table):
      value = {}

   },

   location = {

      observation_location = {80, 0,
                              81, 0,
                              82, 0}

   }

}

for i = 1, 100 do
   for j = 1, 3 do
     observation.operator.value[(i-1)*3 + j] = 0.
   end
end

observation.operator.value[80] = 1.
observation.operator.value[100 + 81] = 1.
observation.operator.value[200 + 82] = 1.


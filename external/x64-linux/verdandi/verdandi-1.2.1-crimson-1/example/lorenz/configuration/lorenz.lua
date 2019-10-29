----------------------------------- MODEL ------------------------------------


lorenz = {

   parameter = {

      Prandtl = 10.,
      Rayleigh = 28.,
      b = 2.6666667

   },

   initial_condition = {

      X = 0.1,
      Y = 0.1,
      Z = 0.1

   },

   time = {

      Delta_t = .001,
      initial_time = 0.,
      final_time = 10.

   },

   output_saver = {

      variable_list = {"X", "Y", "Z"},
      file = "result/%{name}.bin",
      mode_scalar = "binary"

   }

}

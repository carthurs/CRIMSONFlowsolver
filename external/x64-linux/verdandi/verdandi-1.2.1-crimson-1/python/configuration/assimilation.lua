dofile("../example/quadratic_model/configuration/assimilation.lua")
forward.display.show_iteration = false
forward.display.show_time = false
optimal_interpolation.display.show_iteration = false
optimal_interpolation.display.show_time = false
python_model.directory = "../model/"
python_observation_manager.directory = "../observation_manager/"
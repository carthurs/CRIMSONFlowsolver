----------------------------------- MODEL ------------------------------------


petsc_clamped_bar = {

    -- Components of the state vector.
    state = {"displacement", "velocity", "theta_force"},
    -- Reduced state parameters.
    reduced_state = {"theta_force"},

    domain = {

        -- Time step.
        Delta_t = Delta_t_petsc_clamped_bar,
        -- Simulation time.
        final_time = final_time_petsc_clamped_bar,

        bar_length = 1.,
        Nx = 10

    },

    physics = {

        Young_modulus = 1.,
        mass_density = 1.,
        theta_force = {1.5, 1.7},
        theta_stiffness = {1., 1.},
        theta_mass = {1., 1., 1.},
        theta_damp = {1.},
        alpha = 0.001,
        beta = 0.001
    },

    error_statistics = {

        -- Diagonal value of "B".
        state_error_variance = 100,

    },

    output_saver = {

        variable_list = {"disp_0", "velo_0"},
        file = output_directory .. "/%{name}.bin",
        time = "step " .. Delta_t_petsc_clamped_bar * Nskip_save .. " 1.e-6",
        mode = output_mode,
        mode_scalar = output_mode_scalar

   }

}

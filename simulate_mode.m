function simulation = simulate_mode(params, mode)
    time = 0:params.simulation.dt:params.simulation.duration;
    sample_count = numel(time);
    simulation.time = time;
    simulation.mode = mode;
    simulation.position_eci = zeros(3, sample_count);
    simulation.velocity_eci = zeros(3, sample_count);
    simulation.eclipse = false(1, sample_count);
    simulation.true_attitude_body_eci = zeros(3, 3, sample_count);
    simulation.true_rpy = zeros(3, sample_count);
    simulation.true_rate_body = zeros(3, sample_count);

    true_rate_body = zeros(3, 1);
    current_gyro_bias = zeros(3, 1);

    if strcmpi(mode, 'fine')
        simulation.star_tracker_measurement = zeros(3, 3, sample_count);
        simulation.star_tracker_error = zeros(1, sample_count);
        simulation.gyro_measurement = zeros(3, sample_count);
        simulation.gyro_bias = zeros(3, sample_count);
    elseif strcmpi(mode, 'coarse')
        simulation.sun_true_body = zeros(3, sample_count);
        simulation.sun_measured_body = nan(3, sample_count);
        simulation.magnetic_true_body = zeros(3, sample_count);
        simulation.magnetic_measured_body = zeros(3, sample_count);
        simulation.sun_valid = false(1, sample_count);
        simulation.gyro_measurement = zeros(3, sample_count);
        simulation.gyro_bias = zeros(3, sample_count);
    else
        error('Unsupported simulation mode: %s', mode);
    end

    for sample = 1:sample_count
        [simulation.position_eci(:, sample), simulation.velocity_eci(:, sample)] = ...
            propagate_orbit(time(sample), params.orbit);
        sun_eci = sun_vector(time(sample), params);
        attitude_body_eci = attitude_truth(...
            simulation.position_eci(:, sample), simulation.velocity_eci(:, sample), ...
            sun_eci);
        simulation.true_attitude_body_eci(:, :, sample) = attitude_body_eci;
        simulation.true_rpy(:, sample) = dcm_to_rpy(attitude_body_eci);

        if sample > 1
            previous_attitude = simulation.true_attitude_body_eci(:, :, sample - 1);
            C_dot = (attitude_body_eci - previous_attitude) / params.simulation.dt;
            omega_skew = -C_dot * attitude_body_eci';
            true_rate_body = [omega_skew(3, 2); ...
                omega_skew(1, 3); omega_skew(2, 1)];
            simulation.true_rate_body(:, sample) = true_rate_body;
        end

        simulation.eclipse(sample) = check_eclipse(...
            simulation.position_eci(:, sample), sun_eci);

        if strcmpi(mode, 'fine')
            [simulation.star_tracker_measurement(:, :, sample), ...
                simulation.star_tracker_error(sample)] = ...
                simulate_star_tracker(attitude_body_eci, params);
            [simulation.gyro_measurement(:, sample), current_gyro_bias] = ...
                simulate_gyro(true_rate_body, current_gyro_bias, ...
                params.simulation.dt, params);
            simulation.gyro_bias(:, sample) = current_gyro_bias;
        else
            [simulation.sun_measured_body(:, sample), ...
                simulation.sun_true_body(:, sample), simulation.sun_valid(sample)] = ...
                simulate_sun_sensor(sun_eci, attitude_body_eci, ...
                params, simulation.eclipse(sample));
            [simulation.magnetic_measured_body(:, sample), ...
                simulation.magnetic_true_body(:, sample)] = ...
                simulate_magnetometer(simulation.position_eci(:, sample), ...
                attitude_body_eci, params);
            [simulation.gyro_measurement(:, sample), current_gyro_bias] = ...
                simulate_gyro(true_rate_body, current_gyro_bias, ...
                params.simulation.dt, params);
            simulation.gyro_bias(:, sample) = current_gyro_bias;
        end
    end
end
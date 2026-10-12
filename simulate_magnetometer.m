function [measurement_body, true_body] = simulate_magnetometer(...
        position_eci, attitude_body_eci, params)
    magnetic_field_eci = get_magnetic_field(position_eci, params);
    true_body = attitude_body_eci * magnetic_field_eci;
    true_sensor = params.magnetometer.mounting_dcm * true_body;
    noisy_sensor = params.magnetometer.misalignment * true_sensor + ...
        params.magnetometer.bias_vector + ...
        chol(params.magnetometer.covariance, 'lower') * randn(3, 1);
    measurement_body = params.magnetometer.mounting_dcm' * noisy_sensor;
end